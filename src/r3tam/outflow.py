# -*- coding: utf-8 -*-
"""
Created on Fri Feb 18 14:36:06 2022

@author: jgilbert
"""

from r3tam.constants import ACREtoM2, DELT_SEC, FTtoM, AFtoM3, FT3toM3, CFStoKG, Le, \
    cw, CMStoCFS, CFStoAFD, TOL_AF, AFDtoCMS

import numpy as np

import os
import datetime as dt
import pandas as pnd

import copy
from scipy.optimize import linprog

def init_check_target(resmod, outvec):
    
    # placeholder monthly temperature targets for reservoir TCD gate ops
    if 'tempTarg' not in resmod.Outflow.ColumnMap:
        # default temp targets
        resmod.defaultTTarg = True
        try:
            targtemp = resmod.TempTargetData[list(resmod.TempTargetData.keys())[0]]['default']
        except:
            raise Exception("No default temperature target data found - check input *yaml file\nIs there a 'Temperature_Target' section?")
            exit(1)
    else:
        resmod.defaultTTarg = False
        targVar = resmod.Outflow.ColumnMap['tempTarg']
        targVals = resmod.Outflow.DataFrame[targVar]
        numNan = np.sum(np.isnan(targVals))
        if numNan == len(targVals): # then a column existed but no values were provided - go with default
            resmod.defaultTTarg = False
            targtemp = resmod.TempTargetData[list(resmod.TempTargetData.keys())[0]]['default']
        else:
            targtemp = targVals
            
    if type(targtemp)!=pnd.Series:
        targtemp = pnd.Series(list(targtemp.values()), index = pnd.date_range(f"{resmod.SimDates[0].year}-{min(targtemp)}-01",
                                                               f"{resmod.SimDates[0].year}-{max(targtemp)}-31", freq='MS'))
        new = pnd.date_range(targtemp.index[-1], periods=1, freq='M')
        targtemp = targtemp.append(pnd.Series([targtemp[-1]], index=new))
        targtemp = targtemp.resample('d').apply('first').ffill()
    resmod.TempTargetData[list(resmod.TempTargetData.keys())[0]]['Target'] = targtemp
    
    
def outflow3(self, qout, gateDict, log=False, debug=False):
    '''
      trying a more hydraulics-based approach
    '''
    logfp = r'D:\02_Projects\SacTemp\ShastaTempModel\AdpatingUSBRsimpleModel\log\outflow.v3.log'
    if log:
        lf = open(logfp,'a')
        #lf.write("Run date/time: %s\n" %dt.datetime.now().strftime('%Y-%m-%d %H:%M'))
        #lf.write("\n========================================================\n\n")
        lf.write('\n\nDate: %s\n' %self.TimeStep.strftime('%Y-%m-%d'))
        lf.write('Gate      NumberOpen      ResLayer      Outflow      Temp\n')
    
    
    fclimit = 18500. * CFStoAFD  # assuming flows above 18500 cfs are passed through river release/spillways reather than TCD
    qoutTCD = min(qout, fclimit)
    qoutFC = qout-qoutTCD
    totOutFCQ = 0.
    totOutFCE = 0.
    totOutQ = 0.
    totOutE = 0.

    topLyr = 0
    
    tw_elev = 575. # ????: Assuming a default tailwater elevation of 575 ft, based on ratign curve diagram, Chart A-3 Army Corps of Engineers, 1977 (Shasta Standing Operating Procedures)
    for lj in range(0, self.nLyrs):
        if self.Layers[lj].Vol >0:
            topLyr = lj
    topElev = self.Layers[topLyr].Vol/self.Layers[topLyr].MaxVol * \
                (self.Layers[topLyr].MaxElev - self.Layers[topLyr].MinElev) + \
                self.Layers[topLyr].MinElev
    
    # get heads at centerline of each gate level
    glvlHdVels = {}
    for g in gateDict:
        ctrel = self.Outlets[g].CtrElevFt
        B = self.Outlets[g].GateHeight_m
        dlel = self.WSE - tw_elev   # from W2, withdrawal.f90 line 50; see W2 v3 User's Manual section on gates and control file specification
        dlel_m = dlel*FTtoM
        qrel = (self.Outlets[g].A1GT*(dlel_m**self.Outlets[g].B1GT)*B**self.Outlets[g].G1GT)*gateDict[g]
        vel = None
        #dlel = max(0., topElev - ctrel)  # if topElev is below the ctrelevation of the gate level, this will be 0
        #vel = np.sqrt(2.*32.2 * hd)
        #qrel = vel * gateDict[g]  # using number of gates open as a proxy for area in flow calc
        glvlHdVels[g] = [dlel, vel, qrel]
        if debug:
            print("Gate Level %s - Hd:  %0.2f  Vel:  %0.2f ft/s   Qrel:  %0.2f"  %(g, hd, vel, qrel))
    
    gfracs = {} # fraction of qoutTCD that should come from each level
    
    totqrel = 0.
    for g in glvlHdVels:
        totqrel += glvlHdVels[g][2]
    for g in glvlHdVels:
        #if glvlHdVels[g][2]>0:
        gfracs[g] = glvlHdVels[g][2]/totqrel
    
    # do flood control accounting - assuming flood control water comes
    # from the top layer at the top layer temperature
    tmpFCvol = qoutFC
    thisLyr = topLyr
    while tmpFCvol >0.:
        ioutfc = min(self.Layers[thisLyr].Vol, tmpFCvol)
        self.Layers[thisLyr].Vol -= ioutfc
        totOutFCQ += ioutfc
        totOutFCE += ioutfc * self.Layers[thisLyr].Temp
        self.Layers[thisLyr].TotE = max(0., self.Layers[thisLyr].Vol*self.Layers[thisLyr].Temp)
        tmpFCvol -= ioutfc  # decrement the remaining flood control amount
        thisLyr -= 1     # decrement the layer from which to remove the flood control volume
    
    for g in gfracs:
        tmpqoutTCD = gfracs[g]*qoutTCD
        #print(gfracs[g], qoutTCD, tmpqoutTCD)
        resLyrsLtoS = sorted(self.Outlets[g].ResLayers, reverse=True)
        # now distribute to layers
        n = 0
        numNonEmpty = sum(n+1 for i in resLyrsLtoS if self.Layers[i].Vol>0.)
        if numNonEmpty>0:
            lyrFrac = 1/numNonEmpty
        else:
            print("Layers empty for prescribed releases at %s" %self.Outlets[g].Name)
            lyrFrac = 0.

        remainq = 0.
        for l in resLyrsLtoS:
            if self.Layers[l].Vol>0.:
                lyrq = lyrFrac * tmpqoutTCD + remainq
                remainq = 0.
                ioutq = min(lyrq, self.Layers[l].Vol)
                if debug:
                    print("Gate: %s   Layer: %s    Volume: %0.2f   Temp: %0.2f" %(self.Outlets[g].Name, l, ioutq, self.Layers[l].Temp))
                remainq += max(0., lyrq - ioutq)
                self.Layers[l].Vol -= ioutq
                self.Layers[l].TotE = max(0., self.Layers[l].Vol*self.Layers[l].Temp)
                if ~np.isnan(self.Layers[l].Temp):
                    outTotE = ioutq * self.Layers[l].Temp
                #outletQTemp[gi].append([ioutq, outTotE])
                totOutQ += ioutq
                totOutE += outTotE
                if log:
                    lf.write('%d-%s      %s       %d       %0.2f        %0.2f\n' %(g, self.Outlets[g].Name, gateDict[g], l, ioutq, self.Layers[l].Temp ))
        
        if remainq > 0.:
            print("WARNING: Excess outflow not allocated at outlet %s" %(g))
            if log:
                lf.write("WARNING: Excess outflow of %0.2f AF not allocated at outlet %s\n\n" %(remainq, g))
            
    totOutQ = totOutQ + totOutFCQ
    totOutE = totOutE + totOutFCE
    # compile the outlfow and weighted outflow temp
    if abs(totOutQ - qout)> TOL_AF:
        print("WARNING: Calculated outlfow doesn't match prescribed outflow!")
        print("WARNING:\t\t Calculated: %0.2f AF   vs   Prescribed:  %0.2f AF"  %(totOutQ, qout))
        if log:
            lf.write("WARNING: Calculated outlfow doesn't match prescribed outflow!\n")
            lf.write("WARNING: Calculated: %0.2f AF   vs   Prescribed:  %0.2f AF\n\n"  %(totOutQ, qout))
            
    outTemp = totOutE/totOutQ
    if log:
        lf.write("Total outflow: %0.2f AF\n" %totOutQ)
        lf.write("Outflow temperature: %0.2f deg F\n\n" %outTemp)
    if log:
        lf.close()
    return([totOutQ, outTemp])
        
def calc_leakage(self, qout, topElev, gateDict):
    '''
        based on RMA (2003) method and the RAFT (2018) implementation
         - calculate leakage based on penstock flows using factors for 
           each leakge level
           
         - 2020-06-24 rev: include gate opening configuration so that leakage
                         at a level coinciding with an open gate is turned
                         off
    '''
    fclimit = 18750. * CFStoAFD  # assuming flows above 18500 cfs are passed through river release/spillways reather than TCD
    qoutPEN = min(qout, fclimit)  # amount not through river release/spillway = penstock/through TCD
    qoutFC = qout-qoutPEN
    maxElev = self.Layers[self.nLyrs-1].MaxElev  #999. #
    qleakDict = {}
    qleakTot = 0.
    opengates = [g for g,v in gateDict.items() if v>0]
    sumfacs = sum([self.Leakages[x].LeakageFactor for x in self.Leakages])
    for qli, qlv in self.Leakages.items():
        
        if qlv.OutletExclude not in opengates:
            #elev_wght = 2.*(topElev - qlv.BotElevFt)/(maxElev - qlv.BotElevFt)
            elev_wght = 0.33*np.exp(np.sqrt(2*32.2*(topElev-qlv.BotElevFt))/np.sqrt(2.*32.2*(maxElev - qlv.BotElevFt)))
            tmpfac = min(0.33,(qlv.LeakageFactor * elev_wght))
            #elev_wght = 1.
            if qlv.BotElevFt < topElev:  # only leakage if there's water at that elevation
                if qli == 6:
                    elev_wght = elev_wght*1.
                qleaki = min(qoutPEN, qoutPEN * tmpfac)  # qlv.ElevWeight
                qleakDict[qli] = qleaki
                qleakTot += qleaki
    if qleakTot > qoutPEN:
        print("Calculated leakage is greater than penstock release - somethings not right!!!!")
        return
    return([qoutPEN, qoutFC, qleakDict, qleakTot])


def calc_leakageMDWE(resmod, qoutPEN, topElev, gateDict, bypass_frac=0):
    '''
        based on Mike Deas' draft 2020 Shasta Reservoir modeling report for 
        the revised Shasta/Keswick CE-QUAL-W2 models
        
        features: 
            * Leakage fraction depends on reservoir level
            * No leakage through levels where a gate is open
            * 
    
    '''

    
    ## the next step is then calculating what the total and zonal leakage fractions
    # are based on 1) penstock (TCD) flow, 2) water surface elevation, and
    # 3) which gates are open/closed
    
    # first check if we are using pre- or post 2010 formulation
    if resmod.TimeStepDate >=dt.datetime(2010, 10, 1, 0, 0):
        pre2010 = False
    else:
        pre2010 = True
    
    mFrac = 0.2
    relFracs = {}
    if pre2010:
        if topElev >= 1000. : # water surface elevation greater than or equal to 1000 ft
            tmpFracX = mFrac
            relFracs[1] = .1309
            relFracs[2] = .0805 + ((5-gateDict[2])/5.)*0.1165  # zone 2 leakage affected by whether middle gates are open
            relFracs[3] = 0.0934 + ((2-gateDict[0])/2.)*0.0331  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0103 + ((5 - gateDict[1])/5.)*0.1001 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5]= 0.0384
            relFracs[6] = 0.0179
            relFracs[7] = 0.3112
            relFracs[8] = 0.0677
            
        elif (topElev < 1000.) & (topElev >= 945.) : # water surface elevation between 945 and 1000.
            tmpFracX = mFrac - (mFrac*0.1309)*(1-(topElev-945.)/(1000.-945.))
            relFracs[1] = 0.1309* (topElev-945.)/(1000.-945.)
            relFracs[2] = 0.0805 + ((5-gateDict[2])/5.)*0.1165  # zone 2 leakage affected by whether middle gates are open
            relFracs[3] = 0.0934 + ((2-gateDict[0])/2.)*0.0331  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0103 + ((5 - gateDict[1])/5.)*0.1001 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0384
            relFracs[6] = 0.0179
            relFracs[7] = 0.3112
            relFracs[8] = 0.0677                
            
        elif (topElev < 945.) & (topElev >= 900.) : # water surface elevation between 945 and 1000.
            tmpFracX = mFrac - (mFrac*0.1309) - (mFrac*0.1970)*(1-(topElev-900.)/(945.-900.))
            tmpFracX = max(0., tmpFracX) # just in case the calc above comes up with a negative number
            relFracs[1] = 0.
            relFracs[2] = (0.0805 + ((5-gateDict[2])/5.)*0.1165)*((topElev-900)/(945.-900.))  # zone 2 leakage affected by whether middle gates are open
            relFracs[3] = 0.0934 + ((2-gateDict[0])/2.)*0.0331  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0103 + ((5 - gateDict[1])/5.)*0.1001 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0384
            relFracs[6] = 0.0179
            relFracs[7] = 0.3112
            relFracs[8] = 0.0677
        
        elif (topElev < 900.) & (topElev >= 831.) : # water surface elevation between 945 and 1000.
            tmpFracX = mFrac - (mFrac*0.1309) - (mFrac*0.1970) - \
                        (mFrac*0.1265)*(1-(topElev-831.)/(900.-831.))
            tmpFracX = max(0., tmpFracX) # just in case the calc above comes up with a negative number
            relFracs[1] = 0.
            relFracs[2] = 0.
            relFracs[3] = (0.0934 + ((2-gateDict[0])/2.)*0.0331)*(topElev-831.)/(900.-831.)  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0103 + ((5 - gateDict[1])/5.)*0.1001 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0384
            relFracs[6] = 0.0179
            relFracs[7] = 0.3112
            relFracs[8] = 0.0677

        else: # this case was not in Mike Deas' report - including just in case
            tmpFracX = 0.10  # approximate extrapolation from zone above
            relFracs[1] = 0.
            relFracs[2] = 0.
            relFracs[3] = 0.
            relFracs[4] = 0.0103 + ((5 - gateDict[1])/5.)*0.1001 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0384
            relFracs[6] = 0.0179
            relFracs[7] = 0.3112
            relFracs[8] = 0.0677
                        
    else:  # post-2010
        if topElev >= 1000. : # water surface elevation greater than or equal to 1000 ft
            #print("post-2010, topElev>1000 (%0.2f ft)" %topElev)
            tmpFracX = mFrac
            relFracs[1] = 0.163
            relFracs[2] = 0. # zone 2 leakage assumed zero after 2010 modifications
            relFracs[3] = 0.1163 + ((2-gateDict[0])/2.)*0.0412  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0128 + ((5 - gateDict[1])/5.)*0.1247 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0478
            relFracs[6] = 0.0223
            relFracs[7] = 0.3876
            relFracs[8] = 0.0844
            
        elif (topElev < 1000.) & (topElev >= 945.) : # water surface elevation between 945 and 1000.
            
            #print("post-2010, topElev<1000 & topElev>945 (%0.2f ft)" %topElev)
            
            tmpFracX = mFrac - (mFrac*0.1630)*(1-(topElev-945.)/(1000.-945.))
            tmpFracX = max(0., tmpFracX) # just in case the calc above comes up with a negative number
            relFracs[1] = 0.1630*((topElev-945.)/(1000.-945.))
            relFracs[2] = 0. # zone 2 leakage assumed zero after 2010 modifications
            relFracs[3] = 0.1163 + ((2-gateDict[0])/2.)*0.0412  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0128 + ((5 - gateDict[1])/5.)*0.1247 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0478
            relFracs[6] = 0.0223
            relFracs[7] = 0.3876
            relFracs[8] = 0.0844              
            
        elif (topElev < 945.) & (topElev >= 900.) : # water surface elevation between 945 and 1000.
            
            #print("post-2010, topElev<945 & topElev>900 (%0.2f ft)" %topElev)
            
            tmpFracX = mFrac - (mFrac*0.1630) - (mFrac*0.)*(1-(topElev-945.)/(945.-900.))
            tmpFracX = max(0., tmpFracX) # just in case the calc above comes up with a negative number
            
            relFracs[1] = 0.
            relFracs[2] = 0.
            relFracs[3] = 0.1163 + ((2-gateDict[0])/2.)*0.0412  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0128 + ((5 - gateDict[1])/5.)*0.1247 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0478
            relFracs[6] = 0.0223
            relFracs[7] = 0.3876
            relFracs[8] = 0.0844 
        
        elif (topElev < 900.) & (topElev >= 831.) : # water surface elevation between 945 and 1000.
            
            #print("post-2010, topElev<900 & topElev>831 (%0.2f ft)" %topElev)
            
            tmpFracX = mFrac - (mFrac*0.1630) - (mFrac*0.) - \
                        (mFrac*0.1575)*(1-(topElev-831.)/(900.-831.))
            tmpFracX = max(0., tmpFracX) # just in case the calc above comes up with a negative number
            
            relFracs[1] = 0.
            relFracs[2] = 0.
            relFracs[3] = 0.1163 + ((2-gateDict[0])/2.)*0.0412*(topElev-831.)/(900.-831.)  # zone 3 leakage affected by side gates
            relFracs[4] = 0.0128 + ((5 - gateDict[1])/5.)*0.1247 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0478
            relFracs[6] = 0.0223
            relFracs[7] = 0.3876
            relFracs[8] = 0.0844 

        else: # this case was not in Mike Deas' report - including just in case
            
            #print("post-2010, topElev<831  (%0.2f ft)" %topElev)
            tmpFracX = 0.10 # approximate extrapolation from zone above
            relFracs[1] = 0.
            relFracs[2] = 0.
            relFracs[3] = 0.
            relFracs[4] = 0.0128 + ((5 - gateDict[1])/5.)*0.1247 # zone 4 leakage affected by lower (PRG) gates
            relFracs[5] = 0.0478
            relFracs[6] = 0.0223
            relFracs[7] = 0.3876
            relFracs[8] = 0.0844 

    ### we've now calculated the relative fraction of leakage by zone [z#_rfrac] and an
    # initial estimate of total leakage (as a fraction of TCD outflow) [tmpFracX]
    # now we can calculate the actual leakages
       
    relFracSum = sum(relFracs.values())
    fracX = tmpFracX * relFracSum
    
    # total leakage
    qLeakTot = fracX * qoutPEN

    
    if qLeakTot > qoutPEN:
        print("Calculated leakage is greater than total penstock flow..something is wrong!")
        return
    
    # now apply the leakages by zone
    qleakDict = {}
    if len(resmod.Leakages) < 8:
        print("Fewer leakage zones defined in input file (%d) than needed (8) for special Shasta TCD leakage calcs" %len(resmod.Leakages))
        return
    for qli, qlv in resmod.Leakages.items():
        if relFracs[qli] > 0:
            qleakDict[qli] = (relFracs[qli]/relFracSum)*qLeakTot
    
    ## debugging
    totcalclkg = 0
    for qli, qlv in qleakDict.items():
        totcalclkg += qlv
    if abs(totcalclkg-qLeakTot) > 0.01:
        print('Calculated leakage: %0.2f' %totcalclkg)
        print('Intended leakage: %0.2f' %qLeakTot)
        print('Penstock flow: %0.2f' %qoutPEN)
   
    
    return([qleakDict, qLeakTot])
    

def alloc_to_outlet_types(resmod, qout, bypassFrac=0.):
    """
    Distributes outflow to available outlet types
    
    Function to distribute a given outflow volume for this time step (assuming
    units of acre-feet per day for outflow) among the three outlet types: penstocks,
    river outlets, and spillway. If the `bypassFrac` is greater than 0, that portion
    of the release will be assigned to the river outlets and the remainder is
    assigned to the penstocks. Any amount that exceeds the capacity or prescribed
    release through any of these outlet types is assumed to exit the dam as spill.
    
    Parameters
    ----------
    resmod : ResTemp reservoir object
        DESCRIPTION.
    qout : float
        Total outflow volume to be released in this time step [acre-feet/day]
    bypassFrac : float, default=0
        Fraction of total outflow that should bypass the penstocks/hydropower 
        and therefore be allocated to river outlet releases
        
    Returns
    -------
    release_allocs : list of float
        List containing the outflow volumes [acre-feet/day] allocated to (in order)
        the penstocks, river outlets, and spillway

    """

    if resmod.PenstockLimit_Units.upper() == 'CFS':  # TODO: Check units - convert to SI system
        penLim = resmod.PenstockLimit*CFStoAFD
    else:
        penLim = resmod.PenstockLimit

    if resmod.RiverOutleLimit_Units.upper() == 'CFS':  # TODO: Check units - convert to SI system
        rivLim = resmod.RiverOutletLimit*CFStoAFD
    else:
        rivLim = resmod.RiverOutletLimit
        
    # is water high enough to access both river and penstock outlets?
    if resmod.WSE > resmod.PenstockElevation + 10: # ???? Trying to ensure an outflow for all time steps
        pen_access = True
    else:
        pen_access = False
        
    if resmod.WSE > resmod.RiverOutlets[min(resmod.RiverOutlets)].MinElev:
        riv_access = True
    else:
        riv_access = False

    if not pen_access and riv_access:
        # then no access to penstocks - have to use river outlets
        bypassFrac = 1 
    
    
    if bypassFrac > 0.:  # then we're intentionally running water through river outlet instead of penstocks/TCD
        qoutPEN = min((1-bypassFrac)*qout, penLim)

        qoutRIV = min(bypassFrac * qout, rivLim)

    else:
        # is release (qout) greater than the penstock capacity?
        qoutPEN = min(qout, penLim)
        #qoutRIV = 0
        qoutRIV = min(qout-qoutPEN, rivLim)

    # the rest that can't go through the penstocks or river outlets is assuemd to go over the spillway
    qoutSPILL = max(0., qout - qoutPEN - qoutRIV)

    release_allocs = [qoutPEN, qoutRIV, qoutSPILL]
    return(release_allocs)

        
def distPointSinks(resmod,topElev, outletID, nPtSinks,**kwargs):
    
    if 'top_elev_override' in kwargs:
        thisTopElev = kwargs['top_elev_override']
        thisBotElev = thisTopElev - resmod.Outlets[outletID].GateHeight_m/FTtoM
    else:
        thisTopElev = resmod.Outlets[outletID].TopElevFt        
        thisBotElev = resmod.Outlets[outletID].BotElevFt
    
    ptSinkElevs = []
    ctrreldist_ft = min((topElev-thisBotElev)/2.,(thisTopElev-thisBotElev)/2.)
    ctrelev = thisBotElev+ctrreldist_ft
    if nPtSinks == 1:
        # put point sink in center of gate opening range
        hd = topElev - ctrelev
        ptSinkElevs.append([ctrelev,np.sqrt(2*32.2*hd)])
    elif nPtSinks ==2:
        # put one point sink in the middle, one at the bottom
        ptSinkElevs.append([ctrelev, np.sqrt(2*32.2*(topElev - ctrelev))])
        ptSinkElevs.append([thisBotElev, np.sqrt(2*32.2*(topElev-thisBotElev))])
    elif nPtSinks >=3:
        # this distributes the sinks evenly across the vertical height of the 
        # opening, with the lowest sink assigned at the bottom elevation (thus
        # the denominator term: `(nPtSinks-1)`)
        delz = (ctrreldist_ft*2.)/(nPtSinks-1) 
        #delz = (ctrreldist_ft*2.)/(nPtSinks-1)/2 # test: put point sinks only in the bottom half of the opening

        for npt in range(nPtSinks):
            elev = thisBotElev+npt*delz
            # if outletID==3:
            #     print(f"Point sink elev: {elev} ft")
            if elev<topElev: # only assign sink if below water surface elevation
                ptSinkElevs.append([elev, np.sqrt(2*32.2*(topElev-elev))])    
                    
    tmpsum = np.sum([j for i,j in ptSinkElevs])
    fracs = [j/tmpsum for i,j in ptSinkElevs] 
    
    return([ptSinkElevs,fracs])

def tcd_alloc(resmod, gateDict, qoutPEN, qleakTot,  nPtSinks, lp_opt=False,
              outqLyr={}, qleakDict={}, rivDict={},
              **kwargs):
    '''
        given total penstock flow and leakage, allocate
        non-leakage outflow to gate levels based on prescribed operations
            -based on RAFT approach (Daniels et al 2018, NMFS Tech Memo)
            
        outqLyrDist={}, qleakDict={}, rivDict={} are only used if lp_opt=True
        so that the LP opt process can estiamte outflow temperatures
    '''
    # set the debug release level (0,1,2)
    debug = resmod.Debug['Release']
    
    qoutTCD = qoutPEN - qleakTot
    gateQ = {k:[] for k,v in resmod.Outlets.items()}  # dictionary of outflows at each gate
    totGatesOpen = sum(gateDict.values())
    
    maxGateLev = max(gateDict)
    minGateLev = min(gateDict)
    
    # is this the first use of side gates in isolation?
    if resmod.TimeStep==0:
        prevGate = gateDict
    else:
        #prevGate = resmod.PrevGateDict # sgateDict #
        prevGate = copy.deepcopy(resmod.Operations.GateOps[resmod.SimDates[resmod.TimeStep-1]])
                                 
    if (gateDict[0]>0): # and (sum([prevGate[3],prevGate[2],prevGate[1]])<=0):
        # this is the first time the side gates are used 
        resmod.Operations.sidegate_only_cntr += 1
        adj_sdg_ptsinks = True
        
        # side gate top point sink elevation adjustment
        sdgadjlevel_max = 832 #824 #832
        Ladj = sdgadjlevel_max - resmod.Outlets[0].TopElevFt
        kadj = 0.4 # steepness
        x0 = 45 #90 # midpoint (days) of logistic curve
        
        v = Ladj/(1+np.exp(-1*kadj*(resmod.Operations.sidegate_only_cntr-x0)))
        sdgadjlevel = sdgadjlevel_max-v
        
        # sdgadjtime = 15
        # sdgidx = min(resmod.Operations.sidegate_only_cntr, sdgadjtime)
        # sdgadjlevel = min(resmod.WSE, (sdgadjlevel_max +
        #                                (resmod.Outlets[0].TopElevFt-sdgadjlevel_max)
        #                                *sdgidx/sdgadjtime))
        if debug>1:
            print(f"\nAdj sidegate level: {sdgadjlevel} ft")
            print(f"Gate dict: 3: {gateDict[3]} | 2: {gateDict[2]} | 1: {gateDict[1]} | 0: {gateDict[0]}")
    else:
        resmod.Operations.sidegate_only_cntr = 0
        adj_sdg_ptsinks = False
        

    # testing lowering of middle gate point sinks after being open for long 
    # periods (>60 days)
    if resmod.TimeStepDate.timetuple().tm_yday==1 or resmod._time==0: # first day of the year, reset counter
        resmod.Operations.midgatedays = 0
    else:
        if gateDict[2] >0:
            resmod.Operations.midgatedays += 1
        else:
            resmod.Operations.midgatedays = 0
    #gate_history = resmod.Operations.GateOps
    if resmod.Operations.midgatedays > 60:
        adj_mid_ptsinks = True
        midadjlevel = resmod.Outlets[2].CtrElevFt
    else:
        adj_mid_ptsinks = False

    
    # get heads at centerline of each gate level
    #tw_elev = 586. # feet; from x-sect drawing of penstock 4, p78 in pdf of (United States Army Corps of Engineers, 1977)
    glvlHdVels = {}
    for g in gateDict:
        
        if gateDict[g]<=0:
            glvlHdVels[g] = [ 0., 0., 0.]
            continue
        
        if g>minGateLev: # not at bottom gate level
            g_below = g-1
        else:
            g_below = -1
        
        gate_flow_onoff = 1
        if resmod.CheckMinHead:
            if g_below<0 or gateDict[g_below] ==0:  # no gates open on level below
                if resmod.WSE >= resmod.Outlets[g].MinHead + resmod.Outlets[g].BotElevFt:
                    gate_flow_onoff = 1  # no gates open below, but meeting head requirement
                else:
                    gate_flow_onoff = 0 # no gates open below and NOT meeting head reqmt - flow should be 0 thru gate
                    # update gate dict for this time step accordingly
                    gateDict[g] = 0
            else:
                if gateDict[g_below]>0: # there's a gate open below
                    if resmod.WSE >= resmod.Outlets[g].MinHead + resmod.Outlets[g].BotElevFt: 
                        # there's water covering outlet gate and there's a gate open below
                        gate_flow_onoff = 1  
                    else:
                        gate_flow_onoff = 0 
                        # update gate dict for this time step accordingly
                        gateDict[g] = 0
                else:
                    gate_flow_onoff = 0  # there are no gates open below
                    # update gate dict for this time step accordingly
                    gateDict[g] = 0
                
            
        
        #ctrel = resmod.Outlets[g].CtrElevFt
        B = min((resmod.WSE-resmod.Outlets[g].BotElevFt)*FTtoM, resmod.Outlets[g].GateHeight_m)
        B = max(B, 0.0)
        #dlel = resmod.WSE - 819.5 #2021-11-23: changing to penstock centerline elev per new reading of W2 logic; tw_elev   # from W2, withdrawal.f90 line 50; see W2 v3 User's Manual section on gates and control file specification
        #dlel_m = dlel*FTtoM
        dlel_m = max(0.0, resmod.WSE*FTtoM - max(resmod.Outlets[g].BotElevFt*FTtoM,resmod.PenstockElevation_m))
        #dlel_m = resmod.WSE*FTtoM - resmod.PenstockElevation_m
        
        qrel = gate_flow_onoff*(resmod.Outlets[g].A1GT*
                                (dlel_m**(resmod.Outlets[g].B1GT))*
                                B**resmod.Outlets[g].G1GT)*gateDict[g]
        #qrel = min(0, qrel)
        #print(f'qrel = {qrel}\n')
        
        # reduces fraction coming from an upper gate level if WSE is below top of gate
        if (B < resmod.Outlets[g].GateHeight_m) and (g_below>-1): # gate level is only partially submerged
            qrel_frac = gateDict[g]/(gateDict[g]+gateDict[g_below])     
            qrel = qrel*qrel_frac
        
        vel = None
        #dlel = max(0., topElev - ctrel)  # if topElev is below the ctrelevation of the gate level, this will be 0
        #vel = np.sqrt(2.*32.2 * hd)
        #qrel = vel * gateDict[g]  # using number of gates open as a proxy for area in flow calc
        glvlHdVels[g] = [dlel_m, vel, qrel]
        if debug>1:
            print("Gate Level %s - Hd:  %0.2f meters Vel: -none- ft/s   Qrel:  %0.2f"  %(g, dlel_m, qrel))

    # calculate the fraction of flow through each gate level
    tot_qrel = sum([v[2] for v in glvlHdVels.values()])
    for g in gateDict:
        #print(f"gate {g}: {gateDict[g]}")
        if tot_qrel==0:
            glvlHdVels[g].append(0)
        else:
            glvlHdVels[g].append(glvlHdVels[g][2]/tot_qrel)
    
    # find which selective withdrawal ports, g, are currently open 
    openg = [g for g in glvlHdVels if glvlHdVels[g][3]>0]
    openg.sort(reverse=False)

    len_openg = len(openg) # number of open gate levels
    
    fracs_by_level = {}
    
    if len_openg==1: # just flow through a single gate level
        if adj_sdg_ptsinks:
            if debug>1:
                print("Adjusting side gate point sink level")
            # make an adjustment to top point sink elevation for side gate
            ptSinkElevs, fracs = distPointSinks(resmod, resmod.WSE, openg[0], 
                                                nPtSinks, top_elev_override=sdgadjlevel)
        elif adj_mid_ptsinks:
            if debug>1:
                print(f"{resmod.TimeStepDate.strftime('%Y-%m-%d: ')}: Adjusting middle gate point sink elevations")
            ptSinkElevs, fracs = distPointSinks(resmod,  resmod.WSE, openg[0],
                                               nPtSinks, top_elev_override=midadjlevel)
        else:
            ptSinkElevs, fracs = distPointSinks(resmod,  resmod.WSE, openg[0], nPtSinks)

        gateQ[openg[0]] = [[qoutTCD*f for f in fracs], ptSinkElevs]
        
        fracs_by_level[openg[0]] = [1] #one gate level open - all flow through this
        
    elif len_openg==2: # two gate levels are open
        if openg[0]==2 and adj_mid_ptsinks:
            ptSinkElevsUpp, fracsUpp = distPointSinks(resmod, resmod.WSE, openg[0], 
                                                      nPtSinks, top_elev_override=midadjlevel)
            ptSinkElevsLow, fracsLow = distPointSinks(resmod,  resmod.WSE, openg[1], nPtSinks)
        elif openg[1]==2 and adj_mid_ptsinks:
            ptSinkElevsUpp, fracsUpp = distPointSinks(resmod,  resmod.WSE, openg[0], 
                                                      nPtSinks)
            ptSinkElevsLow, fracsLow = distPointSinks(resmod,  resmod.WSE, openg[1], 
                                                      nPtSinks,
                                                      top_elev_override=midadjlevel)
        
        elif openg[1]==0 and adj_sdg_ptsinks:
            ptSinkElevsUpp, fracsUpp = distPointSinks(resmod,  resmod.WSE, openg[0], 
                                                      nPtSinks)
            ptSinkElevsLow, fracsLow = distPointSinks(resmod,  resmod.WSE, openg[1], 
                                                      nPtSinks,
                                                      top_elev_override=sdgadjlevel)
        
        else:
            ptSinkElevsUpp, fracsUpp = distPointSinks(resmod, resmod.WSE, openg[0], nPtSinks)
            ptSinkElevsLow, fracsLow = distPointSinks(resmod, resmod.WSE, openg[1], nPtSinks)
        
        quppfrac= glvlHdVels[openg[0]][3] # with 2 levels open, relative flow fraction from 
                                    # top can be used to distribute between the two levels
        
        gateQ[openg[0]] = [[qoutTCD*f*quppfrac for f in fracsUpp], ptSinkElevsUpp]
        gateQ[openg[1]] = [[qoutTCD*f*(1-quppfrac) for f in fracsLow], ptSinkElevsLow]
        
        fracs_by_level[openg[0]] = [fracsUpp] #fraction of flow through top level
        fracs_by_level[openg[1]] = [fracsLow] # fraction of flow through lower level
    
    elif len_openg==3: # rare instance where 3 levels are open at once
        ptSinkElevsUpp, fracsUpp = distPointSinks(resmod,  resmod.WSE, openg[0], nPtSinks)
        ptSinkElevsMid, fracsMid = distPointSinks(resmod,  resmod.WSE, openg[1], nPtSinks)
        ptSinkElevsLow, fracsLow = distPointSinks(resmod,  resmod.WSE, openg[2], nPtSinks)
        
        gateQ[openg[0]] = [[qoutTCD*f*glvlHdVels[openg[0]][3] for f in fracsUpp], 
                       ptSinkElevsUpp]
        gateQ[openg[1]] = [[qoutTCD*f*glvlHdVels[openg[1]][3] for f in fracsMid], 
                       ptSinkElevsMid]
        gateQ[openg[2]] = [[qoutTCD*f*glvlHdVels[openg[2]][3] for f in fracsLow], 
                       ptSinkElevsLow]
        
        fracs_by_level[openg[0]] = [fracsUpp] #fraction of flow through top level
        fracs_by_level[openg[1]] = [fracsMid] # fraction of flow through middle level
        fracs_by_level[openg[2]] = [fracsLow]
        
    elif len_openg==0:
        print(f"\n--{resmod.TimeStepDate}: No open TCD gate this time step")
    else:
        # if condition not met by one of the above, just use the first (lowest)
        # gate level
        ptSinkElevs, fracs = distPointSinks(resmod, resmod.WSE, openg[0], nPtSinks)
        gateQ[openg[0]] = [[qoutTCD*f for f in fracs], ptSinkElevs] 
        
        fracs_by_level[openg[0]] = [1]
        
    # revise the distribution of flows according to a LP-opt solution
    # considering open gates, point sink elevations, and target temp
    
    if lp_opt==True:
        if 'temp_targ' in kwargs:
            temp_targ = kwargs['temp_targ']
            if temp_targ >=99:
                return(gateQ)
            #print(f'\n\t\t...target is: {temp_targ}')
        else:
            raise BaseException("Temperature target ('temp_targ') value not provided to 'tcd_alloc' function")
            
        
        # set the constraints on the outlets - assume if the gate is open, it can have
        # some minimum fraction (0.01) up to 1 of the flow
        target_bounds = [] # constraints for each gate & pt sink
        gate_temps = [] # temperatures at each gate & pt sink
        ps_elevs = [] # point sink elevations
        ps_temps = [] # temperatures at each point sink
        
        # get current temperature profile for easy interpolation lookup
        modElevs = [v.CtrElev for l, v in resmod.Layers.items()]
        these_temps = []
        for l, v in resmod.Layers.items():
            these_temps.append(v.Temp)
        temp_prof = pnd.DataFrame(index=modElevs, data=these_temps)
        
        #print(f"There are {len(gateQ[3][0])} point sinks at the upper gate")
        for g in gateQ:
            if len(gateQ[g])>0: # if this gate level is not null (i.e. it is active)
                for ps in gateQ[g][1]: # iterate through point sink elevatoins
                    target_bounds.append([0.01, 1])
                    ps_elevs.append(ps[0])
            # if gateDict[g] ==0:
            #     target_bounds[g] == [0., 0.0]
            # else:
            #     target_bounds[g] = [0.01,1]
            
        #print("Point sink elevations:")
        #print(ps_elevs)
        #print(modElevs)
        #print(these_temps)
        ps_temps = np.interp(ps_elevs, modElevs, these_temps) #.reverse() needed because np.interp needs values in ascending order
        
        #print("Interpolated temps:")
        #print(ps_temps)

        # build the inputs to the lp problem
        lhs_eq = [[0]+len(ps_temps)*[1],  
                   [0]+[pt for pt in ps_temps]]
        # print("A_EQ:")
        # print(lhs_eq)
        
        rhs = [1.0, temp_targ]
        clist = [-1*temp_targ]+[pt for pt in ps_temps]
        # print("c: ")
        # print(clist)
        
        bounds = [[0,1]] + target_bounds
        
        # print("Bounds:")
        # print(bounds)
        
        # print("RHS")
        # print(rhs)
        
        # do the LP optimzation
        opt_ = linprog(c=clist,A_eq=lhs_eq, b_eq=rhs,
                       bounds=bounds, method='highs-ds' )
        
        #print(resmod.TimeStepDate.isoformat() + "_" + str(opt_.x)) #@debug

        if type(opt_.x)!=type(None): # there is a solution given constraints
            counter=1
            for g in gateQ:
                if len(gateQ[g])>0:
                    new_fracs = []
                    for ps in gateQ[g][1]:
                        new_fracs.append(opt_.x[counter]) # skip 0-index value in opt_.x data because that's a constant
                        counter+=1
                    gateQ[g][0] = [qoutTCD*nf for nf in new_fracs]

        else: # if opt_x==None, then there is no solution given constraints - either too warm or too cool to achieve
        # in this case, default to existing distribution (i.e. don't need to do anythin)
            #pass
            
            #print(f"{resmod.TimeStepDate.isoformat()}\t-- No initial lp-opt solution found...") #@debug
            for g in gateQ:
                if len(gateQ[g])>0: # if this gate level is not null (i.e. it is active)
                    cumul_frac = 0.0
                    fracs = []
                    num_pt_sinks = len(gateQ[g][1])
                    #fracs = [1/num_pt_sinks for nps in range(num_pt_sinks)]
                    if rhs[1]>max(ps_temps):  #target is warmer than warmest point sink
                        # put max amount of water through highest point sinks
                        for ips, ps in enumerate(gateQ[g][1]): # iterate through point sink elevatoins
                            
                            if ips < num_pt_sinks-1:
                                fracs.append(0.01) # TODO: Change so it can keep track of user-specified different minfracs
                                cumul_frac += 0.01
                                #print(f'{resmod.TimeStepDate.isoformat()}\tAdjusting WARM gate {g} point sink {ips} at elevation {ps[0]} to 0.01...') #@debug
                            else:
                                fracs.append( 1 - cumul_frac)
                                #print(f'{resmod.TimeStepDate.isoformat()}\tAdjusting WARM gate {g} point sink {ips} at elevation {ps[0]} to {1 - cumul_frac}...') #@debug
                    
                    else: # target is colder than coldest point sink
                        # put max amount through lowest point sinks
                       
                        for ips, ps in enumerate(gateQ[g][1]): # iterate through point sink elevatoins
                            #print(f't\tAdjusting COOL gate {g} point sink {ips} at elevation {ps}...')
                            if ips == 0: #num_pt_sinks-1:
                                fracs.append(1-(num_pt_sinks-1)*0.01) # TODO: Change so it can keep track of user-specified different minfracs
                                
                                #print(f'{resmod.TimeStepDate.isoformat()}\tAdjusting COOL gate {g} point sink {ips} at elevation {ps[0]} to {1-(num_pt_sinks-1)*0.01}...') #@debug
                            else:
                                fracs.append( 0.01) 
                                #print(f'{resmod.TimeStepDate.isoformat()}\tAdjusting COOL gate {g} point sink {ips} at elevation {ps[0]} to 0.01...')  #@debug
                                
                            # make sure total fraction = 1  # TODO - make this better
                            tot_fracs = sum(fracs)
                            if abs(tot_fracs - 1.0) > 0.00001:
                                fracs = [f* 1.0/tot_fracs for f in fracs]
                                
                    
                    #print(glvlHdVels[g][3])  #@debug
                    # print("\n\t\t---Adj Est Fracs:")  #@debug
                    # print(fracs)  #@debug
                    gateQ[g][0] = [qoutTCD*f*glvlHdVels[g][3] for f in fracs]
                    
        # check the distribution fo flows by layer
        outqLyrDist_ = get_outflows_by_layer(resmod, gateQ, outqLyr,
                                            qleakDict, rivDict, dz=0)
        
        # get the temperature of just leakage components so this can be used
        # to adjust TCD target
        totQ_ = 0.
        tmplq_ = {}
        for l in resmod.Layers:
            tmplq_[l]=0.
        for rri, rrv in outqLyrDist_.items():
            if rri=='Release_by_Outlet_Type':
                continue
            #print(rri, rrv)
            if 'l' in rri: # a leakage component
                for rli, rlv in rrv[0].items(): # loop through the layer assignments for this outlet flow
                    if np.isnan(rlv):
                        print("Calc'd outflow for %s is %0.2f" %(rri, rlv))
                        return
                    tmplq_[rli] += rlv
                    totQ_ += rlv
        [lkg_totOutQ, lkg_outTemp, lkg_totOutE] = outflow_dist_forGateSelect(resmod,tmplq_)
        
        # update target to take into account leakage
        temp_targ_TCD = (temp_targ*(lkg_totOutQ + qoutTCD)-(lkg_totOutE))/qoutTCD
        #print(f"\n\t\t\tOld and New temp target for TCD: {temp_targ} -> {temp_targ_TCD}") #@debug
        
        rhs2 = [1.0, temp_targ_TCD]
        
        new_ps_temps = []
        openg = [g for g in gateQ if gateQ[g]!=[]]            
        for og in openg: #iterate through each open gate level
            ptsnk_cntr=1
            for qstr, [outletEl,sqrtgh] in list(zip(gateQ[og][0],gateQ[og][1])): # iterate through point sinks for each gate
            
                outByLyr = outqLyrDist_[str(og)+'g'+str(ptsnk_cntr)][0] #get the layer withdrawals for this point sink
                
                totq_ = 0.0
                tote_ = 0.0
                for l in outByLyr: # iterate through layers in wd for this point sink
                    totq_ += outByLyr[l]
                    tote_ += outByLyr[l]*resmod.Layers[l].Temp
                ps_temp_ = tote_/totq_
                new_ps_temps.append(ps_temp_)
                ptsnk_cntr+=1
        
        # update the LHS of the lp input
        lhs_eq2 = [[0]+len(new_ps_temps)*[1],  
                   [0]+[pt for pt in new_ps_temps]]

        #print(new_ps_temps)
        
        # 20250423 - attempting to deal with missed targets - adjust target
        # so that it's within the range of point sinks
        # if temp_targ_TCD < min(new_ps_temps): # target is too cold
        #     temp_targ_TCD2 = min(new_ps_temps) + 0.15
        # elif temp_targ_TCD > max(new_ps_temps): # target is too warm
        #     temp_targ_TCD2 = max(new_ps_temps) - 0.15
        # else:
        #     temp_targ_TCD2 = temp_targ_TCD
            
        # rhs2 = [1.0, temp_targ_TCD2]
        
        opt_2 = linprog(c=clist,A_eq=lhs_eq2, b_eq=rhs2,
                       bounds=bounds, method='highs-ds' )        
        
        if type(opt_2.x)!=type(None): # there is a solution given constraints
            #print("\n\t\t***Updated point sink flow allocations*****")    #@debug
            counter=1
            for g in gateQ:
                if len(gateQ[g])>0:
                    new_fracs = []
                    for ps in gateQ[g][1]:
                        new_fracs.append(opt_2.x[counter]) # skip 0-index value in opt_.x data because that's a constant
                        counter+=1
                    gateQ[g][0] = [qoutTCD*nf for nf in new_fracs]

        else:
            
            #print(f"Could not find lp-opt solution on {resmod.TimeStepDate.isoformat()}") #@debug
            #print(new_ps_temps)
            
            if resmod.TimeStepDate-dt.timedelta(1) in resmod.PointSinkFracs:

                # use the last set of gateQ instead
                prev_psfrac = resmod.PointSinkFracs[resmod.TimeStepDate-dt.timedelta(1)]
                
                # check that prev gate levels and this gate levels match
                prev_gate_levs = [g for g in prev_psfrac if len(prev_psfrac[g])>0]
                this_gate_levs =  [g for g in gateQ if gateQ[g]!=[]] #[g for g in gateQ[g] if len(gateQ[g])>0]
                #print(gateQ)
                if sorted(prev_gate_levs)==sorted(this_gate_levs):
                    for g in gateQ:
                        if g in this_gate_levs:
                            new_fracs = []
                            counter=0
                            for ps in gateQ[g][1]:
                                #print(f'counter: {counter}')
                                #print(prev_psfrac)
                                if len(prev_psfrac[g])<counter:
                                    new_fracs.append(prev_psfrac[g][counter]) 
                                counter+=1
                                
                            gateQ[g][0] = [qoutTCD*nf for nf in new_fracs]
            else:
                # nothing else we can do for now
                # default back to existing gate ops
                pass
        
        # save the point sink levels fractions for debugging
        save_ps_fracs = {k:[] for k,v in resmod.Outlets.items()}
        for g in gateQ:
            if len(gateQ[g])>0:
                save_ps_fracs[g] = [ff/qoutTCD for ff in gateQ[g][0]]
        resmod.PointSinkFracs[resmod.TimeStepDate] = save_ps_fracs
        
    return(gateQ)
    
# def tcd_alloc_lp(resmod, target_temp, gateDict, qoutPEN, qleakTot, nPtSinks):            
#     """
#     Calculate distribution to TCD gates/ports according to a tailbay target 
#     temperature, using an linear programming formulation

#     This is adapted from the method implemented in USBR's WTMP, described in 
#     the 2023 Model Development report
    
#     Parameters
#     ----------
#     resmod : TYPE
#         DESCRIPTION.
#     target_temp : TYPE
#         DESCRIPTION.
#     gateDict : TYPE
#         DESCRIPTION.
#     qoutPEN : TYPE
#         DESCRIPTION.
#     qleakTot : TYPE
#         DESCRIPTION.
#     nPtSinks : TYPE
#         DESCRIPTION.

#     Returns
#     -------
#     gateQ   : dict
#         A dictionary of flow amounts (values) for each open gate/port (keys)

#     """
    
#     # set the debug release level (0,1,2)
#     debug = resmod.Debug['Release']
    
#     qoutTCD = qoutPEN - qleakTot
#     gateQ = {k:[] for k,v in resmod.Outlets.items()}  # dictionary of outflows at each gate
#     totGatesOpen = sum(gateDict.values())
    
#     maxGateLev = max(gateDict)
#     minGateLev = min(gateDict)
    
#     ### NOTE: Shasta-specific code - adjusting side gate elevation vvvvv
#     # is this the first use of side gates in isolation?
#     if resmod.TimeStep==0:
#         prevGate = gateDict
#     else:
#         prevGate = resmod.PrevGateDict # sgateDict #resmod.Operations.GateOps[resmod.SimDates[resmod.TimeStep-1]]
        
#     if (gateDict[0]>0): # and (sum([prevGate[3],prevGate[2],prevGate[1]])<=0):
#         # this is the first time the side gates are used 
#         resmod.Operations.sidegate_only_cntr += 1
#         adj_sdg_ptsinks = True
        
#         # side gate top point sink elevation adjustment
#         sdgadjlevel_max = 832 
#         Ladj = sdgadjlevel_max - resmod.Outlets[0].TopElevFt
#         kadj = 0.4 # steepness
#         x0 = 45 #90 # midpoint (days) of logistic curve
        
#         v = Ladj/(1+np.exp(-1*kadj*(resmod.Operations.sidegate_only_cntr-x0)))
#         sdgadjlevel = sdgadjlevel_max-v
        
#         # sdgadjtime = 15
#         # sdgidx = min(resmod.Operations.sidegate_only_cntr, sdgadjtime)
#         # sdgadjlevel = min(resmod.WSE, (sdgadjlevel_max +
#         #                                (resmod.Outlets[0].TopElevFt-sdgadjlevel_max)
#         #                                *sdgidx/sdgadjtime))
#         if debug>1:
#             print(f"\nAdj sidegate level: {sdgadjlevel} ft")
#             print(f"Gate dict: 3: {gateDict[3]} | 2: {gateDict[2]} | 1: {gateDict[1]} | 0: {gateDict[0]}")
#     else:
#         resmod.Operations.sidegate_only_cntr = 0
#         adj_sdg_ptsinks = False
#     ### NOTE: Shasta-specific code - adjusting side gate elevation ^^^^^
        
#     ### NOTE: Shasta-specific code - adjusting middle gate elevation vvvvv
#     # testing lowering of middle gate point sinks after being open for long 
#     # periods (>60 days)
#     if resmod.TimeStepDate.timetuple().tm_yday==1 or resmod._time==0: # first day of the year, reset counter
#         resmod.Operations.midgatedays = 0
#     else:
#         if gateDict[2] >0:
#             resmod.Operations.midgatedays += 1
#         else:
#             resmod.Operations.midgatedays = 0
#     #gate_history = resmod.Operations.GateOps
#     if resmod.Operations.midgatedays > 60:
#         adj_mid_ptsinks = True
#         midadjlevel = resmod.Outlets[2].CtrElevFt
#     else:
#         adj_mid_ptsinks = False
#     ### NOTE: Shasta-specific code - adjusting middle gate elevation ^^^^^
    
#     # get heads at centerline of each gate level
#     #tw_elev = 586. # feet; from x-sect drawing of penstock 4, p78 in pdf of (United States Army Corps of Engineers, 1977)
#     glvlHdVels = {}
#     for g in gateDict:
        
#         # set head, velocity, and relative flow amount to zero (0) if gate isn't open
#         if gateDict[g]<=0:
#             glvlHdVels[g] = [ 0., 0., 0.]
#             continue
        
#         if g>minGateLev: # not at bottom gate level
#             g_below = g-1
#         else:
#             g_below = -1
        
#         gate_flow_onoff = 1
        
#         if resmod.CheckMinHead:
#             # deactivate gates if prescribed minimum head requirements aren't met
#             if g_below<0 or gateDict[g_below] ==0:  # no gates open on level below
#                 if resmod.WSE >= resmod.Outlets[g].MinHead + resmod.Outlets[g].BotElevFt:
#                     gate_flow_onoff = 1  # no gates open below, but meeting head requirement
#                 else:
#                     gate_flow_onoff = 0 # no gates open below and NOT meeting head reqmt - flow should be 0 thru gate
#                     # update gate dict for this time step accordingly
#                     gateDict[g] = 0
#             else:
#                 if gateDict[g_below]>0: # there's a gate open below
#                     if resmod.WSE >= resmod.Outlets[g].MinHead + resmod.Outlets[g].BotElevFt: 
#                         # there's water covering outlet gate and there's a gate open below
#                         gate_flow_onoff = 1  
#                     else:
#                         gate_flow_onoff = 0 
#                         # update gate dict for this time step accordingly
#                         gateDict[g] = 0
#                 else:
#                     gate_flow_onoff = 0  # there are no gates open below
#                     # update gate dict for this time step accordingly
#                     gateDict[g] = 0
                
            
#         # calculate the head above the outlet/get bottom elevation; see if it's
#         # more than the gate height/opening
#         B = min((resmod.WSE-resmod.Outlets[g].BotElevFt)*FTtoM, resmod.Outlets[g].GateHeight_m)
#         B = max(B, 0.0)
        
#         # calculate the head difference from the penstock elevation
#         dlel_m = max(0.0, resmod.WSE*FTtoM - max(resmod.Outlets[g].BotElevFt*FTtoM,resmod.PenstockElevation_m))

        
        
#         qrel = gate_flow_onoff*(resmod.Outlets[g].A1GT*
#                                 (dlel_m**(resmod.Outlets[g].B1GT))*
#                                 B**resmod.Outlets[g].G1GT)*gateDict[g]
    
#         glvlHdVels[g] = [dlel_m, None, qrel] # head diff; velocity (None, for now), relative q (will be changed later)
        
        
#     # set the constraints on the outlets - assume if the gate is open, it can have
#     # some minimum fraction (0.01) up to 1 of the flow
#     target_bounds = [] # constraints for each gate
#     gate_
#     for g in gateDict:
#         if gateDict[g] ==0:
#             target_bounds[g] == [0., 0.0]
#         else:
#             target_bounds[g] = [0.01,1]
    
#     # go through gates again, set the bounds & get the temps at each gate level
#     # assume gate_dict goes from top->down
#     for g in gate_dict:
        
        
#     return(gateQ)
            
def selective_withdrawal(resmod, qout, bypass_frac, rivDict,
                         nPtSinks=1, dz=1,
                         useMDWElkg=True, logging=False,**kwargs):
    '''
        an adaptation of the CE-QUAL-W2 selective withdrawal mechanism that 
        calculates a withdrawal envelope for purpopses of distribution of outflow to
        multiple structures; 
        
        Base gate configurations are set in resmod.GateDict
        Testing/optional gate configs are set by providing a value for gateDictOpt argument
    '''
    
    # set the debug release level (0,1,2)
    debug = resmod.Debug['Release']
    
    # pick the gate config to use
    if 'gateDictOpt' in kwargs:
        gateDict = copy.deepcopy(kwargs['gateDictOpt'])
    else:
        gateDict = resmod.GateDict
    

        
    totOutQ = 0.
    totOutE = 0.
    
    resmod.getWSE() # check WSE and top layer in case it needs updating
    
    topLyr = resmod.TopLyr
    topElev = resmod.WSE  # get the surface elevation calc'd at the beginning of this time step
    
    # Step 1. Allocate to outlet types
    [qoutPEN, qoutRIV, qoutSPILL] = alloc_to_outlet_types(resmod, qout, bypass_frac)

    # dz = vertical resolution, in meters
    if dz > 0:       
        # do calculation of withdrawal envelope on a finer vertical grid
        lelev_m = resmod.Layers[0].MinElev*FTtoM
        uelev_m = resmod.WSE*FTtoM
        grid = np.arange(lelev_m, uelev_m+dz/2, dz)
        #grid=np.append(grid, uelev_m)
        temp1 = [resmod.Layers[k].Temp for k in resmod.Layers]  #temperature on original grid
        rho1 = [resmod.Layers[k].Rho for k in resmod.Layers]  # density on original grid
        elev1 = [resmod.Layers[k].CtrElev*FTtoM for k in resmod.Layers]  # elevation of original grid
        elev1 = elev1 + [resmod.WSE]
        rho1 = rho1 + [rho1[-1]]
        temp1 = temp1 + [temp1[-1]]
        # interpolate original grid values to refined grid
        rho2 = np.interp(grid, elev1, rho1)
        if np.isnan(rho2[-1]):
            rho2[-1] = rho2[-2]
        temp2 = np.interp(grid, elev1, temp1)
        resmod.FineGrid = grid
        resmod.RhoFG = rho2
        resmod.TempFG = temp2
        resmod.LyrMapFG = np.round(np.interp(grid, elev1,list(range(resmod.nLyrs+1))))

    else:
        resmod.FineGrid = []
        
    # Step 2. Check if leakages exist, if so, then calculate
    if len(resmod.Leakages)>0 and qoutPEN>0:
        if useMDWElkg:
            [qleakDict, qleakTot] = calc_leakageMDWE(resmod, qoutPEN,
                                                     topElev, gateDict,
                                                     bypass_frac=bypass_frac)
        else:
            [qoutPEN, qoutFC, qleakDict, qleakTot] = calc_leakage(qoutPEN, topElev, gateDict)
    else:
        qleakDict = {}
        qleakTot = 0.0
    
    if resmod.Debug['Release']>0:
        print(f"\nQOUTPEN, QOUTRIV, QOUTFC, QLEAKTOT =  {qoutPEN:0.2f}, {qoutRIV:0.2f}, {qoutSPILL:0.2f}, {qleakTot:0.2f}") #%(qoutPEN,qoutFC, qleakTot))
    
    outqLyrDist = {}  # dictionary to hold all outflows as distributed to layers by gate or leakage zone
    
    outqLyrDist['Release_by_Outlet_Type'] = [qoutPEN, qoutRIV, qoutSPILL]
    
    # Step 3. now calculate the bulk flow through selective withdrawal levels
    # set whether to use LP optimization for tcd mixing, set target temp
    if resmod.LP_Opt_Blending:
        targ_col_name = resmod.Outflow.ColumnMap['tempTarg']
        # need to set tailwater target temp
        target_temp_tw = resmod.Outflow.DataFrame.loc[resmod.TimeStepDate,
                                                      targ_col_name]
        gateQ = tcd_alloc(resmod, gateDict, qoutPEN, qleakTot, nPtSinks, 
                          lp_opt=True, temp_targ = target_temp_tw, 
                          outqLyr=outqLyrDist, qleakDict=qleakDict,
                          rivDict=rivDict) 
    else:
        gateQ = tcd_alloc(resmod, gateDict, qoutPEN, qleakTot, nPtSinks, lp_opt=False) #removed topElev from args - use resmod.WSE  for up-to-date data
    
    #print("\nFinished TCD alloc - gateQ is:")
    #print(gateQ)
   # print('-----------------------------------')
    
    # find which selective withdrawal ports, g, are currently open
    openg = [g for g in gateQ if gateQ[g]!=[]]
    
    # outlet-layer assignments - layer with bottom elevation at or above 
    # outlet centerline elevation is assigned to that outlet level
    for og in openg:
        ptsnk_cntr=1
        for qstr, [outletEl,sqrtgh] in list(zip(gateQ[og][0],gateQ[og][1])):
            if debug >=1:
                print("Open Gate: %d - Outlet Elevation: %0.1f" %(og, outletEl))  #type(outletEl), 
        
            if qstr >0:
                # just loop over the layers between the bottom
                # and central assignements for this outlet - take
                # the highest one with water
                oLyr = resmod.Outlets[og].MinLayer
                for l in range(resmod.Outlets[og].MinLayer, resmod.Outlets[og].MaxLayer+1): #CtrLayer+1):
                    if (outletEl <= resmod.Layers[l].MaxElev) & \
                       (outletEl >= resmod.Layers[l].MinElev) & \
                       (resmod.Layers[l].Vol > 0.):
                        oLyr = l
                        
                if (len(gateQ[og][0])==1) & (oLyr == resmod.Outlets[og].MinLayer):
                    outletEl = resmod.Outlets[og].BotElevFt

                if debug>=2:
                    print("Calling wd_env for gates with parameters:   %d, %0.2f ft, %0.2f AF, %d, %0.2f ft" %(oLyr, outletEl, qstr, topLyr, topElev))
                
                if dz>0:
                    [outByLyr, totOutThisGate] = wd_env2(resmod, outletEl, qstr, topLyr, grid, rho2, debug=debug)
                    #print(totOutThisGate)
                else:
                    [outByLyr, totOutThisGate] = wd_env(resmod, oLyr, outletEl, qstr,topLyr, topElev, debug=debug) 
                
                if abs(totOutThisGate - qstr)>TOL_AF:
                    print("WARNING!! Calculated outflow at gate %s and point sink %s does not match what was prescribed!" %og)
                    print("WARNING!! Calculated outflow: %0.4f  --- Prescribed outflow %0.4f" %(totOutThisGate, qstr))
                outqLyrDist[str(og)+'g'+str(ptsnk_cntr)] = [outByLyr]  # adding 'g' to key to indicate gate
                ptsnk_cntr +=1
    
    # now do the outflow distributions for the leakage zones
    for ol in qleakDict:
        #lkgElev = resmod.Leakages[ol].CtrElev # @jmg 20230512 testing this difference in elev assignement
        lkgElev = resmod.Leakages[ol].BotElevFt
        qstr = qleakDict[ol]
        
        # just loop over the layers between the bottom
        # and central assignements for this outlet - take
        # the highest one with water
        lkgLyr = 9999
        for l in range(resmod.Leakages[ol].MinLayer, resmod.Leakages[ol].CtrLayer+1):
            if topElev > resmod.Layers[l].CtrElev: #self.Layers[l].Vol > 0.:
                lkgLyr = l
            else:
                lkgLyr = resmod.Leakages[ol].MinLayer
        if lkgLyr == 9999:
            print("Can't continue with selective withdrawal - no appropriate\nlayer assignment for leakage zone %d could be found." %ol)
            return(None)
        if lkgLyr == resmod.Leakages[ol].MinLayer:
            lkgElev = resmod.Leakages[ol].BotElevFt
            
            
        # now do the outflow allocation to layers based on the outlet level and flow
        if debug>=2:
            print("Calling wd_env with parameters:   %d, %0.2f ft, %0.2f AF, %d, %0.2f ft" %(lkgLyr, lkgElev, qstr, topLyr, topElev))
        
        if dz>0:
            [outByLyr, totOutThisLeak] = wd_env2(resmod, lkgElev, qstr, topLyr, grid, rho2, debug=debug)
        else:
            [outByLyr, totOutThisLeak] = wd_env(resmod, lkgLyr, lkgElev, qstr,topLyr, topElev, debug=debug) 
        
        if abs(totOutThisLeak-qstr)>TOL_AF: # != qstr:
            print("WARNING!! Calculated outflow at gate %s does not match what was prescribed!" %og)
            print("WARNING!! Calculated outflow: %0.2f  --- Prescribed outflow %0.2f" %(totOutThisLeak, qstr))
        outqLyrDist[str(ol)+'l'] = [outByLyr]  # adding 'l' to key to indicate leakage
        
    # do the outflow distributions for river outlets
    tmpFC = 0
    if qoutRIV >0.:
        tmpFC = qoutRIV
        # which river outlets are open?
        
        for ri,rv in resmod.RiverOutlets.items(): # assuming this is ordered from the top down
            if (tmpFC >0) and (resmod.WSE>rv.MinElev+1.) and rivDict[rv.ID]>0:
                rivoutqi = min(tmpFC, rv.Capacity_CFS*CFStoAFD)
                tmpFC = tmpFC - rivoutqi
                
                #print("Releasing %0.2f AF at river outlet %d" %(rivoutqi, ri))
                
                
                for l in range(rv.MinLayer, rv.CtrLayer+1):
                    if resmod.Layers[l].Vol > 0.:
                        rivoLyr = l

                if rivoLyr == rv.MinLayer:
                    rivoEl = rv.MinElev
                else:
                    #rivoEl = rv.CtrElev
                    #rivoEl = (0.5*rv.MinElev+0.5*rv.CtrElev)+rv.MinElev
                    rivoEl = rv.MinElev
                # rivoEl = rv.MinElev
                # rivoLyr = rv.MinLayer
                #print("River out layer: %d   River out elevation: %0.2f" %(rivoLyr, rivoEl))
                # now do the outflow allocation to layers based on river outlet level and flow
                if dz>0.:
                    [outByLyr, totOutThisRiv] = wd_env2(resmod, rivoEl, rivoutqi, topLyr,grid, rho2, debug=debug)
                else:
                    [outByLyr, totOutThisRiv] = wd_env(resmod, rivoLyr, rivoEl, rivoutqi, topLyr, topElev, debug=debug)
                outqLyrDist[str(ri)+'r'] = [outByLyr] # adding 'r' to indicate flow trhough river outlets
    #else:
   #     outqLyrDist['r'] = [None]
        
    # adding outflow distribution for flood release (over spillway/flood gates)
    if (qoutSPILL+tmpFC) >0:
        outqLyrDist['sp'] = [{topLyr: qoutSPILL+tmpFC}]
        # if tmpFC > 0: # if there is required release above river outlet capacity, assume it goes over spillway
        #     outqLyrDist['sp'] = [{topLyr: tmpFC}]


        
    if logging:
        with open(os.path.join(resmod.ProjDir, resmod.RunName +'_OutflowDebugLog.csv'), 'a') as lf:
            #ol = '\nTimeStep %s -- %s\n'  %(self.TimeStep, self.TimeStepDate)
            #lf.write(ol)
            #lf.write('============================================================\n')
            ol = 'TimeStep,DateTime,OutflowLoc,' + ','.join(['L%s' %li for li in range(resmod.nLyrs)])
            lf.write(ol + '\n')
            for ok, ov in outqLyrDist.items():
                
                lyrlist = []
                #print(ov)
                for ni in range(resmod.nLyrs):
                    if ni in ov[0]:
                        lyrlist.append(ov[0][ni])
                    else:
                        lyrlist.append(0.0)
                ol = str(resmod.TimeStep) + ',' + \
                    dt.datetime.strftime(resmod.TimeStepDate,'%Y%m%d') +\
                    ',' + ok + ',' + ','.join([str(nn) for nn in lyrlist]) + '\n'
                lf.write(ol)
                
    return(outqLyrDist)


def get_outflows_by_layer(resmod, gateQ, outqLyrDist,qleakDict,rivDict, dz=1):
    # outqLyrDist needs to be passed through from selective_withdrawal function 
    # to tcd_alloc or whatever function is calling get_outflows_by_layer
    # same with qleakDict and rivDict
    
    # set the debug release level (0,1,2)
    debug = resmod.Debug['Release']
    
    topLyr = resmod.TopLyr
    topElev = resmod.WSE  # get the surface elevation calc'd at the beginning of this time step
    
    qoutRIV = outqLyrDist['Release_by_Outlet_Type'][1]
    qoutPEN = outqLyrDist['Release_by_Outlet_Type'][0]
    qoutSPILL = outqLyrDist['Release_by_Outlet_Type'][2]
    
    # dz = vertical resolution, in meters
    if dz > 0:       
        # do calculation of withdrawal envelope on a finer vertical grid
        lelev_m = resmod.Layers[0].MinElev*FTtoM
        uelev_m = resmod.WSE*FTtoM
        grid = np.arange(lelev_m, uelev_m+dz/2, dz)
        #grid=np.append(grid, uelev_m)
        temp1 = [resmod.Layers[k].Temp for k in resmod.Layers]  #temperature on original grid
        rho1 = [resmod.Layers[k].Rho for k in resmod.Layers]  # density on original grid
        elev1 = [resmod.Layers[k].CtrElev*FTtoM for k in resmod.Layers]  # elevation of original grid
        elev1 = elev1 + [resmod.WSE]
        rho1 = rho1 + [rho1[-1]]
        temp1 = temp1 + [temp1[-1]]
        # interpolate original grid values to refined grid
        rho2 = np.interp(grid, elev1, rho1)
        if np.isnan(rho2[-1]):
            rho2[-1] = rho2[-2]
        temp2 = np.interp(grid, elev1, temp1)
        resmod.FineGrid = grid
        resmod.RhoFG = rho2
        resmod.TempFG = temp2
        resmod.LyrMapFG = np.round(np.interp(grid, elev1,list(range(resmod.nLyrs+1))))

    else:
        resmod.FineGrid = []
    
    
    # find which selective withdrawal ports, g, are currently open
    openg = [g for g in gateQ if gateQ[g]!=[]]
    
    # outlet-layer assignments - layer with bottom elevation at or above 
    # outlet centerline elevation is assigned to that outlet level
    for og in openg:
        ptsnk_cntr=1
        for qstr, [outletEl,sqrtgh] in list(zip(gateQ[og][0],gateQ[og][1])):
            if debug >=1:
                print("Open Gate: %d - Outlet Elevation: %0.1f" %(og, outletEl))  #type(outletEl), 
        
            if qstr >0:
                # just loop over the layers between the bottom
                # and central assignements for this outlet - take
                # the highest one with water
                oLyr = resmod.Outlets[og].MinLayer
                for l in range(resmod.Outlets[og].MinLayer, resmod.Outlets[og].MaxLayer+1): #CtrLayer+1):
                    if (outletEl <= resmod.Layers[l].MaxElev) & \
                       (outletEl >= resmod.Layers[l].MinElev) & \
                       (resmod.Layers[l].Vol > 0.):
                        oLyr = l
                        
                if (len(gateQ[og][0])==1) & (oLyr == resmod.Outlets[og].MinLayer):
                    outletEl = resmod.Outlets[og].BotElevFt

                if debug>=2:
                    print("Calling wd_env for gates with parameters:   %d, %0.2f ft, %0.2f AF, %d, %0.2f ft" %(oLyr, outletEl, qstr, topLyr, topElev))
                
                if dz>0:
                    [outByLyr, totOutThisGate] = wd_env2(resmod, outletEl, qstr, topLyr, grid, rho2, debug=debug)
                    #print(totOutThisGate)
                else:
                    [outByLyr, totOutThisGate] = wd_env(resmod, oLyr, outletEl, qstr,topLyr, topElev, debug=debug) 
                
                if abs(totOutThisGate - qstr)>TOL_AF:
                    print("WARNING!! Calculated outflow at gate %s and point sink %s does not match what was prescribed!" %og)
                    print("WARNING!! Calculated outflow: %0.4f  --- Prescribed outflow %0.4f" %(totOutThisGate, qstr))
                outqLyrDist[str(og)+'g'+str(ptsnk_cntr)] = [outByLyr]  # adding 'g' to key to indicate gate
                ptsnk_cntr +=1
    
    # now do the outflow distributions for the leakage zones
    for ol in qleakDict:
        #lkgElev = resmod.Leakages[ol].CtrElev # @jmg 20230512 testing this difference in elev assignement
        lkgElev = resmod.Leakages[ol].BotElevFt
        qstr = qleakDict[ol]
        
        # just loop over the layers between the bottom
        # and central assignements for this outlet - take
        # the highest one with water
        lkgLyr = 9999
        for l in range(resmod.Leakages[ol].MinLayer, resmod.Leakages[ol].CtrLayer+1):
            if topElev > resmod.Layers[l].CtrElev: #self.Layers[l].Vol > 0.:
                lkgLyr = l
            else:
                lkgLyr = resmod.Leakages[ol].MinLayer
        if lkgLyr == 9999:
            print("Can't continue with selective withdrawal - no appropriate\nlayer assignment for leakage zone %d could be found." %ol)
            return(None)
        if lkgLyr == resmod.Leakages[ol].MinLayer:
            lkgElev = resmod.Leakages[ol].BotElevFt
            
            
        # now do the outflow allocation to layers based on the outlet level and flow
        if debug>=2:
            print("Calling wd_env with parameters:   %d, %0.2f ft, %0.2f AF, %d, %0.2f ft" %(lkgLyr, lkgElev, qstr, topLyr, topElev))
        
        if dz>0:
            [outByLyr, totOutThisLeak] = wd_env2(resmod, lkgElev, qstr, topLyr, grid, rho2, debug=debug)
        else:
            [outByLyr, totOutThisLeak] = wd_env(resmod, lkgLyr, lkgElev, qstr,topLyr, topElev, debug=debug) 
        
        if abs(totOutThisLeak-qstr)>TOL_AF: # != qstr:
            print("WARNING!! Calculated outflow at gate %s does not match what was prescribed!" %og)
            print("WARNING!! Calculated outflow: %0.2f  --- Prescribed outflow %0.2f" %(totOutThisLeak, qstr))
        outqLyrDist[str(ol)+'l'] = [outByLyr]  # adding 'l' to key to indicate leakage
        
    # do the outflow distributions for river outlets
    tmpFC = 0
    if qoutRIV >0.:
        tmpFC = qoutRIV
        # which river outlets are open?
        
        for ri,rv in resmod.RiverOutlets.items(): # assuming this is ordered from the top down
            if (tmpFC >0) and (resmod.WSE>rv.MinElev+1.) and rivDict[rv.ID]>0:
                rivoutqi = min(tmpFC, rv.Capacity_CFS*CFStoAFD)
                tmpFC = tmpFC - rivoutqi
                
                #print("Releasing %0.2f AF at river outlet %d" %(rivoutqi, ri))
                
                
                for l in range(rv.MinLayer, rv.CtrLayer+1):
                    if resmod.Layers[l].Vol > 0.:
                        rivoLyr = l

                if rivoLyr == rv.MinLayer:
                    rivoEl = rv.MinElev
                else:
                    #rivoEl = rv.CtrElev
                    #rivoEl = (0.5*rv.MinElev+0.5*rv.CtrElev)+rv.MinElev
                    rivoEl = rv.MinElev
                # rivoEl = rv.MinElev
                # rivoLyr = rv.MinLayer
                #print("River out layer: %d   River out elevation: %0.2f" %(rivoLyr, rivoEl))
                # now do the outflow allocation to layers based on river outlet level and flow
                if dz>0.:
                    [outByLyr, totOutThisRiv] = wd_env2(resmod, rivoEl, rivoutqi, topLyr,grid, rho2, debug=debug)
                else:
                    [outByLyr, totOutThisRiv] = wd_env(resmod, rivoLyr, rivoEl, rivoutqi, topLyr, topElev, debug=debug)
                outqLyrDist[str(ri)+'r'] = [outByLyr] # adding 'r' to indicate flow trhough river outlets
    #else:
   #     outqLyrDist['r'] = [None]
        
    # adding outflow distribution for flood release (over spillway/flood gates)
    if (qoutSPILL+tmpFC) >0:
        outqLyrDist['sp'] = [{topLyr: qoutSPILL+tmpFC}]
        # if tmpFC > 0: # if there is required release above river outlet capacity, assume it goes over spillway
        #     outqLyrDist['sp'] = [{topLyr: tmpFC}]

    # retunr the updated dictionary of outflow by layers
    return(outqLyrDist)

def wd_env2(self, outletEl, qstr, topLyr, grid, rho2, debug=0):
    

    
    qstr_cms = qstr*AFDtoCMS
    outletEl_m = outletEl*FTtoM
    outLyr = np.abs(grid-outletEl_m).argmin()
    #print("outLyr: %s - density on this layer: %0.2f" %(outLyr, self.RhoFG[outLyr]))
    while np.isnan(self.RhoFG[outLyr]):
        outLyr -= 1
    topElev_m = self.WSE*FTtoM
    
    thisKTSW = np.abs(grid-self.KTSW_elev_m).argmin()
    thisKBSW = np.abs(grid-self.KBSW_elev_m).argmin()
    thisTopLyr = np.abs(grid - topElev_m).argmin()
    
    while np.isnan(self.RhoFG[thisTopLyr]): #move to next lower layer if top lyr selected is dry
        thisTopLyr -= 1
        
    KTOP = int(min(thisKTSW, thisTopLyr))
    KBOT = int(max(0, thisKBSW))
    
    ratio = (outletEl_m - grid[KBOT])/(topElev_m - grid[KBOT])
    coef = 1.
    if ((ratio<0.1) or (ratio>0.9)): coef=2.
    
    if debug>=2:
        print("DEBUG:: WITHDRAWAL ENVELOPE CALCULATION")
        print("########################################")
        print("DEBUG:: KTOP:  %d     KBOT:  %d" %(KTOP, KBOT))
        print("DEBUG:: Ratio: %0.3f, Coefficient %0.2f" %(ratio, coef))
    
     
    # find the avg density frequency (buoyancy frequency?) of region above structure
    if KTOP > outLyr:
        for k in range(outLyr+1, KTOP+1): #???? Exclude the outlet layer form calculation? # have to do KTOP+1 to make sure the layer indexed to KTOP gets inlcuded
            #print("layer: %s" %k)
            dtop = max(0.001, (grid[k] - outletEl_m))  # max() here to deal with instances were a layer minelev might be very close to the outlet elevation
            #w2usc92.for version --> #dft = np.sqrt((abs(self.Layers[k].Rho - self.Layers[outLyr].Rho)/dtop)*(9.81/self.Layers[outLyr].Rho))
            dft = np.sqrt( (abs(rho2[k] - rho2[outLyr]))/(dtop*rho2[outLyr] + 1e-10)*9.81)
            dft = max(dft, 1.0e-10)
            if np.isnan(dft):
                print("dft calculated as nan; component values are:")
                print("  Layer %s rho: %0.2f" %(k, rho2[k]))
                print("  Outlet layer %s rho: %0.2f" %(outLyr, rho2[outLyr]))
                print("  dtop:  %0.2f" %dtop)
                print("  denominator:  %0.2f\n" %(dtop*rho2[outLyr] + 1e-10))                    
            
            # half-height of upper withdrawal zone
            # assuming a point sink
            ztop = ((coef*qstr_cms)/dft)**0.333333
            if dtop >= ztop: 
                KTOP = k
                break
    else:
        dtop = 0.
        dft = 1.e-10
        ztop = 0.
        KTOP = outLyr
    
    eltop = outletEl_m + ztop
    if eltop < topElev_m:
        rdt = abs(rho2[outLyr] - rho2[KTOP])
    elif outletEl_m == topElev_m:
        rdt = 1.e-10
    else:
        rdt = abs(rho2[outLyr] - rho2[thisTopLyr])*ztop/(topElev_m-outletEl_m)
    
    rdt = max(rdt, 1e-10)
    
    if debug>=2:
        print("DEBUG:: dtop:  %0.4f   dft:  %0.4f   ztop: %0.4f   New KTOP: %d" %(dtop, dft, ztop, KTOP))
        print("DEBUG:: rdt: %0.4f" %(rdt))
    
    # find the avg density frequency (buoyancy frequency?) of region below structure
    for k in range(outLyr-1, KBOT-1, -1):
        dbot = max(0.000001,(outletEl_m - grid[k]))
        dfb = np.sqrt((abs(rho2[k] - rho2[outLyr])/dbot) * (9.81/rho2[outLyr]))
        dfb = max(dfb, 1.0e-10)
        
        # and the withdrawal zone half-heihgt
        # assuming point sink
        zbot = ((coef*qstr_cms)/dfb)**0.333333
        if dbot >= zbot:
            KBOT = k
            break
        
    # refrence density calculations
    elbot = outletEl_m - zbot
    if debug >=2:
        print("ELBOT = %0.2f m" %elbot)
        print("KBOT = %d\n\n" %KBOT)
    if elbot > (grid[thisKBSW+1]):
        rdb = abs(rho2[outLyr] - rho2[KBOT])
    else:
        rdb = abs(rho2[outLyr] - rho2[thisKBSW])*zbot/(outletEl_m - grid[thisKBSW+1])
    rdb = max(rdb, 1.e-10)
    
    
    # velocity profile
    vt = 0.0
    dendif = max(rdt, rdb)
    vnorm = {}
    for k in range(KTOP, KBOT-1, -1):
        #jo = k - KTOP  # should make this a zero-indexed set
        vnorm[k] = 1. - ((rho2[k] - rho2[outLyr])/dendif)**2
        vnorm[k] = max(vnorm[k], 0.0)
        vt += vnorm[k]
    
    #print(vnorm)
    # now finally apportion the flow to the layers
    totOutQ = 0.
    remainq = 0.
    lyrOutDict ={}
    for k in range(KTOP, KBOT-1, -1):
        
        lyrq = vnorm[k]/vt * qstr + remainq
        lyrOutDict[k] = lyrq
        totOutQ += lyrq
      
    return([lyrOutDict, totOutQ])

def wd_env(self, outLyr, outletEl, qstr, topLyr, topElev, debug=0):
    '''
        do the withdrawal envelope calculations for a given gate/leakage/outlet
        layer, elevation, and outflow
        return: outflow that should be removed from each layer
    '''
    
    qstr_cms = qstr*AFDtoCMS
    outletEl_m = outletEl*FTtoM
    topElev_m = topElev*FTtoM
    
    #qstr = outletQ # following CE-QUAL nomenclature of 'structure outflow'
    KTOP = int(min(self.KTSW, topLyr))
    KBOT = int(max(0, self.KBSW))
    
    ratio = (outletEl - self.Layers[KBOT].MinElev)/(topElev - self.Layers[KBOT].MinElev)
    coef = 1.
    if ((ratio<0.1) or (ratio>0.9)): coef = 2.
    
    if debug>=2:
        print("DEBUG:: WITHDRAWAL ENVELOPE CALCULATION")
        print("########################################")
        print("DEBUG:: KTOP:  %d     KBOT:  %d" %(KTOP, KBOT))
        print("DEGUB:: Ratio: %0.3f, Coefficient %0.2f" %(ratio, coef))
        
    # find the avg density frequency (buoyancy frequency?) of region above structure
    if KTOP > outLyr:
        for k in range(outLyr+1, KTOP+1): # Exclude the outlet layer form calculation # have to do KTOP+1 to make sure the layer indexed to KTOP gets inlcuded
            #print("layer: %s" %k)
            dtop = max(0.001, (self.Layers[k].MinElev - outletEl)*FTtoM)  # max() here to deal with instances were a layer minelev might be very close to the outlet elevation
            #w2usc92.for version --> #dft = np.sqrt((abs(self.Layers[k].Rho - self.Layers[outLyr].Rho)/dtop)*(9.81/self.Layers[outLyr].Rho))
            dft = np.sqrt( (abs(self.Layers[k].Rho - self.Layers[outLyr].Rho))/(dtop*self.Layers[outLyr].Rho + 1e-10)*9.81)
            dft = max(dft, 1.0e-10)
            if np.isnan(dft):
                print("dft calculated as nan; component values are:")
                print("  Layer %s rho: %0.2f" %(k, self.Layers[k].Rho))
                print("  Outlet layer %s rho: %0.2f" %(outLyr, self.Layers[k].Rho))
                print("  dtop:  %0.2f" %dtop)
                print("  denominator:  %0.2f" %(dtop*self.Layers[outLyr].Rho + 1e-10))                    
            
            # half-height of upper withdrawal zone
            # assuming a point sink
            ztop = ((coef*qstr_cms)/dft)**0.333333
            if dtop >= ztop: 
                KTOP = k
                break
    else:
        dtop = 0.
        dft = 1.e-10
        ztop = 0.
        KTOP = outLyr
    
    eltop = outletEl_m + ztop
    if eltop < topElev_m:
        rdt = abs(self.Layers[outLyr].Rho - self.Layers[KTOP].Rho)
    elif outletEl_m == topElev_m:
        rdt = 1.e-10
    else:
        rdt = abs(self.Layers[outLyr].Rho - self.Layers[topLyr].Rho)*ztop/(topElev_m-outletEl_m)
    
    rdt = max(rdt, 1e-10)
    
    if debug>=2:
        print("DEBUG:: dtop:  %0.4f   dft:  %0.4f   ztop: %0.4f   New KTOP: %d" %(dtop, dft, ztop, KTOP))
        print("DEGUB:: rdt: %0.4f" %(rdt))
    
    # find the avg density frequency (buoyancy frequency?) of region below structure
    for k in range(outLyr-1, KBOT-1, -1):
        dbot = max(0.000001,(outletEl - self.Layers[k].MinElev))*FTtoM
        dfb = np.sqrt((abs(self.Layers[k].Rho - self.Layers[outLyr].Rho)/dbot) * (9.81/self.Layers[outLyr].Rho))
        dfb = max(dfb, 1.0e-10)
        
        # and the withdrawal zone half-heihgt
        # assuming point sink
        zbot = ((coef*qstr_cms)/dfb)**0.333333
        if dbot >= zbot:
            KBOT = k
            break
        
    # refrence density calculations
    elbot = outletEl_m - zbot
    if debug >=2:
        print("ELBOT = %0.2f m" %elbot)
        print("KBOT = %d" %KBOT)
    if elbot > (self.Layers[self.KBSW+1].MinElev)*FTtoM:
        rdb = abs(self.Layers[outLyr].Rho - self.Layers[KBOT].Rho)
    else:
        rdb = abs(self.Layers[outLyr].Rho - self.Layers[self.KBSW].Rho)*zbot/(outletEl_m - (self.Layers[self.KBSW+1].MinElev)*FTtoM)
    rdb = max(rdb, 1.e-10)
    
    
    # velocity profile
    vt = 0.0
    #dendif = max(rdt, rdb)
    vnorm = {}
    for k in range(KTOP, KBOT-1, -1):
        if k > outLyr:  # ensure 
            densdiff = rdt
        else:
            densdiff = rdb
        vnorm[k] = 1. - ((self.Layers[k].Rho - self.Layers[outLyr].Rho)/densdiff)**2
        vnorm[k] = max(vnorm[k], 0.0)
        vt += vnorm[k]
    
    #print(vnorm)
    # now finally apportion the flow to the layers
    totOutQ = 0.
    #totOutE = 0.
    remainq = 0.
    lyrOutDict ={}
    for k in range(KTOP, KBOT-1, -1):
        
        lyrq = vnorm[k]/vt * qstr + remainq
        lyrOutDict[k] = lyrq
        #remainq += max(0., lyrq - ioutq)
        #self.Layers[k].Vol -= ioutq

        totOutQ += lyrq
        #totOutE += outTotE    
   
    return([lyrOutDict, totOutQ])

def gate_level_opts(self):
    # this function just gets a list of gate levels that
    # can be used based on pool elevation
    wse = self.WSE
    
    level_opts = []
    for gl in self.Outlets:
        if wse >= self.Outlets[gl].MinHead + self.Outlets[gl].BotElevFt:
            if gl < 99: # ignore 'dummy' gate zones if included (like for Shasta side gate leakage in Mike Deas' version)
                level_opts.append(gl)
    return(level_opts)  
    
def tcd_gate_open_opts(self, gate_level_opts):
    gate_options = []
    top_gate = max([o for o in self.Outlets if o < 99])
    # if len(gate_level_opts)>3:
    #     gate_level_opts = gate_level_opts[0:3]
    #     #gate_level_opts = gate_level_opts[1:4]
    
    # case where only the side gates are available - happens in CalSim scenarios with extreme drawdown
    if gate_level_opts[0] == 0:
        gatedict = {k:0 for k in self.Outlets if k <99}
        gatedict[0] = self.Outlets[0].NumGates
        gate_options.append(gatedict)
        
    for gl in gate_level_opts:
        if (gl < top_gate) & (gl+1 in gate_level_opts):
            for ug in range(self.Outlets[gl+1].NumGates, -1, -1):
                for lg in range(self.Outlets[gl].NumGates+1):
                
                    if (gl==0) & (lg == self.Outlets[gl].NumGates):
                        # both side gates are open, can have 0-5 of lower gates open
                        gatedict = {k:0 for k in self.Outlets if k <99}
                        gatedict[gl] = lg
                        gatedict[gl+1] = ug
                        gate_options.append(gatedict)
                    elif (gl==0) & (lg < self.Outlets[gl].NumGates):
                        if ug >3:
                            gatedict = {k:0 for k in self.Outlets if k <99}
                            gatedict[gl] = lg
                            gatedict[gl+1] = ug
                            gate_options.append(gatedict)
                    else:
                        if ug + lg >= 5:
                            gatedict = {k:0 for k in self.Outlets if k <99}
                            gatedict[gl] = lg
                            gatedict[gl+1] = ug
                            gate_options.append(gatedict)
        else: # move to next set of gat level options so that 'ug' above can reference upper gate level
            pass
        # else:  # case where only the side gates are available - happens in CalSim scenarios with extreme drawdown
        #     gatedict = {k:0 for k in self.Outlets if k <99}
        #     gatedict[0] = self.Outlets[0].NumGates
        #     gate_options.append(gatedict)
            
    return(gate_options)
#            else:
#                for lg in range(self.Outlets[gl].NumGates+1):

# def gate_options_increment(self, gate_level_opts, prev_gate_dict, outflow,
#                            ttarg, upp_tol, low_tol, nPtSinks=3,
                           
    
def gate_open_opts_incr(self, gate_level_opts, prev_gate_dict, outflow, 
                        targetTemp, 
                        tolerance, low_tolerance, nPtSinks=3, 
                        debugRelease=0, targetLoc='CCR', **kwargs):
    tsIdx = self.TimeStep
    #prevTemp = self.ReleaseTemps[tsIdx-1][0]
    
    kwk_to_ccr = [[1.96289502697935,0.961415095287775], #jan
                [2.01640832688523,0.968522054966461],  #feb
                [-6.60876513076067,1.15289820731661],  #mar
                [-2.70351765355408,1.07591746522629],  #apr
                [0.104476828309609,1.01883116933375],  #may
                [-2.72518633430037,1.07378504815437],  #jun
                [-1.21406813576229,1.04358686039431], #jul
                [-1.98045345947309,1.05989647744030],  #aug
                [2.60426157453926,0.971397038050233], # sep
                [4.26697107315805,0.931174409553145],  #oct
                [2.57265131067974,0.950725149278270],  #nov
                [0.958029934165394,0.975945198262967]] # dec
    
    shd_to_kwk = [[0.368330773672733,0.941041887648361],
                  [0.689539667116569,0.915332054936767],
                  [0.0591216338475616,1.02895155121519],
                  [0.316890327385639,1.03125597205824],
                  [0.503408156619173,1.01849618812037],
                  [-0.0291694668760519,1.08599143360539],
                  [-1.44398203686144,1.23322254104147],
                  [1.60237814080206,0.947430588900742],
                  [3.02003271507210,0.793783989913322],
                  [3.38486416855707,0.734750720922214],
                  [1.40827838574995,0.871442808436708],
                  [-0.307906698674068,0.996563159994448]]
                
    
    [outQ, outT, outE] = self.calcOutTemp(outflow, prev_gate_dict,
                                  nPtSinks=nPtSinks, 
                                  debugRelease=debugRelease)
    if targetLoc=='SHD':
        thisTemp1 = outT
    else:
        
        kwk = kwargs['kwk_model']
        kwkdat = kwargs['kwk_data']
        inTemps = [kwkdat['sppt'], outT]
        inQs = [kwkdat['sppq'], outQ]
        
        sumQ, wTemp = kwk.calcFlowWeightTemp(inQs, inTemps)
        
        


        #thiskwk=kwargs['kwk_model']
        if targetLoc =='KWK':
            #coeff = shd_to_kwk[self.TimeStepDate.month-1]
            #thisTemp1 = outT1 * coeff[1] + coeff[0]
            kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                              kwkdat['kairT'], \
                              stoCorrect=kwkdat['extraq'], \
                              calcTT=True, calcSolRad=True, 
                              radWm2=kwkdat['krad'], 
                              final=False)  
            thisTemp = kotemp
        elif targetLoc=='CCR':     
            #coeff_kwk = shd_to_kwk[self.TimeStepDate.month-1]
            coeff_ccr = kwk_to_ccr[self.TimeStepDate.month-1]
            #thisTemp1 = (outT1*coeff_kwk[1] + coeff_kwk[0])*coeff_ccr[1] + coeff_ccr[0]
            
            for fd in range(3):
                kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                              kwkdat['kairT'], \
                              stoCorrect=kwkdat['extraq'], \
                              calcTT=True, calcSolRad=True, 
                              radWm2=kwkdat['krad'], 
                              final=False)              
                thisTemp = kotemp * coeff_ccr[1] + coeff_ccr[0]
            
            if debugRelease>0:
                print(" -- Trying to meet CCR temperature target:")
                print('\t\t   Target: %0.1f deg F   Calc Temp: %0.1f\n' %(targetTemp, thisTemp))
        else:
            # default to Shasta for now
            thisTemp = outT
        
    #tdiff = thisTemp - targetTemp         
    #gate_options = []
    
    openlevels = [k for k,v in prev_gate_dict.items() if v>0]
    posslevels = gate_level_opts
    lowestopen = min(openlevels)
    highestopen = max(openlevels)
    lowestposs = min(posslevels)
    highestposs = max(posslevels)
    
    too_cold = thisTemp < (targetTemp - low_tolerance) # if evalulates as True, then temperature is too low - prefer to open an upper gate
    too_warm = thisTemp > (targetTemp + tolerance) # if evaluates as True, then temp is too warm, prefer to open a lower gate

    if too_warm:
    #if abs(tdiff)>tolerance and (tdiff > 0.):
        # need to open a lower gate
        #lowestopen = min([k for k,v in prev_gate_dict.items() if v>0])
        #topopen = False

        #print("\n======================\nIt's too warm and...")
        if lowestopen==highestopen: # there's only one level open, need to open below
            #print("\tThere's currently only one level open, so we should open a level below...")
            if lowestopen-1 in posslevels:
                #print("\t\tand as luck would have it, we can do that")
                #print("\t\t opening a gate on level %d" %(lowestopen-1))
                new_gate_dict = copy.deepcopy(prev_gate_dict)
                new_gate_dict[lowestopen-1] = new_gate_dict[lowestopen-1]+1
            else:
                #print("\t\tBut we can't open a lower level right now...")
                #print("\t\t opening all the gates on the lowest level we have access to, level %d" %lowestposs )
                new_gate_dict = copy.deepcopy(prev_gate_dict)
                new_gate_dict[lowestposs] = self.Outlets[lowestposs].NumGates
        else: # then there are already 2 levels open - try opening lower or closing upper
            #print("\tThere are two levels already open; we'll have to try closing an upper or opening a lower...")
            new_gate_dict = copy.deepcopy(prev_gate_dict)
            if prev_gate_dict[lowestopen] < self.Outlets[lowestopen].NumGates:
                #print("\t\t Opening a lower level gate: setting level %d to %d" %(lowestopen, prev_gate_dict[lowestopen]+1))
                new_gate_dict[lowestopen] = prev_gate_dict[lowestopen]+1
            else:
                new_gate_dict[highestopen] = max(0,prev_gate_dict[highestopen]-1)
    
    elif too_cold:
        #print("\n==============\nIt's too cold and..")
    #elif abs(tdiff)>tolerance and (tdiff < 0.):
        # need to open a higher gate
        #highestopen = max([k for k,v in prev_gate_dict.items() if v>0])
        if lowestopen==highestopen: # there's only one level open, open above
            #print("\tThere's currently only one level open, so we should open a level above...")
            if highestopen+1 in posslevels: # can we open a level above?
                #print("\t\t opening a gate on level above (%d)" %(highestopen+1))
                new_gate_dict = copy.deepcopy(prev_gate_dict)
                new_gate_dict[highestopen+1] = new_gate_dict[highestopen+1]+1
            else:
                #print("\t\tThere are no more levels above to open. Setting all gates on level %s to open" %highestposs)
                new_gate_dict = copy.deepcopy(prev_gate_dict)
                new_gate_dict[highestposs] = self.Outlets[highestposs].NumGates
        else: # there are 2 levels already open - try closing a lower or opening an upper
            #print("\tThere are two levels already open; we'll have to try closing a lower or opening an upper...")
            new_gate_dict = copy.deepcopy(prev_gate_dict)
            if prev_gate_dict[highestopen] < self.Outlets[highestopen].NumGates:
                new_gate_dict[highestopen] = prev_gate_dict[highestopen]+1
            else:
                new_gate_dict[lowestopen] = max(0, prev_gate_dict[lowestopen]-1)
                
#            if prev_gate_dict[highestopen] < self.Outlets[highestopen].NumGates:
#                new_gate_dict = copy.deepcopy(prev_gate_dict)
#                new_gate_dict[highestopen] = new_gate_dict[highestopen]+1
#                gate_options.append(new_gate_dict)
#            elif prev_gate_dict[highestopen-1] > 0:  # or close a lower gate
#                new_gate_dict= copy.deepcopy(prev_gate_dict)
#                new_gate_dict[highestopen-1] = new_gate_dict[highestopen-1]-1
#                gate_options.append(new_gate_dict)
#            else:
#                if highestopen+1 in gate_level_opts:
#                    new_gate_dict= copy.deepcopy(prev_gate_dict)
#                    new_gate_dict[highestopen-1] = 0
#                    new_gate_dict[highestopen+1] = new_gate_dict[highestopen+1]+1
#                else:
#                    new_gate_dict= copy.deepcopy(prev_gate_dict)
    else:
        new_gate_dict = prev_gate_dict
        
    # check the result using the new_gate_dict - if the temperature is too warm
    # and the amount by which it exceeds the tolerance is greater than the
    # amount by which it was too cold originally, then revert back to the original
    # configuration - trying to reduce oscillations
    [outQ1, outT1, outE1] = self.calcOutTemp(outflow, new_gate_dict,
                          nPtSinks=nPtSinks, 
                          debugRelease=debugRelease)
    
    if targetLoc=='SHD':
        thisTemp1 = outT1
    else:
        
        kwk = kwargs['kwk_model']
        kwkdat = kwargs['kwk_data']
        inTemps = [kwkdat['sppt'], outT1]
        inQs = [kwkdat['sppq'], outQ1]
        
        sumQ, wTemp = kwk.calcFlowWeightTemp(inQs, inTemps)
        


        #thiskwk=kwargs['kwk_model']
        if targetLoc =='KWK':
            #coeff = shd_to_kwk[self.TimeStepDate.month-1]
            #thisTemp1 = outT1 * coeff[1] + coeff[0]
            
            kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                              kwkdat['kairT'], \
                              stoCorrect=kwkdat['extraq'], \
                              calcTT=True, calcSolRad=True, 
                              radWm2=kwkdat['krad'], 
                              final=False)  
            
            thisTemp1 = kotemp
        elif targetLoc=='CCR':     
            #coeff_kwk = shd_to_kwk[self.TimeStepDate.month-1]
            coeff_ccr = kwk_to_ccr[self.TimeStepDate.month-1]
            #thisTemp1 = (outT1*coeff_kwk[1] + coeff_kwk[0])*coeff_ccr[1] + coeff_ccr[0]
            
            for fd in range(3):
                kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                                                  kwkdat['kairT'], \
                                                  stoCorrect=kwkdat['extraq'], \
                                                  calcTT=True, calcSolRad=True, 
                                                  radWm2=kwkdat['krad'], 
                                                  final=False)  
                thisTemp1 = kotemp * coeff_ccr[1] + coeff_ccr[0]
                
            if debugRelease>0:
                print(" -- Solution to meet CCR temperature target:")
                print('\t\t   Target: %0.1f deg F   Calc Temp: %0.1f\n' %(targetTemp, thisTemp1))
        else:
            # default to Shasta for now
            thisTemp1 = outT1
            
    too_warm_amt = max(0.,thisTemp1- targetTemp+tolerance)
    too_cold_amt = max(0., targetTemp-low_tolerance - thisTemp1)
    
    if (too_warm_amt > too_cold_amt) and (self.TimeStepDate.month < 1):
        print("------ reverting to previous gate config")
        new_gate_dict = prev_gate_dict
    
    return(new_gate_dict, thisTemp1)
        

            
def gate_select3(resmod, targetTemp, outflow, gate_options, tolerance=1., 
                 lower_temp_tolerance=1., nPtSinks=3, debugRelease=0,
                 targetLoc='CCR', **kwargs):
    
    ts_temp_tol = tolerance # tolerance (in deg F) for divergence from target temperature
                          # for changing gate config between time steps; if abs(diff between
                          # targetTemp and previous time steps release) is > `ts_temp_tol`,
                          # then will need to make a gate change
                          
    low_temp_tol = lower_temp_tolerance  # for having non-symmetric tolerance around the target
                                        # e.g. make it less likely to be below target early in season so that it can meet target later in season
                                        
    tsIdx = self.TimeStep
    
    if ['bypassFrac'] in kwargs:
        bypassFrac = kwargs['bypassFrac']
    else:
        bypassFrac = 0.
    
    # coefficients for linear model estimation
    kwk_to_ccr = [[1.96289502697935,0.961415095287775], #jan
                [2.01640832688523,0.968522054966461],  #feb
                [-6.60876513076067,1.15289820731661],  #mar
                [-2.70351765355408,1.07591746522629],  #apr
                [0.104476828309609,1.01883116933375],  #may
                [-2.72518633430037,1.07378504815437],  #jun
                [-1.21406813576229,1.04358686039431], #jul
                [-1.98045345947309,1.05989647744030],  #aug
                [2.60426157453926,0.971397038050233], # sep
                [4.26697107315805,0.931174409553145],  #oct
                [2.57265131067974,0.950725149278270],  #nov
                [0.958029934165394,0.975945198262967]] # dec
    
    shd_to_kwk = [[0.368330773672733,0.941041887648361],
                  [0.689539667116569,0.915332054936767],
                  [0.0591216338475616,1.02895155121519],
                  [0.316890327385639,1.03125597205824],
                  [0.503408156619173,1.01849618812037],
                  [-0.0291694668760519,1.08599143360539],
                  [-1.44398203686144,1.23322254104147],
                  [1.60237814080206,0.947430588900742],
                  [3.02003271507210,0.793783989913322],
                  [3.38486416855707,0.734750720922214],
                  [1.40827838574995,0.871442808436708],
                  [-0.307906698674068,0.996563159994448]]
    
#        gate_dict = self.Operations.GateOps[self.TimeStepDate-dt.timedelta(1)]

    if targetLoc=='SHD':
        prevTemp = np.mean([r[0] for r in resmod.ReleaseTemps[max(0, tsIdx-5):max(0, tsIdx-1)]])
        prevTempError = [r[0]-r[1] for r in resmod.ReleaseTemps[max(0, tsIdx-5):max(0, tsIdx-1)]]
    
    if targetLoc=='CCR':
        ccr_temp = kwargs['ccr_temp']
        #prevTemp = np.mean([r[0] for r in ccr_temp[max(0, tsIdx-3):max(0, tsIdx-1)]])
        #prevTempError = [r[0]-r[1] for r in ccr_temp[max(0, tsIdx-3):max(0, tsIdx-1)]]        
        prevTemp = ccr_temp[-1][0] #tsIdx-2][0]
        prevTempError = [r[0]-r[1] for r in ccr_temp[max(0, tsIdx-4):max(0, tsIdx-1)]] #ccr_temp[tsIdx-2][0] - ccr_temp[tsIdx-2][1] 
        
    increasing_errors = all(i < j for i, j in zip(prevTempError, prevTempError[1:]))
    tdiff = prevTemp - targetTemp

    best_error = 99
    best_gates = {}
         
    print("calculated tdiff = %0.2f" %tdiff)
    print("calculated tolerance level = %0.2f" %(1*ts_temp_tol))
    
    if targetTemp==99:
        gate_dict = gate_options[0]
        
    elif (resmod.TimeStep == 1) or (abs(tdiff) > 1*ts_temp_tol) or increasing_errors:
        if len(gate_options)<1:
            print("NO GATE OPTIONS!!")
            exit
        for gc in gate_options:
            if debugRelease>0:
                print("\n\n%s -- Trying gate option: " %resmod.TimeStepDate)
                print("      Upper:   %d" %gc[3])
                print("      Middle:  %d" %gc[2])
                print("      Lower:   %d" %gc[1])
                print("      Side:    %d" %gc[0])
            
            [outQ, outT, outE] = calcOutTemp(outflow, gc,
                                                  bypassFrac=bypassFrac,  #adding functionality to force bypass of TCD
                                                  nPtSinks=nPtSinks, 
                                                  debugRelease=debugRelease)
            
            if targetLoc=='SHD':
                thisTemp = outT
            else:
                kwk = kwargs['kwk_model']
                kwkdat = kwargs['kwk_data']
                inTemps = [kwkdat['sppt'], outT]
                inQs = [kwkdat['sppq'], outQ]
                
                sumQ, wTemp = kwk.calcFlowWeightTemp(inQs, inTemps)
                
 
                
                if targetLoc =='KWK':
                    #coeff = shd_to_kwk[self.TimeStepDate.month-1]
                    #thisTemp = outT * coeff[1] + coeff[0]
                    kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                                                      kwkdat['kairT'], \
                                                      stoCorrect=kwkdat['extraq'], \
                                                      calcTT=True, calcSolRad=True, 
                                                      radWm2=kwkdat['krad'], 
                                                      final=False) 
                    thisTemp = kotemp

                elif targetLoc=='CCR':     
                    #coeff_kwk = shd_to_kwk[self.TimeStepDate.month-1]
                    coeff_ccr = kwk_to_ccr[self.TimeStepDate.month-1]
                    #thisTemp = (outT*coeff_kwk[1] + coeff_kwk[0])*coeff_ccr[1] + coeff_ccr[0]
                    for fd in range(3):
                        kotemp, ksto = kwk.calcMixTempImp(sumQ, wTemp,kwkdat['kout'],
                                                          kwkdat['kairT'], \
                                                          stoCorrect=kwkdat['extraq'], \
                                                          calcTT=True, calcSolRad=True, 
                                                          radWm2=kwkdat['krad'], 
                                                          final=False) 
                    thisTemp = kotemp * coeff_ccr[1] + coeff_ccr[0]
                else:
                    # default to Shasta for now
                    thisTemp = outT
            
            #this_tdiff = outT - targetTemp
            this_tdiff = thisTemp - targetTemp
            if abs(this_tdiff) < abs(best_error):
                best_error = this_tdiff
                best_gates = gc
            
            if debugRelease >0:
                print("       Target: %0.2f - %0.2f - %0.2f deg F" %(targetTemp+ts_temp_tol, targetTemp, targetTemp-low_temp_tol))
                print("       ResultingTemp:  %0.2f deg F" %thisTemp)
            too_cold = thisTemp < (targetTemp - low_temp_tol) # if evalulates as True, then temperature is too low - prefer to open an upper gate
            too_warm = thisTemp > (targetTemp + ts_temp_tol) # if evaluates as True, then temp is too warm, prefer to open a lower gate
            
            gate_dict = gc
            #if this_tdiff < 0.:
            if (not too_warm) & (not too_cold):
                print("found an option!")
                break
        
        try:
            if gate_dict != best_gates:
                gate_dict = best_gates
        except:
            print("Gate search failed - repeating previous step's configuration")
            gate_dict = self.Operations.GateOps[self.TimeStepDate-dt.timedelta(1)]

    else:
        print("repeating last time step's gate ops")
        gate_dict = self.Operations.GateOps[self.TimeStepDate-dt.timedelta(1)]
        
    return(gate_dict)
        
    
def gate_select2(self, targetTemp, gatePriorities,outflow, tolerance=1., nPtSinks=3, debugRelease=0):
    
    ts_temp_tol = tolerance # tolerance (in deg F) for divergence from target temperature
                      # for changing gate config between time steps; if abs(diff between
                      # targetTemp and previous time steps release) is > `ts_temp_tol`,
                      # then will need to make a gate change
                      
    tsIdx = self.TimeStep
    prevTemp = np.mean([r[0] for r in self.ReleaseTemps[max(0, tsIdx-5):max(0, tsIdx-1)]])
    tdiff = prevTemp - targetTemp
    
    if (self.TimeStep == 1) or (tdiff > 4*ts_temp_tol):
        # if this is the first step, call the standalone `gate_select` function
        # to find a combination of gate openings that will meet the target temperature
        print("tdiff = %0.2f" %tdiff)
        newGateConfig, tdiff2, newOutT = self.gate_select(targetTemp, 
                                                         gatePriorities, 
                                                         outflow, nPtSinks=3, 
                                                         debugRelease=debugRelease)
        print("newGateConfig passed back from gate_select function:")
        print(newGateConfig)
    else:
        # make a decision based on what the previous temperature was, and if 
        # it's within a certain threshold of the target, don't make a change
        # if the previous temperature is beyond this tolerance, make a small
        # gate change (i.e. one or two adjacent gates open) in the direction
        # of the target (open upper gates to warm, lower gates to cool)

        
        # get the gate configuration from the previous step, as we'll need
        # to know what the options are for changing or for keeping the same
        prevDate = self.TimeStepDate - dt.timedelta(1)
        prevGateConfig = copy.deepcopy(self.Operations.GateOps[prevDate])
        newGateConfig = copy.deepcopy(self.Operations.GateOps[prevDate])
        newGateConfig1 = copy.deepcopy(self.Operations.GateOps[prevDate])
        newGateConfig2 = copy.deepcopy(self.Operations.GateOps[prevDate])
        newGateConfig3 = copy.deepcopy(self.Operations.GateOps[prevDate])
        newGateConfig4 = copy.deepcopy(self.Operations.GateOps[prevDate])
        topLevel = max([k for k,v in prevGateConfig.items() if v>0]) # get the top open gate
        botLevel = min([k for k,v in prevGateConfig.items() if v>0])
        if topLevel<max(prevGateConfig.keys()):
            nextUpLevel = topLevel+1
        else:
            nextUpLevel = topLevel
        if botLevel>min(prevGateConfig.keys()):
            nextDnLevel = botLevel-1
        else:
            nextDnLevel = botLevel
            
        nextUpLevelElev = self.Outlets[nextUpLevel].MinHead+self.Outlets[nextUpLevel].BotElevFt
        
        if abs(tdiff)> ts_temp_tol: # then make a gate change
            if tdiff > 0.: # water is too warm, open a gate below
                if botLevel==topLevel: # in the case where only one level is open, set the 
                    botLevel= nextDnLevel #max(0, topLevel-1)  # existing level as the top and make botLevel the level below
                if prevGateConfig[botLevel] < self.Outlets[botLevel].NumGates:
                    # there's a lower gate to be opened
                    newGateConfig[botLevel] = prevGateConfig[botLevel]+1
                elif prevGateConfig[topLevel] > 0: # there are gates above that can be closed
                    newGateConfig[topLevel] = prevGateConfig[topLevel] -1
                else:
                    # for the initial level pair, there are no more gates above to be closed, and no
                    # more gates below to be opened, so we have to move
                    # down a level, if possible
                    if botLevel>0:
                        newGateConfig[botLevel-1] = prevGateConfig[botLevel-1] + 1
                    else:
                        # we're already on the bottom level, set all gates to open
                        newGateConfig[botLevel] = self.Outlets[botLevel].NumGates
            elif (tdiff < 0.) and (self.WSE > nextUpLevelElev): # water is too cool or the reservoir pool has risen enough to 
                tempOpts = []
                if botLevel==topLevel: 
                    topLevel = nextUpLevel
                if prevGateConfig[topLevel] < self.Outlets[topLevel].NumGates:
                    # there's an upper gate to be opened
                    newGateConfig1[topLevel] = prevGateConfig[topLevel] + 1
                    [outQ, outT_upOpen, outE] = self.calcOutTemp(outflow, 
                                                                 newGateConfig1, 
                                                                 nPtSinks=nPtSinks,
                                                                 debugRelease=debugRelease)
                    tempOpts.append([outT_upOpen, outT_upOpen-targetTemp, newGateConfig1])
                    
                if prevGateConfig[botLevel] > 0: # there are bottom gates that can be closed
                    
                    newGateConfig2[botLevel] = prevGateConfig[botLevel] -1
                    if (botLevel==0) & (newGateConfig2[botLevel]<2):
                        newGateConfig2[topLevel] = self.Outlets[topLevel].NumGates
                    elif sum(newGateConfig2.values()) < 5:
                        newGateConfig2[topLevel] = self.Outlets[topLevel].NumGates
                    else:
                        pass
                    [outQ, outT_dnClose, outE] = self.calcOutTemp(outflow, 
                                                                 newGateConfig2, 
                                                                 nPtSinks=nPtSinks,
                                                                 debugRelease=debugRelease)
                    tempOpts.append([outT_dnClose,outT_dnClose-targetTemp, newGateConfig2])
                    
                # try both open up and close down
                if (prevGateConfig[botLevel]>0) & (prevGateConfig[topLevel] < self.Outlets[topLevel].NumGates):
                    newGateConfig3[botLevel] = prevGateConfig[botLevel]-1
                    newGateConfig3[topLevel] = prevGateConfig[topLevel]+1
                    [outQ, outT_updn, outE] = self.calcOutTemp(outflow, 
                                                                 newGateConfig3, 
                                                                 nPtSinks=nPtSinks,
                                                                 debugRelease=debugRelease)
                    tempOpts.append([outT_updn, outT_updn-targetTemp, newGateConfig3])
                    
                if (prevGateConfig[topLevel] == self.Outlets[topLevel].NumGates) & \
                    (prevGateConfig[botLevel]==0):
                    if topLevel < max(prevGateConfig.keys()):
                        newGateConfig4[topLevel+1] = prevGateConfig[topLevel+1]+1
                    else:
                        # we're on the top level, so open up all gates
                        newGateConfig4[topLevel] = self.Outlets[topLevel].NumGates
                    [outQ, outT_upNextLev, outE] = self.calcOutTemp(outflow, 
                                                                 newGateConfig4, 
                                                                 nPtSinks=nPtSinks,
                                                                 debugRelease=debugRelease)     
                    tempOpts.append([outT_upNextLev, outT_upNextLev-targTemp, newGateConfig4])
                
                bestopt = copy.deepcopy(prevGateConfig)
                best_prev_tdiff = tdiff
                #new_best_abstdiff = min([abs(c[1]) for c in tmpOpts])
                #new_best_tdiff = min([])
                for t,d, c in tempOpts:
                    if abs(d) < abs(best_prev_tdiff):
                        if d < 0:
                            best_prev_tdiff = d
                            bestopt = c
                        elif abs(d) < ts_temp_tol:
                            best_prev_tdiff = d
                            bestopt = c
                        else:
                            pass
                newGateConfig = bestopt
            else:
                newGateConfig = prevGateConfig
        else:
            newGateConfig = prevGateConfig
        
    return(newGateConfig)
        
        
def gate_select(self, targetTemp, gatePriorities, outflow, nPtSinks=3, debugRelease=0 ):
    '''
        targetTemp: downstream (tailwater) temperature target to try to hit
                    using blending from one or more gate levels
        gatePriorities: dictionary of priorities given to gate levels, negative 
                    number indicates flow must be allocated to this outlet, 
                    numbers 0 or greater indicate level of decreasing priority
                    (0 = highest priority, flow through this first; 1 = next priority, etc)
    '''
    
    # TODO: add check to see which gates are under water and meeting minimum 
    # head criteria (e.g. 35 ft for Shasta)
    
    # find the top layer and top elevation
    
    gateDict = {}
    for g in gatePriorities.keys():
        gateDict[g] = 0
    
#        topLyr = 0
#        for lj in range(0, self.nLyrs):
#            if self.Layers[lj].Vol >0:
#                topLyr = lj
#        topElev = self.Layers[topLyr].volToElevInterp()        
    #self.getWSE()  #<-- should be called as part of run sequence
    topElev = self.WSE
    
    newGatePriorities = {}
    newGatePriorities[0] = []
    for gpi, gpv in gatePriorities.items():
        if topElev > self.Outlets[gpi].MinHead+self.Outlets[gpi].BotElevFt:
            oLyr = self.Outlets[gpi].MinLayer
            newGatePriorities[gpv*1] = [gpi,self.Layers[oLyr].MinElev,self.Layers[oLyr+1].Temp]
        else:
            newGatePriorities[gpv*0].append([gpi,-9999, -9999.])
    print(newGatePriorities)
    
    botPort = topPort = -9999
    lowestPort = min([g for g in gatePriorities.keys() if gatePriorities[g]>0]) #min(gatePriorities.keys())
    highestPort = max([g for g in gatePriorities.keys() if gatePriorities[g]>0]) # max(gatePriorities.keys())
    lowestPortElev = self.Outlets[lowestPort].BotElevFt
    highestPortElev = self.Outlets[highestPort].BotElevFt
    highestUsablePort = max([newGatePriorities[p][0] for p in newGatePriorities if p>0])
    highestUsablePortElev = self.Outlets[highestUsablePort].BotElevFt
    topPortTemp = self.Layers[self.Outlets[highestPort].MinLayer].Temp #newGatePriorities[highestPort][2]
    botPortTemp = self.Layers[self.Outlets[lowestPort].MinLayer].Temp #newGatePriorities[lowestPort][2]
    
    
    for pri in [p for p in sorted(newGatePriorities.keys(), reverse=False) if p >0]:
        #print(pri)
        if newGatePriorities[pri][2] >= targetTemp:  # if temp at gate is > target
            topPort = newGatePriorities[pri][0]
            topPortTemp = newGatePriorities[pri][2]
            
        if newGatePriorities[pri][2] < targetTemp:  # if temp at gate is less than target
            botPort = newGatePriorities[pri][0]
            botPortTemp = newGatePriorities[pri][2]
            
        if newGatePriorities[pri][1] <= lowestPortElev:  # set the lowest port to equal the lowest (elevation-wise) prioritized gate
            lowestPort = newGatePriorities[pri][0]
            
        if newGatePriorities[pri][1] <= highestUsablePortElev:
            highestPort = newGatePriorities[pri][0]
        
    topPort = highestUsablePort
    print("top port: %d" %topPort)
    print("bot port: %d" %botPort)
    # if botPort still equals -9999, that measn there is no prioritized
    # port meeting the head requirements with temp colder than the target
    # i guess in this instance we go with all flow to the lowest elevation 
    # prioritized outlet?
    if botPort == -9999:
        gateDict[lowestPort] = self.Outlets[lowestPort].NumGates
        [outQ, outT, outE] = self.calcOutTemp(outflow, gateDict, nPtSinks, debugRelease)
        #return([gateDict, outQ, outT, outE])
        print("returning without optimizing - all bottom gates should be open")
        return((gateDict, targetTemp - outT, outT))
    else:
        if (topPort >=0) and (topPort - botPort > 1): # force the port pairs to be adjacent
            topPort = botPort +1
        if topPort==botPort:
            botPort = max(0, topPort-1)
            

    # if topPort still equals -9999, that means there is no prioritized
    # port meeting the head requirements with a temperature warmer than the
    # target, ergo there's no water that can be mixed to warm a release 
    # temperature up to the target - just release from the highest ports
#        if topPort == -9999:
#            gateDict[highestUsablePort] = self.Outlets[highestUsablePort].NumGates
#            [outQ, outT, outE] = self.calcOutTemp(outflow, gateDict, nPtSinks, debugRelease)
#            #return([gateDict, outQ, outT, outE])
#            print("returning without optimizing - all upper gates should be open")
#            return((gateDict, targetTemp - outT, outT))
    
    topFrac = abs((targetTemp - botPortTemp)/(topPortTemp-botPortTemp+TOL_TEMP))
    #botFrac = 1- topFrac
    
    print("For target temp: %0.2f deg F, trying gate levels %s and %s" %(targetTemp, topPort, botPort))
    gateDict[topPort] = self.Outlets[topPort].NumGates
    gateDict[botPort] = np.floor(min(self.Outlets[botPort].NumGates, 
                        max(self.Outlets[botPort].NumGates, 
                            (self.Outlets[topPort].BalanceFactor/topFrac -1)*gateDict[topPort])))
    
    options = {}
    p=0
    minGatesOpen = min(self.Outlets[botPort].NumGates, self.Outlets[topPort].NumGates)
    for u in range(0, self.Outlets[topPort].NumGates+1):
        for l in range(self.Outlets[botPort].NumGates,-1,-1):
            if (u==0) and (l ==0):
                pass
            elif u+l < minGatesOpen:
                pass
            else:
                options[u,l] = p #.append([u, l])
                p += 1
    
    opt2 = {v:k for k,v in options.items()}
    print(options)
    t1 = options.pop((int(gateDict[topPort]), int(gateDict[botPort])))  # this is the 'temperature' ranking of shutter configurations - higher # = expected warmere temps
    t2 = opt2.pop(t1)  # this is the shutter config corresponding to the temperature ranking given by t1
    twarm = options[opt2[p-1]]
    tcold = options[opt2[0]]
    
    [outQ, outT, outE] = self.calcOutTemp(outflow, gateDict, nPtSinks=nPtSinks, debugRelease=debugRelease)
    
    tdiff = targetTemp - outT
    prevdiff = tdiff
#        trials = []
#        trials.append([gateDict, outT])
    
    bestsofar = (gateDict, tdiff, outT)
    allelsefails = (gateDict, tdiff, outT)
    possible_configs = []
    #possible_configs.append(bestsofar)
    for i in range(len(options)):
        if abs(tdiff) < 0.05:
            possible_configs.append(bestsofar)
            print("found solution in %s iterations" %(i+1))
            break
        
        prev_opt_diff = twarm - t1 #options([gateDict[topPort], gateDict[botPort]])
        prev_opt_diff_cold = t1 - tcold
        
        if (tdiff > 0) and (prevdiff < 0.): # current proposed release temp is too cold, previous was too warm
            # so see if there's an option in between
            if prev_opt_diff > 0:  # there's a gate combination between what was just tried and the warmest end member
                next_option_idx2 = t1 + max(1, int(np.ceil(prev_opt_diff/2)))
            else:
                # no more options to choose from - go with the best so far
                print("no more options to choose from - go with the best so far")
                possible_configs.append(bestsofar)
                break
        elif (tdiff > 0.) and (prevdiff > 0.):  # current proposed release is too cold, and the previous one was too
            # go to the warmest option
            next_option_idx2 = max(opt2)
            twarm = next_option_idx2
            
        elif (tdiff < 0.) and (prevdiff > 0.): # current release is too warm, previous was too cold
            # see if there's an option in between
            if prev_opt_diff_cold > 0:
                next_option_idx2 = t1 - max(1, int(np.ceil(prev_opt_diff_cold/2)))
            else:
                # no more options to choose from - go with the best so far
                print("no more options to choose from - go with the best so far")
                break
        elif (tdiff < 0.) and (prevdiff < 0.):
            # go to coldest option
            next_option_idx2 = min(opt2)
            tcold = next_option_idx2
        try:
            testU, testL = opt2.pop(next_option_idx2)
            t2 = options.pop((testU, testL))
        except:
            # if for some reason the 'next_option_idx2' doesn't work, just pop the next item from the working list
            t2,(testU, testL) = opt2.popitem()
            #t2 = options.pop((testU, testL))
        
        gateDict[topPort] = testU
        gateDict[botPort] = testL
        [newOutQ, newOutT, newOutE] = self.calcOutTemp(outflow, gateDict, nPtSinks=nPtSinks, debugRelease=debugRelease)
        
        if debugRelease>0:
            print("\n=============================")
            print("Tested gate config: ")
            print(gateDict)
            print("estimated outflow temperature is: %0.2f" %newOutT)
            print("=============================\n")
        prevdiff = tdiff
        tdiff = targetTemp - newOutT
        if (abs(tdiff) < abs(bestsofar[1])): # & (tdiff>0.) :
            
            if debugRelease>0:
                print("--Updating best option to:")
                print(gateDict)
                print("******************************\n")
            bestsofar = (copy.deepcopy(gateDict), tdiff, newOutT)
            possible_configs.append(copy.deepcopy(bestsofar))
            
    min_abs_tdiff = min([abs(td[1]) for td in possible_configs])
    
    print(possible_configs)
    for pc in possible_configs:
        if (pc[1] > 0) & (self.TimeStepDate.month in [4,5,6,7,8,9]):
            bestsofar = pc
        elif pc[1]==min_abs_tdiff:
            bestsofar = pc
        else:
            bestsofar= possible_configs[0]
#        for i in range(50):
#            if abs(tdiff) < 0.05:
#                print("found solution in %s iterations" %(i+1))
#                break
##            if i < 2:
##                td
#            if tdiff >0: # blended temp is too cold, but more upper gates can be opened or lower gates closed
#                if gateDict[topPort] == self.Outlets[topPort].NumGates: # all upper gates are opened, will need to try closing lower
#                    gateDict[botPort] = max(0, gateDict[botPort]-1)
#                else:   # can open an upper gate
#                    gateDict[topPort] = min(self.Outlets[topPort].NumGates, gateDict[topPort]+1)
#            if tdiff < 0: # blended temp is too warm, but lower gates can be opened or upper gates closed
#                if gateDict[botPort] == self.Outlets[botPort].NumGates:  # all bottom gates are open, will need to close upper
#                    gateDict[topPort] = max(0, gateDict[topPort]-1)
#                else: # can open a lower gate
#                    gateDict[botPort] = min(self.Outlets[botPort].NumGates, gateDict[botPort]+1)
#            [newOutQ, newOutT, newOutE] = self.calcOutTemp(outflow, gateDict, nPtSinks, debugRelease)
#            tdiff = targetTemp - newOutT
#            outT = newOutT
#            trials.append([gateDict, outT])
#            print("\nIteration %s:\nTried configuration Level %s: %s gates, Level %s: %s gates\n ----Resulting temp: %0.2f" %(i, topPort, gateDict[topPort], botPort, gateDict[botPort], outT))
#        print(outQ, outT)
    
    print("this is the best so far:")
    print(bestsofar)
    if type(bestsofar[0])==type(None):
        print("didn't find a better solution..going with fallback:")
        print(allelsefails)        
        bestsofar=allelsefails
    
    return(bestsofar) #[gateDict, outQ, outT, outE])

def gateIncrement(self, targetTemp, outflow, nPtSinks=3, debugRelease=0):
    '''
        instead of doing a full search for a gate config, try out options
        that include only an increment from current configuration - 
        open a single lower gate/close a single upper gate, etc
    '''
    prevDate = self.TimeStepDate - dt.timedelta(1)
    thisGateConfig = copy.deepcopy(self.Operations.GateOps[prevDate])
    
    options = []
    for g in thisGateConfig:
        if thisGateConfig[g] >0:  # if some gates are open on one level
            if g-1 in thisGateConfig:  # and we're not at the bottom of the stack of gates
                if thisGateConfig[g-1] < self.Outlets[g].NumGates: # and the gate below as at least one more gate to open
                    tmp = copy.deepcopy(thisGateConfig)
                    starting_point = copy.deepcopy(tmp[g])
                    for tg in range(starting_point, -1, -1): # incrementally close upper gates
                        if g-1==0:
                            minopengates = 2
                        else:
                            minopengates=5
                        tmp[g] = tg
                        # making sure there are at least 5 total gates open, 
                        # (2 if using side gates), incremental open lower gates
                        for bg in range(max(0, minopengates-tg), self.Outlets[g-1].NumGates+1): 
                            print(tg, bg)
                            tmp[g-1] = bg
                            [newOutQ, newOutT, newOutE] = self.calcOutTemp(outflow, copy.deepcopy(tmp), nPtSinks, debugRelease)
                            #print(tmp, newOutT)
                            if newOutT <= targetTemp:
                                options.append([copy.deepcopy(tmp), newOutT])
                                return(options)
            else:
                print("side gates open already...can't go much colder")
                tmp = copy.deepcopy(thisGateConfig)
                if thisGateConfig[g] < self.Outlets[g].NumGates:
                    tmp[g] = self.Outlets[g].NumGates
                    [newOutQ, newOutT, newOutE] = self.calcOutTemp(outflow, copy.copy(tmp), nPtSinks, debugRelease)
                    if newOutT <= targetTemp:
                        options.append([copy.deepcopy(tmp), newOutT])
                        return(options)
                else:
                    [newOutQ, newOutT, newOutE] = self.calcOutTemp(outflow, copy.copy(tmp), nPtSinks, debugRelease)
                    if newOutT <= targetTemp:
                        options.append([copy.deepcopy(tmp), newOutT])
                        return(options)
                    
    options.append([copy.deepcopy(tmp),-9999])  # if nothing else worked so far, have to go with last option then
    return(options)
                    
            

def calcOutTemp(self, outflow, rivDict, nPtSinks=3, bypassFrac=0, **kwargs):
    '''
    kwarg: gateDictOpt - use this to override the resmod.GateDict (used for 
                                                                   testing options)
    '''
    
    # pick the gate config to use
    if 'gateDictOpt' in kwargs:
        opt_gateDict = copy.deepcopy(kwargs['gateDictOpt'])
        # releases
        outQLyrDist = selective_withdrawal(self, outflow, bypassFrac, rivDict,
                                           nPtSinks=nPtSinks, logging=False,
                                           dz=0, gateDictOpt = opt_gateDict)
    else:
        base_gateDict = self.GateDict
        # releases
        outQLyrDist = selective_withdrawal(self, outflow, bypassFrac,
                                           rivDict,
                                           nPtSinks=nPtSinks, logging=False,
                                           dz=0)
    

    check_outflows = outQLyrDist.pop('Release_by_Outlet_Type')
    #print(check_outflows)
    # create a dictionary to sum outflows assigned to each layer
    totQ = 0.
    tmplq = {}
    for l in self.Layers:
        tmplq[l]=0.
    for ri, rv in outQLyrDist.items():
        #print(rv)
        for rli, rlv in rv[0].items(): # loop through the layer assignments for this outlet flow
            if np.isnan(rlv):
                print("Calc'd outflow for %s is %0.2f" %(ri, rlv))
                return
            tmplq[rli] += rlv
            totQ += rlv
    
    if abs(totQ - outflow)>0.01:
        print("Calculated outflow (%0.4f af) does not match intended outflow (%0.4f af)!" %(totQ, outflow))
        
    [totQ2, outT, outE] = outflow_dist_forGateSelect(self, tmplq) 
    
    return([totQ2, outT, outE])
    
def outflow_dist(self, qoutLyr):
    '''
     this function takes a dictionary (key: layer index, value: outflow in AF)
     and does the accounting to remove water from the reservoir and calculate
     the effective outflow temperature
    '''
    totOutE = 0.
    totOutQ = 0.
    rmndr = 0.
    for l, qo in qoutLyr.items():
        if qo > 0.:
            iqo = min(self.Layers[l].Vol, qo+rmndr)
            rmndr = qo - iqo
            self.Layers[l].Vol -= iqo  # remove volume from layer
            tmpOutE = iqo* self.Layers[l].Temp  # calc the vol*temp for outflow from layer
            totOutE += tmpOutE  # accumulate vol*temp in total outflow
            totOutQ += iqo    # accumulate volume in total outflow
            self.Layers[l].TotE = self.Layers[l].Vol*self.Layers[l].Temp
         
    outTemp = totOutE/totOutQ
    return([totOutQ, outTemp, totOutE])
    

def outflow_check_and_dist(self,outQLyrDist, assignedOutQ ):    

    #print(outQLyrDist)
    # create a dictionary to sum outflows assigned to each layer
    totQ = 0.
    tmplq = {}
    for l in range(max(self.nLyrs,len(self.FineGrid))): #self.Layers:
        tmplq[l]=0.
        
    tcdg = {f"{tg}g":0 for tg in self.Outlets}
    rivg = {f"{r}r":0 for r in self.RiverOutlets}
    lkgg = {f"{lk}l":0 for lk in self.Leakages}
    spg = {"sp":0}
    q_by_outlet = {**tcdg,**rivg,**lkgg,**spg} # dict to organize sum of outflows by outlet level
    
    for ri, rv in outQLyrDist.items():

        for rli, rlv in rv[0].items(): # loop through the layer assignments for this outlet flow
            if np.isnan(rlv):
                print("Calc'd outflow for %s is %0.2f" %(ri, rlv))
                return
            tmplq[rli] += rlv
            q_by_outlet[ri[0:2]] += rlv
            totQ += rlv
    
    # q_by_outlet2 = {}
    # for ps, pv in q_by_outlet.items():
    #     if 'g' in ps:
    #         k = ps[0:2]
    #     else:
    #         k = ps
    #     q_by_outlet2[k] += pv
    
    if abs(totQ - assignedOutQ)>0.01:
        print("Calculated outflow (%0.4f af) does not match intended outflow (%0.4f af)!" %(totQ, assignedOutQ))
        return        

    # now do the accounting by reservoir computational layer, mapping fine grid
    # back to course grid if needed
    totOutE = 0.
    totOutQ = 0.
    rmndr = 0.

    # check if we need to map a refined grid back to the computational grid
    if len(tmplq)> self.nLyrs:
        
        # get aggregated volume by computation layers
        tmplqCG = {}
        tmpleCG = {}
        for cl in range(self.nLyrs):
            if self.Layers[cl].Vol > 0.:
                thisqt = [(tmplq[ki], self.TempFG[ki]) for ki in np.where(self.LyrMapFG==cl)[0] if ~np.isnan(self.TempFG[ki])]
                tmplqCG[cl] = sum([t[0] for t in thisqt])
                tmpleCG[cl] = sum([t[0]*t[1] for t in thisqt])
        #print(tmplqCG)
        #print(tmpleCG)
            
        #for l, qo in tmplqCG.items():
        for l in sorted(list(tmplqCG.keys()), reverse=True):
            qo = tmplqCG[l]
            if qo>0.:
                iqo = min(self.Layers[l].Vol, qo+rmndr)
                rmndr = qo-iqo
                self.Layers[l].Vol -= iqo # remove volume from this layer
                tmpOutE = (iqo/qo)*tmpleCG[l]   # calc the vol*temp for outflow from the layer, adjusting for any difference in withdrawal
                if np.isnan(tmpOutE):
                    print("Layer %s, outflow: %0.2f, q*E: %0.2f" %(l, iqo, tmpOutE))
                totOutE += tmpOutE  # accumulate vol*temp in total outlfow
                totOutQ += iqo   # accumulate voluem in total outflow
                self.Layers[l].TotE -= tmpOutE # self.Layers[l].Vol*self.Layers[l].Temp # update vol*temp in layer
                self.Layers[l].Temp = self.Layers[l].TotE/self.Layers[l].Vol
                #print("Layer %s - New TotE: %0.2f - New Vol: %0.2f - New Temp: %0.2f" %(l, self.Layers[l].TotE,self.Layers[l].Vol, self.Layers[l].Temp))
    else:
        this_botlayer = 0 # lowermost layer here that has water being removed from it
        for l, qo in tmplq.items():
            if qo > 0.:
                iqo = min(self.Layers[l].Vol, qo+rmndr)
                rmndr = qo - iqo
                self.Layers[l].Vol -= iqo  # remove volume from layer
                tmpOutE = iqo* self.Layers[l].Temp  # calc the vol*temp for outflow from layer
                totOutE += tmpOutE  # accumulate vol*temp in total outflow
                totOutQ += iqo    # accumulate volume in total outflow
                self.Layers[l].TotE = self.Layers[l].Vol*self.Layers[l].Temp
                this_botlayer = l

        if rmndr >0:
            #print(f"\n\t\tAAACCCCCKKKKKKK!!!! Missing some water that should be removed: {rmndr}")
            if self.Layers[this_botlayer-1].Vol > rmndr:
                self.Layers[this_botlayer-1].Vol -= rmndr
                tmpOutE = rmndr*self.Layers[this_botlayer-1].Temp
                totOutE += tmpOutE
                totOutQ += rmndr
                self.Layers[this_botlayer-1].TotE = self.Layers[this_botlayer-1].Vol*self.Layers[this_botlayer-1].Temp
        
    # check again if outflow balance is correct
    if abs(totOutQ - assignedOutQ) > 0.0001:
        raise ValueError(f"Specified outflow {assignedOutQ} and calcualted outflow {totOutQ} don't match!!")
    if totOutQ>0:
        outTemp = totOutE/totOutQ
    else:
        outTemp = np.nan
    return([totOutQ, outTemp, totOutE, q_by_outlet])
  
def outflow_dist_forGateSelect(resmod, qoutLyr):
    '''
     version of outflow_dist meant for iterative calculation using gate 
     selection processes - i.e. no water is removed from reservoir but outflow
     temperature is calculated
     
     
     this function takes a dictionary (key: layer index, value: outflow in AF)
     and does the accounting to remove water from the reservoir and calculate
     the effective outflow temperature
    '''
    totOutE = 0.
    totOutQ = 0.
    rmndr = 0.
    for l, qo in qoutLyr.items():
        if qo > 0.:
            iqo = min(resmod.Layers[l].Vol, qo+rmndr)
            rmndr = qo - iqo
            #self.Layers[l].Vol -= iqo  # remove volume from layer
            tmpOutE = iqo* resmod.Layers[l].Temp  # calc the vol*temp for outflow from layer
            totOutE += tmpOutE  # accumulate vol*temp in total outflow
            totOutQ += iqo    # accumulate volume in total outflow
            #self.Layers[l].TotE = self.Layers[l].Vol*self.Layers[l].Temp
         
    outTemp = totOutE/totOutQ
    return([totOutQ, outTemp, totOutE])             
         
  
def extra_outflow_dist(resmod, extra_outflow_vol):
    '''
    Remove any additional water from the reservoir not going through controlled
    outlets (e.g. for direct use and/or delivery). 
    
    This functionality exists to ensure mass balance in cases where there are 
    direct diversions or use from a reservoir pool that affect the volume in
    the reservoir over time. It is assumed that the extraction of water occurs
    in the top-most active layer of the reservoir (i.e. from the surface).
    This could be adjusted in the future if needed. 
    The function updates the reservoir volume in the ResTem object

    Parameters
    ----------
    resmod : Object
        ResTemp reservoir simulation object.

    extra_outflow_vol: float
        Volume (in consistent model units [e.g. acre-feet]) of additional water
        to be removed from the reservoir; ass

    Returns
    -------
    None.

    '''
    topLyr = 0
    for l in resmod.Layers.keys():
        if resmod.Layers[l].Vol > 0:  # this is the top layer
            topLyr = l
            
    # check if the extraction is greater than the layer volume
    remain_out_vol = extra_outflow_vol
    incr_idx = 0 # index to increment down from the top layer if we need to exhaust one layer and go to the next
    while remain_out_vol > 0.0:
        thislyrvol = resmod.Layers[topLyr-incr_idx].Vol
        thislyr_remove = min(thislyrvol, remain_out_vol)
        resmod.Layers[topLyr-incr_idx].Vol = thislyrvol - thislyr_remove
        # update the total energy value as well
        resmod.Layers[topLyr-incr_idx].TotE = resmod.Layers[topLyr-incr_idx].Vol * resmod.Layers[topLyr-incr_idx].Temp
        remain_out_vol = remain_out_vol - thislyr_remove
        incr_idx +=1
            
        
#%%################  Removed code and/or scratch ##############################


# From tcd_alloc (shortened to more concise and generalizable loop)
#     if gateDict[3]>0:  # upper gates open
        
#         if gateDict[2]==0: # no middle gates open - just upper
            
#             # ???? 2020-06-26 - attempting to account for weird 2010 summer
#             # behavior when only upper gates are open and simulated release
#             # temperatures are high compared to actuals - indicating (in real life even)
#             # that much of the flow is coming from deeper
#             if resmod.TimeStepDate.month in [99]: #<--- ??? 2020-08-05: deactivating for testing Mike Deas' pre/post2010 logic ; [6, 7, 8, 9]:  
#                 ptSinkElevs3, fracs3 = distPointSinks(resmod, topElev, 3, nPtSinks)
#                 ptSinkElevs2, fracs2 = distPointSinks(resmod, topElev, 2, nPtSinks)
                
#                 quppfrac = 0.333 #min(1.,self.Outlets[3].BalanceFactor*0.5)
#                 gateQ[3] = [[qoutTCD*f*quppfrac for f in fracs3], ptSinkElevs3]
#                 gateQ[2] = [[qoutTCD*(1-quppfrac)*f for f in fracs2], ptSinkElevs2]
#             else:
#                 # just flow through the upper gates
#                 ptSinkElevs3, fracs3 = distPointSinks(resmod, topElev, 3, nPtSinks)
#                 #ptSinkElevs3 = [[1021., np.sqrt(32.2*(topElev-1021.))], [950., np.sqrt(32.2*(topElev-1000.))]]
#                 #fracs3 = [0.1, 0.9]
#                 gateQ[3] = [[qoutTCD*f for f in fracs3], ptSinkElevs3]
            

#         else:        # also a middle gate is open
#             #hdUpp = max(0., topElev - self.Outlets[3].CtrElevFt)          
#             #hdMid = max(0., topElev - self.Outlets[2].CtrElevFt)
#             ## calculate relative fraction based on assumption that flow
#             ## is proportional to hydrostatic head
#             #sqrthdratio = np.sqrt(hdUpp/hdMid)
#             #a1fac = gateDict[3]*sqrthdratio
            
#             ptSinkElevs3, fracs3 = distPointSinks(resmod, topElev, 3, nPtSinks)
#             ptSinkElevs2, fracs2 = distPointSinks(resmod, topElev, 2, nPtSinks)
#             # trying out simple assumption (i think i did the math right??)
#             # that the fraction of flow to the top gate is 1/(numLowerGates+1)
#             if topElev > resmod.Outlets[3].BotElevFt:
#                 #quppfrac = min(1.,self.Outlets[3].BalanceFactor*(1/(gateDict[2]+1)))
#                 #uppfracopts = [0.9, 0.8, 0.5, 0.3, 0.1]
# #                    if gateDict[2]>1:
# #                        #quppfrac =  uppfracopts(gateDict[3]) 
# #                        quppfrac = min(1., self.Outlets[3].BalanceFactor*(gateDict[3]/totGatesOpen))
# #                    else:
# #                        quppfrac = 0.9
#                 quppfrac = glvlHdVels[3][3]
#             else:
#                 quppfrac = 0.
#             #elev_wght = max(0.,(topElev - self.Outlets[3].CtrElevFt)/100.)
#             #quppfrac = min(1.,(gateDict[3]/totGatesOpen)*elev_wght*self.Outlets[3].BalanceFactor)
#             gateQ[3] = [[qoutTCD *quppfrac*f for f in fracs3], ptSinkElevs3]
#             gateQ[2] = [[qoutTCD*(1-quppfrac)*f for f in fracs2], ptSinkElevs2]
            
#     elif gateDict[2]>0: # assumes top gates are closed/dry and only pressure relief gates can be also open
#         if gateDict[1]==0: # pressure relief gates are closed - just middle
#             ptSinkElevs2, fracs2 = distPointSinks(resmod, topElev, 2, nPtSinks)
#             #ptSinkElevs1, fracs1 = self.distPointSinks(topElev, 1, nPtSinks)
            
#             quppfrac = 1
# #                if self.TimeStep.month in [6, 7, 8, 9]:
# #                    quppfrac = 0.5
#             gateQ[2] = [[qoutTCD * f*quppfrac for f in fracs2],ptSinkElevs2]
#             #gateQ[1] = [[qoutTCD* f * (1-quppfrac) for f in fracs1], ptSinkElevs1]
#         else:    # also at least one PRG is open
#             #quppfrac = min(1., self.Outlets[2].BalanceFactor*(1/(gateDict[1]+1)))
#             #quppfrac = min(1., self.Outlets[2].BalanceFactor*(gateDict[2]/totGatesOpen))
#             quppfrac = glvlHdVels[2][3]
            
#             ptSinkElevs2, fracs2 = distPointSinks(resmod,topElev, 2, nPtSinks)
#             ptSinkElevs1, fracs1 = distPointSinks(resmod,topElev, 1, nPtSinks)
#             #quppfrac = min(1.,(gateDict[2]/totGatesOpen)*self.Outlets[2].BalanceFactor)
#             gateQ[2] = [[qoutTCD*quppfrac*f for f in fracs2],ptSinkElevs2]
#             gateQ[1] = [[qoutTCD*(1-quppfrac)*f for f in fracs1], ptSinkElevs1]
        
#     elif gateDict[1]>0: # assumes top and middle gates are closed and only side gates can also be open
#         if gateDict[0] ==0: #side gates are not open
#             ptSinkElevs, fracs = distPointSinks(resmod,topElev, 1, nPtSinks)
#             gateQ[1] = [[qoutTCD*f for f in fracs], ptSinkElevs]
        
#         else:   #at least one side gate is also open
#             ptSinkElevs1, fracs1 = distPointSinks(resmod,topElev, 1, nPtSinks)
#             ptSinkElevs0, fracs0 = distPointSinks(resmod,topElev, 0, nPtSinks)
            
#             #quppfrac = min(1., self.Outlets[1].BalanceFactor*1/(gateDict[0]+1))  #???? should this be adjusted to account for fewer gates/area at side gate?
#             #quppfrac = min(1.,(gateDict[1]/totGatesOpen)*self.Outlets[1].BalanceFactor)
#             #quppfrac = min(1., self.Outlets[1].BalanceFactor*(gateDict[1]/totGatesOpen))
#             quppfrac = glvlHdVels[1][3]
            
#             gateQ[1] = [[qoutTCD*quppfrac*f for f in fracs1], ptSinkElevs1]
#             gateQ[0] = [[qoutTCD*(1-quppfrac)*f for f in fracs0], ptSinkElevs0]
        
#     else:
#         ptSinkElevs0, fracs0 = distPointSinks(resmod,topElev, 0, nPtSinks)
#         ptSinkElevsD, fracsD = distPointSinks(resmod,topElev, 99, 1) #<-- lower level extra outlet only activated when sidegates are opend
        
#         ptSinkELevs0 = [730.,735., 740. ]
#         # just the side gates are open
#         gateQ[0] = [[qoutTCD*f*1 for f in fracs0], ptSinkElevs0]
#         gateQ[99] = [[qoutTCD*0.0], ptSinkElevsD] #<-- per Mike Deas Shasta draft report, page 72, 35% of tcd outflow goes through this added outlet
        
#     return(gateQ)
    