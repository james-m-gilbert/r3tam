# -*- coding: utf-8 -*-
"""
LongTemp: Longitudinally-Mixed Reservoir Model

Intended for simulations of shallow, well-mixed reservoirs with short (days)
travel/residence times typical of tailbay or reregulating reservoirs downtream
of larger dams and reservoirs.

@author: jgilbert
"""

import os
from r3tam.constants import ACREtoM2,  FTtoM, AFtoM3, FT3toM3, CFStoKG, Le, \
    cw, CMStoCFS, CFStoAFD, BOWEN_CONSTANT, TWOPI, HALFPI
import os, sys
import copy
import pandas as pnd
import numpy as np
from collections import OrderedDict as Odict
import time as time
import yaml
import r3tam.restemp_util as util
import r3tam.surface as surface

import pickle
import datetime as dt

class simspec:
    def __init__(self):
        self.StartDate_str = ''
        self.EndDate_str = ''
        self.DELT_DAY = 1
        self.DELT_SEC = self.DELT_DAY*86400.
        self.InitMethod = ''
        self.IsHindcast = True
        self.CalcGateOps = False
        self.GateChgFreq = 5 # how frequently can a gate change be made, in days
        self.CalcEvaporation = False  # calculate evaporation based on air temp, surface temp and wind?
        self.CalcReleaseSched = False
        self.Reinitialize = False
        self.ReinitMonth = 0
        self.ReinitDay = 0
        
    @property
    def StartDate(self):
        return(dt.datetime.strptime(self.StartDate_str, '%m/%d/%Y'))
        
    @property
    def EndDate(self):
        return(dt.datetime.strptime(self.EndDate_str, '%m/%d/%Y'))

class inputTS:
    
    def __init__(self):
        self.DataFrame = None
        self.ColumnMap = None

class observations:
    
    def __init__(self):
        self.OutflowDF = None
        self.OutflowColMap = None
        self.ProfilesDF = None
        self.ProfilesColMap = None

class LongTemp:
    
    def __init__(self):
        self.K = 0.01
        self.TTconst = 2.1  # in time units of model
        self.CriticalDepth = 30  # ft
        self.CriticalDepthFactor = 1.
        self.DiurnalTempOption = False
        self.MaxVol = np.nan # m3
        self.MinVol = 0.  # m3
        self.Area = 2476000.   #m2
        self.Area_units = 'm2'
        #self.Outlets = Odict()  # TCD outlets
        self.FacilityName = ''
        self.Storage = 0.
        self.WSE = 0.
        self.Rho = 999. # water density, kg/m3
        self.Temp = 55. # water temperature, deg F
        self.IntermedTemp = 55. # water temperature calculated before finalizing time step
        self.TotE = self.Storage*self.Temp
        self.Elevation_Area = None
        self.Elevation_Volume = None
        self.ReleaseTemps = []
        self.ReleaseVolumes = []
        self.StorageTemps = []
        self.StorageRecord = []
        self.Accretions = []
        self.Evaporation= []
        self.TribInflow = []
        self.MainInflow = []
        self.AreaRecord = []
        self.StageRecord = []
        self.LogFilePath = ''
        self.TimeStep = -1
        self._time = -1
        self.TimeStepDate = dt.datetime(2000, 10, 1, 0, 0)
        #self.ReservoirData = res_data()
        self.SimulationSpecs = simspec()
        self.Inflow = inputTS()
        self.Outflow = inputTS()
        self.Met = inputTS()
        self.Observations = observations()
        self.CalcMetric = True
        self.WindHeight_m = 2.
        self.Wind_Func_Params = {}
        self.SurfHeat = {'Evap_Wm2': [], 'Evap_m3':[], 'Evap_J':[], 'SolRad_Wm2': [], 'SolRad_J':[]}
        self.Debug = {'Release':0, 'General':0}
        self.DistrInflowTemps = {1: 5., 2: 5, 3: 5, 4: 6.0, 5: 6.5, 6: 7., 7:7.25, 8:7.5,
                                 9: 7.25, 10: 6.6, 11: 5.5, 12: 5}
        
    def __repr__(self):
        return(f'{self.__class__.__name__} - Longitudinal reservoir: {self.FacilityName}')
    
    @classmethod
    def initialize_longmod(cls, config_fp, **kwargs):
        '''
        Initialize the simulation options and temperature in the well-mixed reservoir.

        Parameters
        ----------
        config_fp : str
            Full filepath pointing to the LongTemp configuration file
        **kwargs : TYPE
            Options for LongTemp model setup

        Returns
        -------
        longmod : LongTemp
            An initialized LongTemp model object.

        '''
        
        print(config_fp)
        
        if config_fp[-3:]=='pkl':
            import pickle
            with open(config_fp,'rb') as infile:
                longmod = pickle.load(infile)
            runName = os.path.basename(config_fp)[0:-4] # assume runname is base filename
        
        if config_fp[-4:]=='yaml':
            inputsDict = util.getInputsDict(config_fp)
            #runName = inputsDict['RunName']
            projDir = os.path.dirname(config_fp)
            print("working directory set as: %s" %projDir)
            longmod = cls()  # instantiate reservoir model object
            
            util.setupSim(longmod, inputsDict, projDir, simType='long')
             
            # get the elevation-area-capacity curves
            elevareafp = os.path.join(projDir,inputsDict['ElevAreaTable']['Values'][0] )
            elevarea_units = inputsDict['ElevAreaTable']['Units']
            if not os.path.exists(elevareafp):
                raise FileNotFoundError(f"Couldn't find specified elev-area table {elevareafp}")
            
            elevarea =  util.getElevAreaTable(elevareafp, elev_col=0, area_col=1)
            if elevarea_units.upper() in ['FT','FEET']:
                # assumes elev-area curve is in ft-acres units - convert to m and m2
                area_m2 = [_*ACREtoM2 for _ in elevarea.Area]
                elev_m = [_*FTtoM for _ in elevarea.Elevation]
            longmod.Elevation_Area =pnd.DataFrame([elev_m, area_m2]).T
            longmod.Elevation_Area.columns=elevarea.columns
            
            # get the elevation-area-capacity curves
            elevvolfp = os.path.join(projDir,inputsDict['ElevStorTable']['Values'][0] )
            elevvol_units = inputsDict['ElevStorTable']['Units']
            if not os.path.exists(elevareafp):
                raise FileNotFoundError(f"Couldn't find specified elev-volume table {elevvolfp}")
                
            elevvol = util.getElevStorTable(elevvolfp, elev_col=0, stor_col=1)

            if elevvol_units.upper() in ['FT','FEET']:
                # assumes elev-area curve is in ft-acres units - convert to m and m2
                vol_m3 = [_*AFtoM3 for _ in elevvol.Storage]
                elev_m = [_*FTtoM for _ in elevvol.Elevation]
            longmod.Elevation_Volume =pnd.DataFrame([elev_m, vol_m3]).T
            longmod.Elevation_Volume.columns=elevvol.columns
            
            longmod.MaxVol = np.max(longmod.Elevation_Volume['Storage'])
            # ???? don't ahve any information on KWK outlets yet - not sure if needed
            #util.setupOutlets(kwk, inputsDict) 

            util.getTimeSeries(longmod, inputsDict)
            
            # check that dew point is provide if evaporation calculation is turned on
            if (longmod.SimulationSpecs.CalcEvaporation) and \
                ('dewpt' not in longmod.Met.ColumnMap):
                
                print("\n\nERROR: Evaporation calculation turned on, but no dew point time series provided!" )
                print("===================================================================================")
                return(None)
                
            
            util.getObservations(longmod, inputsDict)
            
            # ???? no operations data defined - again, not sure how much this is needed/will be needed
            #util.getOpsData(kwk, inputsDict)

            
            util.setLongParams(longmod, inputsDict)   
            
            if 'Calibration Mode' in inputsDict:
                print("Setting up calibration information")
                util.setupCalib(longmod, inputsDict)
        
        
        initDateStr = longmod.SimulationSpecs.StartDate_str
        initStoDate = longmod.SimulationSpecs.StartDate - dt.timedelta(1) # use init storage from end of prior day
        initStoDateStr = initStoDate.strftime('%Y-%m-%d')
        print("Initializing longitudinal reservoir on %s" %initDateStr)
        
        if longmod.SimulationSpecs.InitMethod =='from_obs':
            try:
                stoVar = longmod.Observations.OutflowColMap['storage']
                initVol = longmod.Observations.OutflowDF.loc[initStoDateStr,  stoVar]
                
                if stoVar.split('_')[-1].upper() in ['AF','ACRE-FEET','ACFT','ACREFEET']:
                    initVol = initVol*AFtoM3
                
                initTempVar = longmod.Observations.OutflowColMap['temperature'] #['stoTemp']
                initTemp = longmod.Observations.OutflowDF.loc[initDateStr,  initTempVar]
                
                if initTempVar.split('_')[-1].upper() in ['DEGF','DEG_F','F','DEGREE_F','DEGREEF', 'FAHRENHEIT']:
                    initTemp = (initTemp-32.)/1.8
                
            except:
                print("Couldn't initialize from the provided observation data")
                print(f" -- Check that storage for the date {initStoDateStr} is provided in observation file")
                return(BaseException("Ending setup"))
            #    return(longmod)
        else:

            # expects a set storage and temperature value from dictionary
            initVol = inputsDict['Initial Condition']['Storage']
            initVol_units = inputsDict['Initial Condition']['StorageUnits']
            
            if initVol_units.upper() in ['AF','ACRE-FEET','ACFT','ACREFEET']:
                initVol = initVol*AFtoM3
            
            
            initTemp_ = inputsDict['Initial Condition']['Temperature']
            
            initTemp_units = inputsDict['Initial Condition']['Temperature_Units']
            if initTemp_units.upper() in ['DEGF','DEG_F','F','DEGREE_F','DEGREEF', 'FAHRENHEIT']:
                initTemp = (initTemp_-32.)/1.8    
                
            # if longmod.Debug['General']>1:
            #     print("\n***************************************")
            #     print(f"\t Initializing from manual inputs: \n\t Storage: {initVol} AF  & Temperature: {initTemp_} {initTemp_units}")
            #     print("\n***************************************")
                
        longmod.Storage = initVol
        longmod.Temp = initTemp
        longmod.StorageRecord.append(initVol)
        longmod.StorageTemps.append(initTemp)
        longmod.Rho = longmod.calcRho(initTemp)
        longmod.TotE = initTemp*initVol
        longmod.update_WSE()
        longmod.update_area()
        longmod.AreaRecord.append(longmod.Area)
        longmod.StageRecord.append(longmod.WSE)
        longmod.MaxVol = np.max(longmod.Elevation_Volume['Storage'])
        
        if 'showInit' in kwargs:
            showInit=kwargs['showInit']
        else:
            showInit=False
        if showInit:
            if longmod.CalcMetric:
                print("Setting %s initial condition to: " %longmod.FacilityName)
                print("    Storage:   %0.2f cubic meters" %initVol)
                print("    Temperature:   %0.1f deg C" %initTemp)
            else:                
                print("Setting %s initial condition to: " %longmod.FacilityName)
                print("    Storage:   %0.2f acre-feet" %initVol)
                print("    Temperature:   %0.1f deg F" %initTemp)
        return(longmod)
    
    #@classmethod
    def set_inflow_temps(longmod, flow_temp_dict, set_date=None,
                         temp_units='DEG_F', flow_units='AF'):
        
        if set_date==None:
            thisdt = longmod.TimeStepDate+dt.timedelta(1)
        else:
            thisdt = set_date
        
        invec = longmod.Inflow.DataFrame.loc[thisdt]
        colmp = longmod.Inflow.ColumnMap
        # flow_temp_dict contains variable naems as keys, numbers to be set as values
        if 'inflow_final' in flow_temp_dict:
            #print(f"\t--updating {longmod.FacilityName} inflow to {flow_temp_dict['inflow_final']}")
            if flow_units=='AF':
                invec.loc[colmp['inflow_final']] = flow_temp_dict['inflow_final']*AFtoM3
            else:
                invec.loc[colmp['inflow_final']] = flow_temp_dict['inflow_final']
            
        if 'inflowTemp_final' in flow_temp_dict:
            #print(f"\t--updating {longmod.FacilityName} inflow temperature to {flow_temp_dict['inflowTemp_final']}")
            flow_temp_val = flow_temp_dict['inflowTemp_final']
            if np.isnan(flow_temp_val):
                flow_temp_val = 0.
            if temp_units=='DEG_F':
                invec.loc[[colmp['inflowTemp_final']]] = (flow_temp_val-32.)/1.8
            else:
                invec.loc[[colmp['inflowTemp_final']]] = flow_temp_val
        longmod.Inflow.DataFrame.loc[thisdt] = invec
        
    def advance_longtemp(longmod, final=True):
        
        if final:
            longmod._time +=1
            longmod.TimeStep += 1 #resmod.TimeStep
            longmod.TimeStepDate = longmod.SimDates[longmod.TimeStep]
            time_index = longmod._time 
            
            met_vec = longmod.Met.DataFrame.loc[longmod.TimeStepDate]
            inflow_vec = longmod.Inflow.DataFrame.loc[longmod.TimeStepDate]
            release_vec = longmod.Outflow.DataFrame.loc[longmod.TimeStepDate]  
            longmod.DayOfYear = longmod.TimeStepDate.timetuple().tm_yday
            
        else:
            tmp_time = longmod._time + 1
            tmp_TimeStep = longmod.TimeStep + 1
            tmp_TimeStepDate = longmod.SimDates[tmp_TimeStep]
            time_index = tmp_time
            
            met_vec = longmod.Met.DataFrame.loc[tmp_TimeStepDate]
            inflow_vec = longmod.Inflow.DataFrame.loc[tmp_TimeStepDate]
            release_vec = longmod.Outflow.DataFrame.loc[tmp_TimeStepDate]  
            longmod.DayOfYear = tmp_TimeStepDate.timetuple().tm_yday
        
        # check if needs to be reinitialized
        if longmod.SimulationSpecs.Reinitialize and longmod.TimeStep > 1:
            if (longmod.SimulationSpecs.ReinitMonth == longmod.TimeStepDate.month) & \
               (longmod.SimulationSpecs.ReinitDay == longmod.TimeStepDate.day):  # then reinitialize
                longmod.reintialize()
        
        #airT = met_vec[longmod.Met.ColumnMap['airTemp']]
        this_release = release_vec[longmod.Outflow.ColumnMap['outflow_final']]
        solrad = met_vec[longmod.Met.ColumnMap['solRad']]
        
        if type(longmod.Inflow.ColumnMap['tribInflowTemp'])==list:
            trib_temps = [inflow_vec[v] for v in longmod.Inflow.ColumnMap['tribInflowTemp_final']]
        else:
            trib_temps = [inflow_vec[longmod.Inflow.ColumnMap['tribInflowTemp_final']]]
            
        upTemp = inflow_vec[longmod.Inflow.ColumnMap['inflowTemp_final']]
        #inTemps = [sppt, shaobsT]
            
        #inQs = [extraKESinflow, totQ2]
        
        if type(longmod.Inflow.ColumnMap['tribInflow'])==list:
            trib_vols = [inflow_vec[v] for v in longmod.Inflow.ColumnMap['tribInflow_final']]
        else:
            trib_vols = [inflow_vec[longmod.Inflow.ColumnMap['tribInflow_final']]]
        
        upVol = inflow_vec[longmod.Inflow.ColumnMap['inflow_final']]        

        #inQs = [sppq, totQ2]

        sumQ, wTemp = longmod.calcFlowWeightTemp(trib_vols+[upVol], 
                                                 trib_temps+[upTemp])
        #print(f"Total inflow: {sumQ}   Inflow temp: {wTemp}")
        #print(sumQ/AFtoM3)
        #print(f"trib vol: {np.sum(trib_vols)/AFtoM3} AF | shasta vol: {upVol/AFtoM3} AF")
        
        if longmod.Inflow.ColumnMap['extraInflow_final'] !='':
            xtraInQ = inflow_vec[longmod.Inflow.ColumnMap['extraInflow_final']]
            
        else:
            if type(release_vec[longmod.Outflow.ColumnMap['storage']]) in [np.float64, float]:
                this_storage = release_vec[longmod.Outflow.ColumnMap['storage_obs']]
                xtraInQ = longmod.calcWatBal_inflow(sumQ, this_release, this_storage,
                                                    time_index, inUnits='M3', 
                                                    outUnits='M3', stoUnits='M3')
                if longmod.Debug['General']>1:
                    print(f"Calculated extra inflow of {xtraInQ:0.2f}")
                
            else:
                xtraInQ = 0.
                
        if final: # or len(longmod.Accretions)<(time_index-1):
            longmod.Accretions.append(xtraInQ)
            longmod.MainInflow.append(upVol)
            longmod.TribInflow.append(sum(trib_vols))
        else:
            if time_index==0 and len(longmod.Accretions)==0:
                longmod.Accretions.append(xtraInQ)
                longmod.MainInflow.append(upVol)
                longmod.TribInflow.append(sum(trib_vols))
            else:
                longmod.Accretions[time_index-1] = xtraInQ
                longmod.MainInflow[time_index-1] =upVol
                longmod.TribInflow[time_index-1] = sum(trib_vols)
                
        
        [final_temp, new_sto] = longmod.calcMixTempImp(sumQ, wTemp,this_release,
                                                       met_vec, time_index, 
                                                       stoCorrect=xtraInQ, 
                                                       calcTT=True, calcSolRad=True, 
                                                       radWm2=solrad, final=final)
        longmod.IntermedTemp = final_temp  
        longmod.IntermedRelease = this_release
        # longmod.calcMixTempImp(sumQ, wTemp,thisRelease,airT, stoCorrect=extraKESinflow, \
        #                    calcTT=True, calcSolRad=True, radWm2=solrad, final=True)
            
    def calcWatBal_inflow(self, inflow, outflow, storage, time_index, 
                          inUnits='cfs',outUnits='cfs', stoUnits='AF'):
        '''
            in cases where we know partial inflows, outflows, and storage,
            possibly from different sources, but inflows don't account for
            chagne in storage - this attempts to calc the missing water
        '''

        if inUnits.upper() in ['CFS','FT3/S','CU_FT_PER_SEC']:
            # convert to volume (acre-feet)
            invol = inflow* self.SimulationSpecs.DELT_SEC/43560.*AFtoM3
        elif inUnits.upper() in ['ACRE-FEET', 'ACFT','AC-FT', 'AF']:
            invol = inflow*AFtoM3
        elif inUnits.upper() in ['M3', 'CUBIC_METERS']:
            invol = inflow
        else:
            raise(BaseException('Not a recognized inflow unit, "%s"' %inUnits))
            return
    
        if outUnits.upper() in ['CFS','FT3/S','CU_FT_PER_SEC']:
            # convert to volume (acre-feet)
            outvol = outflow* self.SimulationSpecs.DELT_SEC/43560.*AFtoM3
        elif outUnits.upper() in ['ACRE-FEET', 'ACFT','AC-FT', 'AF']:
            outvol = outflow*AFtoM3
        elif outUnits.upper() in ['M3', 'CUBIC_METERS']:
            outvol = outflow
        else:
            raise(BaseException('Not a recognized inflow unit, "%s"' %outUnits))
            return    

        if stoUnits.upper() in ['AF','ACRE-FEET','AC-FT']:
            # conver to m3
            storage = storage*AFtoM3
        # if stoUnits.upper() in ['AF','ACRE-FEET']:
        #     # convert to volume (acre-feet)
        #     pass
        # else:
        #     raise(BaseException('Not a recognized inflow unit, "%s"' %stoUnits))
        #     return
        
        prevSto = self.StorageRecord[time_index-1]
        stoA = prevSto + invol - outvol # what storage would be without any adjustment
        balA = storage - stoA  # what adjustment is needed to make calc'd storage match prescribed?
        
        return(balA)
            
    def calcFlowWeightTemp(self, inQs, inTs):
        sumQ = sum(inQs)
        wTemp = sum([q*t for q,t in zip(inQs, inTs)])/sumQ
        if self.Debug['General']>1:
            print(self.TimeStepDate, f"Main inflow vol: {sumQ/AFtoM3} AF; Main inflow temp: {wTemp} deg C")
        return([sumQ,wTemp])
        
    def calcMixTempImp(longmod, inVol, inT,outVol, met_vec, time_index, stoCorrect=0, \
                       calcTT=False, calcSolRad=False,radWm2=0.,final=True):
        '''
            assume units are in acre-feet (volume) and cfs (flow)
        '''
        # if calcTT:
        #     outQ = outVol/longmod.SimulationSpecs.DELT_SEC # outVol*43560./longmod.SimulationSpecs.DELT_SEC # convert to cfs from AF
            
        #     # deal with situations where volume goes to zero/very low (as in CalSim runs)
        #     if longmod.Storage < longmod.MaxVol*0.333:
        #         thisTT = longmod.TTconst
        #     else:
        #         thisTT = (longmod.Storage / outQ)/86400.  # divid storage in m3 by outflow in m3/s, divide 86400 to get tt in days
        #     #thisTT = (longmod.Storage*43560. / outQ)/86400.  # divide storage in cubic feet by outflow in cfs, get travel time in seconds, convert to days
        # else:
        #     thisTT = longmod.TTconst
        
        # calculate an initial estimate of air-water warming 
        airT = met_vec[longmod.Met.ColumnMap['airTemp_final']]
        air_temp_delta1 = (airT - longmod.Temp)*longmod.C1
        
        if calcSolRad:
            #adjTemp = self.solrad3(radWm2, seasonal=True)
            solrad_temp_delta = longmod.solrad_w2(radWm2, final=final) #TODO: need to adjsut so that solrad isn't applied every iteration
            #adjTemp = longmod.Temp
        else:
            solrad_temp_delta = 0. #adjTemp = longmod.StorageTemps[time_index-1] #[self.TimeStep-1]
            
        if longmod.SimulationSpecs.CalcEvaporation:
            evap_temp_delta, evap_vol = longmod.calc_evap(met_vec, 
                                                          solrad_temp_delta, 
                                                          air_temp_delta1,
                                                          final=final)
            #delTemp_evap = -1*delTemp_evap
        else:

            evap_col = longmod.Met.ColumnMap['evapVol']
            evap_vol = met_vec[longmod.Met.ColumnMap['evapVol']]
            if evap_col.split("_")[-1].upper() in ['AF', 'ACFT','ACRE_FEET', 'AC-FT', 'AC_FT']:
                evap_vol = evap_vol*AFtoM3
            
            [evap_temp_delta, _] = longmod.calc_latentheat_from_evap(evap_vol) 

        
        
        #print(f"accretions: {stoCorrect/AFtoM3}")
        
        # update the storage to account for inflow addition and evaporation removal
        initSto = longmod.StorageRecord[time_index]
        stoA =  initSto + inVol + stoCorrect #-evap_vol - outVol - stoCorrect #[self.TimeStep-1] + inVol
        if stoCorrect>0:
            #distTemp = 7.5 #longmod.DistrInflowTemps[longmod.TimeStepDate.month] # # inT
            if longmod.TimeStepDate.month in [11,12,1,2]:
                distTemp = 7.5
            else:
                distTemp = longmod.Inflow.DataFrame.loc[longmod.TimeStepDate, 
                                                    longmod.Inflow.ColumnMap['inflowTemp_final']]
        else:
            #distTemp = longmod.Temp
            distTemp = longmod.Inflow.DataFrame.loc[longmod.TimeStepDate, 
                                                    longmod.Inflow.ColumnMap['inflowTemp_final']]
        tempA = (inVol*inT + initSto*longmod.Temp + stoCorrect*distTemp )/stoA  #+ max(stoCorrect,0.)*longmod.Temp
        

        # print(f"\nInvol: {inVol}")
        # print(f"InTemp: {inT}")
        # print(f"InitTemp: {longmod.Temp}")
        # print(f"Init storage: {longmod.StorageRecord[time_index-1]}")
        # print(f"Intermed storage: {stoA}")
        # print(f"Intermed Temp A: {tempA}")
        
        air_temp_delta2 = (airT - tempA)*longmod.C1
        
        tempB = tempA + solrad_temp_delta - (evap_temp_delta*longmod.C3)
        
        air_temp_delta3 = (airT - tempB )*longmod.C1
        
        if longmod.DiurnalTempOption:
            air_temp_delta = (air_temp_delta1+air_temp_delta2+air_temp_delta3)/3
        else:
            air_temp_delta = air_temp_delta1 
        
        # print(f"\nAir Temp delta 1: {air_temp_delta1}")
        # print(f"AirTemp delta2: {air_temp_delta2}")
        # print(f"AirTemp delta: {air_temp_delta}")
        # print(f"Solrad Temp Delta: {solrad_temp_delta}")
        # print(f"Evap temp delta: {evap_temp_delta}")


        if calcTT:
            outQ = outVol/longmod.SimulationSpecs.DELT_SEC # outVol*43560./longmod.SimulationSpecs.DELT_SEC # convert to cfs from AF
            
            # deal with situations where volume goes to zero/very low (as in CalSim runs)
            if longmod.Storage < longmod.MaxVol*0.1:
                thisTT = longmod.TTconst
            else:
                thisTT = (longmod.Storage / outQ)/86400.  # divid storage in m3 by outflow in m3/s, divide 86400 to get tt in days
            #thisTT = (longmod.Storage*43560. / outQ)/86400.  # divide storage in cubic feet by outflow in cfs, get travel time in seconds, convert to days
        else:
            thisTT = longmod.TTconst

        
        final_temp = tempB + (air_temp_delta)*longmod.K *thisTT
        
        prev_temp = longmod.StorageTemps[longmod.TimeStep] #ReleaseTemps[longmod.TimeStep]
        diff_temp = final_temp - prev_temp
        

        #tempA = (inVol*inT + self.StorageRecord[self.TimeStep-1]*adjTemp)/stoA
        #tempA = (inVol*inT + longmod.StorageRecord[time_index-1]*adjTemp)/stoA
        #tempB = tempA + (airT - tempA)*longmod.K*thisTT 
        
        
        #newSto = stoA - outVol + stoCorrect
        #newSto = max(16000., min(23500., newSto))
        newSto = min(longmod.MaxVol, stoA -evap_vol - outVol) # + stoCorrect)
        spill = stoA -evap_vol - outVol  - newSto #+ stoCorrect


        if newSto < 0.2*longmod.MaxVol:
            final_temp = tempA + (air_temp_delta)*longmod.K *thisTT # replaces existing calc when storage low

        if abs(diff_temp) > 6.0 or final_temp < 3.0 or final_temp > 1.5*tempA:
            final_temp = 0.5*tempA + 0.5*prev_temp

        if spill>0:
            newSto = newSto - spill
            print(f"\t**WARNING: Non-Zero Spill ({spill}) Calculated")
        #tempC = tempB #+ max(0., stoCorrect)*inT/(stoA-outVol)
        
        if final:
            longmod.ReleaseTemps.append(final_temp) #tempC)
            longmod.ReleaseVolumes.append(outVol+spill)
            longmod.StorageTemps.append(final_temp) #
            longmod.StorageRecord.append(newSto)
            longmod.Evaporation.append(evap_vol)
            #longmod.Accretions.append(stoCorrect) # handled in advance_longtemp

            longmod.Storage = newSto
            longmod.Temp = final_temp #
            longmod.TotE = final_temp*newSto #tempC*newSto
            
            # update WSE and area
            longmod.update_WSE()
            longmod.update_area() 
            
            return([final_temp, newSto])
        else:
            return([final_temp, newSto])
        
    def reintialize(longmod, initVol=None, initTemp=None, **kwargs):
        init_stor_date = longmod.TimeStepDate - dt.timedelta(1)
        init_stor_date_str = init_stor_date.strftime('%Y-%m-%d')
        
        if 'reset_timestep' in kwargs:
            if kwargs['reset_timestep']:
                longmod.TimeStep = -1
                longmod._time = -1
        
        if 'sim_dates' in kwargs:
            longmod.SimDates = kwargs['sim_dates']
            longmod.TimeStepDate = longmod.SimDates[0]
            
        if longmod.SimulationSpecs.InitMethod =='from_obs':
            try:
                stoVar = longmod.Observations.OutflowColMap['storage']
                initVol = longmod.Observations.OutflowDF.loc[init_stor_date_str,  stoVar]
                
                if stoVar.split('_')[-1].upper() in ['AF','ACRE-FEET','ACFT','ACREFEET']:
                    initVol = initVol*AFtoM3
                
                initTempVar = longmod.Observations.OutflowColMap['temperature'] #['stoTemp']
                initTemp = longmod.Observations.OutflowDF.loc[init_stor_date_str,  initTempVar]
                
                if initTempVar.split('_')[-1].upper() in ['DEGF','DEG_F','F','DEGREE_F','DEGREEF', 'FAHRENHEIT']:
                    initTemp = (initTemp-32.)/1.8
                
            except:
                print("Couldn't initialize from the provided observation data")
                print(f" -- Check that storage for the date {init_stor_date_str} is provided in observation file")
                return(BaseException("Ending setup"))
        else:
            # assume initilization is from specified values
            initVol = initVol
            initTemp = initTemp

        print(f"Reinitializing longitudinal model for {longmod.FacilityName} on {init_stor_date_str} ")
        longmod.Storage = initVol
        longmod.StorageRecord[longmod.TimeStep] = initVol
        longmod.Temp = initTemp

    def initialize(self, storageAF, tempDegF):
        
        self.Storage = storageAF
        self.StorageTemps.append(tempDegF)
        
    def solrad3(self, radWm2, seasonal=False):
        # UPDATE: modified to use radiation directly as W/m2 - to get the same
        # results as with solrad2, coefficient C2 will need to be smaller 
        # with this method
        # heat exchange from solar radiation
        # assumes distribution factors have been calculated previously
        # based on a prescribed critical depth
        # convert solar radiation in W/m2 into cal/m2/day

        if seasonal or (seasonal>0.):
            doy = pnd.to_datetime(self.TimeStepDate).dayofyear
            declination_angle = np.arcsin(0.39795 * np.cos((0.98563* (doy-173))*np.pi/180.))
            frac_rad = (180*declination_angle/np.pi/23.5+1.)/2.
            if seasonal >0:
                frac_rad = frac_rad * (1. - seasonal)+ seasonal
        else:
            frac_rad = 1.
        
        critDepth_m = self.CriticalDepth*FTtoM
        
        tc = (self.StorageTemps[self.TimeStep-1]-32.)*5/9.
        thisRho = 999.85050 + 0.06001*tc - 0.007917*tc**2 + 4.1256e-5*tc**3
        depFac = 0.18
        convfac = 86400. * 1/cw * 1/thisRho
        # #Area_m2 = 2476000. #<-- estimated from quick trace in QGIS
        # if self.Area_units.upper() not in ['M2']:
        #     Area_m2 = self.Area*ACREtoM2
        # else:
        #     Area_m2 = self.Area 
        Area_m2 = self.Area #<-- assumes area is already in m2
        
        E2 = critDepth_m*depFac* frac_rad*self.C2 * radWm2 * convfac * Area_m2/1000. * self.SimulationSpecs.DELT_DAY * self.CriticalDepthFactor #  since radiation is energy per time unit    
        TotE = E2*1.8*1000./AFtoM3 # convert from thousand m3 * degC to ac-ft*degF
        NewTemp = (TotE + self.Storage*self.StorageTemps[self.TimeStep-1])/self.Storage
        
        return(NewTemp)
    
    
    def solrad_w2(longmod, radWm2, final=False):
        """
    
        Parameters
        ----------
        longmod : Object
            LongTemp reservoir simulation object
            
        radWm2 : float
            Solar radiation flux, W/m2
    
        Returns
        -------
        LongTemp object with profile warmed by effects of solar radiation 
    
        Notes
        -----
        This function applies the solar radiation provided in the ``metvec`` [units: :math:`W/m^2`]
        to the surface and for each underlying layer until solar radiation has been
        fully absorbed. 
    
    
        """
        
        
               
        albedo = longmod.MeanAlbedo
        if longmod.SeasonalAlbedo:
            if longmod.Latitude > 5.: # per GLM (Hipsey et al 2017 for the General Lake Model (GLM 2.4)), if latituate greater than 5 North, then northern hemisphere
                albedo = longmod.MeanAlbedo + longmod.AlbedoAmplitude*np.sin((TWOPI*longmod.DayOfYear/365.)+HALFPI)
            elif longmod.Latitude < -5: #similarly, check if in souther hemisphere
                albedo = longmod.MeanAlbedo + longmod.AlbedoAmplitude*np.sin((TWOPI*longmod.DayOfYear/365.)-HALFPI)
            else:
                pass # near the equator - assume albedo is constant at mean albedo
                
        # net shortwave entering the water is incident minus the reflected amount
        sw_reflect = radWm2*albedo 
        sw_water = radWm2 - sw_reflect
        
        # sediment (at lake bottom) reflectivity
        # NOTE(20220831): assumed for now that reflectivity is 1 and that all solar
        #                 radiation incident to sediment is reflected back and absorbed
        #                 by the layer through which it just passed
        sed_reflect = 1.
        
        # a certain amount of sw radiation is absorbed at the waters surface - 
        # controlled by user-specified fraction (sw_beta)
        # do the calc for the top layer first - then loop through lower layers until
        # the radiation making it through the layer is zero
        topThick_m = 0.6
        topThick_ft = topThick_m/FTtoM

        # # Area_m2 = 2476000. #<-- estimated from quick trace in QGIS
        # if longmod.Area_units.upper() not in ['M2']:
        #     Area_m2 = longmod.Area*ACREtoM2
        # else:
        #     Area_m2 = longmod.Area 
        Area_m2 = longmod.Area #<-- assumes area already in m2
        
        convfac =1/cw * 1/longmod.Rho
        solrad_J = longmod.C2*(sw_water*Area_m2)*longmod.SimulationSpecs.DELT_SEC# W/m2*area_m2-> W=> J/s * Sec/time step = J 
        volE = solrad_J*convfac # J * 1/kg/m3 * 1/J/kg*C => m3*C (volume * temp change) 
        
        # ???? There's an issue when storage gets low, solrad temperature influence
        # increases substantially; in real life, Keswick storage does not go to/near 
        # zero, but we can surmise that under such conditions warming would be 
        # similar to that in a river; how best to approximate this condition?
        # version below is 
        if longmod.StorageRecord[longmod._time] < 0.1*longmod.MaxVol:
            solrad_J = solrad_J*0.25
            volE = solrad_J*convfac
        
        if longmod.CalcMetric: #then storage is in m3 and no conversion is necessry
            solrad_temp_delta = volE/longmod.Storage #<-- [m3*degC]/[m3] ==> units of degC
        else:
            solrad_temp_delta = 1.8*volE/(longmod.Storage*AFtoM3) #<-- [m3*degC]/[AF*m3/AF] ==> units of degC*1.8 ==> deg F

            
        #longmod.TotE += volE*1.8/AFtoM3 # convert from m3 * degC to ac-ft*degF #TODO - units check
        #longmod.Temp = longmod.TotE/longmod.Storage
        if final:
            longmod.SurfHeat['SolRad_Wm2'].append(sw_water)
            longmod.SurfHeat['SolRad_J'].append(solrad_J)
        else:
            if len(longmod.SurfHeat['SolRad_Wm2'])==0:
                longmod.SurfHeat['SolRad_Wm2'].append(sw_water)
                longmod.SurfHeat['SolRad_J'].append(solrad_J)
            else:
                #overwrite the value for this time step
                longmod.SurfHeat['SolRad_Wm2'][longmod.TimeStep] = sw_water
                longmod.SurfHeat['SolRad_J'][longmod.TimeStep] = solrad_J
       
        return(solrad_temp_delta)

    def calc_evap(longmod, met_vec, solrad_temp_delta, air_temp_delta,
                  TairUnits='degF', debug=0, final=False): # ws, Tdew, Tair, TairUnits='degF', debug=0):
        '''
            calculate evaporation from open water
            - adapted from ce-qual-w2 v4.2, heat-exchange-f90
            
            ws: wind speed, m/s
            Tdew: dew point temperture in deg c
            Tair: air temperature (ostensibly at 2 m), deg C(?)
        '''
        
        ws =  met_vec[longmod.Met.ColumnMap['wind']]
        Tdew = met_vec[longmod.Met.ColumnMap['dewpt']]
        Tair = met_vec[longmod.Met.ColumnMap['airTemp_final']]
        TairUnits = longmod.Met.ColumnMap['airTemp_final_units']
        

        
        if longmod.CalcMetric:
            longmod.calc_ea_es(met_vec,solrad_temp_delta, air_temp_delta,
                               TsurfUnits='DEG_C' )
        else:               
            longmod.calc_ea_es(met_vec,solrad_temp_delta,air_temp_delta,
                               TsurfUnits='DEG_F')
        # ea_mmHg calculated with call to calc_ea_es prior to calc_evap
        ea_mmHg = longmod.Met.VapPress['ea_mmHg']
        es_mmHg = longmod.Met.VapPress['es_mmHg']
        #ea_mmHg = np.exp(2.3026*((7.5*Tdew)/(Tdew+237.3)+0.6609))  # from CE-QUAL-W2 user's manual
        #ea_hPa = 6.11*10**((7.5*Tdew)/(237.3+Tdew))   # from: https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf
         
        # area should already be in units of m2
        Area_m2 = longmod.Area        
        
    
        fw = surface.wind_func(longmod, met_vec, method='W2')
        
        evap_Wm2 = max(0., fw*(es_mmHg - ea_mmHg))
        
        evapJ = evap_Wm2*Area_m2*longmod.SimulationSpecs.DELT_SEC
        evap_vol_m3 = evap_Wm2 * longmod.SimulationSpecs.DELT_SEC/Le * 1/1000. * Area_m2
        
        
        if evap_Wm2 >1500.:
            print(f'{longmod.TimeStepDate.strftime("%Y-%m-%d")}--High evaporation calculated:')
            print(f'   surface vapor pressure, mmHg: {es_mmHg: 0.4f}')
            print(f'   air vapor pressure, mmHg: {ea_mmHg: 0.4f}')
            #print(f'   tairv: {tairv: 0.4f}')
            #print(f'   dtv: {dtv: 0.4f}')
            #print(f'   dtvl: {dtvl: 0.4f}')
            print(f'   wind func: {fw: 0.4f}')
            print(f'   evap, W/m2: {evap_Wm2: 0.4f}')
            print(f'   evap volume (AF): {evap_vol_m3/AFtoM3}')
            #evap_Wm2 = 800.
            
       
        if final:
            longmod.SurfHeat['Evap_Wm2'].append(evap_Wm2)
            longmod.SurfHeat['Evap_J'].append(evapJ)
            longmod.SurfHeat['Evap_m3'].append(evap_vol_m3)  
            evap_vol = evap_vol_m3
        else:
            evap_vol = evap_vol_m3 
                        
        #evap_vol = longmod.SurfHeat['Evap_m3'][longmod._time]
        if not longmod.CalcMetric:
            # not CalcMetric means units are in AF-degF
            evap_vol = evap_vol/ AFtoM3  

        
        if longmod.CalcMetric:
            masswater = longmod.Storage*longmod.Rho
        else:
            masswater = longmod.Storage * AFtoM3*longmod.Rho #Revised from: 1000. # assumign density = 1000 kg/m3
        
        evap_temp_delta = evapJ/(masswater*cw)
        
        if not longmod.CalcMetric:
            evap_temp_delta = 1.8*evap_temp_delta
        #evapTempFDelta = 1.8*(evapJ/(masswater*cw))
        
        if debug > 0:
            if longmod.CalcMetric:
                print("Evaporation effect on top layer: %0.1f degrees C" %evap_temp_delta)
            else:
                print("Evaporation effect on top layer: %0.1f degrees F" %evap_temp_delta)
        
        return([evap_temp_delta, evap_vol])


    def calc_ea_es(longmod, met_vec,solrad_temp_delta, air_temp_delta, 
                   TsurfUnits='DEG_F'):
        """
        
    
        Parameters
        ----------
        resmod : Object
                 ResTemp reservoir simulation object
            
        met_vec : pandas dataframe row - indexable by column name
            Listing of meteorological variable values for this time step
            
        solrad_temp_delta: float
            Incremental warming (in deg C) caused by absorption of solar radiation
            applied to water prior to calculating vapor pressure to account
            for warming of the surface
            
        TsurfUnits : string, optional
            Indicate the units of the water surface temperature. The default is 'degF'.
    
        Returns
        -------
        None.
    
        """
        Tdew = met_vec[longmod.Met.ColumnMap['dewpt']]
        if TsurfUnits == 'DEG_F':
            Tsurf = (longmod.Temp-32.)*5/9.
        else:
            Tsurf = longmod.Temp

        Tsurf += solrad_temp_delta+air_temp_delta
            
        ea_mmHg = np.exp(2.3026*((7.5*Tdew)/(Tdew+237.3)+0.6609))  # from CE-QUAL-W2 user's manual
        
        if Tsurf < 0.:  # if the surface happens to be frozen
            es_mmHg = np.exp(2.3026*((9.5*Tsurf)/(Tsurf+237.3)+0.6609))
        else:
            es_mmHg = np.exp(2.3026*((7.5*Tsurf)/(Tsurf+237.3)+0.6609))
        longmod.Met.VapPress = {'ea_mmHg': ea_mmHg, 'es_mmHg': es_mmHg}            

    def remove_evap(longmod, method='simple',debug=0, final=True):
        """
        TODO: update docstring for longtemp
        
        Function to remove evaporated water from the reservoir/lake surface and
        to distribute the associated heat exchange through the surface layers.
        Based on the formulation described in [[1]_]
        
        Parameters
        ----------
        resmod : Object
                 ResTemp reservoir simulation object
        method : string, optional
            Switch to pick different methods of distributing the latent heat
            energy associated with evaporation. As of Aug 2022, only the 'simple'
            method is implemented - other entries will The default is 'simple'.
        debug : TYPE, optional
            Sets the 'debug' level for writing out process information - a higher 
            number means more information will be written to the console, and run
            times will be slower. The default is 0.


        Returns
        -------
        ResTemp object with evaporated water removed and profile temperature updated
        by effects of latent heat exchange (evaporation) 

        Notes
        -----
        This function approximates the heat exchange that results when water evaporates
        from the surface of a water body, based on the formulation
        described in [[1]_]. The function assumes that evaporation volume has already
        been calculated and is stored in the ``resmod.SurfHeat['Evap_m3']`` variable.
        This volume is removed from the top layer. If the evaporated volume is greater
        than the top layer volume, the remainder is removed from the next lower layer.
        
        The latent heat energy associated with the evaporated volume (calculated in function
        ``calc_latentheat_from_evap``) is then converted to an equivalent temperature
        change by dividing by layer water mass and specific heat capacity.
        This temperature change is applied to the layers within the critical
        depth using recalculated (to account for layer thickness changes from evaporation
        removal) distribution factors, adjusted by user-specified parameter **C3**.


        Following [[1]_], the formulation used for this function is:
            
            .. math::
                
                E_{LatentHeat} = \sum_{l=top}^{l=extinct} frac_l* C_3* (L_e * m_{evap})/(m_l *c_w)*S_l
                
        Where
        
        :math:`E_{LatentHeat} =` energy transfer through the surface mixing zone 
        assocaited with evaporation, of dimension volume*temperature [units: acre-feet*°F or m3*°C]
        
        :math:`l =` layer index, where top denotes to the top (active) layer in the 
        reservoir and extinct denotes the layer intersected by the extinction depth
            
        :math:`frac_l =` Energy distribution fraction assigned to layer l. According 
        to the Reclamation temperature model formulation [1], this fraction “decreases 
        linearly with depth from 1 at the surface to 0 at the bottom [of the specified 
        depth of energy penetration]”
            
        :math:`C_1 =` User-specified calibration parameter that adjusts the magnitude 
        or efficiency of heat transfer resulting from evaporation; permissible values 
        in range 0-1 [no units]
            
        :math:`L_e =` Latent heat of vaporization for water, approximated as constant 
        2.453 x 106 J/kg at 20°C [units: J/kg]
            
        :math:`m_{evap}=` Mass of evaporated water – either from user-specified 
        inputs or calculated by the model (see ``calc_evap`` function); [units: kg]
            
        :math:`m_{l}=` Mass of water in layer l, calculated as the product of volume 
        and density for layer l [units: kg]
        
        :math:`S_l=` Volume of water in layer *l* [units: acre-feet or :math:`m^3`]
        
        :math:`TCF=` Temperature conversion factor: A value of 1.8 is used to convert 
        from the temperature change (in K or equivalently °C) represented by the 
        :math:`(L_e*m_{evap})/(m_l*c_w)` term to °F (assumes :math:`S_l` is in acre-feet); 
        otherwise use a value of 1
                
        .. [1] Rowell, J. H. (1990). U.S. Bureau of Reclamation Monthly Temperature 
               Model Sacramento River Basin: Draft Report. U.S. Bureau of Reclamation, 
               Mid-Pacific Region.

        """

        # remove evaporated water from top layer volume
        
        evap_vol_AF = longmod.SurfHeat['Evap_m3'][longmod._time] / AFtoM3  # _time index minus one to acct for stored surfheat data being zero-indexed
        evapJ = longmod.SurfHeat['Evap_J'][longmod._time]   

        if final:
            longmod.Storage -= evap_vol_AF 
        
        #TODO: re-write this for metric, check units
        masswater = longmod.Storage * AFtoM3*longmod.Rho #Revised from: 1000. # assumign density = 1000 kg/m3
        evapTempFDelta = 1.8*(evapJ/(masswater*cw))
        if debug > 0:
            print("Evaporation effect on top layer: %0.1f degrees F" %evapTempFDelta)
        
        return([evapTempDelta, evapVol])
        
        
        # if method=='simple':
        #     # first version is for applying latent heat exchange only to top layer
        #     #evapTotE = resmod.C3 * resmod.Layers[topLyr].Vol * evapTempFDelta  # incremental Kelvin is equivalent to incremental Celsius
        #     #resmod.Layers[topLyr].TotE -= evapTotE  # assuming that latent heat is removed from the remaining volume in top layer
        #     #resmod.Layers[topLyr].Temp = resmod.Layers[topLyr].TotE/resmod.Layers[topLyr].Vol
        #     for lyr, depfac1, depfac2 in resmod.CritDepthDist:
        #         thisLayerMass = resmod.Layers[lyr].Vol*AFtoM3*resmod.Layers[lyr].Rho
        #         evapTempFDelta = 1.8*(evapJ/(thisLayerMass*cw))
        #         evapTotE = depfac1*resmod.C3 * resmod.Layers[lyr].Vol * evapTempFDelta  # TODO: change to Celsius - incremental Kelvin is equivalent to incremental Celsius
        #         resmod.Layers[lyr].TotE -= evapTotE  # assuming that latent heat is removed from the remaining volume in top layer
        #         resmod.Layers[lyr].Temp = resmod.Layers[lyr].TotE/resmod.Layers[lyr].Vol
        # else:
        #     # under energy budget approach, heat is distributed through top layers
        #     # according to net of all energy terms
        #     raise Exception("Method other than 'simple' not yet implemented")  #pass
           

        
    def calc_latentheat_from_evap(longmod,evap_vol_m3): #self, evapCFS, debug=False):
        
    #   from constants.py
    #    Le = 2.453e6   # latent heat of vaporization, J/kg, at 20deg C
    #    cw = 4182.  # specific heat of water at 20 deg C, in J/kg/K
    #   cfs2kg = 1/CMStoCFS * 86400.*1000. # convert from cfs to cubic meters/sec, seconds to days, to kg using assumed density of 1000. kg/m3
        # evapCFS = met_vec[resmod.Met.ColumnMap['evap']]
        # evapJ = evapCFS*CFStoKG *Le
        # evapAF = evapCFS*CFStoAFD
        # evap_vol_m3 = evapAF*AFtoM3
        evapJ = evap_vol_m3*longmod.Rho *Le
        
        topArea_m2 = longmod.Area #area should be in m2
        
        evap_Wm2 = evapJ/(topArea_m2*longmod.SimulationSpecs.DELT_SEC)
        
        if longmod.CalcMetric:
            masswater = longmod.Storage*longmod.Rho
        else:
            masswater = longmod.Storage * AFtoM3*longmod.Rho #Revised from: 1000. # assumign density = 1000 kg/m3
        
        evap_temp_delta = evapJ/(masswater*cw)
        
        longmod.SurfHeat['Evap_Wm2'].append(evap_Wm2)
        longmod.SurfHeat['Evap_J'].append(evapJ)
        longmod.SurfHeat['Evap_m3'].append(evap_vol_m3)
        
        return([evap_temp_delta, evap_vol_m3])
            
    def calcRho(self, tc):
        # given temp in C (tc), calulate density
        rho = 999.85050 + 0.06001*tc - 0.007917*tc**2 + 4.1256e-5*tc**3
        return(rho)
    
    def update_WSE(self):
        # given a volume, calculate the state or water surface elevaiton
        elevs = self.Elevation_Volume['Elevation']
        vols = self.Elevation_Volume['Storage']
        self.WSE = np.interp(self.Storage, vols, elevs)
        
        if self.TimeStep <=0 or self.TimeStep>len(self.StageRecord)-1:
            self.StageRecord.append(self.WSE)
        else:
            self.StageRecord[self.TimeStep] = self.WSE
        
    def update_area(self):
        
        self.update_WSE() # update WSE first
        # given a WSE, calcualte the surface area
        elevs = self.Elevation_Area['Elevation']
        areas = self.Elevation_Area['Area']
        
        self.Area = np.interp(self.WSE, elevs, areas)
        
        if self.TimeStep <=0 or self.TimeStep>len(self.AreaRecord)-1:
            self.AreaRecord.append(self.Area)
        else:
            self.AreaRecord[self.TimeStep] = self.Area
        
    def calc_release_temp_errors(self, year=[]):
        yrs = list(set([i.year for i in self.SimDates]))
        
        if len(year)>0:
            yrs = [y for y in year if y in yrs]
                   
        error_dict = {}

        simT = pnd.DataFrame(self.ReleaseTemps, index=self.SimDates, 
                             columns=['SimReleaseT_degC'])
        outTcol = self.Observations.OutflowColMap['obs_outtemp_final']
        obsT = pnd.DataFrame(self.Observations.OutflowDF.loc[self.SimDates[0:len(simT)],
                                                             outTcol[0]])

        for y in yrs:
             
            sim_data = simT.loc[str(y),:] #kwk.calc_errors_df[y][ivar]['sim']
            obs_data = obsT.loc[str(y), :] #resmod.calc_errors_df[y][ivar]['obs']
            
            monerr = {m:[] for m in range(1,13)}  # blank dictionary to hold monthly error calcs
            monobs = {m:[] for m in range(1,13)} # blank dictionary for the observations for purposes of calculating NSE
            monsim = {m:[] for m in range(1,13)} # blank dictionary for the simulation data
            month_stats = {}
            for m in range(1,13):
                thismidx = [ix for ix in sim_data.index if ix.month==m]

                # calculate for both outflow volume and temperature
                # this_sim_vol = sim_data.loc[thismidx,sim_data.columns[0]]
                # this_obs_vol = obs_data.loc[thismidx,obs_data.columns[0]]
                
                this_sim_temp = sim_data.loc[thismidx, sim_data.columns[0]]
                this_obs_temp = obs_data.loc[thismidx, obs_data.columns[0]]
                
                err_temp = [s-o for s,o in zip(this_sim_temp, this_obs_temp)]
                sumsqerr = np.sum([er*er for er in err_temp])
                
                nse = 1 - sumsqerr/(np.nansum((this_obs_temp-np.nanmean(this_obs_temp))**2))
                
                mean_bias = np.sum(err_temp)/len(err_temp)
                mae = np.sum([abs(i) for i in err_temp])/len(err_temp)
                
                rmse = np.sqrt(sumsqerr/len(this_obs_temp))
                corrmat = np.corrcoef(this_obs_temp,pnd.to_numeric(this_sim_temp))
                corrxy = corrmat[0,1]
                r2 = corrxy**2
                
                month_stats[m] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                                 'mae': mae}
                

            # calc on an annual basis too
            err_temp = [s-o for s,o in zip(sim_data.iloc[:,0], obs_data.iloc[:,0])]
            sumsqerr = np.sum([er*er for er in err_temp])
            
            nse = 1 - sumsqerr/(np.nansum((obs_data.iloc[:,0]-np.nanmean(obs_data.iloc[:,0]))**2))
            
            mean_bias = np.sum(err_temp)/len(err_temp)
            mae = np.sum([abs(i) for i in err_temp])/len(err_temp)
            
            rmse = np.sqrt(sumsqerr/len(obs_data.iloc[:,0]))
            corrmat = np.corrcoef(obs_data.iloc[:,0],pnd.to_numeric(sim_data.iloc[:,0]))
            corrxy = corrmat[0,1]
            r2 = corrxy**2
            month_stats['annual'] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                             'mae': mae}
            error_dict[y] = month_stats
        
        return(error_dict)
    
    def finalize(self):
        '''
        Creates a dataframe with simulated output variables. Converts all flow/storage
        units to acre-feet volume

        Returns
        -------
        Dataframe with compiled time series

        '''

        columns = ['ReleaseTemps_degC', 'ReleaseTemps_degF',
                   'ReleaseVolumes', 'StorageTemps',
                  'Storage', 'Accretions', 'Evaporation','TribInflow',
                  'MainInflow', 'Area_ac'] #TODO - add area and stage, 'Area','Stage']

        index = pnd.date_range(self.SimulationSpecs.StartDate-dt.timedelta(1), 
                               self.TimeStepDate)
        
        # prepend [np.nan] for variables that aren't initialized so the 
        # dataframe is aligned with those that are
        releasetemps_degC = [np.nan] + [k for k in self.ReleaseTemps] #[np.nan] + 
        releasetemps_degF = [np.nan] + [k*1.8+32 for k in self.ReleaseTemps]
        releasevol_af = [np.nan] + [k/AFtoM3 for k in self.ReleaseVolumes]
        storage_af = [k/AFtoM3 for k in self.StorageRecord]
        accret_af = [k/AFtoM3 for k in self.Accretions]
        evap_af = [np.nan] + [k/AFtoM3 for k in self.Evaporation]
        trib_af = [k/AFtoM3 for k in self.TribInflow]
        mainin_af = [k/AFtoM3 for k in self.MainInflow]
        mainin_af = mainin_af[0:-1]
        trib_af = trib_af[0:-1]
        accret_af = accret_af[0:-1]
        
        area_ac = [k/ACREtoM2 for k in self.AreaRecord]
        #area_ac = area_ac[1:]
        
        dat_ = [releasetemps_degC, releasetemps_degF, releasevol_af, 
                self.StorageTemps, storage_af, accret_af, evap_af, trib_af, 
                mainin_af, area_ac]
        df_ = pnd.DataFrame(data=dat_, columns=index)
        df_ = df_.T
        df_.columns =columns
        
        return(df_)