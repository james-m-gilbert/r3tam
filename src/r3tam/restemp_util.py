# -*- coding: utf-8 -*-
"""
utilities for helping set up reservoir model object

Created on Wed Jul  1 16:12:50 2020

@author: jgilbert
"""
import os, sys
#import openpyxl as xl
import pandas as pnd
import numpy as np
import datetime as dt


#import restemp as srt
#import r3tam.restemp as rt
from r3tam import restemp as rt

from r3tam.constants import *

def setupRes_ElevSpec(resObj, resInputs):
    """
    set up reservoir layers based on a an ascending list of layer top elevations
    
    :param resObj: reservoir object for which an elevation list has been added
    
    :returns resObj: returned reservoir object has areas and volumes interpolated
    based on area-elevation-capacity table
    
    """
    
#    import yaml
#
#    with open(yamlPath) as stream:
#        resInputs = yaml.load(stream, Loader=yaml.SafeLoader)
    
    lyrElevUnits = resInputs['LayerElevations']['Units']   
    lyrElevVals = resInputs['LayerElevations']['Values']

    elevList = []
    
    fp = lyrElevVals[0]
    if not os.path.isabs(fp):  #if the provided filepath is not absolute (fully back to the drive or root), then assume it is relative to teh project directory
        fp = os.path.join(resObj.ProjDir, fp)
    
    if os.path.exists(fp): # layer elevations specified in a file

        #print9os.path.splite
#        if os.path.splitext(fp)[1] in ['.xls', '.xlsm', '.xlsx', '.xlsb']:
#            wb = xl.open(fp, read_only=True, data_only=True)
#            sht = wb[lyrElevVals[1]]
#            rng1 = lyrElevVals[2] #.split(':') , rng2
#            for row in xl.utils.rows_from_range(rng1): #rowssht.iter_rows(rng1,rng2):
#                for cell in row:
#                    r,c = xl.utils.coordinate_to_tuple(cell)
#                    elevList.append(sht.cell(r,c).value)
        if os.path.splitext(fp)[1] in ['.txt', '.dat']:
            tmpdf = pnd.read_table(fp, sep='\s+')
            elevList = tmpdf.iloc[:,1].to_list()
                    

        if os.path.splitext(fp)[1] in ['.csv']:
            tmpdf = pnd.read_table(fp, sep=',')
            elevList = tmpdf.iloc[:,1].to_list()
                
        if resObj.CalcMetric and lyrElevUnits.upper() in ['FT','FEET']:
            # convert layer elevations to meters
            elevList = [_*FTtoM for _ in elevList]
            
        if not resObj.CalcMetric and lyrElevUnits.upper() in ['M','METERS','METER']:
            # conver layer elevations from meters to feet
            elevList = [_*MtoFT for _ in elevList]
    else:
        print("Error: Couldnt' find layer elevation file %s" %fp)
        return
        #elevList = lyrElevVals # assuming this is parsed as a list of values if not a file path 
    
    #get elevation area
    elevAreaUnits = resInputs['ElevAreaTable']['Units']
    elevAreaVals = resInputs['ElevAreaTable']['Values']
    

    fp = elevAreaVals[0]
    if not os.path.isabs(fp):  #if the provided filepath is not absolute (fully back to the drive or root), then assume it is relative to teh project directory
        fp = os.path.join(resObj.ProjDir, fp)
    
    if os.path.exists(fp):

        if os.path.splitext(fp)[1] in ['.dat','.txt']:
            storElevArea = getElevAreaTable(fp, elev_col=0, area_col=1)
            
            if resObj.CalcMetric and elevAreaUnits.upper() in ['FT','FEET']:
                # convert elevation to m and area to m2
                storElevArea['Elevation'] = [_*FTtoM for _ in storElevArea.Elevation]
                storElevArea['Area'] = [_*ACREtoM2 for _ in storElevArea.Area]
                
            if not resObj.CalcMetric and elevAreaUnits.upper() in ['M','METERS','METER']:
                # conver layer elevations from meters to feet and areas from m2 to acres
                storElevArea['Elevation'] = [_*MtoFT for _ in storElevArea.Elevation]
                storElevArea['Area'] = [_/ACREtoM2 for _ in storElevArea.Area]

    else:
        print("Error: Couldnt' find elevation-area file %s" %fp)
        return
        #storElevArea = pnd.DataFrame(elevAreaVals, columns=['Elevation_ft','Storage_AF','Area_ac']) # assuming this is parsed as a list of values if not a file path        
    
    # get elevation storage
    elevStorUnits = resInputs['ElevStorTable']['Units']
    elevStorVals = resInputs['ElevStorTable']['Values']
    
    fp = elevStorVals[0]
    if not os.path.isabs(fp):  #if the provided filepath is not absolute (fully back to the drive or root), then assume it is relative to teh project directory
        fp = os.path.join(resObj.ProjDir, fp)
    
    if os.path.exists(fp):

        if os.path.splitext(fp)[1] in ['.dat','.txt']:
            storElev = getElevStorTable(fp, elev_col=0, stor_col=1)
            
            if resObj.CalcMetric and elevStorUnits.upper() in ['FT','FEET']:
                # convert elevation to m and area to m2
                storElev['Elevation'] = [_*FTtoM for _ in storElev.Elevation]
                storElev['Storage'] = [_*AFtoM3 for _ in storElev.Storage]
                
            if not resObj.CalcMetric and elevAreaUnits.upper() in ['M','METERS','METER']:
                # conver layer elevations from meters to feet and areas from m2 to acres
                storElevArea['Elevation'] = [_*MtoFT for _ in storElevArea.Elevation]
                storElevArea['Area'] = [_/AFtoM3 for _ in storElevArea.Area]

    else:        
        print("Error: Couldnt' find elevation-storage file %s" %fp)
        return
        storElev = pnd.DataFrame(elevStorVals, columns=['Elevation','Storage']) # assuming this is parsed as a list of values if not a file path        
        
    resObj.ReservoirData.ElevArea = storElevArea
    resObj.ReservoirData.ElevStorage = storElev
    # now calculate areas and incremental volumes
    # assume first elevation is bottom of reservoir, bottom of first layer
    #botElev = elevList[0]
    prevcumulVol = 0.
    for i,ev in enumerate(elevList[1:]):
        # all elevations will be adjsuted to the correct units system above
        botElev = elevList[i]
        il = rt.lyr()
        il.Index = i
        topElev = ev
        b = topElev - botElev
        ctrelev = botElev + 0.5*b
        il.Thick = b
        il.MinElev = botElev
        il.MaxElev = topElev
        il.CtrElev = ctrelev
        il.GateHeight = b 
        #il.GateHeight_m = b*FTtoM #<--no longer necessary
        # now interpolate the area
        il.MaxArea = np.interp(topElev, storElevArea['Elevation'], storElevArea['Area'])
        if resObj.CalcMetric:
            il.MaxAreaUnits = 'm2'
        else:
            il.MaxAreaUnits = 'acre'

        cumulVolAtElev = np.interp(topElev,storElev['Elevation'],storElev['Storage'] )
        il.MaxVol = cumulVolAtElev - prevcumulVol
        
        if resObj.CalcMetric:
            il.MaxVolUnits = 'm3'
        else:
            il.MaxVolUnits = 'acre-ft'
        #il.MaxVol_m3 = il.MaxVol*AFtoM3
        prevcumulVol = cumulVolAtElev
        resObj.Layers[i] = il        
    
    
    return(elevList)
    
def setupOutlets(resObj, inputs):

    strctrs = inputs
    for ro in strctrs['RiverOutlets']:
        rout = rt.riv_out()
        rout.ID = ro
        rodata = strctrs['RiverOutlets'][ro]
        cap = float(rodata['Capacity'])
        capUnits = rodata['CapacityUnits']
        if not resObj.CalcMetric and capUnits.lower() in ['cms','m3/s','cubic meters per second','m3_s','cubic_meters_per_second']:
            cap = cap*CMStoCFS
        elif resObj.CalcMetric and capUnits.lower()  in ['cfs','ft3/s', 'cubic feet per second', 'ft3_s','cubic_feet_per_second']:
            cap = cap/CMStoCFS
        # elif capUnits.lower() in ['cfs','ft3/s', 'cubic feet per second']:
        #     cap = cap
        else:
            print("Unknown units for river outlet capacity - assuming CFS")
            if resObj.CalcMetric:
                cap = cap/CMStoCFS
                capUnits = 'CMS'
        rout.Capacity = cap
        
        [te_ft, be_ft, ce_ft] = evalElevations(resObj, rout, rodata, 'river outlet',ro)
        rout.MinElev = be_ft
        rout.MaxElev = te_ft
        rout.CtrElev = ce_ft
        rout.CapacityUnits = capUnits
        
        resObj.RiverOutlets[ro] = rout
        
    # set up leakage zones
    # TODO: do a check for existence/definition of leakage zones in yaml dictionary    
    for lk in strctrs['LeakageZones']:
        lkz = rt.leakage()
        lkz.Zone = lk
        
        lkdata = strctrs['LeakageZones'][lk]
        lkz.LeakageFactor = lkdata['LeakageFactor']
        lkz.OutletExclude = lkdata['OutletExclude']
        
        [te_ft, be_ft, ce_ft] = evalElevations(resObj, lkz, lkdata, 'leakage zone',lk)
        
        lkz.TopElevFt = te_ft
        lkz.BotElevFt = be_ft
        lkz.CtrElev = ce_ft
        
        resObj.Leakages[lk] = lkz

    print("\n\n-----------------------------------------------")
    print("   Reading in PENSTOCK configuration       ")        
    for pnsk in strctrs['Penstocks']:
        print(f"\t\tVariable:   {pnsk}")
        if pnsk.lower() == 'penstockelevation':
            resObj.PenstockElevation = strctrs['Penstocks'][pnsk]
            #continue
            
        if pnsk.lower() == 'penstockelevunits':
            resObj.PenstockElevUnits = strctrs['Penstocks'][pnsk]
            #continue
            
        if resObj.PenstockElevUnits.upper() in ['FT','FEET']:
            resObj.PenstockElevation_m = resObj.PenstockElevation*FTtoM
            #continue
        else:
            resObj.PenstockElevation_m = resObj.PenstockElevation
            #continue    
        
        if pnsk.lower() == 'penstockreleaselimit':
            
            resObj.PenstockLimit = strctrs['Penstocks'][pnsk]
            print(f"\t\t\tSetting penstock release limit to {resObj.PenstockLimit}")
            #continue
        
        if pnsk.lower() == 'penstockreleaselimit_units':
            resObj.PenstockLimit_Units = strctrs['Penstocks'][pnsk]
            
    print("------------------DONE---------------------------\n")
    # set up selective withdrawal levels
    # TODO: do a check for existence/definition of selective withdrawal levels
    #       in yaml dictionary    
    for swd in strctrs['Selective Withdrawal']:

        if type(swd)==str:
            if swd.lower() == 'checkminhead':
                resObj.CheckMinHead = strctrs['Selective Withdrawal'][swd]
                continue
            else:
                resObj.CheckMinHead = False
                continue
            
        swdl = rt.outlet()
        
        swddata = strctrs['Selective Withdrawal'][swd]
        swdl.Name = swddata['Name']
        swdl.NumGates = swddata['NumberOfGates']
        if 'BalanceFactor' in swddata:
            swdl.BalanceFactor = swddata['BalanceFactor']
            swdl.IsolationFactor = swddata['IsolationFactor']
        else:
            swdl.BalanceFactor = 1
            swdl.IsolationFactor = 1
        swdl.MinHead = swddata['MinHead']
        swdl.MinHdUnits = swddata['MinHeadUnits'].upper()
        swdl.A1GT = swddata['A1GT']
        swdl.B1GT = swddata['B1GT']
        swdl.G1GT = swddata['G1GT']
        swdl.PointSinks = swddata['NumPointSinks']
        [te_ft, be_ft, ce_ft] = evalElevations(resObj, swdl, swddata, 'selective withdrawal', swd)
        
        swdl.TopElevFt = te_ft
        swdl.BotElevFt = be_ft
        swdl.CtrElevFt = ce_ft
        swdl.GateHeight_m = (te_ft-be_ft)*FTtoM
        resObj.Outlets[swd] = swdl
        resObj.GateDict[swd] = 0

def getElevStorTable(fp, elev_col = 0, stor_col=1):
    '''
        function to retrieve elevation-storage paired data table
        assumes fp is path to a ascii text file with two,
        whitespace delimited columns, first row column headers
    '''        
    tmpdf = pnd.read_table(fp, sep='\s+')
    if elev_col==0:
        hdr = ['Elevation', 'Storage']
    else:
        hdr = ['Storage', 'Elevation']
    tmpdf.columns = hdr
    return(tmpdf)
    
def getElevAreaTable(fp, elev_col = 0, area_col=1):
    '''
        function to retrieve elevation-area paired data table
        assumes fp is path to a ascii text file with two,
        whitespace delimited columns, first row column headers
    '''        
    tmpdf = pnd.read_table(fp, sep='\s+')
    if elev_col==0:
        hdr = ['Elevation', 'Area']
    else:
        hdr = ['Area', 'Elevation']
    tmpdf.columns = hdr
    return(tmpdf)    
    
#def getStorageAreaTable(xlFP, xlSht, startRow, startCol, maxRow):
#    wb = xl.open(xlFP, read_only=True, data_only=True)
#    
#    shtIns = wb[xlSht]
#    
#    storElev = []
#    for r in range(startRow, maxRow+1):
#        rd = []
#        for c in range(startCol, startCol+3):
#            rd.append(shtIns.cell(r,c).value)
#        storElev.append(rd)
#    wb.close()
#    storElevDF = pnd.DataFrame(storElev, columns=['Elevation_ft','Storage_AF','Area_ac'])
#    return(storElevDF)
#
#def getStorageElevTable(xlFP, xlSht, startRow, startCol):
#    wb = xl.open(xlFP, read_only=True, data_only=True)
#    
#    shtIns = wb[xlSht]
#    
#    storElev = []
#    for r in range(startRow, shtIns.max_row+1):
#        rd = []
#        for c in range(startCol, startCol+2):
#            rd.append(shtIns.cell(r,c).value)
#        storElev.append(rd)
#    wb.close()
#    storElevDF = pnd.DataFrame(storElev, columns=['Elevation_ft','Storage_AF'])
#    return(storElevDF)
        
def evalElevations(resObj, outletObj, dictObj, outletType, idx):
    from sympy.parsing.sympy_parser import parse_expr
    # top elevation - check for expression to 
    rodata = dictObj
    te = ''  # top elevation
    be = ''  # bottom elevation
    ce = ''  #center elevation
    if 'TopElev' in rodata.keys():
        if any(op in str(rodata['TopElev']) for op in ['+','-','*','/']):
            te = float(parse_expr(rodata['TopElev']))
        else:
            te = rodata['TopElev']
        if rodata['TopElevUnits'].lower() not in ACCEPTABLE_LENGTH_UNITS:
            print("No units assigned to top elevation for %s outlet %s" %(outletType, idx))
            print(" --- Assuming units of FEET")
            if resObj.CalcMetric:
                te = te*FTtoM
        else:
            if not resObj.CalcMetric and rodata['TopElevUnits'].lower() in ['m', 'meter','meters']:
                te = te*MtoFT
                
            if resObj.CalcMetric and rodata['TopElevUnits'].lower() in ['ft','feet']:
                te = te*FTtoM

    
    if 'BottomElev' in rodata.keys():
        if any(op in str(rodata['BottomElev']) for op in ['+','-','*','/']):
            be = float(parse_expr(rodata['BottomElev']))
        else:
            be = rodata['BottomElev']
        if rodata['BottomElevUnits'].lower() not in ACCEPTABLE_LENGTH_UNITS:
            print("No units assigned to bottom elevation for %s outlet %s" %(outletType, idx))
            print(" --- Assuming units of FEET")
            if resObj.CalcMetric:
                be = be*FTtoM
        else:
            if not resObj.CalcMetric and rodata['BottomElevUnits'].lower() in ['m', 'meter','meters']:
                be = be*MtoFT

            if resObj.CalcMetric and rodata['BottomElevUnits'].lower() in ['ft','feet']:
                be = be*FTtoM
    
    if 'CenterElev' in rodata.keys():
        if any(op in str(rodata['CenterElev']) for op in ['+','-','*','/']):
            ce = float(parse_expr(rodata['CenterElev']))
        else:
            ce = rodata['CenterElev']
        if rodata['CtrElevUnits'].lower() not in ACCEPTABLE_LENGTH_UNITS:
            print("No units assigned to center elevation for %s outlet %s" %(outletType, idx))
            print(" --- Assuming units of FEET")
            if resObj.CalcMetric:
                ce = ce*FTtoM
        else:
            if not resObj.CalcMetric and rodata['CtrElevUnits'].lower() in ['m', 'meter','meters']:
                ce = ce*MtoFT
            if resObj.CalcMetric and rodata['CtrElevUnits'].lower() in ['ft','feet']:
                ce = ce*FTtoM
    else:
        ce = 0.5*(te + be)
    
    return([te, be, ce])


def setupSim(resObj, resInputs, projDir, simType='vert'):

        
    startDateStr = resInputs['Timing']['StartDate']
    endDateStr = resInputs['Timing']['EndDate']
    delt_days = resInputs['Timing']['TimeStepDays']
    
    if 'Debug' in resInputs['Simulation Mode']:
        print("Assigning Debug levels")
        if 'Release' in resInputs['Simulation Mode']['Debug']:
            resObj.Debug['Release'] = resInputs['Simulation Mode']['Debug']['Release']
        else:
            resObj.Debug['Release'] = int(resInputs['Simulation Mode']['Debug'])
            
        if 'General' in resInputs['Simulation Mode']['Debug']:
            resObj.Debug['General'] = resInputs['Simulation Mode']['Debug']['General']
        else:
            resObj.Debug['General'] = int(resInputs['Simulation Mode']['Debug'])
    
    
    resObj.SimDates = list(pnd.date_range(startDateStr, endDateStr, freq='d'))  #TODO: make this adjustable in case we want to do something different from daily (3-hourly or 2-daily, for exmaple)
    resObj.RunName = resInputs['RunName']
    resObj.ProjDir = projDir
    resObj.LogFilePath = os.path.join(projDir, f'log{dt.datetime.today().isoformat(timespec="minutes")}.log')
    
    if 'SimulationUnits' in resInputs:
        simunits = resInputs['SimulationUnits']
        if simunits.upper() in ['FT','ACFT','AFDEGF','DEGF','F']:
            resObj.CalcMetric = False
            print("\n-------------------------------------------")
            print("Model calculations will be done in units of:\n\tACRE-FEET, ACRES, and DEG F")
            print("-------------------------------------------")
        elif simunits.upper() in ['M3DEGC','M3','M','DEGC','C']:
            resObj.CalcMetric = True
            print("\n-------------------------------------------")
            print("Model calculations will be done in units of:\n\tCUBIC METERS, SQ METERS, and DEG C")
            print("-------------------------------------------")
        else:
            resObj.CalcMetric = True
            print("\n-------------------------------------------")
            print("No units specifications for calculationfound in config file")
            print(f"\tEx: SimulationUnits: 'm3degC' in {resObj.ConfigFP}")
            print("Model calculations will be done in units of:\n\tCUBIC METERS, SQ METERS, and DEG C")
            print("-------------------------------------------")
            
    
    if 'FacilityName' in resInputs:
        resObj.FacilityName = resInputs['FacilityName']
        
    locdata = resInputs['Location']
    if 'Description' in locdata:
        resObj.LocationDescription = locdata['Description']
    resObj.Latitude = locdata['Latitude']
    resObj.Longitude = locdata['Longitude']
    
    resObj.SimulationSpecs.StartDate_str = startDateStr
    resObj.SimulationSpecs.EndDate_str = endDateStr
    resObj.SimulationSpecs.DELT_DAY = delt_days
    resObj.TimeStepDate = resObj.SimDates[0]
    
    if simType.upper()=='VERT':
        resObj.Surface = resInputs['SurfaceMethod'].lower()
    
    resObj.SimulationSpecs.InitMethod = resInputs['Initial Condition']['Method']

    if str(resInputs['Simulation Mode']['Hindcast']).upper() in ['T', 'TRUE', 'Y','YES']:
        resObj.SimulationSpecs.IsHindcast = True
    else:
        resObj.SimulationSpecs.IsHindcast = False
    
    if simType.upper() in ['VERT','VERTICAL','STRAT','STRATIFIED']:    
        if str(resInputs['Simulation Mode']['SimGateOps']).upper() in ['T','TRUE','Y','YES']:
            resObj.SimulationSpecs.CalcGateOps = True
            resObj.SimulationSpecs.GateChgFreq = int(resInputs['Simulation Mode']['GateChgFreq'])
        else:
            resObj.SimulationSpecs.CalcGateOps = False
        
    if str(resInputs['Simulation Mode']['SimReleaseSchedule']).upper() in ['T','TRUE','Y','YES']:
        resObj.SimulationSpecs.CalcReleaseSched = True
    else:
        resObj.SimulationSpecs.CalcReleaseSched = False
        
    if str(resInputs['Simulation Mode']['CalcEvap']).upper() in ['T', 'TRUE','Y','YES']:
        resObj.SimulationSpecs.CalcEvaporation = True
    else:
        resObj.SimulationSpecs.CalcEvaporation = False
        
    if 'ReinitMonth' in resInputs['Simulation Mode']:
        resObj.SimulationSpecs.Reinitialize = True
        resObj.SimulationSpecs.ReinitMonth = resInputs['Simulation Mode']['ReinitMonth']
        if 'ReinitDay' in resInputs['Simulation Mode']:
            resObj.SimulationSpecs.ReinitDay = resInputs['Simulation Mode']['ReinitDay']
        else:
            resObj.SimulationSpecs.ReinitDay = 1
            
    # add temperature target data if present
    if 'Temperature_Target' in resInputs:
        ttarg_dat = resInputs['Temperature_Target']
        ttargLoc = ttarg_dat['Location']
        ttargDef = ttarg_dat['Default']
        if 'Tolerance' in ttarg_dat:
            ttargTol = ttarg_dat['Tolerance']
        else:
            ttargTol = {m:2 for m in range(1,13)}
        if 'Lower_Tolerance' in ttarg_dat:
            ttargLowTol = ttarg_dat['Lower_Tolerance']
        else:
            ttargLowTol = {m:0 for m in range(1, 13)}
        resObj.TempTargetData = {}
        resObj.TempTargetData[ttargLoc] = {'default': ttargDef, 'tol': ttargTol, 'tol_low': ttargLowTol}

    #if 

def setupCalib(resObj, resInputs):
    
    calib_info = resInputs['Calibration Mode']
    calib_method = calib_info['Method']
    calib_bounds_input = calib_info['Bounds']
    
    if '.' in calib_bounds_input: # assume it's a file then
        fp = os.path.join(resObj.ProjDir, 'inputs',calib_bounds_input)
        if os.path.exists(fp):
            tmp = getInputsDict(fp)
            calib_bounds_dict = tmp['Bounds']
        else:
            print("Error: couldn't find calibration bounds file %s" %fp)
            exit
    
    
    resObj.Calibration.Method = calib_method
    resObj.Calibration.MaxIter = int(calib_info['MaxIter'])
    resObj.Calibration.FATOL = float(calib_info['FATOL'])
    resObj.Calibration.XATOL = float(calib_info['XATOL'])
    resObj.Calibration.Display = bool(calib_info['DISP'])
    
    # set up the initial parameter values (assuming we're starting with whatever
    # value is provided in the original *.yaml file variables) and the boudns
    # on those variables - bounds are assumed to be a list as [lower, upper]
    
    # first, get teh variables we're changing from the bounds dictionary keys
    calib_vars = list(calib_bounds_dict.keys())
    initParamDict = {}
    boundsDict = {}
    
    for par in calib_vars:
        # assign parameter values to reservoir object
        #parsp = par.split(':')
        parsp= par.split('@')
        
        print(par, parsp)
        
        if parsp[0].upper() =='M': # main parameters
            print(par)
            #initParamDict[par] = resObj.__dict__[parsp[1]]
            #boundsDict[par] = calib_bounds_dict[par]
            resObj.Calibration.InitGuess[par] = resObj.__dict__[parsp[1]]
            resObj.Calibration.Bounds[par] = calib_bounds_dict[par]

#        elif parsp[0].upper()=='LKG':
#            sha.Leakages[parsp[1]].LeakageFactor= val
            
        elif parsp[0].upper()=='SWD':   #selective withdrawal variables
            #initParamDict[par] = resObj.Outlets[parsp[1]].__dict__[parsp[2]]
            #boundsDict[par] = calib_bounds_dict[par]
            gateID = int(parsp[1])
            resObj.Calibration.InitGuess[par] = resObj.Outlets[gateID].__dict__[parsp[2]]
            resObj.Calibration.Bounds[par] = calib_bounds_dict[par]
        else:
            pass
   
    #return([initParamDict, boundsDict])
    
def initVolTemp(resObj, initVol, initTempProf, profileTempUnits='degF', profileElevUnits = 'ft'):
    tmpVol = initVol
    
    # # adjust the input profile to the correct units here
    # if resObj.CalcMetric and profileTempUnits.upper() in ['DEG_F','DEGF','F','FAHRENHEIT']:
    #     if profileElevUnits.upper() in ['FT','FEET']:
    #         initTempProfElevs = [_*FTtoM for _ in initTempProf.index]
    #     else:
    #         initTempProfElevs = [_ for _ in initTempProf.index]
            
    #     initTemps = [DEGF_to_DEGC(_) for _ in initTempProf.iloc[:,0]]

    # if not resObj.CalcMetric and profileTempUnits.upper() in ['DEG_C','DEGC','C', 'CELSIUS']:
    #     if profileElevUnits.upper() in ['M','METER', 'METERS']:
    #         initTempProfElevs = [_*MtoFT for _ in initTempProf.index]
    #     else:
    #         initTempProfElevs = [_ for _ in initTempProf.index]
            
    #     initTemps = [DEGC_to_DEGF(_) for _ in initTempProf.iloc[:,0]]
    initTempProfElevs = [_ for _ in initTempProf.index]    
    initTemps = [_ for _ in initTempProf.iloc[:,0]]
    
    for l,v in resObj.Layers.items():
        maxv = v.MaxVol
        if tmpVol - maxv >=0:
            v.Vol = maxv
            tmpVol -= maxv
        else:
            v.Vol = tmpVol
            tmpVol -= tmpVol
        if v.Vol == 0.:
            v.Temp = np.NaN
        else:
            low = np.NaN
            # if resObj.CalcMetric and profileTempUnits.upper() in ['DEGF','DEG_F','F']: #tempUnits =='degC':
            #     up = np.nanmax(initTempProf)*1.8+32.
            #     low = np.nanmin(initTempProf)*1.8 +32.
            # else:
            #     up = np.nanmax(initTempProf)
            #     low = np.nanmin(initTempProf)
            
            up = np.nanmax(initTempProf)
            low = np.nanmin(initTempProf)
            
            initTPelevs = initTempProfElevs #initTempProf.index
            initTPtemps = initTemps #initTempProf.values
            
            if min(initTPelevs) < v.CtrElev:
                lowElev = max(i for i in initTPelevs if i < v.CtrElev)
            else:
                lowElev = min(initTPelevs)
            if max(initTPelevs) > v.CtrElev:
                upElev = min(i for i in initTPelevs if i > v.CtrElev)
            else:
                upElev = max(initTPelevs) #sha.Layers[sha.nLyrs-1].CtrElev
            
            upTemp = initTemps[initTempProfElevs.index(upElev)]
            lowTemp = initTemps[initTempProfElevs.index(lowElev)]
            # if tempUnits=='degC':
            #     upTemp = initTempProf.loc[upElev]*1.8+32.
            #     lowTemp= initTempProf.loc[lowElev]*1.8+32.
            # else:
            #     upTemp = initTempProf.loc[upElev]
            #     lowTemp = initTempProf.loc[lowElev]
            
            if type(upTemp) !=float:
                if len(upTemp)>1: # in case there are redundant readings at same elevation
                    upTemp = upTemp[0] #.iloc[0].iloc[0]
            # else:
            #     upTemp = upTemp.iloc[0]
            
            if type(lowTemp) != float:
                if len(lowTemp)>1: # in case there are redundant readings at same elevation
                    lowTemp = lowTemp[0] #.iloc[0].iloc[0]
            # else:
            #     lowTemp = lowTemp.iloc[0]
                


            # if type(upTemp)==type(pnd.Series):
            #     print("here!")
            #     upTemp = upTemp.values[0]
                        
            
            #print(upElev, upTemp)
            #print("==============")
            #print(lowElev, lowTemp)
            
            
            if np.isnan(upTemp) & np.isnan(lowTemp):
                # if tempUnits=='degC':
                #     upTemp = np.nanmax(initTPtemps)*1.8+32.
                # else:
                #     upTemp = np.nanmax(initTPtemps)
                upTemp = np.nanmax(initTPtemps)
                lowTemp = upTemp
                
            elif np.isnan(upTemp) & ~np.isnan(lowTemp):
                # if tempUnits=='degC':
                #     upTemp = np.nanmax(initTPtemps)*1.8+32.
                # else:
                #     upTemp = np.nanmax(initTPtemps)
                upTemp = np.nanmax(initTPtemps)
                
            elif ~np.isnan(upTemp) & np.isnan(lowTemp):
                # if tempUnits=='degC':
                #     lowTemp = np.nanmin(initTPtemps)*1.8+32.
                # else:
                #     lowTemp = np.nanmin(initTPtemps)
                lowTemp = np.nanmin(initTPtemps)
                
            else:
                pass
    
            if upElev==lowElev:
                v.Temp = upTemp
            else:
                v.Temp = lowTemp + (v.CtrElev - lowElev)*(upTemp - lowTemp)/(upElev - lowElev)

        resObj.Layers[l].updTotE()

def setParams(resObj, inputs):
         
    params = inputs['Parameters']
    resObj.C1 = params['C1']
    resObj.C2 = params['C2']
    resObj.C3 = params['C3']
    resObj.C4 = params['C4']
    resObj.C5 = params['C5']
    resObj.CriticalDepth = params['CriticalDepth']
    resObj.CriticalDepthUnits = params['CriticalDepthUnits']
    resObj.MeanAlbedo = params['MeanAlbedo']
    resObj.ExtinctionCoeff = params['ExtinctionCoeff']
    
    if 'DiffuseDist' in params:
        dd = params['DiffuseDist']
        dd_units = params['DiffuseDistUnits']
    else:
        # default
        dd = 30.
        dd_units = 'ft'
        
    if resObj.CalcMetric and dd_units.upper() in ['FT','FEET']:
        dd = dd*FTtoM
        dd_units = 'm'
        
    if not resObj.CalcMetric and dd_units.upper() in ['M','METERS','METER']:
        dd = dd*MtoFT
        dd_units = 'ft'
    
    resObj.DIFFUSE_DIST = dd
    resObj.DIFFUSE_DIST_UNITS = dd_units 
    
    if 'WindMethod' in params:
        resObj.WindMethod = params['WindMethod'].upper()
    else:
        resObj.WindMethod = 'W2'

    if 'WindAlpha' in params:
        resObj.WindAlpha = params['WindAlpha']
    else:
        resObj.WindAlpha = 9.2 #default, based on CE-QUAL-W2 User's Manual Part 2, p84 (June 2021)
        
    if 'WindBeta' in params:
        resObj.WindBeta = params['WindBeta']
    else:
        resObj.WindBeta = 0.46 #default, based on CE-QUAL-W2 User's Manual Part 2, p84 (June 2021)

    if 'WindGamma' in params:
        resObj.WindGamma = params['WindGamma']
    else:
        resObj.WindGamma = 2. #default, based on CE-QUAL-W2 User's Manual Part 2, p84 (June 2021)
                
    if 'WindHeight_m' in params:
        resObj.WindHeight_m = params['WindHeight_m']
    else:
        # assume wind height is already at or adjusted to 2 m 
        resObj.WindHeight_m = 2.
                
    if 'SeasonalAlbedo' in params:
        resObj.SeasonalAlbedo = params['SeasonalAlbedo']
        if 'AlbedoAmplitude' in params:
            resObj.AlbedoAmplitude = params['AlbedoAmplitude']
        
    
    if 'ElevDependentCD' in params:
        tmp = params['ElevDependentCD']
        resObj.ElevDependentCD_k = float(tmp['k'])
        resObj.ElevDependentCD_A = float(tmp['A'])
        resObj.ElevDependentCD_x0 = float(tmp['x0'])
    
    if 'SeasonalRad' in params:
        resObj.SeasonalRad = float(params['SeasonalRad'])
    else:
        resObj.SeasonalRad = 1.
        
    if 'BlendingMethod' in params:
        if params['BlendingMethod'] ==1:
            resObj.LP_Opt_Blending = True
            print("\n\t\t---Setting blending option to LP Optimization\n")
        elif params['BlendingMethod'] ==0:
            resObj.LP_Opt_Blending = False
        else:
            resObj.LP_Opt_Blending = False
    
def setLongParams(resObj, inputs):
    params= inputs['Parameters']
    resObj.K = params['K']
    resObj.TTconst = params['TT_const']
    #resObj.CriticalDepth = params['CriticalDepth']
    #resObj.CriticalDepthUnits = params['CriticalDepthUnits']
    if 'C1' in params:
        resObj.C1 = params['C1']
    else:
        resObj.C1 = 1.
        
    if 'C2' in params:
        resObj.C2 = params['C2']
    if 'C3' in params:
        resObj.C3 = params['C3'] #evap coeff
    else:
        resObj.C3 = 1.
        
    # if 'SeasonalRad' in params:
    #     resObj.SeasonalRad = float(params['SeasonalRad'])
    # else:
    #     resObj.SeasonalRad = 1.
        
    resObj.MeanAlbedo = params['MeanAlbedo']
    
    if 'Wind_Func_Params' in params:
        wfp = params['Wind_Func_Params']
        for k,v in wfp.items():
            resObj.Wind_Func_Params[k.lower()] = v
            
    #resObj.ExtinctionCoeff = params['ExtinctionCoeff']
    if 'SeasonalAlbedo' in params:
        resObj.SeasonalAlbedo = params['SeasonalAlbedo']
        if 'AlbedoAmplitude' in params:
            resObj.AlbedoAmplitude = params['AlbedoAmplitude']
            
def getInputsDict(yamlPath):
    import yaml
    
    with open(yamlPath) as stream:
        inputs = yaml.load(stream, Loader=yaml.SafeLoader)    
        
    return(inputs)
    
def getTimeSeries(resmod, inputs):
    
    tsInfo = inputs['TimeSeriesInputs']
    thisdelt = resmod.SimulationSpecs.DELT_SEC # number of seconds per time step

    
    if 'Inflow' in  tsInfo:
        fdata = tsInfo['Inflow']
        fp = fdata['Path']
        
        if not os.path.isabs(fp):  # if not abs filepath, assume relative to project directory
            fp = os.path.join(resmod.ProjDir, fp)
            
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Couldnt' find inflow file {fp}")
            
        hdrDict = fdata['Header']
        
        df = pnd.read_table(fp, sep=',', index_col=0)
        try:
            df.index = pnd.to_datetime(df.index, format="%Y-%m-%d")
        except:
            try:
                df.index = pnd.to_datetime(df.index, format="%m/%d/%Y") # parse_dates=True)
            except:
                raise BaseException(f"Couldn't parse the date format in inflow file {fp}\n\t--Should be either 'Y-m-d' or 'm/d/Y'")
        
        try:
            inflowunits  = hdrDict['inflow'].split('_')[-1]
        except:
            raise BaseException("Inflow column name not formatted with units..cannot continue")
        
                
        if resmod.CalcMetric and inflowunits.upper() in ['AF','ACRE-FEET','AC-FT','AC_FT','ACRE_FEET','ACFT']:
            # convert to cubic meters per day
            df['InflowTotal'] = [_*AFDtoCMS*thisdelt for _ in df[hdrDict['inflow']]]
            hdrDict['inflow_total_units'] = 'cubic_meters_per_day'
                             
        elif resmod.CalcMetric and inflowunits.upper() in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
            # convert to cubic meters per day
            df['InflowTotal'] = [_*CFStoCMS*thisdelt for _ in df[hdrDict['inflow']]]
            hdrDict['inflow_total_units'] = 'cubic_meters_per_day'   
            
        elif not resmod.CalcMetric and inflowunits.upper() in ['CMS','M3_S', 'CUBIC_METERS_PER_SECOND']:
            df['InflowTotal'] = [_*CMStoAFD for _ in df[hdrDict['inflow']]] 
            hdrDict['inflow_total_units'] = 'acre_feet_per_day'
                            
        elif not resmod.CalcMetric and inflowunits.upper() in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
            df['InflowTotal'] = [_*CFStoAFD for _ in df[hdrDict['inflow']]]
            hdrDict['inflow_total_units'] = 'acre_feet_per_day'
        
        else:
            df['InflowTotal'] = df[hdrDict['inflow']]
            hdrDict['inflow_total_units'] = inflowunits
            
        hdrDict['inflow_final'] = 'InflowTotal' 
   
        try:
            inflowTunits  = hdrDict['inflowTemp'].split('_')[-1]
        except:
            raise BaseException("Inflow column name not formatted with units..cannot continue")
        
                
        if resmod.CalcMetric and inflowTunits.upper() in ['DEG_F', 'DEGF', 'F']:
            # convert to deg C
            df['InflowTemp_degC'] = [DEGF_to_DEGC(_) for _ in df[hdrDict['inflowTemp']]]
            hdrDict['inflowTemp_final'] = 'InflowTemp_degC'   
                                                 
        elif not resmod.CalcMetric and inflowTunits.upper() in ['DEG_C','DEGC','C']:
            df['InflowTemp_degF'] = [DEGC_to_DEGF(_) for _ in df[hdrDict['inflowTemp']]]   
            hdrDict['inflowTemp_final'] = 'InflowTemp_degF'
        else:
            if hdrDict['inflowTemp']=='' or hdrDict['inflowTemp']==None:
                df[f'InflowTemp_{inflowTunits}'] = [np.nan]*len(df)
            else:
                df[f'InflowTemp_{inflowTunits}'] = df[hdrDict['inflowTemp']]          
            hdrDict['inflowTemp_final'] = f'InflowTemp_{inflowTunits}'
        
        
        if 'tribInflow' in hdrDict:
            tribInflowUnits = hdrDict['tribInflow'].split('_')[-1].upper()
            if resmod.CalcMetric and tribInflowUnits in ['AF','ACRE-FEET','AC-FT','ACFT']:
                df['tribInflow_final'] = [_*AFDtoCMS*thisdelt for _ in df[hdrDict['tribInflow']]]
                hdrDict['tribInflow_final_units'] = 'cubic_meters_per_day'
                hdrDict['tribInflow_final'] = 'tribInflow_final'
            elif resmod.CalcMetric and tribInflowUnits in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
                df['tribInflow_final'] = [_*CFStoCMS*thisdelt for _ in df[hdrDict['tribInflow']]]
                hdrDict['tribInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
                hdrDict['tribInflow_final'] = 'tribInflow_final'
            elif resmod.CalcMetric and tribInflowUnits in ['M3','M3_D','CUBIC_METERS_PER_DAY']:
                df['tribInflow_final'] = [_ for _ in df[hdrDict['tribInflow']]]
                hdrDict['tribInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
                hdrDict['tribInflow_final'] = 'tribInflow_final'
            else:
                df['tribInflow_final'] = [_ for _ in df[hdrDict['tribInflow']]]
                hdrDict['tribInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
                
        if 'tribInflowTemp' in hdrDict:
            tribInflowTUnits = hdrDict['tribInflowTemp'].split('_')[-1].upper()
            if resmod.CalcMetric and tribInflowTUnits.upper() in ['DEG_F', 'DEGF', 'F']:
                # convert to deg C
                df['tribInflowTemp_degC'] = [DEGF_to_DEGC(_) for _ in df[hdrDict['tribInflowTemp']]]
                hdrDict['tribInflowTemp_final'] = 'tribInflowTemp_degC'   
                                                     
            elif not resmod.CalcMetric and inflowTunits.upper() in ['DEG_C','DEGC','C']:
                df['tribInflowTemp_degF'] = [DEGC_to_DEGF(_) for _ in df[hdrDict['tribInflowTemp']]]   
                hdrDict['tribInflowTemp_final'] = 'tribInflowTemp_degF'
            else:
                df[f'tribInflowTemp_{inflowTunits}'] = df[hdrDict['tribInflowTemp']]          
                hdrDict['tribInflowTemp_final'] = f'InflowTemp_{inflowTunits}'
        
        resmod.Inflow.DataFrame= df
        resmod.Inflow.ColumnMap = hdrDict       
    
        if 'extraInflow' in hdrDict:
            extra_inflow_units = hdrDict['extraInflow'].split('_')[-1].upper()
            if resmod.CalcMetric and extra_inflow_units in ['AF','ACRE-FEET','AC-FT','ACFT']:
                # convert to cubic meters per day
                hdrDict['extraInflow_final_units'] = 'cubic_meters_per_day'
                df['extraInflow_final'] = [_*AFDtoCMS*thisdelt for _ in df[hdrDict['extraInflow']]]
                hdrDict['extraInflow_final'] = 'extraInflow_final'
            elif resmod.CalcMetric and extra_inflow_units in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
                # convert to cubic meters per day
                hdrDict['extraInflow_final_units'] = 'cubic_meters_per_day'
                df['extraInflow_final'] = [_*CFStoCMS*thisdelt for _ in df[hdrDict['extraInflow']]]
                hdrDict['extraInflow_final'] = 'extraInflow_final'
            elif resmod.CalcMetric and extra_inflow_units in ['M3/S','M3_S','CUBIC_METERS_PER_SECOND','CMS']:
                df['extraInflow_final'] = [_*thisdelt for _ in df[hdrDict['extraInflow']]]
                hdrDict['extraInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
                hdrDict['extraInflow_final'] = 'extraInflow_final'
            elif resmod.CalcMetric and extra_inflow_units in ['M3','M3_D','CUBIC_METERS_PER_DAY']:
                df['extraInflow_final'] = [_ for _ in df[hdrDict['extraInflow']]]
                hdrDict['extraInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
                hdrDict['extraInflow_final'] = 'extraInflow_final'
            else:
                df['tribInflow_final'] = [_ for _ in df[hdrDict['tribInflow']]]
                hdrDict['tribInflow_final_units'] = 'CUBIC_METERS_PER_DAY'
        else:
            df['extraInflow_final'] = [0 for _ in df.index]
            hdrDict['extraInflow_final_units'] = 'cubic_meters_per_day'
            hdrDict['extraInflow_final'] = 'extraInflow_final'
        
    if 'Outflow' in  tsInfo:
        fdata = tsInfo['Outflow']
        fp = fdata['Path']
        
        if not os.path.isabs(fp):  # if not abs filepath, assume relative to project directory
            fp = os.path.join(resmod.ProjDir, fp)
        
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Couldnt' find outflow file {fp}")
        
        hdrDict = fdata['Header']
        
        df = pnd.read_table(fp, sep=',', index_col=0)
        try:
            df.index = pnd.to_datetime(df.index, format="%Y-%m-%d")
        except:
            try:
                df.index = pnd.to_datetime(df.index, format="%m/%d/%Y") # parse_dates=True)
            except:
                raise BaseException(f"Couldn't parse the date format in outflow file {fp}\n\t--Should be either 'Y-m-d' or 'm/d/Y'")
        
        
        try:
            outflowunits  = hdrDict['outflow'].split('_')[-1]
        except:
            raise BaseException("Outflow column name not formatted with units..cannot continue")
        
        if resmod.CalcMetric and outflowunits.upper() in ['AF','ACRE-FEET','AC-FT','AC_FT','ACRE_FEET','ACFT']:
            # convert to cubic meters per day
            df['OutflowTotal'] = [_*AFDtoCMS*thisdelt for _ in df[hdrDict['outflow']]]
            hdrDict['outflow_final_units'] = 'cubic_meters_per_day' 
            
        elif resmod.CalcMetric and outflowunits.upper() in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
            # convert to cubic meters per day
            df['OutflowTotal'] = [_*CFStoCMS*thisdelt for _ in df[hdrDict['outflow']]]
            hdrDict['outflow_final_units'] = 'cubic_meters_per_day'
            
        elif not resmod.CalcMetric and inflowunits.upper() in ['CMS','M3_S', 'CUBIC_METERS_PER_SECOND']:
            df['OutflowTotal'] = [_*CMStoAFD for _ in df[hdrDict['inflow']]] 
            hdrDict['outflow_final_units'] = 'acre_feet_per_day'
                            
        elif not resmod.CalcMetric and inflowunits.upper() in ['CFS','CUBIC_FEET_PER_SECOND','FT3_S']:
            df['OutflowTotal'] = [_*CFStoAFD for _ in df[hdrDict['inflow']]]
            hdrDict['outflow_final_units'] = 'acre_feet_per_day'
        else:
            df['OutflowTotal'] = df[hdrDict['outflow']]
            hdrDict['outflow_final_units'] = outflowunits
            
        if 'rivOutFlow' in hdrDict:
            if hdrDict['rivOutFlow'] == '' or hdrDict['rivOutFlow']!=None:
                df['RiverOutlet_Flow'] = [0 for _ in df[hdrDict['outflow']]]
                hdrDict['rivOutFlow'] = 'RiverOutlet_Flow'
                
            else:
                df['RiverOutlet_Flow'] = df[hdrDict['rivOutFlow']]
        else:
            df['RiverOutlet_Flow'] = [0 for _ in df[hdrDict['outflow']]]
            hdrDict['rivOutFlow'] = 'RiverOutlet_Flow'
            
        hdrDict['rivoutlet_final_units'] = outflowunits
        
        if 'storage' in hdrDict:
            storage_units = hdrDict['storage'].split('_')[-1].upper()
            if storage_units in  ['AF','ACRE-FEET','AC-FT','AC_FT','ACRE_FEET','ACFT']:
                df['Observed_Storage'] = [_*AFtoM3 for _ in df[hdrDict['storage']]]
                hdrDict['storage_obs'] = 'Observed_Storage'
            else:
                df['Observed_Storage'] = [_ for _ in df[hdrDict['storage']]]
                hdrDict['storage_obs'] = 'Observed_Storage'
                
        
        
        if 'extraOutflow' in hdrDict: # in cases where there are specified withdrawals from reservoir (as with CalSim)
            extraOutflow_units = hdrDict['extraOutflow'].split("_")[-1].upper()
            if extraOutflow_units not in ['AF','ACRE-FEET','AC-FT','AC_FT','ACRE_FEET','ACFT']:
                # assuming extra outflow is provided in acre-feet, if not, throw
                # an error
                raise BaseException("Extra Outflow not defined in units of acre-feet. Try again")
            
            df['DirectDiversions'] = [_ for _ in df[hdrDict['extraOutflow']]]
            hdrDict['extraOutflow'] = 'DirectDiversions'
        else:
            df['DirectDiversions'] = [0 for _ in df.index]
            hdrDict['extraOutflow'] = 'DirectDiversions'
        
        print("assigning outflow data")
        hdrDict['outflow_final'] = 'OutflowTotal'
        resmod.Outflow.DataFrame= df
        resmod.Outflow.ColumnMap = hdrDict    
    
    if 'Meteorology' in  tsInfo:
        fdata = tsInfo['Meteorology']
        fp = fdata['Path']
        
        if not os.path.isabs(fp):  # if not abs filepath, assume relative to project directory
            fp = os.path.join(resmod.ProjDir, fp)
        
        if not os.path.exists(fp):
            raise FileNotFoundError(f"Couldnt' find meteorology input file {fp}")
        
        hdrDict = fdata['Header']
        
        df = pnd.read_table(fp, sep=',', index_col=0)        
        try:
            df.index = pnd.to_datetime(df.index, format="%Y-%m-%d")
        except:
            try:
                df.index = pnd.to_datetime(df.index, format="%m/%d/%Y") # parse_dates=True)
            except:
                raise BaseException(f"Couldn't parse the date format in meteorology file {fp}\n" + \
                                    "\t--Should be either 'Y-m-d' or 'm/d/Y'")
    
        
        # convert air temp to celsius if calcs are in metric
        try:
            airtempunits  = hdrDict['airTemp'].split('_')[-1]
        except:
            raise BaseException("Outflow column name not formatted with units..cannot continue")
        
        hdrDict['airTemp_final_units'] = 'DEG_F' if airtempunits.upper() in ['DEGF','DEG_F','F','FAHRENHEIT'] else 'DEG_C'
        if resmod.CalcMetric and airtempunits.upper() in ['DEGF','DEG_F','F','FAHRENHEIT']:
            df['AirTemp'] = [DEGF_to_DEGC(_) for _ in df[hdrDict['airTemp']]]
            hdrDict['airTemp_final'] = 'AirTemp'
            hdrDict['airTemp_final_units'] = 'DEG_C'
        
        elif not resmod.CalcMetric and airtempunits.upper() in ['DEGC','DEG_C','C', 'CELSIUS']:
            df['AirTemp'] = [DEGC_to_DEGF(_) for _ in df[hdrDict['airTemp']]]
            hdrDict['airTemp_final_units'] = 'DEG_F'
            hdrDict['airTemp_final'] = 'AirTemp'
        
        else:
            hdrDict['airTemp_final'] = 'AirTemp'
            df['AirTemp'] = df[hdrDict['airTemp']]
        
        if 'evap' in hdrDict:
            try:
                evapunits = hdrDict['evap'].split('_')[-1]
            except:
                raise BaseException(f'Could not determine evaporation units from key {hdrDict["evap"]}')
                
            if evapunits.upper() in ['AF', 'AC-FT','ACREFEET']:
                newcol = '_'.join(hdrDict['evap'].split('_')[0:-1]+['CFS'])
                df[newcol] = [_/CFStoAFD for _ in df[hdrDict['evap']]]
                hdrDict['evap'] = newcol 
        
        resmod.Met.DataFrame= df
        resmod.Met.ColumnMap = hdrDict 

    # finally, check if the windspeed needs adjusting to equivalent 2-m height
    if resmod.WindHeight_m != 2.:
        print(f"\n*Adjusting windspeed to equivalent 2-m height from {resmod.WindHeight_m} m\n")
        # adjust windspeed to it's 2-m height equivalent
        # based on CE-QUAL-W2's windspeed height adjustment method 
        # Source: CE-QUAL-W2 User's Manual Part 2, p84-85 (June 2021)
        z0a = 0.001 # roughness height m - for winds < 5 mph or < 2.3 m/s
        z0b = 0.005 # roughtness height, m - for winds > 5mph or >2.3 m/s
        
        # add columns for roughness height and correction factor
        windcol = resmod.Met.ColumnMap['wind']
        resmod.Met.DataFrame['Wind_Rough_Height_m'] = [z0a if i < 2.3 else z0b for i in resmod.Met.DataFrame[windcol]]
        resmod.Met.DataFrame['Wind_Height_Adj_Factor'] = [np.log(2/i)/np.log(resmod.WindHeight_m/i) for i in resmod.Met.DataFrame['Wind_Rough_Height_m']]
        resmod.Met.DataFrame['Wind_2m'] = resmod.Met.DataFrame.apply(lambda row: row[windcol]*row['Wind_Height_Adj_Factor'], axis=1)
    else:
        resmod.Met.DataFrame['Wind_2m'] = resmod.Met.DataFrame[resmod.Met.ColumnMap['wind']]
        
    resmod.Met.ColumnMap['wind2m'] = 'Wind_2m'

def getObservations(resObj, inputs):
    thisdelt = resObj.SimulationSpecs.DELT_SEC # number of seconds per time step
    obsInfo = inputs['Observations']

    if 'Outflow' in obsInfo:
        outflowInfo = obsInfo['Outflow']
        if 'File' in outflowInfo:
            fdata = outflowInfo['File']
            fp = fdata['Path']
            
            if not os.path.isabs(fp):  # if not abs filepath, assume relative to project directory
                fp = os.path.join(resObj.ProjDir, fp)
            
            df = pnd.read_table(fp, sep=',', index_col=0)
            try:
                df.index = pnd.to_datetime(df.index, format="%Y-%m-%d") 
            except:
                try:
                    df.index = pnd.to_datetime(df.index, format="%m/%d/%Y")
                except:
                    raise BaseException(f"Couldn't parse dates in outflow file {fp}")
            df.index = pnd.to_datetime(df.index, format="%m/%d/%Y") # parse_dates=True)
            if 'Header' in fdata:
                hdrDict = fdata['Header']
                # infer units of data
                try:
                    flowunits = hdrDict['flow'].split('_')[-1]
                except:
                    raise BaseException("Couldn't determine units on flow observations")
                
                try:
                    outtempunits = hdrDict['temperature'].split('_')[-1]
                except:
                    raise BaseException("Couldn't determine units on outflow temperature observations")
           
                try:
                    storunits = hdrDict['storage'].split('_')[-1]
                except:
                    raise BaseException("Couldn't determine units on storage observations")
                    
            else:
                print("No header information provided for outflows")
                print("Assuming first column is flow, and second column is temperature")
                hdrDict = {'flow': df.columns[0], 'temperature': df.columns[1]}
                
            # make adjustements for outflow volume units
            if resObj.CalcMetric and flowunits.upper() in ['AF','ACRE-FEET','ACFT']:
                df['ObsOutflow'] = df[hdrDict['flow']]*AFDtoCMS*thisdelt
                
            elif resObj.CalcMetric and flowunits.upper() in ['CFS','FT3/S', 
                                                             'CUBIC FEET PER SECOND', 
                                                             'FT3_S',
                                                             'CUBIC_FEET_PER_SECOND']:
                df['ObsOutflow'] = df[hdrDict['flow']]*CFStoCMS*thisdelt
                
            elif not resObj.CalcMetric and flowunits.upper() in ['CMS','M3/S','CUBIC METERS PER SECOND',
                                                                 'M3_S',
                                                                 'CUBIC_METERS_PER_SECOND']:
                df['ObsOutflow'] = df[hdrDict['flow']]*CMStoAFD*resmod.SimulationSpecs.DELT_DAY 
                
            else:
                df['ObsOutflow'] = df[hdrDict['flow']]
            
            hdrDict['obs_flow_final'] = ['ObsOutflow']


            # make adjustements for outflow temperautre units
            if resObj.CalcMetric and outtempunits.upper() in ['DEG_F','DEGF','F','FAHRENHEIT']:
                df['ObsOutTemp'] = [DEGF_to_DEGC(_) for _ in df[hdrDict['temperature']]]
            elif not resObj.CalcMetric and outtempunits.upper() in ['DEG_C','DEGC','C','CELSIUS']:
                df['ObsOutTemp'] = [DEGC_to_DEGF(_) for _ in df[hdrDict['temperature']]]
            else:
                if hdrDict['temperature']=='' or hdrDict['temperature']==None:
                    df['ObsOutTemp'] = [np.nan]*len(df)
                else:
                    df['ObsOutTemp'] = df[hdrDict['temperature']]
            hdrDict['obs_outtemp_final'] = ['ObsOutTemp']
            
            # finally make adjustements for storage - used in assigning initial condition
            if resObj.CalcMetric and storunits.upper() in ['AF','ACRE-FEET','ACFT']:
                df['ObsStorage'] = [_*AFtoM3 for _ in df[hdrDict['storage']]]
            elif resObj.CalcMetric and storunits.upper() in ['M3','CUBIC_METERS', 'CUBIC METERS']:
                df['ObsStorage'] = [_/AFtoM3 for _ in df[hdrDict['storage']]]
            else:
                df['ObsStorage'] = df[hdrDict['storage']]
            hdrDict['obs_storage_final'] = 'ObsStorage'
                
            resObj.Observations.OutflowDF = df
            resObj.Observations.OutflowColMap = hdrDict
        else:
            print("ERROR::Not currently equiped to deal with anything other than a file here")
            print("-----exiting-----")
            return
    else:
        print("Outflow observations not specified correctly in input source file")
        print("No outflow data added to model")
    
    if 'Profiles' in obsInfo:
        profInfo = obsInfo['Profiles']
        if 'File' in profInfo:
            fdata = profInfo['File']
            fp = fdata['Path']
            
            if not os.path.isabs(fp):  # if not abs filepath, assume relative to project directory
                fp = os.path.join(resObj.ProjDir, fp)
            
            
            if 'Header' in fdata:
                hdrDict = fdata['Header']
                if 'index_cols' in hdrDict:
                    index_col = hdrDict['index_cols']
                    hdrDict = {k:v for k,v in hdrDict.items() if k != 'index_cols'}
                    #df = pnd.read_table(fp, sep=',', index_col=index_col, parse_dates=True)
                    df = pnd.read_table(fp, sep=',')
                    try:
                        df.iloc[:,index_col[0]] = pnd.to_datetime(df.iloc[:,index_col[0]], format="%Y-%m-%d")
                    except:
                        try:
                            df.iloc[:,index_col[0]] = pnd.to_datetime(df.iloc[:,index_col[0]], format="%m/%d/%Y")
                        except:
                            raise BaseException(f"Couldn't parse dates in the profile data file {fp}")
                            
                    colnames = df.columns[index_col]
                    print(df.head())
                    df.set_index([x for x in colnames], inplace=True)
                    #df.index = pnd.to_datetime(df.index, format="%m/%d/%Y") # parse_dates=True)
                else:
                    df = pnd.read_table(fp, sep=',', index_col=0, parse_dates=True)
                    
                # adjust for units change
                elevunits = hdrDict['elevation'].split('_')[-1]
                hdrDict['elevations_final'] = hdrDict['elevation']
                newdf = df.copy(deep=True)
                if resObj.CalcMetric and elevunits.upper() in ['FT','FEET']:
                    newdf.index.set_levels(newdf.index.levels[1]*FTtoM, level=1, inplace=True)
                    newdf.index.set_names('Elev_m',level=1, inplace=True)
                    hdrDict['elevations_final'] = 'Elev_m'
                    
                if not resObj.CalcMetric and elevunits.upper() in ['M','METER','METERS']:
                    newdf.index.set_levels(newdf.index.levels[1]*MtoFT, level=1, inplace=True)
                    newdf.index.set_names('Elev_ft',level=1, inplace=True)  
                    hdrDict['elevations_final'] = 'Elev_ft'
                
                profunits = hdrDict['temperature'].split('_')[-1]
                hdrDict['temperature_final'] = hdrDict['temperature']
                if resObj.CalcMetric and profunits.upper() in ['DEG_F','DEGF', 'F','FAHRENHEIT']:
                    newdf[hdrDict['temperature']] = [DEGF_to_DEGC(_) for _ in df[hdrDict['temperature']]]
                    newdf.rename(columns={hdrDict['temperature']:'Temp_degC'}, inplace=True)
                    hdrDict['temperature_final'] = 'Temp_degC'
                    
                if not resObj.CalcMetric and profunits.upper() in ['DEG_C','DEGC', 'C','CELSIUS']:
                    newdf[hdrDict['temperature']] = [DEGC_to_DEGF(_) for _ in df[hdrDict['temperature']]]
                    newdf.rename(columns={hdrDict['temperature']:'Temp_degF'}, inplace=True)
                    hdrDict['temperature_final'] = 'Temp_degF'
                    
                
            elif 'Column' in fdata:
                #TODO: 20221005: implement units checks and conversions for thi form of profile handling
                hdrDict = {'column': fdata['Column']}
                if 'Index' in fdata:
                    hdrDict['index'] = fdata['Index']
                    
                if fdata['Column']=='elevations':
                    # assume index is datetime
                    df = pnd.read_table(fp, sep=',', index_col=0, parse_dates=True)
                    df.columns = df.columns.astype(int)
                else:
                    tmp = pnd.read_table(fp, sep=',')
                    ncols = len(tmp.columns)
                    datecols = [x for x in list(range(1,ncols))]
                    df = pnd.read_table(fp, sep=',', index_col=0, parse_dates=datecols)
            else:
                print("No header information provided for temperature profiles")
                print("Assuming columns are elevations and index is date")
                hdrDict = {'column': 'elevations', 'index': 'datetime'}
            
            resObj.Observations.ProfilesDF = newdf
            resObj.Observations.ProfilesColMap = hdrDict
            
def getOpsData(resObj, inputs):
    
    opsdata = inputs['Operations']
    
    # check for flood curve info
    if 'Flood_Curve' in opsdata:
        fcdata = opsdata['Flood_Curve']
        for v in ['smin', 'smax', 'xmin', 'xmax']:
            if v in fcdata:
                resObj.Operations.FloodCurve[v] = fcdata[v]
            else:
                print("ERROR!: Couldn't find `%s` under `Flood_Curve` in input yaml file" %v)
        
    # check for minimum flows
    if 'Min_Flows' in opsdata:
        mfdata = opsdata['Min_Flows']
        for m, v in mfdata.items():
            resObj.Operations.MinFlows[m] = v
    
    # assuming both historical_storage_by_date & historical_release_by_date are included
    # 2021-11-22 @jmg: adding 'if' statement here to check so that it's not an error if not included
    if 'Historical_Storage_by_Date' in opsdata:
        hist_sto = opsdata['Historical_Storage_by_Date']
        hist_sto_data = {}
        for k,v in hist_sto.items():
            fp = v['Path']
            if not os.path.isabs(fp):
                fp = os.path.join(resObj.ProjDir, fp)
            hist_sto_data[k] = pnd.read_table(fp, sep=',', index_col=v['Index'], header=0)  #v['Path']
    else:
        hist_sto_data = {'None': np.nan}
        
    if 'Historical_Release_by_Date' in opsdata:        
        hist_rel = opsdata['Historical_Release_by_Date']
        hist_rel_data = {}
        for k,v in hist_rel.items():
            fp = v['Path']
            if not os.path.isabs(fp):
                fp = os.path.join(resObj.ProjDir, fp)
            
            hist_rel_data[k] = pnd.read_table(fp, sep=',', index_col=v['Index'], header=0)  #v['Path']
    else:
        hist_rel_data = {'None': np.nan}
        
    resObj.Operations.HistStorage = hist_sto_data
    resObj.Operations.HistRelease = hist_rel_data
    
def percentileofscore(a, score, kind='rank'):
    """
    Compute the percentile rank of a score relative to a list of scores.
    A `percentileofscore` of, for example, 80% means that 80% of the
    scores in `a` are below the given score. In the case of gaps or
    ties, the exact definition depends on the optional keyword, `kind`.
    
    Parameters
    ----------
    a: array_like
        Array of scores to which `score` is compared.
    score : int or float
        Score that is compared to the elements in `a`.
    kind : {'rank', 'weak', 'strict', 'mean'}, optional
        Specifies the interpretation of the resulting score.
        The following options are available (default is 'rank'):
          * 'rank': Average percentage ranking of score.  In case of multiple
            matches, average the percentage rankings of all matching scores.
          * 'weak': This kind corresponds to the definition of a cumulative
            distribution function.  A percentileofscore of 80% means that 80%
            of values are less than or equal to the provided score.
          * 'strict': Similar to "weak", except that only values that are
            strictly less than the given score are counted.
          * 'mean': The average of the "weak" and "strict" scores, often used
            in testing.  See https://en.wikipedia.org/wiki/Percentile_rank
    Returns
    -------
    pcos : float
        Percentile-position of score (0-100) relative to `a`.
    
    See Also
    --------
    numpy.percentile
    
    Examples
    --------
    Three-quarters of the given values lie below a given score:
    >>> from scipy import stats
    >>> stats.percentileofscore([1, 2, 3, 4], 3)
    75.0
    With multiple matches, note how the scores of the two matches, 0.6
    and 0.8 respectively, are averaged:
    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3)
    70.0
    Only 2/5 values are strictly less than 3:
    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='strict')
    40.0
    But 4/5 values are less than or equal to 3:
    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='weak')
    80.0
    The average between the weak and the strict scores is:
    >>> stats.percentileofscore([1, 2, 3, 3, 4], 3, kind='mean')
    60.0
    """
    if np.isnan(score):
        return np.nan
    a = np.asarray(a)
    n = len(a)
    if n == 0:
        return 100.0

    if kind == 'rank':
        left = np.count_nonzero(a < score)
        right = np.count_nonzero(a <= score)
        pct = (right + left + (1 if right > left else 0)) * 50.0/n
        return pct
    elif kind == 'strict':
        return np.count_nonzero(a < score) / n * 100
    elif kind == 'weak':
        return np.count_nonzero(a <= score) / n * 100
    elif kind == 'mean':
        pct = (np.count_nonzero(a < score) + np.count_nonzero(a <= score)) / n * 50
        return pct
    else:
        raise ValueError("kind can only be 'rank', 'strict', 'weak' or 'mean'")
        
        
def plot_init_storage_and_profile(resObj):
    import matplotlib.pyplot as plt
    import seaborn as sns
    
    reservoir_shape = [1*(1+0.08)**x for x in range(0,100)]
    reservoir_shape = [(x-min(reservoir_shape))/(max(reservoir_shape)-min(reservoir_shape)) for x in reservoir_shape]
    
    lyrTopElevs = [_.MaxElev for _ in resObj.Layers.values()]
    lyrBotElevs = [_.MinElev for _ in resObj.Layers.values()]
    lyrTemps = [_.Temp for _ in resObj.Layers.values()]
    
    minelev = np.min(lyrBotElevs)
    maxelev = np.max(lyrTopElevs)
    
    mintemp = np.nanmin(lyrTemps)
    maxtemp = np.nanmax(lyrTemps)
    temprng = maxtemp-mintemp
    
    reservoir_shape_rescl = [minelev + (maxelev-minelev)*x for x in reservoir_shape]
    
    # for l, v in resObj.Layers.items():
    #     print("  %d       %0.4f     %0.1f      %0.1f       %0.1f" %(l, v.Temp, v.TotE, v.Vol, v.MaxVol)  )

    resObj.getWSE()
    
    with sns.plotting_context('notebook',font_scale=1.1):
        fig, ax = plt.subplots(1,1, figsize=(12,8))
        ax.plot(range(0,100), reservoir_shape_rescl, c='k', lw=2)
        plt.xlim((0,100))
        for l in lyrTopElevs:
            xmax = np.interp(l, reservoir_shape_rescl, range(100))
            plt.axhline(xmin=0, xmax = xmax/100, y=l, ls='--', c='0.2', lw=0.5)
            
        # add stage
        ax.axhline(resObj.WSE, lw=2, c='darkblue')
        
        # profile
        ax2 = ax.twiny()
        if resObj.CalcMetric:
            ax2.set_xlim((max(0,mintemp-0.1*temprng), maxtemp+0.5*temprng))
        else:
            ax2.set_xlim((max(32,mintemp-0.1*temprng), maxtemp+0.5*temprng))
        ax2.plot(lyrTemps, lyrTopElevs, c='darkred')
        
def plot_layers_storage_area(resObj):
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.patches import Rectangle, Circle, Ellipse
    import matplotlib.ticker as tkr
    
    reservoir_shape = [1*(1+0.08)**x for x in range(0,100)]
    reservoir_shape = [(x-min(reservoir_shape))/(max(reservoir_shape)-min(reservoir_shape)) for x in reservoir_shape]
    
    lyrTopElevs = [_.MaxElev for _ in resObj.Layers.values()]
    lyrBotElevs = [_.MinElev for _ in resObj.Layers.values()]
    lyrAreas = [_.MaxArea for _ in resObj.Layers.values()]
    lyrVols = [_.MaxVol for _ in resObj.Layers.values()]
    
    minelev = np.min(lyrBotElevs)
    maxelev = np.max(lyrTopElevs)
    minVol = np.min(lyrVols)
    maxVol = np.max(lyrVols)
    minAreas = np.min(lyrAreas)
    maxAreas = np.max(lyrAreas)
    
    #mintemp = np.nanmin(lyrTemps)
    #maxtemp = np.nanmax(lyrTemps)
    #temprng = maxtemp-mintemp
    
    reservoir_volume_rescl = [minVol + (maxVol-minVol)*100]
    reservoir_shape_rescl = [minelev + (maxelev-minelev)*x for x in reservoir_shape]
    xmax = maxVol *1.1
    # for l, v in resObj.Layers.items():
    #     print("  %d       %0.4f     %0.1f      %0.1f       %0.1f" %(l, v.Temp, v.TotE, v.Vol, v.MaxVol)  )

    #resObj.getWSE()
    
    with sns.plotting_context('notebook',font_scale=1.1):
        fig, ax = plt.subplots(1,1, figsize=(12,8))
        #ax.plot(range(0,100), reservoir_shape_rescl, c='k', lw=2)
        plt.xlim((-10000,xmax))
        plt.ylim((minelev, maxelev))
        aspectratio = (maxVol+10000)/(maxelev-minelev)
        lyridx = 0
        for lb,lt,lv in zip(lyrBotElevs, lyrTopElevs, lyrVols):
            #print(lb, lt, lv)
            # make a rectangle that represents the volume for a given layer
            thisRectangle = Rectangle((0,lb), lv,(lt-lb), color='darkblue', alpha=0.5)
            ax.add_patch(thisRectangle)
            #xmax = np.interp(l, reservoir_shape_rescl, range(100))
            plt.axhline(xmin=0, xmax = lv, y=lt, ls='--', c='0.2', lw=0.5)

            ax.annotate(f'{lyridx}',(xmax*.95,(lt+lb)/2), xycoords='data',
                        xytext=(xmax*0.95,lb), textcoords='data', 
                        fontsize=10)
            lyridx+=1
            
        # add TCD, penstock, and river outlet elevations
        for i, swdo in resObj.Outlets.items(): #selective withdrawal outlets
            #print(swdo.BotElevFt)
            swdRect = Rectangle((-7500, swdo.BotElevFt), 5000, 
                                swdo.TopElevFt-swdo.BotElevFt, edgecolor='k', 
                                facecolor='0.6', alpha=0.5)
            ax.add_patch(swdRect)
            
        # add the penstock location (centerline)
        penRect= Rectangle((-10000,resObj.PenstockElevation-2),10000.,
                           4., color='k')
        ax.add_patch(penRect)
        
        # add the river outlet elevations
        for i, rivo in resObj.RiverOutlets.items():
            rivCir = Ellipse((-5000, rivo.CtrElev), rivo.DiamFt/2*aspectratio, 
                             rivo.DiamFt/2, edgecolor='darkred', 
                                facecolor='darkred', alpha=0.9)
            ax.add_patch(rivCir)
            
        ax.xaxis.set_major_locator(tkr.MultipleLocator(50000.))
        ax.xaxis.set_minor_locator(tkr.MultipleLocator(5000.))
        ax.yaxis.set_minor_locator(tkr.MultipleLocator(10.))     
        
        ax.set_ylabel('Elevation (feet)')
        ax.set_xlabel('Model layer storage volume (acre-feet)')
        
        
        # add a custom legend
        ax.add_patch(Rectangle((1.02, 0.7), width=0.02, height=0.075, clip_on=False,
                               transform=ax.transAxes,edgecolor='k', 
                               facecolor='0.6', alpha=0.5) )
        ax.annotate('TCD\nGates\nLocations', (1.03, 0.72), xycoords='axes fraction',
                    xytext=(1.045, 0.70), textcoords='axes fraction',
                    fontsize=12)
        ax.add_patch(Rectangle((1.014, 0.62), width=0.03, height=0.01, clip_on=False,
                               transform=ax.transAxes, color='k'))
        ax.annotate('Penstock\nLocation', (1.04, 0.62), xycoords='axes fraction',
                    xytext=(1.047, 0.60), textcoords='axes fraction', fontsize=12)
        
        # # add stage
        # ax.axhline(resObj.WSE, lw=2, c='darkblue')
        
        # # profile
        # ax2 = ax.twiny()
        # if resObj.CalcMetric:
        #     ax2.set_xlim((max(0,mintemp-0.1*temprng), maxtemp+0.5*temprng))
        # else:
        #     ax2.set_xlim((max(32,mintemp-0.1*temprng), maxtemp+0.5*temprng))
        # ax2.plot(lyrTemps, lyrTopElevs, c='darkred')