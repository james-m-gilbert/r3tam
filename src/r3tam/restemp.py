# -*- coding: utf-8 -*-
"""
Originally created on Mon May  4 2020

@author: jgilbert
"""

import os
import sys
import copy
import pandas as pnd
import datetime as dt
import numpy as np
from collections import OrderedDict as Odict
import time as time
import yaml
import r3tam.restemp_util as util

import pickle
import datetime as dt

import r3tam.surface as surface
import r3tam.outflow as outflow

from r3tam.constants import DENSE_PREC, DENSE_TOL, FTtoM, cw, AFtoM3


class res_data:
    """
    Represents data associated with a reservoir.

    Currently this includes the elevation-area relationship (as a dataframe),
    and the elevation-storage relationship (also as as dataframe).

    Attributes
    ---------
    ElevArea : pandas.DataFrame
        dataframe of elevation-area relationship; index is elevation, value 
        column is area at that elevation

    ElevStorage : pandas.DataFrame
        dataframe of elevation-reservoir storage relationship; index is 
        elevation, value column is the corresponding volume or reservoir storage

    Methods
    -------
    None.
    """

    def __init__(self):
        """
        Initializes the attributes for the res_data object.

        Creates empty variables to be populated when the model is initialized
        with user inputs.

        Parameters
        ----------
        None.

        Returns
        -------
        None.

        """
        self.ElevArea = pnd.DataFrame()
        self.ElevStorage = pnd.DataFrame()


class simspec:
    """
    A class representing the control or configuration options for a simulation.


    Properties
    ----------
    StartDate : datetime.datetime
        Simulation start date, as a datetime object, automatically updated
        based on model initilization
    EndDate : datetime.datetime
        Simulation end date, as a datetime object, automatically updated
        based on model initilization

    """

    def __init__(self):
        """
        Initializes values for simspec object.

        Initialized variables include things like start and end dates, 
        initialization methods, gate change parameters for dynamic operations, 
        units specification, reinitialization options, and whether or not to 
        calculate evaporation.


        Attributes
        ----------
        StartDate_str : str
            string representation of the simulation start date
        EndDate_str : str
            string representation of the simulation end date
        DELT_DAY : int, default 1
            simulation time step in days
        DELT_SEC : int
            simulation time step in seconds, automatically calculated from DELT_DAY
            as DELT_DAY*86400
        InitMethod: str
            User-specified method for initializing reservoir temperature profiles; 
            As of 2023-03-01, this should be 'from_obs' if using observed 
            temperature profiles as input or 'None' if passing in specified
            profiles from an outside program.
        IsHindcast : bool
            True or False to indicate if the simulation is a hindcast - meaning
            all gate operations are fully specified and no dynamic operations
            are requried; If set to False, the gate operations logic will be
            used.
        CalcGateOps : bool
            True or False to indicate if the gate operations will be calculated
            to meet a downstream temperature targer rather than using specified
            gate operations as input. If True, all gate operations must be fully
            specified in the model inputs. If False,the gate operations logic 
            will be employed to meet a user-specified downstream target.
        GateChgFreq : int
            Value in days that sets the minimum interval at which a selective
            withdrawal gate change operation can be made when dynamic gate 
            operations are being simulated (i.e. when CalcGateOps=True)
        DaysSinceLastGateChange : int, default 0
            Counter to track the number of days since the last selective 
            withdrawal gate change was made. Initialized at zero.
        CalcEvaporation : bool, default False
            Flag to turn on dynamic evaporation calculation (set to True). If
            set to False, the provided input must included an appropriate 
            evaporation time series.
        CalcReleaseSched : bool, default False
            Flag to turn on dynamic release schedule calculation. This function
            is currently in development and should not be used.
        Reinitialize : bool, default False
            Flag to set whether the temperature profile should be reinitialized
            at some point in the simulation. This may be required if running
            multi-year simulations but it is desired to set a specified initial 
            condition at the beginning of each year. If set to True, the 
            ReinitMonth and ReinitDay variables must also be set.
        ReinitMonth : int, default 0
            If Reinitialize is set to True, the month at which the temperature
            profile will be reset. Must be accompanied by a ReinitDay value
        ReinitDay : int, default 0
            If Reinitialize is set to True, the day at which the temperature
            profile will be reset. Must be accompanied by a ReinitMonth value
        Profile_Temp_Units : str, default 'degF'
            Indicates the temperature units of the profile data provided as 
            input to the model.
        Simulation_Temp_Units : str, default 'degF'
            Indicates the temperature units used in the simulation. 

            **NOTE:  Code was originally written for calculations in degF, and 
            functionality for degC has not been fully evaluated at this time.
            It is strongly recommended that `degF` be used**
        """

        self.StartDate_str = ''
        self.EndDate_str = ''
        self.DELT_DAY = 1
        self.DELT_SEC = self.DELT_DAY*86400.
        self.InitMethod = ''
        # self.IsHindcast = True #<-- Deprecated in favor of CalcGateOps
        self.CalcGateOps = False
        self.GateChgFreq = 5  # how frequently can a gate change be made, in days
        # how long has it been since a gate change was last made?
        self.DaysSinceLastGateChange = 0
        # calculate evaporation based on air temp, surface temp and wind?
        self.CalcEvaporation = False
        self.CalcReleaseSched = False
        self.Reinitialize = False
        self.ReinitMonth = 0
        self.ReinitDay = 0
        self.Profile_Temp_Units = 'degF'
        self.Simulation_Temp_Units = 'degF'

    @property
    def StartDate(self):
        """
        Represent simulation start date as datetime object

        Returns
        -------
        StartDate : datetime.datetime object
            Simulation start date
        """
        return(dt.datetime.strptime(self.StartDate_str, '%m/%d/%Y'))

    @property
    def EndDate(self):
        """
        Represent simulation end date as datetime object

        Returns
        -------
        EndDate : datetime.datetime object
            Simulation end date
        """
        return(dt.datetime.strptime(self.EndDate_str, '%m/%d/%Y'))

# ???? Used in early development - keep?
# class calibspec:
#     def __init__(self):
#         self.Method = ''
#         self.MaxIter = 2000
#         self.Display = True
#         self.FATOL = 0.00001
#         self.XATOL = 0.00001
#         self.InitGuess = {}
#         self.Bounds = {}


class inputTS:
    """
    A class representing generic input time series datasets.

    Used for organizing input data like meteorology, gate operations, release
    schedules, and others. 
    """

    def __init__(self):
        """
        Initialize empty input time series data structures

        Returns
        -------
        None.

        """
        self.DataFrame = None
        self.ColumnMap = None


class observations:
    """
    A class representing observation datasets for reservoirs

    Used for organizing observed temperature profile and release data. 
    """

    def __init__(self):
        """
        Initialize empty observation data structures

        Returns
        -------
        None.
        """
        self.OutflowDF = None
        self.OutflowColMap = None
        self.ProfilesDF = None
        self.ProfilesColMap = None


class operations:
    """
    A class representing reservoir operations time series

    Used for organizing release and gate operations data. Originally intended
    to also include data for dynamically calculating release schedules, but is
    now primarily used to store gate operations.    
    """

    def __init__(self):
        """
        Initialize empty operations data structures.

        Returns
        -------
        None.

        """
        self.HistStorage = None
        self.HistRelease = None
        self.FloodCurve = {'smin': 0., 'smax': 0., 'xmin': 0., 'xmax': 0.}
        self.MinFlows = {}
        self.GateOps = {}
        self.sidegate_only_cntr = 0


class Res(object):
    """
    The ResTemp model object.

    The core component of the ResTemp model, this object includes the functions
    to initialize and execute a simulation.    
    """

    def __init__(self):
        """
        Initialize all the properties and variables needed for a ResTemp simulation.

        Properties
        ----------
        C1 : float
            Parameter C1 to adjust air-water temperature exchange efficiency, 
            values should be between 0 and 1 (inclusive)
        C2 : float
            Parameter C2 to adjust solar radiation energy transfer efficiency, 
            values should be between 0 and 1 (inclusive)
        C3 : float
            Parameter C3 to adjust latent heat energy transfer efficiency, 
            values should be between 0 and 1 (inclusive)
        C4 : float
            Parameter C4 to adjust inflow mixing efficiency, values should be 
            between 0 and 1 (inclusive)
        C5 : float
            Parameter C5 to adjust internal diffusion efficiency, values should 
            be between 0 and 1 (inclusive)
        Cd : float
            Wind drag coefficient
        SeasonalAlbedo : bool
            Flag to indicate whether seasonal-variable albedo calculations are 
            turned on
        MeanAlbedo : float, default 0.08
            Mean albedo around which seasonal variation is calculated; default
            value of 0.08 comes from Hipsey et al 2017 & Hamilton & Schaldow 1997
        AlbedoAmplitude : float, default 0.02
            Magnitude of maximum variation in seasonally-varying albedo. Default
            value comes from Hipsey et al 2017 & Hamilton & Schaldow 1997
        CriticalDepth: float
            Distance down from water surface through which surface energy transfer
            occurs
        CriticalDepthFactor: float
            Implemented for time-variable critical depth testing - NOT USED
        Layers : collections.OrderedDict
            Ordered dictionary in which layer information is stored
        MaxVol : float, default 4442000.0
            Maximum storage volume of reservoir. Default is for Shasta. [acre-feet]
        MinVol : float, default 115700.0
            Minimum storage volume of reservoir. Default is for Shasta. [acre-feet]
        Outlets : collections.OrderedDict
            Ordered dictionary in which selective withdrawal outlet information 
            is stored
        PenstockLimit : float, default 0.0
            Maximum release capacity through all penstocks, combined
        PenstockLimit_Units : str, default 'cfs'
            Units of the PenstockLimit variable
        PenstockElevation : float, default 0
            Elevation of the penstock
        PenstockElevUnits : str, default 'ft'
            Units of value set for PenstockElevation (should be 'ft' or 'm')
        PenstockElevation_m : float, default 0.0
            Penstock elevation in units of meters
        RiverOutletLimit : float, default 75000.
            Total maximum capacity of river outlets.
        RiverOutletLimit_Units : str, default 'cfs'
            Units of value set in RiverOutletLimit
        SpillwayOutletLimit : float, default 75000.
            Total maximum capacity of the spillway.
        SpilwayOutletLimit_Units : str, default 'cfs'
            Units of value set in SpillwayOutletLimit.        
        LogFilePath : str
            Full file path to location of a log file - used for debugging
        TimeStep : int
            Integer indicating the time step - starts at -1 such that the first
            time step becomes 0 when advanced.
        _time : int
            Secondary time step counter - used internally
        TimeStepDate : datetime.datetime
            Datetime object tracking the date of the time step
        Leakages : collections.OrderedDictionary 
            Ordered dictionary in which discrete leakage areas or elevations
            are stored
        RiverOutlets : collections.OrderedDictionary 
            Ordered dictionary in which river outlet information is stored
        ReservoirData : restemp.res_data
            A reservoir data (res_data) class object for storing static
            properties of the reservoir (elevaiton-storage-area curves)
        SimulationSpecs : restemp.simspec
            A simulation spec (simspec) object for storing configuration information
            for a run (dates, run options, etc)
        Inflow : restemp.inputTS
            An input time series class object for defining inflow time series
        Outflow : restemp.inputTS
            An input time series class object for defining reservoir release
            time series
        Met : restemp.inputTS
            An input time series class object for defining meteorology time series
        Observations : restemp.observations
            An observations class object for defining profile and release 
            temperature observations
        Operations : restemp.operations
            An operations class object for defining and storing gate operations
            and release data
        Evaporation : list
            A list for storing evaporation values calculated during simulation
        TimeStepList : list
            A list for storing the successfully completed time steps
        ProjDir : str, default ''
            Full file path to the directory of this simulation inputs
        SeasonalRad : float, default 0.9
            NOT CURRENTLY USED; Originally implemented for testing purposes
            Factor for adjusting seasonal radiation efficiency; 
        Wind_Func_Params : dict
            Dictionary for storing wind function parameters.
        Storage : float
            The reservoir storage value for the current time step. This value
            is overwritten each time step.
        WSE : float
            Water surface elevation, calculated multiple times each time step
            based on the storage mass balance; default units of 'ft'
        Simulation_Results : dict
            Dictionary storing lists of values calculated for each output variable.
            Includes storage, release temperautre, temperature profiles, 
            evaporation, releases by outlet type, and target temperatures used.
        SurfHeat : dict
            Dictionary used internally for storing components fo the surface
            heat budget
        Debug : dict
            Dictionary storing the debug settings for General and Release 
            processes. Default values for each category are 0. Setting either
            to an integer value of 1 or higher will increase the number and
            detail of debug information produced during a simulation.
            
        LP_Opt_Blending: boolean
            Flag to indicate whether blending across open gates is determined
            via a linear program optimization scheme; this requires a tailwater
            target temperature be provided; this is set in the input file under
            the outflow section

        Returns
        -------
        None.

        """

        self.C1 = 0.5
        self.C2 = 0.5
        self.C3 = 0.5
        self.C4 = 0.5
        self.C5 = 0.5
        self.Cd = ''  # wind drag coefficient
        self.SeasonalAlbedo = False
        self.MeanAlbedo = 0.08  # based on Hipsey et al 2017 & Hamilton & Schaldow 1997
        self.AlbedoAmplitude = 0.02  # based on Hipsey et al 2017 & Hamilton & Schaldow 1997
        self.CriticalDepth = 30  # ft
        self.CriticalDepthFactor = 1.
        self.Layers = Odict()
        self.MaxVol = 4552000.  # ac-ft
        self.MinVol = 115700.
        self.Outlets = Odict()  # TCD outlets
        self.PenstockLimit = 0.
        self.PenstockLimit_Units = 'cfs'
        self.PenstockElevation = 0.
        self.PenstockElevUnits = 'ft'
        self.PenstockElevation_m = 0.
        self.RiverOutletLimit = 75000.
        self.RiverOutleLimit_Units = 'cfs'
        self.SpillwayOutletLimit = 300000.
        self.SpilwayOutletLimit_Units = 'cfs'
        # self.TempProf = [] # TODO - is this needed anymore?
        self.LogFilePath = ''
        self.TimeStep = -1
        self._time = -1
        self.TimeStepDate = dt.datetime(2000, 10, 1, 0, 0)
        self.Leakages = Odict()
        self.RiverOutlets = Odict()  # direct river release (bypassing power plant) outlets
        self.ReservoirData = res_data()
        self.SimulationSpecs = simspec()
        self.Inflow = inputTS()
        self.Outflow = inputTS()
        self.Met = inputTS()
        self.Observations = observations()
        self.Operations = operations()
        self.Evaporation = []
        self.Extra_Outflow = [] # track additional extracted water
        self.TimeStepList = []
        self.ProjDir = ''
        # self.Calibration = calibspec() #??? keep this option in the code?
        self.ElevDependentCD_k = 0.
        self.ElevDependentCD_A = 0.
        self.ElevDependentCD_x0 = 0.
        self.SeasonalRad = 0.9
        self.Wind_Func_Params = {}

        self.LP_Opt_Blending = False # <-- use linear program optimization to blend selective withdrawals to meet target
                                    # if true, need to provide a target, even if gate operations are set
            
        self.GateDict = {} # configuration of selective withdrawal gates for current time step
        self.PrevGateDict = {} #config of selective withdrawal gates for previous time step
        self.PointSinkFracs = {}  # fraction of flow through point sink levels, by time step

        # state variabls and storing simulated results
        self.Storage = 0.
        self.WSE = 0.
        self.Simulation_Results = {'ReleaseTemp': [], 'ReleaseVol': [],
                                   'ReleaseTempTarget': [],
                                   'TempTargetTol': [],
                                   'Storage': [], 'Stage': [],
                                   'Profiles': [],
                                   'StorageDF': None,
                                   'ReleaseDF': None,
                                   'ProfilesDF': None,
                                   'Release_by_Outlet': {},
                                   'Release_by_Outlet_DF': None}

        self.SurfHeat = {'Evap_Wm2': [], 'Evap_m3': [], 'Evap_J': []}
        self.Debug = {'Release': 0, 'General': 0}

    @property
    def nLyrs(self):
        return(len(self.Layers))

    @classmethod
    def initialize_model(cls, config_fp, **kwargs):
        """ 
        Create and initialize a ResTemp model object from information in a config file


        Parameters
        ----------

        config_fp : string
              Full file path to configuration file - can either be a YAML
              specification or a previously-saved *.pkl file
        **kwargs 
            profTable : specify a profile table to be used in setting the initial
            temperature profile

            initSto : Value to use in setting the initial reservoir storage

            profile_temp_units : "degF" or "degC" to indicate which units the
            provided temperature profile is in

        Raises
        ------
        FileNotFoundError
            - If configuration filepath specified by `config_fp` is not found

        Returns
        --------
        resmod
            A new instance of a ResTemp reservoir temperature model   
        """

        print(f"Loading configuration file: {config_fp}")
        if not os.path.exists(config_fp):
            raise FileNotFoundError(config_fp)
            return

        if os.path.splitext(config_fp)[1] in ['.yaml', '.yml']:
            with open(config_fp) as stream:
                inputs = yaml.load(stream, Loader=yaml.SafeLoader)

            #projDir = os.path.dirname(os.path.dirname(config_fp))
            projDir = os.path.dirname(config_fp)
            runName = inputs['RunName']

            print("working directory set as: %s" % projDir)

            resmod = cls()  # instantiate reservoir model object

            resmod.ConfigFP = config_fp

            util.setupSim(resmod, inputs, projDir)

            util.setParams(resmod, inputs)

            util.setupRes_ElevSpec(resmod, inputs)

            util.setupOutlets(resmod, inputs)

            util.getTimeSeries(resmod, inputs)

            util.getObservations(resmod, inputs)

            if 'Operations' in inputs:
                util.getOpsData(resmod, inputs)

            initDateStr = resmod.SimulationSpecs.StartDate.strftime('%Y-%m-%d')

            resmod.KBSW = 1
            resmod.KTSW = resmod.nLyrs - 1

            # assign gates to layers
            resmod.assignResLayers()

            # set up the diffusion adjacency information
            resmod.getDiffNeighbors()

            if 'profTable' in kwargs:
                initTempProf = kwargs['profTable']
                initTempProf.rename(
                    columns={"Var1": "Temp_degC"}, inplace=True)
                print(initTempProf)

            if 'initSto' in kwargs:
                initVol = kwargs['initSto']

            if 'profile_temp_units' in kwargs:
                tempUnits = kwargs['profile_temp_units']
                resmod.SimulationSpecs.Profile_Temp_Units = tempUnits

            if resmod.SimulationSpecs.InitMethod == 'from_obs':
                outProfiles = resmod.Observations.ProfilesDF
                try:
                    if resmod.Observations.ProfilesColMap['column'] == 'elevations':
                        initTempProf = outProfiles.iloc[outProfiles.index.get_loc(
                            initDateStr, method='nearest')]
                    else:
                        outProfiles = outProfiles.T
                        initTempProf = outProfiles.iloc[outProfiles.index.get_loc(
                            initDateStr, method='nearest')]

                    # stoVar = resmod.Observations.OutflowColMap['obs_storage_final']  # storage moved to the outflow category

                    # initVol = resmod.Outflow.DataFrame.loc[initDateStr,  stoVar] # storage moved to outflow category
                    #initVol = resmod.Observations.OutflowDF.loc[initDateStr,  stoVar]

                except Exception:
                    print("assigning profile assuming 'tidy' format")
                    try:
                        # try assigning profiles assuming a long/tidy format (indexes are date and elevation, values temperature)
                        #idx1 = outProfiles.index.get_level_values(
                        #    0).unique().get_loc(initDateStr, method='nearest')
                        
                        idx1 = outProfiles.index.get_level_values(0).unique().get_indexer([initDateStr], method='nearest')
                        
                        profDat = outProfiles.index.get_level_values(0).unique()[
                            idx1]

                        initTempProf = outProfiles.loc[profDat, :]

                        # stoVar = resmod.Observations.OutflowColMap['obs_storage_final']  # storage moved to the outflow category

                        # initVol = resmod.Outflow.DataFrame.loc[initDateStr,  stoVar] # storage moved to outflow category
                        #initVol = resmod.Observations.OutflowDF.loc[initDateStr,  stoVar]

                    except Exception:
                        print(
                            "Couldn't initialize from array-style or long-style profile tables")
                        # return

                finally:

                    # storage moved to the outflow category
                    stoVar = resmod.Observations.OutflowColMap['obs_storage_final']

                    # initVol = resmod.Outflow.DataFrame.loc[initDateStr,  stoVar] # storage moved to outflow category
                    initVol = resmod.Observations.OutflowDF.loc[initDateStr,  stoVar]

            else:
                print("Trying to initialize from data passed in from external source")

            util.initVolTemp(resmod, initVol, initTempProf,
                             profileTempUnits=tempUnits, profileElevUnits='ft')
            resmod.updateRho()
            resmod.resVolume()

            # save for future use/refernce
            pfp = os.path.join(projDir, runName+'.pkl')

            with open(pfp, 'wb') as of:
                pickle.dump(resmod, of)

            printLayTempVolTotE(resmod)
            return(resmod)

    def advance_restemp(resmod):
        """
        Calculate vertical temperature profile based on meteorology and inflows

        Advances the ResTemp (Res) model one time step considering only meteorological
        boundary conditions and inflows. This determines the temperature profile
        prior to selective withdrawal and releases. The `advance_swd` function
        performs the selective withdrawal and finalizes a time step.


        Parameters
        ----------
        resmod : Res model object
            An instance of the ResTemp model that has been initialized with
            all required information - temperature profile, volume, and 
            boundary condition time series.

        Returns
        -------
        None.

        """
        resmod._time += 1
        resmod.TimeStep += 1  # resmod.TimeStep
        resmod.TimeStepDate = resmod.SimDates[resmod.TimeStep]
        resmod.DayOfYear = resmod.TimeStepDate.timetuple().tm_yday

        # check if needs to be reinitialized
        if resmod.SimulationSpecs.Reinitialize and resmod.TimeStep > 1:
            if (resmod.SimulationSpecs.ReinitMonth == resmod.TimeStepDate.month) & \
               (resmod.SimulationSpecs.ReinitDay == resmod.TimeStepDate.day):  # then reinitialize
                resmod.reinitialize(return_obs_prof=False,
                                    show_init_profile=False)

        met_vec = resmod.Met.DataFrame.loc[resmod.TimeStepDate]
        inflow_vec = resmod.Inflow.DataFrame.loc[resmod.TimeStepDate]
        release_vec = resmod.Outflow.DataFrame.loc[resmod.TimeStepDate]

        # start by ensuring we have the current water surface elevation (WSE)
        # also adds 'TopLyr' index to resmod object
        resmod.getWSE()

        # do surface interaction first - can select between simple
        # or energy budget
        if resmod.Surface == 'simple':

            # set the fractions of top layers over which heat exchange will occur
            surface.distDepth(resmod)

            # calc vapor pressure - puts result in resmod.Met.VaporPress dictionary
            surface.calc_ea_es(resmod, met_vec)

            # calc wind function
            surface.wind_func(resmod, met_vec, method='W2', TairUnits='degF')

            # do wind mixing
            surface.wind_mixing(resmod, met_vec, DensDiffTol=0.001, debug=0)

            # solar radiation
            surface.simple_solrad(resmod, met_vec)
            #surface.solrad_w2(resmod, met_vec)
            #surface.solrad_linear(resmod, met_vec)

            # air-water interchange
            surface.air2water(resmod, met_vec)

            # do evaporation
            if resmod.SimulationSpecs.CalcEvaporation:
                surface.calc_evap(resmod, met_vec)
            else:
                if met_vec[resmod.Met.ColumnMap['evap']] > 0.:
                    # calculate the energy component of evaporation from prescribed evap amount
                    surface.calc_latentheat_from_evap(
                        resmod, met_vec, debug=False)
                else:
                    resmod.SurfHeat['Evap_m3'].append(0.)
                    resmod.SurfHeat['Evap_Wm2'].append(0.)
                    resmod.SurfHeat['Evap_J'].append(0.)
            # now remove the evaporated water from the surface and update temps
            surface.remove_evap(resmod, method='simple',
                                debug=resmod.Debug['General'])

            # precipitation
            surface.precipDepth(resmod, met_vec)

        else:

            # calc vapor pressure - puts result in resmod.Met.VaporPress dictionary
            surface.calc_ea_es(resmod, met_vec)

            # calc wind function
            fw = surface.wind_func(
                resmod, met_vec, method='W2', TairUnits='degF')

            # calc evap
            surface.calc_evap(resmod, met_vec)

            # sensible heat/conductive flux
            surface.calc_sensible()

            raise NotImplementedError(
                "Energy balance method not implemented yet..stay tuned!")
            
        
        # in cases where there are additional extractions (direct deliveries or
        # use from the reservoir), update the volume here
        

        # now do the internal mixing
        thisinflow = inflow_vec[resmod.Inflow.ColumnMap['inflow']]
        thisintemp = inflow_vec[resmod.Inflow.ColumnMap['inflowTemp']]

        # distribute inflows, do mixing
        resmod.distInflowMix(thisinflow, thisintemp)
        
        # remove direct diversions/depletions
        outflow.extra_outflow_dist(resmod, release_vec[resmod.Outflow.ColumnMap['extraOutflow']])

        # distribute volumes
        # TODO: change to consistent debug level instead of true/false
        resmod.distVol(debug=False)

        # apply diffusion
        resmod.diffuseDist()

        # update densities
        resmod.updateRho()

        # do convective mixing check - do less dense layers underly denser ones?
        stabcheck = resmod.checkDensityProfile()
        iterLimit = 1000
        iters = 0
        while (stabcheck > 0) & (iters < iterLimit):
            resmod.checkStable()
            stabcheck = resmod.checkDensityProfile()
            iterLimit += 1

        resmod.resVolume()
        resmod.getWSE()

    def advance_swd(resmod, final=False, return_vals=False, **kwargs):
        '''
        Advance Selective WithDrawal

        Do release/outflow for an assigned gate operation
        if final=True, then remove water from reservoir and recalculate
        thermal profile
        otherwise just calculate what the outflow temperature would be

        Parameters
        ----------
        resmod : Res model object
            An instance of the ResTemp model that has been initialized with
            all required information - temperature profile, volume, and 
            boundary condition time series.
        final : bool, default False
            Flag indicating whether this function call should finalize the
            time step and calculate the profile with releases from a specified
            gate configuration. If True, values are stored in the internal 
            restemp.Res.SimulationResults dictionary
        return_vals : bool, default False
            Flag indicating whether the result of a selective withdrawal calculation
            should be returned from this function call. 
        **kwargs
            check_gate_calcs
            coupling_model
            gate_dict
            temp_target
            riv_dict
            bypass_fraction

        Returns
        -------
        [totQ, outT, outE] : list
            If `return_vals` set to True, returns a list of three values: 
            totQ is the total release volume; outT is the temperature of that
            relased volume; and outE is the thermal energy of that release
            (product of volume and temperature)

        '''
        outflowvec = resmod.Outflow.DataFrame.loc[resmod.TimeStepDate]

        if 'check_gate_calcs' in kwargs:
            check_gate_calcs = kwargs['check_gate_calcs']

        # check if gate operations need to be calculated or they can just be used as assigned
        if resmod.SimulationSpecs.CalcGateOps and not final:
            # # assume you've passed in the simulated temperature at teh target location
            # # from the last time step as kwarg 'prev_sim_targ_temp'
            # prev_targ_temp = kwargs['prev_sim_targ_temp']
            # also assume you;ve passed in thetemperature target time series
            # and the model coupling object
            cplmod = kwargs['coupling_model']
            # ##ttarg_ts = kwargs['temp_target_timeseries']

            spec_targ_col = cplmod.Temperature_Target.Target_TS.ColumnMap['ttarg_specified']
            def_targ_col = cplmod.Temperature_Target.Target_TS.ColumnMap['ttarg_default']
            if spec_targ_col != 'Blank_Spec_Target':
                use_targ_col = spec_targ_col
            else:
                use_targ_col = def_targ_col

            this_date = resmod.TimeStepDate
            #print(f"date at this point is: {this_date}")

            this_ttarg = cplmod.Temperature_Target.Target_TS.DataFrame.loc[this_date, use_targ_col]
            #print(f"temp targ this time step: {this_ttarg}")
            rivDict = {k: int(
                outflowvec[v]) for k, v in resmod.Outflow.ColumnMap['rivOutlets'].items()}

            prev_targ_sim_list = []  # simulated values from last few time steps at target location
            prev_targ_list = []     # target temperature value at target
            prev_temp_error = []
            for prevoffset in range(resmod.SimulationSpecs.GateChgFreq*2, 0, -1): #[5, 4, 3, 2, 1]:
                prev_date = this_date-dt.timedelta(prevoffset)
                if prev_date > cplmod.SimDates[0]:
                    prev_targ_sim = cplmod.Temperature_Target.Sim_TS[prev_date]

                    prev_targ_sim_list.append(prev_targ_sim)

                    prev_targ = cplmod.Temperature_Target.Target_TS.DataFrame.loc[
                        prev_date, use_targ_col]
                    prev_targ_list.append(prev_targ)
                    prev_temp_error.append(prev_targ_sim-prev_targ)
                else:
                    prev_targ_sim = this_ttarg
                    prev_temp_error.append(0)

            #print("\n previous time step errors:")
            #print(f"\t\t {prev_temp_error}\n")
            increasing_errors = all(i < j for i, j in zip(
                prev_temp_error, prev_temp_error[1:]))
            
            if this_date.month>8:
                increasing_errors2 = sum([abs(iv) for iv in prev_temp_error[0:4]])>1.25
            else:
                increasing_errors2 = sum([abs(iv) for iv in prev_temp_error])>3.25 #len(prev_temp_error)*0.2
            
            if resmod.TimeStep <= 0:
                tdiff = 0
            else:
                tdiff = prev_temp_error[-1]

            #print(f"temp difference: {tdiff}")

            this_pos_temp_tol = cplmod.Temperature_Target.Target_TS.DataFrame.loc[this_date,
                                                                                  'Temp_Target_upper_tol_degC']
            this_neg_temp_tol = cplmod.Temperature_Target.Target_TS.DataFrame.loc[this_date,
                                                                                  'Temp_Target_lower_tol_degC']
            targ_data = [this_ttarg, this_pos_temp_tol, this_neg_temp_tol]
            # is the downstream temperature currently warmer or colder than the target?
            # dsTerr = ds_ccr_temps[-1] - thisTempTarg  # positive val = warmer than target

            if resmod.WSE < resmod.PenstockElevation + 10.0: #trying to ensure an outflow for each time step
                # water level is too low to use TCD and penstocks - have to use lowest
                # river outlet
                rivDict = {k: 0 for k in resmod.RiverOutlets}
                rivDict[min(resmod.RiverOutlets)] = 1

                thisOutFlowTot = outflowvec[resmod.Outflow.ColumnMap['outflow']]

                gateDict = {k: 0 for k,
                            v in resmod.Outflow.ColumnMap['gates'].items()}
                
                resmod.GateDict= gateDict
                
                thisRivOut = thisOutFlowTot

                bypassFrac = 1

                if check_gate_calcs:
                    return(['done', resmod.GateDict, targ_data, tdiff, rivDict, bypassFrac])
            else:
                
                thisOutFlowTot = outflowvec[resmod.Outflow.ColumnMap['outflow']]
                thisRivOut = outflowvec[resmod.Outflow.ColumnMap['rivOutFlow']]
                if thisRivOut > 0.:
                    bypassFrac = thisRivOut/thisOutFlowTot
                    
                    # send river outlet water first out through the top-most
                    # outlet, then incrementally down
                    water_to_dist = thisRivOut
                    
                    for ro in resmod.RiverOutlets: # from top down
                        this_capacity_vol = resmod.RiverOutlets[ro].Capacity*resmod.SimulationSpecs.DELT_SEC/43560. #conver to AF volume
                        if resmod.WSE > resmod.RiverOutlets[ro].MinElev:
                            if water_to_dist >0:
                                water_to_dist = water_to_dist - resmod.RiverOutlets[ro].Capacity
                                rivDict[ro] = 1
                        
                            else:
                                rivDict[ro] = 0   
                    
                else:
                    bypassFrac = 0.
                    
            # check if it's been long enough between gate changes
            gate_change_ready = False
            if resmod.SimulationSpecs.DaysSinceLastGateChange >= resmod.SimulationSpecs.GateChgFreq:
                gate_change_ready = True

            need_gate_change = False
            too_warm = prev_targ_sim > (this_ttarg + this_pos_temp_tol)
            too_cool = prev_targ_sim - (this_ttarg - this_neg_temp_tol) < 0
            if (too_warm or too_cool) and gate_change_ready:
                #print(f"check gate ops determine it was TooWarm: {too_warm} and TooCool: {too_cool}")
                need_gate_change = True
                
            if resmod.SimulationSpecs.DaysSinceLastGateChange>6:
                need_gate_change = True

            # get available gate levels and gates that can be opened
            gate_levs = outflow.gate_level_opts(resmod)
            # print(gate_levs)
            gate_options = outflow.tcd_gate_open_opts(resmod, gate_levs)
            
            # add a check for elevation of gate levels - can't use previous
            # gates if level drops below min head threshold
            prev_gate_levs = [k for k,v in resmod.PrevGateDict.items() if v>0]
            if len(set(prev_gate_levs) & set(gate_levs)) >0: # <-- checks if using previous gate levels is allowed under current conditions (is there overlap in the sets)
                pass
            else:
                need_gate_change = True
                gate_change_ready = True

            if this_ttarg == 99 or this_ttarg==99*1.8+32:
                if gate_options == [] and gate_levs==[0]: #  use lowest gate outlet
                    gate_dict = {k: 0 for k,v in resmod.Outlets.items()}
                    gate_dict[0] = resmod.Outlets[0].NumGates
                else:
                    gate_dict = gate_options[0]
                if check_gate_calcs:
                    return(['spec', gate_dict, targ_data, tdiff, rivDict, bypassFrac])
            
            # elif (increasing_errors2) and gate_change_ready:
            #     if check_gate_calcs:
                    
            #         return(['full', gate_options, targ_data, tdiff, rivDict, bypassFrac])
            
            elif gate_change_ready and need_gate_change:
                # try incremental gate search
                if check_gate_calcs:
                    if increasing_errors2 or increasing_errors:
                        if resmod.TimeStepDate.month<10: #was 9
                            return(['full', gate_options, targ_data, tdiff, rivDict, bypassFrac])
                        else:
                            return(['incr', gate_options, targ_data, tdiff, rivDict, bypassFrac])
                    else:
                        return(['incr', gate_options, targ_data, tdiff, rivDict, bypassFrac])

            elif (resmod.TimeStep == 0): 
                    #(increasing_errors and increasing_errors2):
                # (abs(tdiff) > 1*this_pos_temp_tol) or \

                # first time step or trending to higher temperature errors
                if check_gate_calcs:
                    
                    return(['full', gate_options, targ_data, tdiff, rivDict, bypassFrac])

            else:
                # do the previous steps gate options
                if check_gate_calcs:
                    if sum(rivDict.values())==0 and \
                        sum(resmod.Operations.GateOps[prev_date].values())==0:     
                        return['full',gate_options, targ_data, tdiff, rivDict, bypassFrac]
                    
                    return(['prev', resmod.Operations.GateOps[prev_date],
                           targ_data, tdiff, rivDict, bypassFrac])

        elif resmod.SimulationSpecs.CalcGateOps and final:
            # expects specified gate_dict
            gateDict = resmod.GateDict #kwargs['gate_dict']
            thisOutFlowTot = outflowvec[resmod.Outflow.ColumnMap['outflow']]
            thisRivOut = outflowvec[resmod.Outflow.ColumnMap['rivOutFlow']]

            if 'riv_dict' in kwargs:
                rivDict = kwargs['riv_dict']
            else:
                rivDict = {k: int(
                    outflowvec[v]) for k, v in resmod.Outflow.ColumnMap['rivOutlets'].items()}

            if thisRivOut > 0.:
                bypassFrac = thisRivOut/thisOutFlowTot
            else:
                bypassFrac = 0.

            targ_data = kwargs['temp_target']

        else:
            thisOutFlowTot = outflowvec[resmod.Outflow.ColumnMap['outflow']]

            gateDict = {
                k: int(outflowvec[v]) for k, v in resmod.Outflow.ColumnMap['gates'].items()}
            
            resmod.GateDict = gateDict

            thisRivOut = outflowvec[resmod.Outflow.ColumnMap['rivOutFlow']]

            rivDict = {k: int(
                outflowvec[v]) for k, v in resmod.Outflow.ColumnMap['rivOutlets'].items()}

            if thisRivOut > 0.:
                bypassFrac = thisRivOut/thisOutFlowTot
            else:
                bypassFrac = 0.

            if resmod.LP_Opt_Blending:
                targ_col_name = resmod.Outflow.ColumnMap['tempTarg']
                # need to set tailwater target temp
                target_temp_tw = resmod.Outflow.DataFrame.loc[resmod.TimeStepDate,
                                                              targ_col_name]
            elif 'tempTarg' in resmod.Outflow.ColumnMap:
                targ_col_name = resmod.Outflow.ColumnMap['tempTarg']
                # need to set tailwater target temp
                target_temp_tw = resmod.Outflow.DataFrame.loc[resmod.TimeStepDate,
                                                              targ_col_name]
            else:
                target_temp_tw = -9999
                
            thisTempTarg = target_temp_tw
            temptol = -9999
            targ_data = [thisTempTarg, temptol, temptol]

        if final:  # removal of water from reservoir and calculate release temp

            if 'gate_dict' in kwargs:
                gateDict = kwargs['gate_dict']
                resmod.GateDict = gateDict # overwrites any previous gate config


            if 'riv_dict' in kwargs:
                rivDict = kwargs['riv_dict']

            if 'bypass_frac' in kwargs:
                bypassFrac = kwargs['bypass_frac']

            outQLyrDist = outflow.selective_withdrawal(resmod, thisOutFlowTot, # removed gateDict arg - now an attr of resmod
                                                       bypassFrac, rivDict,
                                                       nPtSinks=3, logging=False,
                                                       dz=0)
            #print(f'finalizing temp')
            check_outflows = outQLyrDist.pop('Release_by_Outlet_Type')
            # print(check_outflows)
            if resmod.Debug['Release'] > 1:
                print(f"\n======================================")
                print(f"  Outflow distribution by layer: ")
                print(outQLyrDist, thisOutFlowTot)
                print(f"======================================\n")
                
            [totQ2, outT, outE, q_by_outlet] = outflow.outflow_check_and_dist(
                resmod, outQLyrDist, thisOutFlowTot)
            # outflowvec[resmod.Outflow.ColumnMap['outflow']])

            if abs(totQ2-thisOutFlowTot)>0.0001:
                raise BaseException(f"Calculated outflow {totQ2} doesn't match specified input {thisOutFlowTot}")

            resmod.Simulation_Results['ReleaseVol'].append(
                totQ2)  # 1) total outflow as calculated,
            resmod.Simulation_Results['ReleaseTemp'].append(
                outT)  # 2) mean outflow temp,
            resmod.Simulation_Results['ReleaseTempTarget'].append(
                targ_data[0])  # 3) temperature target (nan if not used)
            resmod.Simulation_Results['TempTargetTol'].append(
                [targ_data[0]+targ_data[1], targ_data[0]-targ_data[2]])  # 4) temp target tolerance (nan if not used)
            # 5) volume release by outlet level and type
            resmod.Simulation_Results['Release_by_Outlet'][resmod.TimeStepDate] = q_by_outlet

            resmod.distVol()  # distribute water after removing via selective withdrawal

            # apply diffusion
            resmod.diffuseDist()

            # update densities
            resmod.updateRho()

            # do convective mixing check - do less dense layers underly denser ones?
            stabcheck = resmod.checkDensityProfile()
            iterLimit = 1000
            iters = 0
            while (stabcheck > 0) & (iters < iterLimit):
                resmod.checkStable()
                stabcheck = resmod.checkDensityProfile()
                iterLimit += 1

            resmod.resVolume()
            resmod.getWSE()
            resmod.Simulation_Results['Storage'].append(resmod.Storage)
            resmod.Simulation_Results['Stage'].append(resmod.WSE)
            resmod.Simulation_Results['Profiles'].append(
                getTempProfile(resmod))
            resmod.TimeStepList.append(resmod.TimeStepDate)

            resmod.Operations.GateOps[resmod.TimeStepDate] = resmod.GateDict #gateDict
            resmod.PrevGateDict = resmod.GateDict

            if return_vals:
                return([totQ2, outT, outE])

        else:

            [totQ2, outT, outE] = outflow.calcOutTemp(resmod, thisOutFlowTot,
                                                      gateDict, rivDict, nPtSinks=3,
                                                      bypassFrac=bypassFrac)
            #[totQ2, outT, outE] = outflow.outflow_dist_forGateSelect(resmod, outQLyrDist)
            return([totQ2, outT, outE])

    def finalize(self):
        '''
        Organize results into dataframes for easier plotting, export to csv, etc

        Dataframes are stored to *DF keys within the restemp.Res.SimulationResults
        data object.

        Returns
        -------
        None

        '''
        if self.CalcMetric:
            stor_units = 'm3'
            temp_units = 'degC'
            flow_units = 'm3'
        else:
            stor_units = 'AF'
            temp_units = 'degF'
            flow_units = 'AF'

        sim_StoragesDF = pnd.DataFrame(self.Simulation_Results['Storage'],
                                       columns=[f'Sim_Storage_{stor_units}'],
                                       index=self.TimeStepList)

        sim_ReleasesDF = pnd.DataFrame([self.Simulation_Results['ReleaseVol'],
                                        self.Simulation_Results['ReleaseTemp'],
                                        self.Simulation_Results['ReleaseTempTarget'],
                                        self.Simulation_Results['TempTargetTol']]).T
        sim_ReleasesDF.index = self.TimeStepList

        sim_ReleasesDF.columns = [f'Sim_Release_{flow_units}',
                                  f'Sim_Release_Temp_{temp_units}',
                                  f'Sim_Temp_Target_{temp_units}',
                                  f'Sim_Target_Tolerance_{temp_units}']

        modElevs = [v.CtrElev for l, v in self.Layers.items()]
        sim_ProfilesDF = pnd.DataFrame(self.Simulation_Results['Profiles'],
                                       index=self.TimeStepList, columns=modElevs)
        self.Simulation_Results['ProfilesDF'] = sim_ProfilesDF
        self.Simulation_Results['ReleaseDF'] = sim_ReleasesDF
        self.Simulation_Results['StorageDF'] = sim_StoragesDF

        self.Simulation_Results['Release_by_Outlet_DF'] = \
            pnd.DataFrame.from_dict(self.Simulation_Results['Release_by_Outlet'],
                                    orient='index')

        if self.SimulationSpecs.CalcEvaporation:
            sim_evap = pnd.DataFrame(self.Evaporation, index=self.TimeStepList)
            self.Simulation_Results['EvapDF'] = sim_evap
        else:
            sim_evap = pnd.DataFrame(
                [np.nan]*len(self.TimeStepList), index=self.TimeStepList)
            self.Simulation_Results['EvapDF'] = sim_evap

    def reinitialize(resmod, return_obs_prof=False, show_init_profile=False, **kwargs):
        """
        Re-initialize the vertical temperature profile mid-simulation

        For use when doing year-by-year validation and for seasonal analysis;        

        Parameters
        ----------
        resmod : Res model object
            An instance of the ResTemp model that has been initialized with
            all required information - temperature profile, volume, and 
            boundary condition time series.
        return_obs_prof : bool, optional, default False
            Return a pandas.DataFrame to the console of the profile used to 
            reintialize the model. Useful for debugging.
        show_init_profile : bool, optional, default False
            Return a pandas.DataFrame to the console of the profile used to 
            intialize the model. Useful for debugging.
        **kwargs : 
            sim_dates : pandas.DataFrame.index
                A datetime index of simulation dates to replace existing values
                Useful if ResTemp model is coupled with another model that needs
                to reset or define simulation specifications.

        Returns
        -------
        None.

        """

        initDateStr = resmod.TimeStepDate.strftime('%Y-%m-%d')
        initDate_daybefore = (resmod.TimeStepDate-dt.timedelta(1))
        initDateStr_daybefore = initDate_daybefore.strftime('%Y-%m-%d')
        print(f"Reinitializing profile on {initDateStr}")

        if 'sim_dates' in kwargs:
            resmod.SimDates = kwargs['sim_dates']

        if resmod.SimulationSpecs.InitMethod == 'from_obs':
            outProfiles = resmod.Observations.ProfilesDF
            try:
                if resmod.Observations.ProfilesColMap['column'] == 'elevations':
                    initTempProf = outProfiles.iloc[outProfiles.index.get_loc(initDateStr,
                                                                              method='nearest')]
                else:
                    outProfiles = outProfiles.T
                    initTempProf = outProfiles.iloc[outProfiles.index.get_loc(initDateStr,
                                                                              method='nearest')]

            except Exception:
                print("assigning profile assuming 'tidy' format")
                try:
                    # try assigning profiles assuming a long/tidy format
                    # (indexes are date and elevation, values temperature)
                    
                    # commented code below required for earlier version of pandas
                    #idx1 = outProfiles.index.get_level_values(0).unique().get_loc(initDateStr,
                    #                                                              method='nearest')
                    
                    idx1 = outProfiles.index.get_level_values(0).unique().get_indexer([initDateStr], method='nearest')
                    
                    profDat = outProfiles.index.get_level_values(0).unique()[
                        idx1]

                    initTempProf = outProfiles.loc[profDat, :]

                except Exception:
                    print(
                        "Couldn't initialize from array-style or long-style profile tables")

            finally:

                # storage moved to the outflow category
                stoVar = resmod.Observations.OutflowColMap['obs_storage_final']

                if initDate_daybefore in resmod.Observations.OutflowDF.index:
                    initVol = resmod.Observations.OutflowDF.loc[initDateStr_daybefore,  stoVar]
                else:
                    initVol = resmod.Observations.OutflowDF.loc[initDateStr,  stoVar]

        else:
            print("Trying to initialize from data passed in from outside")

        util.initVolTemp(resmod, initVol, initTempProf,
                         profileTempUnits=resmod.SimulationSpecs.Profile_Temp_Units,
                         profileElevUnits='ft')
        resmod.updateRho()
        resmod.resVolume()

        if show_init_profile:
            printLayTempVolTotE(resmod)

        if return_obs_prof:
            return(initTempProf)

    def distVol(self, debug=False):
        """
        Distribute a change in volume through the reservoir layers.

        This function "fills in" water removed via releases and propagates
        extra water to overlying layers

        Parameters
        ----------
        debug : bool, optional
            Debug flag to indicate if extra information about the process should
            be printed to the console. The default is False.

        Returns
        -------
        None.

        """

        for l in self.Layers:  # n, ex in enumerate(excessVol):
            self.resVolume()
            # if debug:
            #    print("Layer: %d  Total Storage:  %0.1f" %(l, self.Storage))
            ex = self.Layers[l].Vol - self.Layers[l].MaxVol

            if abs(ex) > 0.:  # there is an imbalance that needs to be redistributed to/from other layers

                # this layer is paritally empty and water from an overlying layer should be brought down, closest first
                if (ex < 0.) & (l+1 < self.nLyrs):

                    volAbv = 0.
                    topLyr = 0
                    for lj in range(l+1, self.nLyrs):
                        if self.Layers[lj].Vol > 0:
                            topLyr = lj
                        volAbv += self.Layers[lj].Vol

                    # limit how much of the deficit can be addressed to how much water is available in layers above
                    tmpex = min(abs(ex), volAbv)

                    for li in range(l+1, topLyr+1):
                        if tmpex > 0:
                            # volume available in layer above
                            tmpvl = min(self.Layers[li].Vol, tmpex)
                            if tmpvl == 0:
                                #print("layer %s has no more water left" %li)
                                continue
                            tmpvlE = tmpvl*self.Layers[li].Temp

                            # add volume and temp/E to layer being filled
                            self.Layers[l].Vol += tmpvl
                            self.Layers[l].TotE += tmpvlE

                            # remove volume and E from layer being drained
                            self.Layers[li].Vol -= tmpvl
                            self.Layers[li].TotE -= tmpvlE

                            # update the temporary accounting variable
                            tmpex -= tmpvl

                            # update temperatures
                            if self.Layers[l].Vol > 0.:
                                self.Layers[l].Temp = self.Layers[l].TotE / \
                                    self.Layers[l].Vol
                            else:
                                self.Layers[l].Temp = np.nan

                            if self.Layers[li].Vol > 0.:
                                self.Layers[li].Temp = self.Layers[li].TotE / \
                                    self.Layers[li].Vol
                            else:
                                self.Layers[li].Temp = np.nan

                        if np.isnan(self.Layers[l].Temp):
                            print(
                                "Bringing water down to layer %s from layer %s gives temp of NaN" % (l, li))

                # this layer has excess water in it that should be distributed to overlying layers
                elif (ex > 0.) & (l+1 < self.nLyrs):

                    tmpvl = abs(ex)
                    tmpvlE = tmpvl*self.Layers[l].Temp

                    # remove volume from this layer
                    self.Layers[l].Vol -= tmpvl
                    self.Layers[l].TotE -= tmpvlE
                    self.Layers[l].Temp = self.Layers[l].TotE / \
                        self.Layers[l].Vol

                    if self.Layers[l+1].Vol == 0.:
                        self.Layers[l+1].TotE = 0.

                    self.Layers[l+1].Vol += tmpvl
                    self.Layers[l+1].TotE += tmpvlE
                    self.Layers[l+1].Temp = self.Layers[l +
                                                        1].TotE/self.Layers[l+1].Vol

                    if self.Layers[l].Temp == np.nan:
                        print("Moving water up to layer %s gives temp of NaN" % l)

                else:
                    # this is the top layer - if it's short, its short, if its overfull, there's spill
                    if ex > 0.:
                        print("spill!!! %0.2f" % ex)

    def distInflowMix(self, inQ, inT, log=True):
        """
        Distribute reservoir inflow to vertical layers according to density.

        Parameters
        ----------
        inQ : float
            Total reservoir inflow volume for the time step. Default units
            of acre-feet
        inT : float
            Temperature of the inflow volume. Default units deg F
        log : bool, optional
            Flag to indicate if details of the inflow mixing should be written
            to the log file. The default is True.

        Returns
        -------
        None.

        """
        tmpQ = inQ
        self.updateRho()  # update density in layers
        inTC = self.FtoC(inT)  # assuming inT is in F, convert to C
        inRho = self.calcRho(inTC)  # calculate density
        inE = inQ*inT
        
        if self.WSE < 940:
            mixfac = 1.1 #1.5
        else:
            mixfac = 1.0
        if log:
            lf = open(self.LogFilePath, 'a')
            lf.write('\n\nDate: %s\n' % self.TimeStepDate.strftime('%Y-%m-%d'))
            lf.write(
                '\nInflow T (deg F): %0.2f\nInflow Q (ac-ft): %0.2f\nInflow Density: %0.2f kg/m3' % (inT, inQ, inRho))
            lf.write('\n=====================================================\n')
            lf.write(
                'Lyr     LyrTemp    LyrVol_AF     WghtAvgT      InflowTempDiff   LayerTempDiff      TLprime      TIprime\n')
        for l in reversed(self.Layers.keys()):
            if tmpQ > 0:
                # inflow temperature equals or is higher than current layer, put water here
                if (inRho <= self.Layers[l].Rho) & (self.Layers[l].Vol > 0.):
                    #print("Inflow settled at layer %s" %l)
                    self.Layers[l].Vol = self.Layers[l].Vol + tmpQ
                    self.Layers[l].TotE = self.Layers[l].TotE + inE
                    self.Layers[l].Temp = self.Layers[l].TotE / \
                        self.Layers[l].Vol
                    tmpQ = 0.
                    if log:
                        ls = 'Layer density: %0.2f kg/m3\n' % self.Layers[l].Rho
                        ls += 'Inflow density: %0.2f kg/m3\n' % inRho
                        ls += 'Inflow mixing directly with layer: %d\n' % (l)
                        lf.write(ls)
                elif l == 0:  # made it to the bottom layer - if the water is this dense, put it here
                    self.Layers[l].Vol = self.Layers[l].Vol + tmpQ
                    self.Layers[l].TotE = self.Layers[l].TotE + inE
                    self.Layers[l].Temp = self.Layers[l].TotE / \
                        self.Layers[l].Vol
                    tmpQ = 0.
                else:
                    #print("inflow mixing with layer %s" %l)
                    # water is denser than current layer, heat exchange with layer before sinking
                    if self.Layers[l].Vol > 0.:
                        wavgT = (
                            self.Layers[l].Vol*self.Layers[l].Temp + inE)/(self.Layers[l].Vol + tmpQ)
                        TLprime = self.Layers[l].Temp + \
                            min(1, mixfac*self.C4)*(wavgT - self.Layers[l].Temp)
                        TIprime = inT + min(1, mixfac*self.C4)*(wavgT - inT)
                        if log:
                            ls = '%d      %0.3f       %0.3f    %0.3f       %0.3f       %0.3f       %0.3f      %0.3f\n' % (
                                l, self.Layers[l].Temp, self.Layers[l].Vol, wavgT, (wavgT-inT), (wavgT - self.Layers[l].Temp), TLprime, TIprime)
                            lf.write(ls)
                        #print("Inflow temp change: %0.1f degF" %(wavgT - inT))

                        # update layer temperature
                        self.Layers[l].Temp = TLprime
                        self.Layers[l].TotE = self.Layers[l].Temp * \
                            self.Layers[l].Vol

                        # update inflow temperature & density
                        inT = TIprime
                        inE = inQ*inT
                        inTC = self.FtoC(inT)
                        inRho = self.calcRho(inTC)
        if log:
            lf.close()

    def diffuse(self):
        """
        Implement diffusion mixing calculations based on vertical temperature
        profile and layer adjacency alone.

        This function replaced by `getDiffNeighbors` and `diffuseDist`

        Returns
        -------
        None.

        """
        updT = {}
        for l, v in self.Layers.items():
            if v.Vol == 0:  # no water in layer - nothing to diffuse
                updT[l] = np.nan
                continue
            if l == 0:  # bottom layer
                avT = (self.Layers[l+1].Vol*self.Layers[l+1].Temp +
                       v.Vol*v.Temp)/(self.Layers[l+1].Vol + v.Vol)
                updT[l] = v.Temp + self.C5*(avT - v.Temp)
            # top layer
            elif (l == self.nLyrs-1) or (self.Layers[l+1].Vol == 0.):
                avT = (self.Layers[l-1].Vol*self.Layers[l-1].Temp +
                       v.Vol*v.Temp)/(self.Layers[l-1].Vol + v.Vol)
                updT[l] = v.Temp + self.C5*(avT - v.Temp)
            else:
                avT = (self.Layers[l+1].Vol*self.Layers[l+1].Temp + v.Vol*v.Temp + self.Layers[l -
                       1].Vol*self.Layers[l-1].Temp)/(self.Layers[l+1].Vol + v.Vol + self.Layers[l-1].Vol)
                updT[l] = v.Temp + self.C5*(avT - v.Temp)

        # now update temperatures in each layer
        for l, v in self.Layers.items():
            if v.Vol > 0.:
                v.Temp = updT[l]
                v.TempC = self.FtoC(updT[l])
                v.Rho = self.calcRho(v.TempC)
                v.TotE = v.Temp*v.Vol

    def getDiffNeighbors(self):
        """
        Determine neighboring layers used in diffusion calculations.

        Determine the number of neighboring layers to be used in calculating 
        temperature differences driving diffusion. The global parameter 
        `DIFFUSE_DIST` and layer thicknesses determine the number of neighbors
        a particular layer has.

        Returns
        -------
        None.

        """
        for l, v in self.Layers.items():
            diffElevUp = v.CtrElev + self.DIFFUSE_DIST
            diffElevDn = v.CtrElev - self.DIFFUSE_DIST
            v.DiffLyrs = []
            if l == 0:  # bottom layer, only diffuse neighbors upwards
                for l2, v2 in self.Layers.items():
                    if v2.CtrElev < diffElevUp:
                        v.DiffLyrs.append(l2)
            elif l == self.nLyrs:  # top layer only diffuse neighbors down
                for l2, v2 in self.Layers.items():
                    if v2.CtrElev > diffElevDn:
                        v.DiffLyrs.append(l2)
            else:
                for l2, v2 in self.Layers.items():
                    if (v2.CtrElev > diffElevDn) & (v2.CtrElev < diffElevUp):
                        v.DiffLyrs.append(l2)

    def diffuseDist(self):
        """
        Calculate the diffusion across neighboring layers as determined by `getDiffNeighbors`

        Returns
        -------
        None.

        """
        updT = {}
        diffac = 1.
        if self.WSE < 940.0:
            diffac = 5.0 #7.5 
            
        for l, v in self.Layers.items():
            if v.Vol == 0:  # no water in layer - nothing to diffuse
                updT[l] = np.nan
                continue
            else:
                tmpVol = 0.
                tmpVolTemp = 0.
                for l2 in v.DiffLyrs:
                    if self.Layers[l2].Vol > 0.:
                        tmpVol += self.Layers[l2].Vol
                        tmpVolTemp += self.Layers[l2].Vol * \
                            self.Layers[l2].Temp
                avT = tmpVolTemp/tmpVol
                updT[l] = v.Temp + min(1, diffac*self.C5)*(avT - v.Temp)

        # now update temperatures in each layer
        for l, v in self.Layers.items():
            if v.Vol > 0.:
                v.Temp = updT[l]
                v.TempC = self.FtoC(updT[l])
                v.Rho = self.calcRho(v.TempC)
                v.TotE = v.Temp * v.Vol
            else:
                v.Temp = np.nan
                v.TempC = np.nan
                v.Rho = np.nan
                v.TotE = 0.

    def resVolume(self):
        """
        Update the current total reservoir storage variable

        Returns
        -------
        None.

        """
        totVol = 0.
        for l, v in self.Layers.items():
            totVol += v.Vol
        self.Storage = totVol

    def checkDensityProfile(self):
        """
        Iterate through layers to identify density instabilities.

        An instability will exist where a denser layer overlies a less dense
        layer. This will cause the Layer.Stable property to be marked as False
        and layers will be mixed by function `checkStable`.

        Returns
        -------
        None.

        """

        # go from top down - if a layer above is denser than one below,
        # set a value of 1, else 0; return the sum of 1's and 0's
        # if all = 0, then stable
        for l in reversed(self.Layers.keys()):
            self.Layers[l].Rho = round(self.Layers[l].Rho, DENSE_PREC)
        cnt = 0
        for l in reversed(self.Layers.keys()):
            if l == 0:  # at the bottom layer
                continue
            # layer above is denser than the one below
            elif (self.Layers[l].Rho - self.Layers[l-1].Rho) > DENSE_TOL:
                self.Layers.Stable = False
                cnt += 1
            else:
                # less dense layer above a denser layer
                self.Layers.Stable = True
        return(cnt)

    def checkStable(self):
        """
        Mix unstable layers identified by `checkDensityProfile`

        Where a denser layer overlies a less dense layer, the two are mixed.
        This mixing proceeds from the top down until now instabilities remain.

        Returns
        -------
        None.

        """
        # go from top down - if a layer above is denser than the one below,
        # mix completely
        for l in reversed(self.Layers.keys()):
            if l == 0:  # at the bottom layer
                continue
            # layer above is denser than the one below
            elif self.Layers[l].Rho > self.Layers[l-1].Rho:
                totVol = self.Layers[l].Vol + self.Layers[l-1].Vol
                totE = self.Layers[l].TotE + self.Layers[l-1].TotE
                totT = totE/totVol
                self.Layers[l].Temp = totT
                self.Layers[l-1].Temp = totT
                self.Layers[l].TotE = self.Layers[l].Temp*self.Layers[l].Vol
                self.Layers[l-1].TotE = self.Layers[l-1].Temp * \
                    self.Layers[l-1].Vol
                self.Layers[l].Rho = self.calcRho(
                    self.FtoC(self.Layers[l].Temp))
                self.Layers[l -
                            1].Rho = self.calcRho(self.FtoC(self.Layers[l-1].Temp))
            else:
                # less dense layer above a denser layer
                continue

    def assignResLayers(self, verbose=False):
        """
        Assign layers to each outlet to expedite withdrawal calculations.

        Given the layer elevations and outlet levels, assign layers to outlets
        they either intersect or completely cover.

        Parameters
        ----------
        verbose : bool, optional
            Flag to indicate whether the results of the layer-to-outlet 
            assignment should be printed to the console. The default is False.

        Returns
        -------
        None.

        """

        for outi, out in self.Outlets.items():
            outTopE = out.TopElevFt
            outBotE = out.BotElevFt
            outCtrEl = 0.5*(outTopE-outBotE) + outBotE  # out.CtrElevFt

            for l, v in self.Layers.items():
                if (v.MaxElev > outCtrEl) & (v.MinElev <= outCtrEl):
                    out.CtrLayer = l
                if (v.MaxElev > outBotE) & (v.MinElev <= outBotE):
                    out.MinLayer = l
                if (v.MaxElev > outTopE) & (v.MinElev <= outTopE):
                    out.MaxLayer = l
        if verbose:
            for ri, ro in self.Outlets.items():
                print("Outlet: %s" % ro.Name)
                print("Center layer assignment: %d" % ro.CtrLayer)
                print("Minimum layer assignment: %d" % ro.MinLayer)
                print("Maximum Layer assignment: %d" % ro.MaxLayer)
    #            print("--Assigned Layers:")
    #            for rl in ro.ResLayers:
    #                print("\t\t LayerIdx: %s" %rl)

        if verbose:
            print("\n********* LEAKAGES ******************")
        for outi, out in self.Leakages.items():
            outTopE = out.TopElevFt
            outBotE = out.BotElevFt
            outCtrEl = out.CtrElev

            for l, v in self.Layers.items():
                if (v.MaxElev > outCtrEl) & (v.MinElev <= outCtrEl):
                    out.CtrLayer = l
                if (v.MaxElev > outBotE) & (v.MinElev <= outBotE):
                    out.MinLayer = l

        if verbose:
            for ri, ro in self.Leakages.items():
                print("Leakage zone: %s" % ro.Zone)
                print("Center layer assignment: %d" % ro.CtrLayer)
                print("Minimum layer assignment: %d" % ro.MinLayer)

        if verbose:
            print("\n********* RIVER OUTLETS ******************")
        for outi, out in self.RiverOutlets.items():
            #outTopE = out.TopElevFt
            outBotE = out.MinElev
            outCtrEl = out.CtrElev

            for l, v in self.Layers.items():
                if (v.MaxElev > outCtrEl) & (v.MinElev <= outCtrEl):
                    out.CtrLayer = l
                if (v.MaxElev > outBotE) & (v.MinElev <= outBotE):
                    out.MinLayer = l
        if verbose:
            for ri, ro in self.RiverOutlets.items():
                print("River outlet level: %s" % ro.ID)
                print("Center layer assignment: %d" % ro.CtrLayer)
                print("Minimum layer assignment: %d" % ro.MinLayer)

    def updateRho(self):
        """
        Calculate density for all reservoir layers.

        Updates density for all layers in reservoir model.

        Returns
        -------
        None.

        """
        for l in self.Layers:
            self.Layers[l].TempC = self.FtoC(self.Layers[l].Temp)
            self.Layers[l].Rho = self.calcRho(self.Layers[l].TempC)

    def calcRho(self, tc):
        """
        Calculate density given a temperature in degrees Celsius

        Parameters
        ----------
        tc : float
            Temperature in degrees Celsius

        Returns
        -------
        rho : float
            Density in units kg/m^3

        """
        # given temp in C (tc), calulate density
        rho = 999.85050 + 0.06001*tc - 0.007917*tc**2 + 4.1256e-5*tc**3
        return(rho)

    def FtoC(self, tf):
        """
        Convert temperature in deg F to deg C

        Parameters
        ----------
        tf : float
            temperature in degrees Fahrenheit

        Returns
        -------
        tc : float
            Temperature in degrees Celsius

        """
        tc = (tf-32.)*5/9.
        return(tc)

    def CtoF(self, tc):
        """
        Convert temperature in deg C to deg F

        Parameters
        ----------
        tc : float
            Temperature in degrees Celsius

        Returns
        -------
        tf : float
            temperature in degrees Fahrenheit


        """
        # convert temp in F to C
        tf = tc*9/5. + 32.
        return(tf)

    def getWSE(self):
        """
        Calculate the water surface elevation given reservoir volume and an elevation-storage curve

        Returns
        -------
        None.

        """
        # find top layer
        topLyr = 0
        for l in self.Layers.keys():
            if self.Layers[l].Vol > 0:  # this is the top layer
                topLyr = l
        self.WSE = self.Layers[topLyr].volToElevInterp()
        self.TopLyr = topLyr
        

    def getLayerElevs(self, loc='ctr', getTemp=False):
        """
        Return the elevations of all layers

        Parameters
        ----------
        loc : str, optional
            Indicate which elevation is desired. Options are 'ctr' for center,
            'top' for top of layer, and 'bot' for bottom. The default is 'ctr'.
        getTemp : bool, optional
            Flag to indicate if the layer temperature should be returned
            with the elevation. The default is False.

        Returns
        -------
        lelev : list
            List of layer elevations, ordered from bottom to top
        ltemps : list, optional
            If `getTemp` is true, a list of layer temperatures

        """
        lelevs = []
        ltemps = []  # <- if we want temperatures too
        for l in self.Layers.keys():
            if getTemp:
                ltemps.append(self.Layers[l].Temp)
            if loc == 'top':
                lelevs.append(self.Layers[l].MaxElev)
            elif loc == 'bot':
                lelevs.append(self.Layers[l].MinElev)
            else:
                lelevs.append(
                    0.5*self.Layers[l].MinElev + 0.5*self.Layers[l].MaxElev)
        if getTemp:
            return([lelevs, ltemps])
        else:
            return(lelevs)


class lyr:
    """
    Class object defining reservoir layers and their properties
    """
    def __init__(self):
        self.Index = 0  # numbered from bottom
        self.Thick = 25
        self.MaxVol = 3.e5  # acre-feet
        self.MaxElev = 1067.  # ft
        self.MinElev = 1060.  # ft
        self.Vol = 1e5    # acre-feet
        self.TotE = 0.   # 'energy' - volume*temp
        self.Temp = 50.  # deg F
        self.ExcessVol = self.Vol - self.MaxVol
        self.TempC = 10.
        self.Rho = 1000.  # in kg/m3
        self.MaxArea = 0.  # in acres

    def updTotE(self):
        self.TotE = self.Temp*self.Vol

    def volToElevInterp(self):
        Elev = (self.Vol/self.MaxVol) * \
            (self.MaxElev - self.MinElev) + self.MinElev
        return(Elev)


class outlet:
    """
    Class object defining selective withdrawal outlets and their properties.
    """
    def __init__(self):
        self.Name = 'Upper'
        #self.Capacity = 9999.
        self.TopElevFt = 1067.
        self.BotElevFt = 1060.
        self.NumGates = 5
        self.ResLayers = []
        self.CtrElevFt = self.TopElevFt*0.5 + self.BotElevFt*0.5
        self.CtrLayer = 0
        self.MinLayer = 0
        self.MinHead = 0.
        self.MinHdUnits = 'ft'


class leakage:
    """
    Class object defining discrete leakage locations and properties that affect 
    selective withdrawal.
    """
    def __init__(self, topElevft=1067., botElevft=1000.):
        self.Zone = 1
        self.TopElevFt = topElevft
        self.BotElevFt = botElevft
        self.LeakageFactor = 0.0306
        self.ElevWeight = 1.
        self.CtrElev = 0.
        self.CtrLayer = 0
        self.MinLayer = 0

    def calcCtrElev(self):
        self.CtrElev = 0.5*self.TopElevFt + 0.5*self.BotElevFt


class riv_out:
    """
    Class object defining river outlets and their properties.
    """
    def __init__(self):
        self.ID = 1
        self.CtrElev = 937.
        self.MinElev = 937.
        self.DiamFt = 102./12.
        self.Capacity_CFS = 15000.
        self.CtrLayer = 0
        self.MinLayer = 0


def getTempProfile(res):
    prof = []
    for l, v in res.Layers.items():
        prof.append(v.Temp)
    return(prof)


def getColdWaterVol(res, ref_temp=52.):
    cwp = 0.
    for l, v in res.Layers.items():
        if v.Temp < ref_temp:
            cwp += v.Vol
    return(cwp)


def getTempCProfile(res):
    prof = []
    for l, v in res.Layers.items():
        prof.append(res.FtoC(v.Temp))
    return(prof)


def printLayVol(res):
    print("Layer ----  MaxVol ----  Vol\n========================")
    for l, v in res.Layers.items():
        print("  %d       %0.1f        %0.1f" % (l, v.MaxVol, v.Vol))


def printLayTemp(res):
    print("Layer ----  MaxVol ----  Vol\n========================")
    for l, v in res.Layers.items():
        print("  %d       %0.4f        %0.1f" % (l, v.Temp, v.Vol))


def printLayTempVolTotE(res):
    print("  Layer ---- Temp  ----- TotE ------ Vol ------ MaxVol \n=====================================")
    for l, v in res.Layers.items():
        print("  %d       %0.4f     %0.1f      %0.1f       %0.1f" %
              (l, v.Temp, v.TotE, v.Vol, v.MaxVol))


def interpProfile(obsProfile, simProfile):
    """
    Interpolate the simulated reservoir profile to the observations

    Parameters
    ----------
    obsProfile : pandas.DataFrame
        A dataframe indexed to elevation, with a single column being the 
        temperature at that elevation

    simProfile : pandas.DataFrame
        A dataframe indexed to elevation, with single column being the 
        simulated temperature at that elevation

    Returns
    -------
    intrpSimTemps : pandas.DataFrame
        A dataframe with the simulated temperatures interpolated to the 
        observation locations.

    """

    tmp = np.interp(obsProfile.index, simProfile.index, simProfile.values)

    if len(obsProfile.columns) > 1:  # <-- in this case, assuming it is matrix format for profile data
        intrpSimTemps = pnd.DataFrame(
            tmp, index=obsProfile.index, columns=[obsProfile.name])
    else:
        intrpSimTemps = pnd.DataFrame(
            tmp, index=obsProfile.index, columns=obsProfile.columns)
    return(intrpSimTemps)


def calc_sse(obs, sim):
    """
    Calculate sum of squared errors - assumes both obs and sim are in same units
    
    A convenience function for evaluating simulation outputs.

    Parameters
    ----------
    obs : pandas.DataFrame
        pandas.DataFrame of floats - observational values
    sim : pandas.DataFrame
        DataFrame of floats - simulated values

    Returns
    -------
    sse : float
        Sum of squared errors value

    """

    # flatten the df's to make sure there is just one column and one index
    if len(sim.shape) > 1:
        sim = sim.iloc[:, 0]
    if len(obs.shape) > 1:
        obs = obs.iloc[:, 0]
    sse = np.nansum((sim.values - obs.values)**2.)
    return(sse)
