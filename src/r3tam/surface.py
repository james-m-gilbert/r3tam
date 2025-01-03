# -*- coding: utf-8 -*-
"""
Created on Thu Feb 17 17:43:23 2022

@author: jgilbert
"""

from r3tam.constants import ACREtoM2, FTtoM, AFtoM3, FT3toM3, CFStoKG, Le, \
    cw, CMStoCFS, CFStoAFD, BOWEN_CONSTANT, TWOPI, HALFPI

import numpy as np

def distDepth(resmod): 
    """
    Determine transfer factors or fractions across the critical depth.

    For a  given critical depth, calculate the fractions to be applied to each 
    layer assuming a linear decrease from the top surface to the layer in which 
    the critical depth exists. Updates the factors and fractions stored internally
    in the restemp.Res object and available to other functions.

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
    critDepth = resmod.CriticalDepth* resmod.CriticalDepthFactor
    topLyr = 0
    for l in resmod.Layers.keys():
        if resmod.Layers[l].Vol > 0:  # this is the top layer
            topLyr = l
    
    resmod.getWSE()
    topLyrDepth = resmod.WSE - resmod.Layers[topLyr].MinElev
    #topLyrDepth = resmod.Layers[topLyr].Thick *(resmod.Layers[topLyr].Vol/resmod.Layers[topLyr].MaxVol)

    cumDepth = topLyrDepth
    lyrs = []
    #lyrs.append([topLyr, 1-(topLyrDepth/2./critDepth)])
    depfrac_temp = []
    depfrac_temp.append(1-(topLyrDepth/critDepth))
    n=1
    while cumDepth <= critDepth:
        l = topLyr-n
        thisThick = resmod.Layers[l].Thick
        #calcDepth = min(cumDepth+thisThick/2., critDepth-0.0001)
        calcDepth = min(cumDepth+thisThick, critDepth)
        if (cumDepth < critDepth) & (calcDepth>=critDepth): # critDepth is in this layer - cumDepth is still the bottom of the overlying layer until it's updated, below
            # then extrapolate the last portion of the linear trend (roughly)
            #print("at bottom of critical depth")
            depftmp = (critDepth-cumDepth)/critDepth * depfrac_temp[n-1]
        else:
            depftmp = 1-(calcDepth/critDepth)
        depfrac_temp.append(depftmp)
        cumDepth += thisThick
        
        #depFrac = 1 - (calcDepth/critDepth)
        #lyrs.append([l, depFrac])
        n+=1
    
    # now rescale the depth fractions so they sum to 1
    depFracs = [x/sum(depfrac_temp) for x in depfrac_temp]
    lyrs.append([topLyr, depFracs[0], depfrac_temp[0]])
    for l in range(1,n):
        lyrs.append([topLyr-l, depFracs[l], depfrac_temp[l]])
        
    resmod.CritDepthDist = lyrs
    
def wind_mixing(resmod,met_vec, tke_timeseries=True, DensDiffTol=0.001, debug=0):
    """
    Calculate and apply effect of surface wind on mixing through the water profile.
    
    The turbulent kinetic energy of the wind shear the water surface is compared
    against the resisting force of density stratification - if the TKE is greater
    then mixing of reservoir layers occurs.
    

    Parameters
    ----------
    resmod : TYPE
        DESCRIPTION.
    met_vec : dict
        Dictionary of meteorological values for this time step; keys are variable name,
        as defined in input yaml file under the meteorology time series section.
    DensDiffTol : float, optional
        Difference below which densities are considered equivalent for calculation
        purposes. The default is 0.001.
    debug : int, optional
        Set the debug level - higher values will produce more information about
        the simulation, more text output in the console window, and will slow 
        simulation. The default is 0.

    Returns
    -------
    None.

    """    
    topLyr = resmod.TopLyr
    rho_w = resmod.Layers[topLyr].Rho
    topElev = resmod.WSE    
    
    if tke_timeseries:
        
        tke_term = met_vec[resmod.Met.ColumnMap['tke']]
    
    else:
        ws = met_vec[resmod.Met.ColumnMap['wind2m']] #use the windspeed at 2m height
        
        rho_a = 1.2041	# kg/m3 - assuming air density at 20degC
        
        if isinstance(resmod.Cd, float):
            Cd = resmod.Cd  #Cd = 1.2e-3   # average/nominal drag coefficient - TODO: check if it is important to do Cd(wind speed), i.e. the way CE-QUALW2 does it 
        else:       # calculate according to formualtion in CE-QUAL-W2 v3.6 - p40 - Theory Manual -Wells & Cole 2021
            if ws < 0.5:  
                Cd = 0.01
            elif ws < 4.:
                Cd = 0.0044*ws**(-1.15)
            elif ws < 15.:
                Cd = 0.005*ws**(0.5)
            else:
                Cd = 0.0026 
                
        tau_aw = rho_a*Cd*ws**2
         
        tke_term = np.sqrt((tau_aw**3)/rho_w)*resmod.SimulationSpecs.DELT_SEC
        
               
        
    Area = resmod.Layers[topLyr].MaxArea
    Area_km2 = Area*ACREtoM2/1000./1000.   #TODO: Units check
    Wshltr = 1. - np.exp(-0.3*Area_km2)
    
    if debug>=1:
        print("\nWind sheltering coefficient: %0.4f\n" %Wshltr)
    
    #TKE_init = Wshltr*Area_km2*np.sqrt((tau_aw**3)/rho_w)*resmod.SimulationSpecs.DELT_SEC
    TKE_init = Wshltr*Area_km2*tke_term
    TKE = TKE_init
    
    nl = topLyr -1
    epiVol_m3 = resmod.Layers[topLyr].Vol*AFtoM3  #TODO: Units check
    epi_botElev_m =resmod.Layers[topLyr].MinElev*FTtoM  #TODO: Units check
    epiRho = resmod.Layers[topLyr].Rho
    
    epiThick_m = topElev*FTtoM - epi_botElev_m  #TODO: Units check
    epiMctr_m = epiThick_m/2. #(topElev+epi_botElev)/2.
  
    while TKE >0.:
        rhoDiff = resmod.Layers[nl].Rho - epiRho
        lyrVol_m3 = resmod.Layers[nl].Vol*FT3toM3 #TODO: Units check
        volratio = (epiVol_m3*lyrVol_m3)/(epiVol_m3 + lyrVol_m3)
        lyrElevDiff = epi_botElev_m - resmod.Layers[nl].CtrElev*FTtoM
    
        if debug>=2:
                print("Components of EPOT:")
                print("rhoDiff: %0.4f  volratio: %0.4f   epiThick_m: %0.4f    lyrElevDiff: %0.4f   epiMctr_m: %0.4f" %(rhoDiff, volratio, epiThick_m, lyrElevDiff, epiMctr_m ))
        
        EPOT = 9.81*rhoDiff*volratio*(epiThick_m + lyrElevDiff - epiMctr_m  )
        if debug>=2:
            print("Layer: %d TKE: %0.3f J  EPOT: %0.3f J" %(nl, TKE, EPOT))
            
        if (TKE > EPOT) or (rhoDiff < 0.):   # turbulent wind-mixing energy greater than potential, so mixing of epilimnion with layer below
            if debug>=1:
                print("Wind-driven mixing of epilimnion with layer %d" %nl)
            
            if nl > 2:
                TKE = max(0., TKE-max(EPOT, 0.))
            else:
                TKE = 0.
            
            # mix the layer below eplimnion with the bottom layer of the epilimnion
            totVol = resmod.Layers[nl].Vol + resmod.Layers[nl+1].Vol
            totE = resmod.Layers[nl].TotE + resmod.Layers[nl+1].TotE
            totT = totE/totVol
            resmod.Layers[nl].Temp = totT
            resmod.Layers[nl+1].Temp = totT
            resmod.Layers[nl].TotE = resmod.Layers[nl].Temp*resmod.Layers[nl].Vol
            resmod.Layers[nl+1].TotE = resmod.Layers[nl+1].Temp*resmod.Layers[nl+1].Vol
            resmod.Layers[nl].Rho = resmod.calcRho(resmod.FtoC(resmod.Layers[nl].Temp))
            resmod.Layers[nl+1].Rho = resmod.calcRho(resmod.FtoC(resmod.Layers[nl+1].Temp))
            
            epiVol_m3 += resmod.Layers[nl].Vol*FT3toM3
            epi_botElev_m = resmod.Layers[nl].MinElev*FTtoM
            epiRho = resmod.Layers[nl].Rho  # density of bottom newly-mixed layer
            
            epiThick_m = topElev*FTtoM - epi_botElev_m
            epiMctr_m = epiThick_m/2.
            
            nl = nl -1
        else:
            # if TKE isn't greater than EPOT in first layer comparison,
            # then no further mixing can happen
            TKE = 0.

def wind_func(resmod, met_vec, method='W2', TairUnits='degF'):
    """
    Calculate the wind function value for selected method.
    
    This function is used to calculate the value of the wind function, with the
    choice of two wind function formulations specified with the `method` 
    keyword.

    Parameters
    ----------
    resmod : Object
             ResTemp reservoir simulation object
        
    met_vec : pandas dataframe row - indexable by column name
        Listing of meteorological variable values for this time step
        
    method : string, optional
        Choice of wind function - Ryan-Harleman ('RH') or the CE-QUAL-W2 version ('W2').
        The default is 'W2'.
        
    TairUnits : string, optional
        Specify the units of the air temperature values provided. The default is 'degF'.

    Returns
    -------
    Value of the wind function calculated given the wind speed for this time step 
    and selected wind function type.
    """
    
    ws = met_vec[resmod.Met.ColumnMap['wind2m']]
    
    if method=='RH': # ryan-harleman formulation

        Tair = met_vec[resmod.Met.ColumnMap['airTemp_final']]
        if TairUnits=='degF':
            Tair = (Tair-32.)*5/9.
        Tsurf = (resmod.Layers[resmod.TopLyr].Temp-32.)*5/9.
        # Ryan-Harleman 1972 formulation, per CE-QUAL-W2 source
        tairv = (Tair+273.0)/(1.-0.378*resmod.Met.VapPress['ea_mmHg']/760.)
        dtv = (Tsurf+273.0)/(1.-0.378*resmod.Met.VapPress['es_mmHg']/760.)-tairv
        dtvl = 0.0084*ws**3  #<-- Lake Hefner f(w) - for comparison and fall-back at higher windspeed
        dtv = max(dtvl, dtv)
        fw = 3.59*dtv**0.3333 + 4.26*ws 
    
    elif method=='W2':
        # general formulation, per CE-QUAL-W2 manual
        if len(resmod.Wind_Func_Params)==0:
            awf = 8.5 #9.45  #9.2  # default, CE-QUAL-W2 manual, v3, eq 4-27
            bwf = 0.3 #0.46 # default
            cwf =  2. #2.05 #2  #default
        else:
            awf = resmod.Wind_Func_Params['awf']
            bwf = resmod.Wind_Func_Params['bwf']
            cwf = resmod.Wind_Func_Params['cwf']
            
        fw = awf + bwf * ws**cwf
    else:
        return(None)
    
    return(fw)
  
def calc_ea_es(resmod, met_vec, TsurfUnits='degF'):
    """
    Calculate actual and saturated vapor pressure from water and air temperatures.

    Parameters
    ----------
    resmod : Object
             ResTemp reservoir simulation object
        
    met_vec : pandas dataframe row - indexable by column name
        Listing of meteorological variable values for this time step
        
    TsurfUnits : string, optional
        Indicate the units of the water surface temperature. The default is 'degF'.

    Returns
    -------
    None.

    """
    Tdew = met_vec[resmod.Met.ColumnMap['dewpt']]
    if TsurfUnits == 'degF':
        Tsurf = (resmod.Layers[resmod.TopLyr].Temp-32.)*5/9.
    else:
        Tsurf = resmod.Layers[resmod.TopLyr].Temp
        
    ea_mmHg = np.exp(2.3026*((7.5*Tdew)/(Tdew+237.3)+0.6609))  # from CE-QUAL-W2 user's manual
    
    if Tsurf < 0.:  # if the surface happens to be frozen
        es_mmHg = np.exp(2.3026*((9.5*Tsurf)/(Tsurf+237.3)+0.6609))
    else:
        es_mmHg = np.exp(2.3026*((7.5*Tsurf)/(Tsurf+237.3)+0.6609))
    resmod.Met.VapPress = {'ea_mmHg': ea_mmHg, 'es_mmHg': es_mmHg}

def calc_evap(resmod, met_vec, TairUnits='degF', debug=0): 
    """
    Calculate evaporation from open water using meteorological variables and water temperature.
    
    Method adapted from ce-qual-w2 v4.2, heat-exchange-f90
        
    Parameters
    ----------
    resmod : restemp.Res
        ResTemp reservoir simulation object
    met_vec : pandas.DataFrame.row
        Row of dataframe containing data for the following variables:
            ws: wind speed, m/s
            Tdew: dew point temperture in deg c
            Tair: air temperature (ostensibly at 2 m)
            ea_mmHg : actual vapor pressure
            es_mmHg : saturated vapor pressure
    """
    
    ws =  met_vec[resmod.Met.ColumnMap['wind']]
    Tdew = met_vec[resmod.Met.ColumnMap['dewpt']]
    Tair = met_vec[resmod.Met.ColumnMap['airTemp']]
    
    # ea_mmHg calculated with call to calc_ea_es prior to calc_evap
    ea_mmHg = resmod.Met.VapPress['ea_mmHg']
    es_mmHg = resmod.Met.VapPress['es_mmHg']
    #ea_mmHg = np.exp(2.3026*((7.5*Tdew)/(Tdew+237.3)+0.6609))  # from CE-QUAL-W2 user's manual
    #ea_hPa = 6.11*10**((7.5*Tdew)/(237.3+Tdew))   # from: https://www.weather.gov/media/epz/wxcalc/vaporPressure.pdf

    topArea_m2 = resmod.Layers[resmod.TopLyr].MaxArea * ACREtoM2 #TODO: check units
    
    # do these calcs prior to call to calc_evap through         
    # Tsurf = (resmod.Layers[resmod.TopLyr].Temp-32.)*5/9.
    
    # if Tsurf < 0.:  # if the surface happens to be frozen
    #     es_mmHg = np.exp(2.3026*((9.5*Tsurf)/(Tsurf+237.3)+0.6609)) 
    # else:
    #     es_mmHg = np.exp(2.3026*((7.5*Tsurf)/(Tsurf+237.3)+0.6609)) 
    

    fw = wind_func(resmod, met_vec, method='W2')
    
    evap_Wm2 = max(0., fw*(es_mmHg - ea_mmHg))
    
    if evap_Wm2 >800.:
        print('Very high evaporation calculated!!!')
        print(f'   surface vapor pressure, mmHg: {es_mmHg: 0.4f}')
        print(f'   air vapor pressure, mmHg: {ea_mmHg: 0.4f}')
        #print(f'   tairv: {tairv: 0.4f}')
        #print(f'   dtv: {dtv: 0.4f}')
        #print(f'   dtvl: {dtvl: 0.4f}')
        print(f'   wind func: {fw: 0.4f}')
        print(f'   evap, W/m2: {evap_Wm2: 0.4f}')
        
    evapJ = evap_Wm2*topArea_m2*resmod.SimulationSpecs.DELT_SEC
    evap_vol_m3 = evap_Wm2 * resmod.SimulationSpecs.DELT_SEC/Le * 1/1000. * topArea_m2

    resmod.SurfHeat['Evap_Wm2'].append(evap_Wm2)
    resmod.SurfHeat['Evap_J'].append(evapJ)
    resmod.SurfHeat['Evap_m3'].append(evap_vol_m3)
    
def remove_evap(resmod, method='simple',debug=0):
    """
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
    
    evap_vol_AF = resmod.SurfHeat['Evap_m3'][resmod._time] / AFtoM3  # _time index minus one to acct for stored surfheat data being zero-indexed
    evapJ = resmod.SurfHeat['Evap_J'][resmod._time]
    topLyr = resmod.TopLyr
    
    #TODO: re-write this for metric, check units
    if evap_vol_AF > resmod.Layers[topLyr].Vol:
        remEvapAF = evap_vol_AF - resmod.Layers[topLyr].Vol
        resmod.Layers[topLyr].Vol = 0.
        resmod.Layers[topLyr].TotE = 0.
        resmod.Layers[topLyr].Temp = np.nan
        topLyr = topLyr -1
        resmod.Layers[topLyr].Vol -= remEvapAF
    else:
        resmod.Layers[topLyr].Vol -= evap_vol_AF
   
    # need to recalculate the distribution factors since the top elevation changed
    # likely not a big difference, but for completeness we recalcualte here
    resmod.getWSE()
    distDepth(resmod)

    
    #TODO: re-write this for metric, check units
    massTopLyr = resmod.Layers[topLyr].Vol * AFtoM3*resmod.Layers[topLyr].Rho #Revised from: 1000. # assumign density = 1000 kg/m3
    evapTempFDelta = 1.8*(evapJ/(massTopLyr*cw))
    if debug > 0:
        print("Evaporation effect on top layer: %0.1f degrees F" %evapTempFDelta)
    
    if method=='simple':
        # first version is for applying latent heat exchange only to top layer
        #evapTotE = resmod.C3 * resmod.Layers[topLyr].Vol * evapTempFDelta  # incremental Kelvin is equivalent to incremental Celsius
        #resmod.Layers[topLyr].TotE -= evapTotE  # assuming that latent heat is removed from the remaining volume in top layer
        #resmod.Layers[topLyr].Temp = resmod.Layers[topLyr].TotE/resmod.Layers[topLyr].Vol
        for lyr, depfac1, depfac2 in resmod.CritDepthDist:
            thisLayerMass = resmod.Layers[lyr].Vol*AFtoM3*resmod.Layers[lyr].Rho
            evapTempFDelta = 1.8*(evapJ/(thisLayerMass*cw))
            evapTotE = depfac1*resmod.C3 * resmod.Layers[lyr].Vol * evapTempFDelta  # TODO: change to Celsius - incremental Kelvin is equivalent to incremental Celsius
            resmod.Layers[lyr].TotE -= evapTotE  # assuming that latent heat is removed from the remaining volume in top layer
            resmod.Layers[lyr].Temp = resmod.Layers[lyr].TotE/resmod.Layers[lyr].Vol
    else:
        # under energy budget approach, heat is distributed through top layers
        # according to net of all energy terms
        raise Exception("Method other than 'simple' not yet implemented")  #pass
       
    resmod.Evaporation.append(evap_vol_AF)
    
def calc_latentheat_from_evap(resmod,met_vec, debug=False): 
    """
    Calculate latent heat for a given evaporation amount.

    Parameters
    ----------
    resmod : restemp.Res
        ResTemp reservoir simulation object
    met_vec : TYPE
        DESCRIPTION.
    debug : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """    
    #   from constants.py
    #    Le = 2.453e6   # latent heat of vaporization, J/kg, at 20deg C
    #    cw = 4182.  # specific heat of water at 20 deg C, in J/kg/K
    #   cfs2kg = 1/CMStoCFS * 86400.*1000. # convert from cfs to cubic meters/sec, seconds to days, to kg using assumed density of 1000. kg/m3
    evapCFS = met_vec[resmod.Met.ColumnMap['evap']]
    evapJ = evapCFS*CFStoKG *Le
    evapAF = evapCFS*CFStoAFD
    evap_vol_m3 = evapAF*AFtoM3
    
    topArea_m2 = resmod.Layers[resmod.TopLyr].MaxArea * ACREtoM2 #TODO: check units
    
    evap_Wm2 = evapJ/(topArea_m2*resmod.SimulationSpecs.DELT_SEC)
    
    resmod.SurfHeat['Evap_Wm2'].append(evap_Wm2)
    resmod.SurfHeat['Evap_J'].append(evapJ)
    resmod.SurfHeat['Evap_m3'].append(evap_vol_m3)

    
def air2water(resmod, met_vec): 
    """
    Function to exchange heat between water and the overlying air.

    Parameters
    ----------
    resmod : Object
             ResTemp reservoir simulation object
        
    met_vec : pandas dataframe row - indexable by column name
        Listing of meteorological variable values for this time step


    Returns
    -------
    ResTemp object with profile warmed by effects of heat transfer with overlying atmosphere 

    Notes
    -----
    This function approximates the heat exchange across the air-water interface 
    due to temperature differences between the two media, based on the formulation
    described in [[1]_]. 
    Uses the air temperature variable in ``metvec`` [units: :math:`^oC` or :math:`^oF`].
    Transfer of heat is proportional to the difference in air temperature and 
    water temperature at each layer within the "energy penetration depth" (or
    user-specified extinction depth, [units: m or ft]), with magnitude of heat
    transfer decreasing with depth from a maximum at the surface to zero at the
    critical depth.

    Following [[1]_], the formulation used is:
        
        .. math::
            
            E_{Air-Water} = \sum_{l=top}^{l=extinct} frac_l* C_1* (T_{air} - T_{water_l})*S_l
            
    Where
    
    :math:`E_{Air-Water} =` energy transfer at the air-water interface and mixing 
    zone, of dimension volume*temperature [units: acre-feet*°F or m3*°C]
    
    :math:`l =` layer index, where top denotes to the top (active) layer in the 
    reservoir and extinct denotes the layer intersected by the extinction depth
        
    :math:`frac_l =` Energy distribution fraction assigned to layer l. According 
    to the Reclamation temperature model formulation [1], this fraction “decreases 
    linearly with depth from 1 at the surface to 0 at the bottom [of the specified 
    depth of energy penetration]”
        
    :math:`C_1 =` User-specified calibration parameter that adjusts the degree of heat 
    transfer for air-water interaction processes (conduction, convection, 
    and, in the simplified approach, longwave radiation).
        
    :math:`T_{air} =` Air temperature for this time step, from user-provided input meteorology 
    time series [units: °F or °C]
        
    :math:`T_{water_l}=` Water temperature in layer *l* at beginning of time step [units: °F or °C]
        
    :math:`S_l=` Volume of water in layer *l* [units: acre-feet or :math:`m^3`]
            
    .. [1] Rowell, J. H. (1990). U.S. Bureau of Reclamation Monthly Temperature 
           Model Sacramento River Basin: Draft Report. U.S. Bureau of Reclamation, 
           Mid-Pacific Region.

    """
    
    
    # heat exchange from air to water (vice versa)
    # assume that the distribution factors over the critical depth have 
    # been calculated already
    
    airTemp = met_vec[resmod.Met.ColumnMap['airTemp']]
    
    for lyr, depfac1, depfac2 in resmod.CritDepthDist:
        lyrTemp = resmod.Layers[lyr].Temp
        lyrVol = resmod.Layers[lyr].Vol
        newE = depfac1*resmod.C1*(airTemp - lyrTemp)*lyrVol
        resmod.Layers[lyr].TotE += newE
        resmod.Layers[lyr].Temp = resmod.Layers[lyr].TotE/resmod.Layers[lyr].Vol

def simple_solrad(resmod, met_vec):
    """
    Calculate temperature change due to absorbed incident solar radiation.
    
    This function adapts the Reclaamtion Temperature Model formualtion for 
    warming surface water due to absorbed solar radiation across a surface
    interval using a linear rate of absorption with depth

    Parameters
    ----------
    resmod : restemp.Res
        ResTemp reservoir simulation object
        
    met_vec : pandas dataframe row - indexable by column name
        Listing of meteorological variable values for this time step

    Returns
    -------
    None.
    
    This function applies the solar radiation provided in the `metvec` [units: W/m^2]
    to the surface and for each underlying layer intersected by the `critical depth`
    defined as an input specification.

    See [[1]_] for description of the original formulation.

    .. [1] Rowell, J. H. (1990). U.S. Bureau of Reclamation Monthly Temperature 
               Model Sacramento River Basin: Draft Report. U.S. Bureau of Reclamation, 
               Mid-Pacific Region.

    """

    # UPDATE: modified to use radiation directly as W/m2 - to get the same
    # results as with solrad2, coefficient C2 will need to be smaller 
    # with this method
    # heat exchange from solar radiation
    # assumes distribution factors have been calculated previously
    # based on a prescribed critical depth
    # convert solar radiation in W/m2 into cal/m2/day
    
    radWm2 = met_vec[resmod.Met.ColumnMap['solRad']]
    
    if resmod.SeasonalRad or (resmod.SeasonalRad>0.):
        doy = resmod.TimeStepDate.dayofyear
        declination_angle = np.arcsin(0.39795 * np.cos((0.98563* (doy-173))*np.pi/180.))
        frac_rad = (180*declination_angle/np.pi/23.5+1.)/2.
        if resmod.SeasonalRad >0:
            frac_rad = frac_rad * (1. - resmod.SeasonalRad)+ resmod.SeasonalRad
    else:
        frac_rad = 1.
        
    for lyr, depfac, depfac2 in resmod.CritDepthDist:
        convfac = resmod.SimulationSpecs.DELT_SEC * 1/cw * 1/resmod.Layers[lyr].Rho
        Area_m2 = resmod.Layers[lyr].MaxArea*ACREtoM2
        E2 = depfac* frac_rad*resmod.C2 * radWm2 * convfac * Area_m2/1000. * resmod.SimulationSpecs.DELT_DAY * resmod.CriticalDepthFactor #  since radiation is energy per time unit    
        resmod.Layers[lyr].TotE += E2*1.8*1000./AFtoM3 # convert from thousand m3 * degC to ac-ft*degF #TODO - units check
        resmod.Layers[lyr].Temp = resmod.Layers[lyr].TotE/resmod.Layers[lyr].Vol
        

def solrad_w2(resmod, met_vec):
    """
    Calculate warming of surface waters due to solar radiation according to Beer-Lambert Law.

    Parameters
    ----------
    resmod : restemp.Res
        ResTemp reservoir simulation object
        
    met_vec : pandas dataframe row - indexable by column name
        Listing of meteorological variables for this time step

    Returns
    -------
    ResTemp object with profile warmed by effects of solar radiation 

    Notes
    -----
    This function applies the solar radiation provided in the ``metvec`` [units: :math:`W/m^2`]
    to the surface and for each underlying layer until solar radiation has been
    fully absorbed. 


    """
    
    
    radWm2 = met_vec[resmod.Met.ColumnMap['solRad']] # incident shortwave radiation for this time step, in units W/m2
    
    print(f"Incoming solrad: {radWm2} W/m2")
    albedo = resmod.MeanAlbedo
    if resmod.SeasonalAlbedo:
        if resmod.Latitude > 5.: # per GLM (Hipsey et al 2017 for the General Lake Model (GLM 2.4)), if latituate greater than 5 North, then northern hemisphere
            albedo = resmod.MeanAlbedo + resmod.AlbedoAmplitude*np.sin((TWOPI*resmod.DayOfYear/365.)+HALFPI)
        elif resmod.Latitude < -5: #similarly, check if in souther hemisphere
            albedo = resmod.MeanAlbedo + resmod.AlbedoAmplitude*np.sin((TWOPI*resmod.DayOfYear/365.)-HALFPI)
        else:
            pass # near the equator - assume albedo is constant at mean albedo
            
    # net shortwave entering the water is incident minus the reflected amount
    sw_reflect = radWm2*albedo 
    sw_water = radWm2 - sw_reflect
    print(f"Solrad entering water: {sw_water} W/m2")
    print(f"Solrad reflected: {sw_reflect}")
    # sediment (at lake bottom) reflectivity
    # NOTE(20220831): assumed for now that reflectivity is 1 and that all solar
    #                 radiation incident to sediment is reflected back and absorbed
    #                 by the layer through which it just passed
    sed_reflect = 0.
    
    # a certain amount of sw radiation is absorbed at the waters surface - 
    # controlled by user-specified fraction (sw_beta)
    # do the calc for the top layer first - then loop through lower layers until
    # the radiation making it through the layer is zero
    topThick_ft = (resmod.Layers[resmod.TopLyr].Vol/resmod.Layers[resmod.TopLyr].MaxVol)*resmod.Layers[resmod.TopLyr].Thick
    topThick_m = topThick_ft*FTtoM
    topSROUT = (1-resmod.C2)*sw_water*np.exp(-1*resmod.ExtinctionCoeff*topThick_m) 
    topAbs = sw_water - topSROUT  #amount absorbed in top layer, W/m2
    print(f"Solrad absorbed, top layer: {topAbs}")
    
    convfac =1/cw * 1/resmod.Layers[resmod.TopLyr].Rho
    Area_m2 = resmod.Layers[resmod.TopLyr].MaxArea*ACREtoM2 # 20220821jmg: right now this assumes light leaving the layer is available to the next layer, but in reality some would be absorbed or reflected from the lake bottom area represented by the difference in layer areas
    nextArea_m2 = resmod.Layers[resmod.TopLyr-1].MaxArea*ACREtoM2
    sedArea_m2 = Area_m2 - nextArea_m2 # area of sediment subject to incident solar radiation from layer above
    
    #calculate the solar radiation incident to and reflected back from the 
    # sediment at the bottom of the top layer and at the edge of the next layer down
    topSEDincid = topSROUT  # incident radiation flux (W/m2) to sediment is same as amount leaving top layer
    topSEDreflect = topSEDincid*sed_reflect  # amount reflected back from sediment is assumed to be absorbed into layer
    #topSEDabs = topSEDincid-topSEDreflect  #20220831jmg: commented out because not used
    
    print(f"Solrad absorbed from sediment, top layer: {topSEDreflect}")
    print(f"Sediment area, top layer: {sedArea_m2} m2")
    #20220831 DONE: add adjustment for layer-to-layer area difference ('flank') area reflection and absorption
    topE = (topAbs*Area_m2 + topSEDreflect*sedArea_m2)*resmod.SimulationSpecs.DELT_SEC*convfac # W/m2*area_m2-> W=> J/s * Sec/time step = J* 1/kg/m3 * 1/J/kg*C => m3*C (volume * temp change) 
    resmod.Layers[resmod.TopLyr].TotE += topE*1.8/AFtoM3 # convert from m3 * degC to ac-ft*degF #TODO - units check
    resmod.Layers[resmod.TopLyr].Temp = resmod.Layers[resmod.TopLyr].TotE/resmod.Layers[resmod.TopLyr].Vol
    
    print(f"Incremental temp change, top layer {topE*1.8/AFtoM3/resmod.Layers[resmod.TopLyr].Vol}")
    # now loop through underlying layers until all solar radiation is absorbed
    srouti = topSROUT 
    thisLayer = max(resmod.TopLyr-1,0)
    while (srouti-0.01) > 0. and thisLayer>=0:
        srin = srouti #incoming to this layer is outgoing of layer above
        thisThick_m = resmod.Layers[thisLayer].Thick*FTtoM
        newsrout = srin*np.exp(-1*resmod.ExtinctionCoeff*thisThick_m) # transmitted solar radiation flux, W/m2
        thisAbs = max(0.,srin - newsrout) # absorbed solar radiation flux, W/m2
        
        convfac =1/cw * 1/resmod.Layers[thisLayer].Rho
        Area_m2 = resmod.Layers[thisLayer].MaxArea*ACREtoM2
        nextArea_m2 = resmod.Layers[max(0, thisLayer-1)].MaxArea*ACREtoM2
        sedArea_m2 = Area_m2 - nextArea_m2 # area of sediment subject to incident solar radiation from layer above
        
        #calculate the solar radiation incident to and reflected back from the 
        # sediment at the bottom of this layer and at the edge of the next layer down
        thisSEDincid = newsrout  # incident radiation flux (W/m2) to sediment is same as amount leaving top layer - energy absorbed will be calculated based on sediment area over which this is applied
        thisSEDreflect = thisSEDincid*sed_reflect  # amount reflected back from sediment is assumed to be absorbed into layer
        #topSEDabs = topSEDincid-topSEDreflect  #20220831jmg: commented out because not used
       
        # 20220831 DONE: add adjustment for layer-to-layer area difference ('flank') area reflection and absorption
        topE = (thisAbs*Area_m2 + thisSEDreflect*sedArea_m2)*resmod.SimulationSpecs.DELT_SEC*convfac # W/m2*area_m2-> W=> J/s * Sec/time step = J* 1/kg/m3 * 1/J/kg*C => m3*C (volume * temp change) 
        resmod.Layers[thisLayer].TotE += topE*1.8/AFtoM3 # convert from m3 * degC to ac-ft*degF #TODO - units check
        resmod.Layers[thisLayer].Temp = resmod.Layers[thisLayer].TotE/resmod.Layers[thisLayer].Vol
        
        srouti = max(0,newsrout)
        thisLayer -= 1


        
def precip(resmod, met_vec):
    """
    Add precipitation falling on the reservoir surface to the top layer &
    update the top layer temperature to reflect heat change assuming
    precipitation is added at teh ambient air temperature

    Parameters
    ----------
    resmod : TYPE
        DESCRIPTION.
    met_vec : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # add precipitation to top layer at air temperature
    
    topLyr = resmod.TopLyr 
    airTemp = met_vec[resmod.Met.ColumnMap['airTemp']]
    precipVol = met_vec[resmod.Met.ColumnMap['precip']]
    resmod.Layers[topLyr].Vol += precipVol
    resmod.Layers[topLyr].TotE += precipVol*airTemp
    resmod.Layers[topLyr].Temp = resmod.Layers[topLyr].TotE/resmod.Layers[topLyr].Vol
    

def precipDepth(resmod, met_vec, Area=0.):
    
    # add precipitation to top layer at air temperature
    
    airTemp = met_vec[resmod.Met.ColumnMap['airTemp']]
    precipFT = met_vec[resmod.Met.ColumnMap['precip']]/12. # TODO: check units - it's assuming precip incoming as inches right now - change to meters
    
    topLyr= resmod.TopLyr
    
    if Area>0.:
        precipVol = precipFT*Area
    else:
        precipVol = precipFT*resmod.Layers[topLyr].MaxArea
        
    resmod.Layers[topLyr].Vol += precipVol
    resmod.Layers[topLyr].TotE += precipVol*airTemp
    resmod.Layers[topLyr].Temp = resmod.Layers[topLyr].TotE/resmod.Layers[topLyr].Vol
    
def calc_sensible(resmod, met_vec, fw):
    
    Tsurf = (resmod.Layers[resmod.TopLyr].Temp-32.)*5/9.  #TODO check units/convert to SI
    Tair = (met_vec[resmod.Met.ColumnMap['airTemp']]-32)/1.8
    Q_sens = fw * BOWEN_CONSTANT*(Tsurf - Tair)
    
    