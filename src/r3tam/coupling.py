# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 11:33:47 2022

@author: jgilbert
"""
import pandas as pnd
import numpy as np
import os, sys
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

class temperature_target:
    
    def __init__(self):
        
        self.Target_Model = None
        self.Target_Location = None
        self.Target_TS = inputTS()
        self.Sim_TS = {} # place to save simulated temperature at target location

    def __repr__(self):
        
        tstr = ("\nTemperature target details:")
        tstr += (f"\n\tTarget model:\t\t{self.Target_Model}")
        tstr += (f"\n\tTarget location:\t\t{self.Target_Location}")
        return(tstr)
    
class coupled_models:
    
    def __init__(self):
        
        self.Name = ''
        self.Description = ''
        self.Input_Config_FP = ''
        self.Component_Models = {}
        self.Connections = {}
        self.Simulation_Specs = simspec()
        self.Temperature_Target = temperature_target()
        
        
        
    def initialize(self, yaml_path, models=[]):
        
        self.Input_Config_FP = yaml_path
        inputs = self.read_input_yaml()
        
        startDateStr = inputs['Timing']['StartDate']
        endDateStr = inputs['Timing']['EndDate']
        delt_days = inputs['Timing']['TimeStepDays']
        self.SimDates = list(pnd.date_range(startDateStr, endDateStr, freq='d'))  #TODO: make this adjustable in case we want to do something different from daily (3-hourly or 2-daily, for exmaple)
        self.Simulation_Specs.StartDate_str = startDateStr
        self.Simulation_Specs.EndDate_str = endDateStr
        self.Simulation_Specs.DELT_DAY = delt_days 
        
        
        if 'Models' in inputs:
            for i,cm in enumerate(inputs['Models']):
                self.Component_Models[cm] = models[i]
                models[i].SimDates = self.SimDates # ensure that all models ahve the same simulaiton dates
                
        if 'Temperature_Target' in inputs:
            ttarg_dict = inputs['Temperature_Target']
            ttarg_mod_loc = ttarg_dict['Location'].split(":")
            self.Temperature_Target.Target_Location = ttarg_mod_loc[1]
            self.Temperature_Target.Target_Model = ttarg_mod_loc[0]
            
            if 'Units' in ttarg_dict:
                in_ttarg_units = ttarg_dict['Units']
            else:
                in_ttarg_units = 'degF'
                
            # we will convert the input target values into deg C
            
            targdefcol = 'Temp_Target_default_degC'
            targ_uppcol = 'Temp_Target_upper_tol_degC'
            targ_lowcol = 'Temp_Target_lower_tol_degC'
            ttarg_def_df = pnd.DataFrame(index=self.SimDates, columns=[targdefcol, 
                                        targ_uppcol,targ_lowcol],
                                         dtype=float)
                
            if 'Default' in ttarg_dict:
                # there is a default/fallback temperature target series, by month or day of year
                ttarg_def = ttarg_dict['Default']
                
                #print(ttarg_def)
                if 'd' in str(list(ttarg_def.keys())[0]):
                    # keys indicate day of year, not months
                    ttargdefts = []
                    for sd in self.SimDates:
                        ttargdefts.append(ttarg_def[sd.dayofyear])
                else:
                    ttargdefts = []
                    for sd in self.SimDates:
                        ttargdefts.append(ttarg_def[sd.month])
                
                if in_ttarg_units.upper() in ['DEGF','F','FAHR','FAHRENHEIT']:
                    ttargdefts = [(_-32)/1.8 if _<99 else 99. for _ in ttargdefts]
                    
                ttarg_def_df[targdefcol] = ttargdefts
                
                
            else:
                # if default temp target isnt' defined in the inputs,
                # the time series is filled in with 99's
                ttargdefts = [99.]*len(self.SimDates)
                ttarg_def_df[targdefcol] = ttargdefts
            
            ttarg_hdr = {'ttarg_default':targdefcol}
            
            if 'File' in ttarg_dict:
                ttarg_fp = ttarg_dict['File']['Path']
                if os.path.exists(ttarg_fp):
                    _ttarg = pnd.read_csv(ttarg_fp, index_col=0, parse_dates=True)
                    if len(_ttarg)!= len(self.SimDates):
                        raise BaseException("Length of specified temperature target does not match\n" +\
                                            f" that of the simulation time series. Check {ttarg_fp}")
                else:
                    raise FileNotFoundError(f"Temperature target time series file not found at {ttarg_fp}")
                
                hdr = ttarg_dict['File']['Header']
                ttarg_hdr['ttarg_specified'] = hdr['target']
                
            else:
                _ttarg = pnd.DataFrame(data=[np.nan]*len(self.SimDates),index=self.SimDates,
                                       columns=['Blank_Spec_Target'])
                ttarg_hdr['ttarg_specified'] = 'Blank_Spec_Target'
            
               

            if 'Upper_Tolerance' in ttarg_dict:
                
                ttarg_uptol = ttarg_dict['Upper_Tolerance']
                
                #print(ttarg_def)
                if 'd' in str(list(ttarg_def.keys())[0]):
                    # keys indicate day of year, not months
                    ttarguppts = []
                    for sd in self.SimDates:
                        ttarguppts.append(ttarg_uptol[sd.dayofyear])
                else:
                    ttarguppts = []
                    for sd in self.SimDates:
                        ttarguppts.append(ttarg_uptol[sd.month])
                if in_ttarg_units.upper() in ['DEGF','F','FAHR','FAHRENHEIT']:
                    ttarguppts = [(_)/1.8 if _<99 else 99. for _ in ttarguppts]
                    
                ttarg_def_df[targ_uppcol] = ttarguppts
                
            else:
                # if default temp target isnt' defined in the inputs,
                # the time series is filled in with 99's
                ttarguppts = [99.]*len(self.SimDates)
                ttarg_def_df[targ_uppcol] = ttarguppts
            
            ttarg_hdr['ttarg_upp_tol'] = targ_uppcol
            
            if 'Lower_Tolerance' in ttarg_dict:
                
                ttarg_lowtol = ttarg_dict['Lower_Tolerance']
                
                #print(ttarg_def)
                if 'd' in str(list(ttarg_def.keys())[0]):
                    # keys indicate day of year, not months
                    ttarglowts = []
                    for sd in self.SimDates:
                        ttarglowts.append(ttarg_lowtol[sd.dayofyear])
                else:
                    ttarglowts = []
                    for sd in self.SimDates:
                        ttarglowts.append(ttarg_lowtol[sd.month])
                
                if in_ttarg_units.upper() in ['DEGF','F','FAHR','FAHRENHEIT']:
                    ttarglowts = [_/1.8 if _<99. else 99. for _ in ttarglowts]
                    
                ttarg_def_df[targ_lowcol] = ttarglowts
                
            else:
                # if default temp target isnt' defined in the inputs,
                # the time series is filled in with 99's
                ttarglowts = [99.]*len(self.SimDates)
                ttarg_def_df[targ_lowcol] = ttarglowts
            
            ttarg_hdr['ttarg_low_tol'] = targ_lowcol
            
            
            ttarg_df = pnd.concat([ttarg_def_df, _ttarg], axis=1)
            
            self.Temperature_Target.Target_TS.ColumnMap = ttarg_hdr
            self.Temperature_Target.Target_TS.DataFrame = ttarg_df
            

        
        
        
    def read_input_yaml(self):
        import yaml
        
        with open(self.Input_Config_FP) as stream:
            inputs = yaml.load(stream, Loader=yaml.SafeLoader)    
            
        return(inputs)


def create_temp_target_shaping(target, shoulder, wlen, cntr_mon, cntr_day,
                               pre_ramp_period = 4, post_ramp_period = 2, 
                               **kwargs):
    
    if 'temp_units' in kwargs:
        temp_units = kwargs['temp_units']
    else:
        temp_units = 'degC'
    
    cust_dates = False
    if 'start_doy' in kwargs:
        start_doy=kwargs['start_doy']
        cust_dates=True
    else:
        start_doy = 1
        

        
    if 'end_doy' in kwargs:
        end_doy = kwargs['end_doy']
        cust_dates=True
    else:
        end_doy = 366
        
    if 'year' in kwargs:
        year = kwargs['year']
    else:
        year = 2021
    
    # a convenience function to ensure WY-based sims start on Oct 1
    on_wy = False
    if 'on_wy' in kwargs:
        on_wy = kwargs['on_wy']
        if on_wy:
            # assume year provided in kwargs is for WY
            start_doy = dt.datetime(year-1, 10, 1).timetuple().tm_yday
    
    if 'sim_dates' in kwargs:
        cust_dates=False
        
    if cust_dates:
        if start_doy > end_doy: # then starting in fall, going through winter/spring/summer of next year
            stdt = dt.date(year-1,1,1)+dt.timedelta(start_doy-1)
        else:
            stdt = dt.date(year,1,1)+dt.timedelta(start_doy-1)
        endt = dt.date(year, 1,1)+dt.timedelta(end_doy-1)
        simdates = pnd.date_range(stdt, endt, freq='d')
    else:
        simdates = kwargs['sim_dates']
        
    cntrDate = dt.date(year, cntr_mon, cntr_day)
    
    wdts= pnd.date_range(cntrDate - dt.timedelta(1)*((wlen/2.)*7-1),
                     cntrDate + dt.timedelta(1)*wlen/2.*7) 
    
    #if wdts[-1]> simdates[-1]:
    wdts = [w for w in wdts if w in simdates]
    
    pre_wdts = pnd.date_range(wdts[0]- dt.timedelta(1)*(7*pre_ramp_period-1),wdts[0])
    post_wdts = pnd.date_range(wdts[-1],wdts[-1]+ dt.timedelta(1)*(7*post_ramp_period-1))
    post_wdts = [pw for pw in post_wdts if pw in simdates]
    
    pre_vals = [shoulder-((shoulder-target)/len(pre_wdts))*(i+1) for i in 
                range(len(pre_wdts))]
    post_vals = [target+((shoulder-target)/len(post_wdts))*(i+1) for i in 
                 range(len(post_wdts))]

    if temp_units=='degC':
        targs = pnd.DataFrame([shoulder]*len(simdates), index=simdates,
                               columns=['target_shape_degC'])
        targs.loc[wdts,'target_shape_degC'] = target
        targs.loc[pre_wdts,'target_shape_degC']=pre_vals
        targs.loc[post_wdts,'target_shape_degC'] = post_vals
    else:
        targs = pnd.DataFrame([shoulder]*len(simdates), index=simdates,
                              columns=['target_shape_degF'])
        targs.loc[wdts,'target_shape_degF'] = target
        targs.loc[pre_wdts,'target_shape_degF']=pre_vals
        targs.loc[post_wdts,'target_shape_degF'] = post_vals
            
    return(targs)


def calc_gate_ops(resmod, cpl):
    
    ttarg_dict = cpl['Temperature_Target']
    ttarg_loc = ttarg_dict['Location']
    
    ttarg_default = ttarg_dict['Default']
    
    
       # if calcGateOps:

       #     if defaultTTarg:
       #         thisTempTarg = targtemp[dv.name.month]
       #     else:
       #         thisTempTarg = dv[targVar]
       #         if 'degC' in targVar:
       #             thisTempTarg = dv[targVar]*1.8+32.  # <-- convert to degF
       #         if dv.name.month in [12, 1]: #, 2]:
       #             thisTempTarg = 99. #<-- disable target temp for winter months
           
           
       #     thisPosTempTol = temptol[dv.name.month]
       #     thisNegTempTol = low_temptol[dv.name.month]
       #     # is the downstream temperature currently warmer or colder than the target?
       #     #dsTerr = ds_ccr_temps[-1] - thisTempTarg  # positive val = warmer than target
           
       #     need_gate_change = False
       #     if ds_ccr_temps[-1][0] > (thisTempTarg + thisPosTempTol):
       #         too_warm = True
       #         too_cool = False
       #         need_gate_change=True
       #     elif ds_ccr_temps[-1][0] -(thisTempTarg - thisNegTempTol) < 0: 
       #         too_cool = True
       #         too_warm = False
       #         need_gate_change=True
       #     else:
       #         too_warm = False
       #         too_cool = False
           

       #     # has it been long enough since our last gate change to be sure
       #     # the change has taken effect?
       #     gate_change_ready=False
       #     if days_since_gate_change >= sha.SimulationSpecs.GateChgFreq:
       #         gate_change_ready = True

       #     if ts==1:
       #         gate_options = sha.gate_open_opts(sha.gate_level_opts())
       #         gateDict = sha.gate_select3(thisTempTarg, 
       #                                     dv[outflowVar], gate_options,
       #                                     tolerance=temptol[dv.name.month],
       #                                     lower_temp_tolerance = low_temptol[dv.name.month],
       #                                     debugRelease=debugRelease,
       #                                     targetLoc='CCR', kwk_data=kwk_data, 
       #                                     kwk_model = kwk, ccr_temp=ds_ccr_temps)

       #     elif gate_change_ready and need_gate_change:
       #         gate_options = sha.gate_open_opts(sha.gate_level_opts())

       #         # if it's too warm, does an incremental change to gates help?
       #         prevGate = sha.Operations.GateOps[sha.TimeStepDate - dt.timedelta(1)]
       #         gateDict, newTemp = sha.gate_open_opts_incr(sha.gate_level_opts(),prevGate, 
       #                                            dv[outflowVar], thisTempTarg,
       #                                            temptol[dv.name.month],
       #                                            low_temptol[dv.name.month],
       #                                            debugRelease=debugRelease,
       #                                            targetLoc='CCR', kwk_data=kwk_data,
       #                                            kwk_model = kwk)                
       #         do_full_gate_search=False
       #         if too_warm: # if we were too warm before..
       #             if newTemp - (thisTempTarg + thisPosTempTol) > 0: # and we're still too warm
       #                 do_full_gate_search=True

       #         if too_cool: # if we were too cool before
       #             if ds_ccr_temps[-1][0] - newTemp < 1.5:
       #                 if (newTemp -(thisTempTarg - thisNegTempTol) < 0) & (thisTempTarg !=99):
       #                     do_full_gate_search=True
       #             else:
       #                 gateDict = prevGate

       #         if do_full_gate_search:
       #             try:
       #                 gateDict = sha.gate_select3(thisTempTarg, 
       #                                             dv[outflowVar], gate_options,
       #                                             tolerance=temptol[dv.name.month],
       #                                             lower_temp_tolerance = low_temptol[dv.name.month],
       #                                             debugRelease=debugRelease,
       #                                             targetLoc='CCR', kwk_data=kwk_data, 
       #                                             kwk_model = kwk, ccr_temp=ds_ccr_temps)
       #             except:
       #                 print("Gate search failed - reverting to previous configuration!")
       #                 gateDict = sha.Operations.GateOps[sha.TimeStepDate - dt.timedelta(1)]

       #         days_since_gate_change=0
       #     else:
       #         gateDict = sha.Operations.GateOps[dv.name-dt.timedelta(1)]  
       #         days_since_gate_change +=1

       # else:
       #     gateDict = {k:allgateDict[k][di] for k in allgateDict.keys()}
       
       # prev_gate_lev_opts = sha.gate_level_opts()
       
       # sha.Operations.GateOps[dv.name] = gateDict