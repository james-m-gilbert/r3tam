# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 18:06:29 2022

@author: jgilbert
"""
import os, sys

import matplotlib.pyplot as plt 

sys.path.insert(0, r'D:\02_Projects\SacTemp\SimTemp\R3TAM\src')
from r3tam import restemp as rt
from r3tam import longtemp as lt
from r3tam.river import River
from r3tam import coupling

from r3tam import make_plots as mkp

# import restemp as rt
# import longtemp as lt
# from river import River
# import coupling
# import make_plots as mkp

import numpy as np

import datetime as dt
import copy, os, time



# read in the config files for each model
config_fp01 = r'D:/02_Projects/SacTemp/SimTemp/R3TAM/examples/shasta/20250314_shasta_only_report_0022_withLPblending.yaml'
config_fp02 = r'D:/02_Projects/SacTemp/SimTemp/R3TAM/examples/keswick/kwk_input.v20240729.yaml'
config_fp03 = r'd:/02_Projects/SacTemp/SimTemp/R3TAM/examples/sac_river/uppsac_river_kwk2rdb.yaml'
couple_config_fp = r'd:/02_Projects/SacTemp/SimTemp/R3TAM/examples/shasta_keswick_sac/coupled_shakwkriv_dynblend.yaml'

#%% set up output directories
proj_dir = r'D:\02_Projects\SacTemp\SimTemp\R3TAM\examples\shasta_keswick_sac'
outdir = os.path.join(proj_dir, 'outputs','figs')

save_plots=False
if not os.path.exists(outdir):
    print(f"making output directory: {outdir}")
    os.mkdir(outdir)
stats_dir = os.path.join(proj_dir, 'outputs', 'stats')
if not os.path.exists(stats_dir):
    print(f"making output stats directory: {stats_dir}")
    os.mkdir(stats_dir)

todaystr = dt.date.today().strftime('%Y%m%d')

#%%  now try using gate selection logic to hit targets
#  at selected stream nodes (ccr)
sha = rt.Res.initialize_model(config_fp01, profile_temp_units='degF')
kwk = lt.LongTemp.initialize_longmod(config_fp02, showInit=True)
uppsac = River.initialize_model(config_fp03)

cpl = coupling.coupled_models()
cpl.initialize(couple_config_fp, models=[sha, kwk, uppsac])

sha.SimulationSpecs.IsHindcast = False #<-- switch so that calcGateOps is engaged
sha.SimulationSpecs.CalcGateOps = True
sha.SimulationSpecs.GateChgFreq = 2 #3
sha.Debug['Release'] = 0

btime_all = time.perf_counter()
runtime_dict = {}

all_profile_errors = {} # dictionary to hold profile error calcs
all_opengate_errors = {}
all_sha_rel_error = {}
all_kwk_rel_error = {}
all_ccr_error = {}
all_bsf_error = {}
all_bnd_error = {}
all_rdb_error = {}
side_gate_compare = {}
kwk_obsT = {} # keswick obs releae temperature time series, by year
kwk_simT = {} # keswick sim release temp by year
#%%
btime_all = time.perf_counter()
for yr in range(2015, 2016): # 2023): #[2017]: #, 2015, 2016]:
    
    btime_this = time.perf_counter()
    # reinitialize coupled model
    if yr==2021 or yr==2020: #2010:
        couple_config_fp = r'd:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/coupled_shakwkriv_2021.yaml'
    else:
        #couple_config_fp = r'd:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/coupled_shakwkriv.yaml
        couple_config_fp = couple_config_fp
    sha = rt.Res.initialize_model(config_fp01, profile_temp_units='degF')
    kwk = lt.LongTemp.initialize_longmod(config_fp02, showInit=True)
    uppsac = River.initialize_model(config_fp03)

    cpl = coupling.coupled_models()
    cpl.initialize(couple_config_fp, models=[sha, kwk, uppsac])

    sha.SimulationSpecs.IsHindcast = False #<-- switch so that calcGateOps is engaged
    sha.SimulationSpecs.CalcGateOps = True
    sha.SimulationSpecs.GateChgFreq = 3
    sha.Debug['Release'] = 0
        
    prev_gate_mode=''

    main_targ_temp_degC = 14
    should_temp_degC = (58-32)/1.8
    targts_shape = coupling.create_temp_target_shaping(main_targ_temp_degC,
                                                       should_temp_degC, 
                                                       24, 8, 1, temp_units='degC',
                                                       pre_ramp_period=2, post_ramp_period=2,
                                                       year=yr, start_doy=1, 
                                                       end_doy= dt.date(yr,12,31).timetuple().tm_yday) #end_doy)
    targts_shape.loc[f'{yr}-01-01':f'{yr}-04-15',:] = 99
    targts_shape.loc[f'{yr}-11-01':f'{yr}-12-31',:] = 99
    
    if yr==2022:
        targts_shape = targts_shape.iloc[:dt.date(yr, 9,30).timetuple().tm_yday]

    # create a target at CCR based on a smoothed version of the observations
    ccr_obs = uppsac.Nodes['ccr'].Observations.DataFrame.loc[f'{yr}']
    ccr_targ = ccr_obs.rolling(30, center=True).agg(np.median).bfill().ffill()
    ccr_targ_degC = (ccr_targ['CCR_degF']-32)/1.8
    ccr_targ_degC = targts_shape
    ccr_targ_degC.name = 'CCR_degC'
    
    # set the new target ts as the specified in cpl model
    cpl.Temperature_Target.Target_TS.DataFrame.loc[targts_shape.index,'Spec_Target_degC'] = ccr_targ_degC.values # targts_shape['target_shape_degC']
    cpl.Temperature_Target.Target_TS.ColumnMap['ttarg_specified'] = 'Spec_Target_degC'
    
    cpl.SimDates = targts_shape.index
    sha.TimeStep = -1
    sha.TimeStepDate = cpl.SimDates[0]
    sha.reinitialize(sim_dates=targts_shape.index)

    kwk.TimeStep = -1
    kwk._time = -1
    kwk.TimeStepDate = cpl.SimDates[0]
    kwk.reintialize(sim_dates=targts_shape.index)
    
    uppsac.TimeStep = -1
    uppsac.TimeStepDate = cpl.SimDates[0]
    uppsac.SimDates = targts_shape.index
    
    prev_targ = cpl.Temperature_Target.Target_TS.DataFrame.loc[cpl.SimDates[0],'Spec_Target_degC']

    gate_method = {}
    sha_sac_diff_dict = {}
    for d in cpl.SimDates: #[0:200]: #[0:154]: #123 = May 2]:
        
        print(d)
        
        sha.advance_restemp() #advances the time step counter for the sha model, does inflow and met calculations
        
        # do an initial check of gate options
        [gate_mode, gates, targ_data, prev_err, rivDict, bypass_frac] = sha.advance_swd(final=False, coupling_model=cpl, 
                                                                  check_gate_calcs=True)
        
        # if abs(targ_data[0] - prev_targ)>0.2:
        #     # do full gate search if the target changes a lot from one day to the next
        #     print("****Setting gate search to 'full' due to target change..")
        #     gate_mode = 'full'
        #     prev_targ = targ_data[0]
        #     # get available gate levels and gates that can be opened
        #     gate_levs = rt.outflow.gate_level_opts(sha)
        #     # print(gate_levs)
        #     gates = rt.outflow.tcd_gate_open_opts(sha, gate_levs)
        
        # if d.month==9 and prev_err > 0.25  and targ_data[0] == main_targ_temp_degC:
        #     # do full gate search if the target changes a lot from one day to the next
        #     print("****Setting gate search to 'full' due to warm water..")
        #     gate_mode = 'full'
        #     prev_targ = targ_data[0]
        #     # get available gate levels and gates that can be opened
        #     gate_levs = rt.outflow.gate_level_opts(sha)
        #     # print(gate_levs)
        #     gates = rt.outflow.tcd_gate_open_opts(sha, gate_levs)
        
        gate_method[d] = gate_mode
        #print(f"gate mode: {gate_mode} - previous gate mode: {prev_gate_mode}")
        #if gate_mode=='full' and prev_gate_mode=='full':
            #print(f"\n{d.strftime('%Y-%m-%d')}:  two full gate modes in a row")
            #print(f"  Days since last gate change: {sha.SimulationSpecs.DaysSinceLastGateChange}")
            #break
        thisOutFlowTot = sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['outflow']]
        thisRivOut = sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['rivOutFlow']]

        if thisRivOut>0.:
            bypassFrac = thisRivOut/thisOutFlowTot
        else:
            bypassFrac=0.
        
        prev_date = d-dt.timedelta(1)
        if prev_date in sha.Operations.GateOps:
            prev_gate_dict = sha.Operations.GateOps[prev_date]
        else:
            prev_gate_dict = {3:0, 2:0,1:0, 0:2}
        
        # do an initial assessment - are temperatures form last time step's gate
        # config warm/cold compared to target?
        # get shasta temps using previous gate config
        [q_o, t_o, e_o ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
                                                 rivDict,
                                                 bypassFrac=bypassFrac,
                                                 gateDictOpt=prev_gate_dict)
        #init_sha_sac_diff = (targ_data[0]*1.8+32)-t_o
        if len(uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp)==0:
            init_sha_sac_diff = 1
        else:
            init_sha_sac_diff = (uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]*1.8+32) - t_o
            
        # calc the running 5-day average of difference between shasta and target
        if len(uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp) < 6:
            running_sha_sac_diff = 2.0
        else:
            last5_sha_degF = sha.Simulation_Results['ReleaseTemp'][-18:-4]
            last5_sac_degC = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-14:]
            last5_sac_degF = [y*1.8+32 for y in last5_sac_degC]
            running_sha_sac_diff = np.mean([y-x for x,y in zip(last5_sha_degF, last5_sac_degF)])
        
        sha_sac_diff_dict[d] = running_sha_sac_diff
        # monthly adjustments for shasta target
        
        # make adjustments according to an interpolated trend that goes up starting
        # in may and back down starting in sept
        may1_adj = 2.25
        sep1_adj = 4.25
        sep15_adj = 4.5
        sep30_adj = 4.75
        nov1_adj = 1.5
        may1_sep1_days = 123
        may1_sep15_days = 138
        may1_sep30_days = 153
        nov1_sep1_days = 61
        nov1_sep15_days = 46
        nov1_sep30_days = 31
        
        # if d.month in [5,6,7,8] or (d.month==9 and d.day <=30):
        #     sha_targ_adj_1 = may1_adj + (sep30_adj - may1_adj)*((d - dt.datetime(d.year, 5,1)).days/may1_sep30_days)
        #     sha_targ_adj_ = (targ_data[0]*1.8+32) - sha_targ_adj_1
        # elif d.month==10: # or (d.month==9 and d.day>15): # [9,10]:
        #     sha_targ_adj_1 = sep30_adj + (nov1_adj - sep30_adj)*((d - dt.datetime(d.year, 9,30)).days/nov1_sep30_days)
        #     sha_targ_adj_ = (targ_data[0]*1.8+32) - sha_targ_adj_1
        # else:
        #     sha_targ_adj_ = (targ_data[0]*1.8+32)- running_sha_sac_diff
        
        if d.month in [5]:
            sha_targ_adj_ = -0.5
        elif d.month in [6,7]:
            sha_targ_adj_ = -0.75
        elif d.month in [8]:
            sha_targ_adj_ = -0.8 #1.0
        elif d.month in [9,10]:
            sha_targ_adj_ = -0.9 #-1.0 #25
        else:
            sha_targ_adj_ = 0.0
            
        prev_err_targ_adj_  = 0.0
        if prev_err > 0.1 and prev_err < 0.33:
            prev_err_targ_adj_ = prev_err
        prev_sha_targ_ = sha.Outflow.DataFrame.loc[d- dt.timedelta(1), sha.Outflow.ColumnMap['tempTarg']]
        
        this_sha_targ_ = (targ_data[0]*1.8+32)- running_sha_sac_diff + sha_targ_adj_
        #this_sha_targ_ = sha_targ_adj_
        
        if d.month in [9,10] and targ_data[0] == main_targ_temp_degC:
            set_sha_targ = np.nanmin([prev_sha_targ_, this_sha_targ_])
        elif d.month in [6,7,8] and targ_data[0] == main_targ_temp_degC:
            set_sha_targ = np.nanmin([prev_sha_targ_+0.1, this_sha_targ_])
        else:
            set_sha_targ = this_sha_targ_
        #sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = set_sha_targ - prev_err_targ_adj_ #
        sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = (targ_data[0]*1.8+32)- running_sha_sac_diff + sha_targ_adj_
        
        if gate_mode=='full':
            print("doing full gate selection...")
            # do full gate selection
            
            # turn off lp-opt for now (remember to turn it back on once gates are selected!)
            #sha.LP_Opt_Blending = True
            
            best_error = 99
            best_neg_error = 99
            best_gates = prev_gate_dict
            
            prev_kwk_intrmd = 0
            for gc in gates: # iterate through all gate options
                #init_sha_targ = (targ_data[0]*1.8+32)-init_sha_sac_diff
                #sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = init_sha_targ
                
                [q_o, t_o, e_o ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
                                                         rivDict,
                                                         bypassFrac=bypassFrac,
                                                         gateDictOpt=gc)
                #print(t_o)
                # kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
                # kwk.set_inflow_temps(kwk_in_dict,set_date=d,
                #                      temp_units='DEG_F', flow_units='AF')
                
                # kwk.advance_longtemp(final=False)
                
                # # if abs(kwk.IntermedTemp - prev_kwk_intrmd) < 0.001:
                # #     raise BaseException("kwk temp inconsistencies")
                # # prev_kwk_intrmd = kwk.IntermedTemp
                
                # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
                #                           kwk.IntermedRelease/86400, set_date=d)
                # uppsac.advance_rivmod(final=False)
                # this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
                
                # #sha_sac_tdiff = t_o - (this_temp*1.8+32.0)
                # #sha_adj_targ = (targ_data[0]*1.8+32)-sha_sac_tdiff
                
                # # current actual difference between shasta and target, 
                # # pos value = warmer river than sha; neg value = shas warmer than river
                # sha_sac_tdiff = (this_temp*1.8+32.0) - t_o
                # #difference between current temp at target loc and target temp
                # sac_targ_diff = (targ_data[0]-this_temp)*1.8
                
                # sha_adj_targ = (targ_data[0]*1.8+32)-sha_sac_tdiff
                # #sha_adj_targ = t_o + sac_targ_diff # targ_data[0]*1.8+32 - sha_sac_tdiff 
                
                # # update the temperature target at Shasta tailwater
                # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ
                
                # set adjusted shasta tailwater target
                # sha_tw_col = sha.Outflow.ColumnMap['tempTarg']
                # sha.Outflow.DataFrame.loc[d, sha_tw_col] = 12.0 # sha_adj_targ
                
                # # now re-run with the adjusted sha tw targ
                # [q_o2, t_o2, e_o2 ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
                #                                          rivDict,
                #                                          bypassFrac=bypassFrac,
                #                                          gateDictOpt=gc)
                
                # kwk_in_dict2 = {'inflow_final': q_o2, 'inflowTemp_final': t_o2}
                # kwk.set_inflow_temps(kwk_in_dict2, set_date=d,
                #                      temp_units='DEG_F', flow_units='AF')
                # kwk.advance_longtemp(final=False)
                # print(f"\n\t\tKeswick temp: {kwk.IntermedTemp}")
                
                # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
                #                           kwk.IntermedRelease/86400, set_date=d)
                # uppsac.advance_rivmod(final=False)
                
                # this_temp2 = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
                
                #this_tdiff = this_temp -targ_data[0]
                # 20250423 compare to shasta estimated target, not river
                this_sha_targ = sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']]
                this_tdiff = t_o - this_sha_targ
                
                # if this_tdiff < 0:
                #     if abs(this_tdiff) < abs(best_neg_error):
                #         print(f">>>> added {gc} to best gates for neg error")
                #         best_neg_error = this_tdiff
                #         best_gates = gc
                
                print(f"Gate config: {gc} - Diff: {this_tdiff}")
                if abs(this_tdiff) < abs(best_error):
                    print(f">>>> added {gc} to best gates")
                    best_error = this_tdiff
                    best_gates = gc
                    
                if sha.Debug['Release'] >0:
                    print("       Target: %0.2f - %0.2f - %0.2f deg F" %(targ_data[0]+targ_data[1], targ_data[0], targ_data[0]-targ_data[2]))
                    print("       ResultingTemp:  %0.2f deg F" %this_temp2)
                    
                #too_cold = this_temp2 < (targ_data[0] - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
                #too_warm = this_temp2 > (targ_data[0] + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate
                
                too_cold = t_o < (this_sha_targ - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
                too_warm = t_o > (this_sha_targ + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate
                
                
                gate_dict = gc
                #if this_tdiff < 0.:
                if (not too_warm) & (not too_cold):
                    print("found an option!")                                       
                    break
                
            try:
                if gate_dict != best_gates:
                    print("found an option from 'best_gates'..")
                    gate_dict = best_gates
                else:
                    print("using last gate option as best...")
            except:
                print("Gate search failed - repeating previous step's configuration")
                gate_dict = prev_gate_dict
                
            
            prop_gates = gate_dict 
            openlevels = [k for k,v in prop_gates.items() if v>0]
            lowestopen = min(openlevels)
            highestopen = max(openlevels)
            lowestposs = min(posslevels)
            highestposs = max(posslevels)
            
            # check that if only a single gate level is open that all 5 gates are
            # open
            if lowestopen==highestopen:
                if lowestopen>0 and sum(prop_gates)<5:
                    print("\n\t\t##### adjusting gate levels under full search1!")
                    prop_gates[lowestopen] = 5
                    new_gate_dict = prev_gate_dict
            else:
                new_gate_dict = prop_gates
            
                
            if gate_dict != prev_gate_dict:
                sha.SimulationSpecs.DaysSinceLastGateChange = 0
            else:
                sha.SimulationSpecs.DaysSinceLastGateChange +=1
                
            sha.LP_Opt_Blending = True
            
            # [q_o3, t_o3, e_o3 ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
            #                                          rivDict,
            #                                          bypassFrac=bypassFrac,
            #                                          gateDictOpt=gate_dict)
            
            # kwk_in_dict3 = {'inflow_final': q_o3, 'inflowTemp_final': t_o3}
            # kwk.set_inflow_temps(kwk_in_dict3, set_date=d,
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
                        
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            # this_temp3 = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            # # set final adj shasta target
            # sha_sac_tdiff3 = (this_temp3*1.8+32.0) - t_o3
               
            # sha_adj_targ3 = (targ_data[0]*1.8+32)-sha_sac_tdiff3
            # #sha_adj_targ = t_o + sac_targ_diff # targ_data[0]*1.8+32 - sha_sac_tdiff 
            
            # # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ3
            
            
        elif gate_mode=='incr':
            print("\t---doing incremental gate search")
            
            sha.LP_Opt_Blending = False
            
            # # use the data from the initial advance_swd run at top (using last
            # # time step's gate configuration)
            # kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
            # kwk.set_inflow_temps(kwk_in_dict,set_date=d,
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
            
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            
            #this_temp_incr = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            this_temp_incr = t_o # using temp from initial advance_swd run at top
            
            # current actual difference between shasta and target, 
            # pos value = warmer river than sha; neg value = shas warmer than river
            # sha_sac_tdiff = (this_temp_incr*1.8+32.0) - t_o
            # #difference between current temp at target loc and target temp
            # sac_targ_diff = (targ_data[0]-this_temp)*1.8
            
            # sha_adj_targ_incr = (targ_data[0]*1.8+32)-sha_sac_tdiff
            #sha_adj_targ = t_o + sac_targ_diff # targ_data[0]*1.8+32 - sha_sac_tdiff 
            
            # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ_incr
            
            # do incrementatl gate search
            prev_gates = sha.Operations.GateOps[d-dt.timedelta(1)]
            openlevels = [k for k,v in prev_gates.items() if v>0]
            posslevels = rt.outflow.gate_level_opts(sha)
            lowestopen = min(openlevels)
            highestopen = max(openlevels)
            lowestposs = min(posslevels)
            highestposs = max(posslevels)
            #too_cold = this_temp_incr < (targ_data[0] - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
            #too_warm = this_temp_incr > (targ_data[0] + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate
            # adjusting to reflect evaluation at shasta rather than river
            this_sha_targ = sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']]
            too_cold = this_temp_incr < (this_sha_targ - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
            too_warm = this_temp_incr > (this_sha_targ + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate

            if too_warm:
            
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
                        new_gate_dict[lowestposs] = sha.Outlets[lowestposs].NumGates
                else: # then there are already 2 levels open - try opening lower or closing upper
                    #print("\tThere are two levels already open; we'll have to try closing an upper or opening a lower...")
                    new_gate_dict = copy.deepcopy(prev_gate_dict)
                    if prev_gate_dict[lowestopen] < sha.Outlets[lowestopen].NumGates:
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
                        new_gate_dict[highestposs] = sha.Outlets[highestposs].NumGates
                else: # there are 2 levels already open - try closing a lower or opening an upper
                    #print("\tThere are two levels already open; we'll have to try closing a lower or opening an upper...")
                    new_gate_dict = copy.deepcopy(prev_gate_dict)
                    if prev_gate_dict[highestopen] < sha.Outlets[highestopen].NumGates:
                        new_gate_dict[highestopen] = prev_gate_dict[highestopen]+1
                    else:
                        new_gate_dict[lowestopen] = max(0, prev_gate_dict[lowestopen]-1)
            else:
                # continuing with previous gate configuration
                # check that if only a single gate level is open that all 5 gates are
                # open
                if lowestopen==highestopen:
                    print(f"###### only one gate open - lowest open: {lowestopen} & sum prev gates: {sum(prev_gates.values())}")
                    if lowestopen>0 and sum(prev_gates.values())<5:
                        print(f'############# setting gate level {lowestopen} to 5')
                        prop_gates = copy.deepcopy(prev_gates)
                        prop_gates[lowestopen] = 5
                        new_gate_dict = prop_gates #prev_gate_dict
                    else:
                        new_gate_dict = prev_gates
                else:
                    new_gate_dict = prev_gates #prev_gate_dict
                
            # # check the result using the new_gate_dict - if the temperature is too warm
            # # and the amount by which it exceeds the tolerance is greater than the
            # # amount by which it was too cold originally, then revert back to the original
            # # # configuration - trying to reduce oscillations
            [q_o2, t_o2, e_o2] = rt.outflow.calcOutTemp(sha,thisOutFlowTot, 
                                              rivDict,
                                              bypassFrac=bypassFrac,
                                              nPtSinks=3, gateDictOpt = new_gate_dict)
            
            # # run through downstream models
            # kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
            # kwk.set_inflow_temps(kwk_in_dict, set_date=d, 
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
            
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            # this_temp_incr2 = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            # sha_sac_tdiff = (this_temp_incr2*1.8+32.0) - t_o
            # sha_adj_targ_incr2 = (targ_data[0]*1.8+32)-sha_sac_tdiff
            
            # # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ_incr2
            
            
            # too_warm_amt = max(0.,this_temp_incr2- targ_data[0] + targ_data[1])
            # too_cold_amt = max(0., targ_data[0] -targ_data[2]-this_temp_incr2)
            too_warm_amt = max(0.,t_o2 - this_sha_targ + targ_data[1])
            too_cold_amt = max(0., this_sha_targ -targ_data[2]-t_o2)
            
            if (too_warm_amt < too_cold_amt): # and (sha.TimeStepDate.month < 1):
                print("------ reverting to previous gate config")
                new_gate_dict = prev_gate_dict

            gate_dict = new_gate_dict 
            
            # ensure when only one gate level is activet that it 
            # isn't sittng with a single isolated gate open
            prop_gates = copy.deepcopy(new_gate_dict)
            openlevels = [k for k,v in prop_gates.items() if v>0]
            posslevels = rt.outflow.gate_level_opts(sha)
            lowestopen = min(openlevels)
            highestopen = max(openlevels)
            lowestposs = min(posslevels)
            highestposs = max(posslevels)
            
            # check that if only a single gate level is open that all 5 gates are
            # open
            if lowestopen==highestopen:
                if lowestopen>0 and sum(prop_gates.values())<5:
                    print("\n\t\t##### adjusting gate levels under incremental search1!")
                    prop_gates[lowestopen] = 5
                    gate_dict = prop_gates
                else:
                    gate_dict = prop_gates
            else:
                gate_dict = prop_gates
            
            
            
            if new_gate_dict != prev_gate_dict:
                sha.SimulationSpecs.DaysSinceLastGateChange = 0
            else:
                sha.SimulationSpecs.DaysSinceLastGateChange +=1
                
            sha.LP_Opt_Blending = True
            
            
        elif gate_mode=='spec':
            # gates = gate_dict - can finalize and move to next step
            gate_dict = gates
            sha.SimulationSpecs.DaysSinceLastGateChange +=1
            
            # get shasta temps using previous gate config
            # [q_o, t_o, e_o ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
            #                                          rivDict,
            #                                          bypassFrac=bypassFrac,
            #                                          gateDictOpt=gate_dict)
            # # propagate through downstream to get an approx correction to Shasta
            # # tailwater target - necessary for lp optimized blending;
            # kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
            # kwk.set_inflow_temps(kwk_in_dict,set_date=d,
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
            
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            # this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            
            # # current actual difference between shasta and target, 
            # # pos value = warmer river than sha; neg value = shas warmer than river
            # sha_sac_tdiff = (this_temp*1.8+32.0) - t_o
            # #difference between current temp at target loc and target temp
            # sac_targ_diff = (targ_data[0]-this_temp)*1.8
            
            # sha_adj_targ = (targ_data[0]*1.8+32)-sha_sac_tdiff
            # #sha_adj_targ = t_o + sac_targ_diff #targ_data[0]*1.8+32 - sha_sac_tdiff 
            
            # # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ
            
            
        elif gate_mode=='prev':
            print('using previous time step gate configuration..')
            
            prop_gates = copy.deepcopy(sha.Operations.GateOps[d-dt.timedelta(1)])
            openlevels = [k for k,v in prop_gates.items() if v>0]
            lowestopen = min(openlevels)
            highestopen = max(openlevels)

            # check that if only a single gate level is open that all 5 gates are
            # open
            if lowestopen==highestopen:
                if lowestopen>0 and sum(prop_gates.values())<5:
                    prop_gates[lowestopen] = 5
                    new_gate_dict = prop_gates
            else:
                new_gate_dict = prop_gates
            
            gate_dict= new_gate_dict #gates
            sha.SimulationSpecs.DaysSinceLastGateChange +=1

            # # get shasta temps using previous gate config
            # [q_prevgate, t_prevgate, e_prevgate ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
            #                                          rivDict,
            #                                          bypassFrac=bypassFrac,
            #                                          gateDictOpt=gate_dict)
            # # propagate through downstream to get an approx correction to Shasta
            # # tailwater target - necessary for lp optimized blending; 
            # kwk_in_dict = {'inflow_final': q_prevgate, 'inflowTemp_final': t_prevgate}
            # kwk.set_inflow_temps(kwk_in_dict,set_date=d,
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
            
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            # this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            
            # # current actual difference between shasta and target, 
            # # pos value = warmer river than sha; neg value = shas warmer than river
            # sha_sac_tdiff = (this_temp*1.8+32.0) - t_o
            # #difference between current temp at target loc and target temp
            # sac_targ_diff = (targ_data[0]-this_temp)*1.8
            
            # sha_adj_targ = (targ_data[0]*1.8+32)-sha_sac_tdiff
            # #sha_adj_targ = t_o + sac_targ_diff #targ_data[0]*1.8+32 - sha_sac_tdiff 
            
            # prev_sha_targ = sha.Outflow.DataFrame.loc[d-dt.timedelta(1), sha.Outflow.ColumnMap['tempTarg']]
            # # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ

            
        else:
            gate_dict = prev_gate_dict
            sha.SimulationSpecs.DaysSinceLastGateChange +=1
            
            # # get shasta temps using previous gate config
            # [q_o, t_o, e_o ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,
            #                                          rivDict,
            #                                          bypassFrac=bypassFrac,
            #                                          gateDictOpt=gate_dict)
            # # propagate through downstream to get an approx correction to Shasta
            # # tailwater target - necessary for lp optimized blending; 
            # kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
            # kwk.set_inflow_temps(kwk_in_dict,set_date=d,
            #                      temp_units='DEG_F', flow_units='AF')
            # kwk.advance_longtemp(final=False)
            
            # uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
            #                           kwk.IntermedRelease/86400, set_date=d)
            # uppsac.advance_rivmod(final=False)
            # this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            
            # # current actual difference between shasta and target, 
            # # pos value = warmer river than sha; neg value = shas warmer than river
            # sha_sac_tdiff = (this_temp*1.8+32.0) - t_o
            # #difference between current temp at target loc and target temp
            # sac_targ_diff = (targ_data[0]-this_temp)*1.8
            
            # sha_adj_targ = (targ_data[0]*1.8+32)-sha_sac_tdiff
            # #sha_adj_targ = t_o + sac_targ_diff #targ_data[0]*1.8+32 - sha_sac_tdiff 
            
            # # update the temperature target at Shasta tailwater
            # sha.Outflow.DataFrame.loc[d, sha.Outflow.ColumnMap['tempTarg']] = sha_adj_targ
        
        #sha.LP_Opt_Blending = False
        prev_gate_mode = gate_mode
        # set the gate dict
        sha.Operations.GateOps[d] = gate_dict

        if gate_mode=='spec':
            [totQ2, outT, outE] = sha.advance_swd(final=True, gate_dict = gate_dict,
                                                  riv_dict = rivDict,
                                                  bypass_frac = bypass_frac,
                                                  temp_target = targ_data, return_vals=True )            
        else:
            [totQ2, outT, outE] = sha.advance_swd(final=True, gate_dict = gate_dict,
                                                  temp_target = targ_data, return_vals=True )
        
        kwk_in_dict = {'inflow_final': totQ2, 'inflowTemp_final': outT}
        kwk.set_inflow_temps(kwk_in_dict, set_date=d, temp_units='DEG_F', flow_units='AF')
        kwk.advance_longtemp(final=True)
        
        uppsac.set_node_flow_temp('kwk', kwk.Temp, 
                                  kwk.ReleaseVolumes[-1]/86400, set_date=d)
        uppsac.advance_rivmod(final=True)
        
    
        cpl.Temperature_Target.Sim_TS[d] = uppsac.Nodes[cpl.Temperature_Target.Target_Location.lower()].OutflowTemp[-1]
        
    etime_this = time.perf_counter()
    runtime_dict[yr] = etime_this - btime_this
    
#%%  quick and dirty plots

    ff =uppsac.plot_sim_obs(cpl.Temperature_Target.Target_Location,
                            plot_units='degC', year=yr,
                        temp_targ = cpl, add_tol=True, plot_off=False)
    ff.show()

#%% check the shasta target

fig, ax = plt.subplots(1,1, figsize=(9,7))
sha_targ_ts = sha.Outflow.DataFrame.loc[f"{yr}", 'Tw_target_degF']
sha_targ_ts[sha_targ_ts>=99] = np.nan
ax.plot(sha_targ_ts)
plt.show()
#%%
    sha.finalize()
    # gather the data and organize for calculating error metrics
    mkp.get_sim_obs_data_this_run(sha)
    
    rfig = mkp.plotReleasesCompare(sha, sha.Simulation_Results['ReleaseDF'], 
                                   obs_label = 'Obs Tailwater',
                                   other_temp = 'temperature_other',
                                   other_label = 'TCD Wt Avg',
                                   viewSave='view',
                                   on_wy=False, predicted_gates=True)
#%%    
    # plot results
    # [1] - river temperature at target
    plt.ioff()

    if save_plots:
        # do the plots in at CCR          
        ff =uppsac.plot_sim_obs(cpl.Temperature_Target.Target_Location,
                                plot_units='degC', year=yr,
                            temp_targ = cpl, add_tol=True, plot_off=False)
        ccrofp = os.path.join(outdir, f'{todaystr}DRAFT_CCR_RivTempFcstTarg__{yr}_degC.png')
        plt.savefig(ccrofp, dpi=600)
        plt.close()
    
        # plot the rest of the river nodes
        for n in uppsac.Nodes:
            if n.lower() != 'ccr':
                if uppsac.Nodes[n].UpNode>=0:
                    f = uppsac.plot_sim_obs(n, plot_units='degC', year=yr, plot_off=True)
                    f.tight_layout()
                    ofn = f'{todaystr}DRAFT_{n.upper()}_RivTempFcstTarg__{yr}_degC.png'
                    f.savefig(os.path.join(outdir, ofn), dpi=600)
                    plt.close()
                    
    sha.finalize()
    # gather the data and organize for calculating error metrics
    mkp.get_sim_obs_data_this_run(sha)

    error_at_open_gates = mkp.getValsAtOpenGateLevels(sha)
    error_profiles = mkp.calcErrorMetrics_profiles(sha)
    all_profile_errors[yr] = error_profiles[yr]
    all_opengate_errors[yr] = error_at_open_gates[yr]
    
    side_gate_compare[yr] = mkp.side_gate_dates(sha)[yr]
    
    if save_plots:
        rfig = mkp.plotReleasesCompare(sha, sha.Simulation_Results['ReleaseDF'], 
                                       obs_label = 'Obs Tailwater',
                                       other_temp = 'temperature_other',
                                       other_label = 'TCD Wt Avg',
                                       viewSave='view',
                                       on_wy=False, predicted_gates=True)
        shrelfp = os.path.join(outdir, f'{todaystr}DRAFT_ShastaRelease_FcstTarg_{yr}.png')
        plt.savefig(shrelfp, dpi=600)
        plt.close()
        
        rfig2 = mkp.plotProfilesCompare2(sha, sha.Simulation_Results['ProfilesDF'],
                                         viewSave='view',on_wy=False, sim_gates=True,
                                         show_stats = False) #, predicted_gates=True)
        shproffp = os.path.join(outdir, f'{todaystr}DRAFT_ShastaProfiles_FcstTarg_{yr}.png')
        plt.savefig(shproffp, dpi=600)
        plt.close()
    
        f = mkp.plotStorageEvapCompare(sha, viewSave='view')
        evapstofp = os.path.join(outdir, f'{todaystr}DRAFT_ShastaEvapStorage_FcstTarg_{yr}.png')
        plt.savefig(evapstofp, dpi=600)
        plt.close()
        
    #calculate the monthly error stats on temperature releases
    error_rel = mkp.calcErrorMetrics_timeseries(sha, ivar='rel', monthly=True)
    all_sha_rel_error[yr] = error_rel[yr]
    simT = mkp.pnd.DataFrame([_ for _ in kwk.StorageTemps[1:]], index=kwk.SimDates[0:len(kwk.StorageTemps[1:])])
    outTcol = kwk.Observations.OutflowColMap['obs_outtemp_final']
    kwk_obsT[yr] = mkp.pnd.DataFrame(kwk.Observations.OutflowDF.loc[kwk.SimDates[0:len(simT)],outTcol[0]])
    kwk_simT[yr] = mkp.pnd.DataFrame([_ for _ in kwk.StorageTemps[1:]], index=kwk.SimDates[0:len(kwk.StorageTemps[1:])])
    
    error_kwk_rel = kwk.calc_release_temp_errors(year=[yr])
    all_kwk_rel_error[yr] = error_kwk_rel[yr]
    
    # post-process the river model data
    df = uppsac.combine_to_df()
    ccr_err = uppsac.calc_error_metrics(df, 'ccr') #, monthly=True)      
    all_ccr_error[yr] = ccr_err[yr]
    #ccr_ann_err = uppsac.calc_error_metrics(df, 'ccr', monthly=False)   

    bsf_err = uppsac.calc_error_metrics(df, 'bsf')      
    all_bsf_error[yr] = bsf_err[yr]
    #bsf_ann_err = uppsac.calc_error_metrics(df, 'bsf', monthly=False)  

    rdb_err = uppsac.calc_error_metrics(df, 'rdb')      
    all_rdb_error[yr] = rdb_err
    #rdb_ann_err = uppsac.calc_error_metrics(df, 'rdb', monthly=False)    
    
    bnd_err = uppsac.calc_error_metrics(df, 'bnd')
    all_bnd_error[yr] = bnd_err


    
etime_all = time.perf_counter()
print(f"Took {etime_all-btime_all} seconds to run all years and make plots")    
plt.ion()

#%% write out the stats
draft = True

mkp.profile_error_to_csv(all_profile_errors, sha, draft=draft, proj_dir=proj_dir)
mkp.timeseries_error_to_csv(sha,all_sha_rel_error, model_name='', proj_dir=proj_dir)
mkp.timeseries_error_to_csv(kwk, all_kwk_rel_error, proj_dir=proj_dir, model_name='Kwk',
                            run_name='FcstTarg')
# gather the river error stats inot a big long-form data frame
# - columns: year, month, nse_degF, rmse_degF, mae_degF, mbias_degF
for stn,error_dict in zip(['ccr','bsf','bnd','rdb'],
                          [all_ccr_error, all_bsf_error, all_bnd_error, 
                           all_rdb_error]):
    
    alldat = []
    for y in error_dict.keys():
        ydat = error_dict[y]
        if list(ydat.keys())[0] == y:
            ydat = error_dict[y][y]
        for m in ydat.keys():
            mdat = ydat[m]
            valdat = [y, m, mdat['nse'], mdat['rmse'], mdat['mean_bias'],
                      mdat['mae'],mdat['r2'],mdat['rmse']*1.8,
                      mdat['mean_bias']*1.8,mdat['mae']*1.8 ]
            alldat.append(valdat)
            
    alldf = mkp.pnd.DataFrame(alldat, columns=['Year','Month','NSE_degC','RMSE_degC',
                                           'MeanBias_degC','MAE_degC','R2','RMSE_degF',
                                           'MeanBias_degF','MAE_degF'])  
    if draft:
        draftstr = 'DRAFT'
    else:
        draftstr = ''
        
    outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_FcstTarg_{stn}_error_stats.csv')
    alldf.to_csv(outfp, header=True)
    
# write out side gate comparison
sdg_df = mkp.pnd.DataFrame(index=list(side_gate_compare.keys()),data=np.empty((len(side_gate_compare),4)),
                       columns=['FirstSDG_Obs','FirstSDG_Sim','FullSDG_Obs','FullSDG_Sim'])
for yk,yv in side_gate_compare.items():
    sdg_df.loc[yk,:] = [yv['first_sdg']['obs'], yv['first_sdg']['sim'],
                        yv['full_sdg']['obs'], yv['full_sdg']['sim']]
outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_FcstTarg_SideGateCompare.csv')
sdg_df.to_csv(outfp, header=True)
    
#%% still ahve to plot the keswick release time series

#obsT = pnd.DataFrame(kwk.Observations.OutflowDF.loc[kwk.SimDates[0:len(simT)],outTcol[0]])
def c_to_f(c):
    return(c*1.8+32)
def f_to_c(f):
    return((f-32)/1.8)
def c_to_f_diff(c):
    return(c*1.8)
def f_to_c_diff(f):
    return(f/1.8)
def af_to_mcm(afv):
    mcmv = afv/1e6*43560.*(12**3)*(2.54**3)/(100**3)
    return(mcmv)
def mcm_to_af(mcmv):
    mafv = mcmv*(100**3)/(2.54**3)/(12**3)/43560
    return(mafv*1e6)
def cfs_to_cms(cfs):
    return(cfs*lt.FT3toM3)
def cms_to_cfs(cms):
    return(cms/lt.FT3toM3)

years = list(range(2000, 2023)) + ['all']
for y in years: #range(2000,2023):
    if y=='all':
        obsT = mkp.pnd.concat([kwk_obsT[iy] for iy in range(2000,2023)])
        simT = mkp.pnd.concat([kwk_simT[iy] for iy in range(2000,2023)])
    else:
        obsT= kwk_obsT[y]
        simT = kwk_simT[y]
    with mkp.sns.plotting_context('notebook', font_scale=1.2):
        
        import matplotlib.gridspec as gridspec
        
        fig2 = plt.figure(figsize=(12,8))
        spec2 = gridspec.GridSpec(ncols=3, nrows=9, figure=fig2)
        
        f2_ax1 = fig2.add_subplot(spec2[0:6,:])
        f2_ax2 = fig2.add_subplot(spec2[7:,:])
        
        #fig, ax = plt.subplots(1,1, figsize=(10,6))
        if y=='all':
            f2_ax1.plot(obsT, c='#3b4d7a', label='Observed') 
            f2_ax1.plot(simT, c='#ed6002', label='Simulated') 
            f2_ax1.xaxis.set_major_locator(mkp.mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
            f2_ax1.xaxis.set_minor_locator(mkp.mdates.MonthLocator(bymonth=[6]))
            f2_ax2.xaxis.set_major_locator(mkp.mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
            f2_ax2.xaxis.set_minor_locator(mkp.mdates.MonthLocator(bymonth=[6]))
            f2_ax1.tick_params(axis='x', labelsize=11 )
            f2_ax2.tick_params(axis='x', labelsize=11 )
            f2_ax1.xaxis.set_major_formatter(mkp.years_fmt)
            f2_ax2.xaxis.set_major_formatter(mkp.years_fmt)   
            f2_ax1.xaxis.set_tick_params(rotation=75)
            f2_ax2.xaxis.set_tick_params(rotation=75)
        else:
            f2_ax1.plot(obsT.loc[str(y)], c='#3b4d7a', label='Observed') #'2022'])
            f2_ax1.plot(simT.loc[str(y)], c='#ed6002', label='Simulated') #'2022'])
            f2_ax1.xaxis.set_major_locator(mkp.months)
            f2_ax1.xaxis.set_major_formatter(mkp.monyr_fmt)
            f2_ax2.xaxis.set_major_locator(mkp.months)
            f2_ax2.xaxis.set_major_formatter(mkp.monyr_fmt)
        
        #f2_ax1.set_xlabel('Date')

        f2_ax1.legend(loc='best', frameon=False)
        
        f2_ax1.yaxis.set_minor_locator(mkp.AutoMinorLocator())
        f2_ax1.set_ylabel('$^oC$', fontsize=14)
        if y=='all':
            f2_ax1.set_title(f'Keswick Release Temperatures', 
                         fontweight='bold', fontsize=16)
        else:   
            f2_ax1.set_title(f'Keswick Release Temperatures - {y}', 
                         fontweight='bold', fontsize=16)

        plt.subplots_adjust(left=0.19, bottom=0.15)
        ax2 = f2_ax1.secondary_yaxis( -0.10, functions=(c_to_f, f_to_c))
        ax2.yaxis.set_major_locator(mkp.AutoLocator()) #tkr.MultipleLocator(1))
        ax2.yaxis.set_minor_locator(mkp.AutoMinorLocator()) #tkr.MultipleLocator(0.2))
        ax2.set_ylabel('$^oF$',fontsize=14)
        
        
        
        f2_ax1.annotate('Keswick Release Temperatures', (0.1, 0.43), xycoords='figure fraction',
                      xytext=(0.04, 0.425), textcoords='figure fraction', 
                      fontsize=14, fontweight='bold', rotation=90) 
        
        if y=='all':
            f2_ax2.plot(simT.index, 
                        simT.iloc[:,0]-obsT.iloc[:,0] , 
                        label='Diff',c='k')     
        else:
            f2_ax2.plot(simT.loc[str(y)].index, 
                        simT.loc[str(y)].values-obsT.loc[str(y)].values , 
                        label='Diff',c='k')
            
        f2_ax2.axhline(0,ls='-',c='0.3', lw=0.8, zorder=0)
        f2_ax2.axhline(1, ls='--',c='0.5',lw=0.3,zorder=0)
        f2_ax2.axhline(-1, ls='--', c='0.5',lw=0.3, zorder=0)
        ax2b = f2_ax2.secondary_yaxis( -0.10, functions=(c_to_f_diff, f_to_c_diff))
        f2_ax2.yaxis.set_major_locator(mkp.MultipleLocator(1))
        # f2_ax1.yaxis.set_major_locator(MultipleLocator(1))
        # f2_ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
        
        ax2b.set_ylabel(f'Diff(Sim-Obs)\n$^oF$', fontsize=12,
                        fontweight='bold')
        f2_ax2.set_ylabel('$^oC$', fontsize=12)
    
        mkp.sns.despine()

        if y=='all':
            plt.savefig(os.path.join(outdir, 
                                     f'{todaystr}DRAFT_Fcst_Keswick_ReleaseTemps_allyears.png'),
                        dpi=600) # kwargs)
        else:
            plt.savefig(os.path.join(outdir, 
                                     f'{todaystr}DRAFT_Fcst_Keswick_ReleaseTemps_{y}.png'),
                        dpi=600) 
        plt.close()
plt.ion()

#%%        

# sha.finalize()
# rfig = mkp.plotReleasesCompare(sha, sha.Simulation_Results['ReleaseDF'],
#                                obs_label = 'Obs Tailwater',
#                                other_temp = 'temperature_other',
#                                other_label = 'TCD Wt Avg', viewSave='view',
#                                on_wy=False, predicted_gates=True)

#%%
rfig2 = mkp.plotProfilesCompare2(sha, sha.Simulation_Results['ProfilesDF'],
                                 viewSave='view',on_wy=False) #, predicted_gates=True)
 #%%

targ_sim_df = pnd.DataFrame.from_dict(cpl.Temperature_Target.Sim_TS, orient='index')
targ_targ_df = cpl.Temperature_Target.Target_TS.DataFrame.loc[targ_sim_df.index]
#plt.plot(targ_sim_df)
#plt.plot(targ_targ_df.Spec_Target)
#plt.plot(targ_sim_df.index, kwk.ReleaseTemps)
plt.plot(uppsac.SimDates,uppsac.Nodes['ccr'].OutflowTemp[0:len(uppsac.SimDates)])
obst = (uppsac.Nodes['ccr'].Observations.DataFrame.loc[str(yr)]-32)/1.8
plt.plot(obst, c='k')
