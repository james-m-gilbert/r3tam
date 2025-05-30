# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 18:06:29 2022

@author: jgilbert
"""

import matplotlib.pyplot as plt 

import restemp as rt
import longtemp as lt
from river import River
import coupling

import datetime as dt
import copy

# read in the config files for each model
config_fp01 = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/shasta_only/shasta_only_hist.yaml'
config_fp02 = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/keswick_only/kwk_input.v2.yaml'
config_fp03 = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/river_example/uppsac_river_kwk2ccr.yaml'
couple_config_fp = r'd:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/coupled_shakwkriv.yaml'
#%%
sha = rt.Res.initialize_model(config_fp01, profile_temp_units='degF')
kwk = lt.LongTemp.initialize_longmod(config_fp02, showInit=True)
uppsac = River.initialize_model(config_fp03)
sha.SeasonalRad = 0. 


cpl = coupling.coupled_models()
cpl.initialize(couple_config_fp, models=[sha, kwk, uppsac])
#%%
targts_shape = coupling.create_temp_target_shaping(12, 13.3, 24, 8, 1, temp_units='degC',
                                                   pre_ramp_period=2, post_ramp_period=2,
                                                   year=2014, start_doy=1, end_doy=365)


for d in cpl.SimDates:
    kwk_vec = kwk.Inflow.DataFrame.loc[d]
    
    sha.advance_restemp()
    #[totQ2, outT, outE] = sha.advance_swd(final=True)
    [totQ2, outT, outE] = sha.advance_swd(final=True, return_vals=True)
    kwk_in_dict = {'inflow_final': totQ2, 'inflowTemp_final': outT}
    kwk.set_inflow_temps(kwk_in_dict, temp_units='DEG_F', flow_units='AF')
    kwk.advance_longtemp(final=True)
    
    uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, kwk.IntermedRelease/86400)
    uppsac.advance_rivmod(final=True)
    
sha.finalize() #<-- put results in dataframes
#%% make some plots

obs_shaout_temp = sha.Observations.OutflowDF['Tw_Temp_degF']
sim_shaout_temp = sha.Simulation_Results['ReleaseDF'].loc[:,'Sim_Release_Temp_degF']

fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(obs_shaout_temp, label='Obs Temp')
ax.plot(sim_shaout_temp, label='Sim Temp')

#%%
check_kwkin = kwk.Inflow.DataFrame['InflowTotal']/lt.AFtoM3
check_shaout = sha.Simulation_Results['ReleaseDF'].loc[:, 'Sim_Release_AF']

check_kwkin_temp = kwk.Inflow.DataFrame['InflowTemp_degC']*1.8+32
check_shaout_temp = sha.Simulation_Results['ReleaseDF'].loc[:,'Sim_Release_Temp_degF']

plt.plot(check_kwkin_temp)
plt.plot(check_shaout_temp)

#%% plot keswick

sim_kwk =  pnd.DataFrame([_ for _ in kwk.StorageTemps[1:]], index=kwk.SimDates[0:len(kwk.StorageTemps[1:])])
obs_kwk = kwk.Observations.OutflowDF.loc[sim_kwk.index, 'ObsOutTemp']

fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(obs_kwk, label='Obs Temp')
ax.plot(sim_kwk, label='Sim Temp')
#%%

sim_ccr = uppsac.Nodes['ccr'].OutflowTemp
obs_ccr = uppsac.Nodes['ccr'].Observations.DataFrame.loc[cpl.SimDates, 'qc1']
sim_ccr_degF = [_*1.8+32 for _ in sim_ccr]
fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(obs_ccr, label='Obs Temp')
ax.plot(cpl.SimDates, sim_ccr_degF, label='Sim Temp')


#%%  now try using gate selection logic to hit targets
#  at selected stream nodes (ccr)
sha = rt.Res.initialize_model(config_fp01, profile_temp_units='degF')
kwk = lt.LongTemp.initialize_longmod(config_fp02, showInit=True)
uppsac = River.initialize_model(config_fp03)

cpl = coupling.coupled_models()
cpl.initialize(couple_config_fp, models=[sha, kwk, uppsac])

sha.SimulationSpecs.IsHindcast = False #<-- switch so that calcGateOps is engaged
sha.SimulationSpecs.CalcGateOps = True
sha.Debug['Release'] = 0

for yr in [2014]: #, 2015, 2016]:
    
    targts_shape = coupling.create_temp_target_shaping(13.3, 13.3, 24, 8, 1, temp_units='degC',
                                                       pre_ramp_period=2, post_ramp_period=2,
                                                       year=yr, start_doy=1, end_doy=365)

    # set the new target ts as the specified in cpl model
    cpl.Temperature_Target.Target_TS.DataFrame.loc[targts_shape.index,'Spec_Target'] = targts_shape['target_shape_degC']
    cpl.Temperature_Target.Target_TS.ColumnMap['ttarg_specified'] = 'Spec_Target'
    
    cpl.SimDates = targts_shape.index
    sha.TimeStep = -1
    sha.TimeStepDate = cpl.SimDates[0]
    sha.reinitialize(sim_dates=targts_shape.index)

    kwk.TimeStep = -1
    kwk.TimeStepDate = cpl.SimDates[0]
    kwk.reintialize(sim_dates=targts_shape.index)
    
    uppsac.TimeStep = -1
    uppsac.TimeStepDate = cpl.SimDates[0]
    uppsac.SimDates = targts_shape.index
    
    
    for d in cpl.SimDates:
        sha.advance_restemp() #advances the time step counter for the sha model, does inflow and met calculations
        
        # do an initial check of gate options
        [gate_mode, gates, targ_data, prev_err, rivDict, bypass_frac] = sha.advance_swd(final=False, coupling_model=cpl, 
                                                                  check_gate_calcs=True)
        print(f"gate mode: {gate_mode}")
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
            
        if gate_mode=='full':
            # do full gate selection
            best_error = 99
            best_gates = prev_gate_dict
            
            for gc in gates: # iterate through all gate options
                [q_o, t_o, e_o ] =rt.outflow.calcOutTemp(sha, thisOutFlowTot,gc,
                                                         rivDict,
                                                         bypassFrac=bypassFrac)
                
                kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
                kwk.set_inflow_temps(kwk_in_dict,set_date=d,
                                     temp_units='DEG_F', flow_units='AF')
                kwk.advance_longtemp(final=False)
                
                uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
                                          kwk.IntermedRelease/86400, set_date=d)
                uppsac.advance_rivmod(final=False)
                this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
                this_tdiff = this_temp -targ_data[0]
                
                if abs(this_tdiff) < abs(best_error):
                    print(f">>>> added {gc} to best gates")
                    best_error = this_tdiff
                    best_gates = gc
                    
                if sha.Debug['Release'] >0:
                    print("       Target: %0.2f - %0.2f - %0.2f deg F" %(targ_data[0]+targ_data[1], targ_data[0], targ_data[0]-targ_data[2]))
                    print("       ResultingTemp:  %0.2f deg F" %this_temp)
                too_cold = this_temp < (targ_data[0] - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
                too_warm = this_temp > (targ_data[0] + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate
                
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
                gate_dict = prev_gate_dict
                
                
            if gate_dict != prev_gate_dict:
                sha.SimulationSpecs.DaysSinceLastGateChange = 0
            else:
                sha.SimulationSpecs.DaysSinceLastGateChange +=1
            
        elif gate_mode=='incr':
            print("\n\n\t---doing incremental gate search")
            # do incrementatl gate search
            prev_gates = sha.Operations.GateOps[d-dt.timedelta(1)]
            openlevels = [k for k,v in prev_gates.items() if v>0]
            posslevels = rt.outflow.gate_level_opts(sha)
            lowestopen = min(openlevels)
            highestopen = max(openlevels)
            lowestposs = min(posslevels)
            highestposs = max(posslevels)
            too_cold = this_temp < (targ_data[0] - targ_data[2]) # if evalulates as True, then temperature is too low - prefer to open an upper gate
            too_warm = this_temp > (targ_data[0] + targ_data[1]) # if evaluates as True, then temp is too warm, prefer to open a lower gate

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
                new_gate_dict = prev_gate_dict
                
            # check the result using the new_gate_dict - if the temperature is too warm
            # and the amount by which it exceeds the tolerance is greater than the
            # amount by which it was too cold originally, then revert back to the original
            # configuration - trying to reduce oscillations
            [q_o, t_o, e_o] = rt.outflow.calcOutTemp(sha,thisOutFlowTot, new_gate_dict,
                                              rivDict,
                                              bypassFrac=bypassFrac,
                                              nPtSinks=3)
            
            # run through downstream models
            kwk_in_dict = {'inflow_final': q_o, 'inflowTemp_final': t_o}
            kwk.set_inflow_temps(kwk_in_dict, set_date=d, 
                                 temp_units='DEG_F', flow_units='AF')
            kwk.advance_longtemp(final=False)
            
            uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
                                      kwk.IntermedRelease/86400, set_date=d)
            uppsac.advance_rivmod(final=False)
            this_temp = uppsac.Nodes[cpl.Temperature_Target.Target_Location].OutflowTemp[-1]
            
            too_warm_amt = max(0.,this_temp- targ_data[0] + targ_data[1])
            too_cold_amt = max(0., targ_data[0] -targ_data[2]-this_temp)
            
            if (too_warm_amt > too_cold_amt): # and (sha.TimeStepDate.month < 1):
                print("------ reverting to previous gate config")
                new_gate_dict = prev_gate_dict

            gate_dict = new_gate_dict 
            
            if new_gate_dict != prev_gate_dict:
                sha.SimulationSpecs.DaysSinceLastGateChange = 0
            else:
                sha.SimulationSpecs.DaysSinceLastGateChange +=1
            
        elif gate_mode=='spec':
            # gates = gate_dict - can finalize and move to next step
            gate_dict = gates
            sha.SimulationSpecs.DaysSinceLastGateChange +=1
            
        elif gate_mode=='prev':
            gate_dict= gates
            sha.SimulationSpecs.DaysSinceLastGateChange +=1
        else:
            gate_dict = prev_gate_dict
            sha.SimulationSpecs.DaysSinceLastGateChange +=1
        
        # set the gate dict
        sha.Operations.GateOps[d] = gate_dict
        #sha.Outflow.ColumnMap
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
        
        
#%%        
import make_plots as mkp
sha.finalize()
rfig = mkp.plotReleasesCompare(sha, sha.Simulation_Results['ReleaseDF'], viewSave='view',
                               on_wy=False, predicted_gates=True)

#%%

targ_sim_df = pnd.DataFrame.from_dict(cpl.Temperature_Target.Sim_TS, orient='index')
targ_targ_df = cpl.Temperature_Target.Target_TS.DataFrame.loc[targ_sim_df.index]
plt.plot(targ_sim_df)
plt.plot(targ_targ_df.Spec_Target)
plt.plot(targ_sim_df.index, kwk.ReleaseTemps)
plt.plot(targ_sim_df.index,uppsac.Nodes['ccr'].OutflowTemp)
