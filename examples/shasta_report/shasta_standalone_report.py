# -*- coding: utf-8 -*-
"""
Shasta-only ResTemp simulation using R3TAM package

This file runs simulations and creates figures and statistical summaries that
are reported in the 2026 R3TAM technical memorandum (report)

@author: James Gilbert (james.gilbert@ucsc.edu)
"""

# include some common libraries
import sys, os
import pandas as pd
import numpy as np
import time


# include the R3TAM package components
from r3tam import restemp as rt
from r3tam import longtemp as lt
from r3tam.river import River
from r3tam import coupling

from r3tam import make_plots as mkp

#%%

# make sure the working directory for this script is set in the 
# 'examples/shasta' folder

# set the configuration filepath
config_fp = './shasta_standalone_report.yaml' #<-- if working directory not set, change this to full system filepath

# write out figures and stats
make_plots = True
write_stats = True

#%% initialize and run simulation 
btime = time.perf_counter() 
resmod = rt.Res.initialize_model(config_fp, profile_temp_units='degF')

resmod.Debug['Release'] = 0
resmod.Debug['General'] = 0
resmod.SeasSolRad = 1.
resmod.Debug['Release'] = 0
for d in resmod.SimDates:
    #print(d)
    resmod.advance_restemp()

    resmod.advance_swd(final=True)

etime = time.perf_counter() 
print(f"took {etime-btime} seconds")

resmod.finalize()

#%% The remaining code writes output and figures for the report

# gather the data and organize for calculating error metrics
mkp.get_sim_obs_data_this_run(resmod)

error_at_open_gates = mkp.getValsAtOpenGateLevels(resmod)
error_profiles = mkp.calcErrorMetrics_profiles(resmod)

error_rel = mkp.calcErrorMetrics_timeseries(resmod, ivar='rel', monthly=True)

# write releae and profile error stats
if write_stats:
    mkp.timeseries_error_to_csv(resmod,error_rel)
    
    mkp.profile_error_to_csv(error_profiles, resmod, draft=True)
    
# make plots for observered-simulated comparisons
obs_rel_std = np.nanstd(resmod.Observations.OutflowDF.loc[:,resmod.Observations.OutflowColMap['obs_outtemp_final']])
print(f"Half std dev of obs release temperatures is: {obs_rel_std/2:0.2f}")
if make_plots:
    rfig = mkp.plotReleasesCompare(resmod,
                                   resmod.Simulation_Results['ReleaseDF'],
                                   viewSave='save',
                                   #target_ts = 'Sim_Temp_Target_degF',
                                   #target_tol = 'Sim_Target_Tolerance_degF',
                                   obs_label = 'Obs Tailwater',
                                   other_temp = 'temperature_other',
                                   other_label = 'TCD Wt Avg',
                                   on_wy=False)    
    
    f= mkp.plotProfilesCompare2(resmod, 
                               resmod.Simulation_Results['ProfilesDF'], 
                               viewSave='save', show_stats=False,
                               on_wy=False)
    
    f = mkp.plotStorageEvapCompare(resmod, viewSave='save')

#%%  plot the scatter of observed and simulated values at open gate elevations
import datetime as dt
figoutDir = os.path.join(resmod.ProjDir, 'outputs', 'figs')
mkp.plt.ioff()
todaystr = f'{dt.date.today().strftime("%Y%m%d")}'
for y in error_at_open_gates:
    erropyr = error_at_open_gates[y]

    with mkp.sns.plotting_context('notebook', font_scale=1.4):
        fig, ax = mkp.plt.subplots(1,1, figsize=(12,9))
        for m in erropyr:
            if m in [11,12,1,2]:
                thisc = 'darkblue'
                label='Winter'
            elif m in [3,4,5]:
                thisc = 'darkgreen'
                label='Spring'
            elif m in [6,7,8]:
                thisc = 'orange'
                label='Summer'
            else:
                thisc = 'red'
                label='Fall'
            thism = erropyr[m]
            thissim = thism['sim_data']
            thisobs = thism['obs_data']
            for vo, vs in zip(thisobs,thissim):
                ax.plot(vo, vs, ls='None', marker='o', 
                        markerfacecolor=thisc,ms=10,
                        markeredgecolor='None', label=label)
                
        minx, maxx = ax.get_xlim()
        ax.plot([minx, maxx], [minx, maxx], ls='--', c='k', lw=1.5, zorder=0)
        
        ax.xaxis.set_minor_locator(mkp.AutoMinorLocator())
        ax.yaxis.set_minor_locator(mkp.AutoMinorLocator())
        ax.set_xlabel('Observed Profile Temperature, $^oF$')
        ax.set_ylabel('Simulated Profile Temperature, $^oF$')
        ax.set_title(f'Profile Temperature Comparison at Open Gate Levels\n{y}',
                     fontweight='bold', fontsize=20)
        h1,l1 = ax.get_legend_handles_labels()
        l2 = list(set(l1))
        h2 = [h1[l1.index(i)] for i in set(l1)]
        mkp.plt.legend(h2, l2, loc='best')
        mkp.sns.despine()
        
        of = os.path.join(figoutDir, f'{todaystr}{resmod.RunName}_ProfileAtOpenGates_{y}.png')
        mkp.plt.savefig(of, dpi=450)
        mkp.plt.close()
        
mkp.plt.ion()    
    

#%% a version of the error-at-gate elevations plot, but for all years combined

with mkp.sns.plotting_context('notebook', font_scale=1.4):
    fig, ax = mkp.plt.subplots(2,2, figsize=(12,10))
    seasons = ['Winter', 'Spring', 'Summer', 'Fall']
    for y in error_at_open_gates:
        erropyr = error_at_open_gates[y]
        axidxs = [[0,0], [0,1], [1,0],[1,1]]
        for m in erropyr:
            if m in [11,12,1,2]:
                thisc = 'darkblue'
                label='Winter'
                axidx = 0
            elif m in [3,4,5]:
                thisc = 'darkgreen'
                label='Spring'
                axidx = 1
            elif m in [6,7,8]:
                thisc = 'orange'
                label='Summer'
                axidx = 2
            else:
                thisc = 'red'
                label='Fall'
                axidx = 3
                
            thism = erropyr[m]
            thissim = thism['sim_data']
            thisobs = thism['obs_data']
            for vo, vs in zip(thisobs,thissim):
                ax[axidxs[axidx][0], axidxs[axidx][1]].plot(vo, vs, ls='None', marker='o', 
                        markerfacecolor=thisc,ms=5,alpha=0.7,
                        markeredgecolor='None', label=label)
                
        
        for tx in [0,1,2,3]:
            
            minx, maxx = ax[axidxs[tx][0], axidxs[tx][1]].get_xlim()
            if tx==0:
                minx, maxx = 40, 65
            elif tx==1:
                minx, maxx = 40, 70
            elif tx==2:
                minx, maxx = 40, 85
            else:
                minx, maxx = 45, 70
            ax[axidxs[tx][0], axidxs[tx][1]].plot([minx, maxx], [minx, maxx], ls='--', c='k', lw=1.5, zorder=0)
        
            ax[axidxs[tx][0], axidxs[tx][1]].xaxis.set_minor_locator(mkp.AutoMinorLocator())
            ax[axidxs[tx][0], axidxs[tx][1]].yaxis.set_minor_locator(mkp.AutoMinorLocator())
            
            ax[axidxs[tx][0], axidxs[tx][1]].set_xlabel('Observed Profile Temperature, $^oF$', fontsize=14)
            ax[axidxs[tx][0], axidxs[tx][1]].set_ylabel('Simulated Profile Temperature, $^oF$', fontsize=14)
            ax[axidxs[tx][0], axidxs[tx][1]].set_title(f'{seasons[tx]}')
            
    mkp.plt.subplots_adjust(hspace=0.3)
    mkp.plt.suptitle(f'Profile Temperature Comparison at Open Gate Levels\nAll Years',
                 fontweight='bold', fontsize=20)
    #h1,l1 = ax.get_legend_handles_labels()
    #l2 = list(set(l1))
    #h2 = [h1[l1.index(i)] for i in set(l1)]
    #mkp.plt.legend(h2, l2, loc='best')
    mkp.sns.despine()
        
    of = os.path.join(figoutDir, f'{todaystr}{resmod.RunName}_ProfileAtOpenGates_AllYears.png')
    mkp.plt.savefig(of, dpi=450)
    mkp.plt.close()
        
mkp.plt.ion()   
