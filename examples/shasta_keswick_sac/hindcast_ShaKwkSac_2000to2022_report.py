# -*- coding: utf-8 -*-
"""
Coupled Shasta-Keswick-Sacramento River Hindcast simulation

This simulation corresponds to Section 3.4.1 in the 2026 R3TAM report 

@author: James Gilbert (james.gilbert@ucsc.edu)
"""

import matplotlib.pyplot as plt 

# include the R3TAM package components
from r3tam import restemp as rt
from r3tam import longtemp as lt
from r3tam.river import River
from r3tam import coupling

from r3tam import make_plots as mkp

import datetime as dt
import copy
import time
import numpy as np
import os

#%% Set the configuration file locations
# read in the config files for each model
curdir = os.path.abspath(os.curdir)
config_fp01 = os.path.join(curdir, r'../shasta_report/shasta_standalone_report.yaml')
config_fp02 = os.path.join(curdir, r'../keswick/kwk_input.v20240729.yaml')
config_fp03 = os.path.join(curdir, r'../sac_river/uppsac_river_kwk2rdb.yaml')
couple_config_fp = os.path.join(curdir, r'coupled_shakwkriv.yaml')
#%% Initialize the component models and the coupling scheme
sha = rt.Res.initialize_model(config_fp01, profile_temp_units='degF')
kwk = lt.LongTemp.initialize_longmod(config_fp02, showInit=True)
uppsac = River.initialize_model(config_fp03)
sha.SeasonalRad = 0. 

cpl = coupling.coupled_models()
cpl.initialize(couple_config_fp, models=[sha, kwk, uppsac])
#%%
# targts_shape = coupling.create_temp_target_shaping(12, 13.3, 24, 8, 1, temp_units='degC',
#                                                    pre_ramp_period=2, post_ramp_period=2,
#                                                    year=2014, start_doy=1, end_doy=365)

btime = time.perf_counter() 
for d in cpl.SimDates:
    kwk_vec = kwk.Inflow.DataFrame.loc[d]
    
    sha.advance_restemp()
    #[totQ2, outT, outE] = sha.advance_swd(final=True)
    [totQ2, outT, outE] = sha.advance_swd(final=True, return_vals=True)
    # if np.isnan(outT) or (outT==0): # couple fo instances in early 2001 when releases are ~0 - so no temp - have to give some useful number here
    #     outT = 48.0 
    kwk_in_dict = {'inflow_final': totQ2, 'inflowTemp_final': outT}
    kwk.set_inflow_temps(kwk_in_dict, temp_units='DEG_F', flow_units='AF')
    kwk.advance_longtemp(final=True)
    
    uppsac.set_node_flow_temp('kwk', kwk.IntermedTemp, 
                              kwk.IntermedRelease/86400, set_date=d)
    uppsac.advance_rivmod(final=True)
    
etime =  time.perf_counter() 
print(f"took {etime-btime} seconds") 
sha.finalize() #<-- put shasta results in dataframes
#%% collect, extract, organize data
import seaborn as sns
import pandas as pnd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as tkr
import datetime as dt

years = mdates.YearLocator()   # every year
fiveyears = mdates.YearLocator(5)
months = mdates.MonthLocator()  # every month
pltdays = mdates.DayLocator((5,10,15,20, 25))
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%b')
monyr_fmt = mdates.DateFormatter('%b\n%Y')

draft = True
#simT = pnd.DataFrame([_*1.8+32 for _ in kwk.StorageTemps[1:]], index=kwk.SimDates)
simT = pnd.DataFrame([_ for _ in kwk.StorageTemps[1:]], index=kwk.SimDates[0:len(kwk.StorageTemps[1:])])
outTcol = kwk.Observations.OutflowColMap['obs_outtemp_final']
obsT = pnd.DataFrame(kwk.Observations.OutflowDF.loc[kwk.SimDates[0:len(simT)],outTcol[0]])

proj_dir = curdir
outdir = os.path.join(proj_dir, 'outputs') 
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
figdir = os.path.join(outdir, 'figs')
if not os.path.exists(figdir):
    os.mkdri(figdir)
    
todaystr = dt.date.today().strftime('%Y%m%d')

#%% make Keswick release temperature comparison plots
import seaborn as sns

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


plt.ioff()
years = list(range(2000, 2023)) + ['all']
for y in years: #range(2000,2023):
    with sns.plotting_context('notebook', font_scale=1.2):
        
        import matplotlib.gridspec as gridspec
        
        fig2 = plt.figure(figsize=(12,8))
        spec2 = gridspec.GridSpec(ncols=3, nrows=9, figure=fig2)
        
        f2_ax1 = fig2.add_subplot(spec2[0:6,:])
        f2_ax2 = fig2.add_subplot(spec2[7:,:])
        
        #fig, ax = plt.subplots(1,1, figsize=(10,6))
        if y=='all':
            f2_ax1.plot(obsT, c='#3b4d7a', label='Observed') 
            f2_ax1.plot(simT, c='#ed6002', label='Simulated') 
            f2_ax1.xaxis.set_major_locator(mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
            f2_ax1.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[6]))
            f2_ax2.xaxis.set_major_locator(mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
            f2_ax2.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[6]))
            f2_ax1.tick_params(axis='x', labelsize=11 )
            f2_ax2.tick_params(axis='x', labelsize=11 )
            f2_ax1.xaxis.set_major_formatter(years_fmt)
            f2_ax2.xaxis.set_major_formatter(years_fmt)   
            f2_ax1.xaxis.set_tick_params(rotation=75)
            f2_ax2.xaxis.set_tick_params(rotation=75)
        else:
            f2_ax1.plot(obsT.loc[str(y)], c='#3b4d7a', label='Observed') #'2022'])
            f2_ax1.plot(simT.loc[str(y)], c='#ed6002', label='Simulated') #'2022'])
            f2_ax1.xaxis.set_major_locator(months)
            f2_ax1.xaxis.set_major_formatter(monyr_fmt)
            f2_ax2.xaxis.set_major_locator(months)
            f2_ax2.xaxis.set_major_formatter(monyr_fmt)
        
        #f2_ax1.set_xlabel('Date')

        f2_ax1.legend(loc='best', frameon=False)
        
        f2_ax1.yaxis.set_minor_locator(tkr.AutoMinorLocator())
        f2_ax1.set_ylabel('$^oC$', fontsize=14)
        if y=='all':
            f2_ax1.set_title(f'Keswick Release Temperatures', 
                         fontweight='bold', fontsize=16)
        else:   
            f2_ax1.set_title(f'Keswick Release Temperatures - {y}', 
                         fontweight='bold', fontsize=16)

        plt.subplots_adjust(left=0.19, bottom=0.15)
        ax2 = f2_ax1.secondary_yaxis( -0.10, functions=(c_to_f, f_to_c))
        ax2.yaxis.set_major_locator(tkr.AutoLocator()) #tkr.MultipleLocator(1))
        ax2.yaxis.set_minor_locator(tkr.AutoMinorLocator()) #tkr.MultipleLocator(0.2))
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
        f2_ax2.yaxis.set_major_locator(tkr.MultipleLocator(1))
        # f2_ax1.yaxis.set_major_locator(MultipleLocator(1))
        # f2_ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
        
        ax2b.set_ylabel(f'Diff(Sim-Obs)\n$^oF$', fontsize=12,
                        fontweight='bold')
        f2_ax2.set_ylabel('$^oC$', fontsize=12)
    
        sns.despine()

        if y=='all':
            plt.savefig(os.path.join(figdir, 
                                     f'{todaystr}DRAFT_Keswick_ReleaseTempsHindcast_allyears.png'),
                        dpi=600) # kwargs)
        else:
            plt.savefig(os.path.join(figdir, 
                                     f'{todaystr}DRAFT_Keswick_ReleaseTempsHindcast_{y}.png'),
                        dpi=600) 
        plt.close()
plt.ion()

#%% calculate error metrics for release temepratures

error_dict = {}
yrs = list(set([i.year for i in kwk.SimDates]))

simT = pnd.DataFrame(kwk.ReleaseTemps, index=kwk.SimDates, columns=['SimReleaseT_degC'])

for y in yrs:
     
    sim_data = simT.loc[str(y),:] #kwk.calc_errors_df[y][ivar]['sim']
    obs_data = obsT.loc[str(y), :] #resmod.calc_errors_df[y][ivar]['obs']
    
    monerr = {m:[] for m in range(1,13)}  # blank dictionary to hold monthly error calcs
    monobs = {m:[] for m in range(1,13)} # blank dictionary for the observations for purposes of calculating NSE
    monsim = {m:[] for m in range(1,13)} # blank dictionary for the simulation data
    month_stats = {}
    for m in range(1,13):
        thismidx = [ix for ix in sim_data.index if ix.month==m]
        
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
    
run_name = 'KWKcoupledhindcast'

stats_dir = os.path.join(outdir, 'stats')
if not os.path.exists(stats_dir):
    os.mkdir(stats_dir)
    
todaystr = f'{dt.date.today().strftime("%Y%m%d")}'

# gather the error stats inot a big long-form data frame
# - columns: year, month, nse_degF, rmse_degF, mae_degF, mbias_degF
alldat = []
for y in error_dict.keys():
    ydat = error_dict[y]
    for m in ydat.keys():
        mdat = ydat[m]
        valdat = [y, m, mdat['nse'], mdat['rmse'], mdat['mean_bias'],
                  mdat['mae'],mdat['r2'],mdat['rmse']*1.8,
                  mdat['mean_bias']*1.8,mdat['mae']*1.8 ]
        alldat.append(valdat)
        
alldf = pnd.DataFrame(alldat, columns=['Year','Month','NSE_degC','RMSE_degC',
                                       'MeanBias_degC','MAE_degC','R2','RMSE_degF',
                                       'MeanBias_degF','MAE_degF'])

if draft:
    draftstr = 'DRAFT'
else:
    draftstr = ''
    
outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_{run_name}_release_temp_error_stats.csv')
alldf.to_csv(outfp, header=True)


#%% make sac river plots
df = uppsac.combine_to_df()

plt.ioff()
yrstoplot = list(range(2000, 2023)) + ['']
for y in yrstoplot: #range(2000, 2023):
    for n in uppsac.Nodes:
        if uppsac.Nodes[n].UpNode>=0:
            if y=='':
                f=uppsac.plot_sim_obs(n, plot_units='degC')
                f.tight_layout()
                ofn = f'{todaystr}{draftstr}_{n.upper()}_RivTempHindcast_allyrs.png'
            else:
                f = uppsac.plot_sim_obs(n, plot_units='degC', year=y)
                f.tight_layout()
                ofn = f'{todaystr}{draftstr}_{n.upper()}_RivTempHindcast_{y}.png'
            f.savefig(os.path.join(figdir, ofn), dpi=600)
            plt.close()
plt.ion()      

#%%  write out error metrics

run_name = 'UPPSACcoupledhindcast'

ccr_err = uppsac.calc_error_metrics(df, 'ccr') #, monthly=True)      
#ccr_ann_err = uppsac.calc_error_metrics(df, 'ccr', monthly=False)   

bsf_err = uppsac.calc_error_metrics(df, 'bsf')      
#bsf_ann_err = uppsac.calc_error_metrics(df, 'bsf', monthly=False)  

rdb_err = uppsac.calc_error_metrics(df, 'rdb')      
#rdb_ann_err = uppsac.calc_error_metrics(df, 'rdb', monthly=False)    

bnd_err = uppsac.calc_error_metrics(df, 'bnd')

todaystr = f'{dt.date.today().strftime("%Y%m%d")}'
draft = True
# gather the error stats inot a big long-form data frame
# - columns: year, month, nse_degF, rmse_degF, mae_degF, mbias_degF
for stn,error_dict in zip(['ccr','bsf','bnd','rdb'],[ccr_err, bsf_err, bnd_err, rdb_err]):
    
    alldat = []
    for y in ccr_err.keys():
        ydat = error_dict[y]
        for m in ydat.keys():
            mdat = ydat[m]
            valdat = [y, m, mdat['nse'], mdat['rmse'], mdat['mean_bias'],
                      mdat['mae'],mdat['r2'],mdat['rmse']*1.8,
                      mdat['mean_bias']*1.8,mdat['mae']*1.8 ]
            alldat.append(valdat)
            
    alldf = pnd.DataFrame(alldat, columns=['Year','Month','NSE_degC','RMSE_degC',
                                           'MeanBias_degC','MAE_degC','R2','RMSE_degF',
                                           'MeanBias_degF','MAE_degF'])
    
    if draft:
        draftstr = 'DRAFT'
    else:
        draftstr = ''
        
    outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_{run_name}_{stn}_error_stats.csv')
    alldf.to_csv(outfp, header=True)


#%% make some plots

obs_shaout_temp = sha.Observations.OutflowDF['Tw_Temp_degF']
sim_shaout_temp = sha.Simulation_Results['ReleaseDF'].loc[:,'Sim_Release_Temp_degF']

fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(obs_shaout_temp, label='Obs Temp')
ax.plot(sim_shaout_temp, label='Sim Temp')

#%% check that keswick inflow equals shasta outflow
check_kwkin = kwk.Inflow.DataFrame['InflowTotal']/lt.AFtoM3
check_shaout = sha.Simulation_Results['ReleaseDF'].loc[:, 'Sim_Release_AF']

check_kwkin_temp = kwk.Inflow.DataFrame['InflowTemp_degC']*1.8+32
check_shaout_temp = sha.Simulation_Results['ReleaseDF'].loc[:,'Sim_Release_Temp_degF']

plt.plot(check_kwkin_temp)
plt.plot(check_shaout_temp)



