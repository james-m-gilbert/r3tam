# -*- coding: utf-8 -*-
"""
Running Keswick *LongTemp* model in isolation using historical upstream inputs
and releases.
The inputs and results of this simulation are reflected in the 2026 
R3TAM report, section 3.2.

@author: jgilbert
"""

# testing keswick model
from r3tam import longtemp as lt
import time 
import pandas as pnd
import os

inpKWKfp = r'./kwk_input.v20240729.yaml' #<-- make sure to set working directory to examples/keswick folder

btime = time.perf_counter()
kwk= lt.LongTemp.initialize_longmod(inpKWKfp, showInit=True) 
kwk.ProjDir = os.path.abspath(os.curdir)
etime = time.perf_counter()
print(f'Took {etime-btime} seconds to initialize')
#kwk.advance_longtemp(final=False)

kwk.Debug['General'] = 0
btime = time.perf_counter()
for t in kwk.SimDates:
    kwk.advance_longtemp(final=True)
    if kwk.Storage <=0:
        print("\n**************************")
        print(t)
        break
etime = time.perf_counter()

numyrs = len(list(set([y.year for y in kwk.SimDates])))

print(f'Took {etime-btime} seconds for {numyrs} years, or {(etime-btime)/numyrs} per year')

#%%
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import matplotlib.ticker as tkr
import datetime as dt
import numpy as np

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

outdir = os.path.join(kwk.ProjDir, 'outputs')
if not os.path.exists(outdir):
    os.mkdir(outdir)
    
figdir = os.path.join(outdir, 'figs')
if not os.path.exists(figdir):
    os.mkdir(figdir)
todaystr = dt.date.today().strftime('%Y%m%d')
#%% make release temperature comparison plots
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

#%%
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
                                     f'{todaystr}DRAFT_Keswick_ReleaseTempsSim_allyears.png'),
                        dpi=600) # kwargs)
        else:
            plt.savefig(os.path.join(figdir, 
                                     f'{todaystr}DRAFT_Keswick_ReleaseTempsSim_{y}.png'),
                        dpi=600) 
        plt.close()
plt.ion()
#ax.plot(kwk.Inflow.DataFrame['InflowTemp_degF'])
#%% check evaporation
obsevap_af = kwk.Met.DataFrame.loc[kwk.SimDates, ['Evap_AF']]
                           
evap_sim_af = pnd.DataFrame(index=kwk.SimDates,
                            data=[_/lt.AFtoM3 for _ in kwk.Evaporation],
                            columns=['SimEvap_AF'])

yrs = list(set([i.year for i in kwk.SimDates])) +['all']
plt.ioff()
for y in yrs:

    with sns.plotting_context('notebook', font_scale=1.2):
        fig, ax = plt.subplots(1,1, figsize=(10,6))
        if y=='all':
            ax.plot(obsevap_af.loc[:, 'Evap_AF'],c='#3b4d7a', label='Observed')
            #ax.plot(kwk.SimDates, [_/lt.AFtoM3 for _ in kwk.StorageRecord[1:]])
            ax.plot(evap_sim_af.loc[:,'SimEvap_AF'], c='#ed6002', label='Simulated')
            ax.xaxis.set_major_locator(mdates.YearLocator())
            ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[6]))
            ax.xaxis.set_major_formatter(years_fmt)
            ax.xaxis.set_tick_params(rotation=75)
            
        else:   
            ax.plot(obsevap_af.loc[str(y), 'Evap_AF'],c='#3b4d7a', label='Observed')
            #ax.plot(kwk.SimDates, [_/lt.AFtoM3 for _ in kwk.StorageRecord[1:]])
            ax.plot(evap_sim_af.loc[str(y),'SimEvap_AF'], c='#ed6002', label='Simulated')
            ax.set_xlabel('Date')
            ax.xaxis.set_major_locator(months)
            ax.xaxis.set_major_formatter(monyr_fmt)
            
        ax.yaxis.set_minor_locator(tkr.AutoMinorLocator())
        ax.set_ylabel('acre-feet', fontsize=14)
        if y=='all':
            ax.set_title(f'Keswick Evaporation', 
                         fontweight='bold', fontsize=16)            
        else:
            ax.set_title(f'Keswick Evaporation - {y}', 
                         fontweight='bold', fontsize=16)

        plt.subplots_adjust(left=0.21,right=0.98, bottom=0.15)
        ax2 = ax.secondary_yaxis( -0.1, functions=(af_to_mcm, mcm_to_af))
        #ax2.yaxis.set_major_locator(tkr.MultipleLocator(1))
        ax2.yaxis.set_minor_locator(tkr.MultipleLocator(0.1))
        ax2.set_ylabel('million $m^3$',fontsize=14)
        
        ax.annotate('Keswick Evaporation', (0.1, 0.33), xycoords='figure fraction',
                      xytext=(0.02, 0.33), textcoords='figure fraction', 
                      fontsize=14, fontweight='bold', rotation=90) 
        
        sns.despine()
        plt.legend(loc='best', frameon=False)
        plt.savefig(os.path.join(figdir, 
                                 f'{todaystr}DRAFT_Keswick_EvapSim_{y}.png'),
                    dpi=600) # kwargs)
        plt.close()
plt.ion()

#%%
fig, ax = plt.subplots(1,1, figsize=(10,6))
ax.plot(kwk.Met.DataFrame.loc[kwk.SimDates, 'Evap_AF'], label='Obs Evap')
ax.plot(kwk.SimDates, evap_sim_af, label='Sim Evap')
plt.legend(loc='best')
#plt.plot(kwk.SimDates, kwk.Accretions)
#%% plot storage comparisons


sim_sto_DF = pnd.DataFrame(index=kwk.SimDates, 
                           data=kwk.StorageRecord[1:],
                           columns=['SimSto_m3'])
sim_sto_DF['SimSto_AF'] = sim_sto_DF.SimSto_m3/lt.AFtoM3

obs_sto_DF = kwk.Observations.OutflowDF.loc[:,[kwk.Observations.OutflowColMap['obs_storage_final']]]
obs_sto_DF['ObsSto_AF'] = obs_sto_DF[kwk.Observations.OutflowColMap['obs_storage_final']]/lt.AFtoM3
plt.ioff()
for y in years: #range(2000,2023):
    with sns.plotting_context('notebook', font_scale=1.2):
        fig, ax = plt.subplots(1,1, figsize=(10,6))
        
        if y=='all':
            ax.plot(obs_sto_DF.loc[:, 'ObsSto_AF'],c='#3b4d7a', label='Observed')
            #ax.plot(kwk.SimDates, [_/lt.AFtoM3 for _ in kwk.StorageRecord[1:]])
            ax.plot(sim_sto_DF.loc[:,'SimSto_AF'], c='#ed6002', label='Simulated')
            ax.xaxis.set_major_locator(mdates.YearLocator())
            ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[6]))
            ax.xaxis.set_major_formatter(years_fmt)
            ax.xaxis.set_tick_params(rotation=75)
            
        else:
            ax.plot(obs_sto_DF.loc[str(y), 'ObsSto_AF'],c='#3b4d7a', label='Observed')
            #ax.plot(kwk.SimDates, [_/lt.AFtoM3 for _ in kwk.StorageRecord[1:]])
            ax.plot(sim_sto_DF.loc[str(y),'SimSto_AF'], c='#ed6002', label='Simulated') # kwk.SimDates,kwk.StorageRecord[1:])
            ax.xaxis.set_major_locator(months)
            ax.xaxis.set_major_formatter(monyr_fmt)
            ax.set_xlabel('Date')

        
        ax.yaxis.set_minor_locator(tkr.AutoMinorLocator())
        ax.set_ylabel('acre-feet', fontsize=14)
        ax.set_title(f'Keswick Storage - {y}', 
                     fontweight='bold', fontsize=16)

        plt.subplots_adjust(left=0.19, bottom=0.15)
        ax2 = ax.secondary_yaxis( -0.14, functions=(af_to_mcm,mcm_to_af))
        ax2.yaxis.set_minor_locator(tkr.MultipleLocator(0.2))
        ax2.set_ylabel('million $m^3$',fontsize=14)
        
        ax.annotate('Storage Volume', (0.1, 0.33), xycoords='figure fraction',
                      xytext=(0.01, 0.36), textcoords='figure fraction', 
                      fontsize=14, fontweight='bold', rotation=90) 
        
        sns.despine()
        plt.legend(loc='best', frameon=False)
        plt.savefig(os.path.join(figdir, 
                                 f'{todaystr}DRAFT_Keswick_StorageSimObs_{y}.png'),
                    dpi=600) # kwargs)
        plt.close()
        
#%% plot release comparisons

sim_rel_DF = pnd.DataFrame(index=kwk.SimDates, 
                           data=kwk.ReleaseVolumes,
                           columns=['SimRel_m3'])
sim_rel_DF['SimRel_cfs'] = sim_rel_DF.SimRel_m3/lt.AFtoM3*43560/86400 # /lt.AFtoM3

obs_rel_DF = kwk.Observations.OutflowDF.loc[:,[kwk.Observations.OutflowColMap['flow']]]
obs_rel_DF['ObsRel_cfs'] = obs_rel_DF[kwk.Observations.OutflowColMap['flow']]*43560/86400
plt.ioff()
for y in years: #range(2000,2023):
    with sns.plotting_context('notebook', font_scale=1.2):
        fig, ax = plt.subplots(1,1, figsize=(10,6))
        
        if y=='all':
            ax.plot(obs_rel_DF.loc[:, 'ObsRel_cfs'],lw=3, alpha=0.8,
                    c='#3b4d7a', label='Observed')
            
            ax.plot(sim_rel_DF.loc[:,'SimRel_cfs'], 
                    c='#ed6002', label='Simulated') 
            ax.xaxis.set_major_locator(mdates.YearLocator())
            ax.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[6]))
            ax.xaxis.set_major_formatter(years_fmt)
            ax.xaxis.set_tick_params(rotation=75)
        else:
            ax.plot(obs_rel_DF.loc[str(y), 'ObsRel_cfs'],
                    c='#3b4d7a', label='Observed', lw=3, alpha=0.8)
            #ax.plot(kwk.SimDates, [_/lt.AFtoM3 for _ in kwk.StorageRecord[1:]])
            ax.plot(sim_rel_DF.loc[str(y),'SimRel_cfs'], 
                    c='#ed6002', label='Simulated') # kwk.SimDates,kwk.StorageRecord[1:])
            ax.set_xlabel('Date')
            ax.xaxis.set_major_locator(months)
            ax.xaxis.set_major_formatter(monyr_fmt)
        
        ax.yaxis.set_minor_locator(tkr.AutoMinorLocator())
        ax.set_ylabel('$ft^3/s$', fontsize=14)
        if y=='all':
            titlestr = 'Keswick Total Releases'
        else:
            titlestr = f'Keswick Total Releases - {y}'
        ax.set_title(titlestr, fontweight='bold', fontsize=16)

        plt.subplots_adjust(left=0.25, bottom=0.15, right=0.98)
        ax2 = ax.secondary_yaxis( -0.15, functions=(cfs_to_cms, cms_to_cfs))
        ax2.yaxis.set_minor_locator(tkr.AutoMinorLocator())
        ax2.set_ylabel('$m^3/s$',fontsize=14)
        
        ax.annotate('Total Daily Release Rate', (0.1, 0.33), xycoords='figure fraction',
                      xytext=(0.01, 0.32), textcoords='figure fraction', 
                      fontsize=14, fontweight='bold', rotation=90) 
        
        sns.despine()
        plt.legend(loc='best', frameon=False)
        plt.savefig(os.path.join(figdir, 
                                 f'{todaystr}DRAFT_Keswick_ReleaseSimObs_{y}.png'),
                    dpi=600) # kwargs)
        plt.close()
plt.ion()
#%%
# def af_to_mcm(afv):
#     mcmv = afv/1e6*43560.*(12**3)*(2.54**3)/(100**3)
#     return(mcmv)
# def mcm_to_af(mcmv):
#     mafv = mcmv*(100**3)/(2.54**3)/(12**3)/43560
#     return(mafv*1e6)

# simSto = pnd.DataFrame(kwk.StorageRecord[1:], index=kwk.SimDates)
# obsSto = pnd.DataFrame(kwk.Observations.OutflowDF['KES_Storage_AF'])
# simStoAF = pnd.DataFrame([mcm_to_af(i/1e6) for i in kwk.StorageRecord[1:]], index=kwk.SimDates)
# fig, ax = plt.subplots(1,1, figsize=(11,6))
# ax.plot(obsSto,c='#3b4d7a', label='Observed')
# ax.plot(simStoAF,c='#ed6002', label='Simulated',alpha=0.7 )


#%%

simRel = pnd.Series([mcm_to_af(i/1e6) for i in kwk.Outflow.DataFrame.loc[simT.index,'OutflowTotal']]) # 'KES_Outflow_af']

sha_kwk_diff = pnd.DataFrame(index=kwk.SimDates) 
sha_kwk_diff['SimT'] = simT
sha_kwk_diff['ObsT'] = obsT
sha_kwk_diff['ShaTobs'] =kwk.Inflow.DataFrame.loc[simT.index, 'InflowTemp_degC'] #'SHD_OutflowTemp_degF']
sha_kwk_diff['Kwk_min_Sha'] = sha_kwk_diff.SimT.iloc[0:] - sha_kwk_diff.ShaTobs.iloc[0:] 

sha_kwk_diff['Kwk_min_Sha_OBS'] = sha_kwk_diff.ObsT.iloc[0:] - sha_kwk_diff.ShaTobs.iloc[0:]
sel_months = np.isin(sha_kwk_diff.index.month, [6,7,8,9])

with sns.plotting_context('notebook', font_scale=1.25):
    
    fig, ax = plt.subplots(1,1, figsize=(10,8))
    ax.scatter(simRel[sel_months]*43560/86400, 
                sha_kwk_diff.loc[sel_months,'Kwk_min_Sha_OBS'], label='Obs')
    ax.scatter(simRel[sel_months]*43560/86400, 
            sha_kwk_diff.loc[sel_months,'Kwk_min_Sha'], label='Sim')

    ax.xaxis.set_minor_locator(tkr.AutoMinorLocator())
    ax.yaxis.set_minor_locator(tkr.AutoMinorLocator())
    
    plt.subplots_adjust(left=0.19)
    ax2b = ax.secondary_yaxis( -0.15, functions=(c_to_f_diff, f_to_c_diff))
    ax2b.set_ylabel('Temperature Difference Through Keswick, $^oF$')
    ax2b.yaxis.set_minor_locator(tkr.AutoMinorLocator())
    ax.set_xlabel('Keswick Release, $ft^3/s$')
    ax.set_ylabel('$^oC$')
    ax.set_title('Keswick Warming as a Function of Keswick Releases\nJune-September',
                 fontweight='bold')
    sns.despine()
    plt.legend(loc='best', frameon=False)
 
    plt.savefig(os.path.join(figdir, 
                             f'{todaystr}DRAFT_Keswick_Warming_JunSep_SimObs.png'),
                dpi=600) # kwargs)
    plt.close()
# fig, ax = plt.subplots(1,1, figsize=(10,6))
# ax.plot(obsSto)
# ax.plot(simSto)

#%% spring creek tunnel contributions

spp_obs = kwk.Inflow.DataFrame.loc['2000':'2022', 
                                   [kwk.Inflow.ColumnMap['tribInflow_final'],
                                    kwk.Inflow.ColumnMap['tribInflowTemp_final']]]
spp_obs['Kwk_outflow'] = kwk.Outflow.DataFrame.loc['2000':'2022', 
                                                   [kwk.Outflow.ColumnMap['outflow_final']]]
spp_obs['SPP_frac_outflow'] = spp_obs['tribInflow_final']/spp_obs['Kwk_outflow']


sel_months = np.isin(sha_kwk_diff.index.month, [6,7,8,9])         
dryyrs_summer = np.where(sha_kwk_diff.index.year.isin([2014, 2015, 2021, 2022]) &
                 sha_kwk_diff.index.month.isin([6,7,8,9]))                          
# plt.scatter(spp_obs.loc[sel_months, 'tribInflow_final'], 
#             obsT.loc[sel_months,obsT.columns[0]],
#             c=spp_obs.loc[sel_months,'tribInflowTemp_degC'], cmap='RdYlBu_r')

sel_spp_obs = spp_obs.loc[np.isin(spp_obs.index.year, np.arange(2000,2023))] #[2012, 2013,2014, 2015,2016, 2017, 2018, 2019, 2020, 2021, 2022]),:]
sel_spp_obs = sel_spp_obs.loc[np.isin(sel_spp_obs.index.month,[4,5,6]),:]
with sns.plotting_context('notebook', font_scale=1.2):
    fig, ax = plt.subplots(1,2, figsize=(13, 6))
    sc0 = ax[0].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
                obsT.loc[sel_spp_obs.index,obsT.columns[0]],
                c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')
    ax[0].set_xlabel('Fraction of Keswick release from Spring Creek PP', 
                     fontsize=14)
    ax[0].xaxis.set_minor_locator(tkr.MultipleLocator(0.05))
    ax[0].set_ylabel('Keswick Release Temperature, $^oC$')
    
    
    sc1 =ax[1].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
                simT.loc[sel_spp_obs.index,simT.columns[0]],
                c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')
    ax[1].set_xlabel('Fraction of Keswick release from Spring Creek PP', 
                     fontsize=14)
    ax[1].set_ylabel('Keswick Release Temperature, $^oC$')
    ax[1].set_ylim(ax[0].get_ylim())
    sns.despine()
    ax[0].set_title("Observed, Apr-May-Jun", fontstyle='italic' )
    ax[1].set_title("Simulated, Apr-May-Jun", fontstyle='italic' )
    
    plt.colorbar(sc0, ax=ax[0], label='SPP Temperature, $^oC$',
                 pad=0.1)
    plt.colorbar(sc1, ax=ax[1], label='SPP Temperature, $^oC$',
                 pad=0.1)
    
    plt.subplots_adjust(left=0.07,right=0.95, top=0.85, wspace=0.2)
    plt.suptitle('Keswick Release Temperatures as a\nFunction of Spring Creek PP Releases and Temperature',
                 fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(figdir, 
                             f'{todaystr}DRAFT_KeswickTemp_vs_SPPfrac_AprJun_SimObs.png'),
                dpi=600) # kwargs)
    plt.close()

sel_spp_obs = spp_obs.loc[np.isin(spp_obs.index.year, np.arange(2000,2023))] #[2012, 2013,2014, 2015,2016, 2017, 2018, 2019, 2020, 2021, 2022]),:]
sel_spp_obs = sel_spp_obs.loc[np.isin(sel_spp_obs.index.month,[7,8,9]),:]
with sns.plotting_context('notebook', font_scale=1.2):
    fig, ax = plt.subplots(1,2, figsize=(13, 6))
    sc0 = ax[0].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
                obsT.loc[sel_spp_obs.index,obsT.columns[0]],
                c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')
    ax[0].set_xlabel('Fraction of Keswick release from Spring Creek PP', 
                     fontsize=14)
    ax[0].xaxis.set_minor_locator(tkr.MultipleLocator(0.05))
    ax[0].set_ylabel('Keswick Release Temperature, $^oC$')
    
    
    sc1 =ax[1].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
                simT.loc[sel_spp_obs.index,simT.columns[0]],
                c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')
    ax[1].set_xlabel('Fraction of Keswick release from Spring Creek PP', 
                     fontsize=14)
    ax[1].set_ylabel('Keswick Release Temperature, $^oC$')
    ax[1].set_ylim(ax[0].get_ylim())
    sns.despine()
    ax[0].set_title("Observed, Jul-Sep", fontstyle='italic' )
    ax[1].set_title("Simulated, Jul-Sep", fontstyle='italic' )
    
    plt.colorbar(sc0, ax=ax[0], label='SPP Temperature, $^oC$',
                 pad=0.1)
    plt.colorbar(sc1, ax=ax[1], label='SPP Temperature, $^oC$',
                 pad=0.1)
    
    plt.subplots_adjust(left=0.07,right=0.95, top=0.85, wspace=0.2)
    plt.suptitle('Keswick Release Temperatures as a\nFunction of Spring Creek PP Releases and Temperature',
                 fontsize=16, fontweight='bold')
    plt.savefig(os.path.join(figdir, 
                             f'{todaystr}DRAFT_KeswickTemp_vs_SPPfrac_JulSep_SimObs.png'),
                dpi=600) # kwargs)
    plt.close()
# sel_spp_obs = spp_obs.loc[np.isin(spp_obs.index.year, np.arange(2000,2023))] #[2012, 2013,2014, 2015,2016, 2017, 2018, 2019, 2020, 2021, 2022]),:]
# sel_spp_obs = sel_spp_obs.loc[np.isin(sel_spp_obs.index.month,[7,8,9]),:]
# fig, ax = plt.subplots(2,1, figsize=(11, 7))
# ax[0].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
#             obsT.loc[sel_spp_obs.index,obsT.columns[0]],
#             c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')
# ax[1].scatter(sel_spp_obs.loc[:, 'SPP_frac_outflow'], #'tribInflow_final'
#             simT.loc[sel_spp_obs.index,simT.columns[0]],
#             c=sel_spp_obs.loc[:,'tribInflowTemp_degC'], cmap='RdYlBu_r')

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
    
run_name = kwk.RunName 
PROJ_DIR = kwk.ProjDir
outDir = outdir #<-- defined above for figures os.path.join(PROJ_DIR, 'outputs')
if not os.path.exists(outDir):
    os.mkdir(outDir)

stats_dir = os.path.join(PROJ_DIR, 'outputs', 'stats')
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