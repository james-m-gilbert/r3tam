# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 11:39:51 2022

@author: jgilbert
"""
import os
import datetime as dt
import pandas as pnd
import matplotlib.pyplot as plt

from river import River
import time 

#%% test one
config_fp = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/river_example/uppsac_river_kwk2ccr.yaml'
#config_fp = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/river_example/uppsac_river_kwk2bnd.yaml'
uppsac = River.initialize_model(config_fp)

#uppsac.calcRivTemps('ccr')

btime = time.perf_counter()
for d in uppsac.SimDates:
    uppsac.advance_rivmod(final=True)
etime = time.perf_counter()
print(f'Took {etime-btime} seconds')

fig = uppsac.plot_sim_obs('ccr', year=2020, plot_units='degC')



#%% test manual updating of keswick boundary temperatures
uppsac2 = River.initialize_model(config_fp)
kwkobs = uppsac2.Nodes['kwk'].BoundaryConditions.DataFrame.loc[uppsac2.SimDates, 'KWK_WatTemp_degC']
kwkobs_mod = pnd.Series([_+np.random.normal(0, 1.5) for _ in kwkobs], index=kwkobs.index)
kwkobsq = uppsac2.Nodes['kwk'].BoundaryConditions.DataFrame.loc[uppsac2.SimDates, 'KWK_Flow_cms']
btime = time.perf_counter()
for d in uppsac2.SimDates:
    uppsac2.set_node_flow_temp('kwk', kwkobs_mod.loc[d], kwkobsq.loc[d])
    uppsac2.advance_rivmod()
    #print(f'simulated date {d}')
etime = time.perf_counter()

print(f'Took {etime-btime} seconds')

uppsac2.plot_sim_obs('ccr')

#%% test running down to red bluff
config_fp = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/river_example/uppsac_river_kwk2bnd.yaml'

uppsac = River.initialize_model(config_fp)

#uppsac.calcRivTemps('ccr')

btime = time.perf_counter()
for d in uppsac.SimDates:
    uppsac.advance_rivmod(final=True)
etime = time.perf_counter()
print(f'Took {etime-btime} seconds')

#%% post-process and plot
df = uppsac.combine_to_df()

figdir = os.path.join(os.path.dirname(config_fp),'outputs','figs')
plt.ioff()
for y in range(2000, 2023):
    for n in uppsac.Nodes:
        if uppsac.Nodes[n].UpNode>=0:
            f = uppsac.plot_sim_obs(n, plot_units='degC', year=y)
            f.tight_layout()
            ofn = f'DRAFT_SimObsRivTemp_{n.upper()}_{y}.png'
            f.savefig(os.path.join(figdir, ofn), dpi=600)
            plt.close()
plt.ion()          

#%%  write out error metrics


run_name = uppsac.RunName
PROJ_DIR = os.path.dirname(uppsac.ConfigFP)
outDir = os.path.join(PROJ_DIR, 'outputs')
if not os.path.exists(outDir):
    os.mkdir(outDir)

stats_dir = os.path.join(PROJ_DIR, 'outputs', 'stats')
if not os.path.exists(stats_dir):
    os.mkdir(stats_dir)

ccr_err = uppsac.calc_error_metrics(df, 'ccr') #, monthly=True)      
#ccr_ann_err = uppsac.calc_error_metrics(df, 'ccr', monthly=False)   

bsf_err = uppsac.calc_error_metrics(df, 'bsf')      
#bsf_ann_err = uppsac.calc_error_metrics(df, 'bsf', monthly=False)  

rdb_err = uppsac.calc_error_metrics(df, 'rdb')      
#rdb_ann_err = uppsac.calc_error_metrics(df, 'rdb', monthly=False)    

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
#%% make plots by year



# resmod = rt.Res.initialize_model(config_fp, profile_temp_units='degF')
# resmod.Debug['Release'] = 0
# btime = time.perf_counter() #process_time()
# for d in resmod.SimDates:
#     print(d)
#     resmod.advance_restemp()

#     resmod.advance_swd(final=True)

# etime = time.perf_counter() #process_time()
# print(f"took {etime-btime} seconds?")