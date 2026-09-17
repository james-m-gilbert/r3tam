# -*- coding: utf-8 -*-
"""
Simulation of Sacramento River segments using regression-based water temperature
model. This is a standalone simulation using historical observations as inputs.

Inputs and results align with the content described in Section 3.3 of the 2026
R3TAM report.

@author: jgilbert
"""
import os
import datetime as dt
import pandas as pnd
import matplotlib.pyplot as plt

from r3tam.river import River
import time 

#%% test one

config_fp = r'./uppsac_river_kwk2bnd.yaml'
uppsac = River.initialize_model(config_fp)


#uppsac.calcRivTemps('ccr')

btime = time.perf_counter()
for d in uppsac.SimDates:
    uppsac.advance_rivmod(final=True)
etime = time.perf_counter()
print(f'Took {etime-btime} seconds')

#fig = uppsac.plot_sim_obs('ccr', year=2020, plot_units='degC')



#%% post-process and plot

outdir = os.path.join(os.path.abspath(os.curdir), 'outputs')
figdir = os.path.join(outdir,'figs')

if not os.path.exists(outdir):
    os.mkdir(outdir)

if not os.path.exists(figdir):
    os.mkdir(figdir)
todaystr = dt.date.today().strftime('%Y%m%d')

df = uppsac.combine_to_df()


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

f = uppsac.plot_sim_obs('ccr',plot_units='degC')
f.tight_layout()
ofn = f'DRAFT_SimObsRivTemp_CCR_allyears.png'
f.savefig(os.path.join(figdir, ofn), dpi=600)
plt.close()

f = uppsac.plot_sim_obs('bsf',plot_units='degC')
f.tight_layout()
ofn = f'DRAFT_SimObsRivTemp_BSF_allyears.png'
f.savefig(os.path.join(figdir, ofn), dpi=600)
plt.close()

f = uppsac.plot_sim_obs('bnd',plot_units='degC')
f.tight_layout()
ofn = f'DRAFT_SimObsRivTemp_BND_allyears.png'
f.savefig(os.path.join(figdir, ofn), dpi=600)
plt.close()

f = uppsac.plot_sim_obs('rdb',plot_units='degC')
f.tight_layout()
ofn = f'DRAFT_SimObsRivTemp_RDB_allyears.png'
f.savefig(os.path.join(figdir, ofn), dpi=600)
plt.close()
#%%  write out error metrics


run_name = uppsac.RunName
# PROJ_DIR = os.path.dirname(uppsac.ConfigFP)
# outDir = os.path.join(PROJ_DIR, 'outputs')
# if not os.path.exists(outDir):
#     os.mkdir(outDir)

stats_dir = os.path.join(outdir, 'stats')
if not os.path.exists(stats_dir):
    os.mkdir(stats_dir)

ccr_err = uppsac.calc_error_metrics(df, 'ccr') #, monthly=True)      
#ccr_ann_err = uppsac.calc_error_metrics(df, 'ccr', monthly=False)   

bsf_err = uppsac.calc_error_metrics(df, 'bsf')      
#bsf_ann_err = uppsac.calc_error_metrics(df, 'bsf', monthly=False)  

bnd_err = uppsac.calc_error_metrics(df, 'bnd')

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
