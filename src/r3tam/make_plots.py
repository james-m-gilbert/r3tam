# -*- coding: utf-8 -*-
"""
make plots

Created on Mon Jul  6 14:58:55 2020

@author: jgilbert
"""
import os, sys
import pickle
import numpy as np

import datetime as dt
import copy

import pandas as pnd
idx = pnd.IndexSlice

import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('ticks')

from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()

from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                               AutoLocator, AutoMinorLocator)

import matplotlib.dates as mdates
years = mdates.YearLocator()   # every year
fiveyears = mdates.YearLocator(5)
months = mdates.MonthLocator()  # every month
pltdays = mdates.DayLocator((5,10,15,20, 25))
years_fmt = mdates.DateFormatter('%Y')
months_fmt = mdates.DateFormatter('%b')
monyr_fmt = mdates.DateFormatter('%b\n%Y')

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
# def cfs_to_cms(cfs):
#     return(cfs*lt.FT3toM3)
# def cms_to_cfs(cms):
#     return(cms/lt.FT3toM3)


def plotStorageCompare(resObj, viewSave='save', on_wy=False):
    
    runName = resObj.RunName 
    figoutDir = os.path.join(resObj.ProjDir, 'outputs', 'figs')
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)

    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
    
    # input time series
    obsVals = resObj.Observations.OutflowDF
    obsHdr = resObj.Observations.OutflowColMap
    
    obsStorage = obsVals.loc[:,[obsHdr['obs_storage_final']]]
    
    simStoragesDF = resObj.Simulation_Results['StorageDF']
    
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simsto_wys = simStoragesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simsto_wys)
    else:
        yrs = np.unique(simStoragesDF.index.map(lambda x: x.year))
    

    for y in yrs:
        
        if on_wy:
            thisDateList = simStoragesDF[simsto_wys==y].index

        else:
            thisDateList = simStoragesDF[str(y)].index        

        this_obs_sto = obsStorage.loc[thisDateList,:]
        this_sim_sto = simStoragesDF.loc[thisDateList,:]

        # plot storage
        with sns.plotting_context('paper', font_scale=1.5):
            fig, ax = plt.subplots(2,1, figsize=(9,6))
            #ax.plot(inpDF_AFF.index[0:len(storages)], storages, color='orange', label='SimpleModel')
            #ax.plot(inpDF_AFF.index, inpDF_AFF.ShastaStorage_AF, color='k', label='Observed Shasta Storage')
            # ax[0].plot(simStoragesDF, color='orange', label='SimpleModel')
            # ax[0].plot(obsStorage.index, obsStorage.values, color='k', label='Observed Shasta Storage')
        
            # ax[1].plot(simStoragesDF.index, simStoragesDF.iloc[:,0] - \
            #           obsStorage.loc[simStoragesDF.index], 
            #           color='k', label='Storage Error')
                
            ax[0].plot(thisDateList, this_sim_sto, color='orange', label=f'ResTemp: {runName}')
            ax[0].plot(thisDateList, this_obs_sto, color='k', label='Observed Shasta Storage')
        
            ax[1].plot(thisDateList, this_sim_sto.values - this_obs_sto.values,
                      color='k', label='Storage Error')
            if viewSave=='save':
                if on_wy:
                    plt.savefig(os.path.join(figoutDir, f'{runName}_Storage_WY{y}.png'), dpi=300)
                else:
                    plt.savefig(os.path.join(figoutDir, f'{runName}_Storage_{y}.png'), dpi=300)
    
                plt.close()

def plotStorageEvapCompare(resObj, viewSave='save', on_wy=False,**kwargs):
    
    runName = resObj.RunName 
    figoutDir = os.path.join(resObj.ProjDir, 'outputs', 'figs')
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)


    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
    
    # input time series
    obsVals = resObj.Observations.OutflowDF
    obsHdr = resObj.Observations.OutflowColMap
    
    obsStorage = obsVals.loc[:,[obsHdr['obs_storage_final']]]
    
    simStoragesDF = resObj.Simulation_Results['StorageDF']
    
    
    if not resObj.SimulationSpecs.CalcEvaporation:
        print("\n****************************************************")
        print("  Evaporation was provided as an input - nothing to plot")
        print("****************************************************\n")
        
    else:
        if 'evap' not in resObj.Met.ColumnMap:
            print(" --- Don't have evaporation observations to compare to,\n --- try adding a column for evaporation to the input time series")
        else:
            evapCol = resObj.Met.ColumnMap['evap']
        
        evapObs = resObj.Met.DataFrame.loc[:,evapCol]
        evapSim = resObj.Simulation_Results['EvapDF']
    
    
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simsto_wys = simStoragesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simsto_wys)
    else:
        yrs = np.unique(simStoragesDF.index.map(lambda x: x.year))
    
    if 'years' in kwargs:
        selyrs = kwargs['years']
        yrs = [y for y in yrs if y in selyrs]
        
    for y in yrs:
        
        if on_wy:
            thisDateList = simStoragesDF[simsto_wys==y].index

        else:
            thisDateList = simStoragesDF[str(y)].index        

        this_obs_sto = obsStorage.loc[thisDateList,:]
        this_sim_sto = simStoragesDF.loc[thisDateList,:]

        this_obs_evap = evapObs.loc[thisDateList]*86400/43560. # obs evap in cfs - convert to AF
        this_sim_evap = evapSim.loc[thisDateList,:]
            
                
        # plot storage
        with sns.plotting_context('notebook', font_scale=1.2):
            fig, ax = plt.subplots(2,1, figsize=(12,8))
            #ax.plot(inpDF_AFF.index[0:len(storages)], storages, color='orange', label='SimpleModel')
            #ax.plot(inpDF_AFF.index, inpDF_AFF.ShastaStorage_AF, color='k', label='Observed Shasta Storage')
            # ax[0].plot(simStoragesDF, color='orange', label='SimpleModel')
            # ax[0].plot(obsStorage.index, obsStorage.values, color='k', label='Observed Shasta Storage')
        
            # ax[1].plot(simStoragesDF.index, simStoragesDF.iloc[:,0] - \
            #           obsStorage.loc[simStoragesDF.index], 
            #           color='k', label='Storage Error')
                
            ax[0].plot(thisDateList, this_obs_sto/1000., color='k', label='Observed')
            ax[0].plot(thisDateList, this_sim_sto/1000., color='orange', label='Simulated')
            ax[0].legend(loc='best', frameon=False)
            ax[0].set_ylabel('Reservoir Storage\nthousand acre-feet')
            ax[0].xaxis.set_major_locator(months)
            ax[0].xaxis.set_major_formatter(monyr_fmt)
            ax[0].yaxis.set_minor_locator(AutoMinorLocator())
            
            ax[1].plot(this_obs_evap.index, this_obs_evap , color='k', 
                       label='Observed Evap')
            ax[1].plot(this_sim_evap.index, this_sim_evap.loc[:,this_sim_evap.columns[0]], 
                    color='orange', label='Simulated Evap')
            
            ax[1].set_xlabel('Date')
            ax[1].set_ylabel('Evaporation Volume\nacre-feet')
            ax[1].xaxis.set_major_locator(months)
            ax[1].xaxis.set_major_formatter(monyr_fmt)
            ax[1].yaxis.set_minor_locator(AutoMinorLocator())
            ax[1].legend(loc='best', frameon=False)
            
            plt.suptitle(f'Shasta Storage & Evaporation\nSimulated-Observed Comparison, {y}', fontweight='bold')
            
            plt.subplots_adjust(hspace=0.2)
            sns.despine()
            
            if viewSave=='save':
                if on_wy:
                    plt.savefig(os.path.join(figoutDir, f'{runName}_Storage_WY{y}.png'), dpi=300)
                else:
                    plt.savefig(os.path.join(figoutDir, f'{runName}_Storage_{y}.png'), dpi=300)
    
                plt.close()

def plotReleasesCompare(resObj,select_years='all', viewSave='save',
                        on_wy=True, predicted_gates=False, target_ts=None, 
                        target_tol=None,**kwargs):
    
    
    runName = resObj.RunName 
    PROJ_DIR = resObj.ProjDir
    figoutDir = os.path.join(PROJ_DIR, 'outputs', 'figs')
    
    # input time series
    simReleaseDF = resObj.Simulation_Results['ReleaseDF']
    sim_rel_vol_col = [c for c in simReleaseDF.columns if 'Sim_Release' in c]
    release_inputs = resObj.Outflow.DataFrame
    release_inputs_colmap = resObj.Outflow.ColumnMap
    obscm = resObj.Observations.OutflowColMap
    
   
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)
        
    if 'ylimits' in kwargs:
        ylimits = kwargs['ylimits']
    else:
        ylimits = None
        
    if 'extra_file_name_part' in kwargs:
        extrfnpt = kwargs['extra_file_name_part']
    else:
        extrfnpt=''
        
    if 'sim_label' in kwargs:
        sim_label = kwargs['sim_label']
    else:
        sim_label='Model'
    
    if 'obs_label' in kwargs:
        obs_label = kwargs['obs_label']
    else:
        obs_label = 'Obs'
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simrels_wys = simReleaseDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simrels_wys)
    else:
        yrs = np.unique(simReleaseDF.index.map(lambda x: x.year))

    if isinstance(select_years, list):
        yrs = [y for y in yrs if y in select_years]
        
    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
        
    obsRelease = resObj.Observations.OutflowDF
    obsHdr = resObj.Observations.OutflowColMap
    
    gateNames = [n for n in release_inputs_colmap['gates'].values()]
#    allgateDict = {}
#    for k, v in inTScm['gates'].items():
#        allgateDict[k] = inTS[v]
    
    ## plot releases
    #simReleaseT = simReleasesDF['ReleaseTemp_F']
    #simReleaseQ = simReleasesDF['Release_AF']
    obstempcol = obsHdr['temperature']
    
    obsRelease.loc[(obsRelease[obstempcol]<42.),obstempcol] = np.nan
    
    for y in yrs:
        
        if on_wy:
            thisDateList = simReleaseDF[simrels_wys==y].index

        else:
            thisDateList = simReleaseDF.loc[str(y)].index
        
        
        print("Making plot for year %d" %y)
        #thisDateList = simReleasesDF[str(y)].index
        
        obsQ = obsRelease.loc[thisDateList, obsHdr['flow']]
        obsT = obsRelease.loc[thisDateList, obsHdr['temperature']]
        
        if 'other_temp' in kwargs:
            obsT2 = obsRelease.loc[thisDateList, obsHdr['temperature_other']]
            plot_otherT = True
            if 'other_label' in kwargs:
                other_label = kwargs['other_label']
            else:
                other_label = 'other_obs'
        else:
            plot_otherT= False
            
        # make a time series of target temperatures
        addTarg = False
        addTol = False
        #target_ts = resObj.Simulation_Results['ReleaseTempTarget']
        #targTol = resObj.Simulation_Results['TempTargetTol']
        if type(target_ts) == list:
            targTS = pnd.DataFrame(target_ts, index=thisDateList)
            targTS = targTS.replace(-9999.0, np.nan)
            if sum(pnd.isna(targTS))<len(targTS):
                addTarg=True
        if type(target_ts) ==str: # then assume it's a field in teh inTS
            targTS = simReleaseDF.loc[thisDateList, target_ts]
            targTS = targTS.replace(-9999.0, np.nan)
            if sum(pnd.isna(targTS))<len(targTS):
                addTarg=True
        if type(target_tol)==list:
            targTol = pnd.DataFrame(target_tol, index=thisDateList)
            targTol = targTol.replace(-9999.0, np.nan)
            if sum(pnd.isna(targTol))<len(targTol):
                addTol=True

        
        # make a time series of gate changes
        gateTS = release_inputs.loc[thisDateList, gateNames] #['Upp','Mid','PRG','SDG']
        gateTS2 = gateTS.T
        gateTS2.index = [4, 3, 2, 1]
        gateTS2 = gateTS2.sort_index()
        
        # get simulated gate changes if predicted_gates is True and the GateOps data exists
        # in the reservoir object
        if (predicted_gates) & (len(resObj.Operations.GateOps)>0):
            sim_gate_dict = resObj.Operations.GateOps
            sim_gateTS = pnd.DataFrame.from_dict(sim_gate_dict)
            sim_gateTS = sim_gateTS.T.loc[thisDateList,:].T
            sim_gateTS.index = [4,3,2,1]
            sim_gateTS = sim_gateTS.sort_index()
            
            
        with sns.plotting_context('paper', font_scale=1.2):
            ann_x = 0.97
            ann_y = 0.92
            if (predicted_gates) & (len(resObj.Operations.GateOps)>0):
                fig, ax = plt.subplots(nrows=4, ncols=1, figsize=(11, 10), sharex=True)  # add an extra row for the predicted gate ops
            else:
                fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(11, 8), sharex=True)
            
            #ax[0].plot(inpDF_AFF.index[0:len(simReleaseQ)], simReleaseQ, color='orange', label='SimpleModel') 
            #ax[0].plot(inpDF_AFF.index[0:len(simReleaseQ)], inpDF_AFF.Outflow_AF.iloc[0:len(simReleaseQ)], color='0.37', lw=0.7, label='Obs Shasta Release')
        
            #ax[0].plot(simReleasesDF[str(y)]['Release_AF'], color='orange', label='SimpleModel') 
            ax[0].plot(simReleaseDF.loc[thisDateList, sim_rel_vol_col[0]], 
                       color='darkorange', label=sim_label) 
            ax[0].plot(obsQ, color='0.37', lw=0.7, label="Observed")
        
            ax[0].set_title('Release Volume Comparison')
            ax[0].set_ylabel('Daily Mean Release Volume\n(acre-feet)', fontsize=11)
            ax[0].legend(loc='best', frameon=False, fontsize=12)
            ax[0].annotate('a', (ann_x, ann_y), xycoords='axes fraction',
                          xytext=(ann_x, ann_y*0.6), textcoords='axes fraction', 
                          fontsize=18 )
        #    ax[1].plot(inpDF_AFF.index[0:len(simReleaseT)], simReleaseT, color='orange', label='SimpleModel')
        #    ax[1].plot(inpDF_AFF.index[0:len(simReleaseQ)], inpDF_AFF.OutflowTemp_degF.iloc[0:len(simReleaseQ)], color='0.37',lw=0.7, label='CDEC SHD - Daily Mean')
            #ax[1].plot(simReleasesDF[str(y)]['ReleaseTemp_F'], color='orange', label='SimpleModel')
            ax[1].plot(simReleaseDF.loc[thisDateList,sim_rel_vol_col[1]], 
                       color='darkorange', label=sim_label, lw=1.2)
            ax[1].plot(obsT, color='0.37',lw=1, label=obs_label)
            
            if plot_otherT:
                ax[1].plot(obsT2, color='darkred', lw=0.7, ls='--',label=other_label)
            
            if addTarg:
                ax[1].plot(targTS, color='k', lw=0.7, ls='--', label='Target Temperature')
                if addTol:
                    ax[1].fill_between(targTol.index, targTS.iloc[:,0]+targTol.iloc[:,0],
                                      y2=targTS.iloc[:,0]-targTol.iloc[:,0],
                                      facecolor='0.3', alpha=0.3)
                    
                    
            ax2 = ax[1].secondary_yaxis( -0.075, functions=(f_to_c, c_to_f))
            ax2.yaxis.set_major_locator(AutoLocator()) #tkr.MultipleLocator(1))
            ax2.yaxis.set_minor_locator(AutoMinorLocator()) #tkr.MultipleLocator(0.2))
            ax2.set_ylabel('$^oC$',fontsize=14)
            ax[1].set_title('Release/Tailwater Temperature Comparison')
            ax2.set_ylabel('Daily Mean Water Temperature, $^oC$', fontsize=11)
            ax[1].set_ylabel('$^oF$', fontsize=13)
            ax[1].legend(loc='best', frameon=False)
            if ylimits!=None:
                ax[1].set_ylim((ylimits[0],ylimits[1]))
            #ax[1].set_ylim((40,65))
            sns.despine()
            ax[1].annotate('b', (ann_x, ann_y), xycoords='axes fraction',
                          xytext=(ann_x, ann_y), textcoords='axes fraction', 
                          fontsize=18 )
            
            plt.suptitle('Shasta Releases Comparison - %s' %y, fontsize=16, fontweight='bold')
            
            # plot based on https://stackoverflow.com/questions/24163313/how-to-create-an-activity-plot-from-pandas-dataframe-like-the-github-contributi
            dayidx2 = pnd.date_range(start=gateTS2.columns[0], end=gateTS2.columns[-1]+dt.timedelta(1), freq='D')
            gateLevel, dayidx = np.mgrid[:gateTS2.shape[0]+1, :gateTS2.shape[1]+1]
            ax[2].set_aspect(9.3) #"equal")
            ax[2].pcolormesh(dayidx2,gateLevel , gateTS2.values, 
                             cmap="gist_gray_r",linewidth=0.001, 
                             edgecolor="None",vmin=0, vmax=5)
            ax[2].set_yticks([0.5, 1.5, 2.5, 3.5])
            ax[2].set_yticklabels(['SDG','LOW','MID','UPP'])
            ax[2].set_title('Observed TCD Gate Operations')
            ax[2].annotate('c', (ann_x, ann_y), xycoords='axes fraction',
                      xytext=(ann_x, ann_y), textcoords='axes fraction', 
                      fontsize=18 )
            #plt.xlim(0, gateTS2.shape[1])
            #plt.xlabel(dayidx2)
            
            if (predicted_gates) & (len(resObj.Operations.GateOps)>0):
                dayidx3 = pnd.date_range(start=sim_gateTS.columns[0], end=sim_gateTS.columns[-1]+dt.timedelta(1), freq='D')
                gateLevel, dayidx = np.mgrid[:sim_gateTS.shape[0]+1, :sim_gateTS.shape[1]+1]
                ax[3].set_aspect(9.3) #"equal")
                ax[3].pcolormesh(dayidx3,gateLevel , sim_gateTS.values, 
                                  cmap="gist_gray_r",linewidth=0.001, 
                                  edgecolor="None",vmin=0, vmax=5)
                ax[3].set_yticks([0.5, 1.5, 2.5, 3.5])
                ax[3].set_yticklabels(['SDG','LOW','MID','UPP'])
                ax[3].set_title('Simulated TCD Gate Operations')
                ax[3].annotate('d', (ann_x, ann_y), xycoords='axes fraction',
                          xytext=(ann_x, ann_y), textcoords='axes fraction', 
                          fontsize=18 )
            else:
                if predicted_gates:
                    print("Tried to plot simulated gate ops, but data was not available.")
                    print(" - Try loading in from a *.pkl file instead of the *.yaml")
            
            ax[2].xaxis.set_major_locator(months)
            ax[2].xaxis.set_major_formatter(monyr_fmt)
            plt.tight_layout(h_pad=0.1)
            plt.subplots_adjust(hspace=0.2, top=0.92)
            ax[2].xaxis.set_minor_locator(mdates.DayLocator([7,14,21,28]))
            ax[0].yaxis.set_minor_locator(AutoMinorLocator())
            ax[1].yaxis.set_minor_locator(AutoMinorLocator())
            
        
            if viewSave=='save':
            
                plt.savefig(os.path.join(figoutDir, f'{runName}{extrfnpt}_Outflow_{y}.png'), dpi=300)
                plt.close()
                #return([None, None])
            else:
                return([fig, ax])

def plotProfilesCompare(resObj, ProfilesDF, simStoDF='', on_wy=True, viewSave='save'):

    from matplotlib.lines import Line2D
    
    storElevDF = resObj.ReservoirData.ElevStorage
    
    def elevFromStor(s):
        #return(np.interp(s,storElevDF.Storage_AF,storElevDF.Elevation_ft))
        return(np.interp(s,storElevDF.Storage,storElevDF.Elevation))
        
    def storFromElev(el):
        #tmp = np.interp(el,storElevDF.Elevation_ft,storElevDF.Storage_AF)
        tmp = np.interp(el,storElevDF.Elevation,storElevDF.Storage)
        return(tmp/1e6)
    
    runName = resObj.RunName
    PROJ_DIR = resObj.ProjDir
    #viewSave ='save' 

    outDir = os.path.join(PROJ_DIR, 'outputs')
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    
    figoutDir = os.path.join(PROJ_DIR, 'outputs', 'figs')
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)
    
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simprof_wys = simProfilesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simprof_wys)
    else:
        yrs = np.unique(simProfilesDF.index.map(lambda x: x.year))

    
    modElevs = [l.CtrElev for l in resObj.Layers.values()]
    
    # get observations
    obsProfiles = resObj.Observations.ProfilesDF
    obsHdr = resObj.Observations.ProfilesColMap
    obsTempCol = obsHdr['temperature']
    obsElevs = list(obsProfiles.columns)
    
    if on_wy:
        if obsProfiles.index.nlevels>1:
            obs_wy = obsProfiles.index.get_level_values(0).map(lambda x: x.year+1 if x.month>9 else x.year)
        else:
            obs_wy = obsProfiles.index.map(lambda x: x.year+1 if x.month>9 else x.year)

    # input time series - for marking res storage line
    #inTS = resObj.InputTimeSeries.DataFrame
    inTScm = resObj.Observations.OutflowColMap
    
    obsStorage = resObj.Observations.OutflowDF[inTScm['obs_storage_final']] # inTS.loc[:,inTScm['storage']]
    
    
    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
    
    for y in yrs:    
        with sns.plotting_context('paper', font_scale=0.8):
            
            plt.rc("axes", labelweight="light")
            # get dates of measured profiles for this year
            if on_wy:
                obsProfs = obsProfiles[obs_wy==y]
            else:
                obsProfs = obsProfiles.loc[str(y),:]
            
            obsProfDates = obsProfs.index.get_level_values(0).unique()
            simProfsDates = simProfilesDF.index.unique()
            
            commonDates = sorted(list(set(obsProfDates) & set(simProfsDates)))
            
            simProfs = simProfilesDF.loc[commonDates,:]
            
            [errors, dates, nse, rmse, r2] = calcErrorMetrics(obsProfs, simProfs, weights=1)
            
#            errors = [er for er2 in errors for er in er2 if ~np.isnan(er)]
#        
#            sumsqerr = np.sum([er*er for er in errors])
#                
#            #obs = [o1 for o2 in obs for o1 in o2 if ~np.isnan(o1)]
#                
#            nse = 1 - sumsqerr/(np.nansum((obsProfs-np.nanmean(obsProfs))**2))
#            
#            rmse = np.sqrt(sumsqerr/len(obsProfs))
            #r2 = np.corrcoef(obsProfs, simProfs)[0,1]**2
            print(nse, rmse, r2)
            
            ncols = 8
            if len(obsProfDates)>32:
                nrows = 5
    
            elif len(obsProfDates)>24:
                nrows=4
            elif len(obsProfDates)>16:
                nrows=3
            else:
                ncols = 6
                nrows = 3
            
            fhght = nrows*2.5
            
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=False, figsize=(17,fhght))
             
    
            i = 0
            #for idx,d in obsProfs.iterrows():  #<-- this was used when the obs profiles were in matrix form on regular elevation levels
            for idx in commonDates: 
                thisObsProf = obsProfs.loc[idx]#<-- get the profile for this date
                print(idx)
                #thisSimProf = profiles[list(datelist).index(idx)]
                #thisSimProf = profiles[list(inpDF_AFF.index).index(idx)]
                thisSimProf = simProfilesDF.loc[idx,:]
                nr = i//ncols
                nc = i%ncols
                
                # get the observed (historical) reservoir storage and, optionally, the simulated
                thisObsSto = obsStorage.loc[idx]
                thisObsElev = elevFromStor(thisObsSto)
                if type(simStoDF)==pnd.DataFrame:
                    thisSimSto = simStoDF.loc[idx]
                    thisSimElev = elevFromStor(thisSimSto)
                else:
                    thisSimSto= np.nan
                    
                #print("Obs storage: %0.2f TAF, Obs elevation: %0.2f" %(thisObsSto, thisObsElev))
       #->> this was used with matrix/defined levels version of profile data -->>  #ax[nr, nc].plot([x*1.8 + 32 for x in reversed(d.values)], list(reversed(obsElevs)), color='k', label='Observed')
                ax[nr, nc].plot(thisObsProf.loc[:,obsTempCol], thisObsProf.index, color='k', label='Observed')
                    
                ax[nr, nc].plot(thisSimProf.values,[float(x) for x in thisSimProf.index], color='orange', label='SimpleModel')
                ax[nr, nc].axhline(thisObsElev, lw=0.5, c='k', ls=':')
                #ax[nr, nc].set_xlim((40, 80))S
                ax[nr, nc].plot([43], [thisObsElev+5], marker='v', c='k',linestyle='None',mew=0.2, ms=4, fillstyle='none')
                ax[nr, nc].set_title(idx.strftime('%b %d, %Y'))
                if np.max(thisSimProf.values)>80:
                    ax[nr,nc].set_xticks([40, 50, 60, 70, 80, 90])
                else:
                    ax[nr, nc].set_xticks([40, 50, 60, 70, 80])
                ax[nr, nc].xaxis.set_minor_locator(MultipleLocator(5))
                ax[nr,nc].set_ylim((650, 1070))
                ax[nr, nc].set_yticks([650, 700, 750, 800, 850,  900, 950, 1000, 1050, 1067])
                
                secax = ax[nr, nc].secondary_yaxis('right',functions=(storFromElev,elevFromStor))
                secax.yaxis.set_minor_locator(AutoMinorLocator(0.5))
                
                if (i == len(obsProfs)-1) or (nc == ncols-1):
                    secax.set_ylabel('Storage(MAF)', fontsize=10.)
                else:
                    secax.set_yticklabels([])
                    
                if nc ==0: #== 0: # set equal to zero to jsut label left-hand most column
                    ax[nr,nc].set_ylabel('Elevation (ft)', fontweight='light', fontsize=10.)
                
                numHanging = nrows*ncols - len(obsProfs)
                hangingIndex = ncols - numHanging
                if (nr ==nrows-1) or ((nr == nrows-2) and (nc >= hangingIndex)):
                    ax[nr, nc].set_xlabel('Temperature (deg F)', fontsize=10.)
                i +=1
                
    #        for ix in range(nrows): 
    #            ax[nrows-1,ix].set_xlabel('Temperature, deg F')
    #            ax[ix, 0].set_ylabel('Elevation, ft')
    #        ax[0, nrows-1].set_ylabel('Elevation, ft')
            
            for idx in range(i, nrows*ncols):
                nr = idx//ncols
                nc = idx%ncols
                fig.delaxes(ax[nr, nc])
            sns.despine()
            
            plt.suptitle('Temperature Profile Comparison for Year %d' %(y), fontsize=18)
            
            plt.tight_layout()
            
            if nrows==2:
                plt.subplots_adjust(top=0.76,bottom=0.05, hspace=0.41, wspace=0.6)
                txtheight = 0.835
            else:
                plt.subplots_adjust(top=0.83, bottom=0.05, hspace=0.41, wspace=0.6)
                txtheight = 0.875
            
            legelems = [Line2D([0],[0], color='orange', lw=2, label='Simple Model'),
                        Line2D([0],[0], color='k', lw=2, label='Observation'),
                        Line2D([0],[0], color='k', lw=1, ls=':', label='Obs Elevation')]
            fig.legend(handles=legelems, bbox_to_anchor=(0.65, 0.95),ncol=3, frameon=False, fontsize=12)
            plt.figtext(0.4, txtheight, "NSE: %0.2f    RMSE: %0.2f deg F (%0.2f deg C)" %(nse, rmse, rmse/1.8),
                        fontsize=11, fontfamily='monospace')
            
            if viewSave=='save':
                plt.savefig(os.path.join(figoutDir, '%s_Profiles_%s.png' %(runName, y)), dpi=300)
                plt.close()
                #return([None, None])
            else:
                return([fig, ax])


def plotProfilesCompare2(resObj, simStoDF='', on_wy=True, viewSave='save',
                         show_stats=True, **kwargs):

    from matplotlib.lines import Line2D
    
    storElevDF = resObj.ReservoirData.ElevStorage
    
    ProfilesDF = resObj.Simulation_Results['ProfilesDF']
    
    def elevFromStor(s):
        #return(np.interp(s,storElevDF.Storage_AF,storElevDF.Elevation_ft))
        return(np.interp(s,storElevDF.Storage,storElevDF.Elevation))
        
    def storFromElev(el):
        #tmp = np.interp(el,storElevDF.Elevation_ft,storElevDF.Storage_AF)
        tmp = np.interp(el,storElevDF.Elevation,storElevDF.Storage)
        return(tmp/1e6)
    
    runName = resObj.RunName
    PROJ_DIR = resObj.ProjDir
    #viewSave ='save' 

    outDir = os.path.join(PROJ_DIR, 'outputs')
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    
    figoutDir = os.path.join(PROJ_DIR, 'outputs', 'figs')
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)
    
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simprof_wys = ProfilesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simprof_wys)
    else:
        yrs = np.unique(ProfilesDF.index.map(lambda x: x.year))

    
    modElevs = [l.CtrElev for l in resObj.Layers.values()]
    
    # get observations
    obsProfiles = resObj.Observations.ProfilesDF
    obsHdr = resObj.Observations.ProfilesColMap
    obsTempCol = obsHdr['temperature']
    obsElevs = list(obsProfiles.columns)
    
    if on_wy:
        if obsProfiles.index.nlevels>1:
            obs_wy = obsProfiles.index.get_level_values(0).map(lambda x: x.year+1 if x.month>9 else x.year)
        else:
            obs_wy = obsProfiles.index.map(lambda x: x.year+1 if x.month>9 else x.year)

    # input time series - for marking res storage line
    #inTS = resObj.InputTimeSeries.DataFrame
    inTScm = resObj.Observations.OutflowColMap
    
    obsStorage = resObj.Observations.OutflowDF[inTScm['obs_storage_final']] # inTS.loc[:,inTScm['storage']]
    
    release_inputs = resObj.Outflow.DataFrame
    release_inputs_colmap = resObj.Outflow.ColumnMap
    obscm = resObj.Observations.OutflowColMap
    
    gateNames = [n for n in release_inputs_colmap['gates'].values()]
    rivNames = [n for n in release_inputs_colmap['rivOutlets'].values()]

    
    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
    
    for y in yrs:    
        with sns.plotting_context('paper', font_scale=0.8):
            
            plt.rc("axes", labelweight="light")
            # get dates of measured profiles for this year
            if on_wy:
                obsProfs = obsProfiles[obs_wy==y]
            else:
                obsProfs = obsProfiles.loc[str(y),:]
            
            obsProfDates = obsProfs.index.get_level_values(0).unique()
            simProfsDates = ProfilesDF.index.unique()
            
            commonDates = sorted(list(set(obsProfDates) & set(simProfsDates)))
            
            simProfs = ProfilesDF.loc[commonDates,:]
            
            [errors, dates, nse, rmse, r2] = calcErrorMetrics(obsProfs, 
                                                              simProfs, 
                                                              weights=1)
            
            
            # get gate config on selected profile dates
            gateTS = release_inputs.loc[commonDates, gateNames] #['Upp','Mid','PRG','SDG']
            gateTS2 = gateTS.T
            gateTS2.index = list(range(len(gateTS2)-1,-1,-1))
            gateTS2 = gateTS2.sort_index()
            

            rivoTS = release_inputs.loc[commonDates, rivNames] #['Upp','Mid','PRG','SDG']
            rivoTS2 = rivoTS.T
            rivoTS2.index = list(range(len(rivoTS2),0,-1))
            rivoTS2 = rivoTS2.sort_index()
            
            if 'sim_gates' in kwargs:
                sim_gates = pnd.DataFrame.from_dict(resObj.Operations.GateOps, 
                                                    orient='columns') #,columns=gateNames)
                #sim_gates.columns = gateNames
                gateTS2 = sim_gates.loc[:,commonDates]
                
                # rivoTS2 = pnd.DataFrame().reindex_like(rivoTS2)
                
            
#            errors = [er for er2 in errors for er in er2 if ~np.isnan(er)]
#        
#            sumsqerr = np.sum([er*er for er in errors])
#                
#            #obs = [o1 for o2 in obs for o1 in o2 if ~np.isnan(o1)]
#                
#            nse = 1 - sumsqerr/(np.nansum((obsProfs-np.nanmean(obsProfs))**2))
#            
#            rmse = np.sqrt(sumsqerr/len(obsProfs))
            #r2 = np.corrcoef(obsProfs, simProfs)[0,1]**2
            print(nse, rmse, r2)
            
            ncols = 8
            if len(obsProfDates)>32:
                nrows = 5
    
            elif len(obsProfDates)>24:
                nrows=4
            elif len(obsProfDates)>16:
                nrows=3
            else:
                ncols = 6
                nrows = 3
            
            fhght = nrows*2.5
            
            fig, ax = plt.subplots(nrows=nrows, ncols=ncols, sharex=False, sharey=False, figsize=(17,fhght))
             
    
            i = 0
            #for idx,d in obsProfs.iterrows():  #<-- this was used when the obs profiles were in matrix form on regular elevation levels
            for idx in commonDates: 
                thisObsProf = obsProfs.loc[idx]#<-- get the profile for this date
                print(idx)
                #thisSimProf = profiles[list(datelist).index(idx)]
                #thisSimProf = profiles[list(inpDF_AFF.index).index(idx)]
                thisSimProf = ProfilesDF.loc[idx,:]
                nr = i//ncols
                nc = i%ncols
                
                # get the observed (historical) reservoir storage and, optionally, the simulated
                thisObsSto = obsStorage.loc[idx]
                thisObsElev = elevFromStor(thisObsSto)
                if type(simStoDF)==pnd.DataFrame:
                    thisSimSto = simStoDF.loc[idx]
                    thisSimElev = elevFromStor(thisSimSto)
                else:
                    thisSimSto= np.nan
                    
                #print("Obs storage: %0.2f TAF, Obs elevation: %0.2f" %(thisObsSto, thisObsElev))
       #->> this was used with matrix/defined levels version of profile data -->>  #ax[nr, nc].plot([x*1.8 + 32 for x in reversed(d.values)], list(reversed(obsElevs)), color='k', label='Observed')
                ax[nr, nc].plot(thisObsProf.loc[:,obsTempCol], thisObsProf.index, 
                                color='k', label='Observed')
                    
                ax[nr, nc].plot(thisSimProf.values,[float(x) for x in thisSimProf.index], 
                                color='orange', label='ResTemp model')
                ax[nr, nc].axhline(thisObsElev, lw=0.5, c='k', ls=':')
                #ax[nr, nc].set_xlim((40, 80))S
                ax[nr, nc].plot([43], [thisObsElev+5], marker='v', c='k',  #<--just puts a little 'hydrat'-like marker at water surface
                                linestyle='None',mew=0.2, ms=4, 
                                fillstyle='none')
                
                # put the open gates/river outlets on figure
                for gi,g in gateTS2.loc[:,[idx]].iterrows(): # i = gate number, g = gate value
                    if g.iloc[0]>0: # there's a gate open on this level
                        gb = resObj.Outlets[gi].BotElevFt
                        gt = resObj.Outlets[gi].TopElevFt
                        ax[nr,nc].axhspan(gt, gb,0,1,color='#a87dbd', alpha=0.2,
                                          zorder=0)
                                          # np.max(thisObsProf.loc[:,obsTempCol]),
                                          # color='#a87dbd', alpha=0.7)
                for ri,r in rivoTS2.loc[:,[idx]].iterrows():
                    if r.iloc[0]>0: # there's a river outlet open
                        rctr = resObj.RiverOutlets[ri].CtrElev
                        ax[nr,nc].axhline(rctr,0,1,c='0.6', lw=2.5, ls='-',zorder=0)
                                          #40, 
                                          #np.max(thisObsProf.loc[:,obsTempCol]),
                                          #c='k', lw=1.5, ls='-')
                        
                
                ax[nr, nc].set_title(idx.strftime('%b %d, %Y'))
                if np.max(thisSimProf.values)>80:
                    ax[nr,nc].set_xticks([40, 50, 60, 70, 80, 90])
                else:
                    ax[nr, nc].set_xticks([40, 50, 60, 70, 80])
                ax[nr, nc].xaxis.set_minor_locator(MultipleLocator(5))
                ax[nr,nc].set_ylim((650, 1070))
                ax[nr, nc].set_yticks([650, 700, 750, 800, 850,  900, 950, 1000, 1050, 1067])
                
                secax = ax[nr, nc].secondary_yaxis('right',functions=(storFromElev,elevFromStor))
                secax.yaxis.set_minor_locator(AutoMinorLocator(0.5))
                
                if (i == len(obsProfs)-1) or (nc == ncols-1):
                    secax.set_ylabel('Storage(MAF)', fontsize=10.)
                else:
                    secax.set_yticklabels([])
                    
                if nc ==0: #== 0: # set equal to zero to jsut label left-hand most column
                    ax[nr,nc].set_ylabel('Elevation (ft)', fontweight='light', fontsize=10.)
                
                numHanging = nrows*ncols - len(obsProfs)
                hangingIndex = ncols - numHanging
                if (nr ==nrows-1) or ((nr == nrows-2) and (nc >= hangingIndex)):
                    ax[nr, nc].set_xlabel('Temperature (deg F)', fontsize=10.)
                i +=1
                
    #        for ix in range(nrows): 
    #            ax[nrows-1,ix].set_xlabel('Temperature, deg F')
    #            ax[ix, 0].set_ylabel('Elevation, ft')
    #        ax[0, nrows-1].set_ylabel('Elevation, ft')
            
            for idx in range(i, nrows*ncols):
                nr = idx//ncols
                nc = idx%ncols
                fig.delaxes(ax[nr, nc])
            sns.despine()
            
            plt.suptitle('Temperature Profile Comparison for Year %d' %(y), fontsize=18)
            
            plt.tight_layout()
            
            if nrows==2:
                if show_stats:
                    top=0.76
                else:
                    top=0.87
                plt.subplots_adjust(top=top,bottom=0.05, hspace=0.41, wspace=0.6)
                txtheight = 0.835
            else:
                if show_stats:
                    top=0.83
                else:
                    top=0.87
                plt.subplots_adjust(top=top, bottom=0.05, hspace=0.41, wspace=0.6)
                txtheight = 0.875
            
            legelems = [Line2D([0],[0], color='orange', lw=2, label='ResTemp Model'),
                        Line2D([0],[0], color='k', lw=2, label='Observation'),
                        Line2D([0],[0], color='k', lw=1, ls=':', label='Obs Elevation')]
            fig.legend(handles=legelems, bbox_to_anchor=(0.65, 0.95),ncol=3, frameon=False, fontsize=12)
            
            if show_stats:
                plt.figtext(0.4, txtheight, "NSE: %0.2f    RMSE: %0.2f deg F (%0.2f deg C)" %(nse, rmse, rmse/1.8),
                            fontsize=11, fontfamily='monospace')
            
            if viewSave=='save':
                plt.savefig(os.path.join(figoutDir, '%s_Profiles_%s.png' %(runName, y)), dpi=300)
                plt.close()
                #return([None, None])
            else:
                return([fig, ax])

def get_sim_obs_data_this_run(resmod, **kwargs):
    '''
    

    Parameters
    ----------
    resmod : TYPE
        DESCRIPTION.
        
    

    Returns
    -------
    resmod: Object
        Added dataframes to the 'calc_errors_df' dictionary - for use in 
        calculating error metrics or plotting

    '''
    if 'dates' in kwargs:
        sdat, edat = kwargs['dates']
    else:
        sdat = resmod.SimulationSpecs.StartDate
        edat = resmod.SimulationSpecs.EndDate
    if 'on_wy' in kwargs:
        on_wy = kwargs['on_wy']
    else:
        on_wy = False
       
    simProfilesDF = resmod.Simulation_Results['ProfilesDF']
    simRelease = resmod.Simulation_Results['ReleaseDF']
    simStorage = pnd.DataFrame(index=resmod.SimDates,data=resmod.Simulation_Results['Storage'])
    simEvap = resmod.Simulation_Results['EvapDF']
    
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simprof_wys = simProfilesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simprof_wys)
    else:
        yrs = np.unique(simProfilesDF.index.map(lambda x: x.year))

    
    modElevs = [l.CtrElev for l in resmod.Layers.values()]
    
    # get observations
    obsProfiles = resmod.Observations.ProfilesDF
    obsHdr = resmod.Observations.ProfilesColMap
    obsTempCol = obsHdr['temperature']
    obsElevs = list(obsProfiles.columns)    
    
    outflow_colmap = resmod.Observations.OutflowColMap
    #obsReleaseTemp = resmod.Observations.OutflowDF.loc[:,outflow_colmap['obs_outtemp_final']]
    relcols = [outflow_colmap['obs_flow_final'][0],outflow_colmap['obs_outtemp_final'][0]]
    obsRelease = resmod.Observations.OutflowDF.loc[:,relcols]
    obsStorage = resmod.Observations.OutflowDF.loc[:,[outflow_colmap['obs_storage_final']]]
    
    resmod.calc_errors_df = {}
    for y in yrs:
        if on_wy:
            if obsProfiles.index.nlevels>1:
                obs_wy = obsProfiles.index.get_level_values(0).map(lambda x: x.year+1 if x.month>9 else x.year)
            else:
                obs_wy = obsProfiles.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        
        # get dates of measured profiles for this year
        if on_wy:
            obsProfs = obsProfiles[obs_wy==y]
            obsRels = obsRelease[obs_wy==y]
        else:
            obsProfs = obsProfiles.loc[str(y),:]
            obsRels = obsRelease.loc[str(y),:]
        
        obsProfDates = obsProfs.index.get_level_values(0).unique()
        simProfsDates = simProfilesDF.index.unique()
        
        obsDates = obsRels.index.unique()
        simDates = simRelease.index.unique()
        
        commonDates = sorted(list(set(obsProfDates) & set(simProfsDates)))
        
        commonDatesTimeSeries = sorted(list(set(obsDates)& set(simDates)))
        
        simProfs = simProfilesDF.loc[commonDates,:]
        simRel = simRelease.loc[commonDatesTimeSeries,:]
        simStor = simStorage.loc[commonDatesTimeSeries,:]
        simEv = simEvap.loc[commonDatesTimeSeries,:]
        
        obsRel = obsRelease.loc[commonDatesTimeSeries, :]
        obsSto = obsStorage.loc[commonDatesTimeSeries,:]
        
        obsRel = obsRelease.loc[commonDatesTimeSeries,:]
        resmod.calc_errors_df[y] = {'prof': {'sim':simProfs,'obs':obsProfs},
                                    'rel': {'sim': simRel, 'obs': obsRel},
                                    'sto': {'sim': simStor, 'obs':obsSto}}
        

def calcErrorMetrics_profiles(resmod,**kwargs):
    

    error_dict = {}
    yrs = [y for y in resmod.calc_errors_df.keys()]
    
    if 'year' in kwargs:
        yrs = [y for y in yrs if y in kwargs['year']]
    
    for y in yrs:
        
        errors = []
        dateList = []
        obsList = []
        simList = []
        month_list = []
        
        sim_profiles_df = resmod.calc_errors_df[y]['prof']['sim']
        obs_profiles_df = resmod.calc_errors_df[y]['prof']['obs']
        monerr = {m:[] for m in range(1,13)}  # blank dictionary to hold monthly error calcs
        monobs = {m:[] for m in range(1,13)} # blank dictionary for the observations for purposes of calculating NSE
        monsim = {m:[] for m in range(1,13)} # blank dictionary for the simulation data
        monerr['annual'] = []
        monobs['annual'] = []
        monsim['annual'] = []
        
        for didx, vals in sim_profiles_df.iterrows():
            print(f"calculating for date: {didx}")
            this_month = didx.month
            if didx in obs_profiles_df.index:
                
                this_sim_prof = sim_profiles_df.loc[didx]
                if obs_profiles_df.index.nlevels==1:  #<-- then only just the date in the index
                    vals = obs_profiles_df.loc[didx].copy()
                else:
                    vals = obs_profiles_df.loc[didx, :].copy()
                vals = vals.sort_index()
                #print(thisSimProf.shape)
                if 'degF' in vals.columns[0]:
                    valsF = vals
                else:
                    valsF = vals.apply(lambda x: x*1.8+32.)
                
                # calculate avg obs temp over the model layer interval, if
                # the obs profile is at higher resolution than the sim
                avg_vals = {}
                if len(valsF) > len(sim_profiles_df.loc[didx]):
                    for l in resmod.Layers:
                        elevrng = [resmod.Layers[l].MinElev, resmod.Layers[l].MaxElev]
                        ctrelev  = resmod.Layers[l].CtrElev
                        try:
                            vi = valsF.loc[elevrng[0]:elevrng[1]]
                        except:
                            vi = []
                        if len(vi)>0:
                            # get average temp over range
                            avg_vals[ctrelev] = np.nanmean(vi)
                        else:
                            avg_vals[ctrelev] = np.nan

                    avg_vals = {k:v for k,v in avg_vals.items() if ~np.isnan(v)}
                    this_sim_prof = sim_profiles_df.loc[didx]
                    this_err = []
                    this_obs = []
                    this_sim = []
                    for i in this_sim_prof.index:
                        if i in avg_vals:
                            simval = this_sim_prof[i]
                            this_err.append(simval - avg_vals[i])
                            this_obs.append(avg_vals[i])
                            this_sim.append(simval)
                            
                else: # sim has same or more # layers as obs
                    this_sim_prof = sim_profiles_df.loc[didx,:]
                    sim_val = np.interp(valsF.index, this_sim_prof.index, 
                                                this_sim_prof.values)
                    this_err = []
                    this_obs = []
                    this_sim = []
                    for i,v in enumerate(sim_val):
                        this_sim.append(v)
                        this_obs.append(valsF.iloc[i][0])
                        this_err.append(v - valsF.iloc[i][0])
                    
                            
                monerr[this_month].append(this_err)
                monobs[this_month].append(this_obs)
                monsim[this_month].append(this_sim)
                monerr['annual'].append(this_err)
                monobs['annual'].append(this_obs)
                monsim['annual'].append(this_sim)
        
        month_stats = {}
        for m in monerr:
            terr = monerr[m]
            tobs = monobs[m]
            tsim = monsim[m]
            flatlist=[element for sublist in terr for element in sublist]
            flatobs = [element for sublist in tobs for element in sublist]
            flatsim = [element for sublist in tsim for element in sublist]
            sumsqerr = np.nansum([er*er for er in flatlist])
            
            nse = 1 - sumsqerr/(np.nansum((flatobs-np.nanmean(flatobs))**2))
            
            mean_bias = np.nansum(flatlist)/len(flatlist)
            mae = np.nansum([abs(i) for i in flatlist])/len(flatlist)
            
            rmse = np.sqrt(sumsqerr/len(flatobs))
            corrmat = np.corrcoef(flatobs, flatsim)
            corrxy = corrmat[0,1]
            r2 = corrxy**2
            
            month_stats[m] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                             'mae': mae}
            
            
        error_dict[y] = month_stats
        
    return(error_dict)

def getValsAtOpenGateLevels(resmod,**kwargs):
    

    error_dict = {}
    yrs = [y for y in resmod.calc_errors_df.keys()]
    
    if 'year' in kwargs:
        yrs = [y for y in yrs if y in kwargs['year']]
    
    for y in yrs:
        
        errors = []
        dateList = []
        obsList = []
        simList = []
        month_list = []
        
        sim_profiles_df = resmod.calc_errors_df[y]['prof']['sim']
        obs_profiles_df = resmod.calc_errors_df[y]['prof']['obs']
        monerr = {m:[] for m in range(1,13)}  # blank dictionary to hold monthly error calcs
        monobs = {m:[] for m in range(1,13)} # blank dictionary for the observations for purposes of calculating NSE
        monsim = {m:[] for m in range(1,13)} # blank dictionary for the simulation data
        for didx, vals in sim_profiles_df.iterrows():
            print(f"calculating for date: {didx}")
            this_month = didx.month
            
            if didx in obs_profiles_df.index:
                
                this_sim_prof = sim_profiles_df.loc[didx]
                
                if obs_profiles_df.index.nlevels==1:  #<-- then only just the date in the index
                    vals = obs_profiles_df.loc[didx].copy()
                else:
                    vals = obs_profiles_df.loc[didx, :].copy()
                vals = vals.sort_index()
                #print(thisSimProf.shape)
                if 'degF' in vals.columns[0]:
                    valsF = vals
                else:
                    valsF = vals.apply(lambda x: x*1.8+32.)
                
                # get the layers that have an open outlet
                open_outlets = [i for i,v in resmod.Operations.GateOps[didx].items() if v>0]
                open_layers = [list(l) for l in [np.arange(resmod.Outlets[o].MinLayer, \
                                                     resmod.Outlets[o].MaxLayer+1) \
                                           for o in open_outlets]]
                open_elevs =  [(resmod.Outlets[o].BotElevFt, resmod.Outlets[o].TopElevFt)  \
                                           for o in open_outlets]
                open_layers = list(set([el for nl in open_layers for el in nl]))
                # calculate avg obs temp over the model layer interval, if
                # the obs profile is at higher resolution than the sim
                avg_vals = {}
                if len(valsF) > len(sim_profiles_df.loc[didx]):
                    for l in open_layers:
                        elevrng = [resmod.Layers[l].MinElev, resmod.Layers[l].MaxElev]
                        ctrelev  = resmod.Layers[l].CtrElev
                        try:
                            vi = valsF.loc[elevrng[0]:elevrng[1]]
                        except:
                            vi = []
                        if len(vi)>0:
                            # get average temp over range
                            avg_vals[ctrelev] = np.nanmean(vi)
                        else:
                            avg_vals[ctrelev] = np.nan

                    avg_vals = {k:v for k,v in avg_vals.items() if ~np.isnan(v)}
                    this_sim_prof = sim_profiles_df.loc[didx]
                    this_err = []
                    this_obs = []
                    this_sim = []
                    for i in this_sim_prof.index:
                        if i in avg_vals:
                            simval = this_sim_prof[i]
                            this_err.append(simval - avg_vals[i])
                            this_obs.append(avg_vals[i])
                            this_sim.append(simval)
                            
                else: # sim has same or more # layers as obs
                    this_sim_prof = sim_profiles_df.loc[didx,:]
                    this_sim_open = this_sim_prof.iloc[open_layers]
                    this_obs_open = {}
                    for l in open_layers:
                        obs_at_sim_open_lyr = np.interp(resmod.Layers[l].CtrElev, valsF.index, 
                                            valsF[valsF.columns[0]])
                        this_obs_open[resmod.Layers[l].CtrElev] = obs_at_sim_open_lyr

                    # sim_val = np.interp(valsF.index, this_sim_open.index, 
                    #                             this_sim_open.values)
                    this_err = []
                    this_obs = []
                    this_sim = []
                    for i,v in this_sim_open.items():
                        this_sim.append(v)
                        this_obs.append(this_obs_open[i])
                        this_err.append(v - this_obs_open[i])
                            
                monerr[this_month].append(this_err)
                monobs[this_month].append(this_obs)
                monsim[this_month].append(this_sim)
            
        
        month_stats = {}
        for m in monerr:
            terr = monerr[m]
            tobs = monobs[m]
            tsim = monsim[m]
            flatlist=[element for sublist in terr for element in sublist]
            flatobs = [element for sublist in tobs for element in sublist]
            flatsim = [element for sublist in tsim for element in sublist]
            sumsqerr = np.nansum([er*er for er in flatlist])
            
            nse = 1 - sumsqerr/(np.nansum((flatobs-np.nanmean(flatobs))**2))
            
            mean_bias = np.nansum(flatlist)/len(flatlist)
            mae = np.nansum([abs(i) for i in flatlist])/len(flatlist)
            
            rmse = np.sqrt(sumsqerr/len(flatobs))
            corrmat = np.corrcoef(flatobs, flatsim)
            corrxy = corrmat[0,1]
            r2 = corrxy**2
            
            month_stats[m] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                             'mae': mae, 'sim_data': tsim, 'obs_data': tobs, 'sim_minus_obs': terr}
            
            
        error_dict[y] = month_stats
        
    return(error_dict)


def calcErrorMetrics_timeseries(resmod,ivar='rel', monthly=False):
    

    error_dict = {}
    yrs = [y for y in resmod.calc_errors_df.keys()]
    
    for y in yrs:
        

        
        sim_data = resmod.calc_errors_df[y][ivar]['sim']
        obs_data = resmod.calc_errors_df[y][ivar]['obs']
        
        if monthly:
            
            monerr = {m:[] for m in range(1,13)}  # blank dictionary to hold monthly error calcs
            monobs = {m:[] for m in range(1,13)} # blank dictionary for the observations for purposes of calculating NSE
            monsim = {m:[] for m in range(1,13)} # blank dictionary for the simulation data
            month_stats = {}
            for m in range(1,13):
                thismidx = [ix for ix in sim_data.index if ix.month==m]
                if ivar=='rel':
                    # calculate for both outflow volume and temperature
                    this_sim_vol = sim_data.loc[thismidx,sim_data.columns[0]]
                    this_obs_vol = obs_data.loc[thismidx,obs_data.columns[0]]
                    
                    this_sim_temp = sim_data.loc[thismidx, sim_data.columns[1]]
                    this_obs_temp = obs_data.loc[thismidx, obs_data.columns[1]]
                    
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
            
            # calculate on an annual basis as well
            this_sim_vol = sim_data.loc[f'{y}',sim_data.columns[0]]
            this_obs_vol = obs_data.loc[f'{y}',obs_data.columns[0]]
            
            this_sim_temp = sim_data.loc[f'{y}', sim_data.columns[1]]
            this_obs_temp = obs_data.loc[f'{y}', obs_data.columns[1]]
    
            err_temp = [s-o for s,o in zip(this_sim_temp, this_obs_temp)]
            sumsqerr = np.nansum([er*er for er in err_temp])
            
            nse = 1 - sumsqerr/(np.nansum((this_obs_temp-np.nanmean(this_obs_temp))**2))
            
            mean_bias = np.nansum(err_temp)/len(err_temp)
            mae = np.nansum([abs(i) for i in err_temp])/len(err_temp)
            
            rmse = np.sqrt(sumsqerr/len(this_obs_temp))
            corrmat = np.corrcoef(this_obs_temp,pnd.to_numeric(this_sim_temp))
            corrxy = corrmat[0,1]
            r2 = corrxy**2
            month_stats['annual'] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                             'mae': mae}
            error_dict[y] = month_stats
        else:
            if ivar=='rel':
                # calculate for both outflow volume and temperature
                this_sim_vol = sim_data.loc[f'{y}',sim_data.columns[0]]
                this_obs_vol = obs_data.loc[f'{y}',obs_data.columns[0]]
                
                this_sim_temp = sim_data.loc[f'{y}', sim_data.columns[1]]
                this_obs_temp = obs_data.loc[f'{y}', obs_data.columns[1]]
        
                err_temp = [s-o for s,o in zip(this_sim_temp, this_obs_temp)]
                sumsqerr = np.sum([er*er for er in err_temp])
                
                nse = 1 - sumsqerr/(np.nansum((this_obs_temp-np.nanmean(this_obs_temp))**2))
                
                mean_bias = np.sum(err_temp)/len(err_temp)
                mae = np.sum([abs(i) for i in err_temp])/len(err_temp)
                
                rmse = np.sqrt(sumsqerr/len(this_obs_temp))
                corrmat = np.corrcoef(this_obs_temp,pnd.to_numeric(this_sim_temp))
                corrxy = corrmat[0,1]
                r2 = corrxy**2
            error_dict[y] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
                             'mae': mae}
        
    return(error_dict)

    #             #thisSimProf = rt.interpProfile(vals, sim_profiles_df.loc[didx])

    #             obsList.append([i[0] if type(i)==np.ndarray else i for i in valsF.values])
    #             #print(thisSimProf.values)
    #             simList.append([i[0] for i in thisSimProf.values]) #[v for i,v in thisSimProf.iterrows()])
    #             # flatten the df's to make sure there is just one column and one index
    #             if len(thisSimProf.shape)>1:
    #                 thisSimProf = thisSimProf.iloc[:,0]
    #             if len(valsF.shape)>1:
    #                 valsF = valsF.iloc[:,0]
                
    #             errors.append([i for i in (valsF.values - thisSimProf.values)])
    #             month_list.append(this_month)
    #             dateList.append(didx)
    #             monerr[this_month].append(errors[-1])
        
    #     errors = [er for er2 in errors for er in er2 if ~np.isnan(er)]
    #     month_stats = {}
    #     for m in list(set(month_list)):
    #         idcs = [i for i in range(len(month_list)) if month_list[i] == m]
    #         if len(idcs)>0:
    #             sel_errors = []
    #             for i in idcs:
    #                 sel_errors.append([errors[i] for i in idcs])
                
            
        


    #     sumsqerr = np.sum([er*er for er in errors])
            
    #     #obs = [o1 for o2 in obs for o1 in o2 if ~np.isnan(o1)]
    #     simList = [i for isub in simList for i in isub]
    #     obsList = [i for isub in obsList for i in isub]
    #     #print(obsList)    
    #     nse = 1 - sumsqerr/(np.nansum((obsList-np.nanmean(obsList))**2))
        
    #     mean_bias = np.sum(errors)/len(errors)
    #     mae = np.sum([abs(i) for i in errors])/len(errors)
        
    #     rmse = np.sqrt(sumsqerr/len(obsList))
    #     corrmat = np.corrcoef(obsList, simList)
    #     corrxy = corrmat[0,1]
    #     r2 = corrxy**2
        
    #     error_dict[y] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
    #                      'mae': mae}
    # return(error_dict)
            
def calcErrorMetrics(obsProfilesDF, simProfilesDF, weights=1):
    
    
    errors = []
    dateList = []
    obsList = []
    simList = []
    for didx, vals in simProfilesDF.iterrows():

        if didx in obsProfilesDF.index:
            if obsProfilesDF.index.nlevels==1:  #<-- then only just the date in the index
                vals = obsProfilesDF.loc[didx].copy()
            else:
                vals = obsProfilesDF.loc[didx, :].copy()
            thisSimProf = rt.interpProfile(vals, simProfilesDF.loc[didx])
            #print(thisSimProf.shape)
            if 'degF' in vals.columns[0]:
                valsF = vals
            else:
                valsF = vals.apply(lambda x: x*1.8+32.)
            obsList.append([i[0] if type(i)==np.ndarray else i for i in valsF.values])
            #print(thisSimProf.values)
            simList.append([i[0] for i in thisSimProf.values]) #[v for i,v in thisSimProf.iterrows()])
            # flatten the df's to make sure there is just one column and one index
            if len(thisSimProf.shape)>1:
                thisSimProf = thisSimProf.iloc[:,0]
            if len(valsF.shape)>1:
                valsF = valsF.iloc[:,0]
            
            errors.append([weights*i for i in (valsF.values - thisSimProf.values)])
            dateList.append(didx)
    
    errors = [er for er2 in errors for er in er2 if ~np.isnan(er)]

    sumsqerr = np.sum([er*er for er in errors])
        
    #obs = [o1 for o2 in obs for o1 in o2 if ~np.isnan(o1)]
    simList = [i for isub in simList for i in isub]
    obsList = [i for isub in obsList for i in isub]
    #print(obsList)    
    nse = 1 - sumsqerr/(np.nansum((obsList-np.nanmean(obsList))**2))
    
    rmse = np.sqrt(sumsqerr/len(obsList))
    corrmat = np.corrcoef(obsList, simList)
    corrxy = corrmat[0,1]
    r2 = corrxy**2
    
    #print(nse, rmse, r2)
    return([errors, dateList, nse, rmse, r2])

def profile_error_to_csv(error_dict, resmod,draft=True, **kwargs):
    '''
    Writes out the error metrics calculated for simulated and observed temperature
    profiles to a csv. Output location will be the 'outputs/stats' in the model 
    directory
    

    Parameters
    ----------
    error_dict : TYPE
        Error dictionary created from a call to calcErrorMetrics_profiles
        
    resmod : OBJECT
        ResTemp model object
        
    draft : Boolean, optional
        Put 'DRAFT' in output file name. The default is True.

    Returns
    -------
    None.

    '''

    run_name = resmod.RunName 
    PROJ_DIR = resmod.ProjDir
    if 'proj_dir' in kwargs:
        outDir = os.path.join(kwargs['proj_dir'], 'outputs')
    else:
        outDir = os.path.join(PROJ_DIR, 'outputs')
    
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    
    if 'proj_dir' in kwargs:
        stats_dir = os.path.join(kwargs['proj_dir'], 'outputs', 'stats')
    else:
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
                      mdat['mae'],mdat['r2'],mdat['rmse']/1.8,
                      mdat['mean_bias']/1.8,mdat['mae']/1.8 ]
            alldat.append(valdat)
            
    alldf = pnd.DataFrame(alldat, columns=['Year','Month','NSE_degF','RMSE_degF',
                                           'MeanBias_degF','MAE_degF','R2','RMSE_degC',
                                           'MeanBias_degC','MAE_degC'])

    if draft:
        draftstr = 'DRAFT'
    else:
        draftstr = ''
        
    outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_{run_name}_profile_error_stats.csv')
    alldf.to_csv(outfp, header=True)


def timeseries_error_to_csv(modelobj, error_dict, model_name='',
                            draft=True, **kwargs):
    
    run_name = modelobj.RunName
    PROJ_DIR = modelobj.ProjDir
    
    if 'run_name' in kwargs:
        run_name = kwargs['run_name']
       
    
    if 'proj_dir' in kwargs:
        outDir = os.path.join(kwargs['proj_dir'], 'outputs')
    else:
        outDir = os.path.join(PROJ_DIR, 'outputs')
    
    if not os.path.exists(outDir):
        os.mkdir(outDir)
    
    if 'proj_dir' in kwargs:
        stats_dir = os.path.join(kwargs['proj_dir'], 'outputs', 'stats')
    else:
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
        
    outfp = os.path.join(stats_dir, f'{todaystr}{draftstr}_{run_name}_{model_name}_releasetemp_errorstats.csv')
    alldf.to_csv(outfp, header=True)
    
def plot_evap(resObj,viewSave='save', on_wy=False):
    
    runName = resObj.RunName
    
    if not resObj.SimulationSpecs.CalcEvaporation:
        print("\n****************************************************")
        print("  Evaporation was provided as an input - nothing to plot")
        print("****************************************************\n")
        
    else:
        if 'evap' not in resObj.Met.ColumnMap:
            print(" --- Don't have evaporation observations to compare to,\n --- try adding a column for evaporation to the input time series")
        else:
            evapCol = resObj.Met.ColumnMap['evap']
        
        evapObs = resObj.Met.DataFrame.loc[:,evapCol]
        #evapSim = pnd.DataFrame([i[1] for i in resObj.Evaporation], index=resObj.TimeStepList)
        evapSim = resObj.Simulation_Results['EvapDF']
        
        PROJ_DIR = resObj.ProjDir
        
        figoutDir = os.path.join(PROJ_DIR, 'outputs', 'figs')
        if not os.path.exists(figoutDir):
            os.mkdir(figoutDir)
            
        dts = copy.deepcopy(resObj.TimeStepList)
        if on_wy:
            #simProfilesDF = np.copy(simProfilesDF, deep=True)
            sim_wys = [x.year+1 if x.month>9 else x.year for x in dts] #simProfilesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
            yrs = np.unique(sim_wys)
        else:
            yrs = np.unique([x.year for x in dts])
            
        for y in yrs:
            if on_wy:
                dtlist = [d for d in sim_wys if d==y]
            else:
                dtlist = [d for d in dts if d.year==y]
                
            obsEvap = evapObs.loc[dtlist]*86400/43560. # obs evap in cfs - convert to AF
            simEvap = evapSim.loc[dtlist,:]
            
            with sns.plotting_context('paper', font_scale=1.2):
                fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(10, 6), sharex=True)        
    
                ax.plot(obsEvap.index, obsEvap, color='k', label='CDEC Measured Evap')
                ax.plot(simEvap.index, simEvap.loc[:,simEvap.columns[0]], color='orange', label='Simulated Evap')
                ax.set_xlabel('Date')
                ax.set_ylabel('Evaporation Volume, acre-feet')
                ax.xaxis.set_major_locator(months)
                ax.xaxis.set_major_formatter(monyr_fmt)


                if on_wy:
                    plt.suptitle('Evaporation Comparison, WY %s' %y)
                else:
                    plt.suptitle('Evaporation Comparions, %s' %y)
                sns.despine()
                plt.tight_layout(h_pad=0.1)
                
            if viewSave=='save':
                plt.savefig(os.path.join(figoutDir, f'{runName}_{y}_Evap.png'), dpi=300)
                plt.close()
                #return([None, None])
            else:
                return([fig, ax]) 

def side_gate_dates(resObj, days_consec_open=5, doy_start_count=181,**kwargs):
    
    release_inputs = resObj.Outflow.DataFrame
    release_inputs_colmap = resObj.Outflow.ColumnMap
    obscm = resObj.Observations.OutflowColMap
    gateNames = [n for n in release_inputs_colmap['gates'].values()]
    yrs = list(set([i.year for i in resObj.SimDates]))
    
    if 'year' in kwargs:
        years = kwargs['year']
    else:
        years = yrs
    
    side_gate_ops = {}
    for y in years:
        thisDateList = release_inputs.loc[str(y),:].index
        gateTS = release_inputs.loc[thisDateList, gateNames] #['Upp','Mid','PRG','SDG']
        # gateTS2 = gateTS.T
        # gateTS2.index = [4, 3, 2, 1]
        # gateTS2 = gateTS2.sort_index()
        obs_gate_ops = gateTS
        
        sim_gate_dict = resObj.Operations.GateOps
        sim_gateTS = pnd.DataFrame.from_dict(sim_gate_dict, orient="index")
        sim_gateTS.columns = gateNames
        # sim_gateTS = sim_gateTS.T.loc[thisDateList,:].T
        # sim_gateTS.index = [4,3,2,1]
        # sim_gate_ops= sim_gateTS.sort_index()
        

        obs_provisional_first_date = thisDateList[doy_start_count]
        sim_provisional_first_date = thisDateList[doy_start_count]
        
        obs_continuous_sdg_days = 0
        sim_continuous_sdg_days = 0
        
        obs_both_sdg_days = 0
        sim_both_sdg_days = 0
        obs_first_date = None
        sim_first_date = None
        obs_both_date = None
        
        for d in thisDateList:
            # are we far enough in to the year to start checking for side gates?
            if d.timetuple().tm_yday >=doy_start_count:
                this_obs_sdg = obs_gate_ops.loc[d,gateNames[3]]
                this_sim_sdg = sim_gateTS.loc[d, gateNames[3]]
                
                if this_obs_sdg > 0:
                    obs_continuous_sdg_days +=1
                    if obs_continuous_sdg_days==1:
                        obs_provisional_first_date=d
                    if obs_continuous_sdg_days>=days_consec_open:
                        obs_first_date = obs_provisional_first_date
                else:
                    obs_continuous_sdg_days = 0
                    
                if this_sim_sdg > 0:
                    sim_continuous_sdg_days +=1
                    if sim_continuous_sdg_days==1:
                        sim_provisional_first_date=d
                    if sim_continuous_sdg_days>=days_consec_open:
                        sim_first_date = sim_provisional_first_date
                else:
                    sim_continous_sdg_days = 0
                    
                if this_obs_sdg==2:
                    obs_both_sdg_days+=1
                    if obs_both_sdg_days==1:
                        obs_prov_both_date = d
                        
                    if obs_both_sdg_days>=days_consec_open:
                        obs_both_date = obs_prov_both_date
                        
                else:
                    obs_both_sdg_days=0
                    
                if this_sim_sdg==2:
                    sim_both_sdg_days+=1
                    if sim_both_sdg_days==1:
                        sim_prov_both_date = d
                    
                    if sim_both_sdg_days>=days_consec_open:
                        sim_both_date = sim_prov_both_date
                        
                else:
                    sim_both_sdg_days = 0
                    
                if sim_both_sdg_days ==0: # both side gates never opened in simulation
                    sim_both_date = None
                    
        side_gate_ops[y] = {'first_sdg': {'obs': obs_first_date, 'sim': sim_first_date},
                            'full_sdg' : {'obs': obs_both_date, 'sim': sim_both_date}}
        
    return(side_gate_ops)
                    
                    
                
        
    
def analyze_gate_ops(runName, resObj, simReleasesDF):
    PROJ_DIR = resObj.ProjDir
    figoutDir = os.path.join(PROJ_DIR, 'outputs', 'figs')
    if not os.path.exists(figoutDir):
        os.mkdir(figoutDir)
        
    if on_wy:
        #simProfilesDF = np.copy(simProfilesDF, deep=True)
        simrels_wys = simReleasesDF.index.map(lambda x: x.year+1 if x.month>9 else x.year)
        yrs = np.unique(simrels_wys)
    else:
        yrs = np.unique(simReleasesDF.index.map(lambda x: x.year))


    #viewSave ='save'
    
    if viewSave=='save':
        plt.ioff()
    else:
        plt.ion()
    
    obsRelease = resObj.Observations.OutflowDF
    obsHdr = resObj.Observations.OutflowColMap 
    
    # input time series
    inTS = resObj.InputTimeSeries.DataFrame
    inTScm = resObj.InputTimeSeries.ColumnMap
    
    gateNames = [n for n in inTScm['gates'].values()]    
    
    for y in yrs:
        
        if on_wy:
            thisDateList = simReleasesDF[simrels_wys==y].index

        else:
            thisDateList = simReleasesDF[str(y)].index
        
        
        print("Making plot for year %d" %y)
        #thisDateList = simReleasesDF[str(y)].index
        
        obsQ = obsRelease.loc[thisDateList, obsHdr['flow']]
        obsT = obsRelease.loc[thisDateList, obsHdr['temperature']]
        

        # make a time series of gate changes
        gateTS = inTS.loc[thisDateList, gateNames] #['Upp','Mid','PRG','SDG']
        gateTS2 = gateTS.T
        gateTS2.index = [4, 3, 2, 1]
        gateTS2 = gateTS2.sort_index()
        
        # get simulated gate changes if predicted_gates is True and the GateOps data exists
        # in the reservoir object
        if (len(resObj.Operations.GateOps)>0):
            sim_gate_dict = resObj.Operations.GateOps
            sim_gateTS = pnd.DataFrame.from_dict(sim_gate_dict)
            sim_gateTS = sim_gateTS.T.loc[thisDateList,:].T
            sim_gateTS.index = [4,3,2,1]
            sim_gateTS = sim_gateTS.sort_index()
        else:
            print("Simulated gate operations not found!")
            print("Try loading from a *.pkl file")
            break
        

        
        for g in sim_gateTS.index:
            ngates = resObj.Outlets[g-1].NumGates
            confmat = np.zeros((ngates, ngates), dtype=int) 
            this_gate_sim = sim_gateTS.loc[g,:]
            this_gate_obs = gateTS2.loc[g,:]
            
            gateopts = list(np.arange(0, ngates+1))
            cats = ['{0[0]}{0[1]}'.format(tup) for tup in product(gateopts, gateopts)]
            cats_levs = ['{0[0]}{0[1]}'.format(tup) for tup in product([0,1], [0,1])]
            t1 = [str(i) for i in this_gate_obs.values.astype(int)]
            t2 =  [str(i) for i in this_gate_sim.values.astype(int)]
            t1lev = [str(min(1,i)) for i in this_gate_obs.values.astype(int)]
            t2lev = [str(min(1,i)) for i in this_gate_sim.values.astype(int)]
            t3 = [i+j for i,j in zip(t1, t2)]
            t3lev = [k+l for k,l in zip(t1lev, t2lev)]
            tmp =pnd.Categorical((this_gate_obs.values.astype(int).astype(str) + this_gate_sim.values.astype(int).astype(str)),
                            categories=cats).value_counts()
            tmp = pnd.Categorical(t3, categories=cats).value_counts()
            tmp_levs = pnd.Categorical(t3lev, categories=cats_levs).value_counts()
            
#            confmat 
#            acc_gates = 
#            for d in thisDateList:
#                simg = sim_gateTS.loc[g,d]
#                obsg = gateTS2.loc[g,d]
#                if simg==obsg
            
        
    
        
if __name__ == "__main__":
    
    
    sys.path.append(r'D:\02_Projects\SacTemp\SimTemp\ResTempPkg')
    import restemp as rt
    import restemp_util as rt_util
    
    PROJ_DIR = os.getcwd()
    print("making plots - called from location: %s" %PROJ_DIR)
        
    import run_model_v7_matlab as rm
    
    args = sys.argv
    
    #inp1 = sys.argv[1]

    #print(inp1)
       
    if len(args)==1:
        print("Insufficient arguments provided - specifiy a plot type to make:")
        print(" Options: \n\t\t\t- (S)torage\n\t\t\t- (R)elease\n\t\t\t- (P)rofile\n\t\t\t- (A)ll" )
        
    else: 
        objfp = args[1]
        
        if not os.path.exists(objfp):
            objfp = os.path.join(PROJ_DIR, objfp)
        if not os.path.exists(objfp):
            print("ERROR:: couldn't find path specified: %s" %objfp)
            exit
        
        if objfp[-3:]=='pkl':
            with open(objfp,'rb') as infile:
                resObj = pickle.load(infile)
            #runName = os.path.basename(inp1)[0:-4] # assume runname is base filename
        else:
            resObj = rm.get_model(objfp)
        
#        bname = os.path.basename(objfp)
#        runName = bname[0:-4]
#        with open(objfp, 'rb') as obf:
#            resObj = pickle.load(obf)
        runName = resObj.RunName
        
        if resObj.SimulationSpecs.StartDate.month==1:
            on_wy=False
        else:
            on_wy=True
        
        # get the dataframes of simulated results too
        simProfDF = pnd.read_table(os.path.join(PROJ_DIR, 'outputs',runName + '_simProf.csv'),
                                   sep=',', index_col=0, parse_dates=True)
        simRelDF = pnd.read_table(os.path.join(PROJ_DIR, 'outputs',runName + '_simRelease.csv'),
                                  sep=',', index_col=0, parse_dates=True)
        simStorDF = pnd.read_table(os.path.join(PROJ_DIR, 'outputs', runName + '_simStor.csv'),
                                   sep=',', index_col=0, parse_dates=True)
         
        argsUpp = [a.upper() for a in args]
        if 'S' in argsUpp:
            plotStorageCompare(resObj, simStorDF)
         
        if 'R' in argsUpp:
            if 'SIMGATE' in argsUpp:
                predgates = True
            else:
                predgates = False
                
            if 'TEMPTARG' in argsUpp:
                targ_ts = args[argsUpp.index('TEMPTARG')+1] # get the variable for target temperature time series (assuming it's in the input ts) for plotting
            else:
                targ_ts = None
                
            print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
            print(" Plotting releases with the following options:            ")
            print("     Include target temperature:   %s                     " %targ_ts)
            print("     Include simulated gate ops:   %s                     " %predgates)
            plotReleasesCompare(runName, resObj, simRelDF, predicted_gates=predgates,
                                target_ts=targ_ts, on_wy = on_wy)

            
        if 'P' in argsUpp:
            plotProfilesCompare(runName, resObj, simProfDF, simStoDF=simStorDF)
            
        if 'A' in argsUpp:
            plotStorageCompare(runName, resObj, simStorDF)
            
            if 'SIMGATE' in argsUpp:
                predgates = True
            else:
                predgates = False
                
            if 'TEMPTARG' in argsUpp:
                targ_ts = args[argsUpp.index('TEMPTARG')+1] # get the variable for target temperature time series (assuming it's in the input ts) for plotting
            else:
                targ_ts = None
                
            plotReleasesCompare(runName, resObj, simRelDF, predicted_gates=predgates,
                                target_ts=targ_ts, on_wy=on_wy)
            #plotReleasesCompare(runName, resObj, simRelDF)
            
            
            plotProfilesCompare(runName, resObj, simProfDF, on_wy=on_wy)
    

else:
    
       
    print(__name__)
    ## first option here - load res temp model code from local source
    #path_to_src = os.path.join(__loader__.path, 'py_src')
    
    # second option - load res temp model from more general source
    sys.path.insert(0, r'D:\02_Projects\SacTemp\SimTemp\R3TAM\src')
    from r3tam import restemp as rt
    from r3tam import restemp_util as rt_util

    
    #import res_ops as rops
    

    print("Loaded run functions successfully")    