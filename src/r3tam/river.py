# -*- coding: utf-8 -*-
"""
River.py - an empirical regression-based river temperature model

@author: jgilbert
"""
import numpy as np
import os, sys
import yaml
import pandas as pnd
import datetime as dt
import time 

from collections import OrderedDict as Odict
import pymatreader as pymat


def c_to_f(c):
    return(c*1.8+32)
def f_to_c(f):
    return((f-32)/1.8)

def c_to_f_diff(c):
    return(c*1.8)
def f_to_c_diff(f):
    return(f/1.8)

class RiverGrid(object):
    
    def __init__(self):
        
        self.Name=''
        
class inputTS:
    
    def __init__(self):
        self.DataFrame = None
        self.ColumnMap = None   
        
class River(object):
    
    def __init__(self):
        
        self.RunName = ''
        self.TimeStep = -1
        self._time = -1
        self.TimeStepDate = dt.datetime(2000, 10, 1, 0, 0)  #default
        self.Nodes = Odict()
        self.Connectivity = Odict() 
        
    @classmethod  
    def initialize_model(cls, config_fp, **kwargs):
        
        print(f"Loading configuration file: {config_fp}")
        if not os.path.exists(config_fp):
            raise FileNotFoundError(config_fp)
            return
        
        if os.path.splitext(config_fp)[1] in ['.yaml','.yml']:
            with open(config_fp) as stream:
                inputs = yaml.load(stream, Loader=yaml.SafeLoader)  
                
        projDir = os.path.dirname(config_fp)
        runName = inputs['RunName']

        print("working directory set as: %s" %projDir)
        
        rivmod = cls()  # instantiate reservoir model object
        rivmod.RunName = runName
        
        rivmod.ConfigFP = config_fp
        
        startDateStr = inputs['Timing']['StartDate']
        endDateStr = inputs['Timing']['EndDate']
        delt_days = inputs['Timing']['TimeStepDays']
        rivmod.SimDates = list(pnd.date_range(startDateStr, endDateStr, freq='d'))
        rivmod.TimeStepDate = rivmod.SimDates[0]
        
        for nn, nv in inputs['Nodes'].items():
            
            newNode = cls.Node()
            newNode.ID = nv[0]
            newNode.UpNode = nv[1]
            newNode.DownNode = nv[2]
            rivmod.Nodes[nn] = newNode
        
        for nn, nv in rivmod.Nodes.items():
            #upnode =nv.UpNode
            for nn2, nv2 in rivmod.Nodes.items():
                if nv2.ID == nv.UpNode:
                    rivmod.Nodes[nn].UpNodeName = nn2
                    
        
        for coefn, coefv in inputs['Coefficients'].items():
            if '_' in coefn:
                coefn,tunits = coefn.split('_')
                units_flag=True
            else:
                tunits='degc'
                units_flag=False
            if coefn in rivmod.Nodes:
                if units_flag:
                    thiscoefv = inputs['Coefficients'][f"{coefn}_{tunits}"]
                else:
                    thiscoefv = inputs['Coefficients'][coefn]
                intrvl = thiscoefv['interval']
                variables = thiscoefv['variables']
                coeflist = []
                coeff_time_keys = list(thiscoefv.keys())[2:]
                for c in coeff_time_keys: #range(1,13):
                    coeflist.append(thiscoefv[c])
                tdf = pnd.DataFrame(data=coeflist, columns=variables, 
                                    index=coeff_time_keys) #range(1,13) )
            rivmod.Nodes[coefn].Coefficients = tdf
            rivmod.Nodes[coefn].CoeffInterval = intrvl
            rivmod.Nodes[coefn].CoeffUnits = tunits
                    
        # get data for boundary and met conditions
        for bcn,bcv in inputs['BoundaryConditions'].items():
            for nn, nv in rivmod.Nodes.items():
                if nv.UpNode==bcn:
                    assign_to_node = nn # nv.ID
            if 'File' in bcv:
                dfp = os.path.join(projDir,bcv['File']['Path'])
                if os.path.exists(dfp):
                    tmpdf = pnd.read_table(dfp,sep=',', index_col=0, parse_dates=True) #, infer_datetime_format=True)
                else:
                    raise FileNotFoundError(f"Couldn't find specified boundary condition file: {dfp}")
                    
                hdrDict = bcv['File']['Header']
                
                # assign
                rivmod.Nodes[assign_to_node].BoundaryConditions.DataFrame = tmpdf
                rivmod.Nodes[assign_to_node].BoundaryConditions.ColumnMap = hdrDict
            
        # get obs data
        for bcn,bcv in inputs['Observations'].items():
            if 'File' in bcv:
                dfp = os.path.join(projDir,bcv['File']['Path'])
                if os.path.exists(dfp):
                    tmpdf = pnd.read_table(dfp,sep=',', index_col=0, parse_dates=True) #, infer_datetime_format=True)
                else:
                    raise FileNotFoundError(f"Couldn't find specified boundary condition file: {dfp}")
                    
                hdrDict = bcv['File']['Header']

                for nn in hdrDict:
                    theseunits = hdrDict[nn][1]
                    col = hdrDict[nn][0]
                    rivmod.Nodes[nn].Observations.DataFrame=tmpdf.loc[:,[col]]
                    rivmod.Nodes[nn].Observations.ColumnMap = \
                                {nn: col, 'units': theseunits}

                        
        return(rivmod)
    
    def set_node_flow_temp(rivmod, selnode, temp, flow, **kwargs):
        
        bcnode = rivmod.Nodes[selnode].BoundaryConditions
        if 'set_date' in kwargs:
            bcnode.DataFrame.loc[kwargs['set_date'], bcnode.ColumnMap['upflow']] = flow
            bcnode.DataFrame.loc[kwargs['set_date'], bcnode.ColumnMap['uptemp']] = temp
        else:
            bcnode.DataFrame.loc[rivmod.TimeStepDate, bcnode.ColumnMap['upflow']] = flow
            bcnode.DataFrame.loc[rivmod.TimeStepDate, bcnode.ColumnMap['uptemp']] = temp
                
        
    
    def advance_rivmod(rivmod, final=False):

        

        if final:
            rivmod._time +=1
            rivmod.TimeStep += 1 #resmod.TimeStep
            rivmod.TimeStepDate = rivmod.SimDates[rivmod.TimeStep]
            rivmod.DayOfYear = rivmod.TimeStepDate.timetuple().tm_yday
            
            for nn, nv in rivmod.Nodes.items():
                if nv.UpNode>0:
                    rivmod.calcRivTemps(nn, final=final)
                   
        else:
            tmp_time = rivmod._time + 1
            tmp_TimeStep = rivmod.TimeStep + 1
            tmp_TimeStepDate = rivmod.SimDates[tmp_TimeStep]
            time_index = tmp_time
            rivmod.DayOfYear = rivmod.TimeStepDate.timetuple().tm_yday
            
            for nn, nv in rivmod.Nodes.items():
                if nv.UpNode>0:
                    rivmod.calcRivTemps(nn, final=final, set_time_step= tmp_TimeStep)
                   
            
        
    def calcRivTemps(rivmod,selnode, final=False, **kwargs):
        
        if 'set_time_step' in kwargs:
            this_step = kwargs['set_time_step']
            this_date = rivmod.SimDates[this_step]
        else:
            this_step = rivmod.TimeStep
            this_date = rivmod.TimeStepDate 
        thismonth = this_date.month #rivmod.TimeStepDate.month
        
        nodemod = rivmod.Nodes[selnode]
        bcnode = rivmod.Nodes[nodemod.UpNodeName].BoundaryConditions
        if type(bcnode.DataFrame) ==type(None): # upstream node is simulated - get value from OutflowTemp
            bc_vec = None
            thisAirT = np.nan
            thisWatT = rivmod.Nodes[nodemod.UpNodeName].OutflowTemp[this_step]
            thisQ = np.nan
        else:
            bc_vec = bcnode.DataFrame.loc[this_date] #rivmod.TimeStepDate]
            thisAirT = bc_vec[bcnode.ColumnMap['airtemp']]
            thisWatT = bc_vec[bcnode.ColumnMap['uptemp']]            
            thisQ = bc_vec[bcnode.ColumnMap['upflow']]
        
        if nodemod.CoeffInterval.lower() == 'monthly':
            thiscoef = nodemod.Coefficients.loc[thismonth,:]
            thiscoefunits = nodemod.CoeffUnits     
            if len(thiscoef)==4:
                # this is a 4-coeff model
                val = thiscoef.airtemp*thisAirT +thiscoef.uptemp*thisWatT +\
                            thiscoef.upflow*thisQ + thiscoef.intcpt
                if final:
                    nodemod.OutflowTemp.append(val)
                else:
                    if len(nodemod.OutflowTemp)<1:
                        nodemod.OutflowTemp.append(val)
                    else:
                        nodemod.OutflowTemp[this_step] = val

            if len(thiscoef)==2:
                # this is a 2-coeff model
                if thiscoefunits=='degf':
                    thisWatT = thisWatT*1.8+32
                val = thiscoef.uptemp*thisWatT + thiscoef.intcpt
                if thiscoefunits=='degf':
                    val = (val-32.)/1.8
                if final:
                    nodemod.OutflowTemp.append(val)
                else:
                    if len(nodemod.OutflowTemp)<1:
                        nodemod.OutflowTemp.append(val)
                    else:
                        nodemod.OutflowTemp[this_step] = val
                    
    def plot_sim_obs(rivmod, selnode, plot_units='degF', **kwargs):
        import seaborn as sns
        sns.set_style('ticks')
        import matplotlib.pyplot as plt
        from pandas.plotting import register_matplotlib_converters
        register_matplotlib_converters()
        
        from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,
                                       AutoMinorLocator)
        
        import matplotlib.dates as mdates
        years = mdates.YearLocator()   # every year
        fiveyears = mdates.YearLocator(5)
        months = mdates.MonthLocator()  # every month
        pltdays = mdates.DayLocator((5,10,15,20, 25))
        years_fmt = mdates.DateFormatter('%Y')
        months_fmt = mdates.DateFormatter('%b')
        monyr_fmt = mdates.DateFormatter('%b\n%Y')
        
        if plot_units=='degF':
            plun = '$^oF$'
        else:
            plun = '$^oC$'
        
        dts = rivmod.SimDates 
        if len(dts)<len(rivmod.Nodes[selnode].OutflowTemp):
            simtemp = pnd.Series(index=dts,data=rivmod.Nodes[selnode].OutflowTemp[1:])
        elif len(dts)> len(rivmod.Nodes[selnode].OutflowTemp):
            simtemp = pnd.Series(index=dts[1:],data=rivmod.Nodes[selnode].OutflowTemp)
        else:
            simtemp = pnd.Series(index=dts,data=rivmod.Nodes[selnode].OutflowTemp)
        
        colname = rivmod.Nodes[selnode].Observations.ColumnMap[selnode]
        obstemp = rivmod.Nodes[selnode].Observations.DataFrame.loc[dts,colname]
        obstempunits = rivmod.Nodes[selnode].Observations.ColumnMap['units']
        
        if 'year' in kwargs:
            sel_year = str(kwargs['year'])
            simtemp = simtemp.loc[sel_year]
            obstemp = obstemp.loc[sel_year]
            dts = list(simtemp.index)
            
        if 'add_target' in kwargs:
            addTarg = kwargs['add_target']
            
        if obstempunits=='degF':
            if plot_units=='degF':
                # convert sim units to deg F
                simtemp_final = pnd.DataFrame(index=dts,data=[_*1.8+32 for _ in simtemp])
                obstemp_final = obstemp
            else: # plot units in deg C, convert obs to deg C
                obstemp_final = (obstemp-32)/1.8 #[(_-32)/1.8 for _ in obstemp.loc[dts]]
                simtemp_final = pnd.DataFrame(index=dts,data=simtemp)
        else: # obstempunits are in deg C
            if plot_units=='degF':
                # convert obs units to deg C
                obstemp_final = (obstemp-32)/1.8 # [(_-32)/1.8 for _ in obstemp.values]
                simtemp_final = pnd.DataFrame(index=dts,data=simtemp)
            else: # plot units are in deg C, nothing needs converting
                obstemp_final = obstemp
                simtemp_final = pnd.DataFrame(index=dts,data=simtemp)
        
        diff = simtemp_final.iloc[:,0]-obstemp_final.iloc[:]
        
        # check if we're adding a temperature target to the time series
        addTarg = False
        if 'temp_targ' in kwargs:
            addTarg = True
            # argument to temp_targ keyword should be a model coupling object
            # that has the temperture target time series associated with it
            ttarg_data = kwargs['temp_targ'].Temperature_Target.Target_TS.DataFrame.copy(deep=True)
            ttarg_cm = kwargs['temp_targ'].Temperature_Target.Target_TS.ColumnMap
            ttarg_units = ttarg_cm['ttarg_specified'].split('_')[-1]
            ttarg_data[ttarg_data>=99] = np.nan
            # get the upper and lower bounds time series
            if 'year' in kwargs:
                upptol = ttarg_data.loc[sel_year, ttarg_cm['ttarg_upp_tol']]
                lowtol = ttarg_data.loc[sel_year, ttarg_cm['ttarg_low_tol']]
                targ_ts = ttarg_data.loc[sel_year, ttarg_cm['ttarg_specified']]
                upptarg = targ_ts+upptol
                lowtarg = targ_ts-lowtol
            else:
                upptol = ttarg_data.loc[:, ttarg_cm['ttarg_upp_tol']]
                lowtol = ttarg_data.loc[:, ttarg_cm['ttarg_low_tol']]
                targ_ts = ttarg_data.loc[:, ttarg_cm['ttarg_specified']]
                upptarg = targ_ts+upptol
                lowtarg = targ_ts-lowtol
        
        if 'add_tol' in kwargs:
            addTol = True
        else:
            addTol = False
        
        if 'plot_off' in kwargs:
            plt.ioff()
        else:
            plt.ion()
        #pnd.DataFrame(index=dts)
        #diff['sim']
        with sns.plotting_context('notebook', font_scale=1.0):
            import matplotlib.gridspec as gridspec
            
            fig2 = plt.figure(figsize=(10,8))
            spec2 = gridspec.GridSpec(ncols=3, nrows=9, figure=fig2)
            
            f2_ax1 = fig2.add_subplot(spec2[0:6,:])
            f2_ax2 = fig2.add_subplot(spec2[7:,:])
            #fig, ax= plt.subplots(1,1, figsize=(10,8))
            f2_ax1.plot(obstemp_final, label='Obs', lw=2, c='k')
            f2_ax1.plot(simtemp_final, label='Sim', lw=2, c='darkred')
            
            if addTarg:
                if plot_units == 'degF':
                    if ttarg_units == 'degC':
                        targ_ts_plot = targ_ts*1.8+32.
                        upptarg_plot = upptarg*1.8+32.
                        lowtarg_plot = lowtarg*1.8+32
                    else:
                        targ_ts_plot = targ_ts
                        upptarg_plot = upptarg
                        lowtarg_plot = lowtarg
                else:
                    # plot units are in deg C
                    if ttarg_units =='degF':
                        targ_ts_plot = (targ_ts-32)/1.8
                        upptarg_plot = (upptarg-32)/1.8
                        lowtarg_plot = (lowtarg-32)/1.8
                    else:
                        targ_ts_plot = targ_ts
                        upptarg_plot = upptarg
                        lowtarg_plot = lowtarg
                        
                f2_ax1.plot(targ_ts_plot, color='k', lw=0.7, ls='--', label='Target\nTemperature')
                if addTol:
                    f2_ax1.fill_between(upptarg_plot.index, upptarg_plot, 
                                        lowtarg_plot,
                                      facecolor='0.3', alpha=0.3)
            
            #f2_ax1.set_ylabel(f'Water Temperature, {plun}',fontsize=13)
            f2_ax1.legend(bbox_to_anchor=[0.22, 0.92], frameon=False, ncol=1)
            f2_ax2.plot(diff , label='Diff',c='k')
            f2_ax2.axhline(0,ls='-',c='0.3', lw=0.8, zorder=0)
            f2_ax2.axhline(1, ls='--',c='0.5',lw=0.3,zorder=0)
            f2_ax2.axhline(-1, ls='--', c='0.5',lw=0.3, zorder=0)
            
            
            if 'year' in kwargs:
                f2_ax1.xaxis.set_major_locator(months)
                f2_ax1.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=[1,8,15,22]))
                f2_ax2.xaxis.set_major_locator(months)
                f2_ax2.xaxis.set_minor_locator(mdates.DayLocator(bymonthday=[1,8,15,22]))
                # f2_ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1,4,7,10])) #interval=3))
                # f2_ax2.xaxis.set_minor_locator(months)
                f2_ax1.xaxis.set_major_formatter(monyr_fmt)
                f2_ax2.xaxis.set_major_formatter(monyr_fmt)   
            else:
                f2_ax1.xaxis.set_major_locator(mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
                f2_ax1.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7,10]))
                f2_ax2.xaxis.set_major_locator(mdates.YearLocator()) #mdates.MonthLocator(bymonth=[1,4,7,10]))
                f2_ax2.xaxis.set_minor_locator(mdates.MonthLocator(bymonth=[1,4,7,10]))
                f2_ax1.tick_params(axis='x', labelsize=11 )
                f2_ax2.tick_params(axis='x', labelsize=11 )
                f2_ax1.xaxis.set_major_formatter(years_fmt)
                f2_ax2.xaxis.set_major_formatter(years_fmt)   
                f2_ax1.xaxis.set_tick_params(rotation=75)
                f2_ax2.xaxis.set_tick_params(rotation=75)
                
   
            if plot_units=='degC':
                ax2 = f2_ax1.secondary_yaxis( -0.10, functions=(c_to_f, f_to_c))
                ax2.yaxis.set_major_locator(MultipleLocator(2))
                ax2.yaxis.set_minor_locator(MultipleLocator(1))
                ax2.set_ylabel(f'Water Temperature, $^oF$',fontsize=13)
                f2_ax1.set_ylabel('$^oC$')
                ax2b = f2_ax2.secondary_yaxis( -0.10, functions=(c_to_f_diff, f_to_c_diff))

                f2_ax1.yaxis.set_major_locator(MultipleLocator(1))
                f2_ax1.yaxis.set_minor_locator(MultipleLocator(0.5))
                ax2.yaxis.set_minor_locator(AutoMinorLocator())
                
                ax2b.set_ylabel(f'Diff in Water Temp\n(Sim-Obs),$^oF$', fontsize=12)
                ax2b.yaxis.set_minor_locator(AutoMinorLocator())
                f2_ax2.yaxis.set_minor_locator(AutoMinorLocator())
                f2_ax2.set_ylabel(f'{plun}', fontsize=12)
            else:
                f2_ax1.yaxis.set_major_locator(MultipleLocator(2))
                f2_ax1.yaxis.set_minor_locator(MultipleLocator(1))
                f2_ax1.set_ylabel(f'Water Temperature, {plun}',fontsize=13)
                
            plt.subplots_adjust(left=0.19, bottom=0.15)

            #f2_ax1.set_ylabel('Water Temperature, $^oF$',fontsize=13)
            
            sns.despine()
            if 'year' in kwargs:
                plt.suptitle(f'Simulated and Observed Water Temperatures\n{selnode.upper()} - {sel_year} ',
                             fontweight='bold')
            else:
                plt.suptitle(f'Simulated and Observed Water Temperatures\n{selnode.upper()}',
                             fontweight='bold')
            #plt.subplots_adjust(hspace=0.4)
        return(fig2)
            
    def gridify(rivmod):
        # convert results into a grid suitable for temperature and tdm landscap
        # analysis
        ## TODO
        return
    
    def combine_to_df(rivmod):
        # combine results at nodes into a single dataframe
        nn = [n for n in rivmod.Nodes]
        df = pnd.DataFrame(index=rivmod.SimDates) #, columns = rivmod.Nodes)
        for n in rivmod.Nodes:
            print(f"adding node {n} to dataframe ")
            d = rivmod.Nodes[n].OutflowTemp
            if rivmod.Nodes[n].UpNode<0:
                # use boundary condition data
                c = rivmod.Nodes[n].BoundaryConditions.ColumnMap['uptemp']
                d=rivmod.Nodes[n].BoundaryConditions.DataFrame.loc[rivmod.SimDates, c]
                df[n] = d
            else:
                if len(rivmod.Nodes[n].OutflowTemp)>len(df):
                    df[n] = rivmod.Nodes[n].OutflowTemp[1:]
                else:
                    df[n] = rivmod.Nodes[n].OutflowTemp
                
            if rivmod.Nodes[n].Observations.DataFrame is not None:
                obs_units = rivmod.Nodes[n].Observations.ColumnMap['units']
                obs_data = rivmod.Nodes[n].Observations.DataFrame.copy(deep=True)
                if obs_units == 'degF':
                    obs_data[obs_data.columns[0]] = [(i-32)/1.8 for i in obs_data[obs_data.columns[0]]]
                df[f'{n}_obs'] = obs_data.loc[:, obs_data.columns[0]]
                
        return(df)
    
    def calc_error_metrics(rivmod, df, node): #, monthly=False):

        error_dict = {}
        yrs = list(set(df.index.year))
        
        for y in yrs:          

            sim_data = df.loc[str(y), [node]] #['sim']
            obs_data = df.loc[str(y), [f'{node}_obs']]
                                                       
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
            error_dict[y] = month_stats
                
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
            
            # else:
               
            #     this_sim_temp = sim_data.loc[f'{y}', sim_data.columns[0]]
            #     this_obs_temp = obs_data.loc[f'{y}', obs_data.columns[0]]
        
            #     err_temp = [s-o for s,o in zip(this_sim_temp, this_obs_temp)]
            #     sumsqerr = np.sum([er*er for er in err_temp])
                
            #     nse = 1 - sumsqerr/(np.nansum((this_obs_temp-np.nanmean(this_obs_temp))**2))
                
            #     mean_bias = np.sum(err_temp)/len(err_temp)
            #     mae = np.sum([abs(i) for i in err_temp])/len(err_temp)
                
            #     rmse = np.sqrt(sumsqerr/len(this_obs_temp))
            #     corrmat = np.corrcoef(this_obs_temp,pnd.to_numeric(this_sim_temp))
            #     corrxy = corrmat[0,1]
            #     r2 = corrxy**2
            #     error_dict[y] = {'nse': nse, 'rmse': rmse, 'r2':r2, 'mean_bias': mean_bias,
            #                      'mae': mae}
            
        return(error_dict)
                
    class Node:
        
        def __init__(self):
            
            self.ID = ''
            self.UpNode = None
            self.UpNodeName = ''
            self.DownNode = None
            self.Inflow = None
            self.Outflow = None
            self.OutflowTemp = []
            self.Coefficients = None
            self.BoundaryConditions = inputTS()
            self.Observations = inputTS()
            # self.WatTempCoeff = {t:np.nan for t in range(1,13)}
            # self.AirTempCoeff = {t:np.nan for t in range(1,13)}
            # self.FlowCoeff = {t:np.nan for t in range(1,13)}
            # self.Intcpt = {t:np.nan for t in range(1,13)}
        
