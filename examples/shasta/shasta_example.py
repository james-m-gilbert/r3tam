# -*- coding: utf-8 -*-
"""
Created on Fri Mar 14 09:36:07 2025

@author: jgilbert
"""
import sys, os
import pandas as pd
import numpy as np

# if not working from installed package uncomment sys.path.insert command
sys.path.insert(0, r'D:\02_Projects\SacTemp\SimTemp\R3TAM\src')
from r3tam import restemp as rt
from r3tam import longtemp as lt
from r3tam.river import River
from r3tam import coupling

from r3tam import make_plots as mkp

import time

#%%

config_fp = r'D:\02_Projects\SacTemp\SimTemp\R3TAM\examples\shasta\20250314_shasta_only_report_0022_withLPblending.yaml'
#config_fp = r'D:/02_Projects/SacTemp/SimTemp/R3TAM/examples/shasta/20250225_shasta_only_report_0022.yaml'
#config_fp = r'D:/02_Projects/SacTemp/SimTemp/ResTempPkg/example/shasta_only/20230502_shasta_only_report_0022continuous.yaml'

make_plots = True
write_stats = False

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
print(f"took {etime-btime} seconds?")

resmod.finalize()
#%%
simReleases = resmod.Simulation_Results['ReleaseDF']
simProf = resmod.Simulation_Results['ProfilesDF']
mkp.plotReleasesCompare(resmod, simReleases, viewSave='view', on_wy=False) #, target_ts=simReleases.Sim_Temp_Target_degF)
mkp.plotProfilesCompare2(resmod, simProf, viewSave='view', on_wy=False)
