# -*- coding: utf-8 -*-
"""
An example script for running a ResTemp model for Shasta Reservoir

This example runs a historic simulation with prescribed gate operations 
but without optimal blending

@author: jgilbert
"""

# import some standard libraries
import sys, os
import pandas as pd
import numpy as np

# if not working from installed package uncomment `sys.path.insert...` 
#  command on next line, replacing path/to/R3TAM/src with your local filepath
#sys.path.insert(0, r'C:\path\to\R3TAM\src')

# import the R3TAM libraries - since this is just a Shasta reservoir 
# simulation, we really only need the ResTemp module
from r3tam import restemp as rt
#from r3tam import longtemp as lt
#from r3tam.river import River
#from r3tam import coupling

# the `make_plots` library has helpful plots for viewing results
from r3tam import make_plots as mkp

# the `time` module is only for measuring the timing of a simulation
# and is not strictly required
import time

#%% Initialize the model

# Set the path to the run specification/configuration file - `shasta_standalone.yaml`
config_fp = r'D:\02_Projects\SacTemp\SimTemp\R3TAM\examples\shasta\shasta_standalone.yaml'

shasta = rt.Res.initialize_model(config_fp, profile_temp_units='degF')

#%% Run the model

btime = time.perf_counter() 
for d in shasta.SimDates:
    print(d)
    shasta.advance_restemp()

    shasta.advance_swd(final=True)

etime = time.perf_counter() 
print(f"took {etime-btime} seconds")


#%%  Get the results and make some plots!

# `finalize` organizes results into pandas DataFrames for easier access
shasta.finalize()

# you can get many of the time-series results as DataFrames, such as for
# simulated releases and temperature profiles, below
simReleases = shasta.Simulation_Results['ReleaseDF']
simProf = shasta.Simulation_Results['ProfilesDF']

# the `make_plots` module has some convenient plotting options
mkp.plotReleasesCompare(shasta, select_years=[2014, 2015], 
                        viewSave='save', 
                        obs_label = 'Obs Tailwater',
                        other_temp = 'temperature_other',
                        other_label = 'TCD Wt Avg',
                        on_wy=False) 
mkp.plotProfilesCompare2(shasta, simProf, viewSave='save', on_wy=False)
