# -*- coding: utf-8 -*-
"""
Created on Wed Jul  1 16:42:16 2020

@author: jgilbert
"""

import numpy as np

def DEGF_to_DEGC(x):
    return((x-32)/1.8)

def DEGC_to_DEGF(x):
    return(x*1.8+32.)
#%% UNIT CONVERSION FACTORS
ACCEPTABLE_LENGTH_UNITS = ['m', 'meter','meters', 'ft','feet']

FT3toM3 = 12.*12.*12.*2.54*2.54*2.54/100./100./100.
AFtoM3 = 43560.*FT3toM3
MtoFT = 100./2.54/12.
FTtoM = 1/MtoFT
ACREtoM2 = 43560.*12*12*2.54*2.54/100/100.

CMStoCFS = 100.*100.*100./(2.54*2.54*2.54)/(12.*12.*12)
CMStoAFD = CMStoCFS*86400/43560.
CFStoCMS = 1./CMStoCFS
CFStoAFD = 86400/43560.
AFDtoCMS = (43560./86400) * CFStoCMS
ELEV_BUFFER = 0.0001

# for precip: kg/m2/s to feet per day - assumes precip density of 1000 kg/m3
KGpM2StoFTpD = 86400./1000.*MtoFT

# time increments
DELT_DAY = 1.
DELT_SEC = DELT_DAY*86400.

# specific heat of water
cwater = 4182. #J/kg*C

# for solar insolation/rdiation calcs
Wm2ToCalcm2pDy = 86400./((100**2)*4.1868)  # seconds/day / (cm2 per 1m2 * Watts per cal/second)

# volume tolerance
TOL_AF = 0.01
DENSE_TOL = 1e-4
DENSE_PREC = int(abs(np.log10(DENSE_TOL)))+1

# temperature tolerance for finding blending ratios
TOL_TEMP = 1e-6



CFStoKG= 1/CMStoCFS * 86400.*1000. # convert from cfs to cubic meters/sec, seconds to days, to kg using assumed density of 1000. kg/m3
Le = 2.453e6   # latent heat of vaporization, J/kg, at 20deg C
cw = 4182.  # specific heat of water at 20 deg C, in J/kg/K

BOWEN_CONSTANT = 0.47


TWOPI = np.pi*2.
HALFPI = np.pi*0.5