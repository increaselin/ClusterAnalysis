#! /usr/bin/env python
# 
# Program: DAOSpecConsts.py
#
# Creator: Mike Lum
#
# Usage: n/a
#
# Description: Contains constants in accessing the DAOSpec software
#
# Revision History:
#    Date          Vers.    Author        Description
#    2018-04-05    1.0f2    Lum         github format and comment cleanup
#    2017-02-05    1.0f1    Lum         First checked in
#
#
#
from constants import directories as dirs

DAOSpecScriptHead = '# !/bin/bash\n#\n'

DAOSpecCommand = 'exec {0}/daospec'.format(dirs.DAOSpecLoc)

# DAOSpec options list (Use the first two characters to set an option):
#
# (FW)HM - Estimate of spectral resolution (in pixels)
# (OR)der for continuum fit
# (SH)ort wavelength limit - Blue limit of region to measure
# (LO)ng wavelength limit - Red limit   "
# (LE)ft edge of window - Blue limit of window to watch fitting function
# (RI)ght edge of window - Red limit   "
# (WA)tch fit process (0=False, 1=True)
# (BA)d data value - Lower limit of good data
# (RE)sidual core flux - Lower limit of core flux for unsaturated lines
# (VE)locity limit - Lines shifted by more than this amount are ignored
# (MI)nimum radial velocity - Limit Rv calculation to greater than this limit
# (MA)ximum radial velocity - "     "      "        "    less  "    "    "
# (SM)allest equivalent width (mA) - Lower limit for valid eqw measures
# (CR)eate output spectra - Produce fitted continuum and residual spectra files
# (SC)ale FWHM with wavelength
# (FI)x FWHM - Require the FWHM (at the center WL) to be the input value

DAOSpecOpts = '''
CR=0
FI=0
OR=5
SC=1
SM=3
VE=5
WA=0
'''


