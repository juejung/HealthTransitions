#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#%% [0] Intro {{{1
"""
Juergen Jung
Towson University
Email: jjung@towson.edu
01/16/2017
"""

import socket
import numpy as np
from pathlib import Path

#}}}1
#%% [1] Filepath {{{1
hostName = socket.gethostname()
home = str(Path.home())

if hostName ==  'BR962':
    # Windows
    print("Run on Windows ...")
    gMEPSdataDir   = "C:/Dropbox/AAA/Data/MEPS/MEPS_Annuals/"
    gHRSdataDir    = "C:/Dropbox/AAA/Data/HRS/RandVersionN/"
    gCPIdir        = "C:/Dropbox/AAA/Data/Indices/"

    gDatDir   = "C:/Dropbox/AAPapers/JungBagchi/Stata/Data/"
    gStataOutDir = "C:/Dropbox/AAPapers/JungBagchi/Stata/Output/"
    gDoDir    = "C:/Dropbox/AAPapers/JungBagchi/Stata/Do/"
    gTexDir   = "C:/Dropbox/AAPapers/JungBagchi/Tex/Graphs/"
else:
    # Linux
    print("Run on Linux ...")
    gMEPSdataDir   = home + "/Dropbox/AAA/Data/MEPS/MEPS_Annuals/"
    gHRSdataDir    = home + "/Dropbox/AAA/Data/HRS/Rand_1992_2016_V2/"
    gCPIdir        = home + "/Dropbox/AAA/Data/Indices/"

    gDatDir        = home + "/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Data/"
    gStataOutDir   = home + "/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Output/"
    gDoDir         = home + "/Dropbox/AAPapers/Jung/HealthMarkov/Stata/Do/"
    gTexDir        = home + "/Dropbox/AAPapers/Jung/HealthMarkov/Tex/Graphs/"
#end

#}}}1
