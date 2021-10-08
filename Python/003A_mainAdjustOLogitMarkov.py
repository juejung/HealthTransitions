#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# [0] Into {{{1
# -----------------------------------------------------------------------------
"""
Markov Switching Matrices

Case with 6 health states that includes death
6 types 
 - (male/female)
 - (black/non-black)
 - (smoker/non-smoker)
 
J.Jung
jjung@towson.edu

08/07/2020

This is the version with 6 health states, where the 6th state is death.
"""
import os
import importlib
import numpy as np
import scipy.linalg as sp
import pandas as pd
import matplotlib.pyplot as plt

import seaborn as sns
sns.set_style("whitegrid", {'grid.linestyle': ':'})
c2, c1, c4, c3, c5, c8, c9, c6, c10, c7 = sns.color_palette("Set1", 10)

import f_stochroot as myfunc
importlib.reload(myfunc)

# Suppress scientific notation in numpy arrays
np.set_printoptions(suppress=True)

title_size = 24
subtitle_size = 14
tick_size = 16
label_size = 14
legend_size = 12
text_size = 14

print('--- Start Plot ---')
# This is important, so that we can call this script from within a Stata do file
os.chdir('/home/jjung/Dropbox/AAPapers/Jung/HealthMarkov/Python/')
path = os.getcwd()
print(path)


plt.close()
exec(open("001_Para.py").read())

# Data: MEPS-HRS 2000- from: 102_oLogitPlot_MEPS_vs_HRS.do
dataIn = pd.read_csv(gStataOutDir \
         + 'MEPS_HRS_1992-2017_Age_All_Dead_OLogit.csv', delimiter=',')

nrRow = 5
nrCol = 6
nrAge = 16
nrGend = 7  
#}}}1
# [2] Transform 2-year Markov matrices into 1-year Markov matrices {{{1
# -----------------------------------------------------------------------------
temp_mat = np.zeros((nrRow+1, nrCol, nrAge, nrGend))
temp_mat_adj = np.zeros((nrRow+1, nrCol, nrAge, nrGend))
typeList = ['a', 'f', 'm', 'b0', 'b1', 's0', 's1']

for i_age in range(nrAge):

    print(dataIn['agev'][i_age])

    for i_type, type_str in enumerate(typeList):

        # Build markov transition matrix for i_age
        for i_row in range(nrRow):
            for i_col in range(nrCol):

                colName = type_str + '_hrs_' + str(i_row+1) + '_to_pr' + str(i_col+1)
                temp_mat[i_row, i_col, i_age, i_type] = dataIn[colName][i_age]

            #end
        #end
        temp_mat[5,5] = 1
        temp_mat_adj[:,:,i_age, i_type] = myfunc.f_stochroot(temp_mat[:,:,i_age, i_type], 1./2)
    #end
#end

# print(temp_mat[:,:,0])
# print(temp_mat_adj[:,:,0])

# print(temp_mat[:,:,10])
# print(temp_mat_adj[:,:,10])
#}}}1
# [3] Add new adjusted matrices to dataframe {{{1
# -----------------------------------------------------------------------------
# Make new columns for adjusted values
for i_type, type_str in enumerate(typeList):
    for i_row in range(nrRow+1):
        for i_col in range(nrCol):
            colName = type_str + '_hrs_' + str(i_row+1) + '_to_pr' + str(i_col+1) + '_adj'
            dataIn[colName] = 0.0
        #end
    #end
#end

for i_age in range(nrAge):

    print(dataIn['agev'][i_age])
    for i_type, type_str in enumerate(typeList):
        for i_row in range(nrRow):
            for i_col in range(nrCol):
                colName = type_str + '_hrs_' + str(i_row+1) + '_to_pr' + str(i_col+1) + '_adj'
                dataIn.loc[i_age, colName] = temp_mat_adj[i_row, i_col, i_age, i_type]
            #end
        #end
    #end
#end


dataIn.to_csv(gStataOutDir + 'MEPS_HRS_1992-2017_Age_All_Dead_OLogit_Adj.csv')
#}}}1