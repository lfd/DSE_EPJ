import os
from os import listdir
from os.path import isfile, join
from typing import overload
from numpy import empty
from openpyxl import load_workbook
import pandas as pd
import math
import seaborn as sns
import numpy as np
from gmpy2 import mpz, mpfr
import sys
import ast

#get fidelity decrease for all the data retrieved from the compiler
def get_fidelity_decrease(compiler_data_path, oneQ_fid = 0.9982,twoQ_fid = 0.9765):

    fidelity_decrease_all = []
    for i,row in compiler_data_path.iterrows(): 
            no_of_1_q_gates_before = row['gates before'] - row['2qgates before']
            no_of_more_q_gates_before = row['2qgates before']
            no_of_more_q_gates_after = row['2qgates']
            no_of_1_q_gates_after = row['gates'] - row['2qgates']
            fidelity_before = (math.pow(oneQ_fid,no_of_1_q_gates_before)) * (math.pow(twoQ_fid,no_of_more_q_gates_before))
            fidelity_after = (math.pow(oneQ_fid,no_of_1_q_gates_after)) * (math.pow(twoQ_fid,no_of_more_q_gates_after))
            fidelity_decrease = ''
            if fidelity_before !=0 :
                fidelity_decrease = (fidelity_before - fidelity_after)/fidelity_before
            fidelity_decrease_all.append(fidelity_decrease)
    #df = pd.DataFrame(fidelity_decrease_all, columns=['fidelity_decrease'])
    #df.to_excel(curdir + "/fidelity_decrease.xlsx", index = False)
    return fidelity_decrease_all

#get depth overhead for all the data retrieved from the compiler
def get_depth_overhead(compiler_data_path):
    gate_overhead_all = []
    for i,row in compiler_data_path.iterrows():     
         gates_after = row['gates']
         gates_before = row['gates before']
         gate_overhead = ((gates_after-gates_before)/gates_before)
         gate_overhead_all.append(gate_overhead)
    return gate_overhead_all
         

     
#get gate overhead for all the data retrieved from the compiler
def get_gate_overhead(compiler_data_path):
    depth_overhead_all = []
    for i,row in compiler_data_path.iterrows():
        depth_before = row['depth before']
        depth_after = row['depth']
        depth_overhead = ((depth_after-depth_before)/depth_before)
        depth_overhead_all.append(depth_overhead)
    return depth_overhead_all

#get overall cost for all the data retrieved from the compiler
def get_cost_data(compiler_data_path, oneQ_fid = 0.9982,twoQ_fid = 0.9765, K=0.995):
    C_all = []
    for i,row in compiler_data_path.iterrows():
        depth_before = row['depth before']
        depth_after = row['depth']
        no_of_1_q_gates_before = row['gates before'] - row['2qgates before']
        no_of_more_q_gates_before = row['2qgates before']
        no_of_more_q_gates_after = row['2qgates']
        no_of_1_q_gates_after = row['gates'] - row['2qgates']
        C_in = -depth_before*np.log(K) - no_of_1_q_gates_before*np.log(oneQ_fid) - no_of_more_q_gates_before*np.log(twoQ_fid)
        C_out = -depth_after*np.log(K) - no_of_1_q_gates_after*np.log(oneQ_fid) - no_of_more_q_gates_after*np.log(twoQ_fid)
        C = C_in / C_out
        C_all.append(C)
    return C_all

if __name__ == '__main__':
     
    curdir = os.path.dirname(__file__)

    compiler_data_path = os.path.join(curdir,"DSE_results")
    cost_data = get_cost_data(compiler_data_path)
    df = pd.DataFrame(fidelity_decrease_all, columns=['cost_improvement'])
    df.to_excel(curdir + "/FOM_all_benchmarks.xlsx", index = False)
