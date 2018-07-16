# sensitivity.py
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# Run the river cluster analysis with varying
# parameters to check sensitivity
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# FJC
# 16/07/18
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=


import os
from glob import glob
import numpy as np, matplotlib.pyplot as plt
import subprocess

DataDirectory ="/raid/fclubb/river_clusters/model_runs/runs_for_analysis/sensitivity_analyses/slope_window/''

print DataDirectory

subdirs = [x[0] for x in os.walk(DataDirectory)]

for dir in subdirs:
    sw = dir.split("/")[-1]
    print sw
    system_call = 'nohup python plot-river-clusters.py -dir '+str(dir)+' -fname fault_uplift_10m1516269869 -len 50 -sw '+str(sw)+' &'
    print system_call
    # subprocess.call(system_call, shell=True)
