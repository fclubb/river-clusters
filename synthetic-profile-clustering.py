#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  1 10:31:24 2018

@author: bodo
"""
# install TimeSynth (if you want more fancy and irregularly sampled signals) with:
#git clone https://github.com/TimeSynth/TimeSynth.git
#cd TimeSynth
#python setup.py install
#import timesynth as ts

# For larger matrices, use:
#Make sure to install CorrCoef:
#git clone https://github.com/UP-RS-ESP/CorrCoef
#cd CorrCoef
#python setup.py install
#import CorrCoef

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as hac
import pandas as pd

#%% create elevation profiles using simple exponential decay, but vary decay constant
# also plot all profiles
length_of_profile = 1000
nr_of_profiles = 100
start_elevation = 200
noise_level = start_elevation / 20 #noise level at 5% of start_elevation
distance = np.arange(0,length_of_profile,1)
exponential_decay_factor = np.random.randint(low=1, high=100, size=(nr_of_profiles))
elevation_profile_matrix = np.empty((nr_of_profiles, length_of_profile))
print (elevation_profile_matrix.shape)

plt.clf()
for i in range(nr_of_profiles):
    elevation_perf = [start_elevation * np.exp(-t / exponential_decay_factor[i]) for t in distance] # exponential perfect decay
    noise = np.random.uniform(0, noise_level, size = len(elevation_perf)) # add some noise
    elevation_noise = elevation_perf + noise
    elevation_profile_matrix[i,:] = elevation_noise
    plt.plot(distance, elevation_profile_matrix[i,:])

plt.grid()
plt.title('Synthethic River Profiles (n=%02d)'%(nr_of_profiles))
plt.xlabel('Distance along profile')
plt.ylabel('Elevation')

#%% Pearson Correlation Analysis of matrix (Pearson product-moment correlation coefficients)
elevation_profile_matrix_corrcoef = np.corrcoef(elevation_profile_matrix)

# When using larger datasets, it is useful to save
#print('saving %s_corr.npy' % pfix)
#np.save('%s_corr.npy' % pfix, elevation_profile_matrix_corrcoef)

print('arccos of correlation coefficient for threshold')
d = np.arccos(elevation_profile_matrix_corrcoef)


#%% Linkage calculation (using arccos(pearson correlation))
method = 'complete'
print('linkage method=%s'%method)
Z = linkage(d, method=method)

# When using larger datasets, it is useful to save
#print('saving %s_linkage.npy' % pfix)
#np.save('%s_linkage.npy' % pfix, z)

plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram - Complete & Pearson')
plt.xlabel('sample index')
plt.ylabel('distance')
hac.dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()

#%% Linkage calculation (other package, look at method and metric options
Z = hac.linkage(elevation_profile_matrix, method='complete', metric='correlation')

# Plot dendogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram - Complete')
plt.xlabel('sample index')
plt.ylabel('distance')
hac.dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()

#%% Perform clustering based on linkages
# may be useful: https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
def print_clusters(results, elevation_profile_pd):
    # check the results
    s = pd.Series(results)
    clusters = s.unique()

    for c in clusters:
        cluster_indeces = s[s==c].index
        print("Cluster %d number of entries %d" % (c, len(cluster_indeces)))
        elevation_profile_pd.T.iloc[:,cluster_indeces].plot()
        plt.title('Cluster %02d of %02d'%(c, len(clusters)))
        plt.show()

#make sure to use complete method:
Z = linkage(d, method=method)

#create cluster based on nr_of_clusters
nr_of_cluster = 6
clst_maxcluster = fcluster(Z, nr_of_cluster, criterion='maxclust')
elevation_profile_pd = pd.DataFrame(elevation_profile_matrix)
print_clusters(clst_maxcluster, elevation_profile_pd)


#create distance-based cluster
#clst = fcluster(Z, d, criterion='distance')
