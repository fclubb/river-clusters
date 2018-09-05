# Plot some random time series for a presentation
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as hac
import pandas as pd

#%% create elevation profiles using simple exponential decay, but vary decay constant
# also plot all profiles
# length_of_profile = 1000
# # nr_of_profiles = 5
# # start_elevation = 50
# # noise_level = start_elevation / 20 #noise level at 5% of start_elevation
# distance = np.arange(0,length_of_profile,1)
# exponential_decay_factor = np.random.randint(low=100, high=200, size=(nr_of_profiles))
# elevation_profile_matrix = np.empty((nr_of_profiles, length_of_profile))
# print (elevation_profile_matrix.shape)

# create random gradient--distance plots
length_of_profile = 1000
distance = np.arange(0,length_of_profile,1)
nr_of_profiles = 5
start_gradient = 0.8
k = np.random.uniform(low=0.1, high=0.2, size=(nr_of_profiles))
t = np.random.uniform(low=-0.8, high=-0.2, size=(nr_of_profiles))
kp_loc = np.random.randint(low=10, high=(length_of_profile-10), size=(nr_of_profiles))
kp_loc2 = np.random.randint(low=10, high=(length_of_profile-10), size=(nr_of_profiles))
noise_level = start_gradient/40

fig, ax = plt.subplots(nrows=nr_of_profiles,ncols=1,sharex=True,sharey=True, figsize =(5,8))
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
colors = ['b', 'r', 'g', 'k', 'orange']

for i in range(nr_of_profiles):
    #elevation_perf = [start_elevation * np.exp(-t / exponential_decay_factor[i]) for t in distance] # exponential perfect decay
    slope_perf = [start_gradient * k[i] * d**(t[i]) for d in distance]
    noise = np.random.uniform(0, noise_level, size = len(slope_perf)) # add some noise
    print(kp_loc[i])
    noise[kp_loc[i]-5:kp_loc[i]+5] = noise[kp_loc[i]-5:kp_loc[i]+5] + 0.02
    noise[kp_loc2[i]-5:kp_loc2[i]+5] = noise[kp_loc2[i]-5:kp_loc2[i]+5] + 0.02
    #elevation_noise = elevation_perf + noise
    gradient = slope_perf+noise
    ax[i].plot(distance, gradient, c=colors[i])

#plt.title('River profiles', fontsize=16)
plt.xlabel('Distance along profile (m)',fontsize=16)
plt.ylabel('Gradient',fontsize=16)

#plt.show()
plt.savefig('/home/clubb/Bitbucket/hrt_workshop/figures/random_river_profiles.png', dpi=300, transparent=True)
plt.clf()

#%% Pearson Correlation Analysis of matrix (Pearson product-moment correlation coefficients)
# elevation_profile_matrix_corrcoef = np.corrcoef(elevation_profile_matrix)
#
# # When using larger datasets, it is useful to save
# #print('saving %s_corr.npy' % pfix)
# #np.save('%s_corr.npy' % pfix, elevation_profile_matrix_corrcoef)
#
# print('arccos of correlation coefficient for threshold')
# d = np.arccos(elevation_profile_matrix_corrcoef)
#
# #%% Linkage calculation (other package, look at method and metric options
# Z = hac.linkage(elevation_profile_matrix, method='ward')
#
# # Plot dendogram
# plt.figure(figsize=(25, 10))
# plt.title('Hierarchical Clustering Dendrogram')
# plt.xlabel('sample index')
# plt.ylabel('distance')
# hac.dendrogram(
#     Z,
#     leaf_rotation=90.,  # rotates the x axis labels
#     leaf_font_size=8.,  # font size for the x axis labels
# )
# plt.show()
