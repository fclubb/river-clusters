# Plot some random time series for a presentation
import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as hac
import pandas as pd

#%% create elevation profiles using simple exponential decay, but vary decay constant
# also plot all profiles
length_of_profile = 1000
nr_of_profiles = 5
start_elevation = 200
noise_level = start_elevation / 20 #noise level at 5% of start_elevation
distance = np.arange(0,length_of_profile,1)
exponential_decay_factor = np.random.randint(low=1, high=100, size=(nr_of_profiles))
elevation_profile_matrix = np.empty((nr_of_profiles, length_of_profile))
print (elevation_profile_matrix.shape)

fig, ax = plt.subplots(nrows=nr_of_profiles,ncols=1,sharex=True,sharey=True, figsize =(5,8))
fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')
colors = ['b', 'r', 'g', 'k', 'orange']

for i in range(nr_of_profiles):
    elevation_perf = [start_elevation * np.exp(-t / exponential_decay_factor[i]) for t in distance] # exponential perfect decay
    noise = np.random.uniform(0, noise_level, size = len(elevation_perf)) # add some noise
    elevation_noise = elevation_perf + noise
    elevation_profile_matrix[i,:] = elevation_noise
    ax[i].plot(distance, elevation_profile_matrix[i,:], c=colors[i])

plt.title('Synthetic River Profiles', fontsize=16)
plt.xlabel('Distance along profile',fontsize=16)
plt.ylabel('Elevation',fontsize=16)

#plt.show()
plt.savefig('/home/s0923330/Bitbucket/HRT_workshop/figures/random_river_profiles.png', dpi=300)
plt.clf()

#%% Pearson Correlation Analysis of matrix (Pearson product-moment correlation coefficients)
elevation_profile_matrix_corrcoef = np.corrcoef(elevation_profile_matrix)

# When using larger datasets, it is useful to save
#print('saving %s_corr.npy' % pfix)
#np.save('%s_corr.npy' % pfix, elevation_profile_matrix_corrcoef)

print('arccos of correlation coefficient for threshold')
d = np.arccos(elevation_profile_matrix_corrcoef)

#%% Linkage calculation (other package, look at method and metric options
Z = hac.linkage(elevation_profile_matrix, method='ward')

# Plot dendogram
plt.figure(figsize=(25, 10))
plt.title('Hierarchical Clustering Dendrogram')
plt.xlabel('sample index')
plt.ylabel('distance')
hac.dendrogram(
    Z,
    leaf_rotation=90.,  # rotates the x axis labels
    leaf_font_size=8.,  # font size for the x axis labels
)
plt.show()
