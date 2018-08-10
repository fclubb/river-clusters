# create some fake slope area plots to illustrate the clustering

import numpy as np
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import fcluster
import pandas as pd
import math


def slope_from_drainage_area(areas, ks, theta):
    slopes = [ks*math.pow(x, -theta) for x in areas]
    #ks * math.pow(areas, theta)
    return slopes

#fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(5,6), sharey=True)

# make a big subplot to allow sharing of axis labels
#fig.add_subplot(111, frameon=False)
# hide tick and tick label of the big axes
#plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

areas = np.linspace(1000,1000000,1000)
theta = 0.5
ks_values = [1, 10, 100]

for i,ks in enumerate(ks_values):
    rand_areas = np.random.choice(areas,size=500)
    slopes = slope_from_drainage_area(rand_areas, ks, theta)
    # ax[i].scatter(areas,slopes)
    # ax[i].set_xscale('log')
    # ax[i].set_yscale('log')
    plt.scatter(rand_areas,slopes)

# set axis labels
plt.xlabel('Drainage area (m$^2$)')
plt.ylabel('Gradient', labelpad=15)
plt.ylim(0.0001,10)
plt.xscale('log')
plt.yscale('log')
plt.show()
