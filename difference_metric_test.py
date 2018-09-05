# test of the normalised difference metric on synthetic
# data

import numpy as np
import matplotlib.pyplot as plt

def difference(X,Y):
    n = len(X)
    # now calculate the difference
    numerator = X - Y
    denominator = X + Y
    div = numerator/denominator
    norm = np.linalg.norm(div)
    diff = norm/np.sqrt(n)
    return div, diff

# these are the two profiles
x_slope = -0.3
y_slope = -0.8
x_intercept= 15
y_intercept= 15
dist = np.linspace(1, 100, 100)
# X = x_slope*dist + x_intercept
# Y = y_slope*dist + y_intercept

# create two different power laws
X = x_intercept * dist**(x_slope)
Y = y_intercept * dist**(y_slope)
n = len(X)

div, diff = difference(X,Y)
# plot the data
fig, ax = plt.subplots(nrows = 2, ncols =1)
ax = ax.ravel()
ax[0].set_title('Different slopes, d = {}'.format(diff))
ax[0].scatter(dist, X, label = 'X', c= 'b',s=2)
ax[0].scatter(dist, Y, label = 'Y', c= 'orange',s=2)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
# plot (X-Y)/(X+Y)
ax2 = ax[0].twinx()
ax2.plot(dist,div, label='diff', c = 'k', ls='--')
#
#make another profile which is parallel to Y but offset by 5 for testing
Z = (y_intercept + 5) * dist**y_slope
div, diff = difference(Z,Y)
#
ax[1].set_title('Different intercepts, d = {}'.format(diff))
ax[1].scatter(dist, Z, label = 'Z', c='purple',s=2)
ax[1].scatter(dist, Y, label = 'Y', c='orange',s=2)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
# plot (X-Y)/(X+Y)
ax2 = ax[1].twinx()
ax2.plot(dist,div, label='diff', c = 'k', ls='--')

# add a legend
ax[0].legend(loc='upper right')
plt.show()

# now subsample Y and Z. does the difference change?
plt.clf()
fig, ax = plt.subplots(nrows = 2, ncols =1)
ax = ax.ravel()
ax[0].set_title('Different intercepts, all data: d = {}'.format(diff))
ax[0].scatter(dist, Z, c= 'b',s=2)
ax[0].scatter(dist, Y, c= 'orange',s=2)
ax[0].set_xscale('log')
ax[0].set_yscale('log')
# plot (X-Y)/(X+Y)
ax2 = ax[0].twinx()
ax2.plot(dist, div, label='diff', c = 'k', ls='--')

# do the subsampling and calculate the difference
sample_dist = dist[0:-1:5]
sample_Z = Z[0:-1:5]
sample_Y = Y[0:-1:5]
div, diff = difference(sample_Z, sample_Y)


ax[1].set_title('Different intercepts, subsampling: d = {}'.format(diff))
ax[1].scatter(sample_dist, sample_Z, c='b',s=2)
ax[1].scatter(sample_dist, sample_Y, c='orange',s=2)
ax[1].set_xscale('log')
ax[1].set_yscale('log')
# plot (X-Y)/(X+Y)
ax2 = ax[1].twinx()
ax2.plot(sample_dist, div, label='diff', c = 'k', ls='--')
plt.show()
