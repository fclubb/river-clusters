import numpy as np
import matplotlib.pyplot as pl
from scipy.cluster.hierarchy import dendrogram, linkage
np.set_printoptions(threshold='nan')
# number of time series
n = 123

# random time series of random length
ts = []

for i in range(n):
    # random length rlen
    rlen = np.random.randint(100, 300)
    ts.append(np.random.random(rlen))

# correlation coefficients
cc = np.zeros(int((n * (n - 1)) / 2))
print len(cc)

k = 0
for i in range(n):
    for j in range(i+1, n):
        tsi = ts[i]
        tsj = ts[j]
        if len(tsi) > len(tsj):
            tsi = tsi[:len(tsj)]
        else:
            tsj = tsj[:len(tsi)]
        dts = tsi - tsj
        l = 0
        while dts[l] == 0:
            l += 1
        tsi, tsj = tsi[l:], tsj[l:]
        cc[k] = np.corrcoef(tsi, tsj)[0, 1]
        k += 1

print cc
dd = np.arccos(cc)
print len(dd)
ln = linkage(dd, method = 'complete')
dendrogram(ln)
pl.show()
