# river-clusters

This repository contains code for clustering rivers in different catchments based on their long profiles.

## Dependencies

To run this code you need to have Python installed. You must also have the following packages:

* matplotlib
* scipy
* numpy
* pandas


In addition to these packages, you also need to have `CorrCoef` installed, a Python C-extension for memory efficient and multithreaded Pearson product-moment correlation coefficient estimation using OpenMP, created by Aljoscha Rheinwalt. For more details see the [GitHub page](https://github.com/UP-RS-ESP/CorrCoef).

To install:

```
git clone https://github.com/Rheinwalt/CorrCoef.git
cd CorrCoef
sudo python setup.py install
```
## Contact

For more information please contact Fiona Clubb at the University of Potsdam: `clubb@uni-potsdam.de`
