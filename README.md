# river-clusters

This repository contains code for clustering rivers in different catchments based on their long profiles.

## Dependencies

The best way to run the code is to use `conda`. Install Anaconda or Miniconda with Python 3.6: https://conda.io/miniconda.html

First of all clone the repository:
```
git clone https://github.com/fclubb/river-clusters.git
cd river-clusters
```
Then create a new conda environment using the included `environment.yml` file:

```
conda env create -f environment.yml
```
and activate the new environment by:
```
source activate river-clusters
```
You can then run the code from inside this environment using:
```
python cluster-river-profiles.py -h
```
which will bring up a help menu with the options.  The standard format is:
```
python cluster-river-profiles.py -dir </path/to/data/folder/> -fname <DEM_name> -so <stream_order_of_choice>
```
## Contact

For more information please contact Fiona Clubb at the University of Potsdam: `clubb@uni-potsdam.de`
