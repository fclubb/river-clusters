# river-clusters

This repository contains code for clustering rivers in different catchments based on their long profiles. These instructions are for install on a linux machine. For help installing on other platforms, feel free to contact Fiona Clubb (details below).

## Dependencies

In order to run the clustering code we first of all need to extract all of the river profiles from a catchment using LSDTopoTools.  To do this, first clone the LSDTopoTools2 git repository:

```
git clone https://github.com/LSDtopotools/LSDTopoTools2
cd LSDTopoTools2
```
Then install LSDTopoTools2 by running the following within the `LSDTopoTools2` directory:
```
bash lsdtt2_setup.sh
```
This will have installed the river clustering tool onto your machine.

## Download and install

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
Deactivate the environment at any time with
```
source deactivate river-clusters
```
You can also run the code without conda by ensuring that you have all of the dependencies installed e.g. via pip (see list in the `environment.yml` file.

## Look at the example data

We have provided an example run in the directory `example_data`. This DEM is from a landscape evolution model with a simple contrast in erodibility between the N and S half of the raster. To run the river cluster code using `LSDTopoTools` you need to provide a parameter file. You can find the parameter file for the example data at: `example_data/river_clusters.driver`. It should look like this:

```
# This is a driver file for LSDTopoTools
# Any lines with the # symbol in the first row will be ignored

# File information
dem read extension: bil
dem write extension: bil
read path: /home/clubb/GitHub/river-clusters/example_data/
write path: /home/clubb/GitHub/river-clusters/example_data/
read fname: spatial_K
write fname: spatial_K
threshold_contributing_pixels: 1000

# Parameters for DEM processing
raster_is_filled: false
min_slope_for_fill: 0.001

# What basins do you want?
select_basins_by_order: true
basin_order: 5

# What analyses you want to do
print_all_tributaries_to_csv: true

# Parameters to print data
print_filled_raster: true
print_basin_raster: true
print_junctions_to_csv: true
write hillshade: true
print_channels_to_csv: true
print_stream_order_raster: true
print_junction_index_raster: true
```
**NOTE: You MUST edit the read path and the write path to point to the correct directory on your machine!**

This parameter file should be saved in the *SAME DIRECTORY* as your DEM.

The parameter `basin_order` controls which catchments are selected for clustering.  For the example data, we chose 5th order so that we just get the one basin in the domain. You can adjust this based on your specific requirements.

## Run LSDTopoTools to get the river profiles

Go into the directory with the example data and then run the river profile clustering:

```
cd example_data
lsdtt-river-clusters ./ river_clusters.driver
```
This always follows the format `lsdtt-river-clusters <path-to-data> <name-of-DEM>`

This should have produced the file `DEM_name_all_tribs.csv` along with some other rasters (such as the stream network and a hillshade, for example). So for the example data, your directory should now look like this:

```
ls
```

## Run the clustering code

To run the clustering code, go back to the main repository and activate the river clustering environment (if using conda):
```
source activate river-clusters
```
Then run:
```
python cluster-river-profiles.py -h
```
which will bring up a help menu with the options.  The standard format is:
```
python cluster-river-profiles.py -dir </path/to/data/folder/> -fname <DEM_name> -so <stream_order_of_choice>
```
For the example data, this command would look like:
```
python cluster-river-profiles.py -dir ./example_data/ -fname spatial_K -so 1
```

## Contact

For more information please contact Fiona Clubb at the University of Potsdam: `clubb@uni-potsdam.de`
