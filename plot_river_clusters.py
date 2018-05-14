#---------------------------------------------------------------------#
# Clustering of river profiles
# Developed by Fiona Clubb
#              Bodo Bookhagen
#              Aljoscha Rheinwalt
# University of Potsdam
#---------------------------------------------------------------------#

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
from glob import glob
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy import stats
from CorrCoef import Pearson
import math

def read_river_profile_csv(DataDirectory, fname_prefix):
    """
    Function to read in a csv file with the river profile data

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_spaghetti_profiles.csv')
    return df

def find_nearest_idx(array,value):
    """
    Given a value, find the index of the point in the array which is closest
    to that value
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

#---------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
#---------------------------------------------------------------------#

def ResampleProfiles(df, profile_len = 100, step=1):
    """
    This function takes the dataframe of river profiles and creates an array
    that can be used for the time series clustering.  For each profile, the slope
    is assigned to a common distance step, in metres (default=1m)

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): step size that you want in metres, default = 1

    Returns: array of size (n_profiles, profile_len) containing the slopes along each profile

    Author: FJC
    """
    # create new array of regularly spaced differences
    reg_dist = np.arange(0, profile_len, step)

    # find the minimum length that the array can be (profile length/root2)
    min_length = profile_len/(math.sqrt(2))

    # loop through the dataframe and store the data for each profile as an array of
    # slopes and distances
    profiles = []
    source_ids = df['source_id'].unique()
    for i, source in enumerate(source_ids):
        this_df = df[df['source_id'] == source]
        slopes = this_df['slope'].as_matrix()
        distances = this_df['distance_from_source'].as_matrix()
        if (len(slopes) > min_length):
            profiles.append((distances, slopes))

    # now create the 2d array to store the data
    n_profiles = len(profiles)
    data = np.empty((n_profiles, profile_len))

    # loop through the profiles. For each point in the regularly spaced array,
    # find the index of the closest point in the distance array. Then use this to
    # assign the slope to the regularly spaced array
    for i, p in enumerate(profiles):
        reg_slope = []
        for d in reg_dist:
            idx = find_nearest_idx(p[0], d)
            reg_slope.append(p[1][idx])
        data[i] = reg_slope

    return data

def ClusterProfiles(df, profile_len=100, step=1, min_corr=0.5):
    """
    Cluster the profiles based on gradient and distance from source.
    Aggolmerative clustering.

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): the spacing in metres between the data points that you want, default = 1m
        min_corr (float): minimum correlation threshold for clustering

    Author: AR, FJC
    """
    data = ResampleProfiles(df, profile_len, step)

    # we could have a look at the ranks too ..
    # correlations
    cc = Pearson(data)

    # distances
    dd = np.arccos(cc)

    # do agglomerative clustering by stepwise pair matching
    # based on angle between scalar products of time series
    ln = linkage(dd, method = 'complete')

    # define threshold for cluster determination
    thr = np.arccos(min_corr)

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    dendrogram(ln)
    #dendro=dendrogram(
    #    link,
    #    leaf_rotation=90.,
    #    leaf_font_size=8.,
    #)
    plt.axhline(y = thr, color = 'r')
    plt.show()

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print cl.min(), cl.max()

def CalculateSlope(df, slope_window_size):
    """
    This function takes in a dataframe with elevation and distance
    from source data, and adds a column with the slope of each point
    fitted over a certain window size.

    Args:
        df: dataframe
        slope_window_size (int): total number of points used to calculate
        slope (INCLUDES the node of interest)

    Author: FJC
    """
    dist = df['distance_from_source']
    elev = df['elevation']
    slopes = np.empty(len(dist))

    pts_array = np.column_stack((dist,elev))

    slicer = (slope_window_size - 1)/2

    for index, x in enumerate(pts_array):
        # find the rows above and below relating to the window size
        this_slice = pts_array[index-slicer:index+slicer+1]
        if len(this_slice) != 0:
            # now regress this slice
            x = this_slice[:,0]
            y = this_slice[:,1]
            #print x, y
            slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
            slopes[index] = abs(slope)
            #print slope
        else:
            slopes[index] = np.nan

    df['slope'] = slopes

    return df

#---------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#---------------------------------------------------------------------#

def PlotAllProfilesNormalised(DataDirectory, fname_prefix):
    """
    Function to make a plot of all the normalised river profiles

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.85,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the river profile csv
    df = read_river_profile_csv(DataDirectory,fname_prefix)
    print ('Got the dataframe!')

    # normalise the distance by total length of channels upstream
    df['normalised_dist'] = df['distance_from_outlet']/df['total_length_upstream']

    # get the profile IDs
    river_ids = df['id'].unique()

    for id in river_ids:
        print ('This id is: ', id)
        this_df = df[df['id'] == id]
        max_dist = this_df['distance_from_outlet'].max()
        this_df['normalised_dist'] = this_df['distance_from_outlet']/max_dist
        outlet_elev = this_df['elevation'].min()
        this_df['normalised_elev'] = this_df['elevation']/outlet_elev
        ax.plot(this_df['normalised_dist'], this_df['normalised_elev'])

    ax.set_xlabel('Distance from outlet / Total channel length')
    ax.set_ylabel('Elevation normalised by outlet')

    plt.savefig(DataDirectory+fname_prefix+'_normalised_dist.png', dpi=300)

def PlotAllProfilesNormalisedElev(DataDirectory, fname_prefix):
    """
    Function to make a plot of all the normalised river profiles

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the river profile csv
    df = read_river_profile_csv(DataDirectory,fname_prefix)
    print ('Got the dataframe!')

    # normalise the distance by total length of channels upstream
    df['normalised_dist'] = df['distance_from_outlet']/df['total_length_upstream']

    # get the profile IDs
    river_ids = df['id'].unique()

    for id in river_ids:
        print ('This id is: ', id)
        this_df = df[df['id'] == id]
        outlet_elev = this_df['elevation'].min()
        this_df['normalised_elev'] = this_df['elevation']/outlet_elev
        ax.plot(this_df['distance_from_outlet'], this_df['normalised_elev'])

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation normalised by outlet')

    plt.savefig(DataDirectory+fname_prefix+'_normalised_elev.png', dpi=300)
    plt.clf()
    #plt.show()

def PlotAllProfiles(DataDirectory, fname_prefix):
    """
    Function to make a plot of all the river profiles

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # read in the river profile csv
    df = read_river_profile_csv(DataDirectory,fname_prefix)
    print ('Got the dataframe!')

    # get the profile IDs
    river_ids = df['id'].unique()

    for id in river_ids:
        print ('This id is: ', id)
        this_df = df[df['id'] == id]
        ax.plot(this_df['distance_from_outlet'], this_df['elevation'])

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+'_profiles.png', dpi=300)
    plt.clf()
    #plt.show()

def PlotProfilesAllTribuatires(DataDirectory,fname_prefix):
    """
    Function to make individual plots for each basin, with all the tributaries.
    Coloured by source.

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    for fname in glob(DataDirectory+"*_all_tribs*.csv"):

        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        df = pd.read_csv(fname)

        # get the profile IDs
        source_ids = df['source_id'].unique()
        basin_id = df['basin_id'][0]

        for id in source_ids:
            print ('This id is: ', id)
            this_df = df[df['source_id'] == id]
            print (this_df)
            ax.plot(this_df['distance_from_outlet'], this_df['elevation'])

        ax.set_xlabel('Distance from outlet (m)')
        ax.set_ylabel('Elevation (m)')

        plt.savefig(DataDirectory+fname_prefix+'_'+str(basin_id)+'_profiles.png', dpi=300)
        plt.clf()

def PlotProfilesAllSourcesSlope(DataDirectory,fname_prefix,slope_window_size=3):
    """
    Function to make a plot of all the channels coloured by source

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension
        slope_window_size (int): number of points on the channel from which to calculate slope

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    for fname in glob(DataDirectory+"*_all_sources*.csv"):

        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        df = pd.read_csv(fname)

        # get the profile IDs
        source_ids = df['source_id'].unique()

        for id in source_ids:
            print ('This id is: ', id)
            this_df = df[df['source_id'] == id]
            this_df = CalculateSlope(this_df, slope_window_size)
            this_df = this_df[this_df['slope'] != np.nan]
            ax.plot(this_df['distance_from_source'],this_df['slope'], lw=1)

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Gradient')

        plt.savefig(DataDirectory+fname_prefix+'_profiles_sources_slope.png', dpi=300)
        plt.clf()

def PlotProfilesAllSourcesElev(DataDirectory,fname_prefix,slope_window_size=3):
    """
    Function to make a plot of all the channels coloured by source

    Args:
        DataDirectory (str): the data directory
        fname_prefix (str): the filename of the DEM without extension
        slope_window_size (int): number of points on the channel from which to calculate slope

    Author: FJC
    """
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    for fname in glob(DataDirectory+"*_all_sources*.csv"):

        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        df = pd.read_csv(fname)

        # get the profile IDs
        source_ids = df['source_id'].unique()

        for id in source_ids:
            print ('This id is: ', id)
            this_df = df[df['source_id'] == id]
            source_elev = this_df['elevation'].max()
            this_df['normalised_elev'] = this_df['elevation']/source_elev
            ax.plot(this_df['distance_from_source'], this_df['elevation'], lw=1)

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Elevation (m)')

        plt.savefig(DataDirectory+fname_prefix+'_profiles_sources.png', dpi=300)
        plt.clf()

if __name__ == '__main__':

    DataDirectory = '/home/clubb/Data_for_papers/river_spaghetti/Pozo/cat1/'
    fname_prefix = 'Pozo_cat1_UTM11_WGS84_1m'
    # PlotAllProfiles(DataDirectory,fname_prefix)
    # PlotAllProfilesNormalisedElev(DataDirectory,fname_prefix)
    # PlotAllProfilesNormalised(DataDirectory,fname_prefix)
    #PlotProfilesAllTribuatires(DataDirectory,fname_prefix)
    slope_window_size = 25
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_sources1000.csv')
    df = CalculateSlope(df, slope_window_size)
    ClusterProfiles(df, profile_len = 1000, step=1, min_corr = 0.5)
    #PlotProfilesAllSourcesSlope(DataDirectory,fname_prefix, slope_window_size)
    #PlotProfilesAllSourcesElev(DataDirectory,fname_prefix)
