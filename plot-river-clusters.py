#---------------------------------------------------------------------#
# Clustering of river profiles
# Developed by Fiona Clubb
#              Bodo Bookhagen
#              Aljoscha Rheinwalt
# University of Potsdam
#---------------------------------------------------------------------#

# setting backend to run on server
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import rcParams
import matplotlib.cm as cm
from matplotlib.colors import LinearSegmentedColormap
from glob import glob
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster, set_link_color_palette
import scipy.stats as ss
from scipy.ndimage.interpolation import shift
from CorrCoef import Pearson
import math
import LSDPlottingTools as LSDP
from LSDMapFigure.PlottingRaster import MapFigure
import sys
from collections import defaultdict
import os
#import seaborn as sns

#=============================================================================
# This is just a welcome screen that is displayed if no arguments are provided.
#=============================================================================
def print_welcome():

    print("\n\n=======================================================================")
    print("Hello! I'm going to do some river profile clustering for you.")
    print("You will need to tell me which directory to look in.")
    print("Use the -dir flag to define the working directory.")
    print("If you don't do this I will assume the data is in the same directory as this script.")
    print("You also need the -fname flag which will give the prefix of the raster files.")
    print("For help type:")
    print("   python plot_river_clusters.py -h\n")
    print("=======================================================================\n\n ")

#---------------------------------------------------------------------#
# ANALYSIS FUNCTIONS
#---------------------------------------------------------------------#
def find_nearest_idx(array,value):
    """
    Given a value, find the index of the point in the array which is closest
    to that value.
    Author: FJC
    """
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return idx-1
    else:
        return idx

def ProfilesRegularDistance(df, profile_len = 100, step=2, slope_window_size=25):
    """
    This function takes the dataframe of river profiles and creates an array
    that can be used for the time series clustering.  For each profile, the slope
    is assigned to a common distance step, in metres (default=1m)

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): step size that you want in metres. This should be greater than the maximum
        possible spacing between your points (max_spacing = sqrt(2 * DataRes^2)). Default = 2
        slope_window_size (int): window over which slope was calculated

    Returns: array of size (n_profiles, profile_len) containing the slopes along each profile

    Author: FJC
    """
    print("Resampling the profiles to a common distance step of {} m".format(step))
    # find the slope radius. We don't calculate slope for the first or last few nodes, which
    # are below the min required window size. Therefore reg dist needs to be smaller
    slope_radius = ((slope_window_size-1)/2)

    # create new array of regularly spaced differences
    reg_dist = np.arange(slope_radius+step, profile_len-(slope_radius), step)

    # find the minimum length that the array can be (profile length/root2)
    min_length = profile_len/(math.sqrt(2))

    # loop through the dataframe and store the data for each profile as an array of
    # slopes and distances
    profiles = []
    source_ids = df['id'].unique()
    final_sources = []
    for i, source in enumerate(source_ids):
        this_df = df[df['id'] == source]
        this_df = this_df[np.isnan(this_df['slope']) == False]  # remove nans
        slopes = this_df['slope'].as_matrix()
        distances = this_df['distance_from_source'].as_matrix()
        if (len(slopes) >= min_length):
            profiles.append((distances, slopes))
            #thinned_df = thinned_df.append(this_df)
            final_sources.append(source)

    # now create the 2d array to store the data
    n_profiles = len(profiles)
    data = np.empty((n_profiles, len(reg_dist)))

    # create a new dataframe for storing the data about the selected profiles
    thinned_df = pd.DataFrame()

    # loop through the profiles. For each point in the regularly spaced array,
    # find the index of the closest point in the distance array. Then use this to
    # assign the slope to the regularly spaced array
    for i, p in enumerate(profiles):
        reg_slope = []
        for d in reg_dist:
            idx = find_nearest_idx(p[0], d)
            reg_slope.append(p[1][idx])
            # get this distance and append the regular distance to the thinned df
            thinned_df = thinned_df.append(df.loc[(df['id']==final_sources[i]) & (df['distance_from_source'] == p[0][idx])])
        thinned_df.loc[(thinned_df['id'] == final_sources[i]), 'reg_dist'] = reg_dist
        data[i] = reg_slope

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_reg_dist.csv')


    return thinned_df, data

def CrossCorrelateProfiles(y1, y2):
    """
    Cross correlate two profiles to calculate the lag between them, and shift
    y2 to remove the lag.

    Args:
        y1: numpy array of profile 1
        y2: numpy array of profile 2

    Returns:
        shifted_y2: array of the shifted profile.

    Author: FJC
    """
    npts = len(y1)
    # make an array with the lags
    lags = np.arange(-npts + 1, npts)
    ccov = np.correlate(y1 - y1.mean(), y2 - y2.mean(), mode='full')
    ccor = ccov / (npts * y1.std() * y2.std())

    # find the maximum lag
    maxlag = lags[np.argmax(ccor)]
    print ("max correlation at lag {}".format(maxlag))

    # shift y2 so to reduce the lag
    shifted_y2 = shift(y2, maxlag, cval=np.nan)

    return shifted_y2

def ShiftProfiles(df, shift_steps=100):
    """
    Shift the profiles by a certain number of steps backwards and forwards compared
    to a reference profile. Check the correlation of each, and take the max correlation.

        Args:
        df: the dataframe of the profiles with the regular distances
        shift_steps: the number of steps to check for the correlation. Default = 5 steps.

    Author: FJC
    """
    print("Shifting the profiles based on their maximum correlation...")
    # get the source ids of each profile
    sources = df['id'].unique()
    pts_in_profile = len(df[df['id'] == sources[0]])
    data = np.empty((len(sources), pts_in_profile))

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]

        data[i] = this_df['slope']

    y1 = data[0]
    for idx, y in enumerate(data):
        src_df = df[df.id == sources[idx]]
        #shift the profile
        shifted_y2 = CrossCorrelateProfiles(y1, y)
        # append to the dataframe
        df.loc[df.id == sources[idx], 'slope'] = shifted_y2
        plt.plot(src_df['distance_from_source'], shifted_y2, lw=1)


    plt.xlabel('Distance from source (m)')
    plt.ylabel('Gradient (m/m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_shifted.png'), dpi=300)
    plt.clf()

    df.to_csv(DataDirectory+fname_prefix+'_profiles_shifted.csv')

def ClusterProfiles(df, profile_len=100, step=2, min_corr=0.5, method='complete'):
    """
    Cluster the profiles based on gradient and distance from source.
    Aggolmerative clustering, see here for more info:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html#scipy.cluster.hierarchy.linkage

    Args:
        df: pandas dataframe from the river profile csv.
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): the spacing in metres between the data points that you want, default = 1m
        min_corr (float): minimum correlation threshold for clustering
        method (str): clustering method to use, see scipy docs. Can be 'single', 'complete', 'average',
        'weighted', 'centroid', 'median', or 'ward'. Default is 'complete'.

    Author: AR, FJC
    """
    print ("Now I'm going to do some hierarchical clustering...")

    # get the data from the dataframe into the right format for clustering
    sources = df['id'].unique()
    pts_in_profile = len(df[df['id'] == sources[0]])
    data = np.empty((len(sources), pts_in_profile))

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]

        data[i] = this_df['slope']

    # we could have a look at the ranks too ..
    # correlations
    cc = Pearson(data)

    # distances
    dd = np.arccos(cc)

    # do agglomerative clustering by stepwise pair matching
    # based on angle between scalar products of time series
    ln = linkage(dd, method=method)

    # define threshold for cluster determination
    thr = np.arccos(min_corr)

    # compute cluster indices
    cl = fcluster(ln, thr, criterion = 'distance')
    print("I've finished! I found {} clusters for you :)".format(cl.max()))
    print cl

    set_link_color_palette(colors)

    source_ids = df['id'].unique()

    plt.title('Hierarchical Clustering Dendrogram')
    plt.xlabel('sample index')
    plt.ylabel('distance')
    R = dendrogram(ln, color_threshold=1, above_threshold_color='k')

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_dendrogram.png", dpi=300)
    plt.clf()

    for i,id in enumerate(source_ids):
        df.loc[df.id==id, 'cluster_id'] = cl[i]
        # write the colour code for this cluster ID

    return df

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
    ids = df['id'].unique()
    for id in ids:
        dist = df['distance_from_source'][df['id'] == id]
        elev = df['elevation'][df['id'] == id]
        slopes = np.empty(len(dist))

        pts_array = np.column_stack((dist,elev))

        slicer = (slope_window_size - 1)/2

        for index, x in enumerate(pts_array):
            # find the rows above and below relating to the window size
            this_slice = pts_array[index-slicer:index+slicer+1]
            if len(this_slice) == slope_window_size:
                # now regress this slice
                x = this_slice[:,0]
                y = this_slice[:,1]
                #print x, y
                slope, intercept, r_value, p_value, std_err = stats.linregress(x, y)
                slopes[index] = abs(slope)
                #print slope
            else:
                slopes[index] = np.nan

        df.loc[df.id==id, 'slope'] = slopes

    print("Got the slope over a window radius of {} m").format(slope_window_size)

    return df

#---------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#---------------------------------------------------------------------#

def PlotProfilesAllSourcesElev(slope_window_size=3,profile_len=100, step=2, min_corr=0.5):
    """
    Function to make a plot of all the channels coloured by source

    Args:
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
        source_ids = df['id'].unique()

        for id in source_ids:
            print ('This id is: ', id)
            this_df = df[df['id'] == id]
            source_elev = this_df['elevation'].max()
            this_df['normalised_elev'] = this_df['elevation']/source_elev
            ax.plot(this_df['distance_from_source'], this_df['elevation'], lw=1)

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Elevation (m)')

        plt.savefig(DataDirectory+fname_prefix+'_profiles_sources.png', dpi=300)
        plt.clf()

def PlotProfilesByCluster(slope_window_size=3,profile_len=100, step=2, min_corr=0.5, method = 'complete'):
    """
    Function to make plots of the river profiles in each cluster

    Args:
        slope_window_size: window size in metres over which to calculate slope
        profile_len (int): number of data points in each profile (= distance from source)
        step (int): the spacing in metres between the data points that you want, default = 1m
        min_corr (float): minimum correlation threshold for clustering
        method (str): method for performing the clustering

    Author: FJC
    """
    # read in the csv
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_reg_dist.csv')

    # do the clustering
    cluster_df = ClusterProfiles(df, profile_len = profile_len, step=step, min_corr = min_corr, method = method)

    # find the unique clusters for plotting
    clusters = cluster_df['cluster_id'].unique()

    for cl in clusters:
        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        this_df = cluster_df[cluster_df['cluster_id'] == cl]
        cl = int(this_df.iloc[0]['cluster_id'])
        sources = this_df['id'].unique()
        if (len(sources) > 1):
            for idx, src in enumerate(sources):
                src_df = this_df[this_df['id'] == src]
                src_df = src_df[src_df['slope'] != np.nan]
                ax.plot(src_df['reg_dist'], src_df['slope'], lw=1, color=colors[cl-1])
                # save the colour to the cluster dataframe for later plots
                cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = colors[cl-1]
        else:
            ax.plot(this_df['reg_dist'], this_df['slope'], lw=1, color='k')
            # save the colour to the cluster dataframe for later plots
            cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = 'k'

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Gradient')
        ax.set_title('Cluster {}'.format(int(cl)))

        plt.savefig(DataDirectory+fname_prefix+('_profiles_clustered_{}.png').format(int(cl)), dpi=300)
        plt.clf()

    # write the clustered dataframe to csv
    cluster_df.to_csv(DataDirectory+fname_prefix+'_profiles_clustered.csv')

    return cluster_df

def PlotProfileShifting():
    """
    Make a plot of the profiles before and after they were shifted
    for checking

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_reg_dist.csv')
    shifted_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_shifted.csv')

    sources = df['id'].unique()

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    for src in sources:
        this_df = df[df['id'] == src]
        ax.plot(this_df['reg_dist'], this_df['slope'], lw=0.5)

    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Gradient (m/m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_reg_dist.png'), dpi=300)
    plt.clf()

    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    for src in sources:
        this_shifted_df = shifted_df[shifted_df['id'] == src]
        ax.plot(this_shifted_df['reg_dist'], this_shifted_df['slope'], lw=0.5)

    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Gradient (m/m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_shifted.png'), dpi=300)
    plt.clf()


def MakeHillshadePlotClusters():
    """
    Make a shaded relief plot of the raster with the channels coloured by the cluster
    value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

    Args:
        cluster_df: dataframe with the clustered information

    Author: FJC
    """
    cluster_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_clustered.csv')
    # Set up fonts for plots
    label_size = 10
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = ['arial']
    rcParams['font.size'] = label_size

    # set figure sizes based on format
    fig_width_inches = 4.92126

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HillshadeName = fname_prefix+'_hs'+raster_ext

    # create the map figure
    MF = MapFigure(HillshadeName, DataDirectory,coord_type="UTM_km")
    clusters = cluster_df.cluster_id.unique()
    for cl in clusters:
        this_df = cluster_df[cluster_df.cluster_id == cl]
        this_colour = str(this_df.colour.unique()[0])
        ChannelPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_hs_clusters.png', FigFormat='png', Fig_dpi = 300) # Save the figure

def PlotMedianProfiles():
    """
    Make a summary plot showing the median profile for each cluster, both in
    gradient-distance and elevation-distance space.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_clustered.csv')

    # find out some info
    clusters = df.cluster_id.unique()
    sources = df.id.unique()
    dist_array = df[df.id == sources[0]].reg_dist.as_matrix()

    # set up a figure
    fig,ax = plt.subplots(nrows=len(clusters),ncols=1, figsize=(5,6), sharex=True, sharey=True)
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # for each cluster, get the mean gradient for each regular distance
    for i, cl in enumerate(clusters):

        cluster_df = df[df.cluster_id == cl]
        median_gradients = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.median() for x in dist_array])
        lower_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.quantile(0.25) for x in dist_array])
        upper_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].slope.quantile(0.75) for x in dist_array])
        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax[i].plot(dist_array,median_gradients,color=this_colour, lw=1)
        ax[i].fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)

    # set axis labels
    plt.xlabel('Distance from source (m)')
    plt.ylabel('Gradient')

    # save and clear the figure
    plt.savefig(DataDirectory+fname_prefix+('_profiles_median.png'), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # for each cluster, get the mean gradient for each regular distance
    for cl in clusters:
        cluster_df = df[df.cluster_id == cl]
        median_elevs = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.median() for x in dist_array])
        lower_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.25) for x in dist_array])
        upper_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.75) for x in dist_array])
        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax.plot(dist_array,median_elevs,color=this_colour, lw=1)
        ax.fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)

    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_median_elev.png'), dpi=300)
    plt.clf()

def MakeSlopeAreaPlots():
    """
    Function to make a plot of gradient against drainage area for each cluster.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_clustered.csv')

    # find out some info
    clusters = df.cluster_id.unique()
    sources = df.id.unique()

    # set up a figure
    fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # for each cluster, get the median gradient for a SA plot
    for cl in clusters:
        cluster_df = df[df.cluster_id == cl]
        # get the drainage area bins
        max_area = cluster_df.drainage_area.max()
        median_gradients = [cluster_df[cluster_df.drainage_area == x].slope.median() for x in drainage_area]
        # get the colour from the dataframe
        this_colour = str(cluster_df.colour.unique()[0])
        ax.scatter(drainage_area,median_gradients,color=this_colour, lw=1)

    ax.set_xlabel('Drainage area (m$^2$)')
    ax.set_ylabel('Gradient (m/m)')

    plt.savefig(DataDirectory+fname_prefix+('_slope_area_median.png'), dpi=300)
    plt.clf()


    # # for each cluster, get the mean gradient for each regular distance
    # for cl in clusters:
    #     # set up a figure
    #     fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    #     gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    #     ax = fig.add_subplot(gs[5:100,10:95])
    #
    #     cluster_df = df[df.cluster_id == cl]
    #     median_elevs = [cluster_df[cluster_df.reg_dist == x].elevation.median() for x in dist_array]
    #     # get the colour from the dataframe
    #     this_colour = str(cluster_df.colour.unique()[0])
    #     ax.plot(dist_array,median_elevs,color=this_colour, lw=1)
    #
    # ax.set_xlabel('Distance from source (m)')
    # ax.set_ylabel('Elevation (m)')
    #
    # plt.savefig(DataDirectory+fname_prefix+('_profiles_median_elev.png'), dpi=300)
    # plt.clf()

def PlotUniqueChannelsWithLengthScale(profile_len=1000):
    """
    This function takes the csv file with all the sources, and makes a plot
    showing the number of unique channels with each length scale.

    Args:
        profile_len: max length of the profile in metres.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_sources{}.csv'.format(profile_len)))

    # make an array of the profile lengths you want to test. The one in the csv has to be the max
    # one. We'll test 100 m spacings from this.
    step=100
    lengths = np.arange(0+step, profile_len+step, step=step)
    print lengths




if __name__ == '__main__':

    # If there are no arguments, send to the welcome screen
    if not len(sys.argv) > 1:
        full_paramfile = print_welcome()
        sys.exit()

    # Get the arguments
    import argparse
    parser = argparse.ArgumentParser()

    # The location of the data files
    parser.add_argument("-dir", "--base_directory", type=str, help="The base directory with the m/n analysis. If this isn't defined I'll assume it's the same as the current directory.")
    parser.add_argument("-fname", "--fname_prefix", type=str, help="The prefix of your DEM WITHOUT EXTENSION!!! This must be supplied or you will get an error (unless you're running the parallel plotting).")

    # The options for clustering
    parser.add_argument("-len", "--profile_len", type=int, help="The length of the profiles, you should have specified this in the parameter file for the spaghetti code. Default is 1000 m.", default=1000)
    parser.add_argument("-sw", "--slope_window", type=int, help="The window size for calculating the slope based on a regression through an equal number of nodes upstream and downstream of the node of interest. This is the total number of nodes that are used for calculating the slope. For example, a slope window of 25 would fit a regression through 12 nodes upstream and downstream of the node, plus the node itself. The default is 25 nodes.", default=25)
    parser.add_argument("-m", "--method", type=str, help="The method for clustering, see the scipy linkage docs for more information. The default is 'complete'.", default='complete')
    parser.add_argument("-c", "--min_corr", type=float, help="The minimum correlation for defining the clusters. Use a smaller number to get less clusters, and a bigger number to get more clusters (from 0 = no correlation, to 1 = perfect correlation). The default is 0.5.", default=0.5)
    parser.add_argument("-step", "-step", type=int, help="The regular spacing in metres that you want the profiles to have for the clustering. This should be greater than sqrt(2* DataRes^2).  The default is 2 m.", default = 2)
    parser.add_argument("-shift", "--shift_steps", type=int, help="The number of steps that you want to shift the profiles by to prior to clustering. The default is 50. You can get the shift in metres by multiplying the shift steps by the spacing for the regularly spaced profiles (e.g. 50 with a step of 2 will be 100 m for testing the shift). Default = 50 steps", default = 50)

    args = parser.parse_args()

    if not args.fname_prefix:
        print("WARNING! You haven't supplied your DEM name. Please specify this with the flag '-fname'")
        sys.exit()
    else:
        fname_prefix = args.fname_prefix

    # get the base directory
    if args.base_directory:
        DataDirectory = args.base_directory
        # check if you remembered a / at the end of your path_name
        if not DataDirectory.endswith("/"):
            print("You forgot the '/' at the end of the directory, appending...")
            DataDirectory = DataDirectory+"/"
    else:
        print("WARNING! You haven't supplied the data directory. I'm using the current working directory.")
        DataDirectory = os.getcwd()

    # set colour palette: 6 class Dark 2 from http://colorbrewer2.org
    colors = ['#1b9e77', '#d95f02', '#7570b3', '#e7298a', '#66a61e', '#e6ab02']
    cmap_name = 'Dark2'
    dark2 = LinearSegmentedColormap.from_list(cmap_name, colors, N=len(colors))

    # check to see if you have ran the analyses before
    regular_csv = DataDirectory+args.fname_prefix+'_profiles_reg_dist.csv'
    shifted_csv = DataDirectory+args.fname_prefix+'_profiles_shifted.csv'
    clustered_csv = DataDirectory+args.fname_prefix+'_profiles_clustered.csv'
    if not os.path.isfile(clustered_csv):
        if not os.path.isfile(shifted_csv):
            if not os.path.isfile(regular_csv):
                # read in the original csv
                df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_sources{}.csv'.format(args.profile_len))
                # calculate the slope
                df = CalculateSlope(df, args.slope_window)
                # shift the profiles to a common distance frame
                regular_df, data = ProfilesRegularDistance(df, profile_len = args.profile_len, step=args.step, slope_window_size=args.slope_window)
            # shift the profiles to reduce lag
            regular_df = pd.read_csv(regular_csv)
            ShiftProfiles(regular_df, shift_steps=args.shift_steps)
            #PlotProfileShifting()

        # now do the clustering
        #cluster_df = PlotProfilesByCluster(slope_window_size=args.slope_window,profile_len=args.profile_len,step=args.step,method=args.method,min_corr=args.min_corr)

    #PlotMedianProfiles()
    #MakeHillshadePlotClusters()
