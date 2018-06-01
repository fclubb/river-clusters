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
from scipy import stats
from scipy.signal import correlate
from CorrCoef import Pearson
import math
import LSDPlottingTools as LSDP
from LSDMapFigure.PlottingRaster import MapFigure
import sys
from collections import defaultdict
import os

# Set up fonts for plots
label_size = 10
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size

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

def ProfilesRegularDistance(df, profile_len = 1000, step=2, slope_window_size=25):
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
    print("The profile length is {} m, any smaller channels will be removed.".format(profile_len))
    print("Assigning the profiles a common distance step of {} m".format(step))
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
        if not this_df.empty:
            slopes = this_df['slope'].as_matrix()[::-1] #need to reverse these as the distances need to be sorted
            distances = this_df['distance_from_outlet'].as_matrix()[::-1]
            if (distances.max() >= profile_len):
                profiles.append((distances, slopes))
                #thinned_df = thinned_df.append(this_df)
                final_sources.append(source)

    # now create the 2d array to store the data
    n_profiles = len(profiles)
    data = np.empty((n_profiles, len(reg_dist)))

    # create a new dataframe for storing the data about the selected profiles
    thinned_df = pd.DataFrame()

    # make a plot of the gradient vs. distance from source aligned to the outlet with
    # the resampled distance frame
    # set up a figure
    # fig = plt.figure(1, facecolor='white')
    # gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    # ax = fig.add_subplot(gs[5:100,10:95])
    #
    # loop through the profiles. For each point in the regularly spaced array,
    # find the index of the closest point in the distance array. Then use this to
    # assign the slope to the regularly spaced array
    for i, p in enumerate(profiles):
        reg_slope = []
        for d in reg_dist:
            idx = find_nearest_idx(p[0], d)
            reg_slope.append(p[1][idx])
            this_dist = p[0][idx]
            # get this distance and append the regular distance to the thinned df
            thinned_df = thinned_df.append(df.loc[(df['id']==final_sources[i]) & (df['distance_from_outlet'] == p[0][idx])])
        thinned_df.loc[(thinned_df['id'] == final_sources[i]), 'reg_dist'] = reg_dist
        data[i] = reg_slope

        # plot this profile
        # ax.plot(reg_dist[::-1], reg_slope, lw=1)

    # write the thinned_df to output in case we want to reload
    thinned_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.csv')
    #
    # # now save the figure
    # ax.set_xlabel('Distance from source (m)')
    # ax.set_ylabel('Gradient')
    # plt.savefig(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.png', dpi=300)
    # plt.clf()

    return thinned_df, data

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
    #print sources
    pts_in_profile = len(df[df['id'] == sources[0]])
    data = np.empty((len(sources), pts_in_profile))

    for i, src in enumerate(sources):
        this_df = df[df['id'] == src]
        data[i] = this_df['slope']
    #print data


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
    R = dendrogram(ln, color_threshold=thr, above_threshold_color=threshold_color)

    plt.axhline(y = thr, color = 'r', ls = '--')
    plt.savefig(DataDirectory+fname_prefix+"_upstream_dendrogram.png", dpi=300)
    plt.clf()

    for i,id in enumerate(source_ids):
        df.loc[df.id==id, 'cluster_id'] = cl[i]

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

    # make a plot of the gradient vs. distance from source aligned to the outlet
    fig_width_inches = 4.92126

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # get all unique source IDs from the dataframe
    ids = df['id'].unique()
    for id in ids:
        dist = df['distance_from_outlet'][df['id'] == id]
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
        ax.plot(dist[::-1], slopes, lw=1)

    # now save the figure
    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Gradient')
    plt.savefig(DataDirectory+fname_prefix+'_profiles_upstream.png', dpi=300)
    plt.clf()

    print("Got the slope over a window radius of {} m").format(slope_window_size)

    return df

def RemoveProfilesShorterThanThresholdLength(df, profile_len=1000):
    """
    Remove any profiles that are shorter than a threshold length.
    I can't believe I worked out how to do that in one line of code...!
    """
    # find the maximum distance for each id and remove any less than the profile length
    df = df[df.groupby('id')['distance_from_outlet'].transform('max') >= profile_len]
    return df

def RemoveNonUniqueProfiles(df):
    """
    From the regularly spaced distance dataframe, remove any profiles
    which are non-unique (e.g., they have the same nodes as another
    profile)
    """
    # group the dataframe by the source id, and then check if the
    # node column is repeated. If it is then remove the index from the dataframe
    if (len(df['id'].unique()) > 1):
        s = df.groupby('id')['node'].apply(tuple).drop_duplicates().index
        new_df = df.loc[df['id'].isin(s)]
    else:
        new_df = df

    return new_df
#---------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#---------------------------------------------------------------------#

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
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_reg_dist.csv')

    df = RemoveNonUniqueProfiles(df)

    # do the clustering
    cluster_df = ClusterProfiles(df, profile_len = profile_len, step=step, min_corr = min_corr, method = method)

    # find the unique clusters for plotting
    clusters = cluster_df['cluster_id'].unique()
    clusters.sort()

    counter = 0
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
                ax.plot(src_df['reg_dist'].as_matrix()[::-1], src_df['slope'].as_matrix(), lw=1, color=colors[counter])
                # save the colour to the cluster dataframe for later plots
                cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = colors[counter]
            counter +=1
        else:
            ax.plot(this_df['reg_dist'].as_matrix()[::-1], this_df['slope'].as_matrix(), lw=1, color=threshold_color)
            # save the colour to the cluster dataframe for later plots
            cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = threshold_color

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Gradient')
        ax.set_title('Cluster {}'.format(int(cl)))

        plt.savefig(DataDirectory+fname_prefix+('_profiles_upstream_clustered_{}.png').format(int(cl)), dpi=300)
        plt.clf()

    # write the clustered dataframe to csv
    cluster_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')

    return cluster_df

def MakeHillshadePlotClusters():
    """
    Make a shaded relief plot of the raster with the channels coloured by the cluster
    value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

    Args:
        cluster_df: dataframe with the clustered information

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
    cluster_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')


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
        # plot the whole channel network in black
        ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.8',manual_size=2, zorder=1)
        # plot the clustered profiles in the correct colour
        this_df = cluster_df[cluster_df.cluster_id == cl]
        this_colour = str(this_df.colour.unique()[0])
        ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_hs_clusters.png', FigFormat='png', Fig_dpi = 300) # Save the figure

def PlotMedianProfiles():
    """
    Make a summary plot showing the median profile for each cluster, both in
    gradient-distance and elevation-distance space.

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')

    # find out some info
    clusters = df.cluster_id.unique()
    clusters.sort()
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
        ax[i].plot(dist_array[::-1],median_gradients,color=this_colour, lw=1)
        ax[i].fill_between(dist_array[::-1], lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)
        ax[i].text(0.9, 0.8,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes,fontsize=8)

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
        ax.plot(dist_array[::-1],median_elevs,color=this_colour, lw=1)
        ax.fill_between(dist_array[::-1], lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)

    ax.set_xlabel('Distance from source (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+('_profiles_median_elev.png'), dpi=300)
    plt.clf()

def PlotUniqueStreamsWithLength(step=2, slope_window_size=25):
    """
    Function to make a plot of the number of unique channels you get (at least one non-overlapping
    node with different profile lengths that you analyse).

    Author: FJC
    """
    # read in the original csv
    df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_tribs.csv')

    max_length = df['distance_from_outlet'].max()
    len_step=100
    profile_lengths = np.arange(len_step, (max_length/100).astype(int)*100, step=len_step)

    # calculate the slope
    df = CalculateSlope(df, slope_window_size)

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    # for each length, remove any sources that are below the desired length.
    for len in profile_lengths:
        thinned_df = RemoveProfilesShorterThanThresholdLength(df, len)
        reg_df, __ = ProfilesRegularDistance(thinned_df, profile_len = len, step=step, slope_window_size=slope_window_size)
        reg_df = RemoveNonUniqueProfiles(reg_df)
        n_sources = reg_df['id'].unique().size
        print n_sources
        ax.scatter(len, n_sources, s = 25, c='w', edgecolors='k')

    ax.set_xlabel('Profile length (m)')
    ax.set_ylabel('Number of unique channels')

    plt.savefig(DataDirectory+fname_prefix+'_n_channels_with_length.png', dpi=300)
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

    # set colour palette: 8 class Set 1 from http://colorbrewer2.org
    N_colors = 8
    colors = LSDP.colours.list_of_hex_colours(N_colors, 'Dark2')
    threshold_color = '#377eb8'

    # read in the original csv
    df = pd.read_csv(DataDirectory+args.fname_prefix+'_all_tribs.csv')

    # calculate the slope
    df = RemoveProfilesShorterThanThresholdLength(df, args.profile_len)
    df = CalculateSlope(df, args.slope_window)

    # # shift the profiles to a common distance frame
    regular_df, data = ProfilesRegularDistance(df, profile_len = args.profile_len, step=args.step, slope_window_size=args.slope_window)


    # now do the clustering
    cluster_df = PlotProfilesByCluster(slope_window_size=args.slope_window,profile_len=args.profile_len,step=args.step,method=args.method,min_corr=args.min_corr)

    PlotMedianProfiles()
    MakeHillshadePlotClusters()
    #PlotUniqueStreamsWithLength()
