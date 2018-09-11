import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from matplotlib import rcParams
from scipy import stats

# Set up fonts for plots
label_size = 10
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['arial']
rcParams['font.size'] = label_size

def list_of_hex_colours(N, base_cmap):
    """
    Return a list of colors from a colourmap as hex codes

        Arguments:
            cmap: colormap instance, eg. cm.jet.
            N: number of colors.

        Author: FJC
    """
    cmap = cm.get_cmap(base_cmap, N)

    hex_codes = []
    for i in range(cmap.N):
        rgb = cmap(i)[:3] # will return rgba, we take only first 3 so we get rgb
        hex_codes.append(mcolors.rgb2hex(rgb))
    return hex_codes

def lower_p(x):
    """
    Calculate IQR of x
    """
    lower_p = np.percentile(x, 25)
    return lower_p

def upper_p(x):
    """
    Calculate IQR of x
    """
    upper_p = np.percentile(x, 75)
    return upper_p

def bin_slope_area_data(slope, area, nbins=20):
    """
    Perform log binning on the slope-area data
    Args:
        slope: array of slope data
        area: array of area data
    Returns: arrays of LOG median slopes, interquartile range, bin centres (to plot area), bin edges
    """
    # log the data
    log_slope = np.log10(slope)
    log_area = np.log10(area)

    # do log binning
    bin_meds, bin_edges, binnumber = stats.binned_statistic(log_area, log_slope, bins=nbins, statistic='median')
    bin_width = (bin_edges[1] - bin_edges[0])
    bin_centres = bin_edges[1:] - bin_width/2

    # get the IQR
    lower_per, bin_edges, binnumber = stats.binned_statistic(log_area, log_slope, bins=nbins, statistic=lower_p)
    upper_per, bin_edges, binnumber = stats.binned_statistic(log_area, log_slope, bins=nbins, statistic=upper_p)

    return bin_meds, lower_per, upper_per, bin_centres, bin_edges
#---------------------------------------------------------------------#
# PLOTTING FUNCTIONS
#---------------------------------------------------------------------#
def PlotProfilesByCluster(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Function to make plots of the river profiles in each cluster

    Author: FJC
    """

    # find the maximum distance from outlet in each cluster and use this to sort
    # the data
    # idx = cluster_df.groupby(['cluster_id'])['distance_from_outlet'].transform(max) == df['distance_from_outlet']
    # lengths = cluster_df[idx]
    # lengths = lengths.sort_values(by=['distance_from_outlet'])
    # print lengths
    # # find the unique clusters for plotting
    # clusters = lengths['cluster_id'].tolist()
    cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))
    clusters = cluster_df['cluster_id'].unique()
    #clusters.sort()


    #counter = 0
    for cl in clusters:
        # set up a figure
        fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
        gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
        ax = fig.add_subplot(gs[5:100,10:95])

        this_df = cluster_df[cluster_df['cluster_id'] == cl]
        cl = int(this_df.iloc[0]['cluster_id'])
        this_colour = str(this_df.colour.unique()[0])
        sources = this_df['id'].unique()
        if (len(sources) > 1):
            for idx, src in enumerate(sources):
                src_df = this_df[this_df['id'] == src]
                src_df = src_df[src_df['slope'] != np.nan]
                ax.plot(src_df['reg_dist'].values, src_df['slope'].values, lw=1, color=this_colour)
                # save the colour to the cluster dataframe for later plots
                #cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = colors[counter]
            #counter +=1
        # else:
        #     ax.plot(this_df['distance_from_outlet'].values, this_df['elevation'].values, lw=1, color=threshold_color)
            # save the colour to the cluster dataframe for later plots
            #cluster_df.loc[cluster_df.cluster_id==cl, 'colour'] = threshold_color

        ax.set_xlabel('Distance from source (m)')
        ax.set_ylabel('Gradient')
        ax.set_title('Cluster {}'.format(int(cl)))

        plt.savefig(OutDirectory+fname_prefix+('_profiles_SO{}_CL{}.png').format(stream_order, int(cl)), dpi=300)
        plt.clf()

    # write the clustered dataframe to csv
    #cluster_df.to_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv', index=False)

    return cluster_df

def PlotMedianProfiles(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Make a summary plot showing the median profile for each cluster, both in
    gradient-distance and elevation-distance space.

    Author: FJC
    """
    df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))

    # find out some info
    clusters = df.cluster_id.unique()
    clusters.sort()
    sources = df.id.unique()
    max_source = df.loc[df['reg_dist'].idxmax()]['id']
    dist_array = df[df.id == max_source].reg_dist.values

    # set up a figure
    fig,ax = plt.subplots(nrows=1,ncols=1, figsize=(6,4), sharex=True, sharey=True)
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
        ax.grid(color='0.8', linestyle='--', which='both')
        ax.plot(dist_array,median_gradients,color=this_colour, lw=1)
        ax.fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)
        #ax.set_ylim(0,0.4)
        #ax.text(0.9, 0.8,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax.transAxes,fontsize=8)

    # set axis labels
    plt.xlabel('Distance from source (m)', labelpad=0.1, fontsize=14)
    plt.ylabel('Gradient', labelpad=10, fontsize=14)
    plt.subplots_adjust(bottom=0.2, left=0.15)

    # save and clear the figure
    plt.savefig(OutDirectory+fname_prefix+('_profiles_median_SO{}.png'.format(stream_order)), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

    # # set up a figure
    # fig = plt.figure(1, facecolor='white',figsize=(4.92126,3.2))
    # gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    # ax = fig.add_subplot(gs[5:100,10:95])
    #
    # # for each cluster, get the mean gradient for each regular distance
    # for cl in clusters:
    #     cluster_df = df[df.cluster_id == cl]
    #     median_elevs = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.median() for x in dist_array])
    #     lower_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.25) for x in dist_array])
    #     upper_quantile = np.asarray([cluster_df[cluster_df.reg_dist == x].elevation.quantile(0.75) for x in dist_array])
    #     # get the colour from the dataframe
    #     this_colour = str(cluster_df.colour.unique()[0])
    #     ax.grid(color='0.8', linestyle='--', which='both')
    #     ax.plot(dist_array,median_elevs,color=this_colour, lw=1)
    #     ax.fill_between(dist_array, lower_quantile, upper_quantile, facecolor=this_colour, alpha=0.2)
    #
    # ax.set_xlabel('Distance from outlet (m)')
    # ax.set_ylabel('Elevation (m)')
    #
    # plt.savefig(OutDirectory+fname_prefix+('_profiles_median_elev_SO{}.png'.format(stream_order)), dpi=300)
    # plt.clf()

def PlotSlopeAreaAllProfiles(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Make a summary plot showing the S-A plot for each cluster. Includes the data from
    all the way down the channel profile

    Author: FJC
    """
    cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))
    df = pd.read_csv(DataDirectory+fname_prefix+'_slopes.csv')

    # find out some info
    clusters = cluster_df.cluster_id.unique()
    clusters.sort()
    sources = cluster_df.id.unique()

    # set up a figure
    fig,ax = plt.subplots(nrows=len(clusters),ncols=1, figsize=(5,8), sharex=False, sharey=False)
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # for each cluster, get the mean gradient for each regular distance
    for i, cl in enumerate(clusters):

        # find out which sources are in the full df since we want to plot the full profiles for each cluster.
        this_df = cluster_df[cluster_df.cluster_id == cl]
        these_sources = this_df.id.unique()

        all_df =  df[df['id'].isin(these_sources)]

        # calculate the channel steepness
        area = all_df['drainage_area'].values
        med_slopes, lower_per, upper_per, bin_centres, _ = bin_slope_area_data(all_df['slope'], area)
        #print(med_slopes)
        # print(iqr_slopes)

        #print(log_slope)
        gradient, intercept, r, p, std = stats.linregress(bin_centres, med_slopes)
        print(intercept)
        intercept = float(10**intercept)
        print("Steepness index: {}".format(intercept))
        print("concavity: {}".format(gradient))
        x2 = np.linspace(1,area.max(),100)
        y2 = intercept*x2**(gradient)

        # transform binned data into normal for plotting
        med_slopes = np.array([10**x for x in med_slopes])
        med_areas = np.array([10**x for x in bin_centres])
        lower_per = np.array([10**x for x in lower_per])
        upper_per = np.array([10**x for x in upper_per])
        upper_err = upper_per - med_slopes
        lower_err = med_slopes - lower_per

        # get the colour from the dataframe
        this_colour = str(this_df.colour.unique()[0])
        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].scatter(area, all_df['slope'].values, color=this_colour, s=1)
        ax[i].errorbar(med_areas, med_slopes, xerr=None, yerr=[lower_err, upper_err], fmt='o', ms=5, marker='D', mfc='w', mec='k', zorder=3, c='k')
        # ax[i].scatter(med_areas, med_slopes, color='w',zorder=3, s=20, marker='D', edgecolors='k')
        ax[i].plot(x2, y2, "--", c='k')
        ax[i].text(0.15, 0.1,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes,fontsize=12)
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        #ax[i].set_xlim(900,10000)
        ax[i].set_ylim(0.0001, 10)
        ax[i].set_title('$k_s$ = {}; $\\theta$ = {}'.format(round(intercept,4), round(abs(gradient),2)), fontsize=16)


    # set axis labels
    plt.xlabel('Drainage area (m$^2$)', fontsize=14)
    plt.ylabel('Gradient (m/m)', labelpad=15, fontsize=14)
    plt.subplots_adjust(left=0.15, hspace=0.3)

    # save and clear the figure
    plt.savefig(OutDirectory+fname_prefix+('_SA_median_SO{}.png'.format(stream_order)), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

def PlotSlopeArea(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Make a summary plot showing the S-A plot for each cluster. Only the data from
    the profiles that were included in the clustering.

    Author: FJC
    """
    cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))

    # find out some info
    clusters = cluster_df.cluster_id.unique()
    clusters.sort()
    sources = cluster_df.id.unique()

    # set up a figure
    fig,ax = plt.subplots(nrows=len(clusters),ncols=1, figsize=(5,8), sharex=False, sharey=False)
    # make a big subplot to allow sharing of axis labels
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # for each cluster, get the mean gradient for each regular distance
    for i, cl in enumerate(clusters):

        # find out which sources are in the full df since we want to plot the full profiles for each cluster.
        this_df = cluster_df[cluster_df.cluster_id == cl]

        # calculate the channel steepness
        area = this_df['drainage_area'].values
        nbins=10
        med_slopes, lower_per, upper_per, bin_centres, _ = bin_slope_area_data(this_df['slope'], area, nbins)
        # print(iqr_slopes)

        #print(log_slope)
        gradient, intercept, r, p, std = stats.linregress(bin_centres, med_slopes)
        print(intercept)
        intercept = float(10**intercept)
        print("Steepness index: {}".format(intercept))
        print("concavity: {}".format(gradient))
        x2 = np.linspace(1,8000,100)
        y2 = intercept*x2**(gradient)

        # transform binned data into normal for plotting
        med_slopes = np.array([10**x for x in med_slopes])
        med_areas = np.array([10**x for x in bin_centres])
        lower_per = np.array([10**x for x in lower_per])
        upper_per = np.array([10**x for x in upper_per])
        upper_err = upper_per - med_slopes
        lower_err = med_slopes - lower_per

        # get the colour from the dataframe
        this_colour = str(this_df.colour.unique()[0])
        ax[i].grid(color='0.8', linestyle='--', which='both')
        ax[i].scatter(area, this_df['slope'].values, color=this_colour, s=1)
        ax[i].errorbar(med_areas, med_slopes, xerr=None, yerr=[lower_err, upper_err], fmt='o', ms=5, marker='D', mfc='w', mec='k', zorder=3, c='k')
        # ax[i].scatter(med_areas, med_slopes, color='w',zorder=3, s=20, marker='D', edgecolors='k')
        ax[i].plot(x2, y2, "--", c='k')
        ax[i].text(0.15, 0.1,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[i].transAxes,fontsize=12)
        ax[i].set_xscale('log')
        ax[i].set_yscale('log')
        ax[i].set_xlim(900,7000)
        ax[i].set_ylim(0.001, 1)
        ax[i].set_title('$k_s$ = {}; $\\theta$ = {}'.format(round(intercept,4), round(abs(gradient),2)), fontsize=14)


    # set axis labels
    plt.xlabel('Drainage area (m$^2$)', fontsize=14)
    plt.ylabel('Gradient (m/m)', labelpad=15, fontsize=14)
    plt.subplots_adjust(left=0.15, hspace=0.3)

    # save and clear the figure
    plt.savefig(OutDirectory+fname_prefix+('_SA_median_SO{}.png'.format(stream_order)), dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

def PlotSlopeAreaVsChi(DataDirectory, fname_prefix):
    """
    Make a summary plot showing a SA plot and a chi plot for all the channels in the basin

    Author: FJC
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_slopes.csv')

    # find out some info
    sources = df.id.unique()

    # set up a figure
    fig,ax = plt.subplots(nrows=1,ncols=2, figsize=(10,4.5), sharex=False, sharey=False)
    # make a big subplot to allow sharing of axis labels
    #fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axes
    #plt.tick_params(labelcolor='none', top='off', bottom='off', left='off', right='off')

    # calculate the channel steepness
    area = df['drainage_area'].values
    med_slopes, lower_per, upper_per, bin_centres, _ = bin_slope_area_data(df['slope'], area)
    gradient, intercept, r, p, std = stats.linregress(bin_centres, med_slopes)
    intercept = float(10**intercept)
    print("Steepness index: {}".format(intercept))
    print("concavity: {}".format(gradient))
    x2 = np.linspace(1,area.max(),100)
    y2 = intercept*x2**(gradient)

    # transform binned data into normal for plotting
    med_slopes = np.array([10**x for x in med_slopes])
    med_areas = np.array([10**x for x in bin_centres])
    lower_per = np.array([10**x for x in lower_per])
    upper_per = np.array([10**x for x in upper_per])
    upper_err = upper_per - med_slopes
    lower_err = med_slopes - lower_per


    # LEFT - slope area plot
    ax[0].grid(color='0.8', linestyle='--', which='both', zorder=1)
    ax[0].scatter(area, df['slope'].values, color='0.5', s=1, zorder=2)
    ax[0].errorbar(med_areas, med_slopes, xerr=None, yerr=[lower_err, upper_err], fmt='o', ms=5, marker='D', mfc='r', mec='k', zorder=3, c='k')
    ax[0].plot(x2, y2, "--", c='k')
    #ax[0].text(0.15, 0.1,'Cluster {}'.format(int(cl)),horizontalalignment='center',verticalalignment='center',transform = ax[0][i].transAxes,fontsize=12)
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    #ax[0].set_xlim(900,10000)
    ax[0].set_ylim(0.0001, 10)
    ax[0].set_title('$k_s$ = {}; $\\theta$ = {}'.format(round(intercept,4), round(abs(gradient),2)), fontsize=14)


    # set axis labels
    ax[0].set_xlabel('Drainage area (m$^2$)', fontsize=14)
    ax[0].set_ylabel('Gradient (m/m)', labelpad=15, fontsize=14)
    plt.subplots_adjust(wspace=0.3,bottom=0.15)

    # RIGHT - chi plot
    chi_df = pd.read_csv(DataDirectory+fname_prefix+'_MChiSegmented.csv')
    norm = mcolors.Normalize(vmin=chi_df['m_chi'].min(),vmax=30)
    ax[1].grid(color='0.8', linestyle='--', which='both', zorder=1)
    sc = ax[1].scatter(chi_df['chi'], chi_df['elevation'], c=chi_df['m_chi'], cmap=cm.viridis, norm=norm, s=1, zorder=2)
    ax[1].set_xlabel('$\chi$ (m)', fontsize=14)
    ax[1].set_ylabel('Elevation (m)', labelpad=15, fontsize=14)

    # add a colourbar_location
    cax = fig.add_axes([0.58, 0.82, 0.2, 0.02])
    cbar = plt.colorbar(sc,cmap=cm.viridis, orientation='horizontal',cax=cax)
    cbar.set_label('$k_s$', fontsize=10)

    # save and clear the figure
    plt.savefig(DataDirectory+fname_prefix+'_SA_all.png', dpi=300)
    plt.clf()
    plt.cla()
    plt.close()

def PlotUniqueStreamsWithLength(DataDirectory, OutDirectory, fname_prefix, step=2, slope_window_size=25):
    """
    Function to make a plot of the number of unique channels you get (at least one non-overlapping
    node with different profile lengths that you analyse).

    Author: FJC
    """
    # read in the original csv
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')

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
        print (n_sources)
        ax.scatter(len, n_sources, s = 25, c='w', edgecolors='k')

    ax.set_xlabel('Profile length (m)')
    ax.set_ylabel('Number of unique channels')

    plt.savefig(OutDirectory+fname_prefix+'_n_channels_with_length.png', dpi=300)
    plt.clf()

def PlotLongitudinalProfiles(DataDirectory, fname_prefix):
    """
    Just make a simple plot of the river long profiles
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_upstream_clustered.csv')

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    sources = df['id'].unique()

    for src in sources:
        this_df = df[df['id'] == src]
        ax.plot(this_df['distance_from_outlet'], this_df['elevation'])

    ax.set_xlabel('Distance from outlet (m)')
    ax.set_ylabel('Elevation (m)')

    plt.savefig(DataDirectory+fname_prefix+'_long_profiles.png', dpi=300)
    plt.clf()

def PlotTrunkChannel(DataDirectory, fname_prefix):
    """
    Make a simple plot of the longest channel. This is mostly to use for the model runs.
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')

    # set up a figure
    fig,ax = plt.subplots(nrows=1,ncols=1, figsize=(6,4), sharex=True, sharey=True)

    trunk_src = df.loc[df['distance_from_outlet'].idxmax()]['id']

    this_df = df[df['id'] == trunk_src]
    ax.grid(color='0.8', linestyle='--', which='major')
    ax.plot(this_df['distance_from_outlet'], this_df['elevation'], c='k')

    ax.set_xlabel('Distance from outlet (m)', fontsize=14)
    ax.set_ylabel('Elevation (m)', fontsize=14)
    #ax.set_xlim(0,2500)
    ax.set_ylim(0,40)
    plt.subplots_adjust(bottom=0.2, left=0.15)

    plt.savefig(DataDirectory+fname_prefix+'_trunk_profile.png', dpi=300)
    plt.clf()

def PlotElevDistanceTrunkChannel(DataDirectory, fname_prefix, stream_order=1):
    """
    Make a plot of the elevation and slope against the distance from the channel head
    For the paper
    """
    df = pd.read_csv(DataDirectory+fname_prefix+'_slopes.csv')

    # set up a figure
    fig = plt.figure(1, facecolor='white')
    gs = plt.GridSpec(100,100,bottom=0.15,left=0.1,right=0.9,top=0.9)
    ax = fig.add_subplot(gs[5:100,10:95])

    trunk_src = df.loc[df['distance_from_outlet'].idxmax()]['id']

    # fit a linear regression
    this_df = df[df['id'] == trunk_src]
    this_dist = this_df['distance_from_outlet'][::-1][12:37]
    #print(this_dist)
    this_elev = this_df['elevation'][12:37]
    slope, intercept, r, p, std = stats.linregress(this_dist, this_elev)
    #print(slope, intercept)
    new_elev = slope* this_dist + intercept

    reg_df = pd.read_csv(DataDirectory+fname_prefix+'_profiles_SO{}.csv'.format(stream_order))
    print(reg_df['distance_from_outlet'][::-1][1:20])
    print(reg_df['reg_dist'][1:20])

    #plotting
    ax.grid(color='0.8', linestyle='--', which='major',zorder=1)
    ax.scatter(this_df['distance_from_outlet'][::-1][:50], this_df['elevation'][:50], edgecolors='k', facecolor='white', s=30, label=None,zorder=2)
    ax.scatter(this_df['distance_from_outlet'][::-1][12:37], this_df['elevation'][12:37], edgecolors='k', facecolor='0.5', s=30, label=None,zorder=3)
    ax.scatter(this_df['distance_from_outlet'][::-1][24:25], this_df['elevation'][24:25], edgecolors='k', facecolor='red', s=40, label='Node of interest',zorder=4)
    #ax.plot(this_dist-1, new_elev-0.1, c='k', ls='--')
    ax.text(30, 8.2, 'Linear fit, $S = {}$'.format(np.round(abs(slope),4)),fontsize=10)
    ax.legend(loc='upper right', fontsize=10)
    ax.set_xlabel('Distance downstream (m)')
    ax.set_ylabel('Elevation (m)')
    plt.title('Slope window size: $W_s = 25$',fontsize=12)
    plt.annotate(s='', xy=(13,9), xytext=(44,6.8), arrowprops=dict(arrowstyle='<->'))
    #plt.arrow(13,9,30,-2.1,length_includes_head=True, head_width=0.2,shape='right')
    #plt.arrow(13,9,30,-2.1,length_includes_head=True, head_width=0.2,shape='right',head_starts_at_zero=False)
    #ax.set_xlim(0,2500)
    #ax.set_ylim(0,35)

    plt.savefig(DataDirectory+fname_prefix+'_trunk_elev_dist.png', dpi=300)
    plt.clf()
