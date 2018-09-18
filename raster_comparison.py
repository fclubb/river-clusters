# Scripts to compare the results of the clustering to different rasters
# FJC 12/09/18

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from LSDPlottingTools import LSDMap_GDALIO as IO
from LSDPlottingTools import LSDMap_PointTools as PT

def BoxPlotByCluster(DataDirectory, OutDirectory, fname_prefix,  raster_name, stream_order=1):
    """
    Make a boxplot of the results of the clustering compared to the raster specified
    by raster_name
    """
    #df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))

    # read in the raster
    raster_ext = '.bil'
    this_raster = IO.ReadRasterArrayBlocks(DataDirectory+raster_name)
    EPSG_string = IO.GetUTMEPSG(DataDirectory+raster_name)
    NDV, xsize, ysize, GeoT, Projection, DataType = IO.GetGeoInfo(DataDirectory+raster_name)
    CellSize,XMin,XMax,YMin,YMax = IO.GetUTMMaxMin(DataDirectory+raster_name)

    pts = PT.LSDMap_PointData(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order),data_type ='csv')
    easting, northing = pts.GetUTMEastingNorthing(EPSG_string=EPSG_string)
    cluster_id = pts.QueryData('cluster_id', PANDEX=True)
    clusters = list(set(cluster_id))

    # dict for the data
    data = {k: [] for k in clusters}
    for x, (i, j) in enumerate(zip(northing, easting)):
    # convert to rows and cols
        X_coordinate_shifted_origin = j - XMin;
        Y_coordinate_shifted_origin = i - YMin;

        col_point = int(X_coordinate_shifted_origin/CellSize);
        row_point = (ysize - 1) - int(round(Y_coordinate_shifted_origin/CellSize));
        # check for data at this cell
        this_value = this_raster[row_point][col_point]
        if not np.isnan(this_value):
            if this_value < 10:
                # get the cluster id
                data[cluster_id[x]].append(this_value)

    print(data)

    # now make a boxplot
    labels, these_data = [*zip(*data.items())]  # 'transpose' items to parallel key, value lists

    plt.boxplot(these_data)
    plt.xticks(range(1, len(labels) + 1), labels)
    plt.show()

def BoxPlotGradient(DataDirectory, OutDirectory, fname_prefix, stream_order=1):
    """
    Make a boxplot showing the channel gradient stats for each cluster
    """

    # read the csv and get some info
    df = pd.read_csv(OutDirectory+fname_prefix+"_profiles_clustered_SO{}.csv".format(stream_order))
    colors = df['colour'].unique()

    # set props for fliers
    flierprops = dict(marker='o', markerfacecolor='none', markersize=1,
                  linestyle='none', markeredgecolor='k')

    # make the boxplot and return the dict with the boxplot properties
    bp_dict = df.boxplot(column=['slope'], by=['cluster_id'], return_type='both', patch_artist= True,flierprops=flierprops, figsize=(5,5))
    # make the median lines black
    #[[item.set_color('k') for item in bp_dict[key]['medians']] for key in bp_dict.keys()]

    # change the colours based on the cluster ID
    for row_key, (ax,row) in bp_dict.iteritems():
        plt.xlabel('')
        ax.set_title('')
        j=-1 #stupid thing because there are double the number of caps and whiskers compared to boxes
        for i,cp in enumerate(row['caps']):
            if i%2==0:
                j+=1
            cp.set(color=colors[j])
        j=-1
        for i,wh in enumerate(row['whiskers']):
            if i%2==0:
                j+=1
            wh.set_color(colors[j])
        for i,box in enumerate(row['boxes']):
            box.set_facecolor(colors[i])
            box.set_alpha(0.7)
            box.set_edgecolor(colors[i])
        for i,med in enumerate(row['medians']):
            med.set(color=colors[i])
        for i,pt in enumerate(row['fliers']):
            pt.set_markeredgecolor(colors[i])


    ax.grid(color='0.8', linestyle='--', which='major', zorder=1)
    ax.set_xticklabels(['Canada shale', 'Breccia/volcaniclastics'], fontsize=12)
    #print(boxplot)
    plt.suptitle('')
    ax.set_ylabel('Gradient (m/m)', fontsize=14)
    plt.subplots_adjust(left=0.2)
    plt.savefig(OutDirectory+fname_prefix+'_boxplot.png', dpi=300)
    plt.clf()

DataDirectory = '/home/clubb/OneDrive/river_clusters/Pozo/'
OutDirectory = DataDirectory+'threshold_0/'
fname_prefix = 'Pozo_DTM'
raster_name = 'pozo_geol_WGS84_reclass.tif'
# BoxPlotByCluster(DataDirectory, OutDirectory, fname_prefix,  raster_name, stream_order=1)
BoxPlotGradient(DataDirectory,OutDirectory,fname_prefix,stream_order=1)
