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

DataDirectory = '/home/clubb/OneDrive/river_clusters/Pozo/'
OutDirectory = DataDirectory+'threshold_0/'
fname_prefix = 'Pozo_DTM'
raster_name = 'pozo_geol_WGS84_reclass.tif'
BoxPlotByCluster(DataDirectory, OutDirectory, fname_prefix,  raster_name, stream_order=1)
