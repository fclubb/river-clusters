#---------------------------------------------------------------------#
# Raster plotting functions for clustering analysis
# Developed by Fiona Clubb
#              Bodo Bookhagen
#              Aljoscha Rheinwalt
# University of Potsdam
#---------------------------------------------------------------------#

from LSDPlottingTools import LSDMap_VectorTools as VT
from LSDPlottingTools import LSDMap_GDALIO as IO
from LSDPlottingTools import LSDMap_BasicManipulation as BM
from LSDMapFigure import PlottingRaster
import pandas as pd
import matplotlib.pyplot as plt
import os
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar


def PlotElevationWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, cbar_loc='right', custom_cbar_min_max = []):
    """
    Make a plot of the raster with the channels coloured by the cluster
    value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

    Args:
        stream_order: the stream order of the profiles that you are analysing
        cbar_loc: location of the colourbar, can be right, top, bottom, left, or none.
        custom_cbar_min_max: list of [min, max] to recast the raster to for display.

    Author: FJC
    """
    import LSDPlottingTools as LSDP
    from LSDMapFigure.PlottingRaster import MapFigure

    df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
    cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


    # set figure sizes based on format
    fig_width_inches = 4.92126

    # some raster names
    raster_ext = '.bil'
    BackgroundRasterName = fname_prefix+raster_ext
    HSName = fname_prefix+'_hs'+raster_ext

    if not os.path.isfile(DataDirectory+HSName):
        # make a hillshade
        BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

    MF = MapFigure(HSName, DataDirectory,coord_type="UTM",colourbar_location = cbar_loc)
    MF.add_drape_image(BackgroundRasterName,DataDirectory,colourmap = 'gray', alpha=0.8, colorbarlabel = "Elevation (m)",colour_min_max = custom_cbar_min_max)

    # # create the map figure
    # MF = MapFigure(BackgroundRasterName, DataDirectory,coord_type="UTM",colour_min_max = custom_cbar_min_max)

    clusters = cluster_df.cluster_id.unique()
    for cl in clusters:
        # plot the whole channel network in black
        ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
        MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='w',manual_size=1, zorder=2, alpha=1)
        # plot the clustered profiles in the correct colour
        this_df = cluster_df[cluster_df.cluster_id == cl]
        this_colour = str(this_df.colour.unique()[0])
        ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
        MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2)

    MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_elev_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

def PlotHillshadewithClusters(DataDirectory, OutDirectory, fname_prefix,stream_order=1):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.9',manual_size=2, zorder=1, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=3)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_hs_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

def PlotLithologyWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, shapefile_name = 'geol.shp', geol_field = 'geol'):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 4.92126

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            print("Making a hillshade for you")
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")
        res = IO.GetUTMMaxMin(DataDirectory+BackgroundRasterName)[0]

        #rasterise the shapefile
        new_shp, geol_dict = VT.geologic_maps_modify_shapefile(DataDirectory+ shapefile_name, geol_field)
        lith_raster = VT.Rasterize_geologic_maps_pythonic(new_shp, res, geol_field)
        MF.add_drape_image(lith_raster, "", colourmap=plt.cm.jet, alpha=0.4, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='0.9',manual_size=2, zorder=1, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=3)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_lith_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure


def PlotRasterLithologyWithClusters(DataDirectory, OutDirectory, fname_prefix, stream_order=1, geol_raster = 'geol'):
        """
        Make a hillshade of the raster with the channels coloured by the cluster
        value. Rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        df = pd.read_csv(DataDirectory+fname_prefix+'_all_tribs.csv')
        cluster_df = pd.read_csv(OutDirectory+fname_prefix+'_profiles_clustered_SO{}.csv'.format(stream_order))


        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        #geology
        LithName = geol_raster
        print("The geology raster is"+LithName)
        MF.add_drape_image(LithName, DataDirectory, colourmap=plt.cm.jet, alpha=0.5, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        clusters = cluster_df.cluster_id.unique()
        for cl in clusters:
            # plot the whole channel network in black
            ChannelPoints = LSDP.LSDMap_PointData(df, data_type="pandas", PANDEX = True)
            MF.add_point_data(ChannelPoints,show_colourbar="False", unicolor='white',manual_size=1, zorder=2, alpha=0.5)
            # plot the clustered profiles in the correct colour
            this_df = cluster_df[cluster_df.cluster_id == cl]
            this_colour = str(this_df.colour.unique()[0])
            ClusteredPoints = LSDP.LSDMap_PointData(this_df, data_type = "pandas", PANDEX = True)
            MF.add_point_data(ClusteredPoints,show_colourbar="False",zorder=100, unicolor=this_colour,manual_size=2)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = OutDirectory+fname_prefix+'_lith_clusters_SO{}.png'.format(stream_order), FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

def PlotRasterLithology(DataDirectory, fname_prefix, geol_raster = 'geol'):
        """
        Make a hillshade and rasterise a geology shapefile and drape on the top.
        Uses the LSDPlottingTools libraries. https://github.com/LSDtopotools/LSDMappingTools

        Args:
            stream_order: the stream order of the profiles that you are analysing
            shapefile_name: name of the lithology shapefile
            geol_field: the field of the shapefile that has the lithology information

        Author: FJC
        """
        import LSDPlottingTools as LSDP
        from LSDMapFigure.PlottingRaster import MapFigure

        # set figure sizes based on format
        fig_width_inches = 8

        # some raster names
        raster_ext = '.bil'
        BackgroundRasterName = fname_prefix+raster_ext
        HSName = fname_prefix+'_hs'+raster_ext

        if not os.path.isfile(DataDirectory+HSName):
            # make a hillshade
            BM.GetHillshade(DataDirectory+BackgroundRasterName, DataDirectory+HSName)

        # create the map figure
        MF = MapFigure(HSName, DataDirectory,coord_type="UTM")

        #geology
        LithName = geol_raster
        print("The geology raster is"+LithName)
        MF.add_drape_image(LithName, DataDirectory, colourmap=plt.cm.jet, alpha=0.5, show_colourbar = False, discrete_cmap=True, cbar_type=int,mask_value=0)

        MF.save_fig(fig_width_inches = fig_width_inches, FigFileName = DataDirectory+fname_prefix+'_lith.png', FigFormat='png', Fig_dpi = 300, fixed_cbar_characters=6, adjust_cbar_characters=False) # Save the figure

if __name__ == '__main__':

    DataDirectory = '/home/clubb/Data_for_papers/river_clusters/Pozo/'
    fname_prefix  = 'Pozo_DTM_basin_208'
    stream_order = 1
    shp = 'pozo_geol_WGS84.shp'
    field = 'Lithology'
    #PlotLithologyWithClusters(DataDirectory, fname_prefix, stream_order, shp, field)
    PlotRasterLithologyWithClusters(DataDirectory, fname_prefix, stream_order, geol_raster='pozo_geol_WGS84_new_padded')
