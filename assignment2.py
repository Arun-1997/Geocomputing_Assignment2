############# 1 ################
import os
import sys
import json
from osgeo import gdal, ogr,osr
import psycopg2
from psycopg2.extras import RealDictCursor
from matplotlib import pyplot as plt
import geopandas as gpd
from shapely import speedups

##################################
print("#" * 10, " Block 1", "#" * 10)
print()

# Top 10 consumed crops according to FAO stats
crop_list = [
    "bananas",
    "fruits",
    "groundnuts",
    "millets",
    "potatoes",
    "pulses",
    "sorghum",
    "sweet potatoes",
    "vegetables",
    "wheat",
]

food_per_capita = {
    "groundnuts": 0.0411,
    "millets": 0.0222,
    "bananas": 0.0223,
    "vegetables": 0.0874,
    "sweet potatoes": 0.0061,
    "wheat": 0.0650,
    "pulses": 0.0085,
    "fruits": 0.0403,
    "sorghum": 0.1010,
    "potatoes": 0.0112,
}

data_dir = r".\data"
os.chdir(data_dir)
# Providing the file path seperator and the raster file name in a variable
# in case of file name change
path_sep = "/"
raster_file_name = "sdn_crops.tif"
sdn_cropsDs = gdal.Open(raster_file_name)

# Check if the sdn_cropsDS is opened or not. If not opened,
# print an warning message and exit

if sdn_cropsDs is None:
    relative_file_path = path_sep.join([data_dir,raster_file_name])    
    absolute_file_path = os.path.abspath(relative_file_path)
    print(f"The raster file is not opened. Kindly check the path and file name of the raster: \n{absolute_file_path}")
    sys.exit()

# Defining nrbands to get the number of bands from the raster
nrbands = sdn_cropsDs.RasterCount
print("There are " + str(nrbands) + " bands")

# Defining sdn_xSize to get the Raster X size i.e width
sdn_xSize = sdn_cropsDs.RasterXSize

# Defining sdn_xSize to get the Raster Y size i.e height
sdn_ySize = sdn_cropsDs.RasterYSize
print("x size: ", sdn_xSize, " y size: ", sdn_ySize)

rasterProjection = sdn_cropsDs.GetProjection()

# Get the spatial reference from ogr module and set the well known text (wkt)
# to the raster projection
spatialRef = osr.SpatialReference(wkt=rasterProjection)
print("EPSG Code:", int(spatialRef.GetAttrValue("AUTHORITY", 1)))


# Get the No data value for all the bands and assign it to a dictionary
# with key as band number
no_data_value = {i:sdn_cropsDs.GetRasterBand(i).GetNoDataValue()
                                for i in range(1,nrbands+1)}
print("NoData value", json.dumps(no_data_value,indent=1))

# Checking if the No data values for all the bands are
# same using the all() method. True if same, False if not same
same_no_data_value = all(no_data_value.values())
print(f"NoData value same for all bands? {same_no_data_value}")


rasterGeotransform = sdn_cropsDs.GetGeoTransform()
if rasterGeotransform is not None:
    print()
    ulx = rasterGeotransform[0]
    uly = rasterGeotransform[3]

    # The formula to get the lower right x and lower right y coordinates
    # To get lrx (lower right x), Add the ulx to the product of number of columns
    # i.e. RasterXsize and pixel resolution 
    # along the w-e. Also add the RasterYsize to the rotation of x angle (In case there is skewness)
    # The same is done vice versa for the lower right y (lry) 
    lrx = ulx + sdn_xSize*rasterGeotransform[1] + sdn_ySize*rasterGeotransform[2]
    lry = uly + sdn_ySize*rasterGeotransform[5] + sdn_xSize*rasterGeotransform[4]
    extent = [ulx, lry, lrx, uly]

    print("extent:", extent)
    print()

# Running a comprehensive loop to get the minimum and maximum values of all the bands
raster_statistics = {f"Band {i}":{"Minimum":sdn_cropsDs.GetRasterBand(i).GetMinimum(),
                                  "Maximum":sdn_cropsDs.GetRasterBand(i).GetMaximum()}
                                  for i in range(1,nrbands+1)}


print("Minimum and maximum values of the bands : ")
print(json.dumps(raster_statistics,indent=2))

############# 2 ################
print("#" * 10, " Block 2", "#" * 10)

catchment_data_file = "catchments.shp"
city_catchmentsDs = ogr.Open(catchment_data_file)

#Extracting the shapefile layer
input_layer = city_catchmentsDs.GetLayer(0)

#2.1
layerExtents = input_layer.GetExtent()      #Getting the layer extents
print("x_min = %.2f x_max = %.2f y_min = %.2f y_max = %.2f" % (layerExtents[0],
layerExtents[1], layerExtents[2], layerExtents[3]))

#2.2
layerDefinition = input_layer.GetLayerDefn()    #Getting the layer definition for further metadata
fieldCount = layerDefinition.GetFieldCount()
print("Number of fields: " + str(fieldCount))

#Printing the field name
fields = [layerDefinition.GetFieldDefn(i).GetName() for i in range(fieldCount)]
print("Fields: ", fields)

#2.3
layerFeatureNum = input_layer.GetFeatureCount()     #Getting the number of features
print("Number of features: " + str(layerFeatureNum))

#2.4

#Creating a dictionary of areas and object id of all the catchment
areas = {}
for feature in input_layer:
    featureArea = feature.GetGeometryRef().Area()
    areas[int(feature['objectid'])] = featureArea

#Extracting the minimum area and the corresponding object id
small_and_largeCatchments = {
    "small": {"area": min(areas.values()), "objectid": str(min(areas, key=areas.get))},
    "large": {"area": max(areas.values()), "objectid": str(max(areas, key=areas.get))},
}


print(
    "largest catchment: ",
    small_and_largeCatchments["large"]["objectid"],
    " area: ",
    small_and_largeCatchments["large"]["area"],
)

print(
    "smallest catchment: ",
    small_and_largeCatchments["small"]["objectid"],
    " area: ",
    small_and_largeCatchments["small"]["area"],
)


# ############# 3 ################
#
# print("#" * 20, " Block 3", "#" * 20)
# print()
#
# most_important_crop_per_catchment = {}
# for feature in input_layer:
#     sdn_crops_clippedDs = gdal.Warp(
#         cutlineDSName=catchment_data_file,
#         cutlineWhere="objectid = {}".format(objectId),
#         cropToCutline=False,
#     )
#
#     sdn_crops_prodArray = gdarr.DatasetReadAsArray(input_layer)
#     sdn_crops_prodArray[sdn_crops_prodArray > no_data_value] = 0
#
#     most_important_crop_per_catchment[objectId] = [
#         most_important_crop,
#         int(most_important_crop_amount),
#         0,
#     ]
#
#
# ############# 4 ################
# print("#" * 10, " Block 4", "#" * 10)
# print()
#
#
# cred_path = r"..\credentials.json"
#
# # Read the login details from json file for safe connection.
# with open(cred_path, "r") as login:
#     db_con_data = json.load(login)
#
# db_host = db_con_data["host"]
#
# conn_string = "host='%s' port='%d' user='%s' password='%s' dbname='%s'" % (
#     db_host,
#     db_port,
#     db_username,
#     db_pwd,
#     db_name,
# )
#
# with psycopg2.connect(conn_string) as conn:
#     with conn.cursor(cursor_factory=RealDictCursor) as cursor:
#         selectQuery = """
#         """
#         cursor.execute(selectQuery)
#         records = cursor.fetchall()
#         for record in records:
#             most_important_crop_per_catchment[record["objectid"]][2] = int(
#                 sum(crops_cons_total)
#             )
#
#
# print(most_important_crop_per_catchment)
# print()
#
# ############# 5 ################
# print("#" * 10, " Block 5", "#" * 10)
# print()
#
# output_data_file = "sdn_market_dynamics.json"
# driver = ogr.GetDriverByName("GeoJSON")
# output_ds = driver.CreateDataSource(output_data_file)
#
# field_names = ["crop", "prod", "cons"]
#
# for infeature in output_layer:
#     objectId = infeature.GetFieldAsInteger("objectid")
#     if objectId in most_important_crop_per_catchment.keys():
#         outfeature.SetField(
#             field_names[0], most_important_crop_per_catchment[objectId][0]
#         )
#         outfeature.SetGeometry(outfeature.GetGeometryRef())
#         output_layer.AddFeature(outfeature)
#
# input_layer = None
# output_layer = None
# city_catchmentsDs = None
# output_ds = None
#
# ############# 6 ################
# print("#" * 10, " Block 6", "#" * 10)
# print()
#
# speedups.disable()
#
# print(catchment_df)
# print()
#
# catchment_df.plot(column="diff", legend="True")
# plt.show()
#
# ############# END ################
