############# 1 ################
import os
import sys
import json
from osgeo import gdal, ogr, osr, gdal_array as gdarr
import psycopg2
from psycopg2.extras import RealDictCursor
from matplotlib import pyplot as plt
import geopandas as gpd
from shapely import speedups
import numpy as np
try:
    import mapclassify
except ModuleNotFoundError:
    print("Warning! Map Classify not installed. Legend for map might not be properly classified!")
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
    relative_file_path = path_sep.join([data_dir, raster_file_name])
    absolute_file_path = os.path.abspath(relative_file_path)
    print(f"The raster file is not opened. Kindly check the path and file name of the raster: \n{absolute_file_path}")
    sys.exit()

# Defining nrbands to get the number of bands from the raster
nrbands = sdn_cropsDs.RasterCount
print("There are " + str(nrbands) + " bands")

# 1.1 Defining sdn_xSize to get the Raster X size i.e width
sdn_xSize = sdn_cropsDs.RasterXSize

# Defining sdn_xSize to get the Raster Y size i.e height
sdn_ySize = sdn_cropsDs.RasterYSize
print("x size: ", sdn_xSize, " y size: ", sdn_ySize)

rasterProjection = sdn_cropsDs.GetProjection()

# 1.2 Get the spatial reference from ogr module and set the well known text (wkt)
# to the raster projection
spatialRef = osr.SpatialReference(wkt=rasterProjection)
print("EPSG Code:", int(spatialRef.GetAttrValue("AUTHORITY", 1)))

rasterGeotransform = sdn_cropsDs.GetGeoTransform()
if rasterGeotransform is not None:
    print()
    ulx = rasterGeotransform[0]
    uly = rasterGeotransform[3]

    # 1.3 The formula to get the lower right x and lower right y coordinates
    # To get lrx (lower right x), Add the ulx to the product of number of columns
    # i.e. RasterXsize and pixel resolution along the w-e. 
    # Also add the RasterYsize to the rotation of x angle (In case there is skewness)
    # The same is done vice versa for the lower right y (lry) 
    lrx = ulx + sdn_xSize * rasterGeotransform[1] + sdn_ySize * rasterGeotransform[2]
    lry = uly + sdn_ySize * rasterGeotransform[5] + sdn_xSize * rasterGeotransform[4]
    extent = [ulx, lry, lrx, uly]

    print("extent:", extent)
    print()

# 1.4 Get the No data value for all the bands and assign it to a dictionary
# with key as band number
no_data_value = {i: sdn_cropsDs.GetRasterBand(i).GetNoDataValue()
                 for i in range(1, nrbands + 1)}
print("NoData value", json.dumps(no_data_value, indent=1))

# Checking if the No data values for all the bands are
# same using the all() method. True if same, False if not same
no_data_value_list = list(no_data_value.values())
same_no_data_value = all(i == no_data_value_list[0] for i in no_data_value_list)
print(f"NoData value same for all bands? {same_no_data_value}")

# 1.5 Running a comprehensive loop to get the minimum and maximum values of all the bands
raster_statistics = {f"Band {i}": {"Minimum": sdn_cropsDs.GetRasterBand(i).GetMinimum(),
                                   "Maximum": sdn_cropsDs.GetRasterBand(i).GetMaximum()}
                     for i in range(1, nrbands + 1)}

print("Minimum and maximum values of the bands : ")
print(json.dumps(raster_statistics, indent=2))

############# 2 ################
print("#" * 10, " Block 2", "#" * 10)

catchment_data_file = "catchments.shp"
city_catchmentsDs = ogr.Open(catchment_data_file)

# Extracting the shapefile layer
input_layer = city_catchmentsDs.GetLayer(0)

# 2.1
layerExtents = input_layer.GetExtent()  # Getting the layer extents
print("x_min = %.2f x_max = %.2f y_min = %.2f y_max = %.2f" % (layerExtents[0],
                                                               layerExtents[1], layerExtents[2], layerExtents[3]))

# 2.2
layerDefinition = input_layer.GetLayerDefn()  # Getting the layer definition for further metadata
fieldCount = layerDefinition.GetFieldCount()
print("Number of fields: " + str(fieldCount))

# Printing the field name
fields = [layerDefinition.GetFieldDefn(i).GetName() for i in range(fieldCount)]
print("Fields: ", fields)

# 2.3
layerFeatureNum = input_layer.GetFeatureCount()  # Getting the number of features
print("Number of features: " + str(layerFeatureNum))

# 2.4
# Creating a dictionary of areas and object id of all the catchment
areas = {int(feature["objectid"]): feature.GetGeometryRef().GetArea()
         for feature in input_layer}

# Extracting the minimum area and the corresponding object id
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

print("#" * 20, " Block 3", "#" * 20)
print()

most_important_crop_per_catchment = {}
for feature in input_layer:

    objectId = feature.GetFieldAsString("objectid")
    sdn_crops_clippedDs = gdal.Warp("",
                                    sdn_cropsDs, format="Mem",
                                    cutlineDSName=catchment_data_file,
                                    cutlineWhere="objectid = {}".format(objectId),
                                    cropToCutline=False
                                    )

    sdn_crops_prodArray = gdarr.DatasetReadAsArray(sdn_crops_clippedDs)

    # Set the No data value to zero 
    # If noData value is same, for all the bands the same value can be taken 
    # from first band and set to zero
    if same_no_data_value:
        nodata_value = no_data_value.get(1)
        sdn_crops_prodArray[sdn_crops_prodArray == nodata_value] = 0
    # If noData values are different for each band, then loop over them for each band
    # and set zero to those values
    else:
        for i, j in no_data_value.items():
            sdn_crops_prodArray_band = sdn_crops_prodArray[i - 1]
            sdn_crops_prodArray_band[sdn_crops_prodArray_band == j] = 0

    # Getting the Numberof bands, ysize and xsize from the array shape
    nr_bands, ysize, xsize = sdn_crops_prodArray.shape

    # Reshaping the array into a 2 dimensional array where each row 
    # has pixel values of the same band and columns are the same position 
    # of different bands
    # 2D Array (number of bands * number of pixels in raster)
    sdn_crops_prod2dArray = sdn_crops_prodArray.reshape(nr_bands, ysize * xsize)

    # Getting the sum of all the pixels in each band
    sdn_crop_prodSum = np.sum(sdn_crops_prod2dArray, axis=1)

    # Getting the maximum value and its index from the sum of bands 
    # This gives the crop that has the most production for the feature
    most_important_crop_index, most_important_crop_amount = np.argmax(sdn_crop_prodSum), np.max(sdn_crop_prodSum)

    # If maximum value is zero, it means that there are no crops produced in that catchment
    # None is assigned to crop in that case
    most_important_crop = crop_list[most_important_crop_index] if int(most_important_crop_amount) != 0 else None

    most_important_crop_per_catchment[int(objectId)] = [
        most_important_crop,
        int(most_important_crop_amount),
        0,
    ]

# ############# 4 ################
print("#" * 10, " Block 4", "#" * 10)
print()

cred_path = r"../credentials.json"

# Read the login details from json file for safe connection.
with open(cred_path, "r") as login:
    db_con_data = json.load(login)

# Defining the credentials to connect to the server from the json file
db_host = db_con_data["host"]
db_port = db_con_data["port"]
db_username = db_con_data["user"]
db_pwd = db_con_data["pwd"]
db_name = db_con_data["db"]

conn_string = "host='%s' port='%d' user='%s' password='%s' dbname='%s'" % (
    db_host,
    db_port,
    db_username,
    db_pwd,
    db_name,
)

# Connecting to Postgresql driver with the conn_string parameters
# from the json file
with psycopg2.connect(conn_string) as conn:
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        # SQL query to get the cities that correspond to each catchment
        # In the given data, the catchments are service areas for each cities
        # hence, every city will correspond to one catchment. A spatial relation is
        # made to join the city attributes with the catchment 
        selectQuery = """
        select c.id as city_id,cat.id as catchment_id,c.city,c.pop2017,cat.objectid
        from sudan.cities as c
        join sudan.catchments as cat on st_within(c.geom,cat.geom)
        """
        # Run the SQL query
        cursor.execute(selectQuery)

        # Get the records from the query
        records = cursor.fetchall()

        # Loop over the records
        for record in records:
            most_important_crop = most_important_crop_per_catchment[record["objectid"]][0]
            # The crop consumption is the product of food per capita of the crop (provided above) 
            # with the population from the cities table
            pop = float(record["pop2017"])
            # If no crop is important i.e. catchment has None value, returns zero for crop consumption
            crops_cons_total = food_per_capita[most_important_crop] * pop if most_important_crop is not None else 0
            most_important_crop_per_catchment[record["objectid"]][2] = int(crops_cons_total)

print(json.dumps(most_important_crop_per_catchment, indent=3))
print()

############# 5 ################
print("#" * 10, " Block 5", "#" * 10)

data_directory = os.getcwd()
output_data_file = "sdn_market_dynamics.geojson"
file_path_output = os.path.join(data_directory, output_data_file)  # saving in the current directory

if os.path.exists(file_path_output):
    print("File already existed, removing the old version")
    os.remove(file_path_output)  # remove the file if it already exists to create the new file

driver = ogr.GetDriverByName('GeoJSON')
output_ds = driver.CreateDataSource(file_path_output)

input_layer_srs = input_layer.GetSpatialRef()  # using the same crs info as the original layer

output_layer = output_ds.CreateLayer(output_data_file, input_layer_srs, ogr.wkbPolygon)
oft_int = "OFTInteger"
oft_str = "OFTString"
field_names = [("objectid",oft_int),("crop",oft_str), ("prod",oft_int), ("cons",oft_int)]

for val in field_names:
    field_name = ogr.FieldDefn(val[0], getattr(ogr,val[1]))    # type can be changed -- Best way to do it - discuss
    field_name.SetWidth(24)
    output_layer.CreateField(field_name)

for infeature in input_layer:
    objectId = infeature.GetFieldAsInteger("objectid")
    if objectId in most_important_crop_per_catchment.keys():
        outfeature = ogr.Feature(output_layer.GetLayerDefn())
        
        for i in range(len(field_names)):
            fld_attr = field_names[i]
            if fld_attr[0] == "objectid":
                outfeature.SetField(fld_attr[0],objectId)
            else:
                outfeature.SetField(fld_attr[0], most_important_crop_per_catchment[objectId][i-1])
        
        outfeature.SetGeometry(infeature.GetGeometryRef())
        output_layer.CreateFeature(outfeature)
        outfeature = None

input_layer = None
output_layer = None
city_catchmentsDs = None
output_ds = None

# ############# 6 ################
print("#" * 10, " Block 6", "#" * 10)
print()

speedups.disable()
catchment_df = gpd.read_file(file_path_output)


# Add a new column in the geodataframe for the difference
catchment_df["diff"] = catchment_df["prod"] - catchment_df["cons"]
print(catchment_df)

# plotting the map classified by 
# diff value and displaying the chloropleth map
catchment_df.plot(column="diff", legend="True",cmap = "RdYlBu_r",legend_kwds={"label": "Surplus/Deficit",
                                                                                            'orientation': "vertical"})
plt.title("Surplus or Deficit of crops in Sudan")
plt.show()


#
# ############# END ################
