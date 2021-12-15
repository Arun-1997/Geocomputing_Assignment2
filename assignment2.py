############# 1 ################
import os
import json
from osgeo import gdal, ogr
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

sdn_cropsDs = gdal.Open("sdn_crops.tif")

print("There are " + str(nrbands) + " bands")

print("x size: ", sdn_xSize, " y size: ", sdn_ySize)

rasterProjection = sdn_cropsDs.GetProjection()
print("EPSG Code:", int(spatialRef.GetAttrValue("AUTHORITY", 1)))

print("No data value", no_data_value)

rasterGeotransform = sdn_cropsDs.GetGeoTransform()
if rasterGeotransform is not None:
    print()
    ulx = rasterGeotransform[0]
    uly = rasterGeotransform[3]

    extent = [ulx, lry, lrx, uly]

    print("extent:", extent)
    print()

############# 2 ################
print("#" * 10, " Block 2", "#" * 10)
print()

catchment_data_file = "catchments.shp"
city_catchmentsDs = ogr.Open(catchment_data_file)

input_layer = city_catchmentsDs.GetLayer(2)

print("x_min = %.2f x_max = %.2f y_min = %.2f y_max = %.2f" % ())

layerDefinition = input_layer.GetLayerDefn()
print("Number of fields: " + str(fieldCount))

fields = []

print("Fields: ", fields)

print("Number of features: " + str(layerFeatureNum))

small_and_largeCatchments = {
    "small": {"area": 0, "objectid": ""},
    "large": {"area": 0, "objectid": ""},
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

print()

############# 3 ################

print("#" * 20, " Block 3", "#" * 20)
print()

most_important_crop_per_catchment = {}
for feature in input_layer:
    sdn_crops_clippedDs = gdal.Warp(
        cutlineDSName=catchment_data_file,
        cutlineWhere="objectid = {}".format(objectId),
        cropToCutline=False,
    )

    sdn_crops_prodArray = gdarr.DatasetReadAsArray(input_layer)
    sdn_crops_prodArray[sdn_crops_prodArray > no_data_value] = 0

    most_important_crop_per_catchment[objectId] = [
        most_important_crop,
        int(most_important_crop_amount),
        0,
    ]


############# 4 ################
print("#" * 10, " Block 4", "#" * 10)
print()


cred_path = r"..\credentials.json"

# Read the login details from json file for safe connection.
with open(cred_path, "r") as login:
    db_con_data = json.load(login)

db_host = db_con_data["host"]

conn_string = "host='%s' port='%d' user='%s' password='%s' dbname='%s'" % (
    db_host,
    db_port,
    db_username,
    db_pwd,
    db_name,
)

with psycopg2.connect(conn_string) as conn:
    with conn.cursor(cursor_factory=RealDictCursor) as cursor:
        selectQuery = """
        """
        cursor.execute(selectQuery)
        records = cursor.fetchall()
        for record in records:
            most_important_crop_per_catchment[record["objectid"]][2] = int(
                sum(crops_cons_total)
            )


print(most_important_crop_per_catchment)
print()

############# 5 ################
print("#" * 10, " Block 5", "#" * 10)
print()

output_data_file = "sdn_market_dynamics.json"
driver = ogr.GetDriverByName("GeoJSON")
output_ds = driver.CreateDataSource(output_data_file)

field_names = ["crop", "prod", "cons"]

for infeature in output_layer:
    objectId = infeature.GetFieldAsInteger("objectid")
    if objectId in most_important_crop_per_catchment.keys():
        outfeature.SetField(
            field_names[0], most_important_crop_per_catchment[objectId][0]
        )
        outfeature.SetGeometry(outfeature.GetGeometryRef())
        output_layer.AddFeature(outfeature)

input_layer = None
output_layer = None
city_catchmentsDs = None
output_ds = None

############# 6 ################
print("#" * 10, " Block 6", "#" * 10)
print()

speedups.disable()

print(catchment_df)
print()

catchment_df.plot(column="diff", legend="True")
plt.show()

############# END ################
