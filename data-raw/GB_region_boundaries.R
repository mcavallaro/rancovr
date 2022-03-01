
crs = structure(list(input = "4326", wkt = "GEOGCS[\"WGS 84\",\n    DATUM[\"World Geodetic System 1984\",\n        SPHEROID[\"WGS 84\",6378137.0,298.257223563,\n            AUTHORITY[\"EPSG\",\"7030\"]],\n        AUTHORITY[\"EPSG\",\"6326\"]],\n    PRIMEM[\"Greenwich\",0.0,\n        AUTHORITY[\"EPSG\",\"8901\"]],\n    UNIT[\"degree\",0.017453292519943295],\n    AXIS[\"Geodetic longitude\",EAST],\n    AXIS[\"Geodetic latitude\",NORTH],\n    AUTHORITY[\"EPSG\",\"4326\"]]"), class = "crs")
GB_region_boundaries = sf::st_read('data-raw/bdline_essh_gb/Data/Supplementary_Historic_European_Region/historic_european_region_region.shx')
GB_region_boundaries = sf::st_transform(GB_region_boundaries, crs=crs)
save(GB_region_boundaries, file = "data/GB_region_boundaries.RData")
