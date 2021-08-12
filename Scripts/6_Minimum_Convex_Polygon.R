# Credits -----------------------------------------------------------------

# Code created by Luara Tourinho (https://github.com/luaratourinho)
# Changed by Julia Niemeyer
# 04 August 2021


# First install packages
# Then run the library required

library(maps)
library(mapdata)
library(rgdal)
library(raster)
library(rgeos)
library(dplyr)


# Some examples of projections --------------------------------------------

# Creating two objects: (1) with the original projection of your data
# and (2) with an equal area projection
dir.create("./outputs")

## Entrar com a planilha limpa pós spthin
sp <- read.table("./data/03_clean_df_thin_1_BSF.csv",header=TRUE, sep=",")# %>%
# select(species, lon, lat) %>%
#filter(species == "Acestrorhynchus_britskii")
sp_names <- unique(sp$species)
# running for one species
#sp.n = sp_names[[2]]

#Começa o for pras sps

#for (a in 1:length(sp_names)) {

 # message("starting the analysis for ", paste0(sp_names[a]))

# WGS84
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")

# South America Albers Equal Area Conic
crs.albers <-
  CRS("+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60
                  +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs")

# Lambert Azimuthal Equal Area (ideal for Pacific Ocean islands)
#crs.azim <-
 # CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
 #               +datum=WGS84 +units=m +no_defs")



# Building a minimum convex polygon ---------------------------------------


# Read your table of occurrence records

##Planilha de ocorrência pós spthin aqui
occurrence_records <- read.table("./data/03_clean_df_thin_1_BSF.csv",header=TRUE, sep=",") %>%
   filter(species == paste0(sp_names[a])) %>%
  select(species, lon, lat)

# Check the column names for your coordinates and place below within c("","")
coordinates(occurrence_records) <- c("lon", "lat")

# Define the original projection of your occurrence records
proj4string(occurrence_records) <- crs.wgs84

# Reproject to an equal area projection
occurrence_records_albers <-
  spTransform(occurrence_records, crs.albers)

# Minimum convex polygon (mpc) function
mpc <- gConvexHull(occurrence_records_albers)



# Adding a buffer to the mpc -----------------------------------------------


# You can choose a fixed value for a buffer, for example
# If you want to 500km of ratio, use b = 500000, because b is in meters
# If you want to a proportional value as suggested by Barve et al. 2011
# (doi:10.1016/j.ecolmodel.2011.02.011), use the function gArea.
# The function gArea is in m2, so multipling by 0.000001 you can get a area in km2
# For example, if you want to draw a buffer of 20%, in km, around of the mpc


# 20% of the polygon area is 2e-07
# draw a buffer of 2° (111km) around of the minimum convex polygon (adaptado de Barve et al. 2011)
## width -> can be replaced for another buffer (in meters) for species that have large distributions based on the literature
b <- (111000) ##Mudar valor pro valor que decidi
buffer_mpc <- gBuffer(mpc, width = b)



# An example of mpc utility -------------------------------------------------


# If you want to clip a worldclim environmental raster
# (https://www.worldclim.org/), which is in WGS84

# Reproject your polygon to the same projection of your environmental data
polygon_wgs <- spTransform(buffer_mpc, crs.wgs84)
plot(polygon_wgs)

# Read your present AND FUTURE environmental raster selected by correlation
pres_files <- list.files("./Maps/Present", full.names = T, 'tif$|bil$')
head(pres_files)

envi <- stack(pres_files)

# Cut your study area for the same extention and shape of your mpc
present_ly <- crop(envi, polygon_wgs)
present_ly2 <- mask(present_ly, polygon_wgs)

# Plot the results
#plot(present_ly2[[1]])
#plot(occurrence_records, add = T)

fut_files <- list.files("./Maps/Future/rcp45", full.names = T, 'tif$|bil$')
head(fut_files)

envi_fut <- stack(fut_files)

# Cut your study area for the same extention and shape of your mpc
future_ly <- crop(envi_fut, polygon_wgs)
future_ly2 <- mask(future_ly, polygon_wgs)

# Plot the results
#plot(future_ly2[[1]])
#plot(occurrence_records, add = T)



# Save your layers --------------------------------------------------------


dir.create(paste0("./outputs/", sp_names[a])) ##Fazer o loop aqui



writeOGR(
  polygon_wgs,
  dsn = paste0("./outputs/", sp_names[a]), ##SPECIES deve ser o nome da especie
  layer = "polygon_mpc_wgs",
  driver = "ESRI Shapefile",
  overwrite = T)
## TA DANDO ERRO PRA SALVAR
#Error in writeOGR(polygon_wgs, dsn = "./outputs",  :
#obj must be a SpatialPointsDataFrame, SpatialLinesDataFrame or
#SpatialPolygonsDataFrame

dir.create(paste0("./outputs/", sp_names[a], "/Pres_env_crop/"))

#writeRaster(present_ly2,
 #           "./outputs/SPECIES/Pres_env_crop/env_crop", ##SPECIES deve ser o nome da especie
  #          format = "GTiff",
   #         overwrite = T)

writeRaster(present_ly2, filename=paste0("./outputs/", sp_names[a], "/Pres_env_crop/", names(present_ly2)), bylayer=TRUE, format="GTiff")

dir.create(paste0("./outputs/",sp_names[a], "/Fut_env_crop/"))
writeRaster(future_ly2, filename=paste0("./outputs/", sp_names[a], "/Fut_env_crop/", names(future_ly2)), bylayer=TRUE, format="GTiff")

#}
