# Credits -----------------------------------------------------------------

# Code created by Luara Tourinho (https://github.com/luaratourinho)
# Edited by
# Julia Niemeyer
# Date: 08 august 2021


# First install packages
#install.packages("beepr")
# Then run the library required

library(maps)
library(mapdata)
library(rgdal)
library(raster)
library(rgeos)
library(dplyr)
library(beepr)


# Some examples of projections --------------------------------------------

# Creating two objects: (1) with the original projection of your data
# and (2) with an equal area projection
dir.create("./outputs")

file = "./data/03_clean_df_thin_1.csv" ##enter the name of your table

##minimum occurrence records to run analysis
n_min <- 15

## Enter name of folders to read environmental layers
#Present
pres_folder = "./Maps/Present/30s"
pres_files <- list.files(pres_folder, full.names = T, 'tif$|bil$') #don't change this
##standardize names of the variables as in your model [in the order of appearance in head(pres_files)]
head(pres_files)
names_var <- c('bio_15', 'bio_18', 'bio_4', 'bio_5')

##First RCP (name of the folder where your environmental layers are)
RCP1 <- "rcp_45" ##change number according to RCP

##Second RCP
RCP2 <- "rcp_85"


########## END OF ACTION  NEEDED ############

# Reading files -----------------------------------------------------------
sp <- read.table(file, header=TRUE, sep=",")#%>%
#select(species, lon, lat) %>%
#filter(species == "Hoplerythrinus_unitaeniatus")
sp_names <- unique(sp$species)
# running for one species
#sp.n = sp_names[[2]]

# WGS84
crs.wgs84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


# South America Albers Equal Area Conic
crs.albers <-
  CRS("+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0
                                 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

envi <- stack(pres_files)
names(envi) <- names_var
# All America
envi.cut <- raster::crop(envi, c(-160, -28, -60, 90))
#plot(envi.cut[[1]])

fut_files <- list.files(paste0("./Maps/Future/", RCP1), full.names = T, 'tif$|bil$')
#head(fut_files)
envi_fut <- stack(fut_files)
names(envi_fut) <- names_var ##standardize names of the variables as in your model
envi_fut_cut <- crop(envi_fut, c(-160, -28, -60, 90))
#fut_envi_res <- resample(envi_fut, envi.cut, method='bilinear')
#plot(envi_fut_cut[[1]])

fut_files2 <- list.files(paste0("./Maps/Future/", RCP2), full.names = T, 'tif$|bil$')
#head(fut_files2)
envi_fut2 <- stack(fut_files2)
names(envi_fut2) <- names_var  ##standardize names of the variables as in your model

envi_fut_cut2 <- crop(envi_fut2, c(-160, -28, -60, 90))

#length(sp_names) <- 370

for (a in 1:length(sp_names)){
  #
  message("starting the analysis for ", paste0(sp_names[a]))

  sp.n = sp_names[[a]]
  # Lambert Azimuthal Equal Area (ideal for Pacific Ocean islands)
  #crs.azim <-
  # CRS("+proj=laea +lat_0=90 +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84
  #               +datum=WGS84 +units=m +no_defs")



  # Building a minimum convex polygon ---------------------------------------


  # Read your table of occurrence records

  ##Planilha de ocorrência pós spthin aqui
  occurrence_records <- sp %>%
    dplyr::filter(species == paste0(sp_names[a])) %>%
    dplyr::select(species, lon, lat)

  # Number of occurrences to perform pseudoabsence sampling
  #Tive que add o plyr::count aqui pra não confundir com dplyr senão dá erro

  if (nrow(occurrence_records) < n_min){ ##Will not analyze species with less than n_min occurences
    print('species has less than 15 records and will not be analyzed')
    target_dir = paste0("./outputs/", sp_names[a], "/")
    dir.create( target_dir )
    write(format('species has less than 15 records and will not be analyzed'), file=paste(target_dir, 'STOP.txt', sep=""))
    next
  }

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
  #b <- (111000) ##Mudar valor pro valor que decidi
  #buffer_mpc <- gBuffer(mpc, width = b)

  # 20% of the polygon area is 2e-07
  b20 <- (gArea(mpc)*2e-07) ##20%
  b10 <- (gArea(mpc)*1e-07) ##10%
  b5 <- (gArea(mpc)*5e-8) ##5%
  b.list <- list(b20, b10, b5)

  for (i in 1:length(b.list)) {

    skip_to_next <<- FALSE

    # Note that print(b) fails since b doesn't exist
    buffer_mpc <- gBuffer(mpc, width = b.list[[i]])

    tryCatch(polygon_wgs <- spTransform(buffer_mpc, crs.wgs84),
             error = function(e) {skip_to_next <<- TRUE})

    if(skip_to_next == TRUE) {
      print(paste0('buffer of place ',i, ' got error, skipping to next'))
      next
    } else {

      print(paste0('using buffer of place ', i, ' in the list: ', b.list[[i]], ' m2'))
      # polygon_wgs <- spTransform(buffer_mpc, crs.wgs84)
      break
    }}


  # An example of mpc utility -------------------------------------------------


  # If you want to clip a worldclim environmental raster
  # (https://www.worldclim.org/), which is in WGS84

  # Reproject your polygon to the same projection of your environmental data
  polygon_wgs <- spTransform(buffer_mpc, crs.wgs84)
  #plot(polygon_wgs)

  # Cut your study area for the same extention and shape of your mpc
  present_ly <- crop(envi.cut, polygon_wgs)
  present_ly2 <- mask(present_ly, polygon_wgs)

  # Plot the results
  #plot(present_ly2[[4]])
  #plot(occurrence_records, add = T)


  # Cut your study area for the same extention and shape of your mpc
  ##First RCP
  future_ly <- crop(envi_fut_cut, polygon_wgs)
  future_ly2 <- mask(future_ly, polygon_wgs)
  #future_ly2[[3]] <- future_ly2[[3]]/10 #Bio4 tem que ser dividida por 10 ou por 100?

  # Plot the results
  #plot(future_ly2[[4]])
  #plot(occurrence_records, add = T)

  ##Second RCP
  future_ly3 <- crop(envi_fut_cut2, polygon_wgs)
  future_ly4 <- mask(future_ly3, polygon_wgs)
  #future_ly2[[3]] <- future_ly2[[3]]/10 #Bio4 tem que ser dividida por 10 ou por 100?

  # Plot the results
  #plot(future_ly2[[4]])
  #plot(occurrence_records, add = T)


  # Save your layers --------------------------------------------------------


  dir.create(paste0("./outputs/", sp_names[a])) ##Fazer o loop aqui

  #writeOGR(
  # polygon_wgs,
  #dsn = paste0("./outputs/", sp_names[a]), ##SPECIES deve ser o nome da especie
  #layer = "polygon_mpc_wgs",
  #driver = "ESRI Shapefile",
  #overwrite = T)
  ## TA DANDO ERRO PRA SALVAR
  #Error in writeOGR(polygon_wgs, dsn = "./outputs",  :
  #obj must be a SpatialPointsDataFrame, SpatialLinesDataFrame or
  #SpatialPolygonsDataFrame

  dir.create(paste0("./outputs/", sp_names[a], "/Pres_env_crop/"))

  #writeRaster(present_ly2,
  #           "./outputs/SPECIES/Pres_env_crop/env_crop", ##SPECIES deve ser o nome da especie
  #          format = "GTiff",
  #         overwrite = T)
  #names(present_ly2) <- c('bio_15', 'bio_18', 'bio_4', 'bio_5')
  writeRaster(present_ly2, filename=paste0("./outputs/", sp_names[a], "/Pres_env_crop/", names(present_ly2)), bylayer=TRUE, format="GTiff")

  dir.create(paste0("./outputs/",sp_names[a], "/Fut_env_crop/"))
  dir.create(paste0("./outputs/",sp_names[a], "/Fut_env_crop/", RCP1 ,"/"))
  #names(future_ly2) <- c(bio_1, bio_2, bio_3, bio_4) ##name the rasters as in your model
  writeRaster(future_ly2, filename=paste0("./outputs/", sp_names[a], "/Fut_env_crop/", RCP1 ,"/", names(future_ly2)), bylayer=TRUE, format="GTiff")
  dir.create(paste0("./outputs/",sp_names[a], "/Fut_env_crop/", RCP2 ,"/"))
  #names(future_ly4) <- c(bio_1, bio_2, bio_3, bio_4) ##name the rasters as in your model
  writeRaster(future_ly4, filename=paste0("./outputs/", sp_names[a], "/Fut_env_crop/", RCP2 ,"/", names(future_ly4)), bylayer=TRUE, format="GTiff")

  rm(future_ly)
  rm(future_ly2)
  rm(future_ly3)
  rm(future_ly4)
  rm(present_ly)
  rm(present_ly2)
gc()
}

beep(5) ##R will play a tune when this analysis is done
