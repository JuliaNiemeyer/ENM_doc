# Credits ---------------------------

# Original routine by
# Vitor Cavalcante


# Edited by
# Julia Niemeyer & Luara Tourinho

# Date: 04 Aug 2021

# This script is an example for generating pseudoabsence using biomod2 ----


# Required packages

library(biomod2)
library(raster)
library(dplyr)
library(beepr)

#.rs.restartR()

intersect_mask <- function(x){
  values_x <- getValues(x)
  inter_x <- values_x %*% rep(1,nlayers(x))
  mask <- setValues(subset(x,1),values = (inter_x>0))
  return(mask)
}

# Opening occurences ------------------------------------------------------

## Entrar com a planilha limpa pós spthin
file = "./data/03_clean_df_thin_1.csv" ##enter the name of your table

##minimum occurrence records to run analysis
n_min <- 15

########## END OF ACTION  NEEDED ############

# Reading files -----------------------------------------------------------
sp <- read.table(file, header=TRUE, sep=",")#%>%
sp_names <- unique(sp$species)


#length(sp_names) <- 380

for (a in 381:length(sp_names)){
  # message("starting the analysis for ", paste0(sp_names[a]))
  sp <- read.table(file, header=TRUE, sep=",") %>%
    filter(species == paste0(sp_names[a])) %>%
    select(species, lon, lat)

  message("starting the analysis for ", paste0(sp_names[a]))

  if (nrow(sp) < n_min){ ##Will not analyze species with less than 15 occurences
    print('species has less than 15 records and will not be analyzed')
    next
  }

  #nsp <- unique(sp$species)
  #species <- "Acestrorhynchus_britskii"

  #Get only lat and long
  My_target_species <- sp[,2:3]



  # Reading environmental variables -----------------------------------------

  # Read your environmental raster selected by correlation
  raster_files <- list.files(paste0("./outputs/", sp_names[a], "/Pres_env_crop"), full.names = T, 'tif$|bil$')
  #raster_files <- list.files('./Maps/Present', full.names = T, 'tif$|bil$')
  head(raster_files)

  environment <- stack(raster_files)

  ## keep only all cells that are defined for all layers
  environment <- stack(mask(environment, intersect_mask(environment)))

  #names(environment) <- c("Bio15","Bio18","Bio4","Bio5") ##nome das biovariáveis na ordem

  occurrence.resp <- rep(1, length(My_target_species$lon))

  skip_to_next <<- FALSE
  PA.list <- list(50000, 10000) ##it will run for 50km first. For species with restricted distribution, 50km will
  #return an error and the function will run for 10km

  for (i in 1:length(PA.list)) {

    tryCatch(Mymodel <- BIOMOD_FormatingData(
      resp.var = occurrence.resp,
      expl.var = environment,
      resp.xy = My_target_species,
      resp.name = "Occurrence",
      PA.nb.rep = 1,
      PA.nb.absences = length(sp$species),
      PA.strategy = "disk",
      PA.dist.min = PA.list[[i]],
      PA.dist.max = 20000000,
      na.rm = TRUE) , error = function(e) {skip_to_next <<- TRUE})

    if(skip_to_next == TRUE) {
      print(paste0('Species has restricted distribution, trying 10km'))
      next
    } else {

      print(paste0('using minimum distance of ', PA.list[[i]]))
      break
    }}
  #Mymodel

  gc()
  ##Another model
  for (i in 1:length(PA.list)) {

    tryCatch(Mymodel2 <- BIOMOD_FormatingData(
      resp.var = occurrence.resp,
      expl.var = environment,
      resp.xy = My_target_species,
      resp.name = "Occurrence",
      PA.nb.rep = 1,
      PA.nb.absences = length(sp$species)*10, ##Isso tem que entrar no loop
      PA.strategy = "disk",
      PA.dist.min = PA.list[[i]],
      PA.dist.max = 20000000,
      na.rm = TRUE) , error = function(e) {skip_to_next <<- TRUE})

    if(skip_to_next == TRUE) {
      print(paste0('Species has restricted distribution, trying 10km'))
      next
    } else {

      print(paste0('using minimum distance of ', PA.list[[i]]))
      break
    }}
  #Mymodel2



  ## function to get PA dataset
  get_PAtab <- function(bfd){
    dplyr::bind_cols(
      x = bfd@coord[, 1],
      y = bfd@coord[, 2],
      status = bfd@data.species,
      bfd@PA
    )
  }

  ## function to get background mask
  get_mask <- function(bfd){
    bfd@data.mask
  }

  (pres.xy <- get_PAtab(Mymodel) %>%
      filter(status == 1) %>%
      select(x, y))

  (pres.xy2 <- get_PAtab(Mymodel2) %>%
      filter(status == 1) %>%
      select(x, y))


  ## get the coordinates of pseudoabsences
  ## all repetition of pseudoabsences sampling merged
  (pa.all.xy <- get_PAtab(Mymodel) %>%
      filter(is.na(status)) %>%
      select(x, y)) %>%
    distinct()

  (pa.all.xy2 <- get_PAtab(Mymodel2) %>%
      filter(is.na(status)) %>%
      select(x, y)) %>%
    distinct()

  ##isso aqui não farei
  #plot(environment[[1]])
  #points(pa.all.xy, pch = 18)
  #points(pres.xy, pch = 20, col= "red")


  write.csv(pa.all.xy, paste0("./outputs/",sp_names[a], "/pseudoabs1.csv"), row.names = F) ##no loop colocar dentro da pasta de cada espécie
  write.csv(pa.all.xy2, paste0("./outputs/",sp_names[a], "/pseudoabs2.csv"), row.names = F) ##no loop colocar dentro da pasta de cada espécie


  ##esse aqui não farei
  pseudoabs <- pa.all.xy
  #head(pseudoabs)
  #dim(pseudoabs)

  pres = sp
  # Replace using your species name in "Genus_epithet" ##aqui ajeitar para loop
  pres$`species` <- sub(pattern = paste0(sp_names[a]), replacement = "1", x = pres$`species`)
  #tail(pres)
  pseudo_0 <- rep(0,nrow(pseudoabs))
  #pseudoabs$species <- NA
  #pseudoabs$species <- 0
  pseudoabs$species <- pseudo_0
  pseudoabs <- pseudoabs[,c(3,1,2)]
  names(pseudoabs) <-c("species","lon","lat")
  pres_pseudo_table <- rbind(pres,pseudoabs)
  #head(pres_pseudo_table)
  #tail(pres_pseudo_table)
  #dim(pres_pseudo_table)
  names(pres_pseudo_table) <-c("pa","lon","lat")



  pseudoabs2 <- pa.all.xy2
  #head(pseudoabs2)
  #dim(pseudoabs)

  pres = sp
  # Replace using your species name in "Genus_epithet" ##aqui ajeitar para loop
  pres$`species` <- sub(pattern = paste0(sp_names[a]), replacement = "1", x = pres$`species`)
  #tail(pres)
  pseudo_0 <- rep(0,nrow(pseudoabs2))
  pseudoabs2$species <- pseudo_0
  #pseudoabs2$species <- NA
  #pseudoabs2$species <- 0
  pseudoabs2 <- pseudoabs2[,c(3,1,2)]
  names(pseudoabs2) <-c("species","lon","lat")
  pres_pseudo_table2 <- rbind(pres,pseudoabs2)
  #head(pres_pseudo_table2)
  #tail(pres_pseudo_table)
  #dim(pres_pseudo_table)
  names(pres_pseudo_table2) <-c("pa","lon","lat")


  ##aqui mudar no loop para entrar dentro da pasta da espécie
  write.csv(pres_pseudo_table,paste0("./outputs/", sp_names[a],"/pres_pseudoabs.csv"), row.names = F)
  write.csv(pres_pseudo_table2,paste0("./outputs/", sp_names[a],"/pres_pseudoabs2.csv"), row.names = F)
  #pres_pseudo_table <- read.csv("./data/output/pres_pseudoabs.csv")
  #head(pres_pseudo_table)
gc()
}

beep(5)

############# Retirar
# Check pseudo abs points to see if they fall inside area
#coordinates(pres) <- c("lon", "lat")
# Define the original projection of your occurrence records
#proj4string(pres) <- crs.wgs84

#abs <- SpatialPoints(pres, crs.wgs84)
#plot(environment[[1]])
#plot(abs, add = T)
#
