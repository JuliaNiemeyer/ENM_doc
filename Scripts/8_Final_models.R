# Credits ---------------------------

# Script created by
# Bruno M. Carvalho (https://github.com/brunomc-eco)

# Edited by
# Julia Niemeyer & Luara Tourinho (https://github.com/luaratourinho)
# and others

# Date: 08 august 2021


# Required packages
library(parallel)
library(foreach)
library(doParallel)
require(raster)
library(rgeos)
library(rgdal)
library(dismo)
library(rJava)

library(kernlab)
library(randomForest)
library(maptools)


library(SDMTools)
library(dplyr)



indexOf <- function(v,findFor) {

  i=0
  for(i2 in v) {
    i = i + 1
    if (i2==findFor){
      return(i)
    }
  }
  return(0)
}


#setwd("/Mydirectory")



# Settings ----------------------------------------------------------------

# Name of .csv file with species occurrences, columns: 'species', 'lon', 'lat'
##Botar o caminho inteiro pra Ler a planilha do zero
file = "./data/03_clean_df_thin_1_BSF.csv"


# Reading files -----------------------------------------------------------

pres <- read.csv(file) # %>%
  #select(species, lon, lat) #%>%
  #filter(species == "Acestrorhynchus_britskii")
#pres <- as.data.frame(pres)

sp.names <- as.character(unique(pres$species))
# running for one species
sp.n = sp.names[[1]]

#Read predictor variables (i.e. present maps cropped by mcp)
raster_files <- list.files("./outputs/SPECIES/Pres_env_crop/", full.names = T, 'tif$|bil$')
head(raster_files)

predictors <- stack(raster_files)

# Read your future environmental rasters selected by correlation
raster_files2 <- list.files("./Maps/Future/rcp45", full.names = T, 'tif$|bil$')
head(raster_files2)

future_variable <- stack(raster_files2)
#names(future_variable)

# for GLM, RandomForest and SVM
model <- pa ~bio_1+bio_12+bio_5+bio_6
  # partitions (4 = 75% train e 25% test; 5 = 80% train e 20% test; 10 = 90% train e 10% test)
  k = 10
# number of background
bg.pt = 10000
# threshold ('spec_sens' = max. spec + sens)
t.met = 'spec_sens'
# Minimum TSS value for ensemble
tss.lim = 0.7

cont.maps = T # save continuous maps by algorithm
bin.maps = T # save binary maps by algorithm
ens.maps = T # save ensemble maps by algorithm

# ENM ---------------------------------------------------------------------


ini = Sys.time()
pb <- txtProgressBar(min = 1, max = length(sp.names)+1, style = 3)

# Number of occurrences to perform pseudoabsence sampling
#Tive que add o plyr::count aqui pra não confundir com dplyr senão dá erro
lim = plyr::count(pres$species)$freq[indexOf(plyr::count(pres$species)$x, sp.n)]

started_time = Sys.time()
cat( format( started_time, "%a %b %d %X %Y"), '-', 'STARTED', '\n')
cat( format( started_time, "%a %b %d %X %Y"), '-', 'Preparing train and test datasets for', sp.n, 'with ', lim, 'lines...', '\n')

##A pasta final onde tudo vai ser salvo (no loop é dentro da pasta de cada espécie_results talvez)
target_dir = paste("./outputs/SPECIES/results", '/', sep="" )
dir.create( target_dir )

if(file.exists(paste(target_dir, '/STARTED.txt', sep="")))
  stop("You MUST DELETE results folder before continue")

write(format( started_time, "%a %b %d %X %Y"), file=paste(target_dir, 'STARTED.txt', sep=""))

# For using different number of pseudoabsence in each algorithm, for example:
## pseudoausencia = n*10
sp.data <- read.csv(paste("./outputs/SPECIES/pres_pseudoabs2.csv"), header=TRUE, sep=',')

#colocar aqui embaixo  outra planilha depois de rodar pro outro modelo
##pseudoausencia = 100
sp.data2 <- read.csv(paste("./outputs/SPECIES/pres_pseudoabs.csv"), header=TRUE, sep=',')

#sp.data3 <- read.csv(paste('./spdata/', "pres_pseudoabs_10000", '.csv', sep=""), header=TRUE, sep=',')


# For bioclim and Maxent:
pres <- sp.data[sp.data$pa==1,2:3]
abs <- sp.data[sp.data$pa!=1,2:3]
bg <- randomPoints(predictors, bg.pt)
colnames(bg) <- c("lon", "lat")


## Crossvalidation partition in 10 groups (são 10 partições)
set.seed(10)
foldpres <- kfold(pres, 4)
set.seed(10)
foldabs <- kfold(abs, 4)

prestrain <- list()
prestest <- list()
abstrain <- list()
abstest <- list()
for(i in 1:k){
  foldpres <- kfold(pres, 4)
  foldabs <- kfold(abs, 4)
  prestrain[[i]] <- pres[foldpres != 1,]
  prestest[[i]] <- pres[foldpres == 1,]
  abstrain[[i]] <- abs[foldabs != 1,]
  abstest[[i]] <- abs[foldabs == 1,]
}


# For SVM:
train <- list()
pa_train <- list()
predtrain <- list()
testpres <- list()
testabs <- list()

for(i in 1:k){
  train[[i]] <- rbind(prestrain[[i]], abstrain[[i]])
  pa_train[[i]] <- c(rep(1, nrow(prestrain[[i]])), rep(0, nrow(abstrain[[i]])))
  predtrain[[i]] <- extract(predictors, train[[i]])
  predtrain[[i]] <- data.frame(cbind(pa=pa_train[[i]], predtrain[[i]]))
}


#For Random Forest
pres2 <- sp.data2[sp.data2$pa==1,2:3]
abs2 <- sp.data2[sp.data2$pa!=1,2:3]
set.seed(10)
foldpres2 <- kfold(pres2, 4)
set.seed(10)
foldabs2 <- kfold(abs2, 4)

prestrain2 <- list()
prestest2 <- list()
abstrain2 <- list()
abstest2 <- list()
for(i in 1:k){
  foldpres2 <- kfold(pres2, 4)
  foldabs2 <- kfold(abs2, 4)
  prestrain2[[i]] <- pres2[foldpres2 != 1,]
  prestest2[[i]] <- pres2[foldpres2 == 1,]
  abstrain2[[i]] <- abs2[foldabs2 != 1,]
  abstest2[[i]] <- abs2[foldabs2 == 1,]
}

train2 <- list()
pa_train2 <- list()
predtrain2 <- list()
testpres2 <- list()
testabs2 <- list()

for(i in 1:k){
  train2[[i]] <- rbind(prestrain2[[i]], abstrain2[[i]])
  pa_train2[[i]] <- c(rep(1, nrow(prestrain2[[i]])), rep(0, nrow(abstrain2[[i]])))
  predtrain2[[i]] <- extract(predictors, train2[[i]])
  predtrain2[[i]] <- data.frame(cbind(pa=pa_train2[[i]], predtrain2[[i]]))
}


# For GLM:
pres3 <- sp.data[sp.data$pa==1,2:3]
abs3 <- sp.data[sp.data$pa!=1,2:3]
set.seed(10)
foldpres3 <- kfold(pres3, 4)
set.seed(10)
foldabs3 <- kfold(abs3, 4)

prestrain3 <- list()
prestest3 <- list()
abstrain3 <- list()
abstest3 <- list()
for(i in 1:k){
  foldpres3 <- kfold(pres3, 4)
  foldabs3 <- kfold(abs3, 4)
  prestrain3[[i]] <- pres3[foldpres3 != 1,]
  prestest3[[i]] <- pres3[foldpres3 == 1,]
  abstrain3[[i]] <- abs3[foldabs3 != 1,]
  abstest3[[i]] <- abs3[foldabs3 == 1,]
}

train3 <- list()
pa_train3 <- list()
predtrain3 <- list()
testpres3 <- list()
testabs3 <- list()

for(i in 1:k){
  train3[[i]] <- rbind(prestrain3[[i]], abstrain3[[i]])
  pa_train3[[i]] <- c(rep(1, nrow(prestrain3[[i]])), rep(0, nrow(abstrain3[[i]])))
  predtrain3[[i]] <- extract(predictors, train3[[i]])
  predtrain3[[i]] <- data.frame(cbind(pa=pa_train3[[i]], predtrain3[[i]]))
}



# Running models ----------------------------------------------------------



cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Bioclim model for', sp.n, '...', '\n')
bc <- list()
evbc <- list()
bcTSS <- list()
bcAUC <- list()
bckappa <- list()
bcthres <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Bioclim (', i, ') model for', sp.n, '...', '\n')
  bc[[i]] <- bioclim(predictors, prestrain[[i]])
  evbc[[i]] <- evaluate(prestest[[i]], abstest[[i]], bc[[i]], predictors)
  bcTSS[[i]] <- max(evbc[[i]]@TPR + evbc[[i]]@TNR)-1
  bcAUC[[i]] <- evbc[[i]]@auc
  bckappa[[i]] <- max(evbc[[i]]@kappa)
  bcthres[[i]] <- threshold(evbc[[i]], t.met)
}

cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running GLM (logistic regression) for', sp.n, '...', '\n')
gm <- list()
evgm <- list()
gmTSS <- list()
gmAUC <- list()
gmkappa <- list()
gmthres <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running GLM (logistic regression) (', i, ') for', sp.n, '...', '\n')
  gm[[i]] <- glm(model, family=binomial(link="logit"), data=predtrain3[[i]])
  evgm[[i]] <- evaluate(prestest3[[i]], abstest3[[i]], gm[[i]], predictors)
  gmTSS[[i]] <- max(evgm[[i]]@TPR + evgm[[i]]@TNR)-1
  gmAUC[[i]] <- evgm[[i]]@auc
  gmkappa[[i]] <- max(evgm[[i]]@kappa)
  gmthres[[i]] <- threshold(evgm[[i]], t.met)
}

cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Random Forest model for', sp.n, '...', '\n')
rf <- list()
evrf <- list()
rfTSS <- list()
rfAUC <- list()
rfkappa <- list()
rfthres <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Random Forest (', i, ') model for', sp.n, '...', '\n')
  rf[[i]] <- randomForest(model, data=predtrain2[[i]], na.action=na.omit)
  evrf[[i]] <- evaluate(prestest2[[i]], abstest2[[i]], rf[[i]], predictors)
  rfTSS[[i]] <- max(evrf[[i]]@TPR + evrf[[i]]@TNR)-1
  rfAUC[[i]] <- evrf[[i]]@auc
  rfkappa[[i]] <- max(evrf[[i]]@kappa)
  rfthres[[i]] <- threshold(evrf[[i]], t.met)
}

cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Maxent model for', sp.n, '...', '\n')
mx <- list()
evmx <- list()
mxTSS <- list()
mxAUC <- list()
mxkappa <- list()
mxthres <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running Maxent (', i, ') model for', sp.n, '...', '\n')
  mx[[i]] <- maxent(predictors, prestrain[[i]], a=bg)
  evmx[[i]] <- evaluate(prestest[[i]], abstest[[i]], mx[[i]], predictors)
  mxTSS[[i]] <- max(evmx[[i]]@TPR + evmx[[i]]@TNR)-1
  mxAUC[[i]] <- evmx[[i]]@auc
  mxkappa[[i]] <- max(evmx[[i]]@kappa)
  mxthres[[i]] <- threshold(evmx[[i]], t.met)
}

cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running SVM model for', sp.n, '...', '\n')
sv <- list()
evsv <- list()
svTSS <- list()
svAUC <- list()
svkappa <- list()
svthres <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Running SVM (', i, ') model for', sp.n, '...', '\n')
  sv[[i]] <- ksvm(model, data=predtrain[[i]])
  evsv[[i]] <- evaluate(prestest[[i]], abstest[[i]], sv[[i]], predictors)
  svTSS[[i]] <- max(evsv[[i]]@TPR + evsv[[i]]@TNR)-1
  svAUC[[i]] <- evsv[[i]]@auc
  svkappa[[i]] <- max(evsv[[i]]@kappa)
  svthres[[i]] <- threshold(evsv[[i]], t.met)
}



# Validation table --------------------------------------------------------


cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Generating Validation table for models of', sp.n, '...', '\n')
bcTSSval <- unlist(bcTSS)
#bcTSSval <- data.frame(matrix(unlist(bcTSS), nrow=length(bcTSS), byrow=TRUE), stringsAsFactors=FALSE)

gmTSSval <- unlist(gmTSS)
#gmTSSval <- data.frame(matrix(unlist(gmTSS), nrow=length(gmTSS), byrow=TRUE), stringsAsFactors=FALSE)

rfTSSval <- unlist(rfTSS)
#rfTSSval <- data.frame(matrix(unlist(rfTSS), nrow=length(rfTSS), byrow=TRUE,), stringsAsFactors=FALSE)

mxTSSval <- unlist(mxTSS)
#mxTSSval <- data.frame(matrix(unlist(mxTSS), nrow=length(mxTSS), byrow=TRUE), stringsAsFactors=FALSE)


svTSSval <- unlist(svTSS)
bcAUCval <- unlist(bcAUC)
gmAUCval <- unlist(gmAUC)
rfAUCval <- unlist(rfAUC)
mxAUCval <- unlist(mxAUC)
svAUCval <- unlist(svAUC)
bckappaval <- unlist(bckappa)
gmkappaval <- unlist(gmkappa)
rfkappaval <- unlist(rfkappa)
mxkappaval <- unlist(mxkappa)
svkappaval <- unlist(svkappa)
mod.names <- c(rep ('bc', k), rep('gm', k), rep('rf', k), rep('mx', k), rep('sv', k))
mod.sp <- c(rep(sp.n, k*5)) #conferir a planilha

TSS <- c(bcTSSval, gmTSSval, rfTSSval, mxTSSval, svTSSval)
AUC <- c(bcAUCval, gmAUCval, rfAUCval, mxAUCval, svAUCval)
kappa <- c(bckappaval, gmkappaval, rfkappaval, mxkappaval, svkappaval)
Valid <- data.frame(mod.sp, mod.names, TSS, AUC, kappa, stringsAsFactors=FALSE)

##Pasta de resultados na pasta de cada espécie
write.csv(Valid, file = paste( target_dir, '/Valid_', sp.n, '.csv', sep=""))



# Projecting models -------------------------------------------------------



# Projecting Bioclim -------------------------------------------------------

#Drop bad models, project under current and future conditions
cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Bioclim model of', sp.n, '...', '\n')
cur.bc <- list()
cur.bc.bin <- list()
future_variable.bc <- list()
future_variable.bc.bin <- list()


for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Bioclim (', i, ') model of', sp.n, '...', '\n')
  if(bcTSS[[i]] >= tss.lim){
    cur.bc[[i]] <- predict(predictors, bc[[i]])
    cur.bc.bin[[i]] <- cur.bc[[i]] > bcthres[[i]]
    future_variable.bc[[i]] <- predict(future_variable, bc[[i]])
    future_variable.bc.bin[[i]] <- future_variable.bc[[i]] > bcthres[[i]]

  } else {
    cur.bc[[i]] <- NULL
    cur.bc.bin[[i]] <- NULL
    future_variable.bc[[i]] <- NULL
    future_variable.bc.bin[[i]] <- NULL
  }
}


##########
cur.bc <- Filter(Negate(is.null), cur.bc) # remove null rasters
cur.bc.bin <- Filter(Negate(is.null), cur.bc.bin)

#Ensemble of binary models (Majority Rule)
cur.bc.ens <- Reduce('+', cur.bc.bin) #sum of pixels
tval <- unique(cur.bc.ens)
tval <- tval[tval != 0]
tval <- median(tval)
cur.bc.ens.bin <- cur.bc.ens >= tval


##########
future_variable.bc <- Filter(Negate(is.null), future_variable.bc)
future_variable.bc.bin <- Filter(Negate(is.null), future_variable.bc.bin)

future_variable.bc.ens <- Reduce('+', future_variable.bc.bin)
tval <- unique(future_variable.bc.ens)
tval <- tval[tval != 0]
tval <- median(tval)
future_variable.bc.ens.bin <- future_variable.bc.ens >= tval

##########

# Taking out the zero values and normalizing rasters

for(z in 1:length(cur.bc)){
  adeq = cur.bc[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  cur.bc[[z]] <- calc(adeq, adeq_norm)
}



for(z in 1:length(future_variable.bc)){
  adeq = future_variable.bc[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  future_variable.bc[[z]] <- calc(adeq, adeq_norm)
}


#Ensemble of continuos models

x <- stack(cur.bc)
cur.bc.cont <- calc(x, fun = mean)

x <- stack(future_variable.bc)
future_variable.bc.cont <- calc(x, fun = mean)

##obs. isso aqui talve remover - ver oq faz mas talvez não precise (rodar e ver)
#for(i in 1:k){
 # min_max <- range(gcmmg.bc[[i]][], na.rm=T)
  #write.table(min_max, paste("gcmmg.bc", i, sp.n, ".txt",  sep="_"))
#}

#rm(cur.bc.bin, cur.bc.cont, cur.bc.ens)
#rm(future_variable.bc.bin, future_variable.bc.cont, future_variable.bc.ens)
###

# Projecting GLM ----------------------------------------------------------


cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting GLMs of', sp.n, '...', '\n')
cur.gm <- list()
cur.gm.bin <- list()
future_variable.gm <- list()
future_variable.gm.bin <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting GLMs (', i, ') of', sp.n, '...', '\n')
  if(gmTSS[[i]] >= tss.lim){
    cur.gm[[i]] <- predict(predictors, gm[[i]])
    cur.gm.bin[[i]] <- cur.gm[[i]] > gmthres[[i]]
    future_variable.gm[[i]] <- predict(future_variable, gm[[i]])
    future_variable.gm.bin[[i]] <- future_variable.gm[[i]] > gmthres[[i]]

  } else {
    cur.gm[[i]] <- NULL
    cur.gm.bin[[i]] <- NULL
    future_variable.gm[[i]] <- NULL
    future_variable.gm.bin[[i]] <- NULL

  }
}

##########
cur.gm <- Filter(Negate(is.null), cur.gm)
cur.gm.bin <- Filter(Negate(is.null), cur.gm.bin)

cur.gm.ens <- Reduce('+', cur.gm.bin)
tval <- unique(cur.gm.ens)
tval <- tval[tval != 0]
tval <- median(tval)
cur.gm.ens.bin <- cur.gm.ens >= tval

##########
future_variable.gm <- Filter(Negate(is.null), future_variable.gm)
future_variable.gm.bin <- Filter(Negate(is.null), future_variable.gm.bin)

future_variable.gm.ens <- Reduce('+', future_variable.gm.bin)
tval <- unique(future_variable.gm.ens)
tval <- tval[tval != 0]
tval <- median(tval)
future_variable.gm.ens.bin <- future_variable.gm.ens >= tval


##########

for(z in 1:length(cur.gm)){
  adeq = cur.gm[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  cur.gm[[z]] <- calc(adeq, adeq_norm)
}


for(z in 1:length(future_variable.gm)){
  adeq = future_variable.gm[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  future_variable.gm[[z]] <- calc(adeq, adeq_norm)
}


x <- stack(cur.gm)
cur.gm.cont <- calc(x, fun = mean)

x <- stack(future_variable.gm)
future_variable.gm.cont <- calc(x, fun = mean)




# Projecting RF ----------------------------------------------------------



cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Random Forest models of', sp.n, '...', '\n')
cur.rf <- list()
cur.rf.bin <- list()
future_variable.rf <- list()
future_variable.rf.bin <- list()


for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Random Forest (', i, ') models of', sp.n, '...', '\n')
  if(rfTSS[[i]] >= tss.lim){
    cur.rf[[i]] <- predict(predictors, rf[[i]])
    cur.rf.bin[[i]] <- cur.rf[[i]] > rfthres[[i]]

    future_variable.rf[[i]] <- predict(future_variable, rf[[i]])
    future_variable.rf.bin[[i]] <- future_variable.rf[[i]] > rfthres[[i]]
  } else {

    cur.rf[[i]] <- NULL
    cur.rf.bin[[i]] <- NULL
    future_variable.rf[[i]] <- NULL
    future_variable.rf.bin[[i]] <- NULL

  }
}

##########
cur.rf <- Filter(Negate(is.null), cur.rf)
cur.rf.bin <- Filter(Negate(is.null), cur.rf.bin)

cur.rf.ens <- Reduce('+', cur.rf.bin)
tval <- unique(cur.rf.ens)
tval <- tval[tval != 0]
tval <- median(tval)
cur.rf.ens.bin <- cur.rf.ens >= tval


##########
future_variable.rf <- Filter(Negate(is.null), future_variable.rf)
future_variable.rf.bin <- Filter(Negate(is.null), future_variable.rf.bin)

future_variable.rf.ens <- Reduce('+', future_variable.rf.bin)
tval <- unique(future_variable.rf.ens)
tval <- tval[tval != 0]
tval <- median(tval)
future_variable.rf.ens.bin <- future_variable.rf.ens >= tval


##########

for(z in 1:length(cur.rf)){
  adeq = cur.rf[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    cur.rf[[z]] <- calc(adeq, adeq_norm)
  }

}


for(z in 1:length(future_variable.rf)){
  adeq = future_variable.rf[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    future_variable.rf[[z]] <- calc(adeq, adeq_norm)
  }
}

cur.rf <- Filter(Negate(is.null), cur.rf)
future_variable.rf <- Filter(Negate(is.null), future_variable.rf)

x <- stack(cur.rf)
cur.rf.cont <- calc(x, fun = mean)

x <- stack(future_variable.rf)
future_variable.rf.cont <- calc(x, fun = mean)



# Projecting Maxent --------------------------------------------------------



cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Maxent models of', sp.n, '...', '\n')
cur.mx <- list()
cur.mx.bin <- list()
future_variable.mx <- list()
future_variable.mx.bin <- list()

for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting Maxent (', i, ') models of', sp.n, '...', '\n')
  if(mxTSS[[i]] >= tss.lim){
    cur.mx[[i]] <- predict(predictors, mx[[i]])
    cur.mx.bin[[i]] <- cur.mx[[i]] > mxthres[[i]]

    future_variable.mx[[i]] <- predict(future_variable, mx[[i]])
    future_variable.mx.bin[[i]] <- future_variable.mx[[i]] > mxthres[[i]]

  } else {
    cur.mx[[i]] <- NULL
    cur.mx.bin[[i]] <- NULL
    future_variable.mx[[i]] <- NULL
    future_variable.mx.bin[[i]] <- NULL

  }
}


##########
cur.mx <- Filter(Negate(is.null), cur.mx)
cur.mx.bin <- Filter(Negate(is.null), cur.mx.bin)

cur.mx.ens <- Reduce('+', cur.mx.bin)
tval <- unique(cur.mx.ens)
tval <- tval[tval != 0]
tval <- median(tval)
cur.mx.ens.bin <- cur.mx.ens >= tval


##########
future_variable.mx <- Filter(Negate(is.null), future_variable.mx)
future_variable.mx.bin <- Filter(Negate(is.null), future_variable.mx.bin)

future_variable.mx.ens <- Reduce('+', future_variable.mx.bin)
tval <- unique(future_variable.mx.ens)
tval <- tval[tval != 0]
tval <- median(tval)
future_variable.mx.ens.bin <- future_variable.mx.ens >= tval


##########


for(z in 1:length(cur.mx)){
  adeq = cur.mx[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    cur.mx[[z]] <- calc(adeq, adeq_norm)
  }
}

for(z in 1:length(future_variable.mx)){
  adeq = future_variable.mx[[z]]
  if(sum(adeq[], na.rm=T)!=0){
    minimo <- min(adeq[], na.rm=T)
    maximo <- max(adeq[], na.rm=T)
    adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
    future_variable.mx[[z]] <- calc(adeq, adeq_norm)
  }
}

x <- stack(cur.mx)
cur.mx.cont <- calc(x, fun = mean)

x <- stack(future_variable.mx)
future_variable.mx.cont <- calc(x, fun = mean)



# Projecting SVM ----------------------------------------------------------



cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting SVM models of', sp.n, '...', '\n')
cur.sv <- list()
cur.sv.bin <- list()
future_variable.sv <- list()
future_variable.sv.bin <- list()


for(i in 1:k){
  cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting SVM (', i, ') models of', sp.n, '...', '\n')
  if(svTSS[[i]] >= tss.lim){
    cur.sv[[i]] <- predict(predictors, sv[[i]])
    cur.sv.bin[[i]] <- cur.sv[[i]] > svthres[[i]]

    future_variable.sv[[i]] <- predict(future_variable, sv[[i]])
    future_variable.sv.bin[[i]] <- future_variable.sv[[i]] > svthres[[i]]

  } else {
    cur.sv[[i]] <- NULL
    cur.sv.bin[[i]] <- NULL
    future_variable.sv[[i]] <- NULL
    future_variable.sv.bin[[i]] <- NULL

  }
}


##########
cur.sv <- Filter(Negate(is.null), cur.sv)
cur.sv.bin <- Filter(Negate(is.null), cur.sv.bin)

cur.sv.ens <- Reduce('+', cur.sv.bin)
tval <- unique(cur.sv.ens)
tval <- tval[tval != 0]
tval <- median(tval)
cur.sv.ens.bin <- cur.sv.ens >= tval


##########
future_variable.sv <- Filter(Negate(is.null), future_variable.sv)
future_variable.sv.bin <- Filter(Negate(is.null), future_variable.sv.bin)

future_variable.sv.ens <- Reduce('+', future_variable.sv.bin)
tval <- unique(future_variable.sv.ens)
tval <- tval[tval != 0]
tval <- median(tval)
future_variable.sv.ens.bin <- future_variable.sv.ens >= tval


##########

for(z in 1:length(cur.sv)){
  adeq = cur.sv[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  cur.sv[[z]] <- calc(adeq, adeq_norm)
}


for(z in 1:length(future_variable.sv)){
  adeq = future_variable.sv[[z]]
  minimo <- min(adeq[], na.rm=T)
  maximo <- max(adeq[], na.rm=T)
  adeq_norm <- function(x) {(x-minimo)/(maximo-minimo)}
  future_variable.sv[[z]] <- calc(adeq, adeq_norm)
}


x <- stack(cur.sv)
cur.sv.cont <- calc(x, fun = mean)

x <- stack(future_variable.sv)
future_variable.sv.cont <- calc(x, fun = mean)




###################### ensemble current ALL ALGO
# binarios
cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Projecting ensemble models of', sp.n, '...', '\n')

##aqui ver de add um ifelse caso algum algoritmo tenha dado abaixo do TSS, que vai dar erro - ignorar se dar zero (NULL)
#algo.cur <- list(cur.bc.ens.bin, cur.gm.ens.bin, cur.rf.ens.bin, cur.mx.ens.bin, cur.sv.ens.bin)

list <- c(cur.gm.ens.bin, cur.rf.ens.bin, cur.mx.ens.bin, cur.sv.ens.bin)
algo.cur <- list[is.na(list) == FALSE]
#algo.cur <- list(list)
#algo.cur <- list(cur.gm.ens.bin, cur.rf.ens.bin, cur.mx.ens.bin, cur.sv.ens.bin) ##pegar apenas os que não deram NULL (colocar no ifelse)
algo.cur2 <- algo.cur[is.na(algo.cur) == FALSE]
ens.cur <- Reduce('+', algo.cur2)
tval <- unique(ens.cur)
tval <- tval[tval != 0]
tval <- median(tval)
ens.cur.bin <- ens.cur >= tval  ##ensemble presente binario


# Ensemble ----------------------------------------------------------------

# For continuous model, by TSS-weighted average

##Mesmo ifelse de acima
w <- c(bcTSSval, gmTSSval, rfTSSval, mxTSSval, svTSSval) ##ja ignora o que for null
#w_fut <- w[is.null(w)] ##adicionei isso pra retirar qualquer um que for NULL


if (is.null(cur.bc == FALSE)) {
st1 <- stack(cur.bc)
for(i in 1:length(cur.bc)){
  min_max <- range(cur.bc[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.bc", i, sp.n, ".txt",  sep="_"))
}
#rm(cur.bc)
} else {
  st1 <- NULL
  }

if (is.null(cur.gm) == FALSE){
st2 <- stack(cur.gm)
for(i in 1:length(cur.gm)){
  min_max <- range(cur.gm[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.gm", i, sp.n, ".txt",  sep="_"))
}
#rm(cur.gm)
} else {
  st2 <- NULL
}

if (is.null(cur.rf) == FALSE){
st3 <- stack(cur.rf)
for(i in 1:length(cur.rf)){
  min_max <- range(cur.rf[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.rf", i, sp.n, ".txt",  sep="_"))
}
} else {
  st3 <- NULL
}

if (is.null(cur.mx) == FALSE){
st4 <- stack(cur.mx)
for(i in 1:length(cur.mx)){
  min_max <- range(cur.mx[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.mx", i, sp.n, ".txt",  sep="_"))
}
#rm(cur.mx)
} else {
  st4 <- NULL
}

if (is.null(cur.sv) == FALSE){
st5 <- stack(cur.sv)
for(i in 1:length(cur.sv)){
  min_max <- range(cur.sv[[i]][], na.rm=T)
  #write.table(min_max, paste("cur.sv", i, sp.n, ".txt",  sep="_"))
}
#rm(cur.sv)
} else {
  st5 <- NULL
}

stacks <- c(st1, st2, st3, st4, st5) ## ignora o que for null
st <- stack(stacks)

w.sem.os.nulos <- w[w >= tss.lim]
ens2.cur <- weighted.mean(st, w.sem.os.nulos) ##ensamble present continue --resultado final
ens2.cur.sd.w <- sum(w.sem.os.nulos * (st - ens2.cur)^2) ## variancia
ens2.cur.sd.w <- sqrt(ens2.cur.sd.w) ##desvio padrão (incerteza) -- pra caso revisor peça



# Uncertainty for each algo (each partition) -------------------------------------------------------------

wbc <- c(bcTSSval)
wbc.sem.os.nulos <- wbc[wbc >= tss.lim]
if(is.null(wbc.sem.os.nulos)){
 print('bioclim is null') } else {
cur.bc.mean.w <- weighted.mean(st2, wbc.sem.os.nulos)
cur.bc.sd.w <- sum(wbc.sem.os.nulos * (st2 - cur.bc.mean.w)^2)
cur.bc.sd.w <- sqrt(cur.bc.sd.w) }


wgm <- c(gmTSSval)
wgm.sem.os.nulos <- wgm[wgm >= tss.lim]
if(is.null(wgm.sem.os.nulos)){
  print('glm is null') } else {
cur.gm.mean.w <- weighted.mean(st2, wgm.sem.os.nulos)
cur.gm.sd.w <- sum(wgm.sem.os.nulos * (st2 - cur.gm.mean.w)^2)
cur.gm.sd.w <- sqrt(cur.gm.sd.w)
  }

wrf <- c(rfTSSval)
wrf.sem.os.nulos <- wrf[wrf >= tss.lim]
if (is.null(wrf.sem.os.nulos)){
  print('rf is null') } else {
cur.rf.mean.w <- weighted.mean(st3, wrf.sem.os.nulos)
cur.rf.sd.w <- sum(wrf.sem.os.nulos * (st3 - cur.rf.mean.w)^2)
cur.rf.sd.w <- sqrt(cur.rf.sd.w)
  }

wmx <- c(mxTSSval)
wmx.sem.os.nulos <- wmx[wmx >= tss.lim]
if (is.null(wmx.sem.os.nulos)) {
  print('maxent is null')
} else {
cur.mx.mean.w <- weighted.mean(st4, wmx.sem.os.nulos)
cur.mx.sd.w <- sum(wmx.sem.os.nulos * (st4 - cur.mx.mean.w)^2)
cur.mx.sd.w <- sqrt(cur.mx.sd.w)
}

wsv <- c(svTSSval)
wsv.sem.os.nulos <- wsv[wsv >= tss.lim]
if (is.null(wsv.sem.os.nulos)){
  print('svmk is null')
} else {
cur.sv.mean.w <- weighted.mean(st5, wsv.sem.os.nulos)
cur.sv.sd.w <- sum(wsv.sem.os.nulos * (st5 - cur.sv.mean.w)^2)
cur.sv.sd.w <- sqrt(cur.sv.sd.w)
}

# Ensemble (future) -------------------------------------------------------------

# For binary model

# binarios
list_fut <- c(future_variable.bc.ens.bin, future_variable.gm.ens.bin, future_variable.rf.ens.bin, future_variable.mx.ens.bin, future_variable.sv.ens.bin)
algo.future_variable <- list[is.na(list) == FALSE] ##will ignore null rasters
algo.future_variable <- algo.future_variable[is.na(algo.future_variable) == FALSE]
ens.future_variable <- Reduce('+', algo.future_variable)
tval <- unique(ens.future_variable)
tval <- tval[tval != 0]
tval <- median(tval)
ens.future_variable.bin <- ens.future_variable >= tval

# continuous
if (is.null(future_variable.bc)){
  st1 <- NULL
}else {
  st1 <- stack(future_variable.bc)
  for(i in 1:length(future_variable.bc)){
    min_max <- range(future_variable.bc[[i]][], na.rm=T)
    #write.table(min_max, paste("future_variable.bc", i, sp.n, ".txt",  sep="_"))
  }
  }
#rm(future_variable.bc)

if (is.null(future_variable.gm)){
  st2 <- NULL }else {
    st2 <- stack(future_variable.gm)
    for(i in 1:length(future_variable.gm)){
      min_max <- range(future_variable.gm[[i]][], na.rm=T)
      #write.table(min_max, paste("future_variable.gm", i, sp.n, ".txt",  sep="_"))
    }
}
#rm(future_variable.gm)

if (is.null(future_variable.rf)){
  st3 <- NULL
}else {
  st3 <- stack(future_variable.rf)
  for(i in 1:length(future_variable.rf)){
    min_max <- range(future_variable.rf[[i]][], na.rm=T)
    #write.table(min_max, paste("future_variable.rf", i, sp.n, ".txt",  sep="_"))
  }
}
#rm(future_variable.rf)

if (is.null(future_variable.mx)){
  st4 <- NULL
}else {
  st4 <- stack(future_variable.mx)
  for(i in 1:length(future_variable.mx)){
    min_max <- range(future_variable.mx[[i]][], na.rm=T)
    #write.table(min_max, paste("future_variable.mx", i, sp.n, ".txt",  sep="_"))
  }
}
#rm(future_variable.mx)

if (is.na(future_variable.sv)){
  st5 <- NULL
} else {
  st5 <- stack(future_variable.sv)
  for(i in 1:length(future_variable.sv)){
    min_max <- range(future_variable.sv[[i]][], na.rm=T)
    #write.table(min_max, paste("future_variable.sv", i, sp.n, ".txt",  sep="_"))
  }
}
#rm(future_variable.sv)

stacks_fut <- c(st1, st2, st3, st4, st5) ## ignora o que for null
st_fut <- stack(stacks_fut)

w.sem.os.nulos <- w[w >= tss.lim]
ens2.future_variable <- weighted.mean(st_fut, w.sem.os.nulos)
ens.fut.cont <- ens2.future_variable
ens2.fut.sd.w <- sum(w.sem.os.nulos * (st_fut - ens.fut.cont)^2)
ens2.fut.sd.w <- sqrt(ens2.fut.sd.w)


# Uncertainty future ------------------------------------------------------

wbc <- c(bcTSSval)
wbc.sem.os.nulos <- wbc[wbc >= tss.lim]
future_variable.bc.mean.w <- weighted.mean(st2, wbc.sem.os.nulos)
future_variable.bc.sd.w <- sum(wbc.sem.os.nulos * (st2 - future_variable.bc.mean.w)^2)
future_variable.bc.sd.w <- sqrt(future_variable.bc.sd.w)

wgm <- c(gmTSSval)
wgm.sem.os.nulos <- wgm[wgm >= tss.lim]
future_variable.gm.mean.w <- weighted.mean(st2, wgm.sem.os.nulos)
future_variable.gm.sd.w <- sum(wgm.sem.os.nulos * (st2 - future_variable.gm.mean.w)^2)
future_variable.gm.sd.w <- sqrt(future_variable.gm.sd.w)

wrf <- c(rfTSSval)
wrf.sem.os.nulos <- wrf[wrf >= tss.lim]
future_variable.rf.mean.w <- weighted.mean(st3, wrf.sem.os.nulos)
future_variable.rf.sd.w <- sum(wrf.sem.os.nulos * (st3 - future_variable.rf.mean.w)^2)
future_variable.rf.sd.w <- sqrt(future_variable.rf.sd.w)

wmx <- c(mxTSSval)
wmx.sem.os.nulos <- wmx[wmx >= tss.lim]
future_variable.mx.mean.w <- weighted.mean(st4, wmx.sem.os.nulos)
future_variable.mx.sd.w <- sum(wmx.sem.os.nulos * (st4 - future_variable.mx.mean.w)^2)
future_variable.mx.sd.w <- sqrt(future_variable.mx.sd.w)

wsv <- c(svTSSval)
wsv.sem.os.nulos <- wsv[wsv >= tss.lim]
future_variable.sv.mean.w <- weighted.mean(st5, wsv.sem.os.nulos)
future_variable.sv.sd.w <- sum(wsv.sem.os.nulos * (st5 - future_variable.sv.mean.w)^2)
future_variable.sv.sd.w <- sqrt(future_variable.sv.sd.w)

####

ens.fut <- ens.future_variable.bin
tval <- unique(ens.fut)
tval <- tval[tval != 0]
tval <- median(tval)
ens.fut.bin <- ens.fut >= tval


cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Saving ensemble maps of', sp.n, '...', '\n')
writeRaster(ens2.cur, file = paste(target_dir, '/CUR.cont_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(ens2.future_variable, file = paste(target_dir, '/MPI8570.cont_', sp.n, '.asc', sep=""),overwrite=TRUE)

writeRaster(ens.cur.bin, file = paste(target_dir, '/CUR.bin_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(ens.future_variable.bin, file = paste(target_dir, '/MPI8570.bin_', sp.n, '.asc', sep=""),overwrite=TRUE)

#weighted.mean and weighted.sd
writeRaster(ens2.cur.sd.w, file = paste(target_dir, '/ens.cur.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(ens2.fut.sd.w, file = paste(target_dir, '/ens.mpi8570.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.bc.mean.w, file = paste(target_dir, '/cur.bc.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.bc.sd.w, file = paste(target_dir, '/cur.bc.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.gm.mean.w, file = paste(target_dir, '/cur.gm.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(cur.gm.sd.w, file = paste(target_dir, '/cur.gm.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.rf.mean.w, file = paste(target_dir, '/cur.rf.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(cur.rf.sd.w, file = paste(target_dir, '/cur.rf.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.mx.mean.w, file = paste(target_dir, '/cur.mx.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(cur.mx.sd.w, file = paste(target_dir, '/cur.mx.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(cur.sv.mean.w, file = paste(target_dir, '/cur.sv.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(cur.sv.sd.w, file = paste(target_dir, '/cur.sv.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.bc.mean.w, file = paste(target_dir, '/mpi8570.bc.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.bc.sd.w, file = paste(target_dir, '/mpi8570.bc.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.gm.mean.w, file = paste(target_dir, '/mpi8570.gm.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(future_variable.gm.sd.w, file = paste(target_dir, '/mpi8570.gm.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.rf.mean.w, file = paste(target_dir, '/mpi8570.rf.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(future_variable.rf.sd.w, file = paste(target_dir, '/mpi8570.rf.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.mx.mean.w, file = paste(target_dir, '/mpi8570.mx.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(future_variable.mx.sd.w, file = paste(target_dir, '/mpi8570.mx.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
#writeRaster(future_variable.sv.mean.w, file = paste(target_dir, '/mpi8570.sv.mean.w_', sp.n, '.asc', sep=""),overwrite=TRUE)
writeRaster(future_variable.sv.sd.w, file = paste(target_dir, '/mpi8570.sv.sd.w_', sp.n, '.asc', sep=""),overwrite=TRUE)



cat( format( Sys.time(), "%a %b %d %X %Y"), '-', 'Finished train and test datasets for', sp.n, 'with ', lim, 'lines...', '\n')

save.image("./outputs/my_analysis.rData")

finished_time = Sys.time()
cat( format( finished_time, "%a %b %d %X %Y"), '-', 'FINISHED', '\n')
write(format( finished_time, "%a %b %d %X %Y"), file=paste('./outputs/', '/FINISHED.txt', sep=""))
