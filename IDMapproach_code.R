##==================================#
#                                   #
#      IDM approach analysis        #
#   Doñana Nathional Park Analysis  #
#                                   #
##==================================#

####
#  PCOUNT MODEL (ABUNDANCE)  (backward stepwise selection procedure)
####


# Load needed libraries
library(lattice)
library(parallel)
library(Rcpp)
library(unmarked)
# We select our working directory
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/publicacion TFM/escrito/maquetación/datos a compartir")

# Charge data recorded for the analysis
data <- read.csv("IDMapproach_data.csv")
# See the structure of data
head(data)

###### We are going to work with each species separately

#####################
### red deer data ###
#####################

##PREPARE DATA FRAMEWORK##

# we select the data of the species of interest
dataciervo <- data[data$sp == "ciervo",]
head(dataciervo)
# The data structure is as follows:
# Column 2: ID assigned to each of the cameras
# Column 3: represent the different species that we are considering (deer, fallow deer and wild boar)
# Columns 4-28: represent each of the counts (replicates) that we have carried out, 
#               in this case 25, indicating the number of animals detected, being NA 
#               when the camera was not active in that count
# Columns 29-63: predictor variables measured at each site and/or occasion (covariates),

# We are going to save in "y" the data on the detection of the species in each one
# of replicates carried out and in "n" the number of sites visited
y <- dataciervo[,c(3:27)]
n <- nrow(dataciervo)

# We store the covariates in the object "d1516.site"
std_covamb_ciervo <- scale(dataciervo[,c(28:35)])
covamb <- data.frame(scale(dataciervo[,c(28:35)]))#we standardize the covariates
names(covamb) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6")
typeuse <- data.frame(as.factor(dataciervo[,37])) #we consider the dominant soil type as a factor
names(typeuse) <- "typeuse"
year <- data.frame(dataciervo[,36])
names(year) <- "year"
d1516.site <- data.frame(covamb,typeuse,year) #we join all the covariates that characterize each site

# We also store the covariates that vary for each observation (replication), in this case the 
# time of visit, in the object "d1516.obs"
d1516.obs <- data.frame(as.vector(t(dataciervo[,c(38:62)])))
names(d1516.obs) <- "time"

# Finally we put in a single object the history of detections and covariates
# in a special object that unmarked understands, we'll call it "d1516c"
d1516c <- unmarkedFramePCount(y = y, siteCovs = d1516.site, obsCovs = d1516.obs)

# We look the structure and ask for a summary:
head(d1516c)
summary(d1516c)


##MODEL SELECTION##

# We already have the deer data, and we proceed to choose the type of distribution
# that best fits our data
fmPOIS<- pcount(~1 ~1, d1516c, mixture="P",K=150) #model with poisson distribution
fmZIP <- pcount(~1 ~1, d1516c, mixture="ZIP",K=150) #model with negative binomial distribution
fmNB <- pcount(~1 ~1, d1516c, mixture="NB",K=150) #model with zero inflated distribution

# we look at the fit of the model, taking into account the AIC
fmlist<-fitList("POIS" = fmPOIS,"ZIP" = fmZIP,"NB" = fmNB)
modSel(fmlist)
# We see that the NB model has the lowest AIC, with more than two units of difference
# That is why we are left with this type of distribution: negative binomial (NB)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)

#_Observation process_#
am0<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am2<-pcount(~typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am3<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am4<-pcount(~time+typeuse ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am0 = am0,am2 = am2,am3 = am3,am4 = am4)
modSel(fmlist)

# Now we try to remove a second covariate, to see if it improves.
am5<-pcount(~time ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, d1516c, mixture="NB", K = 150)
am6<-pcount(~year ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, d1516c, mixture="NB", K = 150)

fmlist<-fitList(am0 = am0,am2 = am2,am3 = am3,am4 = am4,am5 = am5,am6 = am6)
modSel(fmlist)
# We see that none of the models in this step improve the previous model (m3), for
# what we are left with the m3 model

# With the best model, we move on to the next process

#_Ecological process_#
am.0<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.1<-pcount(~time+year ~dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.2<-pcount(~time+year ~dvera+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.3<-pcount(~time+year ~dvera+dwat+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.4<-pcount(~time+year ~dvera+dwat+v1+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.5<-pcount(~time+year ~dvera+dwat+v1+v2+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.6<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v5+v6+year, mixture="NB", d1516c, K = 150)
am.7<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v6+year, mixture="NB", d1516c, K = 150)
am.8<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+year, mixture="NB", d1516c, K = 150)
am.9<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, mixture="NB", d1516c, K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am.0 = am.0,am.1 = am.1,am.2 = am.2,am.3 = am.3,am.4 = am.4,
                am.5 = am.5,am.6 = am.6,am.7 = am.7,am.8 = am.8,am.9 = am.9)
modSel(fmlist)
# the best model is am6, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_v4" = am.6,  
                    "no_v3" = am.5,
                    "no_v1" = am.3,
                    "no_v2" = am.4,
                    "no_v5" = am.7,
                    "no_v6" = am.8,
                    "no_year" = am.9,
                    "no_dwat" = am.2,
                    "null" = am.0)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn) 

Mod.avg_ciervo <- model.avg(Cand.models, fit = TRUE)
Mod.avg_ciervo

# we can look at a summary of this average model
Sum_Mod.avg_ciervo <- summary(model.avg(Cand.models))
Sum_Mod.avg_ciervo


##MODEL PREDICTION##

library(sp)
library(raster)
library(rgdal)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

# we load some new layers of "year", to observe how the occupancy changes between years. 
nd <- nrow(datadoñana)
vyear1 <- rep(1, times = nd )
vyear0 <- rep(0, times = nd )
datadoñana_df$YEAR1 <- vyear1
datadoñana_df$YEAR0 <- vyear0

# We store each of the covariates in a raster
dvera <- rasterFromXYZ(datadoñana_df[,c("X","Y","DVERA_mean")],
                       crs=NA)
dwat <- rasterFromXYZ(datadoñana_df[,c("X","Y","DWAT_meanm")],
                      crs=NA)
v1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V1_sumsum")],
                    crs=NA)
v2 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V2_sumsum")],
                    crs=NA)
v3 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V3_sumsum")],
                    crs=NA)
v4 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V4_sumsum")],
                    crs=NA)
v5 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V5_sumsum")],
                    crs=NA)
v6 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V6_sumsum")],
                    crs=NA)
year1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR1")],
                       crs=NA)
year0 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR0")],
                       crs=NA)

# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
std_covamb_ciervo

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-15.891329)/14.972568  
dwat.s <- (dwat-6.121166)/4.177509
v1.s <- (v1-20.694915)/31.053488
v2.s <- (v2-32.864407)/38.440934
v3.s <- (v3-26.491525)/36.015980
v4.s <- (v4-3.067797)/15.858061
v5.s <- (v5-4.338983)/14.577510
v6.s <- (v6-10.542373)/23.819235


# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
library(AICcmodavg)
extractX(cand.set = Cand.models, parm.type = "lambda") #we look at which are the 
                                    #variables present in the set of final models

ef1 <- stack(dvera.s,dwat.s,v1.s,v2.s,v3.s,v4.s,v5.s,v6.s,year0) # We join all the raster layers
names(ef1) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6","year")


ef1_df <- data.frame( rasterToPoints( ef1 ) ) #convert the raster multilayer into a data frame

# We make the prediction at DNP level with the average model
Mod.avg_ciervo_pred <- predict(Mod.avg_ciervo, newdata=ef1_df, type="state", full=TRUE, 
                               re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_ciervo_pred$fit))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualize the prediction

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_ciervo_pcount_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###################################################################

######################
### wild boar data ###
######################

##PREPARE DATA FRAMEWORK##

# we select the data of the species of interest
datajabali <- data[data$sp == "jabali",]
head(datajabali)
# The data structure is as follows:
# Column 2: ID assigned to each of the cameras
# Column 3: represent the different species that we are considering (deer, fallow deer and wild boar)
# Columns 4-28: represent each of the counts (replicates) that we have carried out, 
#               in this case 25, indicating the number of animals detected, being NA 
#               when the camera was not active in that count
# Columns 29-63: predictor variables measured at each site and/or occasion (covariates),

# We are going to save in "y" the data on the detection of the species in each one
# of replicates carried out and in "n" the number of sites visited
y <- datajabali[,c(3:27)]
n <- nrow(datajabali)

# We store the covariates in the object "d1516.site"
std_covamb_jabali <- scale(datajabali[,c(28:35)])
covamb <- data.frame(scale(datajabali[,c(28:35)]))#we standardize the covariates
names(covamb) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6")
typeuse <- data.frame(as.factor(datajabali[,37])) #we consider the dominant soil type as a factor
names(typeuse) <- "typeuse"
year <- data.frame(datajabali[,36])
names(year) <- "year"
d1516.site <- data.frame(covamb,typeuse,year) #we join all the covariates that characterize each site

# We also store the covariates that vary for each observation (replication), in this case the 
# time of visit, in the object "d1516.obs"
d1516.obs <- data.frame(as.vector(t(datajabali[,c(38:62)])))
names(d1516.obs) <- "time"

# Finally we put in a single object the history of detections and covariates
# in a special object that unmarked understands, we'll call it "d1516c"
d1516c <- unmarkedFramePCount(y = y, siteCovs = d1516.site, obsCovs = d1516.obs)

# We look the structure and ask for a summary:
head(d1516c)
summary(d1516c)


##MODEL SELECTION##

# We already have the wild boar data, and we proceed to choose the type of distribution
# that best fits our data
fmPOIS<- pcount(~1 ~1, d1516c, mixture="P",K=150) #model with poisson distribution
fmZIP <- pcount(~1 ~1, d1516c, mixture="ZIP",K=150) #model with negative binomial distribution
fmNB <- pcount(~1 ~1, d1516c, mixture="NB",K=150) #model with zero inflated distribution

# we look at the fit of the model, taking into account the AIC
fmlist<-fitList("POIS" = fmPOIS,"ZIP" = fmZIP,"NB" = fmNB)
modSel(fmlist)
# We see that the NB model has the lowest AIC, with more than two units of difference
# That is why we are left with this type of distribution: negative binomial (NB)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)

#_Observation process_#
am0<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am2<-pcount(~typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am3<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am4<-pcount(~time+typeuse ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am0 = am0,am2 = am2,am3 = am3,am4 = am4)
modSel(fmlist)
# We see that none of the models in this step improve the previous model (m0), for
# what we are left with the m0 model

# With the best model, we move on to the next process

#_Ecological process_#
am.0<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.1<-pcount(~time+typeuse+year ~dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.2<-pcount(~time+typeuse+year ~dvera+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.3<-pcount(~time+typeuse+year ~dvera+dwat+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.4<-pcount(~time+typeuse+year ~dvera+dwat+v1+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.5<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.6<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v5+v6+year, mixture="NB", d1516c, K = 150)
am.7<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v6+year, mixture="NB", d1516c, K = 150)
am.8<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+year, mixture="NB", d1516c, K = 150)
am.9<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, mixture="NB", d1516c, K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am.0 = am.0,am.1 = am.1,am.2 = am.2,am.3 = am.3,am.4 = am.4,
                am.5 = am.5,am.6 = am.6,am.7 = am.7,am.8 = am.8,am.9 = am.9)
modSel(fmlist)
# the best model is am6, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_v4" = am.6,  
                    "no_year" = am.9,
                    "null" = am.0)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn) 

Mod.avg_jabali <- model.avg(Cand.models, fit = TRUE)
Mod.avg_jabali

# we can look at a summary of this average model
Sum_Mod.avg_jabali <- summary(model.avg(Cand.models))
Sum_Mod.avg_jabali


##MODEL PREDICTION##

library(sp)
library(raster)
library(rgdal)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

# we load some new layers of "year", to observe how the occupancy changes between years. 
nd <- nrow(datadoñana)
vyear1 <- rep(1, times = nd )
vyear0 <- rep(0, times = nd )
datadoñana_df$YEAR1 <- vyear1
datadoñana_df$YEAR0 <- vyear0

# We store each of the covariates in a raster
dvera <- rasterFromXYZ(datadoñana_df[,c("X","Y","DVERA_mean")],
                       crs=NA)
dwat <- rasterFromXYZ(datadoñana_df[,c("X","Y","DWAT_meanm")],
                      crs=NA)
v1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V1_sumsum")],
                    crs=NA)
v2 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V2_sumsum")],
                    crs=NA)
v3 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V3_sumsum")],
                    crs=NA)
v4 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V4_sumsum")],
                    crs=NA)
v5 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V5_sumsum")],
                    crs=NA)
v6 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V6_sumsum")],
                    crs=NA)
year1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR1")],
                       crs=NA)
year0 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR0")],
                       crs=NA)

# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
std_covamb_jabali

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-15.891329)/14.972568 
dwat.s <- (dwat-6.121166)/4.177509
v1.s <- (v1-20.694915)/31.053488
v2.s <- (v2-32.864407)/38.440934
v3.s <- (v3-26.491525)/36.015980
v4.s <- (v4-3.067797)/15.858061
v5.s <- (v5-4.338983)/14.577510
v6.s <- (v6-10.542373)/23.819235

# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
library(AICcmodavg)
extractX(cand.set = Cand.models, parm.type = "lambda") #we look at which are the 
#variables present in the set of final models

ef1 <- stack(dvera.s,dwat.s,v1.s,v2.s,v3.s,v4.s,v5.s,v6.s,year0)# We join all the raster layers 
names(ef1) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6","year")


ef1_df <- data.frame( rasterToPoints( ef1 ) )#convert the raster multilayer into a data frame 

# We make the prediction at DNP level with the average model
Mod.avg_jabali_pred <- predict(Mod.avg_jabali, newdata=ef1_df, type="state", full=TRUE, 
                               re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_jabali_pred$fit))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualize the prediction

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_jabali_pcount_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###################################################################

########################
### fallow deer data ###
########################


# we select the data of the species of interest
datagamo <- data[data$sp == "gamo",]
head(datagamo)
# The data structure is as follows:
# Column 2: ID assigned to each of the cameras
# Column 3: represent the different species that we are considering (deer, fallow deer and wild boar)
# Columns 4-28: represent each of the counts (replicates) that we have carried out, 
#               in this case 25, indicating the number of animals detected, being NA 
#               when the camera was not active in that count
# Columns 29-63: predictor variables measured at each site and/or occasion (covariates),

# We are going to save in "y" the data on the detection of the species in each one
# of replicates carried out and in "n" the number of sites visited
y <- datagamo[,c(3:27)]
n <- nrow(datagamo)

# We store the covariates in the object "d1516.site"
std_covamb_gamo <- scale(datagamo[,c(28:35)])
covamb <- data.frame(scale(datagamo[,c(28:35)]))#we standardize the covariates
names(covamb) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6")
typeuse <- data.frame(as.factor(datagamo[,37])) #we consider the dominant soil type as a factor
names(typeuse) <- "typeuse"
year <- data.frame(datagamo[,36])
names(year) <- "year"
d1516.site <- data.frame(covamb,typeuse,year) #we join all the covariates that characterize each site

# We also store the covariates that vary for each observation (replication), in this case the 
# time of visit, in the object "d1516.obs"
d1516.obs <- data.frame(as.vector(t(datagamo[,c(38:62)])))
names(d1516.obs) <- "time"

# Finally we put in a single object the history of detections and covariates
# in a special object that unmarked understands, we'll call it "d1516c"
d1516c <- unmarkedFramePCount(y = y, siteCovs = d1516.site, obsCovs = d1516.obs)

# We look the structure and ask for a summary:
head(d1516c)
summary(d1516c)


##MODEL SELECTION##

# We already have the falloww deer data, and we proceed to choose the type of distribution
# that best fits our data
fmPOIS<- pcount(~1 ~1, d1516c, mixture="P",K=150) #model with poisson distribution
fmZIP <- pcount(~1 ~1, d1516c, mixture="ZIP",K=150) #model with negative binomial distribution
fmNB <- pcount(~1 ~1, d1516c, mixture="NB",K=150) #model with zero inflated distribution

# we look at the fit of the model, taking into account the AIC
fmlist<-fitList("POIS" = fmPOIS,"ZIP" = fmZIP,"NB" = fmNB)
modSel(fmlist)
# We see that the NB model has the lowest AIC, with more than two units of difference
# That is why we are left with this type of distribution: negative binomial (NB)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)

#_Observation process_#
am0<-pcount(~time+typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am2<-pcount(~typeuse+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am3<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)
am4<-pcount(~time+typeuse ~dvera+dwat+v1+v2+v3+v4+v5+v6, d1516c, mixture="NB", K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am0 = am0,am2 = am2,am3 = am3,am4 = am4)
modSel(fmlist)

# Now we try to remove a second covariate, to see if it improves.
am5<-pcount(~time ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, d1516c, mixture="NB", K = 150)
am6<-pcount(~year ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, d1516c, mixture="NB", K = 150)

fmlist<-fitList(am0 = am0,am2 = am2,am3 = am3,am4 = am4,am5 = am5,am6 = am6)
modSel(fmlist)
# We see that none of the models in this step improve the previous model (m3), for
# what we are left with the m3 model

# With the best model, we move on to the next process

#_Ecological process_#
am.0<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.1<-pcount(~time+year ~dwat+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.2<-pcount(~time+year ~dvera+v1+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.3<-pcount(~time+year ~dvera+dwat+v2+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.4<-pcount(~time+year ~dvera+dwat+v1+v3+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.5<-pcount(~time+year ~dvera+dwat+v1+v2+v4+v5+v6+year, mixture="NB", d1516c, K = 150)
am.6<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v5+v6+year, mixture="NB", d1516c, K = 150)
am.7<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v6+year, mixture="NB", d1516c, K = 150)
am.8<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+year, mixture="NB", d1516c, K = 150)
am.9<-pcount(~time+year ~dvera+dwat+v1+v2+v3+v4+v5+v6, mixture="NB", d1516c, K = 150)

# we choose the model with the lowest AIC
fmlist<-fitList(am.0 = am.0,am.1 = am.1,am.2 = am.2,am.3 = am.3,am.4 = am.4,
                am.5 = am.5,am.6 = am.6,am.7 = am.7,am.8 = am.8,am.9 = am.9)
modSel(fmlist)
# the best model is am2, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_dwat" = am.2,   #datos de 2015 y 2016
                    "no_v6" = am.8,
                    "no_v3" = am.5,
                    "no_year" = am.9,
                    "no_dvera" = am.1,
                    "no_v2" = am.4,
                    "no_v4" = am.6,
                    "no_v1" = am.3,
                    "null" = am.0)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn) 

Mod.avg_gamo <- model.avg(Cand.models, fit = TRUE)
Mod.avg_gamo

# we can look at a summary of this average model
Sum_Mod.avg_gamo <- summary(model.avg(Cand.models))
Sum_Mod.avg_gamo


##MODEL PREDICTION##

library(sp)
library(raster)
library(rgdal)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

# we load some new layers of "year", to observe how the occupancy changes between years. 
nd <- nrow(datadoñana)
vyear1 <- rep(1, times = nd )
vyear0 <- rep(0, times = nd )
datadoñana_df$YEAR1 <- vyear1
datadoñana_df$YEAR0 <- vyear0

# We store each of the covariates in a raster
dvera <- rasterFromXYZ(datadoñana_df[,c("X","Y","DVERA_mean")],
                       crs=NA)
dwat <- rasterFromXYZ(datadoñana_df[,c("X","Y","DWAT_meanm")],
                      crs=NA)
v1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V1_sumsum")],
                    crs=NA)
v2 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V2_sumsum")],
                    crs=NA)
v3 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V3_sumsum")],
                    crs=NA)
v4 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V4_sumsum")],
                    crs=NA)
v5 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V5_sumsum")],
                    crs=NA)
v6 <- rasterFromXYZ(datadoñana_df[,c("X","Y","V6_sumsum")],
                    crs=NA)
year1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR1")],
                       crs=NA)
year0 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR0")],
                       crs=NA)

# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
std_covamb_jabali

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-15.891329)/14.972568 
dwat.s <- (dwat-6.121166)/4.177509
v1.s <- (v1-20.694915)/31.053488
v2.s <- (v2-32.864407)/38.440934
v3.s <- (v3-26.491525)/36.015980
v4.s <- (v4-3.067797)/15.858061
v5.s <- (v5-4.338983)/14.577510
v6.s <- (v6-10.542373)/23.819235

# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
library(AICcmodavg)
extractX(cand.set = Cand.models, parm.type = "lambda") #we look at which are the 
#variables present in the set of final models

ef1 <- stack(dvera.s,dwat.s,v1.s,v2.s,v3.s,v4.s,v5.s,v6.s,year0)# We join all the raster layers 
names(ef1) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6","year")


ef1_df <- data.frame( rasterToPoints( ef1 ) )#convert the raster multilayer into a data frame 

# We make the prediction at DNP level with the average model
Mod.avg_gamo_pred <- predict(Mod.avg_gamo, newdata=ef1_df, type="state", full=TRUE, 
                             re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_gamo_pred$fit))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualize the prediction

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_gamo_pcount_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###########################END#####################################