##==================================#
#                                   #
#      RSF approach analysis        #
#   Doñana Nathional Park Analysis  #
#                                   #
##==================================#

####
#  GLMMTMB MODEL (HABITAT USE)  (backward stepwise selection procedure)
####


# Load needed libraries
library(adehabitatHS)
library(lubridate)
library(caTools)
library(rgdal)
# We select our working directory
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/publicacion TFM/escrito/maquetación/datos a compartir")

# Charge data recorded for the analysis
locs_env_ok <- read.csv("RSFapproach_data.csv")

###### We are going to work with each species separately

#####################
### red deer data ###
#####################

##PREPARE DATA FRAMEWORK##

# we select the data of the species of interest
locs_hs<-subset(locs_env_ok, SPECIES == 1)
locs_hs$NAME <- as.factor(locs_hs$NAME) #we consider the code of each animal (NAME) as a factor

# Get the identification data
locs_data <- locs_hs[,c(2:4)]
sp_data <- locs_hs[,c(13:15)]

# Standardize the variables at the species level
scale_covar<- scale(locs_hs[,c(5:12)])

# Group the identification data with the standardized variables
locs_hs <- data.frame(locs_data, scale_covar, sp_data)


##MODEL SELECTION##

library(glmmTMB)
library(lmerTest)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)
mm0 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm1 <- glmmTMB(Type~dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm2 <- glmmTMB(Type~dvera+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm3 <- glmmTMB(Type~dvera+dwat+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm4 <- glmmTMB(Type~dvera+dwat+v1+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm5 <- glmmTMB(Type~dvera+dwat+v1+v2+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm6 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm7 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm8 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm9 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))

# we choose the model with the lowest AIC
library(AICcmodavg)
models <- list("mm2" = mm2,
               "mm9" = mm9)

aictab(models)

# We see that most of the models have overparameterization problems. Focusing on
# the models that do not have this problem, we keep the mm9 as a reference, which
# is the one with the best fit, as base models.  

# Now we try to start with mm9 model as colapsed model(mm0).
mm0 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm1 <- glmmTMB(Type~dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm2 <- glmmTMB(Type~dvera+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm3 <- glmmTMB(Type~dvera+dwat+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm4 <- glmmTMB(Type~dvera+dwat+v1+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm5 <- glmmTMB(Type~dvera+dwat+v1+v2+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm6 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm7 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm8 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))

# we choose the model with the lowest AIC
models <- list("mm0" = mm0,
               "mm1" = mm1,
               "mm2" = mm2,
               "mm3" = mm3,
               "mm4" = mm4,
               "mm5" = mm5,
               "mm6" = mm6,
               "mm7" = mm7,
               "mm8" = mm8)

aictab(models)

# the best model is mm1, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_dvera" = mm1,  
                    "no_v6" = mm8,
                    "no_v1" = mm3,
                    "no_v5" = mm7,
                    "full" = mm0)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn)

avg_model_ciervo <- model.avg(Cand.models)
avg_model_ciervo

# we can look at a summary of this average model
average_model_ciervo <- summary(model.avg(Cand.models))
average_model_ciervo

##MODEL PREDICTION##

library(sp)
library(raster)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

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

# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
scale_covar

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-26.093009)/22.998754
dwat.s <- (dwat-7.697859)/5.431122
v1.s <- (v1-21.374431)/35.626516
v2.s <- (v2-35.413868)/42.813352   
v3.s <- (v3-24.678873)/36.913927
v4.s <- (v4-2.553088)/12.942146
v5.s <- (v5-2.373673)/9.737937
v6.s <- (v6-11.935536)/28.346536

COLLAR <- ((v6+1)/(v6+1))+7000

# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
ef1 <- stack(dvera.s,dwat.s,v1.s,v2.s,v3.s,v4.s,v5.s,v6.s,COLLAR)# We join all the raster layers 
names(ef1) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6","COLLAR")

ef1_df <- data.frame( rasterToPoints(ef1) )#convert the raster multilayer into a data frame  

# We make the prediction at DNP level with the average model
Mod.avg_ciervo_pred <- predict(avg_model_ciervo, newdata=ef1_df, type="response", 
                               full=TRUE, re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_ciervo_pred))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualize the prediction

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_ciervo_glmmTMB_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###################################################################

########################
### fallow deer data ###
########################

##PREPARE DATA FRAMEWORK##

# we select the data of the species of interest
locs_hs<-subset(locs_env_ok, SPECIES == 2)
locs_hs$NAME <- as.factor(locs_hs$NAME) #we consider the code of each animal (NAME) as a factor

# Get the identification data
locs_data <- locs_hs[,c(2:4)]
sp_data <- locs_hs[,c(13:15)]

# Standardize the variables at the species level
scale_covar<- scale(locs_hs[,c(5:12)])

# Group the identification data with the standardized variables
locs_hs <- data.frame(locs_data, scale_covar, sp_data)

##MODEL SELECTION##

library(glmmTMB)
library(lmerTest)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)
mm0 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm1 <- glmmTMB(Type~dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm2 <- glmmTMB(Type~dvera+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm3 <- glmmTMB(Type~dvera+dwat+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm4 <- glmmTMB(Type~dvera+dwat+v1+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm5 <- glmmTMB(Type~dvera+dwat+v1+v2+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm6 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm7 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm8 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm9 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))


# we choose the model with the lowest AIC
library(AICcmodavg)
models <- list("mm0" = mm0,
               "mm1" = mm1,
               "mm2" = mm2,
               "mm3" = mm3,
               "mm4" = mm4,
               "mm5" = mm5,
               "mm6" = mm6,
               "mm7" = mm7,
               "mm8" = mm8,
               "mm9" = mm9)

aictab(models)

# the best model is mm6, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_v4" = mm6,     
                    "no_v2" = mm4,
                    "no_year" = mm9,
                    "no_v3" = mm5,
                    "full" = mm0)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn)

avg_model_gamo <- model.avg(Cand.models)
avg_model_gamo

# we can look at a summary of this average model
average_model_gamo <- summary(model.avg(Cand.models))
average_model_gamo


##MODEL PREDICTION##

library(sp)
library(raster)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

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

# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
scale_covar

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-15.802196)/14.05384  
dwat.s <- (dwat-4.035350)/2.63624
v1.s <- (v1-2.734345)/12.77665
v2.s <- (v2-16.582187)/30.50410   
v3.s <- (v3-36.360887)/36.21834
v4.s <- (v4-7.080408)/22.76735
v5.s <- (v5-8.226637)/19.92417
v6.s <- (v6-28.698055)/38.80786

COLLAR <- ((v6+1)/(v6+1))+7000

# we load some new layers of "year", to observe how the occupancy changes between years.
nd <- nrow(datadoñana)
vyear1 <- rep(1, times = nd )
vyear0 <- rep(0, times = nd )
datadoñana_df$YEAR1 <- vyear1

year1 <- rasterFromXYZ(datadoñana_df[,c("X","Y","YEAR1")],
                       crs=NA)

# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
ef1 <- stack(dvera.s,dwat.s,v1.s,v2.s,v3.s,v4.s,v5.s,v6.s,COLLAR)# We join all the raster layers 
names(ef1) <- c("dvera","dwat","v1","v2","v3","v4","v5","v6","COLLAR")

ef1_df <- data.frame( rasterToPoints(ef1) )#convert the raster multilayer into a data frame 

# We make the prediction at DNP level with the average model
Mod.avg_gamo_pred <- predict(avg_model_gamo, newdata=ef1_df, type="response", 
                             full=TRUE, re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_gamo_pred))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualize the prediction

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_gamo_glmmTMB_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###################################################################

######################
### wild boar data ###
######################

##PREPARE DATA FRAMEWORK##

# we select the data of the species of interest
locs_hs<-subset(locs_env_ok, SPECIES == 3) 
locs_hs$NAME <- as.factor(locs_hs$NAME) #we consider the code of each animal (NAME) as a factor

# Get the identification data
locs_data <- locs_hs[,c(2:4)]
sp_data <- locs_hs[,c(13:15)]

# Standardize the variables at the species level
scale_covar<- scale(locs_hs[,c(5:12)])

# Group the identification data with the standardized variables
locs_hs <- data.frame(locs_data, scale_covar, sp_data)

##MODEL SELECTION##

library(glmmTMB)
library(lmerTest)

# In this part, we proceed to the selection of the best fit model or group of models
# by backward stepwise selection procedure (starting with the complex model and removing
# covariates step by step until the fit does not improve)
mm0 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm1 <- glmmTMB(Type~dwat+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm2 <- glmmTMB(Type~dvera+v1+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm3 <- glmmTMB(Type~dvera+dwat+v2+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm4 <- glmmTMB(Type~dvera+dwat+v1+v3+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm5 <- glmmTMB(Type~dvera+dwat+v1+v2+v4+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm6 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v5+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm7 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v6+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm8 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+year+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm9 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))

# we choose the model with the lowest AIC
library(AICcmodavg)
models <- list("mm0" = mm0,
               "mm1" = mm1,
               "mm2" = mm2,
               "mm3" = mm3,
               "mm4" = mm4,
               "mm5" = mm5,
               "mm6" = mm6,
               "mm7" = mm7,
               "mm8" = mm8,
               "mm9" = mm9)

aictab(models)

# the best model is mm9, with more than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# mm9 for the next step.
mm0 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm1 <- glmmTMB(Type~dwat+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm2 <- glmmTMB(Type~dvera+v1+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm3 <- glmmTMB(Type~dvera+dwat+v2+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm4 <- glmmTMB(Type~dvera+dwat+v1+v3+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm5 <- glmmTMB(Type~dvera+dwat+v1+v2+v4+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm6 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v5+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm7 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v6+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))
mm8 <- glmmTMB(Type~dvera+dwat+v1+v2+v3+v4+v5+(1|COLLAR),data=locs_hs, family=binomial(link="logit"))

# we choose the model with the lowest AIC
models <- list("mm0" = mm0,
               "mm1" = mm1,
               "mm2" = mm2,
               "mm3" = mm3,
               "mm4" = mm4,
               "mm5" = mm5,
               "mm6" = mm6,
               "mm7" = mm7,
               "mm8" = mm8)

aictab(models)

# the best model is mm1, with less than two units difference from the
# collapsed model (m0). For this reason, we look to stay with the 
# collapsed model, and all those that have a deltaAIC value lower than 2 units.

# We keep a list of the final models with a name assigned to each one.
Cand.models <- list("no_dvera" = mm1,
                    "no_v5" = mm7,
                    "no_v1" = mm3,
                    "full" = mm0, 
                    "no_v4" = mm6)

# We now look at averaging the set of best-fit models, with the
# modavg() function from MuMIn package
library(MuMIn)

avg_model_jabali <- model.avg(Cand.models)
avg_model_jabali

# we can look at a summary of this average model
average_model_jabali <- summary(model.avg(Cand.models))
average_model_jabali


##MODEL PREDICTION##

library(sp)
library(raster)

# To make the prediction we load the covariate layer of the area we are interested
# in predicting (DOÑANA NATHIONAL PARK)
datadoñana <- readOGR("DNP_covariates.shp")
datadoñana

datadoñana_df <- data.frame(datadoñana)
head(datadoñana_df)

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


# We standardize the variables of the area to be predicted in the same way that
# we have done with the variables of the model
scale_covar

# And with these we standardize the covariates of the map with "scaled:centre" as
# the mean, and "scaled:scale" as the standard deviation
dvera.s <- (dvera-17.891775)/15.480415   
dwat.s <- (dwat-3.554962)/2.057801
v1.s <- (v1-11.693954)/24.704517
v2.s <- (v2-28.849652)/37.692168   
v3.s <- (v3-27.996790)/33.505457
v4.s <- (v4-5.158909)/19.126814
v5.s <- (v5-3.918138)/12.610507
v6.s <- (v6-21.958801)/33.472877

COLLAR <- ((v6+1)/(v6+1))+7000

# We join all the raster layers and make the prediction according to the estimates
#obtained with the best fit model (the average model)
ef1 <- stack(dvera.s,dwat.s,v1.s,v3.s,v4.s,v5.s,v6.s,COLLAR)# We join all the raster layers 
names(ef1) <- c("dvera","dwat","v1","v3","v4","v5","v6","COLLAR")

ef1_df <- data.frame( rasterToPoints(ef1) )#convert the raster multilayer into a data frame

# We make the prediction at DNP level with the average model
Mod.avg_jabali_pred <- predict(avg_model_jabali, newdata=ef1_df, type="response", 
                               full=TRUE, re.form = ~0, allow.new.levels=TRUE)

# Create a data frame with each of the predicted values together with its coordinate value
df_pred <- as.data.frame(cbind(x=ef1_df$x, y=ef1_df$y, pred=Mod.avg_jabali_pred))
r_pred <- rasterFromXYZ(df_pred) # create a raster from this data frame

plot(r_pred, axes=FALSE, col=topo.colors(100)) #visualizamos la predicción

# we save the raster layer in our documents
route_rst<-paste("capas_finales/pred_jabali_glmmTMB_1516.tiff")
writeRaster(r_pred, route_rst,"GTiff", overwrite=TRUE)

###################################################################
###########################END#####################################