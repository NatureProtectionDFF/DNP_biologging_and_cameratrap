##==================================#
#                                   #
#      Predictions agreement        #
#   Doñana Nathional Park Analysis  #
#                                   #
##==================================#

####
#  Cohen's weighted kappa (with quantile reclassified predictions)
####


# Load needed libraries
library(raster)
library(mnormt)
library(psych)
# We select our working directory
setwd("D:/OneDrive - Universidad de Castilla-La Mancha/DOCTORADO UCLM/publicacion TFM/escrito/maquetación/datos a compartir")

###### We are going to work with each species separately

#####################
### red deer data ###
#####################

##DOÑANA NATHIONAL PARK##

# Charge data recorded for the analysis
# We load the rasters of the predictions at the Park level (management areas)
Q_foto <- raster("capas_finales/pred_ciervo_pcount_1516.tif")

Q_tele <- raster("capas_finales/pred_ciervo_glmmTMB_1516.tif")

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the camera trap map. 
quantile(Q_foto)

quan_1516f <- matrix(c(0,  12.1847215, 1,          
                       12.1847215,  22.8368874, 2,
                       22.8368874, 35.4550076, 3,
                       35.4550076, 252, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
foto4 <- reclassify(Q_foto, quan_1516f)

my_col = rev(terrain.colors(n = 4))
plot(foto4, legend = FALSE, col = my_col, main= "ciervo fototrampeo")

# since this is a categorical variable, we plot a map with a categorical legend
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
library(rgdal)

reserva <- readOGR(dsn="datos/FINCAS_DN.shp", layer="FINCAS_DN")
plot(reserva, add = TRUE, border = "red", lwd = 2)
calibra <- readOGR(dsn="datos/zona_calibración_rec.shp", layer="zona_calibración_rec")
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the biologging map.
quantile(Q_tele)

quan_1516t <- matrix(c(0,  0.3563618, 1,        
                       0.3563618,  0.4588813, 2,
                       0.4588813, 0.5313628, 3,
                       0.5313628, 0.70, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
tele4 <- reclassify(Q_tele, quan_1516t)

my_col = rev(terrain.colors(n = 4))
plot(tele4, legend = FALSE, col = my_col, main= "ciervo telemetria")

# since this is a categorical variable, we plot a map with a categorical legend
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
plot(reserva, add = TRUE, border = "red", lwd = 2)
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto4[]),na.omit(tele4[]))) 

##BIOLOGICAL RESERVE##

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Camera Trap prediction
foto_reserva <- crop(mask(foto4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_reserva, legend = FALSE, col = my_col, main= "reserva ciervo fototrampeo")

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Biologging prediction
tele_reserva <- crop(mask(tele4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_reserva, legend = FALSE, col = my_col, main= "reserva ciervo telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto_reserva[]),na.omit(tele_reserva[])))

##CALIBRATION REGION##

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Camera Trap prediction
foto_calibra <- crop(mask(foto4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_calibra, legend = FALSE, col = my_col, main= "zona calibración ciervo fototrampeo")

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Biologging prediction
tele_calibra <- crop(mask(tele4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_calibra, legend = FALSE, col = my_col, main= "zona calibración ciervo telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto_calibra[]),na.omit(tele_calibra[])))

###################################################################
###################################################################

######################
### wild boar data ###
######################

##DOÑANA NATHIONAL PARK##

# Charge data recorded for the analysis
# We load the rasters of the predictions at the Park level (management areas)
Q_foto <- raster("capas_finales/pred_jabali_pcount_1516.tif")

Q_tele <- raster("capas_finales/pred_jabali_glmmTMB_1516.tif")

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the camera trap map.
quantile(Q_foto)

quan_1516f <- matrix(c(0,  11.2000875, 1,          
                       11.2000875,  32.8410339, 2,
                       32.8410339, 67.7764683, 3,
                       67.7764683, 608, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
foto4 <- reclassify(Q_foto, quan_1516f)

# since this is a categorical variable, we plot a map with a categorical legend
my_col = rev(terrain.colors(n = 4))
plot(foto4, legend = FALSE, col = my_col, main= "jabalí fototrampeo")
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
library(rgdal)

reserva <- readOGR(dsn="datos/FINCAS_DN.shp", layer="FINCAS_DN")
plot(reserva, add = TRUE, border = "red", lwd = 2)
calibra <- readOGR(dsn="datos/zona_calibración_rec.shp", layer="zona_calibración_rec")
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the biologging map.
quantile(Q_tele)

quan_1516t <- matrix(c(0,  0.0889961198, 1,        
                       0.0889961198,  0.2194845602, 2,
                       0.2194845602, 0.3711865097, 3,
                       0.3711865097, 0.75, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
tele4 <- reclassify(Q_tele, quan_1516t)

my_col = rev(terrain.colors(n = 4))
plot(tele4, legend = FALSE, col = my_col, main= "jabalí telemetria")

# since this is a categorical variable, we plot a map with a categorical legend
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
plot(reserva, add = TRUE, border = "red", lwd = 2)
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto4[]),na.omit(tele4[])))

##BIOLOGICAL RESERVE##

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Camera Trap prediction
foto_reserva <- crop(mask(foto4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_reserva, legend = FALSE, col = my_col, main= "reserva jabalí fototrampeo")

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Biologging prediction
tele_reserva <- crop(mask(tele4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_reserva, legend = FALSE, col = my_col, main= "reserva jabalí telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
cohen.kappa(x=cbind(na.omit(foto_reserva[]),na.omit(tele_reserva[])))

##CALIBRATION REGION##

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Camera Trap prediction
foto_calibra <- crop(mask(foto4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_calibra, legend = FALSE, col = my_col, main= "zona calibración jabalí fototrampeo")

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Biologging prediction
tele_calibra <- crop(mask(tele4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_calibra, legend = FALSE, col = my_col, main= "zona calibración jabalí telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto_calibra[]),na.omit(tele_calibra[])))

###################################################################
###################################################################

########################
### fallow deer data ###
########################

##DOÑANA NATHIONAL PARK##

# Charge data recorded for the analysis
# We load the rasters of the predictions at the Park level (management areas)
Q_foto <- raster("capas_finales/pred_gamo_pcount_1516.tif")

Q_tele <- raster("capas_finales/pred_gamo_glmmTMB_1516.tif")

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the camera trap map. 
quantile(Q_foto)

quan_1516f <- matrix(c(0,  2.513766, 1,          
                       2.513766,  8.526722, 2,
                       8.526722, 20.04485, 3,
                       20.04485, 1300, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
foto4 <- reclassify(Q_foto, quan_1516f)

my_col = rev(terrain.colors(n = 4))
plot(foto4, legend = FALSE, col = my_col, main= "gamo fototrampeo")

# since this is a categorical variable, we plot a map with a categorical legend
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
library(rgdal)

reserva <- readOGR(dsn="datos/FINCAS_DN.shp", layer="FINCAS_DN")
plot(reserva, add = TRUE, border = "red", lwd = 2)
calibra <- readOGR(dsn="datos/zona_calibración_rec.shp", layer="zona_calibración_rec")
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# We use the quantiles BASED ON THE WHOLE PREDICTED AREA (DNP) for the biologging map.
quantile(Q_tele)

quan_1516t <- matrix(c(0,  0.0906489529, 1,        
                       0.0906489529,  0.2365275249, 2,
                       0.2365275249, 0.4032026827, 3,
                       0.4032026827, 0.80, 4), 4,3, byrow = TRUE)

# We reclassify the map, we categorize it into 4 and we plot it to see what it looks like:
tele4 <- reclassify(Q_tele, quan_1516t)

my_col = rev(terrain.colors(n = 4))
plot(tele4, legend = FALSE, col = my_col, main= "gamo telemetria")

# since this is a categorical variable, we plot a map with a categorical legend
legend(x = 'toprigh', legend = c("low", "medium", "high", "very high"), fill = my_col)

# we incorporate the limitation of the reserve zone and calibration zone
plot(reserva, add = TRUE, border = "red", lwd = 2)
plot(calibra, add = TRUE, border = "blue", lwd = 2)

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto4[]),na.omit(tele4[])))

##BIOLOGICAL RESERVE##

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Camera Trap prediction
foto_reserva <- crop(mask(foto4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_reserva, legend = FALSE, col = my_col, main= "reserva gamo fototrampeo")

# We cut the prediction layer to the BIOLOGICAL RESERVE area, to compare
# predictions at this level for Biologging prediction
tele_reserva <- crop(mask(tele4, reserva), reserva)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_reserva, legend = FALSE, col = my_col, main= "reserva gamo telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto_reserva[]),na.omit(tele_reserva[])))

##CALIBRATION REGION##

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Camera Trap prediction
foto_calibra <- crop(mask(foto4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(foto_calibra, legend = FALSE, col = my_col, main= "zona calibración gamo fototrampeo")

# We cut the prediction layer to the CALIBRATION REGION area, to compare
# predictions at this level for Biologging prediction
tele_calibra <- crop(mask(tele4, calibra), calibra)

# We plot it to see what it looks like:
my_col = rev(terrain.colors(n = 4))
plot(tele_calibra, legend = FALSE, col = my_col, main= "zona calibración gamo telemetría")

# Calculate Cohen's kappa with the psych package. We stay with the weighted value
# Values range from -1 (reverse match) to 1 (full match). A value close to zero 
# would indicate random categorization.
cohen.kappa(x=cbind(na.omit(foto_calibra[]),na.omit(tele_calibra[])))

###################################################################
###########################END#####################################