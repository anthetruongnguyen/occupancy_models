
################################## CAMERA-TRAP WORKFLOW ##########################################
library(camtrapR)
wd <- "D:/IZW/Occupancy training/cameratrap photos/Pu Mat NP"
setwd(wd)

### create raw folders###
createStationFolders(inDir = "raw", 
                     stations = c("stc001", "stc002"),
                     createinDir = TRUE)

### rename photos ###
imageRename(inDir = "raw", 
            outDir = "rename", 
            hasCameraFolders = TRUE,
            keepCameraSubfolders = TRUE,
            createEmptyDirectories = FALSE,
            copyImages = TRUE, 
            writecsv = FALSE)

### create species folders ###
species <- c("Annam partridge", "Annamite striped rabbit", "Asiatic brush tailed porcupine", "Asiatic black bear", "Assamese Macaque",
             "Bar backed partridge", "Bar bellied pitta", "Black drongo","Black Giant Squirrel", "Black hooded laughingthrush", 
             "Blank", 
             "Blue pitta", "Blue rumped pitta", "Blue whistling thrush", "Buff-breasted babbler", 
             "Common palm civet", "Crab eating mongoose", "Crested argus","Crested goshawk", "Crested serpent eagle", 
             "Dark muntjac", 
             "Domestic buffalo", "Domestic cow", "Domestic dog", "Domestic fowl", "Domestic pig", 
             "Emerald dove", "Wild pig", "Fairy pitta", "Ferret badger", 
             "Great eared nightjar", "Grey peacock pheasant",
             "Large scimitar babbler", "Leopard cat", "Lesser necklaced laughingthrush",
             "Lesser oriental chevrotain",
             "Local people", 
             "Malayan night heron", "East Asian porcupine", "Masked palm civet", 
             "Mountian hawk eagle", 
             "Murid", 
             "Monitor lizard",
             "Northern treeshrew", "Orange headed thrush", "Owston's civet", "Pangolin", "Pig tailed macaque", 
             "Puff-throated babbler", "Racket tailed treepie", 
             "Ranger", 
             "Red collared woodpecker","Red junglefowl", "Red muntjac", "Red shanked douc langur", 
             "Researcher", "Retrieval", 
             "Rhesus macaque", "Rufous throated fulvetta", "Rufous cheeked laughingthrush", "Rufous throated partridge", "Rufous-tailed robin", 
             "Sambar", "Scaly thrush", "Serow", 
             "Setting", 
             "Short-tailed Scimitar babbler","Siamese fireback", "Siberian blue robin", 
             "Siberian thrush", "Silver pheasant", "Slaty legged crake", "Spot necked babbler", 
             "Spotted linsang", 
             "Squirrel",
             "Stripe backed Weasel", "Stump tailed macaque",
             "Unidentified animal","Unidentified bat", "Unidentified bird", 
             "Unidentified macaque", "Unidentified mammal", "Unidentified muntjac", "Unidentified pheasant", "unidentified weasel",
             "White breasted waterhen", "White crested laughingthrush", "White rumped sharma", 
             "White winged magpie", "Yellow bellied weasel", "Yellow throated marten", "Impressed tortoise", "Slender tailed treeshrew", "Red cheeked Squirrel", "Hairy footed Squirrel", 
             "Unidentified flying squirrel",
             "Sunda pangolin")

createSpeciesFolders(inDir= "renamed", 
                     hasCameraFolders=TRUE,
                     species, 
                     removeFolders = FALSE)

###append species names###
wd_SBC <- "D:/Dropbox (ScreenForBio)/Field Data/Vietnam 2017-2018/Camera Trap Photos/SBC_grid 2"
setwd(wd_SBC)
wd_images_renamed <- "renamed"
appendSpeciesNames(inDir = wd_images_renamed, 
                   IDfrom = "directory",
                   hasCameraFolders = TRUE,
                   #metadataSpeciesTag,
                   #metadataHierarchyDelimitor = "|",
                   removeNames = FALSE
                   #writecsv = FALSE
)
###record data###

exclude = c("Blank", "Ranger", "Researcher", "Retrieval", "Setting", 
            "Unidentified animal", "Unidentified bird", 
            "Unidentified macaque", "Unidentified muntjac", "Unidentified pheasant",
            "Unidentified flying squirrel")
inDir = "renamed"
outDir = wd_SBC

rc.table <- recordTable(inDir,
                        IDfrom = "directory",
                        cameraID = "filename",
                        camerasIndependent = FALSE,
                        minDeltaTime = 60,
                        deltaTimeComparedTo = "lastRecord",
                        timeZone = "Asia/Saigon",
                        stationCol = "station",
                        exclude,
                        writecsv = TRUE,
                        outDir, 
                        removeDuplicateRecords = TRUE
)

###create detection plots and shapfiles###
library(readxl)
CTtable <- read_excel("datatable_Nui Chua_grid 2.xlsx", sheet = 'station info')
CTtable <- data.frame(CTtable)
recordTable <- read.csv("record_table_60min_deltaT_2019-09-15.csv")
recordTable <- data.frame(recordTable)

#recordTable$Species[recordTable$Species == "Red cheeked Squirrel"] <- "Squirrel" ## repace red-cheeked squirrel as squirrel

det.map <- detectionMaps(CTtable,
                         recordTable = recordTable,
                         Xcol  = "X",
                         Ycol = "Y",
                         stationCol = "station",
                         speciesCol = "Species",
                         richnessPlot = TRUE,
                         speciesPlots = TRUE,
                         addLegend = TRUE,
                         printLabels = TRUE,
                         smallPoints = 0,
                         plotR = TRUE,
                         writePNG = TRUE,
                         plotDirectory = "detection",
                         createPlotDir = TRUE,
                         pngMaxPix = 1000,
                         writeShapefile = TRUE,
                         shapefileName = "detection_map",
                         shapefileDirectory = "shapefile",
                         shapefileProjection = "+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
)

#### camera trap operation###

Cam.op <- cameraOperation(CTtable, 
                          stationCol = "station", 
                          setupCol = "date_setting", 
                          retrievalCo = "date_retrieval", 
                          hasProblems = FALSE,
                          dateFormat = "%Y-%m-%d", 
                          writecsv = FALSE,
                          outDir = wd_SBC
)

### survey report###
report <- surveyReport (recordTable,
                        CTtable,
                        speciesCol           = "Species",
                        stationCol           = "station",
                        setupCol             = "date_setting",
                        retrievalCol         = "date_retrieval",
                        CTDateFormat         = "%Y-%m-%d", 
                        recordDateTimeCol    = "DateTimeOriginal",
                        recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                        sinkpath             = getwd())

### activity ###
acti.sp1 <- activityDensity(recordTable, 
                species = "Common palm civet",
                allSpecies = FALSE,
                speciesCol = "Species",
                recordDateTimeCol = "DateTimeOriginal",
                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                plotR = TRUE, 
                writePNG = FALSE, 
                plotDirectory, 
                createDir = FALSE, 
                pngMaxPix = 1000,
                add.rug = TRUE
)

### activity overlap ###
acti.sp1.sp2 <- activityOverlap(recordTable, 
                                speciesA = "Common palm civet",
                                speciesB = "Ferret badger",
                                speciesCol = "Species",
                                recordDateTimeCol = "DateTimeOriginal",
                                recordDateTimeFormat = "%Y-%m-%d %H:%M:%S",
                                plotR = TRUE, 
                                writePNG = FALSE, 
                                addLegend = TRUE,
                                legendPosition = "topleft",
                                plotDirectory, 
                                createDir = FALSE, 
                                pngMaxPix = 1000,
                                add.rug = TRUE
                               
)


################################## OCCUPANCY MODELS ##############################################

### get covariates from GIS raster ###

library(raster)
library(rgdal)
library(rgeos)
library(ggplot2)
library(reshape2)

# get points from CT table #
CTtable <- read.csv(".csv")
CTtable <- as.data.frame(CTtable)
#CTtable <- as.data.frame(CTtable)
ct_coords <- data.frame(CTtable$UTM_N, CTtable$UTM_E)
ct_points <- SpatialPoints(ct_coords) 
ct_points <- SpatialPointsDataFrame(ct_coords,
                                    data = CTtable[,seq(1, 15)], # add information (get from CTtable) to points 
                                    proj4string=  CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")) 


# load raster #
wd_raster <- "D:/Dropbox (ScreenForBio)/Projects/Vietnam 2017-2019/SBC project/Occupancy testing"
setwd(wd_raster)
elevation <- raster("elevation.tif")
fs <- raster("forest_score.tif")
d_roads <- raster("roads_distance_30m.tif")
#covariates_stacked <- stack(d_roads,d_water,d_villages)
#names(covariates_stacked) <- c("d_roads", "d_water","d_villages")

# extract values/covariates from rasters/stracked raster #
covs <- as.data.frame(extract(covariates_stacked,ct_points))
covs_table <- cbind(CTtable[,seq(1, 2)],   #have no clue why CTtable$station doesnt work
                    d_roads=covs$d_roads,
                    d_water=covs$d_water,
                    d_villages=covs$d_villages)
write.csv(covs_table, file="gis_covariates.csv")

### bulding models ###
# load data
rec.table <- read.csv("recordtable_test.csv")
rec.table$DateTimeOriginal   <- strptime(paste(rec.table$Date, rec.table$Time), format = "%m/%d/%Y %H:%M", tz = "UTC")
colnames(rec.table)[colnames(rec.table) == "station"] <- "Station"
#colnames(rec.table)[colnames(rec.table) == "species"] <- "Species"

CTtable <- read.csv("LGXM_CT_infoDot1.csv")
colnames(CTtable)[colnames(CTtable) == "station"] <- "Station"
#date.format.old <- "%m/%d/%Y"  
#CTtable$Date_set     <- as.Date(CTtable$Date_setup , format = date.format.old)
#CTtable$Date_retrieved <- as.Date(CTtable$Date_retrieval, format = date.format.old)

covs_gis <- read.csv("gis_covariates.csv")

# scale numeric values

covtable <- cbind(d_roads=covs_gis$d_roads,d_water=covs_gis$d_water,d_villages=covs_gis$d_villages)
covtable_scaled <- scale(covtable)
covtable_scaled <- as.data.frame(covtable_scaled)


# camera trap operation
camop2 <- cameraOperation(CTtable, 
                          stationCol = "Station",
                          allCamsOn = TRUE,
                          setupCol = "date_setup", 
                          retrievalCol = "date_retrieval", 
                          hasProblems = FALSE,
                          dateFormat = "%Y-%m-%d",
                          camerasIndependent = TRUE)

#### RUN MODELS ###
detHist1 <- detectionHistory(recordTable       = rec.table,
                             camOp             = camop2,
                             stationCol        = "Station",
                             speciesCol        = "Species",
                             recordDateTimeCol = "DateTimeOriginal",
                             occasionLength    = 10,
                             minActiveDaysPerOccasion = 0, ########  IMPORTANT
                             species           = "Ferret badger",
                             includeEffort     = TRUE, ########  NEED DOUBLE CHECK
                             scaleEffort       = FALSE,
                             day1              = "station",
                             timeZone          = "Asian/Saigon")


uf_occu1 <- unmarkedFrameOccu( y        = detHist1$detection_history,
                               siteCovs = covtable_scaled,
                               obsCovs  = list(effort = as.data.frame(detHist1$effort))
)

m1 <- occu(~ effort ~d_roads, data=uf_occu1)  
m2 <- occu(~ effort ~d_water, data=uf_occu1)
m3 <- occu(~ effort ~d_villages, data=uf_occu1)
m4 <- occu(~ effort ~d_roads+d_water+d_villages, data=uf_occu1)

##### CREATE PREDICTION MAP ####

#load rasters
elevation <- raster("elevation.tif")
fs <- raster("forest_score.tif")
d_roads <- raster("roads_distance_30m.tif")

#get covariates value for prediction (e.g. d_roads)
dat4plot <- data.frame(d_roads=scale(values(d_roads)),
                       d_water=scale(values(d_water)),
                       d_villages=scale(values(d_villages)))

# predict occupancy 
predict_m1 <- predict(m1, newdata=dat4plot, type="state")
#predict_m2 <- predict(m2, newdata=dat4plot, type="state")
#predict_m3 <- predict(m3, newdata=dat4plot, type="state")
#predict_m4 <- predict(m4, newdata=dat4plot, type="state")

# create empty prediction rasters and fill it with predicted values
prediction_raster_m1 <- raster(d_roads)
#prediction_raster_m2 <- raster(d_roads)
#prediction_raster_m3 <- raster(d_roads)
#prediction_raster_m4 <- raster(d_roads)

values(prediction_raster_m1) <- predict_m1$Predicted
#values(prediction_raster_m2) <- predict_m2$Predicted
#values(prediction_raster_m3) <- predict_m3$Predicted
#values(prediction_raster_m4) <- predict_m4$Predicted

#writeRaster(prediction_raster_m1, file="prediction_map_with_effort_d_roads", overwrite=TRUE)
