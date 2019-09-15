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

createSpeciesFolders(inDir= "rename", 
                     hasCameraFolders=TRUE,
                     species, 
                     removeFolders = FALSE)


t <- 0
