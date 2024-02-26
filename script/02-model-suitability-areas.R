#...................................................
#...................................................
# Model species distribution using bioclimatic variables
# ...................................................
# ...................................................
# Packages ####
library("rgdal")
library("raster")
library("dismo")
library("rgeos")
library("sf")
library("geodata")
#...................................................
#...................................................
# Data ####
# the BiodiversityR saves outputs on the current working directory
# get the parent wd to return here if needed
wcpath = "data/wc2.1-global"
outputwd = "output/bioclim"
dir.create(outputwd, showWarnings = FALSE, recursive = TRUE)


vif = read.csv("data/vif-bioclim-selection.csv")
vif = vif$bio

# bioclimatic variables
bio = worldclim_global("bio", res = 5, wcpath)
bio = stack(bio)
names(bio) = gsub("wc2.1_5m_", "", names(bio))
bio = subset(bio, subset = vif)
bio = stack(bio)

ext = extent(bio)
ext@ymin = -62

bio = crop(bio, ext)

# read file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

gcm = model_runs[keep, ]

rownames(gcm) = 1:nrow(gcm)

# passport data
dat = read.csv("data/spam-crop-presence.csv")

sp = unique(dat$crop)

sp

#...................................................
#...................................................
# Run ensemble modelling ####
for (i in seq_along(sp)) {
  
  s = dat$crop == sp[i]
  
  # subset data
  coord = dat[s, ]
  
  message("\n Modelling ", 
          sp[i], "\n Time: ", 
          date(), "\n")
  
  coord = coord[, c("x", "y")]
  
  # create background points over the area
  nbg = ceiling(nrow(coord))
  
  bg = randomPoints(bio[[1]], 
                    n = nbg, 
                    p = coord, 
                    ext = ext, 
                    extf = 1.5)
  
  bg = as.data.frame(bg)

  # plot(bio[[1]])
  # points(coord[,c(1,2)])
  # points(bg, pch = "+")

  
  cat("Using ", nrow(coord), " presence points \n")
  
  group1 = kfold(coord, 5)
  group2 = kfold(bg, 5)
  
  bio_layers = c()
  
  mod = list()
  
  thres = c()
  
  for(k in 1:5) {
    
    pres_train = coord[group1 != k, ]
    pres_test = coord[group1 == k, ]
    
    backg_train = bg[group2 != k, ]
    backg_test = bg[group2 == k, ]
    
    bc = bioclim(bio, pres_train)
  
    mod[[k]] = bc
    
    e = evaluate(pres_test, backg_test, bc, bio)
    
    tr = threshold(e, 'spec_sens')
    
    thres = c(thres, tr)
    
    pb = predict(bio, bc, ext=ext, progress = '')
    pb[pb < tr] = 0
    
    bio_layers = c(bio_layers, pb)
    
  }
  
  bio_layers = stack(bio_layers)
  
  bio_layers = calc(bio_layers, max)
  
  writeRaster(bio_layers, 
              filename = paste0(outputwd, "/", sp[i], "-current.tif"),
              overwrite = TRUE)
  
  ### write raster for each gcm model in RCP45 and RCP85 
  for (j in seq_len(nrow(gcm))) {
    
    cat("Predict future distribution, GCM", toupper(gcm[j,]), "\n")
    
    #load GCM layers
    future = cmip6_world(model = gcm[j, 1],
                         ssp = gcm[j, 2],
                         time = gcm[j, 3],
                         var = "bio",
                         res = 5, 
                         path = wcpath)
    
    names(future) = paste0("bio_", 1:19)
    
    future = subset(future, subset = vif)
  
    future = stack(future)
      
    future = crop(future, ext)
    
    bio_future = c()
    
    for(l in seq_along(mod)) {
      
      m = mod[[l]]
      
      pb = predict(future, m, ext=ext, progress = '')
      
      pb[pb < thres[l]] = 0
      
      bio_future = c(bio_future, pb)
      
    }
    
    bio_future = stack(bio_future)
    
    bio_future = calc(bio_future, max)
    
    writeRaster(bio_future, 
                filename = paste0(outputwd, "/", 
                                  sp[i], "-", 
                                  paste(model_runs[j, 1:3], collapse = "-"), 
                                  ".tif"),
                overwrite = TRUE)
    
    
  }
  
}

message("Done at ", Sys.time())

