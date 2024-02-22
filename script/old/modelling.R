#...................................................
#...................................................
# Model species distribution using bioclimatic variables
# ...................................................
# ...................................................
# Packages ####
library("data.table")
library("rgdal")
library("raster")
library("dismo")
library("rgeos")
library("gbm")
library("glmnet")
library("maxnet")
library("sf")
library("BiodiversityR")
library("PresenceAbsence")
library("maxlike")
library("geodata")
#...................................................
#...................................................
# Data ####
# the BiodiversityR saves outputs on the current working directory
# get the parent wd to return here if needed
parentwd = getwd()
wcpath = "data/wc2.1-global"
outputwd = "output/enm/"
dir.create(outputwd, showWarnings = FALSE, recursive = TRUE)


vif = read.csv("data/vif-bioclim-selection.csv")
vif = vif$bio

# bioclimatic variables
bio = worldclim_global("bio", res = 5, wcpath)
bio = stack(bio)
names(bio) = gsub("wc2.1_5m_", "", names(bio))
bio = subset(bio, subset = vif)
bio = stack(bio)

# define projection and extension
myproj = proj4string(bio)
myext  = extent(bio) 
myres  = res(bio)

# read file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

gcm = model_runs[keep, ]

rownames(gcm) = 1:nrow(gcm)

# passport data
dt = fread("data/spam-crop-presence.csv")

sp = unique(dt$crop)
sp

#...................................................
#...................................................
# In case you need to stop the process here is the 
# point where the algorithm looks for species already 
# processed
filepattern = paste(c("current"), collapse = "|")
nfiles = length(c("current")) * 2

for (i in seq_along(sp)) {
  
  pres = paste0(outputwd, sp[i], "/ensembles/presence/")
  pres = grepl(filepattern, list.files(pres))
  
  suit = paste0(outputwd, sp[i], "/ensembles/suitability/")
  suit = grepl(filepattern, list.files(suit))
  
  # if folder doesnt contain these files then it should be removed
  if (sum(suit) != nfiles | sum(pres) != nfiles) {
    
    unlink(paste0(outputwd, sp[i]), recursive = TRUE)
    
  }
  
}

done = list.dirs(outputwd)[-1]
done = strsplit(done, "/")
done = suppressWarnings(do.call("rbind", done)[,4])
done = unique(done)


# filter and run the model for the remnant 
sp = sp[!sp %in% done] 

#...................................................
#...................................................
# Run ensemble modelling ####
for (i in seq_along(sp)) {
  
  # create a dir for the species and work in that dir
  output = paste0(outputwd, sp[i])
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  #setwd(output)
  
  k = dt$crop == sp[i]
  
  # subset data
  coord = dt[k, ]
  
  message("\n Ensemble modelling for ", 
          sp[i], "\n Time: ", 
          date(), "\n")
  
  coord = coord[, .(x, y)]
  
  coord = as.data.frame(coord)
  
  # coord = coord[coord$x > 0, ]
  # 
  # coord = coord[coord$y < 20, ]
  
  # create a convexHull to limit the model to the 
  # area where the species is actually present
  # calculate largest distance
  largedist = pointDistance(coord, longlat = FALSE)
  largedist = max(largedist, na.rm = TRUE)
  
  # make a convex hull 
  hull = convHull(coord, lonlat = TRUE)
  
  # extent convex hull
  hull = gBuffer(hull@polygons, 
                  width = 0.1 * largedist)
  
  crs(hull) = myproj
  
  # crop bioclim layers to fit this extention
  bio_i = mask(bio, hull)
  
  bio_i = crop(bio_i, hull)
  
  bio_i = stack(bio_i)
  
  ext = extent(bio_i)
  
  # create a buffer with presence points to avoid background
  # points too close to presence
  lonlatb = st_as_sf(coord, 
                      coords = c("x", "y"), 
                      crs = 4326)
  
  # set the buffer around the points
  lonlatb = suppressWarnings(st_buffer(lonlatb,
                                        dist = 1))
  
  lonlatb = st_union(lonlatb)
  
  lonlatb = as_Spatial(lonlatb)
  
  lonlatb = rasterize(lonlatb, bio_i[[1]], field = 1, background = NA)
  
  bg_mask = mask(bio_i[[1]], lonlatb, inverse = TRUE)
  
  # create background points over the area
  nbg = ceiling(nrow(coord))
  bg = randomPoints(bg_mask, 
                     n = nbg, 
                     p = coord, 
                     ext = ext, 
                     extf = 1.25)
  bg = as.data.frame(bg)

  plot(bio_i[[1]])
  points(coord[,c(1,2)])
  points(bg, pch = "+")
  
  # reduce sampling bias removing points within the same grid cell
  r = raster(ext)
  res(r) = myres
  
  coord = as.data.frame(coord)
  
  coord = gridSample(coord, r, n = 1)
  
  cat("Using ", nrow(coord), " presence points \n")
  
  
  group1 <- kfold(coord, 5)
  group2 <- kfold(bg, 5)
  
  bio_layers = c()
  
  for(k in 1:5) {
    
    pres_train <- coord[group1 != k, ]
    pres_test <- coord[group1 == k, ]
    
    backg_train <- bg[group2 != k, ]
    backg_test <- bg[group2 == k, ]
    
    # plot(bio_i[[1]])
    # points(backg_train, pch='-', cex=0.5, col='yellow')
    # points(backg_test, pch='-',  cex=0.5, col='black')
    # points(pres_train, pch= '+', col='green')
    # points(pres_test, pch='+', col='blue')
    
    bc <- bioclim(bio_i, pres_train)
  
    mod[[k]] = bc
    
    tr <- threshold(e, 'spec_sens')
    
    pb <- predict(bio_i, bc, ext=ext, progress='')
    pb[pb < tr] = 0
    
    mod = c(mod, pb)
    
  }
  
  mod = stack(mod)
  
  mod = calc(mod, max)
  
  m = mod
  
  m[m>0] = 1
  
  plot(m)
  
  # par(mfrow=c(1,2))
  # plot(pb, main='Bioclim, raw values')
  # plot(pb > tr, main='presence/absence')
  # points(pres_train, pch='+')
  # 
  
  
  
  
  # Run ensemble modelling
  # step 1: model calibration
  # here the function tests for the best algorithms
  # since the algorithms were previously selected,
  # a 3-fold cross-validation is performed to make sure that all 
  # pass the output.weights threshold
  message("\n Step 1: Calibrating algorithms \n", "Time: ", date(), "\n")
  set.seed(9999)
  enm_step1 = ensemble.calibrate.weights(x = bio_i, 
                                          p = coord, 
                                          a = bg,
                                          k = 3,
                                          layer.drops = NULL,
                                          SINK = TRUE, 
                                          species.name = sp[[i]],
                                          BIOCLIM = 1,
                                          DOMAIN = 1, 
                                          MAXNET = 1,
                                          MAHAL = 0,
                                          GBM = 0,
                                          GAM = 0,
                                          GLM = 0,
                                          SVM = 0,
                                          RPART = 0, 
                                          GBMSTEP = 0,
                                          MAXENT = 0, 
                                          NNET = 0, 
                                          RF = 0, 
                                          EARTH = 0,
                                          GLMSTEP = 0, 
                                          GAMSTEP = 0, 
                                          MGCV = 0, 
                                          MGCVFIX = 0, 
                                          CF = 0, 
                                          FDA = 0,
                                          SVME = 0,
                                          ENSEMBLE.tune = TRUE, 
                                          PROBIT = TRUE,
                                          # see Liu et al (2013) doi:10.1111/jbi.12058
                                          threshold.method = "threshold2013.mean", 
                                          threshold.PresenceAbsence = TRUE, 
                                          ENSEMBLE.min = 0.7)
  
  # step 2: create models that will be used for the raster predictions
  # models with output.weights <0.05 are excluded
  output_weights = enm_step1$output.weights
  output_weights[output_weights < 0.05] = 0
  
  message("Step 2: Model species distribution with selected ENM algorithms \n")
  
  set.seed(9999)
  enm_step2 = ensemble.calibrate.models(x = bio_i, 
                                         p = coord, 
                                         a = bg,
                                         k = 10, 
                                         layer.drops = NULL,
                                         SINK = TRUE, 
                                         species.name = sp[[i]],
                                         models.keep = TRUE,
                                         input.weights = output_weights,
                                         # see Liu et al (2013) doi:10.1111/jbi.12058
                                         threshold.method = "threshold2013.mean", 
                                         threshold.PresenceAbsence = TRUE, 
                                         ENSEMBLE.tune = FALSE, 
                                         ENSEMBLE.min = 0.7,
                                         PROBIT = TRUE,
                                         Yweights = "BIOMOD", 
                                         models.save = FALSE)
  
  
  # save AUCs
  auc = data.frame(auc = enm_step2$AUC.testing)
  auc$model = rownames(auc)
  auc = auc[!grepl("MAHAL|ENSEMBLE", rownames(auc)), ]
  write.csv(auc, file = "auc_testing.csv")
  
  message("Step 3.1: Generate map of current distribution \n")
  #step3: use previously calibrated models to construct consensus layers
  ensemble_current = ensemble.raster(xn = bio_i,
                                      models.list = enm_step2$models,
                                      input.weights = output_weights,
                                      thresholds = enm_step2$models$thresholds,
                                      SINK = TRUE,
                                      RASTER.species.name = sp[[i]], 
                                      RASTER.stack.name = "current")
  
  ### write raster for each gcm model in RCP45 and RCP85 
  for (k in seq_len(nrow(gcm))) {
    cat("Step 3.2: Predict future distribution, GCM", toupper(gcm[k]), "\n")
    
    #load GCM layers
    future = cmip6_world(model = gcm[i, 1],
                         ssp = gcm[i, 2],
                         time = gcm[i, 3],
                         var = "bio",
                         res = res, 
                         path = paste0(parentwd, wcpath))
    
    names(future) = paste0("bio_", 1:19)
    
    future = subset(future, subset = vif)
  
    future = stack(future)
      
    #gcm_model <- crop(gcm_model, mesoam)
    
    #gcm_model <- stack(gcm_model)
    
    ensemble_gcm_model <- ensemble.raster(xn = future,
                                          models.list = enm_step2$models,
                                          input.weights = output_weights,
                                          thresholds = enm_step2$models$thresholds,
                                          SINK = TRUE,
                                          RASTER.species.name = species[i], 
                                          RASTER.stack.name = gcm[k])
    
  }
  
  # remove working files created in the third step
  unlink("models", recursive = TRUE)
  unlink("ensembles/count", recursive = TRUE)
  file.remove(list.files(pattern = "working", full.names = TRUE))
  
  # return to parent wd
  setwd(parentwd)
  
}

message("Done at ", Sys.time())

