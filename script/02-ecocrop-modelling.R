library("geodata")
library("Recocrop")
library("terra")
library("janitor")

outputdir = "output/ecocrop-raw"

dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

# https://reagro.org/blocks/when/ecocrop.html
wcpath = "data/wc2.1-global"

res = 5

# read file with crop parameters
plant_spp = read.csv("data/calib-eco-pars.csv")

plant_spp$name2 = make_clean_names(plant_spp$NAME)

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, ]

# read current climate data
rain = worldclim_global("prec", res = res, wcpath)
names(rain) = paste0("rain", 1:length(names(rain)))

tmin = worldclim_global("tmin", res = res, wcpath)
names(tmin) = paste0("tmin", 1:length(names(tmin)))

tmax = worldclim_global("tmax", res = res, wcpath)
names(tmax) = paste0("tmax", 1:length(names(tmax)))

tavg = mosaic(tmin, tmax, fun = "mean")
names(tavg) = paste0("tavg", 1:length(names(tavg)))


for (p in seq_along(plant_spp$NAME)) {
  
  print(plant_spp$NAME[p])
  
  crop = ecocropPars(plant_spp$NAME[p])
  
  m = ecocrop(crop)
  
  control(m, get_max = TRUE)
  
  plant = predict(m, 
                  tavg = tavg, 
                  prec = rain, 
                  filename = paste0(outputdir, "/", plant_spp$name2[p], "-current.tif"),
                  overwrite = TRUE)
  
  # now run over climate scenarios 
  for (i in seq_len(nrow(model_runs))) {
    
    print(paste(plant_spp$NAME[p], paste(model_runs[i, 1:3], collapse = " - ")))
    
    rain_future = cmip6_world(model = model_runs[i, 1],
                              ssp = model_runs[i, 2],
                              time = model_runs[i, 3],
                              var = "prec",
                              res = res, 
                              path = wcpath)
    names(rain_future) = paste0("rain", 1:length(names(rain_future)))
    
    tmax_future = cmip6_world(model = model_runs[i, 1],
                              ssp = model_runs[i, 2],
                              time = model_runs[i, 3],
                              var = "tmax",
                              res = res, 
                              path = wcpath)
    names(tmax_future) = paste0("tmax", 1:length(names(tmax_future)))
    
    tmin_future = cmip6_world(model = model_runs[i, 1],
                              ssp = model_runs[i, 2],
                              time = model_runs[i, 3],
                              var = "tmin",
                              res = res, 
                              path = wcpath)
    names(tmin_future) = paste0("tmin", 1:length(names(tmin_future)))
    
    tavg_future = mosaic(tmin_future, tmax_future, fun = "mean")
    names(tavg_future) = paste0("tavg", 1:length(names(tavg_future)))
    
    plant_future = predict(m, 
                           tavg = tavg_future, 
                           prec = rain_future,
                           filename = paste0(outputdir, "/", 
                                             plant_spp$name2[p], 
                                             "-", paste(model_runs[i, 1:3], collapse = "-"), 
                                             ".tif"),
                           overwrite = TRUE)
    
  }
  
}


