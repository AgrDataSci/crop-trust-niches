library("BiodiversityR")
library("janitor")
library("geodata")
library("raster")

# https://reagro.org/blocks/when/ecocrop.html
wcpath = "data/wc2.1-global"
outputdir = "output/ecocrop-raw"
dir.create(outputdir, showWarnings = FALSE, recursive = TRUE)

res = 5

# read file with crop parameters
plant_spp = read.csv("data/calib-eco-pars.csv")

plant_spp$name2 = make_clean_names(plant_spp$NAME)

# read file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, ]

rownames(model_runs) = 1:nrow(model_runs)

# read current climate data
current = worldclim_global("bio", res = res, wcpath)
current = stack(current)
names(current) = paste0("bio", 1:19)
current = subset(current, subset = c("bio5", "bio6", "bio12"))
current@title = "base"

#index = grep("cowpea", plant_spp$crop)
# index = 26
# sel = index:nrow(plant_spp)
# 
# plant_spp = plant_spp[sel, ]

# run over species names
for (p in seq_along(plant_spp$NAME)) {
  
  print(plant_spp$name2[p])
  
  # As the raster data correspond to WorldClim version 1,
  # the temperatures need to be multiplied by 10
  temp_thres = c(plant_spp[p, "TMIN"], plant_spp[p, "TOPMN"],
                 plant_spp[p, "TOPMX"], plant_spp[p, "TMAX"])
  
  
  rain_thres = c(plant_spp[p, "RMIN"], plant_spp[p, "ROPMN"],
                 plant_spp[p, "ROPMX"], plant_spp[p, "RMAX"])
  
  
  ecocrop_pars <- ensemble.ecocrop.object(temp.thresholds = temp_thres, 
                                          rain.thresholds = rain_thres,
                                          temp.multiply = 1,
                                          annual.temps = FALSE, 
                                          name = plant_spp$name2[p])
  e = ensemble.ecocrop(current, 
                       ecocrop.object = ecocrop_pars,
                       RASTER.stack.name = "base")
  
  writeRaster(e, 
              filename = paste0(outputdir, "/", plant_spp$name2[p], "-current.tif"),
              overwrite = TRUE)
  
  # now run over climate scenarios 
  for (i in seq_len(nrow(model_runs))) {
    
    print(paste(plant_spp$NAME[p], paste(model_runs[i, 1:3], collapse = " - ")))
    
    future = cmip6_world(model = model_runs[i, 1],
                         ssp = model_runs[i, 2],
                         time = model_runs[i, 3],
                         var = "bio",
                         res = res, 
                         path = wcpath)
    future = stack(future)
    future = subset(future, subset = c("bio05", "bio06", "bio12"))
    future@title = "base"
    
    e_future = ensemble.ecocrop(future, 
                                ecocrop.object = ecocrop_pars,
                                RASTER.stack.name = "base")
    
    writeRaster(e_future, 
                filename = paste0(outputdir, "/", 
                                  plant_spp$name2[p], "-", 
                                  paste(model_runs[i, 1:3], collapse = "-"), 
                                  ".tif"),
                overwrite = TRUE)
    
  }
}
