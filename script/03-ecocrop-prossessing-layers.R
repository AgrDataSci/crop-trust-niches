library("geodata")
library("terra")
library("janitor")

here = "output/ecocrop-raw/"

output = "output/ecocrop/"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# read file with crop parameters
plant_spp = read.csv("data/calib-eco-pars.csv")

plant_spp$name2 = make_clean_names(plant_spp$NAME)

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, c("V2", "V3")]

# ........................................
# ........................................
# Average rasters by SSPs #####
# this will run over GCMs to get the average and over crops 
crop_spp = unique(plant_spp$crop)

gcm = paste(model_runs$V2, model_runs$V3, sep = "-")

for (i in seq_along(plant_spp$name2)) {
  
  index = plant_spp$name2[i]
  
  files = list.files(here, pattern = index, full.names = TRUE)
  
  for (j in seq_along(gcm)) {
    
    f = grep(gcm[j], files)
    
    f = files[f]
    
    r = rast(f)
    
    r = mean(r)
    
    writeRaster(r, filename = paste0(here, index, "-SSP-", gcm[j], ".tif"))
    
  }
  
}

# ........................................
# ........................................
# Max value per SPAM crop ####
# get the max value in each cell combining ecocrop plant species
# into a single SPAM crop 
spam = unique(plant_spp$SPAM_Code)

ssp = c("current", paste0("SSP-", gcm))

for (i in seq_along(spam)) {
  
  index = plant_spp[plant_spp$SPAM_Code == spam[i], "name2"]
  
  index1 = paste(index, collapse = "|")
  
  files = list.files(here, pattern = index1, full.names = TRUE)
  
  for (j in seq_along(ssp)) {
    
    print(paste(spam[i], ssp[j]))
    
    f = grep(ssp[j], files)
    
    f = files[f]
    
    r = rast(f)
    
    r = max(r)
    
    names(r) = gsub(" |-", "", paste(spam[i], ssp[j]))
    
    writeRaster(r, filename = paste0(output, spam[i], "-", ssp[j], ".tif"), overwrite = TRUE)
  }
}

list.files(output)



files = list.files(output, pattern = "maize", full.names = T)

r = rast(files)

plot(r)









