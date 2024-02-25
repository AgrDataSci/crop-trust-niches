library("geodata")
library("terra")
library("janitor")

here = "output/ecocrop-raw"

output = "output/ecocrop"

dir.create(output, recursive = TRUE, showWarnings = FALSE)

# read file with crop parameters
plant_spp = read.csv("data/calib-eco-pars.csv")

plant_spp$name2 = make_clean_names(plant_spp$NAME)

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, c("V2", "V3")]

# # ........................................
# # ........................................
# # Average rasters by SSPs #####
# # this will run over GCMs to get the average and over crops 

spam = unique(plant_spp$SPAM_Code)

gcm = paste(model_runs$V2, model_runs$V3, sep = "-")

gcm = c("current", gcm)

for (i in seq_along(spam)) {

  index = spam[i]

  files = list.files(here, pattern = index, full.names = TRUE)

  for (j in seq_along(gcm)) {

    f = grep(gcm[j], files)

    f = files[f]

    r = rast(f)

    r = max(r)
    
    if (isFALSE(spam[i] %in% c("chic", "cowp", "pige", "bean"))) {
      r[r < 0.07] = 0 
    }
    writeRaster(r, filename = paste0(output,"/", index, "-SSP-", gcm[j], ".tif"),
                overwrite = TRUE)

  }

}

