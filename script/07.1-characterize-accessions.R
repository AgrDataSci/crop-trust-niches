# Per crop and group (landrace group, or CWR genepool) produce a table
# with the average, sigma, CV, q1, q2, q3, q4, min/max, of the 
# historical conditions of the accessions.
library("terra")
library("geodata")

# ....................................
# ....................................
# Input data #####
# read file with crop parameters
plant_spp = read.csv("data/calib-eco-pars.csv")

sel = list.files("data/cwr-processed", pattern = ".csv")
sel = gsub(".csv", "", sel)

crop_name = plant_spp$NAME

crop_name = crop_name[crop_name %in% sel ]

plant_spp = plant_spp[plant_spp$NAME %in% sel, ]

plant_spp = plant_spp[!duplicated(plant_spp$NAME), ]

# read crop wild relatives data
# run over crop files
cwr = list()

for(i in seq_along(crop_name)) {
  cwr[[i]] = read.csv(paste0("data/cwr-processed/", crop_name[i], ".csv"))
}

cwr = do.call(rbind, cwr)

cwr$x = as.numeric(cwr$final_lon)

cwr$y = as.numeric(cwr$final_lat)

cwr$NAME = cwr$Crop_Common_Name

keep = !is.na(cwr$x) & !is.na(cwr$y)

table(keep)

cwr = cwr[keep, ]

cwr = merge(cwr, 
            plant_spp[,c("NAME", "SPAM_Name", "SPAM_Code")],
            by = "NAME", 
            all.x = TRUE)

cwr$GR = "CWR"


table(cwr$SPAM_Name)
table(cwr$NAME)

# ggplot(cwr, aes(x = x, y = y, color = SPAM_Code)) +
#   geom_point()


# read landrace data
# run over crop files
landrace = list()

for(i in seq_along(crop_name)) {
  landrace[[i]] = read.csv(paste0("data/landrace-processed/", crop_name[i], ".csv"))
}

landrace = do.call(rbind, landrace)

names(landrace)

landrace$x = as.numeric(landrace$Longitude)

landrace$y = as.numeric(landrace$Latitude)

landrace$NAME = landrace$Crop_Name

table(landrace$NAME)

keep = !is.na(landrace$x) & !is.na(landrace$y)

table(keep)

landrace = merge(landrace, 
                 plant_spp[,c("NAME", "SPAM_Name", "SPAM_Code")], 
                 by = "NAME",
                 all.x = TRUE)


table(landrace$SPAM_Name)
table(landrace$NAME)

landrace$GR = "Landrace"

dat = rbind(cwr[,c("GR", "NAME", "SPAM_Name", "SPAM_Code","x","y")],
            landrace[,c("GR", "NAME", "SPAM_Name", "SPAM_Code","x","y")])

head(dat)

# ....................................
# ....................................
# Climate variables and SSP #####
wcpath = "data/wc2.1-global"

bio = worldclim_global("bio", res = 5, wcpath)
names(bio) = gsub("wc2.1_5m_", "", names(bio))

model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

gcm = model_runs[keep, ]

gcm = unique(gcm$V1)

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

ssp = model_runs[keep, c("V2", "V3")]

# ....................................
# ....................................
# Extract climate data #####
# read with file with models to run 
# extract bioclim current 
bio_e = extract(bio, dat[,c("x", "y")])

bio_e = bio_e[, -grep("ID", names(bio_e))]

dat2 = cbind(dat, bio_e)

head(dat2)

dat2$ssp = "current"

bio_f = list()

for(i in seq_along(ssp$V2)) {
  
  print(ssp[i, ])
  
  bio_future = list()
  
  # run over gcm layers and extract values
  for (j in seq_along(gcm)) {
    
    b = cmip6_world(model = gcm[j],
                    ssp = ssp[i, 1],
                    time = ssp[i, 2],
                    var = "bio",
                    res = 5, 
                    path = wcpath)
    
    names(b) = paste0("bio", 1:19)
    
    b = extract(b, dat[, c("x", "y")])
    
    bio_future[[j]] = b
  }
  
  # get the mean of gcm
  mean_b = matrix(NA, nrow = nrow(dat), ncol = 19)
  
  for(k in seq_along(1:19)) {
    v = lapply(bio_future, function(x){
      x[,paste0("bio", k)]
    })
    v = as.data.frame(do.call(cbind, v))
    v = rowMeans(v)
    mean_b[, k] = v
  }
  
  mean_b = as.data.frame(mean_b)
  
  names(mean_b) = paste0("bio", 1:19)
  
  mean_b$ssp = paste0(ssp[i, ], collapse = "-")
  
  bio_f[[i]] = cbind(dat, mean_b)
  
}

bio_f = do.call("rbind", bio_f)

names(bio_f) = gsub("bio_", "bio", names(bio_f))

names(dat2) = gsub("bio_", "bio", names(dat2))

names(bio_f)

names(dat2)

dat = rbind(dat2, bio_f)

write.csv(dat, "output/sampled-points-spam-ecocrop/landrace-and-cwr-bioclim-extracted.csv", 
          row.names = FALSE)





