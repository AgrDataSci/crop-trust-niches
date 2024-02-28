# Combine variables using PCA for historical and future climate, and analyze 
# multi-dimensional shifts in climate conditions.
library("terra")
library("geodata")

wcpath = "data/wc2.1-global"

bio = worldclim_global("bio", res = 5, wcpath)
names(bio) = gsub("wc2.1_5m_", "", names(bio))
names(bio) = gsub("_", "", names(bio))
bio

spam_points = read.csv("output/sampled-points-spam-ecocrop/sampled-points-spam-ecocrop.csv")

head(spam_points)

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = unique(plant_spp$SPAM_Code)

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))

gcm = model_runs[keep, ]

gcm = unique(gcm$V1)

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

ssp = model_runs[keep, c("V2", "V3")]

# extract bioclim current 
dat2 = spam_points[spam_points$suitability == 1 & spam_points$ssp == "126-2041-2060", ]

bio_e = extract(bio, dat2[,c("x", "y")])

dat2 = cbind(dat2, bio_e)

dat2 = dat2[,!grepl("ID", names(dat2))]

head(dat2)

dat2$ssp = "current"

# split the data to extract bioclim by ssp
dat = split(spam_points, spam_points$ssp)

check_point = as.character(unlist(lapply(dat, function(x) unique(x$ssp))))

ssp
check_point

for(i in seq_along(ssp$V2)) {
  
  print(check_point[i])
  
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
    
    b = extract(b, dat[[i]][, c("x", "y")])
    
    bio_future[[j]] = b
  }
  
  # get the mean of gcm
  mean_b = matrix(NA, nrow = nrow(dat[[i]]), ncol = 19)
  
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
  
  dat[[i]] = cbind(dat[[i]], mean_b)
  
}

dat = do.call(rbind, dat)

dat = rbind(dat2, dat)

write.csv(dat, "output/sampled-points-spam-ecocrop/bioclim-extracted.csv", row.names = FALSE)

