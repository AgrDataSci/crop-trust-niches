# ..................................
# ..................................
# Draw current production sample from MapSPAM
# using relative harvested area as weight, 
# then remove any pixels that according to the suitability layers become unsuitable under 
# future climates, and add pixels that are projected to become suitable. The added pixels will 
# also be sampled, aiming to maintain the proportion between new and currently cultivated areas, 
# and the sampling weights will be the future projected suitability values. 
# first run: Feb 2024
# Kauê de Sousa CGIAR
# ..................................
# ..................................
library("terra")

suitability_path = "output/change-suitability/"

ecocrop_path = "output/ecocrop/"

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = toupper(unique(plant_spp$SPAM_Code))

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, c("V2", "V3")]

ssp = paste(model_runs$V2, model_runs$V3, sep = "-")

# get ssp layer extention 
myext = ext(rast(paste0(suitability_path, tolower(spam[1]), "-SSP-", ssp[1], ".tif")))

dat = data.frame()

set.seed(2148)
seeds = floor(runif(length(spam), 0, 10000))

for(i in seq_along(spam)) {
  
  for(j in seq_along(ssp)) {
    
    print(paste0(spam[i], " SSP ", ssp[j] ))
    
    # open ecocrop layer
    ecocrop = rast(paste0(ecocrop_path, tolower(spam[i]), "-SSP-", ssp[j], ".tif"))
    
    ecocrop = crop(ecocrop, myext)
    
    # open the map with changes in suitability
    suit = rast(paste0(suitability_path, tolower(spam[i]), "-SSP-", ssp[j], ".tif"))
    
    # open SPAM harvest layer
    s_harv = rast(paste0("data/spam/spam2010V2r0_global_H_", spam[i], "_H.tif"))
    
    s_harv = crop(s_harv, myext)
    
    spam_p = s_harv
    
    # reclassify values to fit in the current SPAM map 
    spam_p[spam_p > 0] = 10
    
    overlay = function(x, y) {(x + y)}
    change_suit = do.call(overlay, list(x = suit, y = spam_p))
    
    change_suit[change_suit==9] = -1
    
    change_suit[change_suit==10] = 1
    
    change_suit[change_suit==11] = 1
    
    change_suit[change_suit==12] = 2
    
    # get sample using harvest area as weight
    set.seed(seeds[i])
    s_sample = spatSample(s_harv, 
                          size = 1000, 
                          method = "weights", 
                          xy = TRUE,
                          na.rm = TRUE)
    
    names(s_sample)[3] = "spam_harvest"
    
    s_sample$crop = tolower(spam[i])
    
    values = extract(change_suit, s_sample[,c(1:2)])
    
    s_sample$suitability = values[,2]
    
    # use new areas as mask as get probabilities from ecocrop
    # to add novel areas
    novel = change_suit
    novel[novel != 2] = NA
    
    novel = mask(ecocrop, novel)
    
    set.seed(seeds[i])
    novel = spatSample(novel, 
                       size = 1000, 
                       method = "weights", 
                       xy = TRUE,
                       na.rm = TRUE)
    
    novel$crop = tolower(spam[i])
    
    novel$suitability = 2
    
    novel = as.data.frame(novel[,-3])
    
    s_sample = as.data.frame(s_sample[,-3])
    
    s_sample = rbind(s_sample, novel)
    
    s_sample = s_sample[s_sample$suitability != -1, ]
    
    s_sample$ssp = ssp[j]
    
    dat = rbind(dat, s_sample)
  
  }
  
}

output = "output/sampled-points-spam-ecocrop/"

dir.create(output, showWarnings = FALSE, recursive = TRUE)

write.csv(dat, paste0(output, "sampled-points-spam-ecocrop.csv"), 
          row.names = FALSE)
