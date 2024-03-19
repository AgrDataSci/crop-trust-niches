# .................................................
# .................................................
# Task 2b. 
# Cluster accessions within species and group (landrace group or CWR genepool) 
# using current climate data. Consider pooling CWR together. Climate groups 
# not to be used in case landrace groups are climate-based.
#
# This analysis was performed with a range of k (number of groups) to assess
# variability at different levels of aggregation (or resolutions of the analysis).
# Since our objective is to understand the range and diversity of climate 
# characteristics mainly under historical climates (i.e., climates likely to 
# represent the range of climate adaptations of the accessions), we only performed
# this analysis using historical climate data. We calculated central tendency, 
# variability, and probability distributions for the various groups identified,
# focusing on those climate variables that most explain overall variation. 
# To facilitate interpretation, we produced summary tables, biplots, cumulative
# probability distribution plots, and force-directed graphs.

library("terra")
library("geodata")
library("tidyverse")
library("magrittr")
library("factoextra")
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

cwr$taxon = cwr$taxon_final

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

landrace$taxon = landrace$crop

dat = rbind(cwr[,c("GR", "NAME", "SPAM_Name", "SPAM_Code", "taxon", "x","y")],
            landrace[,c("GR", "NAME", "SPAM_Name", "SPAM_Code", "taxon", "x","y")])


rm(landrace, cwr)

head(dat)


# ggplot(dat, aes(x = x, y = y, color = GR)) +
#   geom_point()

# remove duplicates 
dat %<>% 
  group_by(GR, NAME, taxon, SPAM_Code) %>% 
  mutate(lonlat = paste0(x, y)) %>% 
  filter(!duplicated(lonlat)) %>% 
  ungroup() 

dat %>% 
  group_by(GR, NAME, taxon, SPAM_Code) %>% 
  summarise(n = length(x))

# ggplot(dat, aes(x = x, y = y, color = GR)) +
#   geom_point()

# ....................................
# ....................................
# Climate variables and SSP #####
wcpath = "data/wc2.1-global"

bio = worldclim_global("bio", res = 5, wcpath)
bionames = paste0("bio", 1:19)
names(bio) = bionames

plot(bio[[1]])

# model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")
# 
# keep = !duplicated(paste0(model_runs$V1, model_runs$V2, model_runs$V3))
# 
# gcm = model_runs[keep, ]
# 
# gcm = unique(gcm$V1)
# 
# keep = !duplicated(paste0(model_runs$V2, model_runs$V3))
# 
# ssp = model_runs[keep, c("V2", "V3")]

# ....................................
# ....................................
# Extract climate data #####
# read with file with models to run 
# extract bioclim current 
bio_e = terra::extract(bio, dat[, c("x", "y")])

bio_e = bio_e[, -grep("ID", names(bio_e))]

dat2 = cbind(dat, bio_e)

head(dat2)

dat2$group = paste(dat$GR, dat$SPAM_Code, sep = " - ")

table(dat2$group)

groups = unique(dat2$group)

i = 3

groups[i]

keep = dat2$group == groups[i]

d = dat2[keep, bionames]

d = scale(d)

d = dist(d)

clust = hclust(d, method = "ward.D2")

plot(clust, cex = 0.5, hang = -1)

clust_sqt = sqrt(clust$height)

plot(clust_sqt)

k = cutree(clust, 20)

table(k)

longlat = dat[keep, c("x", "y")]

longlat$clust = as.factor(k)

ggplot(longlat, aes(x = x, y = y, color = clust)) +
  geom_point()

m = as.matrix(d)

z = fviz_nbclust(m, kmeans, "gap_stat")



