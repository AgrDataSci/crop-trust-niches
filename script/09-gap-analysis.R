# .........................................................
# .........................................................
# Gap analysis assessing the extent to which future climates
# in key production areas are represented in the international 
# collections of the 10 selected crops
# 
# To understand differences between future climates in production areas
# and the climates represented by accessions in international collections, 
# we used the cluster classifications (genetic resources groups, hereafter) 
# developed in Task 2 as the basis, and matched the production areas to
# these groups. The analysis was performed at a chosen climate cluster 
# resolution, seeking to intersect with defined genetic groups (for landraces) 
# and genepools (for CWR). The climate matching was done using the future climate
# of the current production sample (data from Task 1). The matching was done
# based on climate similarity with a threshold (e.g., the centroid of the 
# cluster ± 2 times its standard deviation for at least 80% of the climate 
# variables) to determine when a given crop production site is ‘matched’ or 
# ‘not’ to a genetic resources group. All production sites not matched to a
# genetic resources group were analysed and plotted to understand their climate 
# characteristics and geographic distribution.
#
# Assign crop growing area sample to accession climate pools, using different 
# num of climate pools (i.e., different k cluster numbers). The matching will
# be done based on climate similarity (e.g., Mahalanobis, Gower, or Euclidean distance)
# with a threshold (e.g., the centroid of the cluster ± 2 times its standard deviation
# for at least 80% of the climate variables) to determine when a given crop production
# site is ‘matched’ or ‘not’ to a genetic resources group. Do this for (1) future 
# climate of the current production sample; and (2) future climate of the future
# production sample.
# 
# Output 1: table of accession clusters that more frequently had matching 
# crop (current, future) production area samples, produce table of means, 
# sigma, max/min, PDF plot, maps
# 
# Output 2: analyze current, future production samples that fell outside of
# accession climate pools, plot mean climate conditions, sigma, probability 
# distribution of the samples). Plot maps of where these are.
# first run: Mar 2024
# Kauê de Sousa CGIAR
# .........................................................
# .........................................................
library("terra")
library("tidyverse")
library("magrittr")

clust = read.csv("output/cluster-accessions/cluster-data-landrace-cwr.csv")

names(clust)

bio = c("bio1", "bio5", "bio6", "bio12", "bio16", "bio17")

ssp = c("ssp126", "ssp245", "ssp370", "ssp585")

year = c(2050, 2070)

pathprod = "data/future-prod-areas"

gap = clust %>% 
  group_by(group) %>% 
  summarise(nclust = length(unique(clust)),
            bio1 = paste0(round(min(bio1), 0), " - ", round(max(bio1), 0)),
            bio5 = paste0(round(min(bio5), 0), " - ", round(max(bio5), 0)),
            bio6 = paste0(round(min(bio6), 0), " - ", round(max(bio6), 0)),
            bio12 = paste0(round(min(bio12), 0), " - ", round(max(bio12), 0)),
            bio16 = paste0(round(min(bio16), 0), " - ", round(max(bio16), 0)),
            bio17 = paste0(round(min(bio17), 0), " - ", round(max(bio17), 0))) %>% 
  separate(group, c("Crop", "GR"), sep = " - ")

write_csv(gap, "output/climate-ranges-gr-crops.csv")

# plot points by clusters
pdat = 
  clust[c("group", bio)] %>% 
  pivot_longer(!group , names_to = "bio", values_to = "value") %>% 
  separate(group, c("Crop", "GR"), sep = " - ")

pdat$bio = factor(pdat$bio, levels = bio)

pdat$Crop = factor(pdat$Crop, levels = rev(unique(pdat$Crop)))

pdat %>% 
  group_by(bio) %>% 
  summarise(min = min(value),
            max = max(value))

pdat$value = ifelse(pdat$bio == "bio17" |
                      pdat$bio == "bio12" |
                      pdat$bio == "bio16", log(pdat$value),
                    pdat$value)


pdat$bio = factor(pdat$bio, labels = c("bio1", "bio5", "bio6",
                            "log(bio12)", "log(bio16)", "log(bio17)"))

pdat %>% 
  group_by(bio) %>% 
  summarise(min = min(value),
            max = max(value))

p = ggplot(pdat, aes(y = Crop, x = value, fill = GR)) +
  geom_violin() +
  facet_wrap( ~ bio, scale = "free") +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "",
       y = "")

ggsave("output/overlay-crop-gr-groups-bio-clim.png",
       plot = p,
       width = 20,
       height = 20,
       units = "cm")

# ...................................
# ...................................
ranges = 
  clust %>% 
  group_by(group, clust) %>% 
  summarise(bio1_1 = min(bio1), 
            bio1_2 = max(bio1),
            bio5_1 = min(bio5),
            bio5_2 = max(bio5),
            bio6_1 = min(bio6),
            bio6_2 = max(bio6),
            bio12_1 = min(bio12),
            bio12_2 = max(bio12),
            bio16_1 = min(bio16),
            bio16_2 = max(bio16),
            bio17_1 = min(bio17), 
            bio17_2 = max(bio17)) %>% 
  separate(group, c("Crop", "GR"), sep = " - ", remove = FALSE)

groups = unique(ranges$group)

# .........................
# .........................
# SPAM prod areas 
files = list.files("data/SPAM", pattern = "A.tif", full.names = TRUE)

spam_crop = gsub("data/SPAM/spam2010V2r0_global_H_|_A.tif",
                 "", 
                 files)

spam = rast(files)

spam[spam <= 50] = 0

spam = sum(spam)

spam[spam > 50] = 1

myext = ext(spam)

# .........................
# .........................
# Reclassify over prod areas 

i = 1
j = 2
g = 1
b = 1


for(g in seq_along(groups)) {
  
  prod_group = c()
  
  for (i in seq_along(years)) {
    
    for (j in seq_along(ssp)) {
    
      f = paste0(pathprod, "/", year[i], "/wc2.1_5m_av5_", 
                 ssp[j], "_", bio, "_", year[i], ".tif")
      
      r = rast(f)
      
      names(r) = bio
      
      r = terra::crop(r, myext)
      
      r = mask(r, spam)
      
      # do it for each bio within the crop group
      for(b in seq_along(bio)) {
        
        z = r
        
        # and now do it for each cluster
        class_matrix = ranges[ranges$group == groups[g], ]
        
        class_matrix = class_matrix[, c(paste0(bio[b], c("_1", "_2")), "clust")]
        
        class_matrix = as.matrix(class_matrix)
        
        clusters = unique(class_matrix[,3])
        
        for(k in seq_along(clusters)) {
          
          z[[b]] = classify(z[[b]], class_matrix[k, ])
          
          z[[b]][z[[b]] < 10000] = 0
        }
        
      }
        
    }
    
  }

}








