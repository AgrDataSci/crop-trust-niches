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
library("geodata")

output = "output/prod-areas-gr-clusters/"

dir.create(output, showWarnings = FALSE)

clust = read.csv("output/cluster-accessions/cluster-data-landrace-cwr.csv")

names(clust)

bio = c("bio1", "bio5", "bio6", "bio12", "bio16", "bio17")

ssp = "ssp370" #c("ssp126", "ssp245", "ssp370", "ssp585")

year = c("hist", "2050", "2070")

# this is a path to averaged bioclim variables for future climates
# I was in a hurry when I wrote this script,
# but reproducibility can be achieved with 
# averaging the layers from 
# geodata::cmip6_world(model = gcm[j, 1],
#                      ssp = gcm[j, 2],
#                      time = gcm[j, 3],
#                      var = "bio",
#                      res = 5, 
#                      path = wcpath)

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
year = 2050

for (i in seq_along(year)) {
  
  print(year[i])
  
  #  for (j in seq_along(ssp)) {
  
  if(isTRUE(year[i] == "hist")) {
    r = worldclim_global("bio", res = 5, "data/wc2.1-global")
    
    names(r) = paste0("bio", 1:19)
    
    r = subset(r, subset = bio)

  }
  
  if (isFALSE(year[i] == "hist")) {
    f = paste0(pathprod, "/", year[i], "/wc2.1_5m_av5_", 
               ssp[1], "_", bio, "_", year[i], ".tif")
    
    r = rast(f)
    
    names(r) = bio
  }
  
  r = terra::crop(r, myext)
  
  r = mask(r, spam)
  
  for(g in seq_along(groups)) {
    
    print(groups[g])
    
    clusters = unique(ranges$clust[ranges$group == groups[g]])
    
    prod_clust = c()
    
    for(k in seq_along(clusters)) {
      
      clust_r = r
      
      for(b in seq_along(bio)) {
        
        # and now do it for each cluster
        class_matrix = ranges[ranges$group == groups[g], ]
        
        class_matrix = class_matrix[, c(paste0(bio[b], c("_1", "_2")), "clust")]
        
        class_matrix = as.matrix(class_matrix)
        
        clusters = unique(class_matrix[,3])
        
        clust_r[[b]][clust_r[[b]] >= class_matrix[k, 1] & clust_r[[b]] <= class_matrix[k, 2]] = -1000
        
        clust_r[[b]][clust_r[[b]] != -1000] = 0
        
        clust_r[[b]][clust_r[[b]] == -1000] = 1
        
        
      }
      
      clust_r = sum(clust_r)
      
      clust_r[clust_r < 6] = 0
      
      clust_r[clust_r == 6] = class_matrix[k, 3]
      
      prod_clust = c(prod_clust, clust_r)
      
    }
    
    prod_clust = rast(prod_clust)
    
    plot(prod_clust)
    
    prod_clust = sum(prod_clust)  
    
    prod_clust[prod_clust > max(clusters)] = 0
    
    writeRaster(prod_clust, 
                filename = paste0(output, 
                                  groups[g], "-", 
                                  year[i], "-",
                                  ssp[1],
                                  ".tif"),
                overwrite = TRUE)
    
  }
  
}








