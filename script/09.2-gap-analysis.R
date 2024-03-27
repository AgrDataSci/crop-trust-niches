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
library("sf")
library("patchwork")
library("RColorBrewer")


here = "output/prod-areas-gr-clusters/"

list.files(here)

pathprod = "data/future-prod-areas"

clust = read.csv("output/cluster-accessions/cluster-data-landrace-cwr.csv")

names(clust)

gr = unique(clust$GR)

spam = toupper(unique(clust$SPAM_Code))

crops = unique(clust$NAME)

crops = crops[-7]

bio = c("bio1", "bio5", "bio6", "bio12", "bio16", "bio17")

ssp = "ssp370" #c("ssp126", "ssp245", "ssp370", "ssp585")

year = c("hist", "2050")#, "2070")

adm = st_read("data/world-map/world_borders_adm0.shp")
adm = st_as_sf(adm)
adm = adm[adm$CONTINENT != "Antarctica", ]

nk = 16

mycolors = colorRampPalette(c('#e41a1c','#377eb8','#4daf4a','#984ea3',
                              '#ff7f00','#ffff33','#a65628','#f781bf'))(nk)




# .........................
# .........................
# SPAM prod areas 
spam_files = list.files("data/SPAM", pattern = "A.tif", full.names = TRUE)

result = data.frame()

for(i in seq_along(spam)) {
  
  print(spam[i])
  
  s = spam_files[grep(spam[i], spam_files)]
  
  s = rast(s)
  
  s[s <= 500] = 0
  
  s[s > 500] = 20
  
  myext = ext(s)

  for(j in seq_along(gr)) {
    
    for (k in seq_along(year)) {
      
      # read climate data
      if(isTRUE(year[i] == "hist")) {
        b = worldclim_global("bio", res = 5, "data/wc2.1-global")
        
        names(b) = paste0("bio", 1:19)
        
        b = subset(b, subset = bio)
        
      }
      
      if (isFALSE(year[i] == "hist")) {
        b = paste0(pathprod, "/", year[i], "/wc2.1_5m_av5_", 
                   ssp[1], "_", bio, "_", year[i], ".tif")
        
        b = rast(b)
        
        names(b) = bio
      }
      
      r = paste0(here, crops[i], " - ", gr[j], "-", year[k], "-ssp370.tif")
      
      r = rast(r)
      
      r = terra::crop(r, myext)
      
      r = sum(r, s)
      
      r[r < 20] = 0
      
      dat = as.data.frame(r, xy = TRUE)
      
      dat$sum[dat$sum == 0] = NA
    
      dat = na.omit(dat)
      
      dat$sum = dat$sum - 20
      
      dat2 = cbind(terra::extract(b, dat[,c("x", "y")]))
      
      dat = cbind(dat, dat2)
      
      # get climate ranges
      ranges = dat %>% 
        group_by(sum) %>% 
        summarise(Cluster = unique(sum),
                  bio1 = paste0(round(min(bio1), 0), " -- ", round(max(bio1), 0)),
                  bio5 = paste0(round(min(bio5), 0), " -- ", round(max(bio5), 0)),
                  bio6 = paste0(round(min(bio6), 0), " -- ", round(max(bio6), 0)),
                  bio12 = paste0(round(min(bio12), 0), " -- ", round(max(bio12), 0)),
                  bio16 = paste0(round(min(bio16), 0), " -- ", round(max(bio16), 0)),
                  bio17 = paste0(round(min(bio17), 0), " -- ", round(max(bio17), 0))) %>% 
        select(-sum)
      
      
      val = round(table(dat$sum) / length(dat$sum), 3)
      
      res = data.frame(Crop = crops[i],
                       GR = gr[j], 
                       Cluster = names(val),
                       Match = as.vector(val),
                       Year = year[k])
      
      res = merge(res, ranges, by = "Cluster")
      
      result = rbind(result, res)
      
    }
    
    
  } 
}


result = split(result, result$Year)

dat = result[[1]]
dat$Match2 = result[[2]]$Match

write_csv(dat, "output/gap-analysis-crops.csv")



dat = do.call("rbind", result)

write_csv(dat, "output/gap-analysis-crops-long.csv")

dat$Cluster = as.factor(as.integer(dat$Cluster))

dat$Year = ifelse(dat$Year == "hist", "Historical", dat$Year)

dat$Year = factor(dat$Year, levels = c("Historical", "2050"))

dat$GR

dat %>% 
  filter(GR == "Landrace") %>% 
ggplot(aes(y=Crop, x=Match, fill=Cluster)) +
  geom_bar(stat="identity", position="fill") +
  theme_light() +
  facet_wrap(~ Year) +
  theme_bw() +
  scale_fill_manual(values = mycolors,
                     guide = guide_legend(ncol = 2)) +
labs(title="Landrace",
       x="Proportion of Matching Areas",
       y="")

ggsave("output/gap-landrace.png",
       plot = last_plot(),
       width = 30,
       height = 15,
       units = "cm")


dat %>% 
  filter(GR == "CWR") %>% 
  ggplot(aes(y=Crop, x=Match, fill=Cluster)) +
  geom_bar(stat="identity", position="fill") +
  theme_light() +
  facet_wrap(~ Year) +
  theme_bw() +
  scale_fill_manual(values = mycolors,
                    guide = guide_legend(ncol = 2)) +
  labs(title="Crop wild relatives",
       x="Proportion of Matching Areas",
       y="")

ggsave("output/gap-cwr.png",
       plot = last_plot(),
       width = 30,
       height = 15,
       units = "cm")
