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
library("patchwork")
library("broom")
library("sf")
library("RColorBrewer")
source("script/helper-01-functions.r")
# ....................................
# ....................................
# Input data #####
# read file with crop parameters
output = "output/cluster-accessions/"
dir.create(output, recursive = TRUE, showWarnings = FALSE)

list.files("data")
plant_spp = read.csv("data/calib-eco-pars.csv")

adm = st_read("data/world-map/world_borders_adm0.shp")
adm = st_as_sf(adm)
adm = adm[adm$CONTINENT != "Antarctica", ]

biovif = read.csv("data/vif-bioclim-selection.csv")
biovif = gsub("_", "", biovif$bio)
biovif = paste0("bio", 1:19)[paste0("bio", 1:19) %in% biovif]

sel = list.files("data/cwr-processed", pattern = ".csv")
sel = gsub(".csv", "", sel)

crop_name = plant_spp$NAME

crop_name = crop_name[crop_name %in% sel ]

plant_spp = plant_spp[plant_spp$NAME %in% sel, ]

plant_spp = plant_spp[!duplicated(plant_spp$NAME), ]

# ....................................
# ....................................
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

# ....................................
# ....................................
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

all(biovif %in% names(bio))

# ....................................
# ....................................
# Extract climate data #####
# read with file with models to run 
# extract bioclim current 
bio_e = terra::extract(bio, dat[, c("x", "y")])

bio_e = bio_e[, -grep("ID", names(bio_e))]

dat2 = cbind(dat, bio_e)

head(dat2)

dat2$group = paste(dat$NAME, dat$GR, sep = " - ")

dat2$group = gsub("African rice", "Rice", dat2$group)

table(dat2$group)

groups = sort(unique(dat2$group))

pca_list = list()
map_list = list()
desn_list = list()
boxp_list = list()

for (i in seq_along(groups)){
  
  print(groups[i])
  
  keep = dat2$group == groups[i]
  
  di = dat2[keep, ]
  
  di[di == 0] = NA
  
  di = na.omit(di)
  
  d = scale(di[bionames])
  
  d = dist(d)
  
  hc = hclust(d, method = 'ward.D2')
  # plot(hc, cex = 0.5, hang = -1)
  
  # select optimal cutoff point using lm
  false = FALSE
  nk = 15
  
  mycolors = colorRampPalette(c('#e41a1c','#377eb8','#4daf4a','#984ea3',
                                '#ff7f00','#ffff33','#a65628','#f781bf'))(nk)
  
  while(isFALSE(false)) {
    
    k = cutree(hc, k = nk) 
    di$clust = as.factor(k)
    
    dil = di[c("clust", biovif)] %>% 
      pivot_longer(!clust , names_to = "bio", values_to = "value") %>% 
      group_by(clust, bio)
    
    dil$value = ifelse(dil$value <= 0, dil$value + 100, dil$value)
    
    dil$value = log(dil$value)
    
    mod = lm(value ~ clust, data = dil)
    
    false = all(tidy(mod)[,5] < 0.001)
    
    if(isTRUE(false)) {
      false = TRUE
    }
    
    if(isFALSE(false)) {
      nk = nk - 1
    }
    
    if(isTRUE(nk < 3)){
      false = TRUE
      nk = 3
    }
    
  }
  
  cat("using", nk, "clusters for", groups[i], "\n")
  
  pc = scale(di[bionames])
  
  PCA2 = princomp(pc)
  
  pc_plot = plot_pca(PCA2, labels = di$clust, scale = 10)  +
    scale_color_manual(values = mycolors) +
    #scale_color_brewer(palette = "Set1", direction = -1) +
    labs(title = groups[i])

  pca_list[[i]] = pc_plot
  
  ggsave(gsub(" ", "", paste0(output, groups[i], "-pca-clust.pdf")),
         plot = pc_plot,
         width = 20,
         height = 20,
         units = "cm")
  
  map_clust = ggplot() +
    geom_sf(adm$geometry,
            mapping = aes(), 
            colour = "white", 
            fill = "#828282") +
    geom_jitter(data = di, aes(x = x, y = y, color = clust)) +
    theme_void() +
    scale_color_manual(values = mycolors,
                       guide = guide_legend(ncol = 2)) +
    #scale_color_brewer(palette = "Set1") +
    theme(legend.text = element_text(size = 14),
          legend.position = c(0.11, 0.45),
          plot.title = element_text(hjust = 0.5)) +
    labs(color = groups[i]) 
  
  map_list[[i]] = map_clust
  
  ggsave(gsub(" ", "", paste0(output, groups[i], "-map-clust.pdf")),
         plot = map_clust,
         width = 40,
         height = 20,
         units = "cm")
  
  # take descriptive stats by cluster
  sums = di[c("clust", biovif)] %>% 
    pivot_longer(!clust , names_to = "bio", values_to = "value") %>% 
    group_by(clust, bio) %>% 
    summarise(mean = mean(value),
              median = median(value),
              sd = sd(value),
              min = min(value),
              max = max(value)) %>% 
    ungroup()
  
  write.csv(sums, 
            gsub(" ", "", paste0(output, groups[i], "-summaries.csv")),
            row.names = FALSE)
  
  # plot points by clusters
  boxp = di[c("clust", biovif)] %>% 
    pivot_longer(!clust , names_to = "bio", values_to = "value") %>% 
    ggplot(aes(x = bio, y = value, fill = clust)) +
    geom_boxplot() +
    facet_wrap(bio ~ ., scale = "free") +
    scale_fill_manual(values = mycolors) +
    #scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "",
         y = "Value",
         fill = groups[i])
  
  boxp_list[[i]] = boxp
  
  ggsave(gsub(" ", "", paste0(output, groups[i], "-boxplot-clust.pdf")),
         plot = boxp,
         width = 25,
         height = 25,
         units = "cm")
  
  densp = di[c("clust", biovif)] %>% 
    pivot_longer(!clust , names_to = "bio", values_to = "value") %>% 
    ggplot(aes(x = value, fill = clust)) +
    geom_density(alpha = 0.2) + 
    facet_wrap(bio ~ ., scale = "free") +
    scale_fill_manual(values = mycolors) +
    #scale_fill_brewer(palette = "Set1") +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(x = "",
         y = "Density",
         fill = groups[i])
  
  desn_list[[i]] = densp
  
  ggsave(gsub(" ", "", paste0(output, groups[i], "-densityplot-clust.pdf")),
         plot = densp,
         width = 25,
         height = 25,
         units = "cm")
  
}


groups

grep("CWR", groups)

map_cwr = map_list[[1]] + map_list[[3]] + 
  map_list[[5]] + map_list[[7]] +
  map_list[[9]] + map_list[[11]] +
  map_list[[13]] + map_list[[15]] + 
  map_list[[17]] + map_list[[19]] + 
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-cwr-map.png"),
       plot = map_cwr,
       width = 50,
       height = 55,
       units = "cm",
       dpi = 500)

map_landr = map_list[[2]] + map_list[[4]] + 
  map_list[[6]] + map_list[[8]] +
  map_list[[10]] + map_list[[12]] +
  map_list[[14]] + map_list[[16]] + 
  map_list[[18]] + map_list[[20]] + 
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-landrace-map.png"),
       plot = map_landr,
       width = 50,
       height = 55,
       units = "cm",
       dpi = 500)

pca_cwr = pca_list[[1]] + pca_list[[3]] + 
  pca_list[[5]] + pca_list[[7]] +
  pca_list[[9]] + pca_list[[11]] +
  pca_list[[13]] + pca_list[[15]] + 
  pca_list[[17]] + pca_list[[19]] + 
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-cwr-pca.png"),
       plot = pca_cwr,
       width = 40,
       height = 65,
       units = "cm",
       dpi = 500)

pca_landr = pca_list[[2]] + pca_list[[4]] + 
  pca_list[[6]] + pca_list[[8]] +
  pca_list[[10]] + pca_list[[12]] +
  pca_list[[14]] + pca_list[[16]] + 
  pca_list[[18]] + pca_list[[20]] + 
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-landrace-pca.png"),
       plot = pca_landr,
       width = 40,
       height = 65,
       units = "cm",
       dpi = 500)


groups

bp = boxp_list[[2]] + boxp_list[[4]] + 
  boxp_list[[8]] + boxp_list[[16]] +  
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-landrace-boxplot.png"),
       plot = bp,
       width = 30,
       height = 30,
       units = "cm",
       dpi = 500)


bp = boxp_list[[1]] + boxp_list[[3]] + 
  boxp_list[[7]] + boxp_list[[15]] +  
  plot_layout(ncol = 2)

ggsave(paste0(output, "all-cwr-boxplot.png"),
       plot = bp,
       width = 30,
       height = 30,
       units = "cm",
       dpi = 500)



