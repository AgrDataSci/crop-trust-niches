# .................................................
# .................................................
# Task 1b. 
# Combine variables using PCA for historical and future climate, and analyze 
# multi-dimensional shifts in climate conditions.
# first run: Feb 2024
# Kauê de Sousa CGIAR
# ..................................
# ..................................
library("ggplot2")
library("patchwork")

source("script/helper-01-functions.r")

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = unique(plant_spp$SPAM_Code)

dat = read.csv("output/sampled-points-spam-ecocrop/bioclim-extracted.csv", 
               na.strings = "NA")

dat = na.omit(dat)

dat$crop = factor(dat$crop, levels = sort(spam))

ssp = unique(dat$ssp)

ssp_labs = toupper(ssp)
ssp_labs = gsub("-2041-2060", " 2050", ssp_labs)
ssp_labs = gsub("-2061-2080", " 2070", ssp_labs)
ssp_labs = paste0("SSP-", ssp_labs)
ssp_labs[1] = "Historical" 

plots = list()

for(i in seq_along(ssp)) {
  
  d_i = dat
  
  d_i = d_i[d_i$ssp == ssp[i], ]
  
  #d_i = d_i[d_i$crop == "whea", ]
  
  d_i = na.omit(d_i)
  
  bio_i = d_i[,paste0("bio", c(1, 2, 4, 5, 6, 7, 12, 16, 17, 18, 19))]
  
  pca = scale(bio_i)
  
  PCA2 = princomp(pca)
  
  plots[[i]] = plot_pca(PCA2, labels = d_i$crop, scale = 10) + 
    scale_color_brewer(palette = "RdYlBu", direction = -1) +
    labs(title = ssp_labs[i])
  
}

p = 
  plots[[1]] + plots[[2]] + plots[[3]] +
  plots[[4]] + plots[[5]] + plots[[6]] +
  plots[[7]] + plots[[8]] + plots[[9]] +
  plot_layout(ncol = 3)

ggsave("output/climate-envelope-ecocrop.pdf",
       plot = p,
       width = 35,
       height = 35,
       units = "cm",
       dpi = 600)

ggsave("output/climate-envelope-ecocrop.png",
       plot = p,
       width = 45,
       height = 45,
       units = "cm",
       dpi = 600)

# plots = list()
# 
# sel = c('current', '245-2041-2060', '245-2061-2080', '370-2041-2060', '370-2061-2080')
# 
# ssp = ssp[which(ssp %in% sel)]
# 
# ssp_labs = ssp_labs[which(ssp %in% sel)]

# #i=2
# j=1
# 
# for(j in seq_along(spam)) {
#   
#   #p = list()
#   
#   #for(i in seq_along(ssp)) {
#     
#     #d_i = dat[dat$ssp == ssp[i], ]
#     d_i = dat
#     d_i = d_i[d_i$crop == spam[j], ]
#     
#     d_i = na.omit(d_i)
#     
#     bio_i = d_i[,paste0("bio", c(1, 2, 4, 5, 6, 7, 12, 16, 17, 18, 19))]
#     
#     pca = scale(bio_i)
#     
#     PCA2 = princomp(pca)
#     
#     plots[[j]] = plot_pca(PCA2, labels = d_i$ssp, scale = 3) + 
#       scale_color_brewer(palette = "RdYlBu") +
#       labs(title = spam[j])
#     
#   #}
#   
#   #plots = c(plots, p)
#   
# }
# 
# 
# 
# p = 
#   plots[[1]] + plots[[2]] + plots[[3]] +
#   plots[[4]] + plots[[5]] + plots[[6]] +
#   plots[[7]] + plots[[8]] + plots[[9]] + plots[[10]] +
#   plot_layout(ncol = 2)
# p
# 
# 
# ggsave("output/climate-envelope-ecocrop-by-crop.pdf",
#        plot = p,
#        width = 35,
#        height = 50,
#        units = "cm",
#        dpi = 600)


