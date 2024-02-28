# Combine variables using PCA for historical and future climate, and analyze 
# multi-dimensional shifts in climate conditions.
library("ggplot2")
#library("factoextra")
#library("GDAtools")
library("patchwork")
library("ggfortify")

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
  
  d_i = dat[dat$ssp == ssp[i], ]
  
  d_i = na.omit(d_i)
  
  bio_i = d_i[,paste0("bio", c(1, 2, 4, 5, 6, 7, 12, 16, 17, 18, 19))]
  
  pca = scale(bio_i)
  
  PCA2 = princomp(pca)
  
  if(i == 8) leg = "bottom" else leg = "none"
  
  plots[[i]] = 
    autoplot(PCA2, 
           data = d_i, 
           color = "crop",
           loadings.label = TRUE, 
           loadings = TRUE, 
           loadings.colour = 'grey20',
           loadings.label.size = 5,
           loadings.label.color = "grey20") +
    theme_bw() +
    labs(title = ssp_labs[i]) +
    scale_color_brewer(palette = "Paired") +
    theme(legend.position = leg,
          legend.title = element_blank())
}

p = 
  plots[[1]] + plots[[2]] + plots[[3]] +
  plots[[4]] + plots[[5]] + plots[[6]] +
  plots[[7]] + plots[[8]] + plots[[9]] +
  plot_layout(ncol = 3)
p


ggsave("output/climate-envelope-ecocrop.pdf",
       plot = p,
       width = 35,
       height = 35,
       units = "cm",
       dpi = 600)


# 
# g = fviz_pca_biplot(PCA,
#                     label = "all",
#                     habillage = d_i$crop,
#                     geom = "text",
#                     show.legend=FALSE)
# 
# #plot(g)
# 
# g$data$suitability = factor(d_i$crop)
# g$data$ssp = factor(d_i$ssp)
# 
# g + geom_point(aes(color = suitability,
#                    shape = ssp))
# 


# #bio_i = scale(bio_i)
# 
# crops = as.factor(d_i$crop)
# 
# PCA <- bcPCA(bio_i, crops)
# # categories of class
# p1 = plot(PCA, 
#           choix = "ind", 
#           invisible = "none", 
#           label = "none",
#           col.ind.sup =  c('#9e0142'))
# p1
# 
# # variables in decathlon data
# p2 = plot(PCA, choix = "varcor")
# p2
# # between-class inertia percentage
# p1 + p2


# # bio_i
# # 
# # summary(bio_i)
# # 
