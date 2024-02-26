# Combine variables using PCA for historical and future climate, and analyze 
# multi-dimensional shifts in climate conditions.
library("ggplot2")
library("factoextra")

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = unique(plant_spp$SPAM_Code)

dat = read.csv("output/sampled-points-spam-ecocrop/bioclim-extracted.csv", 
               na.strings = "NA")

dat = na.omit(dat)

head(dat)

i = 1

d_i = dat[dat$crop == spam[i], ]

d_i$id = paste0(d_i$x, d_i$y, d_i$ssp, d_i$suitability)

d_i = d_i[!duplicated(d_i$id), ]

table(d_i$suitability, d_i$ssp)

bio_i = d_i[,paste0("bio_", 1:19)]

bio_i = scale(bio_i)

PCA = princomp(bio_i)

plot(PCA)

g = fviz_pca_biplot(PCA, 
                    label = "all", 
                    habillage = d_i$ssp, 
                    geom="text",
                    show.legend=FALSE) 

plot(g)

g$data$suitability = factor(d_i$suitability)
g$data$ssp = factor(d_i$ssp)

g + geom_point(aes(color = suitability,
                   shape = ssp))

