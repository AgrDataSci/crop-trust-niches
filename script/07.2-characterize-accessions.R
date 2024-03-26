# ..................................
# ..................................
# Per crop and group (landrace group, or CWR genepool) produce a 
# table with the average, sigma, CV, q1, q2, q3, q4, min/max, of the 
# historical conditions of the accessions.
# first run: Mar 2024
# Kauê de Sousa CGIAR
# ..................................
# ..................................
library("ggplot2")
library("patchwork")
library("ggfortify")
library("tidyverse")

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = unique(plant_spp$SPAM_Code)

dat = read.csv("output/sampled-points-spam-ecocrop/landrace-and-cwr-bioclim-extracted.csv", 
               na.strings = "NA")

dat = na.omit(dat)

dat$crop = factor(dat$SPAM_Code, levels = sort(spam))

dat$GR = as.factor(dat$GR)

ssp = unique(dat$ssp)

ssp_labs = toupper(ssp)
ssp_labs = gsub("-2041-2060", " 2050", ssp_labs)
ssp_labs = gsub("-2061-2080", " 2070", ssp_labs)
ssp_labs = paste0("SSP-", ssp_labs)
ssp_labs[1] = "Historical" 

gr = as.character(unique(dat$GR))

crops = sort(unique(dat$NAME))

dat$NAME = factor(dat$NAME, levels = crops)

# ................................
# ................................
# Climate envelopes ####

for (g in seq_along(gr)) {
  
  plots = list()
  
  for(i in seq_along(ssp)) {
    
    d_i = dat[dat$ssp == ssp[i] & dat$GR == gr[g], ]
    
    d_i = na.omit(d_i)
    
    bio_i = d_i[, paste0("bio", 1:19)]#c(1, 2, 4, 5, 6, 7, 12, 16, 17, 18))]
    
    PCA2 = princomp(bio_i)
    
    if(i == 8) leg = "bottom" else leg = "none"
    
    plots[[i]] = 
      autoplot(PCA2, 
               data = d_i, 
               color = "NAME",
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
  
  ggsave(paste0("output/climate-envelope-", gr[g] ,"-species-name.pdf"),
         plot = p,
         width = 35,
         height = 35,
         units = "cm",
         dpi = 600)
  
  ggsave(paste0("output/climate-envelope-", gr[g] ,"-species-name.png"),
         plot = p,
         width = 35,
         height = 35,
         units = "cm",
         dpi = 600)
}



# ................................
# ................................
# Create tables with summarized information ####
group = paste(dat$GR, dat$NAME, dat$ssp, sep = "_")

dat_sum = split(dat, group)

bio_index = paste0("bio", c(1, 2, 4, 5, 6, 7, 12, 16, 17, 18, 19))

bio_sum = lapply(dat_sum, function(Y) {
  
  s = lapply(Y[bio_index], function(x){
    q = round(quantile(x, probs = c(0.25, 0.5, 0.75)), 3)
    
    data.frame(mean = round(mean(x, na.rm = TRUE), 3),
               min = round(min(x, na.rm = TRUE), 3),
               max = round(max(x, na.rm = TRUE), 3),
               sd = round(sd(x, na.rm = TRUE), 3),
               q1 = q[1],
               q2 = q[2],
               q3 = q[3])
  })
  
  s = do.call("rbind", s)
  
  s$bio = rownames(s)
  
  s
  
})

bio_sum = do.call("rbind", bio_sum)

bio_sum

values = strsplit(rownames(bio_sum), "_")

values = do.call(rbind, values)

values[,3] = unlist(lapply(strsplit(values[,3], "[.]"), function(x) x[1]))

values = as.data.frame(values)

names(values) = c("GR", "Crop", "SSP")

bio_sum = cbind(values, bio_sum)

bio_sum$bio = factor(bio_sum$bio, levels = bio_index)

write.csv(bio_sum, "output/summary-table-bio-cwr-landrace.csv", row.names = FALSE)


x = bio_sum[bio_sum$SSP == 'current' & bio_sum$GR == "CWR", ]



