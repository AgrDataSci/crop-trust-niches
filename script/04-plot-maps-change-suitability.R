# ..................................
# ..................................
# Prepare ecocrop maps
# first run: Feb 2024
# Kauê de Sousa CGIAR
# ..................................
# ..................................
library("terra")
library("geodata")
library("sf")
library("tidyverse")
library("patchwork")

output = "output/change-suitability/"

dir.create(output, showWarnings = FALSE, recursive = TRUE)

wcpath = "data/world-map"

layerpath = "output/ecocrop/"

w = world(resolution = 5, level = 0, path = wcpath)

w = st_as_sf(w)     

w = w[w$NAME_0 != "Antarctica", ]

plant_spp = read.csv("data/calib-eco-pars.csv")

spam = plant_spp[,c("SPAM_Name", "SPAM_Code")]

spam_code = unique(spam$SPAM_Code)

spam_name = ClimMobTools:::.title_case(unique(spam$SPAM_Name))

# read with file with models to run 
model_runs = read.csv("data/worldclim-cmip6-model-runs.csv")

keep = !duplicated(paste0(model_runs$V2, model_runs$V3))

model_runs = model_runs[keep, c("V2", "V3")]

gcm = paste(model_runs$V2, model_runs$V3, sep = "-")

for (i in seq_along(spam_code)) {
  
  current = rast(paste0(layerpath, spam_code[i], "-SSP-current.tif"))
  current[current > 0] = 1
  current[current <= 0] = 0
  
  suit_plots = list()
  
  for(j in seq_along(gcm)) {
    
    print(paste0(spam_name[i], "-SSP-", gcm[j]))
    
    future = rast(paste0(layerpath, spam_code[i],
                         "-SSP-", gcm[j], ".tif"))
    future[future > 0] = 2
    future[future <= 0] = 0
    
    #identify change in suitability in RCP scenario
    #future minus current raster
    #change in suitability codes
    #  1 = always suitable
    #  0 = never suitable
    # -1 = no longer suitable
    #  2 = novel areas
    overlay = function(x, y) {(x - y)}
    change_suit = do.call(overlay, list(x = future, y = current))
    change_suit = crop(change_suit, ext(w))
    
    writeRaster(change_suit, 
                filename = paste0(output,spam_code[i], "-SSP-", gcm[j], ".tif"),
                overwrite = TRUE)
    
    r <- as.data.frame(change_suit, xy = TRUE)
    r <- r[!is.na(r[, 3]), ]
    names(r)[3] = "layer"
    r$layer = as.factor(r$layer)
    
    p <- ggplot() +
      geom_tile(r, mapping = aes(x = x, y = y, fill = layer)) +
      geom_sf(w$geometry, mapping = aes(), colour = "black", fill = NA) +
      scale_fill_manual(values = c("#d7191c", "grey95","#abd9e9","#053061")) +
      theme_void() + 
      labs(title = paste0(spam_name[i], " SSP-", gcm[j])) +
      theme(legend.position = "bottom",
            legend.title=element_blank(),
            plot.title = element_text(hjust = 0.5))
    
    suit_plots[[j]] = p
    
  }
  
  final_plot = suit_plots[[1]] + suit_plots[[2]] +
    suit_plots[[3]] + suit_plots[[4]] +
    suit_plots[[5]] + suit_plots[[6]] +
    suit_plots[[7]] + suit_plots[[8]] +  plot_layout(ncol = 2)

  ggsave(paste0(output, spam_code[i], ".png"),
         plot = final_plot,
         width = 30,
         height = 60,
         units = "cm")
  
}

