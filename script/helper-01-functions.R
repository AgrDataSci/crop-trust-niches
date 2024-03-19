library("ggrepel")
library("ggplot2")

# plot pca
plot_pca = function(object, labels = "", scale = 3){
  
  sums = summary(object)
  
  vars = sums$sdev^2
  
  vars = round(vars/sum(vars) * 100, 1)
  
  pcd = as.data.frame(object$scores[,1:2])
  names(pcd) = c("dim1", "dim2")
  pcd$item = labels
  
  
  pcd_text = split(pcd, pcd$item)
  
  pcd_text = lapply(pcd_text, function(x){
    data.frame(item = x$item[1],
               dim1_m = mean(x$dim1),
               dim2_m = mean(x$dim2))
  })
  
  pcd_text = do.call("rbind", pcd_text)
  
  loadings = as.data.frame(object$loadings[1:length(vars), ])
  names(loadings) = paste0("dim", 1:ncol(loadings))
  loadings$traits = rownames(loadings)
  
  pcplot = 
    ggplot(pcd) +
    geom_point(aes(x = dim1, y = dim2, color = item)) +
    geom_segment(data = loadings, aes(x = 0, 
                                      y = 0, 
                                      xend = dim1 * scale,
                                      yend = dim2 * scale),
                 arrow = arrow(length = unit(0.3, "cm"), 
                               type = "open", angle = 25),
                 linewidth = 0.7, color = "grey20") +
    geom_label(data = pcd_text, aes(x = dim1_m,
                                          y = dim2_m,
                                          label = item,
                                          color = item)) +
    geom_label_repel(data = loadings,
                    aes(label = traits,
                        x = dim1 * scale,
                        y = dim2 * scale),
                    box.padding = 0.2,
                    point.padding = 0.3,
                    size = 3.5,
                    color = "grey20", 
                    arrow = arrow(length = unit(0.3, "cm"), 
                                  type = "closed",
                                  angle = 25),
                    force = 4) +
    theme_bw() +
    theme(legend.position = "none") +
    labs(x = paste0("Comp.1 (", vars[1], "%)"),
         y = paste0("Comp.2 (", vars[2], "%)"))
  
  pcplot
}
