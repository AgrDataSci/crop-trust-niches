library("terra")
library("geodata")

files = list.files("data/SPAM", pattern = "A.tif", full.names = TRUE)

spam_crop = gsub("data/SPAM/spam2010V2r0_global_H_|_A.tif",
                "", 
                files)

spam_crop = tolower(spam_crop)

dat = data.frame()

for (i in seq_along(spam_crop)) {
  
  w = rast(files[i])
  
  w[w <= 500] = 0
  
  w2 = as.data.frame(w, xy = TRUE)
  
  names(w2)[3] = "layer" 
  
  w2 = w2[w2$layer >= 500, ]
  
  rownames(w2) = 1:nrow(w2)
  
  set.seed(9816)
  s = sample(1:nrow(w2), 2000)
  
  xy = w2[s, c(1, 2) ]
  
  xy$crop = spam_crop[i]
  
  dat = rbind(dat, xy)
  
}

plot(dat[,c(1:2)], col = factor(unique(dat$crop)))

write.csv(dat, "data/spam-crop-presence.csv", row.names = FALSE)


