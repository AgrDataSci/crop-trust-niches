library("terra")


plant_spp = read.csv("data/calib-eco-pars.csv")

spam = toupper(unique(plant_spp$SPAM_Code))

spam

i=1

list.files("data/spam")

s_area = list.files("data/spam/", pattern = paste0(spam[i], "_A.tif"),
                    full.names = TRUE)
s_area = rast(s_area)

#s_area = subst(s_area, from = 0, to = NA)

plot(s_area)

s_harv = list.files("data/spam/", pattern = paste0(spam[i], "_H.tif"),
                    full.names = TRUE)
s_harv = rast(s_harv)

plot(s_harv)

s_harv
s_area

s_sample = spatSample(s_area, 
                      size = 10000, 
                      method = "random", 
                      as.points = TRUE,
                      na.rm = FALSE)

plot(s_sample)



