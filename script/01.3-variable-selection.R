#...................................................
# Select bioclimatic variables for modelling
# using VIF analysis
# first run: Feb 2024
# Kauê de Sousa CGIAR
#...................................................
#...................................................
# Packages ####
library("data.table")
library("raster")
library("dismo")
library("BiodiversityR")
library("car")
library("geodata")
#...................................................
#...................................................
# Data ####
wcpath = "data/wc2.1-global"

# bioclimatic variables
bio = worldclim_global("bio", res = 5, wcpath)
bio = stack(bio)
names(bio) = gsub("wc2.1_5m_", "", names(bio))
# define projection and extension
myproj = proj4string(bio)
myext  = extent(bio)
myres  = res(bio)

# species acronyms
list.files("data")

# passport data
df = fread("data/spam-crop-presence.csv")

# .......................................
# .......................................
# Set background points ####
xy = df[,c(1:2)]
names(xy) = c("lon", "lat")

xy = unique(xy, by = c("lon", "lat"))

set.seed(123)
bg = randomPoints(bio[[1]], 
                  n = 10000, 
                  ext = myext, 
                  extf = 1.25)

plot(bio[[1]])
points(bg)
#...................................................
#...................................................
# Variable selection with VIF ####
vif = ensemble.VIF(
  x = bio,
  a = xy,
  an = bg,
  VIF.max = 10,
  keep = NULL,
  factors = NULL,
  dummy.vars = NULL
)

result = data.frame(bio = names(vif$VIF.final),
                    vif = vif$VIF.final)

write.csv(result, "data/vif-bioclim-selection.csv", row.names = FALSE)

