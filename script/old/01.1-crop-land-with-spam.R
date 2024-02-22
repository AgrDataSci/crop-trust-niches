library("terra")

path = "data/SPAM"

files = list.files(path, pattern = "A.tif", full.names = TRUE)

files

r = rast(files)

r2 = sum(r)

plot(r2)

r2[r2 > 0] = 3

r2[r2 == 0] = 1

plot(r2)

writeRaster(r2, "data/cropland/global-crop-land.tif")


