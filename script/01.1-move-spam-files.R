library("geodata")

from = "/Users/kauedesousa/local-workflow/geo-raster-shapefile/spam2010v2r0_global_harv_area/"

here = "data/SPAM/"

crops = read.csv("data/mapspam-crop-names.csv")

crops = toupper(crops$spamname)

identifier = paste(crops, collapse = "|")

files = list.files(from, pattern = identifier, full.names = TRUE)

f = grep("H.tif$|A.tif$", files)

f = files[f]

n = gsub(from, "", f)

file.copy(f, paste0(here, n))



