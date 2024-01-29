library("geodata")
library("Recocrop")
library("terra")
library("data.table")

# https://reagro.org/blocks/when/ecocrop.html
wcpath = "data/wc2.1-global"

res = 5

rain = worldclim_global("prec", res = res, wcpath)
names(rain) = paste0("rain", 1:length(names(rain)))

plot(rain)

tmin = worldclim_global("tmin", res = res, wcpath)
names(tmin) = paste0("tmin", 1:length(names(tmin)))

tmax = worldclim_global("tmax", res = res, wcpath)
names(tmax) = paste0("tmax", 1:length(names(tmax)))

tavg = mosaic(tmin, tmax, fun = "mean")
names(tavg) = paste0("tavg", 1:length(names(tavg)))

crop = ecocropPars("Wheat, common")

crop

m = ecocrop(crop)
m
plant = predict(m, prec = rain, tavg = tavg)
p = classify(plant > 0, cbind(0,NA)) * 1:12
pm = median(p, na.rm = TRUE)
hv = pm + m$duration
hv = ifel(hv > 12, hv - 12, hv)
plot(hv)


control(m, which_max=TRUE)
mplant = predict(m, prec = rain, tavg = tavg)
mhv = mplant + m$duration
mhv = ifel(mhv > 12, mhv - 12, mhv)
plot(mhv)

# crop parameters obtained from 
# https://github.com/AramburuMerlos/globcropdiv
# http://doi.org/10.1088/1748-9326/ac62ab
# dpars = fread("data/calib-eco-pars.csv")
# crops = unique(dpars$crop)
# 
# rpar = c("RMIN", "ROPMN", "ROPMX", "RMAX")
# phpar = c("PHMIN", "PHOPMN", "PHOPMX", "PHMAX")
# tpar = c("TMIN","TOPMN","TOPMX","TMAX")
# mrpar = paste0(rpar, "_M")
# 
# 
# 
# for(i in 1:length(crops)){
#    
# spp = which(dpars$crop == crops[i])
# score_spp = paste0("score_", spp) 
#   
#   # run for each species within crop category
#   for(j in 1:length(spp)){
#     # crop species parameters
#     crop_pars = cbind(
#       duration = c(dpars[spp[j], duration], NA, NA, NA),
#       ktmp = dpars[spp[j], KTMP] + c(-1, +1, NA, NA),
#       tavg = unlist(dpars[spp[j], ..tpar], use.names = FALSE),
#       prec = unlist(dpars[spp[j], ..mrpar], use.names = FALSE),
#       ph = unlist(dpars[spp[j], ..phpar], use.names = FALSE),
#       anpr = unlist(dpars[spp[j], ..rpar], use.names = FALSE)
#     )
#     
#     foo = list(name = dpars[spp[j], NAME], 
#                 parameters = crop_pars)
#     # rainfed 
#     m = ecocrop(foo)
#     
#     control(m, get_max = TRUE)
#     
#     plant = predict(m, prec = rain, tavg = tavg)
#     
#     p = classify(plant > 0, cbind(0,NA)) * 1:12
#     pm = median(p, na.rm = TRUE)
#     hv = pm + m$duration
#     hv = ifel(hv > 12, hv - 12, hv)
#     plot(hv)
#     
#     
#     
#     dynamicPredictors(m) = cbind(
#       ktmp = as.vector(t(d[, ..tmin])),
#       tavg = as.vector(t(d[, ..tavg])),
#       prec = as.vector(t(d[, ..prec]))
#     ) 
#     staticPredictors(m) = cbind(
#       ph = d$ph,
#       anpr = d$anpr
#     )
#     d[, rfp:= run(m)]
#     
#     # irrigated
#     m = ecocrop(foo)
#     control(m, get_max=TRUE)
#     dynamicPredictors(m) = cbind(
#       ktmp = as.vector(t(d[, ..tmin])),
#       tavg = as.vector(t(d[, ..tavg]))
#     ) 
#     staticPredictors(m) = cbind(ph = d$ph)
#     d[, irp:= run(m)]
#     
#     # species prediction
#     d[, (score_spp[j]):= rfp * (1 - irri/100) + irp * irri/100]
#     d[, c("rfp", "irp"):= NULL]
#   }
#   # get maximum score for all species within crop i
#   d[, (crops[i]):= do.call(pmax, .SD), .SDcols = score_spp]
#   d[, (score_spp):= NULL]
#   
#   # write values to spatRast and disk
#   r = rast(totcl)
#   v = rep(NA_real_, ncell(r))
#   v[d$cell] = d[[crops[i]]]
#   values(r) = v
#   writeRaster(r, filename = file.path(outdir, paste0(crops[i], ".tif")),
#               overwrite = T, 
#               wopt = list(names = crops[i], filetype = "GTiff",
#                           gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
#   )
# }
# 
