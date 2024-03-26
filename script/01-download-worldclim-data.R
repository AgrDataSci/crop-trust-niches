# ..................................
# ..................................
# Download worldclim data 
# first run: Feb 2024
# Kauê de Sousa CGIAR
slibrary("geodata")

wcpath = "data/wc2.1-global"

models = c("ACCESS-CM2", "EC-Earth3-Veg", 
           "INM-CM5-0", "MPI-ESM1-2-HR", 
           "MRI-ESM2-0")

ssps = c("126", "245", "370", "585")

times = c("2041-2060", "2061-2080")

vars = "bio" #c("tmin", "tmax", "prec")

res = 5

model_runs = matrix(NA, nrow = 1, ncol = 4)
# now run over climate scenarios 
for (i in seq_along(models)) {
  for(j in seq_along(ssps)) {
    for (k in seq_along(times)) {
      for (l in seq_along(vars)) {
        x = cbind(models[i], ssps[j], times[k], vars[l])
        model_runs = rbind(model_runs, x)
}}}}

model_runs = na.omit(model_runs)

write.csv(model_runs, "data/worldclim-cmip6-model-runs.csv", row.names = FALSE)

# start with worldclim current
for (i in seq_along(vars)) {
  worldclim_global(vars[i], res = res, path = wcpath)
}

# now run over climate scenarios 
for (i in seq_len(nrow(model_runs))) {
  print(paste0(i, " in ", nrow(model_runs)))
  cmip6_world(model = model_runs[i, 1],
              ssp = model_runs[i, 2],
              time = model_runs[i, 3],
              var = model_runs[i, 4],
              res = res, 
              path = wcpath)
}



