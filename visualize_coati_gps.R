#Visualize coati GPS data

#-------PARAMETERS---------

#directory where code is stored
codedir <- '~/Dropbox/coati/coati_code/'

#path to the gps data file to use
datafile <- '~/Dropbox/coati/processed/galaxy2021/galaxy2021_gps_level0.RData'

#map file
mapfile <- '~/Dropbox/coati/visualizations/coati_maps.RData'

plotdir <- '~/Dropbox/coati/visualizations'

#-----------SETUP-----------
setwd(codedir)
source('coati_functions.R')

#---------MAIN FUNCTION-------

#load data file
load(datafile)

#load maps
load(mapfile)


for(i in 2:length(idx_days)){
  t0 <- idx_days[i-1]
  tf <- idx_days[i]
  trajectories.movie(lats = lats, 
                   lons = lons, 
                   start.time= t0, 
                   end.time = tf, 
                   tail.time = 600,
                   ind.names = coati_ids$id,
                   plot.legend = T,
                   times = times_1Hz,
                   base.dir = plotdir,
                   map = coati_map17,
                   utm.zone = '17',
                   southern_hemisphere = FALSE,
                   step = 10,
                   show.scale.bar = T,
                   scale.bar.len = 200)
}
