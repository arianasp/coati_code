#Read coati data into R and store in standard format
#Save to a file in outdir called <group_id>_gps_level0.RData
#So far this only saves the high-res data (1Hz), within the specified time window (see below)

#----------PARAMETERS-------------

#may need to be changed

#directory where data is stored
indir <- '~/Dropbox/coati/rawdata/galaxy2021/gps/2021-12-29/' 

#directory where processed data is stored
outdir <- '~/Dropbox/coati/processed/galaxy2021/'

#directory where code is stored
codedir <- '~/Dropbox/coati/coati_code/'

#start hour of 1 Hz recording each day (UTC)
start_1Hz <- 11

#end time of 1 Hz recording each day (UTC)
end_1Hz <- 14

#name of the group (to be used as the name of the output file of processed data)
group_id <- 'galaxy2021'

#----------SETUP---------
setwd(codedir)
source('coati_functions.R')

#---------LOAD AND PROCESS DATA---------------

#get list of files in the directory
files <- list.files(path = indir)

#set working directory
setwd(indir)

#loop over files and read in lat/lon data
dat <- data.frame()
for(i in 1:length(files)){
  f <- files[i]
  tmp <- read.delim(file = f, header=F, sep = ',')
  gps <- tmp[,c('V14','V16','V6','V7')]
  colnames(gps) <- c('date','time','lon','lat')
  gps$id <- gsub('_gps.txt','',f)
  dat <- rbind(dat, gps)
}

dat$date <- gsub( '\\.', '-', dat$dat)

#create datetime objects
dat$datetime <- as.POSIXct(x = paste(dat$date, dat$time), tz = 'UTC', format = "%d-%m-%Y %H:%M:%S")

#remove any dates before 2000
dat <- dat[which(dat$datetime > as.POSIXct('2000-01-01', tz= 'UTC')),]

#remove any 0 lats and lons (these are not real)
dat <- dat[which(dat$lon !=0 & dat$lat !=0),]

#remove unnecessary columns 
dat <- dat[, c('id','datetime','lon','lat')]

#add east and north (UTM) columns
eastsNorths <- latlon.to.utm(cbind(dat$lon, dat$lat), utm.zone = '17', southern_hemisphere = F)
dat$east <- eastsNorths[,1]
dat$north <- eastsNorths[,2]


#get dates in the sample
dates <- sort(unique(date(dat$datetime)))

#remove anything prior to 2000 (there are some random bad fixes in there)
dates <- dates[which(dates > date('2000-01-01'))]

#timeline for high res intervals only
times_1Hz <- c()
idx_days <- c(1)
for(i in 1:length(dates)){
  times_date <- seq(from = as.POSIXct(paste0(dates[i],' ', start_1Hz,':00:00'), tz = 'UTC'),
                    to = as.POSIXct(paste0(dates[i],' ', end_1Hz-1,':59:59'), tz = 'UTC'),
                    by = 1)
  times_1Hz <- c(times_1Hz, as.character(times_date))
  idx_days <- c(idx_days, length(times_1Hz)+1)
}

#ids table - could put more info here
coati_ids <- data.frame(id = unique(dat$id), group_id = rep(group_id,length(unique(dat$id))))

#lats, lons, xs, and ys matrices
N <- nrow(coati_ids)
lats <- lons <- xs <- ys <- matrix(NA, nrow = N, ncol = length(times_1Hz))
for(i in 1:N){
  id <- coati_ids$id[i]
  dat_id <- dat[which(dat$id==id),]
  idxs <- match(as.character(dat_id$datetime),times_1Hz)
  non.nas <- which(!is.na(idxs))
  lats[i,idxs[non.nas]] <- dat_id$lat[non.nas]
  lons[i,idxs[non.nas]] <- dat_id$lon[non.nas]
  xs[i,idxs[non.nas]] <- dat_id$east[non.nas]
  ys[i,idxs[non.nas]] <- dat_id$north[non.nas]
}

save(file = paste0(outdir,group_id,'_gps_level0.RData'), list = c('times_1Hz','xs','ys','lats','lons','coati_ids','idx_days'))



