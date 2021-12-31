#Useful functions for working with coati data

library(circular)
library(dismo)
library(lubridate)
library(gplots)
library(viridis)
library(raster)
library(rgdal)
library(zoo)

#LAT/LON TO UTM CONVERSIONS (AND VICE VERSA)
#Converts a matrix of lons and lats (lons first column, lats second column) to UTM
#Inputs:
#	LonsLats: [N x 2 matrix] of lons (col 1) and lats (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	EastingsCol1: whether eastings should be given in first column of output (default) or not
#Outputs:
#	EastNorths or NorthEasts: [N x 2 matrix] of Eastings and Northings - eastings are first column by default
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='34',southern_hemisphere=TRUE){
  latlons <- data.frame(X=LonsLats[,2],Y=LonsLats[,1])
  non.na.idxs <- which(!is.na(latlons$X) & !is.na(latlons$Y))
  len <- nrow(latlons)
  non.na.latlons <- latlons[non.na.idxs,]
  coordinates(non.na.latlons) <- ~Y + X
  proj4string(non.na.latlons) <- CRS('+proj=longlat +ellps=WGS84 +datum=WGS84')
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  utm <- spTransform(non.na.latlons,CRS(projection.string))
  EastNorths <- matrix(NA,nrow=len,ncol=2)
  EastNorths[non.na.idxs,] <- utm@coords
  if(!EastingsCol1){
    NorthEasts <- EastNorths[,c(2,1)]
    return(NorthEasts)
  } else{
    return(EastNorths)
  }
}

#Converts a matrix of eastings and northings (eastings first column, northings second column) to UTM
#Inputs:
#	EastNorths: [N x 2 matrix] of eastings (col 1) and northings (col2)
#	utm.zone: [numeric or string], by default 34 (this is where the KRR is)
#	southern_hemisphere: [boolean], by default TRUE
#	LonsCol1: whether lons should be given in first column of output (default) or not
#Outputs:
#	LonLats or LatLons: [N x 2 matrix] of longitudes and latitudes - lons are first column by default 
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '34',southern_hemisphere=TRUE){
  utms <- data.frame(X=EastNorths[,1],Y=EastNorths[,2])
  non.na.idxs <- which(!is.na(utms$X) & !is.na(utms$Y))
  len <- nrow(utms)
  if(southern_hemisphere){
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +south',sep='')
  } else{
    projection.string <- paste('+proj=utm +zone=',utm.zone, '+ellps=WGS84 +north',sep='')
  }
  non.na.utms <- SpatialPoints(utms[non.na.idxs,],proj4string=CRS(projection.string))
  lonlat <- spTransform(non.na.utms,CRS('+proj=longlat +ellps=WGS84 +datum=WGS84'))
  LonLats <- matrix(NA,nrow=len,ncol=2)
  LonLats[non.na.idxs,] <- lonlat@coords
  if(!LonsCol1){
    LatLons <- LonLats[,c(2,1)]
    return(LatLons)
  } else{
    return(LonLats)
  }	
}


#Generate a bunch of images (PNGs) showing trajectory data over a given time period (optionally on top of a Google Earth map)
#These can then be ffmpeg'ed together into a video.
#Inputs:
#	lats: [N x T matrix] of latitude coordinates (N individuals, T timesteps)
#	lons: [N x T matrix] of longitude coordinates (N individuals, T timesteps)
#	start.time: [numeric] time index at which to start the video
#	end.time: [numeric] time index at which to end the video
#	tail.time: [numeric] number of previous time steps to plot as a "tail" which trails the point showing the current location
#	base.dir: [string] directory in which to store the folder of outputted images
#	colors: OPTIONAL [vector of strings of length N] to specify different colors for different individuals (if NULL, default is to color individuals randomly using a rainbow palette)
#	map: OPTIONAL [gmap object]. If included, trajectories will be plotted on a Google Earth map. If NULL (default), they will be plotted on a white background.
#	ind.names: [vector of strings] giving names of individuals
#	plot.legend: [boolean] whether to show a legend
#	times: [vector of times] to show on map
#	show.all.past: [boolean] whether to show all the past data in a lighter color 
#	calls: OPTIONAL [data frame] of calls to plot on top of trajectories (defaults to NULL)
#	call.persist.time: OPTIONAL [numeric] how long the calls are shown for (by default, this matches tail time)
#	call.types: OPTIONAL [data frame] of call types and colors / symbols to plot them in (defaults to black circles for all)
#	zoom: OPTIONAL [numeric] from 8 to 21 (21 is fully zoomed in)
# places: either NULL (if no landmark places) or a data frame containing:
#   $name: name of each place
#   $lon, lat: lon and lat coordinates of each place
#   $color: color to plot in
#   $cex: size of the point
#   $pch: pch of the point
# traces: either NULL (if no 'traces' i.e. paths) or a list of n x 2 matrices containing lon/lat coordinates of lines  
#Outputs:
#	a folder of PNG images inside base.dir
trajectories.movie <-function(lats,lons,start.time,end.time,step=1,tail.time=9,base.dir='/Users/astrandb/Desktop',colors=NULL,map=NULL,ind.names=NULL,plot.legend=F,times=NULL,show.all.past=FALSE,show.scale.bar=TRUE,scale.bar.len=1000,utm.zone=36,southern_hemisphere=T,scale.bar.text=NULL,scale.bar.text.offset=0,playback.time=NULL,playback.lat=NULL,playback.lon=NULL,inds=NULL,calls=NULL,call.persist.time=NULL,call.types=NULL,zoom=20,past.step=1,places=NULL, traces = NULL){
  
  #create directory in which to store images
  dir.create(paste(base.dir,'/seq',start.time,'-',end.time,sep=''))
  
  #number of individuals
  N <- dim(lats)[1]
  
  #get colors if not specified
  if(is.null(colors)){
    colors <- rainbow(N)
  }
  
  #get lats and lons to use
  if((start.time - tail.time) > 1){
    currlats = lats[,(start.time-tail.time):end.time]
    currlons = lons[,(start.time-tail.time):end.time]
  } else{
    currlats = lats[,1:end.time]
    currlons = lons[,1:end.time]
  }
  
  
  #get times vector
  if(!is.null(times)){
    if((start.time - tail.time) > 1){
      currtimes <- times[(start.time-tail.time):end.time]
    } else{
      currtimes <- times[1:end.time]
    }
  }
  
  if((start.time - tail.time) < 1){
    min.time <- 1
  } else{
    min.time <- start.time - tail.time
    
  }
  #get map boundaries - depending on your scale you may need to adjust the .0005's
  if(!is.null(map)){
    bb <- attr(map, 'bb')
  }
  
  if(is.null(inds)){
    ind.idxs <- seq(1,N)
  } else{
    ind.idxs <- inds
  }
  
  xmin = quantile(lons[ind.idxs,min.time:end.time],0.0001,na.rm=T)
  xmax = quantile(lons[ind.idxs,min.time:end.time],0.9999,na.rm=T)
  ymin = quantile(lats[ind.idxs,min.time:end.time],0.0001,na.rm=T)
  ymax = quantile(lats[ind.idxs,min.time:end.time],0.9999,na.rm=T)
  
  #get legend info if needed
  if(plot.legend){
    legend.x <- xmin + (xmax-xmin)/100
    legend.y <- ymax - (ymax-ymin)/100
    legend.x <- bb$ll.lon + (bb$ur.lon - bb$ll.lon) / 10
    legend.y <- bb$ur.lat - (bb$ur.lat - bb$ll.lat) / 10
  }
  
  #get time colors
  time.cols.palette <- colorRampPalette(c('black','black','orange','yellow','yellow','yellow','orange','black','black'))
  time.cols <- time.cols.palette(24)
  
  #if scale bar needed, source utm to lat lon, and get easts and norths
  if(show.scale.bar){
    
    #get UTM boundaries
    eastsNorths <- latlon.to.utm(cbind(bb$ll.lon, bb$ll.lat), southern_hemisphere = southern_hemisphere, utm.zone = utm.zone)
    ll.x <- eastsNorths[,1]
    ll.y <- eastsNorths[,2]
    
    eastsNorths <- latlon.to.utm(cbind(bb$ur.lon, bb$ur.lat), southern_hemisphere = southern_hemisphere, utm.zone = utm.zone)
    ur.x <- eastsNorths[,1]
    ur.y <- eastsNorths[,2]
    
    #get scale bar values in UTM
    scalemax.x <- ur.x - (ur.x - ll.x) / 10
    scalemin.x <- scalemax.x - scale.bar.len
    scale.y <- ll.y + (ur.y - ll.y) / 10
    
    #convert to lat/lon
    lonsLats <- utm.to.latlon(cbind(scalemin.x, scale.y), southern_hemisphere = southern_hemisphere, utm.zone = utm.zone)
    scale.minlon <- lonsLats[,1]
    scale.minlat <- lonsLats[,2]
    lonsLats <- utm.to.latlon(cbind(scalemax.x, scale.y), southern_hemisphere = southern_hemisphere, utm.zone = utm.zone)
    scale.maxlon <- lonsLats[,1]
    scale.maxlat <- lonsLats[,2]
    
  }
  
  #idx<-1
  idx <- 1
  prev.lons <- lons[,start.time]
  prev.lats <- lats[,start.time]
  for(i in seq(start.time,end.time,step)){
    
    #get lats and lons for past and present data
    x<-lons[,i]
    y<-lats[,i]
    
    for(j in 1:N){
      if(!is.na(x[j])){
        prev.lons[j] <- x[j]
        prev.lats[j] <- y[j]
      }
    }
    
    #get data for tail
    if(tail.time > 0){
      if((i - tail.time) < 1){
        past.idxs <- 1:i
      } else{
        past.idxs <- (i-tail.time):i
      }
      x.past<-as.matrix(lons[,past.idxs])
      y.past<-as.matrix(lats[,past.idxs])
    }
    
    #get data for translucent tail (all past data)
    if(show.all.past){
      x.past.all <- as.matrix(lons[,seq(1,i,past.step)])
      y.past.all <- as.matrix(lats[,seq(1,i,past.step)])
    }
    
    #make figure
    filename = paste(base.dir,'/seq',start.time,'-',end.time,'/',idx,'.png',sep='')
    png(file=filename,width=8,height=8,units='in',res=300)
    par(mar=c(0,0,0,0))
    
    #initialized plot or map
    if(!is.null(map)){
      plot(map,interpolate=TRUE)
      par(usr = c(bb$ll.lon, bb$ur.lon, bb$ll.lat, bb$ur.lat))
    } else{
      par(bg='black')
      plot(NULL,xlim=c(xmin,xmax),ylim=c(ymin,ymax),xaxt='n',yaxt='n',xlab='',ylab='',bg='black',asp=1)
    }
    
    #plot all past data (if required)
    if(show.all.past){
      for(j in 1:N){
        points(x.past.all[j,],y.past.all[j,],col=adjustcolor(colors[j],alpha.f=0.1),bg=adjustcolor(colors[j],alpha.f=0.1),pch=19,cex=0.05)
      }
      
    }
    
    #plot "tails" (past locations)
    for(j in 1:N){
      lines(x.past[j,],y.past[j,],col=colors[j],lwd=2)
      #points(x.past[j,],y.past[j,],col=colors[j],cex=0.7,pch=19,bg=colors[j])
    }
    
    #add a playback (if needed)
    if(!is.null(playback.time)){
      if(abs(i-playback.time) < 1200){
        if(i < playback.time){
          points(playback.lon,playback.lat,pch=15,cex=2,col='black')
        }
        if(i >= playback.time){
          if(abs(i-playback.time)<10){
            points(playback.lon,playback.lat,pch=15,cex=2,col='yellow')
            points(playback.lon,playback.lat,pch=8,cex=2,col='black')
          } else{
            points(playback.lon,playback.lat,pch=15,cex=2,col='yellow')
          }
        }
      }		
    }
    
    #plot a legend
    if(plot.legend & idx < 60){
      if(!is.null(ind.names) & is.null(call.types)){
        legend(legend.x,legend.y,as.character(ind.names),col=colors,fill=colors,border='white',text.col='white',box.col='white',bty='n')
      }
      if(!is.null(ind.names) & !is.null(call.types)){
        legend('bottomleft',c(as.character(ind.names),'',as.character(call.types$type)),col=c(colors,'black',as.character(call.types$col)),pch=c(rep(19,length(ind.names)+1),call.types$sym),border='white',text.col='white',box.col='white',bty='n')
      }
    }
    
    #plot the time
    if(!is.null(times)){
      #text(x=xmin,y=ymin,labels=times[i],col=time.cols[hour(times[i])+1],cex=1)
      #text(x=xmin+(xmax-xmin)/8,y=ymin + (ymax - ymin)/15,labels=times[i],col='white',cex=1)
      text(x=bb$ll.lon+(bb$ur.lon-bb$ll.lon)/4,y=bb$ll.lat + (bb$ur.lat - bb$ll.lat)/15,labels=times[i],col=time.cols[hour(times[i])+1],cex=2)
    }
    
    #make a scale bar
    if(show.scale.bar){
      #lines(c(scale.minlon,scale.maxlon),c(scale.minlat,scale.maxlat),lwd=2,col='white')
      lines(c(scale.minlon,scale.maxlon),c(scale.minlat,scale.maxlat),lwd=2,col='white')
      #text(x=text.lon,y=text.lat,labels=scale.bar.text,col='white')
      text(x=(scale.minlon + scale.maxlon)/2,y=scale.minlat + (bb$ur.lat - bb$ll.lat)/50,labels=scale.bar.text,col='white')
    }
    
    #plot current locations
    #for(j in 1:N){
    #	points(prev.lons[j],prev.lats[j],col=adjustcolor(colors[j],alpha.f=0.2),bg=adjustcolor(colors[j],alpha.f=0.2),pch=19,cex=1)
    #}
    if(!is.null(map)){
      points(x,y,pch=21,cex=1.5,col='white',bg=colors)
    } else{
      points(x,y,pch=21,cex=.5,col=colors,bg=colors)
    }
    
    #plot calls in the past
    if(!is.null(calls)){
      calls.curr <- calls[which(is.na(calls$nonfoc) & calls$t.idx < i & calls$t.idx >= (i-call.persist.time)),]
      if(nrow(calls.curr)>0){
        calls.lon <- lons[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
        calls.lat <- lats[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
        calls.type.idxs <- match(calls.curr$call.type,call.types$type)
        calls.col = as.character(call.types$col[calls.type.idxs])
        calls.sym = call.types$sym[calls.type.idxs]
        points(calls.lon,calls.lat,col=calls.col,pch=calls.sym,cex=1)
      }
    }
    
    #plot current calls
    if(!is.null(calls)){
      calls.curr <- calls[which(calls$t.idx == i & is.na(calls$nonfoc)),]
      if(nrow(calls.curr)>0){
        calls.lon <- lons[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
        calls.lat <- lats[cbind(calls.curr$ind.idx,calls.curr$t.idx)]
        calls.type.idxs <- match(calls.curr$call.type,call.types$type)
        calls.col = as.character(call.types$col[calls.type.idxs])
        calls.sym = call.types$sym[calls.type.idxs]
        points(calls.lon,calls.lat,col=calls.col,pch=calls.sym,cex=2)
      }
    }
    
    if(!is.null(traces)){
      for(tr in 1:length(traces)){
        xtr <- traces[[tr]][,1]
        ytr <- traces[[tr]][,2]
        lines(xtr,ytr,lwd=2,col='white')
      }
    }
    
    if(!is.null(places)){
      points(places$lon,places$lat,col=places$col,pch=places$pch,cex=places$cex,lwd=places$lwd)
    }
    
    dev.off()
    
    idx<-idx+1
    
  }
  
}
