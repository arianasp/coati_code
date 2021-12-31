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
latlon.to.utm <- function(LonsLats,EastingsCol1 = TRUE,utm.zone='17',southern_hemisphere=FALSE){
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
utm.to.latlon <- function(EastNorths,LonsCol1=TRUE,utm.zone = '17',southern_hemisphere=FALSE){
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
