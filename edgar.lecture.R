#This script presents an introduction to spatial data analysis with R
#It is based on Chapters 2 and 9 of the book:
# Applied Spatial Data Analysis with R
# Roger S. Bivand, Edzer Pebesma and V. Gómez-Rubio
# UseR! Series, Springer, 2nd ed. 2013

# Running this script requires that the code & data bundles from the book 
# be unzipped in the same folder as the script. 
# For instance the zip bundle for chapter 2 is 
# available at: 
# http://www.asdar-book.org/data2ed.php?chapter=2


#set up: clear all  & set wd
rm(list=ls()) 
setwd("C:/Dropbox/253/lat_bundle")

#load required libraries
library(rgdal) #processes GDAL - Geospatial Data Analysis files
library(spdep) #Spatial dependence: weighting schemes, statistics and models

#Read data sets

# Leukemia incidence in 281 census tracts, 8 central NY state counties (1980)
NY8 <- readOGR(".", "NY8_utm18") 
# read Shapefile (spdep package), 
# arguments (1)"." - source: a directory; here current wd 
# (2) layer: shapefile name
# Note: find shapefiles by googling "[location] shapefile"

#city names
cities <- readOGR(".", "NY8cities") 
#locations of 11 inactive hazardous waste sites; TCE: Trichloroethylene
TCE <- readOGR(".", "TCE")

# How is the data stored?
class(NY8) 
getClass("SpatialPolygonsDataFrame")
#it's a SpatialPolygonsDataFrame - a data.frame on polygons

#spatial polygons contain polygons and Spatial* characteristics
getClass("SpatialPolygons")

# a polygon is a sequence of closed lines; point coordinates where the first point equals the last 
getClass("Polygon")
#labpt - label point, centroid of polygon
# a line is an ordered list of coordinates
getClass("Line")

#finally... Spatial is the mother class of all Spatial* classes 
# used in in the sp package
getClass("Spatial")


# What does the data contain?
summary(NY8)
#Note the special slots of this class
#Coordinates
#proj4string
#Spatial data.frame - with coordinates X,Y

# What about cities
summary(cities)

#plot cities & TCE locations
plot(NY8, border="grey60", axes=TRUE)
text(coordinates(cities), labels=as.character(cities$names), font=2, cex=1.5)
plot(NY8, border="grey60", axes=TRUE)
points(TCE, pch=2, cex=1.5)
text(coordinates(TCE), labels=as.character(TCE$name), cex=1.5,
     font=1, pos=c(4,1,4,1,4,4,4,2,3,4,2), offset=0.3)

#plot one of the features - percent age > 65
spplot(NY8, c("PCTAGE65P"))#, col="transparent"

spplot(NY8, c("PCTAGE65P"), col="transparent")

#different plot: new color palette
#color palette creator function
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
#get a range for the values
tr_at <- seq(min(NY8$PCTAGE65P), max(NY8$PCTAGE65P), length.out=20)
#create a color interpolating function taking the required
#number of shades as argument
tr_rds <- rds(20)
#parameters
# at - at which values colors change
# col.regions - specify fill colors 
tr_pl <- spplot(NY8, c("PCTAGE65P"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Age>65", cex=0.8))
plot(tr_pl)

# reads a GAL lattice file into a neighbors list 
NY_nb <- read.gal("NY_nb.gal", region.id=row.names(NY8))

summary(NY_nb) #which states are neighbors?

plot(NY8, border="grey60", axes=TRUE)
plot(NY_nb, coordinates(NY8), pch=19, cex=0.6, add=TRUE)

#Data Analysis:

# linear model with:
# z = log(1000(Y_i+1)/n_i) - transformation of leukemia cases
# exposure to TCE- log inverse distance from nearest site
# percent aged >65, percent owning home
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
summary(nylm)


#plot fit and residual
NY8$lm_fit <- nylm$fit
NY8$lm_residual <- nylm$residuals
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
fit_pl <- spplot(NY8, c("lm_fit"), col="transparent", cex=0.8)
res_pl <- spplot(NY8, c("lm_residual"), col="transparent", cex=0.8)
plot(fit_pl, split=c(1,1,2,1), more=TRUE)
plot(res_pl, split=c(2,1,2,1), more=FALSE)


# Autoregressive model
# generate weight object from neighbor list object
# "B" - generates binary weights
NYlistw<-nb2listw(NY_nb, style = "B")

# fit model (I - lambda* W)(Y- X* beta) = epsilon
nysar<-spautolm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, listw=NYlistw)

summary(nysar)
#note: lambda significantly different from 0

#plot trend and stochastic component
NY8$sar_trend <- nysar$fit$signal_trend
NY8$sar_stochastic <- nysar$fit$signal_stochastic
rds <- colorRampPalette(brewer.pal(8, "RdBu"))
tr_at <- seq(-1, 1.3, length.out=21)
tr_rds <- rds(sum(tr_at >= 0)*2)[-(1:(sum(tr_at >= 0)-sum(tr_at < 0)))]
tr_pl <- spplot(NY8, c("sar_trend"), at=tr_at, col="transparent", col.regions=tr_rds, main=list(label="Trend", cex=0.8))
st_at <- seq(-0.16, 0.39, length.out=21)
st_rds <- rds(sum(st_at >= 0)*2)[-(1:(sum(st_at >= 0)-sum(st_at < 0)))]
st_pl <- spplot(NY8, c("sar_stochastic"), at=st_at, col="transparent", col.regions=st_rds, main=list(label="Stochastic", cex=0.8))
plot(tr_pl, split=c(1,1,2,1), more=TRUE)
plot(st_pl, split=c(2,1,2,1), more=FALSE)


# Ok this is great, but how do I transform my data into spatial format?
#Key data type #1: "SpatialPointsDataFrame"
getClass("SpatialPointsDataFrame")

#read data matrix
CRAN_df <- read.table("CRAN051001a.txt", header=TRUE)
CRAN_mat <- cbind(CRAN_df$long, CRAN_df$lat)
row.names(CRAN_mat) <- 1:nrow(CRAN_mat)
str(CRAN_mat)

#set CRS 
llCRS <- CRS("+proj=longlat +ellps=WGS84")
CRAN_sp <- SpatialPoints(CRAN_mat, proj4string=llCRS)
summary(CRAN_sp)

#if you don't need CRS
llCRS <- CRS(as.character(NA))

CRAN_spdf <- SpatialPointsDataFrame(CRAN_sp, CRAN_df)
summary(CRAN_spdf)

#Key data type #1: Spatial Weights - listw
#listw is an old-style (S3) class, so has no formal definition
# it is based on the Spatial neighbors ("nb") class
vignette("nb", package = "spdep") #read more
