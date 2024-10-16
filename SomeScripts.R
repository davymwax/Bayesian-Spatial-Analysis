####################
## Load Packages ###
####################

library("MASS")
library("maps")
library("fields")
library("sf")
library("akima")
library("spBayes")
library("RColorBrewer")
library("classInt")
library("gstat")
library("sp")
library("tidyverse")
library("tidyr")
library("ggplot2")
library("geoR")
library("kableExtra")
library("nlme")
library("lme4")
library("CARBayes")
library("maptools")
library("spatstat")
library("SpatialEpi")
library("spatialreg")
library("spdep")
library("foreign")
library("splancs")

############################

x.grid <- seq(0, 10, length=30)
y.grid <- seq(0, 10, length=30)
n.grid <- length(x.grid)
coord.pts <- matrix(0, nrow = n.grid*n.grid,
                    ncol = 2)
coord.pts[,1] <- rep(x.grid, each=n.grid)
coord.pts[,2] <- rep(y.grid, n.grid)
locations <- plot(coord.pts[,1], coord.pts[,2], 
                  xlab="X", ylab="Y",
                  main="Locations where the spatial process \n is simulated at",
                  type="p",
                  pch=21,
                  col="gray")

