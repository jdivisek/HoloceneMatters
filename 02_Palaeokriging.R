#####################################################################################
#		     INTERPOLATE POLLEN COMPOSITION GRADINETS                         #
#####################################################################################

library(spdep)
library(raster)
library(gstat)
library(rgdal)
library(dismo)
library(Metrics)
library(vegan)
library(fields)

#source required functions
myfunctions <- list.files(paste(getwd(), "/MyFunctions", sep=""), full.names=TRUE)
for (i in 1:length(myfunctions))
{
  source(myfunctions[i])
}

###READ DATA-------------------------------------------------------------------------------
Alt <- raster(paste(getwd(), "/rasters/Alt.tif", sep=""))
country <- readOGR(paste(getwd(), "/shapes/country_line_WGS84web.shp", sep=""))

###PERIOD1 - LATE GLACIAL--------------------------------------------------------------
period1 <- read.delim("period1.txt", header=T, row.names=1)

period1$Alt <- extract(Alt, period1[, c("x", "y")])

dem.df <- as.data.frame(Alt, xy=TRUE)

temp1.dem <- paleokriging(dat=period1, coord=period1[, c("x", "y")], y = "PCO1b", 
                          x = "Alt", new.coord = dem.df[, c("x", "y")], x.new = dem.df$Alt, 
                          radius=NULL, k=12, do.R2=T)
temp1.dem <- cbind(dem.df[, c("x", "y")], temp1.dem)
head(temp1.dem)

temp1 <- rasterFromXYZ(temp1.dem, res=c(1200, 1200), crs=Alt@crs, digits=5)
plot(temp1, axes=F, main= "Temperateness in Late Glacial")
plot(country, add=T)
points(period1[, c("x", "y")], cex=(decostand(period1$PCO1b, "range")*2)+1)

writeRaster(temp1, filename = paste(getwd(), "/history_rasters/temp1.tif", sep=""), format="GTiff", overwrite=TRUE)
temp1 <- raster(paste(getwd(), "/history_rasters/temp1.tif", sep=""))

###PERIOD2 - HOLOCENE ONSET --------------------------------------------------------------
period2 <- read.delim("period2.txt", header=T, row.names=1)

period2$Alt <- extract(Alt, period2[, c("x", "y")])

dem.df <- as.data.frame(Alt, xy=TRUE)

temp2.dem <- paleokriging(dat=period2, coord=period2[, c("x", "y")], y = "PCO_1", 
                          x = "Alt", new.coord = dem.df[, c("x", "y")], x.new = dem.df$Alt, 
                          radius=NULL, k=12, do.R2=T)
temp2.dem <- cbind(dem.df[, c("x", "y")], temp2.dem)
head(temp2.dem)

temp2 <- rasterFromXYZ(temp2.dem, res=c(1200, 1200), crs=Alt@crs, digits=5)
plot(temp2, axes=F, main= "PCO 1 (temperatnost) - Holocene onset")
plot(country, add=T)
points(period2[, c("x", "y")], cex=(decostand(period2$PCO_1, "range")*2)+1)

writeRaster(temp2, filename = paste(getwd(), "/history_rasters/temp2.tif", sep=""), format="GTiff", overwrite=TRUE)
temp2 <- raster(paste(getwd(), "/history_rasters/temp2.tif", sep=""))

###PERIOD3 - PRE-NEOLITHIC--------------------------------------------------------------
period3 <- read.delim("period3.txt", header=T, row.names=1)

period3$Alt <- extract(Alt, period3[, c("x", "y")])

dem.df <- as.data.frame(Alt, xy=TRUE)

step3.dem <- paleokriging(dat=period3, coord=period3[, c("x", "y")], y = "PCO1b", 
                          x = "Alt", new.coord = dem.df[, c("x", "y")], x.new = dem.df$Alt, 
                          radius=NULL, k=12, do.R2=T)
step3.dem <- cbind(dem.df[, c("x", "y")], step3.dem)
head(step3.dem)

step3 <- rasterFromXYZ(step3.dem, res=c(1200, 1200), crs=Alt@crs, digits=5)
plot(step3, axes=F, main= "PCO 1 (stepnost) - Pre-Neolithic")
plot(country, add=T)
points(period3[, c("x", "y")], cex= (decostand(period3$PCO1b, "range")*2)+1)

writeRaster(step3, filename = paste(getwd(), "/history_rasters/step3.tif", sep=""), format="GTiff", overwrite=TRUE)
step3 <- raster(paste(getwd(), "/history_rasters/step3.tif", sep=""))

###PERIOD4 - LATE GLACIAL--------------------------------------------------------------
period4 <- read.delim("period4_WGS84web.txt", header=T, row.names=1)

period4$Alt <- extract(Alt, period4[, c("x", "y")])

dem.df <- as.data.frame(Alt, xy=TRUE)

step4.dem <- paleokriging(dat=period4, coord=period4[, c("x", "y")], y = "PCO1b", 
                          x = "Alt", new.coord = dem.df[, c("x", "y")], x.new = dem.df$Alt, 
                          radius=NULL, k=12, do.R2=T)
step4.dem <- cbind(dem.df[, c("x", "y")], step4.dem)
head(step4.dem)

step4 <- rasterFromXYZ(step4.dem, res=c(1200, 1200), crs=Alt@crs, digits=5)
plot(step4, axes=F, main= "PCO 1 (stepnost) - Late Neolithic")
plot(country, add=T, lwd=1.5)
points(period4[, c("x", "y")], cex= (decostand(period4$PCO1b, "range")*2)+1)

writeRaster(step4, filename = paste(getwd(), "/history_rasters/step4.tif", sep=""), format="GTiff", overwrite=TRUE)
step4 <- raster(paste(getwd(), "/history_rasters/step4.tif", sep=""))

###PERIOD5--------------------------------------------------------------
period5 <- read.delim("period5_WGS84web.txt", header=T, row.names=1)

period5$Alt <- extract(Alt, period5[, c("x", "y")])

dem.df <- as.data.frame(Alt, xy=TRUE)

taiga5.dem <- paleokriging(dat=period5, coord=period5[, c("x", "y")], y = "PCO_1", 
                           x = "Alt", new.coord = dem.df[, c("x", "y")], x.new = dem.df$Alt, 
                           radius=NULL, k=12, do.R2=T)
taiga5.dem <- cbind(dem.df[, c("x", "y")], taiga5.dem)
head(taiga5.dem)

taiga5 <- rasterFromXYZ(taiga5.dem, res=c(1200, 1200), crs=Alt@crs, digits=5)
plot(taiga5, axes=F, main= "PCO 1 (taigovitost) - Late Prehistory")
plot(country, add=T)
points(period5[, c("x", "y")], cex=(decostand(period5$PCO_1, "range")*2)+1)

writeRaster(taiga5, filename = paste(getwd(), "/history_rasters/taiga5.tif", sep=""), format="GTiff", overwrite=TRUE)
taiga5 <- raster(paste(getwd(), "/history_rasters/taiga5.tif", sep=""))

###MAPS-------------------------------------------------------------------------
windows(8.44, 8.07)#812, 777
# tiff("Paleomaps.tif", 8.44, 8.07, units = "in", res=400, compression = "lzw")
layout(matrix(1:6, ncol=2, nrow=3, byrow = F))

par(mar=c(1, 0.5, 3, 5), oma=c(1.5, 3, 0, 0))

#period1
my.ramp2 <- colorRampPalette(c(rgb(0,100,0, maxColorValue = 255),
                               rgb(245,242,170, maxColorValue = 255)), bias=0.5)

breaks <- c(-0.334155,-0.250606,-0.188895,-0.143314,-0.109647,-0.084780,-0.066413,-0.052846,-0.042826,-0.035424,-0.029958,-0.022556,-0.012536,0.001031,0.019398,0.044265,0.077932,0.123513,0.185224,0.268773,0.381888)
breaks[1] <- round(minValue(temp1),3)-0.001
breaks[length(breaks)] <- round(maxValue(temp1),3)

arg <- list(at=c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), labels=round(c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), 2))
sel <- arg$labels >= 0
arg$labels <- as.character(arg$labels)
arg$labels[sel] <- paste(" ", arg$labels[sel], sep="")
arg$labels

plot(temp1, breaks=breaks, col=rev(my.ramp2(20)), axes=FALSE, main = "Temperateness in the Late Glacial\n15,000-11,700 cal BP", legend=FALSE, cex.main=1.2)
plot(temp1, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.84,0.87,0.05,0.84), col=rev(my.ramp2(20)))
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
points(period1[, c("x", "y")], cex=(decostand(period1$PCO1b, "range")*2)+0.8)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="white", lend=2)
text(2350000, 6050000, "100 km", col="white", cex=1.1, pos=3, font=2)
text(2300000, 6460000, "PL")
text(1860000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

#period 2
my.ramp2 <- colorRampPalette(c(rgb(0,100,0, maxColorValue = 255),
                               rgb(245,242,170, maxColorValue = 255)), bias=0.5)

breaks <- c(-0.211428,-0.145665,-0.098452,-0.064556,-0.040221,-0.022750,-0.010207,-0.001202,0.005264,0.009905,0.013237,0.016570,0.021211,0.027676,0.036682,0.049225,0.066696,0.091031,0.124927,0.172140,0.237902)
breaks[1] <- round(minValue(temp2),3)-0.001
breaks[length(breaks)] <- round(maxValue(temp2),2)

arg <- list(at=c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), labels=round(c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), 2))
sel <- arg$labels >= 0
arg$labels <- as.character(arg$labels)
arg$labels[sel] <- paste(" ", arg$labels[sel], sep="")
arg$labels

plot(temp2, breaks=breaks, col=rev(my.ramp2(20)), axes=FALSE, main = "Temperateness at the Holocene onset\n11,000-10,000 cal BP", legend=FALSE, cex.main=1.2)
plot(temp2, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.84,0.87,0.05,0.84), col=rev(my.ramp2(20)))
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
points(period2[, c("x", "y")], cex=(decostand(period2$PCO_1, "range")*2)+0.8)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="white", lend=2)
text(2350000, 6050000, "100 km", col="white", cex=1.1, pos=3, font=2)
text(2300000, 6460000, "PL")
text(1860000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

#period 3
my.ramp2 <- colorRampPalette(c(rgb(255,255,128, maxColorValue = 255), 
                               rgb(139,69,19, maxColorValue = 255)), bias=1)

breaks <- c(-0.361443,-0.279597,-0.217987,-0.171612,-0.136703,-0.110426,-0.090646,-0.075757,-0.064549,-0.056113,-0.044905,-0.030016,-0.010236,0.016041,0.050950,0.097326,0.158935,0.240782,0.349513,0.493961,0.685858)
breaks[1] <- round(minValue(step3),3)-0.001
breaks[length(breaks)] <- round(maxValue(step3),2)

arg <- list(at=c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), labels=round(c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), 2))
sel <- arg$labels >= 0
arg$labels <- as.character(arg$labels)
arg$labels[sel] <- paste(" ", arg$labels[sel], sep="")
arg$labels

plot(step3, breaks=breaks, col=my.ramp2(20), axes=FALSE, main = "Landscape openness in the Pre-Neolithic\n8,500-7,500 cal BP", legend=FALSE, cex.main=1.2)
plot(step3, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.84,0.87,0.05,0.84), col=my.ramp2(20))
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
points(period3[, c("x", "y")], cex= (decostand(period3$PCO1b, "range")*2)+0.8)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="white", lend=2)
text(2350000, 6050000, "100 km", col="white", cex=1.1, pos=3, font=2)
text(2300000, 6460000, "PL")
text(1860000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

#period 4
my.ramp2 <- colorRampPalette(c(rgb(255,255,128, maxColorValue = 255), 
                               rgb(139,69,19, maxColorValue = 255)), bias=1)


breaks <- c(-0.214518,-0.141825,-0.095241,-0.065388,-0.046258,-0.033998,-0.026142,-0.021108,-0.017881,-0.015814,-0.014489,-0.012422,-0.009195,-0.004161,0.003695,0.015955,0.035085,0.064937,0.111521,0.184215,0.297650)#
breaks[1] <- round(minValue(step4),3)-0.001
breaks[length(breaks)] <- round(maxValue(step4),2)

arg <- list(at=c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), labels=round(c(breaks[c(1,2,4)], 0, breaks[c(18,20,21)]), 2))
# arg$at[length(arg$at)] <- 0.29
sel <- arg$labels >= 0
arg$labels <- as.character(arg$labels)
arg$labels[sel] <- paste(" ", arg$labels[sel], sep="")
arg$labels

plot(step4, breaks=breaks, col=my.ramp2(20), axes=FALSE, main = "Landscape openness in the Late Neolithic\n6,500-5,500 cal BP", legend=FALSE, cex.main=1.2)
plot(step4, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.84,0.87,0.05,0.84), col=my.ramp2(20))
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
points(period4[, c("x", "y")], cex= (decostand(period4$PCO1b, "range")*2)+0.8)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="white", lend=2)
text(2350000, 6050000, "100 km", col="white", cex=1.1, pos=3, font=2)
text(2300000, 6460000, "PL")
text(1860000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

#period 5
my.ramp2 <- colorRampPalette(c("turquoise4", "#BCEBEC",
                               "lightcyan"), bias=1)

breaks <- c(-0.160219,-0.135425,-0.114908,-0.097930,-0.083882,-0.072256,-0.062637,-0.054676,-0.048089,-0.042638,-0.036051,-0.028091,-0.018471,-0.006846,0.007203,0.024180,0.044697,0.069492,0.099455,0.135664,0.179423)
breaks[1] <- round(minValue(taiga5),3)-0.001
breaks[length(breaks)] <- round(maxValue(taiga5),2)

arg <- list(at=c(breaks[c(1,3,6)], 0, breaks[c(18,20,21)]), labels=round(c(breaks[c(1,3,6)], 0, breaks[c(18,20,21)]), 2))
sel <- arg$labels >= 0
arg$labels <- as.character(arg$labels)
arg$labels[sel] <- paste(" ", arg$labels[sel], sep="")
arg$labels

plot(taiga5, breaks=breaks, col=rev(my.ramp2(20)), axes=FALSE, main = "Representation of taiga in the Late Prehistory\n3,000-2,000 cal BP", legend=FALSE, cex.main=1.2)
plot(taiga5, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.84,0.87,0.05,0.84), col=rev(my.ramp2(20)))
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
points(period5[, c("x", "y")], cex= (decostand(period5$PCO_1, "range")*2)+0.8)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="white", lend=2)
text(2350000, 6050000, "100 km", col="white", cex=1.1, pos=3, font=2)
text(2300000, 6460000, "PL")
text(1860000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

# dev.off()
