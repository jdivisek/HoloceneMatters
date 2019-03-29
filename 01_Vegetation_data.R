#################################################################################
#		IMPORT VEGETATION PLOT DATA AND MAP SPECIES RICHNESS		#
#################################################################################


###IMPORT DATA----------------------------------------------------------------
#already selected plots

stijeh <- read.delim("stijeh.sel.f.txt", header=T, row.names = 1) #dark conifereous forests
svetle <- read.delim("svetle.sel.f.txt", header=T, row.names = 1) #light foresst
travniky <- read.delim("travniky.sel.f.txt", header=T, row.names = 1) #semi-dry and steppe grasslands

###MAPS FOR PAPER-------------------------------------------------------------
library(raster)
library(rgdal)
library(berryFunctions)
library(classInt)
library(BAMMtools)

hillshd <- raster(paste(getwd(), "/hillshade/hillshade.tif", sep=""))

shades <- colorRampPalette(c("gray80", "white"))

##Maps of total species richness----------------------------------------------
##Dark conifereous forests----------------------------------------------------

# windows(9.89, 6.55)
tiff("Plotmap_stijeh_tot.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- stijeh.sel.f[, c("tot_rich", "POINT_X", "POINT_Y")]
plot.data <- plot.data[order(plot.data$tot_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Total species richness of dark coniferous forests")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$tot_rich, k=7+1)
clas <- cut(plot.data$tot_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$tot_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1930000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")


dev.off()

##Light forests------------------------------------------------------

# windows(9.89, 6.55)
tiff("Plotmap_svetle_tot.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- svetle.sel.f[, c("tot_rich", "POINT_X", "POINT_Y")]
plot.data <- plot.data[order(plot.data$tot_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Total species richness of light forests")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$tot_rich, k=7+1)
clas <- cut(plot.data$tot_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$tot_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1930000, 6350000, "CZ")
text(2180000, 6248000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

dev.off()

##Grasslands------------------------------------------------------

# windows(9.89, 6.55)
tiff("Plotmap_travniky_tot.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- travniky.sel.f[, c("Pocet_druhu_celkem", "POINT_X", "POINT_Y")]
colnames(plot.data)[1] <- "tot_rich"
plot.data <- plot.data[order(plot.data$tot_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Total species richness of semi-dry and steppe grasslands")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$tot_rich, k=7+1)
clas <- cut(plot.data$tot_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$tot_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1920000, 6350000, "CZ")
text(2180000, 6248000, "SK")
text(2200000, 6050000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

dev.off()

##Maps of diagnostic species richness----------------------------------------------
##Stinne jehlicnate lesy------------------------------------------------------
# windows(9.89, 6.55)
tiff("Plotmap_stijeh_dg.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- stijeh.sel.f[, c("dg_rich", "POINT_X", "POINT_Y")]
plot.data <- plot.data[order(plot.data$dg_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Specialist species richness of dark coniferous forests")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$dg_rich, k=7+1)
clas <- cut(plot.data$dg_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$dg_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1930000, 6350000, "CZ")
text(2100000, 6190000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

dev.off()

##Svetle lesy------------------------------------------------------
# windows(9.89, 6.55)
tiff("Plotmap_svetle_dg.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- svetle.sel.f[, c("dg_rich", "POINT_X", "POINT_Y")]
plot.data <- plot.data[order(plot.data$dg_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Specialist species richness of light forests")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$dg_rich, k=7+1)
clas <- cut(plot.data$dg_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$dg_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1930000, 6350000, "CZ")
text(2180000, 6248000, "SK")
text(2200000, 6070000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

dev.off()

##Travniky------------------------------------------------------
# windows(9.89, 6.55)
tiff("Plotmap_travniky_dg.rich.tif", 9.89, 6.55, units = "in", res = 400, compression = "lzw")
par(mar=c(3,4.5,3,0))

plot.data <- travniky.sel.f[, c("Pocet_vybranych_druhu", "POINT_X", "POINT_Y")]
colnames(plot.data)[1] <- "dg_rich"
plot.data <- plot.data[order(plot.data$dg_rich),]

plot(hillshd, col=shades(100), asp=1, axes=F, legend=FALSE,
     main = "Specialist species richness of semi-dry and steppe grasslands")
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)

my.col <- my.ramp(7)
breaks <- getJenksBreaks(plot.data$dg_rich, k=7+1)
clas <- cut(plot.data$dg_rich, breaks, include.lowest = TRUE, labels=1:(length(breaks)-1))

points(plot.data$POINT_X, plot.data$POINT_Y, pch=21, col="black", cex=0.9, bg=my.col[clas], lwd=0.5)

colPointsLegend(plot.data$dg_rich, colors = my.col, title="Number of species",
                bb = breaks, nlab = length(breaks), at=breaks, labels = breaks, cex=0.8,
                x1 = 0.60, x2 = 0.89, y1 = 0.80, y2 = 0.90, lines = F, bg=NA)

lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)
text(2300000, 6430000, "PL")
text(1920000, 6350000, "CZ")
text(2180000, 6248000, "SK")
text(2200000, 6050000, "HU")
text(1840000, 6150000, "AT")
text(2550000, 6220000, "UA")

dev.off()