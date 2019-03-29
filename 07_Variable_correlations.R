#########################################################################
#			Variable correlations				#
#########################################################################

library(corrplot)

###Light forests-------------------------------------------------------------------------
tiff("Correlation_matrix.tif", 9.62, 9.02, units = "in", res = 400, compression = "lzw")
cor.mat <- cor(svetle[, c(env.vars[-1], hist.vars[-6], "Releve_area")])
colnames(cor.mat) <- c("Terrain ruggedness index", "Heat load index", "Topographic wetness index", "Soil pH", "Mean annual temperature",
                       "Temperature seasonality", "Annual precipitation", "Precipitation seasonality", "Area of forests",
                       "Area of grasslands", "Area of arable land", "Build-up area", "Landscape diversity", "Temperateness in Late Glacial",
                       "Temperateness at Holocene onset", "Landscape openness in Pre-Neolithic", "Landscape openness in Late Neolithic", "Representation of taiga in Late Prehistory",
                       "Plot size")
rownames(cor.mat) <- colnames(cor.mat)


cor.col  <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

corrplot.mixed(cor.mat, lower = "number", upper = "ellipse", tl.pos = "lt",
               tl.cex=0.7, tl.col = "black", number.cex=0.7,
               title="Light forests", upper.col=rev(cor.col(200)), 
               lower.col = rev(cor.col(200)), mar=c(0,0,1.5,0))
dev.off()

###Semi-dry and steppe grasslands--------------------------------------------------------
tiff("Correlation_matrix.tif", 9.62, 9.02, units = "in", res = 400, compression = "lzw")
cor.mat <- cor(travniky[, c(env.vars[-1], hist.vars[-6], "Releve_area")])
colnames(cor.mat) <- c("Terrain ruggedness index", "Heat load index", "Topographic wetness index", "Soil pH", "Mean annual temperature",
                       "Temperature seasonality", "Annual precipitation", "Precipitation seasonality", "Area of forests",
                       "Area of grasslands", "Area of arable land", "Build-up area", "Landscape diversity", "Temperateness in Late Glacial",
                       "Temperateness at Holocene onset", "Landscape openness in Pre-Neolithic", "Landscape openness in Late Neolithic", "Representation of taiga in Late Prehistory",
                       "Plot size")
rownames(cor.mat) <- colnames(cor.mat)


cor.col  <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

corrplot.mixed(cor.mat, lower = "number", upper = "ellipse", tl.pos = "lt",
               tl.cex=0.7, tl.col = "black", number.cex=0.7,
               title="Semi-dry and steppe grasslands", upper.col=rev(cor.col(200)), 
               lower.col = rev(cor.col(200)), mar=c(0,0,1.5,0))
dev.off()

###Dark coniferous forests---------------------------------------------------------------
tiff("Correlation_matrix.tif", 9.62, 9.02, units = "in", res = 400, compression = "lzw")

cor.mat <- cor(stijeh[, c(env.vars[-1], hist.vars[-6], "Releve_area")])
colnames(cor.mat) <- c("Terrain ruggedness index", "Heat load index", "Topographic wetness index", "Soil pH", "Mean annual temperature",
                       "Temperature seasonality", "Annual precipitation", "Precipitation seasonality", "Area of forests",
                       "Area of grasslands", "Area of arable land", "Build-up area", "Landscape diversity", "Temperateness in Late Glacial",
                       "Temperateness at Holocene onset", "Landscape openness in Pre-Neolithic", "Landscape openness in Late Neolithic", "Representation of taiga in Late Prehistory",
                       "Plot size")
rownames(cor.mat) <- colnames(cor.mat)


cor.col  <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                               "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                               "#4393C3", "#2166AC", "#053061"))

corrplot.mixed(cor.mat, lower = "number", upper = "ellipse", tl.pos = "lt",
               tl.cex=0.7, tl.col = "black", number.cex=0.7,
               title="Dark coniferous forests", upper.col=rev(cor.col(200)), 
               lower.col = rev(cor.col(200)), mar=c(0,0,1.5,0))
dev.off()
