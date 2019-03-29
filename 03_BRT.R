######################################################################################
#				MODEL SPECIES RICHNESS				     #
######################################################################################

library(dismo)
library(raster)
library(gbm)
library(viridis)
library(pgirmess)

#import environmental data
r.paths <- list.files(path=paste(getwd(), "/rasters", sep=""), pattern='tif', full.names=TRUE ) 
r.paths <- r.paths[c(19,14,20,16,2,5,3,4,12,13,11,21,7)]

env <- stack(r.paths)

#import historical data
r.paths <- list.files(path=paste(getwd(), "/history_rasters", sep=""), pattern='tif', full.names=TRUE ) 

paleo <- stack(r.paths)

#merge datasets
env.hist <- stack(env, paleo)
names(env.hist)

#import habitat ranges
r.paths <- list.files(path=paste(getwd(), "/range_rasters", sep=""), pattern='tif', full.names=TRUE )
ranges <-stack(r.paths)
plot(ranges)

#colors for maps
my.ramp <- colorRampPalette(c(rgb(40,146,199, maxColorValue = 255), 
                              rgb(250,250,100, maxColorValue = 255),
                              rgb(232,16,20, maxColorValue = 255)))

###DARK CONIFEROUS FORESTS----------------------------------------------------------------
Releve_area <- 400
Releve_area <- as.data.frame(Releve_area)

###all species---------------------------------------------------------------------------
#current environment only
set.seed(1234)
mod1.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1.stijeh$gbm.call$best.trees
hist(mod1.stijeh$residuals)
mod1.stijeh$contributions

cor(stijeh$tot_rich, mod1.stijeh$fitted)^2 #pseudo-R2
mean(mod1.stijeh$residuals * mod1.stijeh$residuals) #MSE
mod1.stijeh$cv.statistics$deviance.mean #minimum cv deviance
mod1.stijeh$cv.statistics$deviance.se #cv deviance se

pred1.stijeh <- predict(env, mod1.stijeh, const=Releve_area, n.trees=mod1.stijeh$gbm.call$best.trees, type="response")
plot(pred1.stijeh, main="Dark coniferous forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred1.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1.stijeh <- mask(pred1.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred1.stijeh <- round(pred1.stijeh)
plot(pred1.stijeh, main="Dark coniferous forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred1.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred1.stijeh.mask.tif", sep=""))
crs(pred1.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod1.stijeh$residuals)
plot(mor.cor)

#historical factors
set.seed(1234)
mod2.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22.stijeh <-  gbm.step(data=stijeh, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "tot_rich", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod2.stijeh$gbm.call$best.trees
hist(mod2.stijeh$residuals)
mod22.stijeh$contributions

cor(stijeh$tot_rich, mod2.stijeh$fitted)^2 #pseudo-R2
mean(mod2.stijeh$residuals * mod2.stijeh$residuals) #MSE
mod2.stijeh$cv.statistics$deviance.mean #min cv deviance
mod2.stijeh$cv.statistics$deviance.se #cv deviance se

pred2.stijeh <- predict(env.hist, mod2.stijeh, const=Releve_area, n.trees=mod2.stijeh$gbm.call$best.trees, type="response")
plot(pred2.stijeh, main="Dark coniferous forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred2.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2.stijeh <- mask(pred2.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred2.stijeh <- round(pred2.stijeh)
plot(pred2.stijeh, main="Dark coniferous forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred2.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred2.stijeh.mask.tif", sep=""))
crs(pred2.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod2.stijeh$residuals)
plot(mor.cor)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1.stijeh$gbm.call$best.trees, mod2.stijeh$gbm.call$best.trees, mod22.stijeh$gbm.call$best.trees)))
y.lim <- c(min(c(mod1.stijeh$cv.values - mod1.stijeh$cv.loss.ses,
                 mod2.stijeh$cv.values - mod2.stijeh$cv.loss.ses,
                 mod22.stijeh$cv.values - mod22.stijeh$cv.loss.ses)),
           max(c(mod1.stijeh$cv.values + mod1.stijeh$cv.loss.ses,
                 mod2.stijeh$cv.values + mod2.stijeh$cv.loss.ses,
                 mod22.stijeh$cv.values + mod22.stijeh$cv.loss.ses)))

plot(mod1.stijeh$trees.fitted, mod1.stijeh$cv.values - mod1.stijeh$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of dark coniferous forests\nd = 5, lr = 0.001")

polygon(c(mod1.stijeh$trees.fitted, rev(mod1.stijeh$trees.fitted)),
        c(mod1.stijeh$cv.values - mod1.stijeh$cv.loss.ses, rev(mod1.stijeh$cv.values + mod1.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1.stijeh$trees.fitted, mod1.stijeh$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22.stijeh$trees.fitted, rev(mod22.stijeh$trees.fitted)),
        c(mod22.stijeh$cv.values - mod22.stijeh$cv.loss.ses, rev(mod22.stijeh$cv.values + mod22.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22.stijeh$trees.fitted, mod22.stijeh$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2.stijeh$trees.fitted, rev(mod2.stijeh$trees.fitted)),
        c(mod2.stijeh$cv.values - mod2.stijeh$cv.loss.ses, rev(mod2.stijeh$cv.values + mod2.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2.stijeh$trees.fitted, mod2.stijeh$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1.stijeh$cv.values, mod2.stijeh$cv.values, mod22.stijeh$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###habitat specialists------------------------------------------------------------------
#current environment only
set.seed(1234)
mod3.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3.stijeh$gbm.call$best.trees
# hist(mod3.stijeh$residuals)
mod3.stijeh$contributions

cor(stijeh$dg_rich, mod3.stijeh$fitted)^2 #pseudo-R2
mean(mod3.stijeh$residuals * mod3.stijeh$residuals) #MSE
mod3.stijeh$cv.statistics$deviance.mean #min cv deviance
mod3.stijeh$cv.statistics$deviance.se #cv deviance se

pred3.stijeh <- predict(env, mod3.stijeh, const=Releve_area, n.trees=mod3.stijeh$gbm.call$best.trees, type="response")
plot(pred3.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred3.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3.stijeh <- mask(pred3.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred3.stijeh <- round(pred3.stijeh)
plot(pred3.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred3.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred3.stijeh.mask.tif", sep=""))
crs(pred3.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod3.stijeh$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod4.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42.stijeh <-  gbm.step(data=stijeh, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "dg_rich", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod4.stijeh$gbm.call$best.trees
hist(mod4.stijeh$residuals)
mod4.stijeh$contributions

cor(stijeh$dg_rich, mod4.stijeh$fitted)^2 #pseudo-R2
mean(mod4.stijeh$residuals * mod4.stijeh$residuals) #MSE
mod4.stijeh$cv.statistics$deviance.mean #min cv deviance
mod4.stijeh$cv.statistics$deviance.se #cv deviance se

pred4.stijeh <- predict(env.hist, mod4.stijeh, const=Releve_area, n.trees=mod4.stijeh$gbm.call$best.trees, type="response")
plot(pred4.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred4.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4.stijeh <- mask(pred4.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred4.stijeh <- round(pred4.stijeh)
plot(pred4.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred4.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred4.stijeh.mask.tif", sep=""))
crs(pred4.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod4.stijeh$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3.stijeh$gbm.call$best.trees, mod4.stijeh$gbm.call$best.trees, mod42.stijeh$gbm.call$best.trees)))
y.lim <- c(min(c(mod3.stijeh$cv.values - mod3.stijeh$cv.loss.ses,
                 mod4.stijeh$cv.values - mod4.stijeh$cv.loss.ses,
                 mod42.stijeh$cv.values - mod42.stijeh$cv.loss.ses)),
           max(c(mod3.stijeh$cv.values + mod3.stijeh$cv.loss.ses,
                 mod4.stijeh$cv.values + mod4.stijeh$cv.loss.ses,
                 mod42.stijeh$cv.values + mod42.stijeh$cv.loss.ses)))

plot(mod3.stijeh$trees.fitted, mod3.stijeh$cv.values - mod3.stijeh$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of dark coniferous forests\nd = 5, lr = 0.001")

polygon(c(mod3.stijeh$trees.fitted, rev(mod3.stijeh$trees.fitted)),
        c(mod3.stijeh$cv.values - mod3.stijeh$cv.loss.ses, rev(mod3.stijeh$cv.values + mod3.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3.stijeh$trees.fitted, mod3.stijeh$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42.stijeh$trees.fitted, rev(mod42.stijeh$trees.fitted)),
        c(mod42.stijeh$cv.values - mod42.stijeh$cv.loss.ses, rev(mod42.stijeh$cv.values + mod42.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42.stijeh$trees.fitted, mod42.stijeh$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4.stijeh$trees.fitted, rev(mod4.stijeh$trees.fitted)),
        c(mod4.stijeh$cv.values - mod4.stijeh$cv.loss.ses, rev(mod4.stijeh$cv.values + mod4.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4.stijeh$trees.fitted, mod4.stijeh$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3.stijeh$cv.values, mod4.stijeh$cv.values, mod42.stijeh$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###LIGHT FORESTS-------------------------------------------------------------------------
###all species---------------------------------------------------------------------------
#current environment only
set.seed(1234)
mod1.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1.svetle$gbm.call$best.trees
hist(mod1.svetle$residuals)
mod1.svetle$contributions

cor(svetle$tot_rich, mod1.svetle$fitted)^2 #pseudo-R2
mean(mod1.svetle$residuals * mod1.svetle$residuals) #MSE
mod1.svetle$cv.statistics$deviance.mean #min cv deviance
mod1.svetle$cv.statistics$deviance.se #cv deviance se

pred1.svetle <- predict(env, mod1.svetle, const=Releve_area, n.trees=mod1.svetle$gbm.call$best.trees, type="response")
plot(pred1.svetle, main="Light forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred1.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1.svetle <- mask(pred1.svetle, ranges[["list_1km"]], maskvalue=0)
pred1.svetle <- round(pred1.svetle)
plot(pred1.svetle, main="Light forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred1.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred1.svetle.mask.tif", sep=""))
crs(pred1.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod1.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod2.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22.svetle <-  gbm.step(data=svetle, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "tot_rich", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod2.svetle$gbm.call$best.trees
hist(mod2.svetle$residuals)
mod2.svetle$contributions

cor(svetle$tot_rich, mod2.svetle$fitted)^2 #pseudo-R2
mean(mod2.svetle$residuals * mod2.svetle$residuals) #MSE
mod2.svetle$cv.statistics$deviance.mean #min cv deviance
mod2.svetle$cv.statistics$deviance.se #cv deviance se

pred2.svetle <- predict(env.hist, mod2.svetle, const=Releve_area, n.trees=mod2.svetle$gbm.call$best.trees, type="response")
plot(pred2.svetle, main="Light forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred2.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2.svetle <- mask(pred2.svetle, ranges[["list_1km"]], maskvalue=0)
pred2.svetle <- round(pred2.svetle)
plot(pred2.svetle, main="Light forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred2.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred2.svetle.mask.tif", sep=""))
crs(pred2.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod2.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1.svetle$gbm.call$best.trees, mod2.svetle$gbm.call$best.trees, mod22.svetle$gbm.call$best.trees)))
y.lim <- c(min(c(mod1.svetle$cv.values - mod1.svetle$cv.loss.ses,
                 mod2.svetle$cv.values - mod2.svetle$cv.loss.ses,
                 mod22.svetle$cv.values - mod22.svetle$cv.loss.ses)),
           max(c(mod1.svetle$cv.values + mod1.svetle$cv.loss.ses,
                 mod2.svetle$cv.values + mod2.svetle$cv.loss.ses,
                 mod22.svetle$cv.values + mod22.svetle$cv.loss.ses)))

plot(mod1.svetle$trees.fitted, mod1.svetle$cv.values - mod1.svetle$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of light forests\nd = 5, lr = 0.001")

polygon(c(mod1.svetle$trees.fitted, rev(mod1.svetle$trees.fitted)),
        c(mod1.svetle$cv.values - mod1.svetle$cv.loss.ses, rev(mod1.svetle$cv.values + mod1.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1.svetle$trees.fitted, mod1.svetle$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22.svetle$trees.fitted, rev(mod22.svetle$trees.fitted)),
        c(mod22.svetle$cv.values - mod22.svetle$cv.loss.ses, rev(mod22.svetle$cv.values + mod22.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22.svetle$trees.fitted, mod22.svetle$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2.svetle$trees.fitted, rev(mod2.svetle$trees.fitted)),
        c(mod2.svetle$cv.values - mod2.svetle$cv.loss.ses, rev(mod2.svetle$cv.values + mod2.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2.svetle$trees.fitted, mod2.svetle$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1.svetle$cv.values, mod2.svetle$cv.values, mod22.svetle$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###diagnostic species------------------------------------------------------------------
#current environment only
set.seed(1234)
mod3.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3.svetle$gbm.call$best.trees
hist(mod3.svetle$residuals)
mod3.svetle$contributions

cor(svetle$dg_rich, mod3.svetle$fitted)^2 #pseudo-R2
mean(mod3.svetle$residuals * mod3.svetle$residuals) #MSE
mod3.svetle$cv.statistics$deviance.mean #min cv deviance
mod3.svetle$cv.statistics$deviance.se #cv deviance se

pred3.svetle <- predict(env, mod3.svetle, const=Releve_area, n.trees=mod3.svetle$gbm.call$best.trees, type="response")
plot(pred3.svetle, main="Light forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred3.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3.svetle <- mask(pred3.svetle, ranges[["list_1km"]], maskvalue=0)
pred3.svetle <- round(pred3.svetle)
plot(pred3.svetle, main="Light forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred3.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred3.svetle.mask.tif", sep=""))
crs(pred3.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod3.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod4.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42.svetle <-  gbm.step(data=svetle, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "dg_rich", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod4.svetle$gbm.call$best.trees
hist(mod4.svetle$residuals)
mod4.svetle$contributions

cor(svetle$dg_rich, mod4.svetle$fitted)^2 #pseudo-R2
mean(mod4.svetle$residuals * mod4.svetle$residuals) #MSE
mod4.svetle$cv.statistics$deviance.mean #min cv deviance
mod4.svetle$cv.statistics$deviance.se #cv deviance se

pred4.svetle <- predict(env.hist, mod4.svetle, const=Releve_area, n.trees=mod4.svetle$gbm.call$best.trees, type="response")
plot(pred4.svetle, main="Light forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred4.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4.svetle <- mask(pred4.svetle, ranges[["list_1km"]], maskvalue=0)
pred4.svetle <- round(pred4.svetle)
plot(pred4.svetle, main="Light forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred4.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred4.svetle.mask.tif", sep=""))
crs(pred4.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod4.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3.svetle$gbm.call$best.trees, mod4.svetle$gbm.call$best.trees, mod42.svetle$gbm.call$best.trees)))
y.lim <- c(min(c(mod3.svetle$cv.values - mod3.svetle$cv.loss.ses,
                 mod4.svetle$cv.values - mod4.svetle$cv.loss.ses,
                 mod42.svetle$cv.values - mod42.svetle$cv.loss.ses)),
           max(c(mod3.svetle$cv.values + mod3.svetle$cv.loss.ses,
                 mod4.svetle$cv.values + mod4.svetle$cv.loss.ses,
                 mod42.svetle$cv.values + mod42.svetle$cv.loss.ses)))

plot(mod3.svetle$trees.fitted, mod3.svetle$cv.values - mod3.svetle$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of light forests\nd = 5, lr = 0.001")

polygon(c(mod3.svetle$trees.fitted, rev(mod3.svetle$trees.fitted)),
        c(mod3.svetle$cv.values - mod3.svetle$cv.loss.ses, rev(mod3.svetle$cv.values + mod3.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3.svetle$trees.fitted, mod3.svetle$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42.svetle$trees.fitted, rev(mod42.svetle$trees.fitted)),
        c(mod42.svetle$cv.values - mod42.svetle$cv.loss.ses, rev(mod42.svetle$cv.values + mod42.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42.svetle$trees.fitted, mod42.svetle$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4.svetle$trees.fitted, rev(mod4.svetle$trees.fitted)),
        c(mod4.svetle$cv.values - mod4.svetle$cv.loss.ses, rev(mod4.svetle$cv.values + mod4.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4.svetle$trees.fitted, mod4.svetle$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3.svetle$cv.values, mod4.svetle$cv.values, mod42.svetle$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###SEMI-DRY AND STEPPE GRASSLANDS------------------------------------------------------
Releve_area <- 25
Releve_area <- as.data.frame(Releve_area)

#all species---------------------------------------------------------------------------
#current environment only

set.seed(1234)
mod1.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], "Releve_area"), 
                          gbm.y = "Pocet_druhu_celkem", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1.travniky$gbm.call$best.trees
hist(mod1.travniky$residuals)
mod1.travniky$contributions

cor(travniky$Pocet_druhu_celkem, mod1.travniky$fitted)^2
mean(mod1.travniky$residuals * mod1.travniky$residuals) #MSE
mod1.travniky$cv.statistics$deviance.mean #min cv deviance
mod1.travniky$cv.statistics$deviance.se #cv deviance se

pred1.travniky <- predict(env, mod1.travniky, const=Releve_area, n.trees=mod1.travniky$gbm.call$best.trees, type="response")
plot(pred1.travniky, main="Travniky (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred1.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1.travniky <- mask(pred1.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred1.travniky <- round(pred1.travniky)
plot(pred1.travniky, main="Travniky (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred1.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred1.travniky.mask.tif", sep=""))
crs(pred1.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod1.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod2.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                          gbm.y = "Pocet_druhu_celkem", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22.travniky <- gbm.step(data=travniky, gbm.x = c(hist.vars[-6], "Releve_area"), 
                           gbm.y = "Pocet_druhu_celkem", family = "poisson",
                           tree.complexity = 5, learning.rate = 0.003, 
                           bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod2.travniky$gbm.call$best.trees
hist(mod2.travniky$residuals)
mod2.travniky$contributions

cor(travniky$Pocet_druhu_celkem, mod2.travniky$fitted)^2 #pseudo-R2
mean(mod2.travniky$residuals * mod2.travniky$residuals) #MSE
mod2.travniky$cv.statistics$deviance.mean #min cv deviance
mod2.travniky$cv.statistics$deviance.se #cv deviance se 

pred2.travniky <- predict(env.hist, mod2.travniky, const=Releve_area, n.trees=mod2.travniky$gbm.call$best.trees, type="response")
plot(pred2.travniky, main="Travniky (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred2.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2.travniky <- mask(pred2.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred2.travniky <- round(pred2.travniky)
plot(pred2.travniky, main="Travniky (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred2.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred2.travniky.mask.tif", sep=""))
crs(pred2.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod2.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1.travniky$gbm.call$best.trees, mod2.travniky$gbm.call$best.trees, mod22.travniky$gbm.call$best.trees)))
y.lim <- c(min(c(mod1.travniky$cv.values - mod1.travniky$cv.loss.ses,
                 mod2.travniky$cv.values - mod2.travniky$cv.loss.ses,
                 mod22.travniky$cv.values - mod22.travniky$cv.loss.ses)),
           max(c(mod1.travniky$cv.values + mod1.travniky$cv.loss.ses,
                 mod2.travniky$cv.values + mod2.travniky$cv.loss.ses,
                 mod22.travniky$cv.values + mod22.travniky$cv.loss.ses)))

plot(mod1.travniky$trees.fitted, mod1.travniky$cv.values - mod1.travniky$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of semi-dry and steppe grasslands\nd = 5, lr = 0.003")

polygon(c(mod1.travniky$trees.fitted, rev(mod1.travniky$trees.fitted)),
        c(mod1.travniky$cv.values - mod1.travniky$cv.loss.ses, rev(mod1.travniky$cv.values + mod1.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1.travniky$trees.fitted, mod1.travniky$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22.travniky$trees.fitted, rev(mod22.travniky$trees.fitted)),
        c(mod22.travniky$cv.values - mod22.travniky$cv.loss.ses, rev(mod22.travniky$cv.values + mod22.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22.travniky$trees.fitted, mod22.travniky$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2.travniky$trees.fitted, rev(mod2.travniky$trees.fitted)),
        c(mod2.travniky$cv.values - mod2.travniky$cv.loss.ses, rev(mod2.travniky$cv.values + mod2.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2.travniky$trees.fitted, mod2.travniky$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1.travniky$cv.values, mod2.travniky$cv.values, mod22.travniky$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

####diagnostic species--------------------------------------------------
#current environment only
set.seed(1234)
mod3.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], "Releve_area"), 
                          gbm.y = "Pocet_vybranych_druhu", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3.travniky$gbm.call$best.trees
hist(mod3.travniky$residuals)
mod3.travniky$contributions

cor(travniky$Pocet_vybranych_druhu, mod3.travniky$fitted)^2 #pseudo-R2
mean(mod3.travniky$residuals * mod3.travniky$residuals) #MSE
mod3.travniky$cv.statistics$deviance.mean #min cv deviance
mod3.travniky$cv.statistics$deviance.se #cv deviance se

pred3.travniky <- predict(env, mod3.travniky, const=Releve_area, n.trees=mod3.travniky$gbm.call$best.trees, type="response")
plot(pred3.travniky, main="Travniky (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred3.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3.travniky <- mask(pred3.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred3.travniky <- round(pred3.travniky)
plot(pred3.travniky, main="Travniky (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred3.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred3.travniky.mask.tif", sep=""))
crs(pred3.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod3.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod4.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                          gbm.y = "Pocet_vybranych_druhu", family = "poisson",
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42.travniky <- gbm.step(data=travniky, gbm.x = c(hist.vars[-6], "Releve_area"), 
                           gbm.y = "Pocet_vybranych_druhu", family = "poisson",
                           tree.complexity = 5, learning.rate = 0.003, 
                           bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod4.travniky$gbm.call$best.trees
hist(mod4.travniky$residuals)
mod4.travniky$contributions

cor(travniky$Pocet_vybranych_druhu, mod4.travniky$fitted)^2 #pseudo-R2
mean(mod4.travniky$residuals * mod4.travniky$residuals) #MSE
mod4.travniky$cv.statistics$deviance.mean #min cv deviance
mod4.travniky$cv.statistics$deviance.se #cv deviance se

pred4.travniky <- predict(env.hist, mod4.travniky, const=Releve_area, n.trees=mod4.travniky$gbm.call$best.trees, type="response")
plot(pred4.travniky, main="Travniky (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred4.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4.travniky <- mask(pred4.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred4.travniky <- round(pred4.travniky)
plot(pred4.travniky, main="Travniky (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred4.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred4.travniky.mask.tif", sep=""))
crs(pred4.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod4.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3.travniky$gbm.call$best.trees, mod4.travniky$gbm.call$best.trees, mod42.travniky$gbm.call$best.trees)))
y.lim <- c(min(c(mod3.travniky$cv.values - mod3.travniky$cv.loss.ses,
                 mod4.travniky$cv.values - mod4.travniky$cv.loss.ses,
                 mod42.travniky$cv.values - mod42.travniky$cv.loss.ses)),
           max(c(mod3.travniky$cv.values + mod3.travniky$cv.loss.ses,
                 mod4.travniky$cv.values + mod4.travniky$cv.loss.ses,
                 mod42.travniky$cv.values + mod42.travniky$cv.loss.ses)))

plot(mod3.travniky$trees.fitted, mod3.travniky$cv.values - mod3.travniky$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of semi-dry and steppe grasslands\nd = 5, lr = 0.003")

polygon(c(mod3.travniky$trees.fitted, rev(mod3.travniky$trees.fitted)),
        c(mod3.travniky$cv.values - mod3.travniky$cv.loss.ses, rev(mod3.travniky$cv.values + mod3.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3.travniky$trees.fitted, mod3.travniky$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42.travniky$trees.fitted, rev(mod42.travniky$trees.fitted)),
        c(mod42.travniky$cv.values - mod42.travniky$cv.loss.ses, rev(mod42.travniky$cv.values + mod42.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42.travniky$trees.fitted, mod42.travniky$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4.travniky$trees.fitted, rev(mod4.travniky$trees.fitted)),
        c(mod4.travniky$cv.values - mod4.travniky$cv.loss.ses, rev(mod4.travniky$cv.values + mod4.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4.travniky$trees.fitted, mod4.travniky$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3.travniky$cv.values, mod4.travniky$cv.values, mod42.travniky$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###extract variable contributions from the models--------------------------------
##all species + current factors
contrib1 <- list(mod1.stijeh$contributions[c(env.vars[-1], "Releve_area"),],
                 mod1.svetle$contributions[c(env.vars[-1], "Releve_area"),],
                 mod1.travniky$contributions[c(env.vars[-1], "Releve_area"),]) 
contrib1 <- do.call(cbind.data.frame, contrib1)[,c(2,4,6)]
colnames(contrib1) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib1, 1)

#all species + current & historical factors
contrib2 <- list(mod2.stijeh$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod2.svetle$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod2.travniky$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),]) 
contrib2 <- do.call(cbind.data.frame, contrib2)[,c(2,4,6)]
colnames(contrib2) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib2, 1)

##diagnostic species + current factors
contrib3 <- list(mod3.stijeh$contributions[c(env.vars[-1], "Releve_area"),],
                 mod3.svetle$contributions[c(env.vars[-1], "Releve_area"),],
                 mod3.travniky$contributions[c(env.vars[-1], "Releve_area"),]) 
contrib3 <- do.call(cbind.data.frame, contrib3)[,c(2,4,6)]
colnames(contrib3) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib3, 1)

#diagnostic species + current & historical factors
contrib4 <- list(mod4.stijeh$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod4.svetle$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod4.travniky$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),]) 
contrib4 <- do.call(cbind.data.frame, contrib4)[,c(2,4,6)]
colnames(contrib4) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib4, 1)

###PARTIAL DEPENDENCE PLOTS---------------------------------------------------------------------
###Dark coniferous forests-----------------------------------------------------
mod1.stijeh$contributions

pdp.stijeh <- list()
pdp.stijeh$mod1 <- list()
pdp.stijeh$mod2 <- list()
pdp.stijeh$mod3 <- list()
pdp.stijeh$mod4 <- list()

brt.models <- list(mod1.stijeh, mod2.stijeh, mod3.stijeh, mod4.stijeh)

for(q in c(1:4))#1:4
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdp.stijeh[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")

  }
  names(pdp.stijeh[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)

par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdp.stijeh$mod4[[i]], type="l", lwd=1.3, main=names(pdp.stijeh$mod4)[i],
       ylab="No. species", col="blue")
  points(pdp.stijeh$mod4[[i]], pch=16)
}

plot(pdp.stijeh$mod4[[i]], type="p", lwd=1.3, main=names(pdp.stijeh$mod4)[i],
     ylab="No. species", col="gray90")

###Light forests-------------------------------------------------------------
mod1.svetle$contributions

pdp.svetle <- list()
pdp.svetle$mod1 <- list()
pdp.svetle$mod2 <- list()
pdp.svetle$mod3 <- list()
pdp.svetle$mod4 <- list()

brt.models <- list(mod1.svetle, mod2.svetle, mod3.svetle, mod4.svetle)

for(q in c(1:4))
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdp.svetle[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")
    
  }
  names(pdp.svetle[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)

par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdp.svetle$mod4[[i]], type="l", lwd=1.3, main=names(pdp.svetle$mod4)[i],
       ylab="No. species", col="blue")
  points(pdp.svetle$mod4[[i]], pch=16)
}

###Semi-dry and steppe grasslands-------------------------------------------------------------
mod1.travniky$contributions

pdp.travniky <- list()
pdp.travniky$mod1 <- list()
pdp.travniky$mod2 <- list()
pdp.travniky$mod3 <- list()
pdp.travniky$mod4 <- list()

brt.models <- list(mod1.travniky, mod2.travniky, mod3.travniky, mod4.travniky)

for(q in c(1:4))
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdp.travniky[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")
    
  }
  names(pdp.travniky[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)

par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdp.travniky$mod4[[i]], type="l", lwd=1.3, main=names(pdp.travniky$mod4)[i],
       ylab="No. species", col="blue")
  points(pdp.travniky$mod4[[i]], pch=16)
}

###PREDICTION MAPS-----------------------------------------------------
library(raster)
library(rgdal)
library(berryFunctions)
library(classInt)
library(BAMMtools)

#all species------------------------------------------------------------
# windows(7.44, 7.64)#812, 777
tiff("Map_TSR.tif", 7.44, 7.64, res = 400, units = "in", compression = "lzw")
layout(matrix(1:6, ncol=2, nrow=3, byrow = F))

par(mar=c(1, 0, 3, 0), oma=c(1.5, 3.5, 0, 2.5))

#Light forests
#ArcGIS breaks
breaks <- c(12.000000, 19.684729, 23.527094, 27.369458, 30.955665, 34.541872, 38.384236, 42.995074, 49.399015, 64.000000)
breaks[1] <- min(c(minValue(pred1.svetle), minValue(pred2.svetle)))
breaks[length(breaks)] <- max(c(maxValue(pred1.svetle), maxValue(pred2.svetle)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred1.svetle, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Light forests (TSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred2.svetle, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Light forests (TSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
#axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred2.svetle, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

#Semi-dry and steppe grasslands
#ArcGIS breaks
breaks <- c(21.000000, 26.587045, 30.230769, 33.388664, 36.303644, 39.461538, 43.591093, 48.935223, 56.708502, 80.271255)
breaks[1] <- min(c(minValue(pred1.travniky), minValue(pred2.travniky)))
breaks[length(breaks)] <-  max(c(maxValue(pred1.travniky), maxValue(pred2.travniky)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred1.travniky, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Semi-dry and steppe grasslands (TSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred2.travniky, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Semi-dry and steppe grasslands (TSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred2.travniky, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

#Dark coniferous forests
#ArcGIS breaks
breaks <- c(9.000000, 16.780392, 20.796078, 24.309804, 27.823529, 31.337255, 34.850980, 38.866667, 43.886275, 73.000000)
breaks[1] <- min(c(minValue(pred1.stijeh), minValue(pred2.stijeh)))
breaks[length(breaks)] <- max(c(maxValue(pred1.stijeh), maxValue(pred2.stijeh)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred1.stijeh, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Dark coniferous forests (TSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred2.stijeh, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Dark coniferous forests (TSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred2.stijeh, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))#image.plot {fields}
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

dev.off()

#diagnostic speces------------------------------------------------------------
# windows(7.44, 7.64)#812, 777
tiff("Map_SSR.tif", 7.44, 7.64, res = 400, units = "in", compression = "lzw")
layout(matrix(1:6, ncol=2, nrow=3, byrow = F))

par(mar=c(1, 0, 3, 0), oma=c(1.5, 3.5, 0, 2.5))

#Light forests
#ArcGIS breaks
breaks <- c(2.000000, 3.937255, 4.980392, 6.992157, 8.929412, 10.941176, 12.952941, 14.964706, 16.976471, 21.000000)
breaks[1] <- min(c(minValue(pred3.svetle), minValue(pred4.svetle)))
breaks[length(breaks)] <- max(c(maxValue(pred3.svetle), maxValue(pred4.svetle)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred3.svetle, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Light forests (SSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred4.svetle, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Light forests (SSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred4.svetle, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

#Semi-dry and steppe grasslands
#ArcGIS breaks
breaks <- c(9.000000, 12.917647, 15.964706, 17.996078, 20.898039, 23.945098, 27.862745, 31.925490, 36.858824, 46.000000)
breaks[1] <- min(c(minValue(pred3.travniky), minValue(pred4.travniky)))
breaks[length(breaks)] <-  max(c(maxValue(pred3.travniky), maxValue(pred4.travniky)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred3.travniky, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Semi-dry and steppe grasslands (SSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred4.travniky, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Semi-dry and steppe grasslands (SSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
# axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
# axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred4.travniky, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

#Dark coniferous forests
#breaks
breaks <- seq(1,10)
breaks[1] <- min(c(minValue(pred3.stijeh), minValue(pred4.stijeh)))
breaks[length(breaks)] <- max(c(maxValue(pred3.stijeh), maxValue(pred4.stijeh)))

arg <- list(at=round(breaks, 0), labels=round(breaks, 0), cex.axis = 0.7, tck=-0.4, mgp=c(3, 0.6, 0))

plot(pred3.stijeh, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Dark coniferous forests (SSR ~ E)", legend=FALSE, cex.main=1.1, asp=1)
axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

plot(pred4.stijeh, breaks=breaks, col=my.ramp(length(breaks)-1), axes=FALSE, main = "Dark coniferous forests (SSR ~ E + H)", legend=FALSE, cex.main=1.1, asp=1)
#axis(2, at=c(6446060.617789, 6274652.112374, 6106625.407736), labels = c("50°N", "49°N", "48°N"), las=1)
axis(1, at=c(2448819.53367, 2226182.6357, 2003545.737729), labels = c("22°E", "20°E", "18°E"), las=1)
plot(country, add=T)
plot(pred4.stijeh, breaks=breaks, axis.args=arg, legend.only=TRUE, smallplot= c(0.94,0.97,0.05,0.84), col=my.ramp(length(breaks)-1))#image.plot {fields}
lines(x=c(2300000, 2400000), y=c(6050000, 6050000), lwd=3, col="black", lend=2)
text(2350000, 6050000, "100 km", col="black", cex=1.1, pos=3, font=2)

dev.off()
