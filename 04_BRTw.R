#########################################################################################
#                BOOSTED REGRESSION TREES WITH PLOT WEIGHTS				#
#########################################################################################


library(dismo)
library(raster)
library(gbm)
library(viridis)
library(pgirmess)

###coordinates of pollen profiles
paleo.coords <- rbind(period1[,c("x", "y")],
                      period2[,c("x", "y")],
                      period3[,c("x", "y")],
                      period4[,c("x", "y")],
                      period5[,c("x", "y")])
head(paleo.coords)
dim(paleo.coords)
plot(paleo.coords)
colnames(paleo.coords) <- c("POINT_X", "POINT_Y")


###DARK CONIFEROUS FORESTS----------------------------------------------------------------
#create weights

d <- as.matrix(dist(rbind(stijeh[, c("POINT_X", "POINT_Y")], paleo.coords)))
dim(d)
d <- d[(nrow(stijeh)+1):nrow(d),1:nrow(stijeh)]
d[1:5,1:5]

w.stijeh <- apply(d, 2, FUN=min)
hist(w.stijeh)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(stijeh$POINT_X, stijeh$POINT_Y, w.stijeh)

w.stijeh <- w.stijeh/max(w.stijeh)
w.stijeh <- (1-w.stijeh)+min(w.stijeh)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(stijeh$POINT_X, stijeh$POINT_Y, w.stijeh)
colPoints(stijeh$POINT_X, stijeh$POINT_Y, w.stijeh^2)

rm(d)


###BRT-----------------------------------------------------------------------------------
Releve_area <- 400
Releve_area <- as.data.frame(Releve_area)

###all species---------------------------------------------------------------------------
#current environment only
set.seed(1234)
mod1w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson", site.weights = w.stijeh,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1w.stijeh$gbm.call$best.trees
hist(mod1w.stijeh$residuals)
mod1w.stijeh$contributions

cor(stijeh$tot_rich, mod1w.stijeh$fitted)^2 #pseudo-R2
mean(mod1w.stijeh$residuals * mod1w.stijeh$residuals) #MSE
mod1w.stijeh$cv.statistics$deviance.mean #minimum cv deviance
mod1w.stijeh$cv.statistics$deviance.se #cv deviance se

pred1w.stijeh <- predict(env, mod1w.stijeh, const=Releve_area, n.trees=mod1w.stijeh$gbm.call$best.trees, type="response")
plot(pred1w.stijeh, main="Dark coniferous forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred1w.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1w.stijeh <- mask(pred1w.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred1w.stijeh <- round(pred1w.stijeh)
plot(pred1w.stijeh, main="Dark coniferous forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred1w.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1w.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred1w.stijeh.mask.tif", sep=""))
crs(pred1w.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod1w.stijeh$residuals)
plot(mor.cor)


#historical factors
set.seed(1234)
mod2w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson", site.weights = w.stijeh,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "tot_rich", family = "poisson", site.weights = w.stijeh,
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod2w.stijeh$gbm.call$best.trees
hist(mod2w.stijeh$residuals)
mod22w.stijeh$contributions

cor(stijeh$tot_rich, mod2w.stijeh$fitted)^2 #pseudo-R2
mean(mod2w.stijeh$residuals * mod2w.stijeh$residuals) #MSE
mod2w.stijeh$cv.statistics$deviance.mean #min cv deviance
mod2w.stijeh$cv.statistics$deviance.se #cv deviance se

pred2w.stijeh <- predict(env.hist, mod2w.stijeh, const=Releve_area, n.trees=mod2w.stijeh$gbm.call$best.trees, type="response")
plot(pred2w.stijeh, main="Dark coniferous forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred2w.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2w.stijeh <- mask(pred2w.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred2w.stijeh <- round(pred2w.stijeh)
plot(pred2w.stijeh, main="Dark coniferous forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred2w.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2w.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred2w.stijeh.mask.tif", sep=""))
crs(pred2w.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod2w.stijeh$residuals)
plot(mor.cor)


###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1w.stijeh$gbm.call$best.trees, mod2w.stijeh$gbm.call$best.trees, mod22w.stijeh$gbm.call$best.trees)))
y.lim <- c(min(c(mod1w.stijeh$cv.values - mod1w.stijeh$cv.loss.ses,
                 mod2w.stijeh$cv.values - mod2w.stijeh$cv.loss.ses,
                 mod22w.stijeh$cv.values - mod22w.stijeh$cv.loss.ses)),
           max(c(mod1w.stijeh$cv.values + mod1w.stijeh$cv.loss.ses,
                 mod2w.stijeh$cv.values + mod2w.stijeh$cv.loss.ses,
                 mod22w.stijeh$cv.values + mod22w.stijeh$cv.loss.ses)))

plot(mod1w.stijeh$trees.fitted, mod1w.stijeh$cv.values - mod1w.stijeh$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of dark coniferous forests\nd = 5, lr = 0.001")

polygon(c(mod1w.stijeh$trees.fitted, rev(mod1w.stijeh$trees.fitted)),
        c(mod1w.stijeh$cv.values - mod1w.stijeh$cv.loss.ses, rev(mod1w.stijeh$cv.values + mod1w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1w.stijeh$trees.fitted, mod1w.stijeh$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22w.stijeh$trees.fitted, rev(mod22w.stijeh$trees.fitted)),
        c(mod22w.stijeh$cv.values - mod22w.stijeh$cv.loss.ses, rev(mod22w.stijeh$cv.values + mod22w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22w.stijeh$trees.fitted, mod22w.stijeh$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2w.stijeh$trees.fitted, rev(mod2w.stijeh$trees.fitted)),
        c(mod2w.stijeh$cv.values - mod2w.stijeh$cv.loss.ses, rev(mod2w.stijeh$cv.values + mod2w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2w.stijeh$trees.fitted, mod2w.stijeh$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1w.stijeh$cv.values, mod2w.stijeh$cv.values, mod22w.stijeh$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")


###habitat specialists------------------------------------------------------------------
#current environment only
set.seed(1234)
mod3w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson", site.weights = w.stijeh,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3w.stijeh$gbm.call$best.trees
# hist(mod3w.stijeh$residuals)
mod3w.stijeh$contributions

cor(stijeh$dg_rich, mod3w.stijeh$fitted)^2 #pseudo-R2
mean(mod3w.stijeh$residuals * mod3w.stijeh$residuals) #MSE
mod3w.stijeh$cv.statistics$deviance.mean #min cv deviance
mod3w.stijeh$cv.statistics$deviance.se #cv deviance se

pred3w.stijeh <- predict(env, mod3w.stijeh, const=Releve_area, n.trees=mod3w.stijeh$gbm.call$best.trees, type="response")
plot(pred3w.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred3w.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3w.stijeh <- mask(pred3w.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred3w.stijeh <- round(pred3w.stijeh)
plot(pred3w.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred3w.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3w.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred3w.stijeh.mask.tif", sep=""))
crs(pred3w.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod3w.stijeh$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod4w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson", site.weights = w.stijeh,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42w.stijeh <-  gbm.step(data=stijeh, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "dg_rich", family = "poisson", site.weights = w.stijeh,
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod4w.stijeh$gbm.call$best.trees
hist(mod4w.stijeh$residuals)
mod4w.stijeh$contributions

cor(stijeh$dg_rich, mod4w.stijeh$fitted)^2 #pseudo-R2
mean(mod4w.stijeh$residuals * mod4w.stijeh$residuals) #MSE
mod4w.stijeh$cv.statistics$deviance.mean #min cv deviance
mod4w.stijeh$cv.statistics$deviance.se #cv deviance se

pred4w.stijeh <- predict(env.hist, mod4w.stijeh, const=Releve_area, n.trees=mod4w.stijeh$gbm.call$best.trees, type="response")
plot(pred4w.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred4w.stijeh.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4w.stijeh <- mask(pred4w.stijeh, ranges[["stijeh_1km"]], maskvalue=0)
pred4w.stijeh <- round(pred4w.stijeh)
plot(pred4w.stijeh, main="Dark coniferous forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4w.stijeh, filename = paste(getwd(), "/brt_rasters/stijeh/pred4w.stijeh.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4w.stijeh <- raster(paste(getwd(), "/brt_rasters/stijeh/pred4w.stijeh.mask.tif", sep=""))
crs(pred4w.stijeh) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(stijeh[, c("POINT_X", "POINT_Y")], mod4w.stijeh$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3w.stijeh$gbm.call$best.trees, mod4w.stijeh$gbm.call$best.trees, mod42w.stijeh$gbm.call$best.trees)))
y.lim <- c(min(c(mod3w.stijeh$cv.values - mod3w.stijeh$cv.loss.ses,
                 mod4w.stijeh$cv.values - mod4w.stijeh$cv.loss.ses,
                 mod42w.stijeh$cv.values - mod42w.stijeh$cv.loss.ses)),
           max(c(mod3w.stijeh$cv.values + mod3w.stijeh$cv.loss.ses,
                 mod4w.stijeh$cv.values + mod4w.stijeh$cv.loss.ses,
                 mod42w.stijeh$cv.values + mod42w.stijeh$cv.loss.ses)))

plot(mod3w.stijeh$trees.fitted, mod3w.stijeh$cv.values - mod3w.stijeh$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of dark coniferous forests\nd = 5, lr = 0.001")

polygon(c(mod3w.stijeh$trees.fitted, rev(mod3w.stijeh$trees.fitted)),
        c(mod3w.stijeh$cv.values - mod3w.stijeh$cv.loss.ses, rev(mod3w.stijeh$cv.values + mod3w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3w.stijeh$trees.fitted, mod3w.stijeh$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42w.stijeh$trees.fitted, rev(mod42w.stijeh$trees.fitted)),
        c(mod42w.stijeh$cv.values - mod42w.stijeh$cv.loss.ses, rev(mod42w.stijeh$cv.values + mod42w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42w.stijeh$trees.fitted, mod42w.stijeh$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4w.stijeh$trees.fitted, rev(mod4w.stijeh$trees.fitted)),
        c(mod4w.stijeh$cv.values - mod4w.stijeh$cv.loss.ses, rev(mod4w.stijeh$cv.values + mod4w.stijeh$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4w.stijeh$trees.fitted, mod4w.stijeh$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4w.stijeh$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3w.stijeh$cv.values, mod4w.stijeh$cv.values, mod42w.stijeh$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")


###LIGHT FORESTS-------------------------------------------------------------------------
#create weights

d <- as.matrix(dist(rbind(svetle[, c("POINT_X", "POINT_Y")], paleo.coords)))
dim(d)
d <- d[(nrow(svetle)+1):nrow(d),1:nrow(svetle)]
d[1:5,1:5]

w.svetle <- apply(d, 2, FUN=min)
hist(w.svetle)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(svetle$POINT_X, svetle$POINT_Y, w.svetle)

w.svetle <- w.svetle/max(w.svetle)
w.svetle <- (1-w.svetle)+min(w.svetle)
range(w.svetle)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(svetle$POINT_X, svetle$POINT_Y, w.svetle)
colPoints(svetle$POINT_X, svetle$POINT_Y, w.svetle^2)


###BRT-----------------------------------------------------------------------------------
Releve_area <- 400
Releve_area <- as.data.frame(Releve_area)

###all species---------------------------------------------------------------------------
#current environment only
set.seed(1234)
mod1w.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson", site.weights = w.svetle,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1w.svetle$gbm.call$best.trees
hist(mod1w.svetle$residuals)
mod1w.svetle$contributions

cor(svetle$tot_rich, mod1w.svetle$fitted)^2 #pseudo-R2
mean(mod1w.svetle$residuals * mod1w.svetle$residuals) #MSE
mod1w.svetle$cv.statistics$deviance.mean #min cv deviance
mod1w.svetle$cv.statistics$deviance.se #cv deviance se

pred1w.svetle <- predict(env, mod1w.svetle, const=Releve_area, n.trees=mod1w.svetle$gbm.call$best.trees, type="response")
plot(pred1w.svetle, main="Light forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred1w.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1w.svetle <- mask(pred1w.svetle, ranges[["list_1km"]], maskvalue=0)
pred1w.svetle <- round(pred1w.svetle)
plot(pred1w.svetle, main="Light forests (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred1w.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1w.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred1w.svetle.mask.tif", sep=""))
crs(pred1w.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod1w.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod2w.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "tot_rich", family = "poisson", site.weights = w.svetle,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22w.svetle <-  gbm.step(data=svetle, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "tot_rich", family = "poisson", site.weights = w.svetle,
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)


mod2w.svetle$gbm.call$best.trees
hist(mod2w.svetle$residuals)
mod2w.svetle$contributions

cor(svetle$tot_rich, mod2w.svetle$fitted)^2 #pseudo-R2
mean(mod2w.svetle$residuals * mod2w.svetle$residuals) #MSE
mod2w.svetle$cv.statistics$deviance.mean #min cv deviance
mod2w.svetle$cv.statistics$deviance.se #cv deviance se

pred2w.svetle <- predict(env.hist, mod2w.svetle, const=Releve_area, n.trees=mod2w.svetle$gbm.call$best.trees, type="response")
plot(pred2w.svetle, main="Light forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred2w.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2w.svetle <- mask(pred2w.svetle, ranges[["list_1km"]], maskvalue=0)
pred2w.svetle <- round(pred2w.svetle)
plot(pred2w.svetle, main="Light forests (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred2w.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2w.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred2w.svetle.mask.tif", sep=""))
crs(pred2w.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod2w.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)


###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1w.svetle$gbm.call$best.trees, mod2w.svetle$gbm.call$best.trees, mod22w.svetle$gbm.call$best.trees)))
y.lim <- c(min(c(mod1w.svetle$cv.values - mod1w.svetle$cv.loss.ses,
                 mod2w.svetle$cv.values - mod2w.svetle$cv.loss.ses,
                 mod22w.svetle$cv.values - mod22w.svetle$cv.loss.ses)),
           max(c(mod1w.svetle$cv.values + mod1w.svetle$cv.loss.ses,
                 mod2w.svetle$cv.values + mod2w.svetle$cv.loss.ses,
                 mod22w.svetle$cv.values + mod22w.svetle$cv.loss.ses)))

plot(mod1w.svetle$trees.fitted, mod1w.svetle$cv.values - mod1w.svetle$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of light forests\nd = 5, lr = 0.001")

polygon(c(mod1w.svetle$trees.fitted, rev(mod1w.svetle$trees.fitted)),
        c(mod1w.svetle$cv.values - mod1w.svetle$cv.loss.ses, rev(mod1w.svetle$cv.values + mod1w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1w.svetle$trees.fitted, mod1w.svetle$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22w.svetle$trees.fitted, rev(mod22w.svetle$trees.fitted)),
        c(mod22w.svetle$cv.values - mod22w.svetle$cv.loss.ses, rev(mod22w.svetle$cv.values + mod22w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22w.svetle$trees.fitted, mod22w.svetle$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2w.svetle$trees.fitted, rev(mod2w.svetle$trees.fitted)),
        c(mod2w.svetle$cv.values - mod2w.svetle$cv.loss.ses, rev(mod2w.svetle$cv.values + mod2w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2w.svetle$trees.fitted, mod2w.svetle$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1w.svetle$cv.values, mod2w.svetle$cv.values, mod22w.svetle$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")


###diagnostic species------------------------------------------------------------------
#current environment only
set.seed(1234)
mod3w.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson", site.weights = w.svetle,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3w.svetle$gbm.call$best.trees
hist(mod3w.svetle$residuals)
mod3w.svetle$contributions

cor(svetle$dg_rich, mod3w.svetle$fitted)^2 #pseudo-R2
mean(mod3w.svetle$residuals * mod3w.svetle$residuals) #MSE
mod3w.svetle$cv.statistics$deviance.mean #min cv deviance
mod3w.svetle$cv.statistics$deviance.se #cv deviance se

pred3w.svetle <- predict(env, mod3w.svetle, const=Releve_area, n.trees=mod3w.svetle$gbm.call$best.trees, type="response")
plot(pred3w.svetle, main="Light forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred3w.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3w.svetle <- mask(pred3w.svetle, ranges[["list_1km"]], maskvalue=0)
pred3w.svetle <- round(pred3w.svetle)
plot(pred3w.svetle, main="Light forests (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred3w.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3w.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred3w.svetle.mask.tif", sep=""))
crs(pred3w.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod3w.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)


#historical factors
set.seed(1234)
mod4w.svetle <-  gbm.step(data=svetle, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                         gbm.y = "dg_rich", family = "poisson", site.weights = w.svetle,
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42w.svetle <-  gbm.step(data=svetle, gbm.x = c(hist.vars[-6], "Releve_area"), 
                          gbm.y = "dg_rich", family = "poisson", site.weights = w.svetle,
                          tree.complexity = 5, learning.rate = 0.001, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod4w.svetle$gbm.call$best.trees
hist(mod4w.svetle$residuals)
mod4w.svetle$contributions

cor(svetle$dg_rich, mod4w.svetle$fitted)^2 #pseudo-R2
mean(mod4w.svetle$residuals * mod4w.svetle$residuals) #MSE
mod4w.svetle$cv.statistics$deviance.mean #min cv deviance
mod4w.svetle$cv.statistics$deviance.se #cv deviance se

pred4w.svetle <- predict(env.hist, mod4w.svetle, const=Releve_area, n.trees=mod4w.svetle$gbm.call$best.trees, type="response")
plot(pred4w.svetle, main="Light forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred4w.svetle.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4w.svetle <- mask(pred4w.svetle, ranges[["list_1km"]], maskvalue=0)
pred4w.svetle <- round(pred4w.svetle)
plot(pred4w.svetle, main="Light forests (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4w.svetle, filename = paste(getwd(), "/brt_rasters/svetle/pred4w.svetle.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4w.svetle <- raster(paste(getwd(), "/brt_rasters/svetle/pred4w.svetle.mask.tif", sep=""))
crs(pred4w.svetle) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(svetle[, c("POINT_X", "POINT_Y")], mod4w.svetle$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3w.svetle$gbm.call$best.trees, mod4w.svetle$gbm.call$best.trees, mod42w.svetle$gbm.call$best.trees)))
y.lim <- c(min(c(mod3w.svetle$cv.values - mod3w.svetle$cv.loss.ses,
                 mod4w.svetle$cv.values - mod4w.svetle$cv.loss.ses,
                 modh.svetle$cv.values - mod42w.svetle$cv.loss.ses)),
           max(c(mod3w.svetle$cv.values + mod3w.svetle$cv.loss.ses,
                 mod4w.svetle$cv.values + mod4w.svetle$cv.loss.ses,
                 mod42w.svetle$cv.values + mod42w.svetle$cv.loss.ses)))

plot(mod3w.svetle$trees.fitted, mod3w.svetle$cv.values - mod3w.svetle$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of light forests\nd = 5, lr = 0.001")

polygon(c(mod3w.svetle$trees.fitted, rev(mod3w.svetle$trees.fitted)),
        c(mod3w.svetle$cv.values - mod3w.svetle$cv.loss.ses, rev(mod3w.svetle$cv.values + mod3w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3w.svetle$trees.fitted, mod3w.svetle$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42w.svetle$trees.fitted, rev(mod42w.svetle$trees.fitted)),
        c(mod42w.svetle$cv.values - mod42w.svetle$cv.loss.ses, rev(mod42w.svetle$cv.values + mod42w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42w.svetle$trees.fitted, mod42w.svetle$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4w.svetle$trees.fitted, rev(mod4w.svetle$trees.fitted)),
        c(mod4w.svetle$cv.values - mod4w.svetle$cv.loss.ses, rev(mod4w.svetle$cv.values + mod4w.svetle$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4w.svetle$trees.fitted, mod4w.svetle$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4w.svetle$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3w.svetle$cv.values, mod4w.svetle$cv.values, mod42w.svetle$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###SEMI-DRY AND STEPPE GRASSLANDS------------------------------------------------------
#create weights

d <- as.matrix(dist(rbind(travniky[, c("POINT_X", "POINT_Y")], paleo.coords)))
dim(d)
d <- d[(nrow(travniky)+1):nrow(d),1:nrow(travniky)]
d[1:5,1:5]

w.travniky <- apply(d, 2, FUN=min)
hist(w.travniky)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(travniky$POINT_X, travniky$POINT_Y, w.travniky)

w.travniky <- w.travniky/max(w.travniky)
w.travniky <- (1-w.travniky)+min(w.travniky)
range(w.travniky)

plot(paleo.coords, pch=16, col="red", cex=1.5)
colPoints(travniky$POINT_X, travniky$POINT_Y, w.travniky)
colPoints(travniky$POINT_X, travniky$POINT_Y, w.travniky^2)


###BRT---------------------------------------------------------------------------------
Releve_area <- 25
Releve_area <- as.data.frame(Releve_area)

#all species---------------------------------------------------------------------------
#current environment only

set.seed(1234)
mod1w.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], "Releve_area"), 
                          gbm.y = "Pocet_druhu_celkem", family = "poisson", site.weights = w.travniky,
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod1w.travniky$gbm.call$best.trees
hist(mod1w.travniky$residuals)
mod1w.travniky$contributions

cor(travniky$Pocet_druhu_celkem, mod1w.travniky$fitted)^2
mean(mod1w.travniky$residuals * mod1w.travniky$residuals) #MSE
mod1w.travniky$cv.statistics$deviance.mean #min cv deviance
mod1w.travniky$cv.statistics$deviance.se #cv deviance se

pred1w.travniky <- predict(env, mod1w.travniky, const=Releve_area, n.trees=mod1w.travniky$gbm.call$best.trees, type="response")
plot(pred1w.travniky, main="Travniky (No. species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred1w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred1w.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred1w.travniky <- mask(pred1w.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred1w.travniky <- round(pred1w.travniky)
plot(pred1w.travniky, main="Travniky (No. species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred1w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred1w.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred1w.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred1w.travniky.mask.tif", sep=""))
crs(pred1w.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod1w.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod2w.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                          gbm.y = "Pocet_druhu_celkem", family = "poisson", site.weights = w.travniky,
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod22w.travniky <- gbm.step(data=travniky, gbm.x = c(hist.vars[-6], "Releve_area"), 
                           gbm.y = "Pocet_druhu_celkem", family = "poisson", site.weights = w.travniky,
                           tree.complexity = 5, learning.rate = 0.003, 
                           bag.fraction = 0.5, step.size=100, max.trees = 30000)


mod2w.travniky$gbm.call$best.trees
hist(mod2w.travniky$residuals)
mod2w.travniky$contributions

cor(travniky$Pocet_druhu_celkem, mod2w.travniky$fitted)^2 #pseudo-R2
mean(mod2w.travniky$residuals * mod2w.travniky$residuals) #MSE
mod2w.travniky$cv.statistics$deviance.mean #min cv deviance
mod2w.travniky$cv.statistics$deviance.se #cv deviance se 

pred2w.travniky <- predict(env.hist, mod2w.travniky, const=Releve_area, n.trees=mod2w.travniky$gbm.call$best.trees, type="response")
plot(pred2w.travniky, main="Travniky (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred2w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred2w.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred2w.travniky <- mask(pred2w.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred2w.travniky <- round(pred2w.travniky)
plot(pred2w.travniky, main="Travniky (No. species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred2w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred2w.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred2w.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred2w.travniky.mask.tif", sep=""))
crs(pred2w.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod2w.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod1w.travniky$gbm.call$best.trees, mod2w.travniky$gbm.call$best.trees, mod22w.travniky$gbm.call$best.trees)))
y.lim <- c(min(c(mod1w.travniky$cv.values - mod1w.travniky$cv.loss.ses,
                 mod2w.travniky$cv.values - mod2w.travniky$cv.loss.ses,
                 mod22w.travniky$cv.values - mod22w.travniky$cv.loss.ses)),
           max(c(mod1w.travniky$cv.values + mod1w.travniky$cv.loss.ses,
                 mod2w.travniky$cv.values + mod2w.travniky$cv.loss.ses,
                 mod22w.travniky$cv.values + mod22w.travniky$cv.loss.ses)))

plot(mod1w.travniky$trees.fitted, mod1w.travniky$cv.values - mod1w.travniky$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Total species richness of semi-dry and steppe grasslands\nd = 5, lr = 0.003")

polygon(c(mod1w.travniky$trees.fitted, rev(mod1w.travniky$trees.fitted)),
        c(mod1w.travniky$cv.values - mod1w.travniky$cv.loss.ses, rev(mod1w.travniky$cv.values + mod1w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod1w.travniky$trees.fitted, mod1w.travniky$cv.values, lwd=2, col=my.col[3])

polygon(c(mod22w.travniky$trees.fitted, rev(mod22w.travniky$trees.fitted)),
        c(mod22w.travniky$cv.values - mod22w.travniky$cv.loss.ses, rev(mod22w.travniky$cv.values + mod22w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod22w.travniky$trees.fitted, mod22w.travniky$cv.values, lwd=2, col=my.col[2])

polygon(c(mod2w.travniky$trees.fitted, rev(mod2w.travniky$trees.fitted)),
        c(mod2w.travniky$cv.values - mod2w.travniky$cv.loss.ses, rev(mod2w.travniky$cv.values + mod2w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod2w.travniky$trees.fitted, mod2w.travniky$cv.values, lwd=2, col=my.col[1])

lines(rep(mod1w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod22w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod2w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod1w.travniky$cv.values, mod2w.travniky$cv.values, mod22w.travniky$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

####diagnostic species--------------------------------------------------
#current environment only
set.seed(1234)
mod3w.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], "Releve_area"), 
                          gbm.y = "Pocet_vybranych_druhu", family = "poisson", site.weights = w.travniky,
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

mod3w.travniky$gbm.call$best.trees
hist(mod3w.travniky$residuals)
mod3w.travniky$contributions

cor(travniky$Pocet_vybranych_druhu, mod3w.travniky$fitted)^2 #pseudo-R2
mean(mod3w.travniky$residuals * mod3w.travniky$residuals) #MSE
mod3w.travniky$cv.statistics$deviance.mean #min cv deviance
mod3w.travniky$cv.statistics$deviance.se #cv deviance se

pred3w.travniky <- predict(env, mod3w.travniky, const=Releve_area, n.trees=mod3w.travniky$gbm.call$best.trees, type="response")
plot(pred3w.travniky, main="Travniky (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred3w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred3w.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred3w.travniky <- mask(pred3w.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred3w.travniky <- round(pred3w.travniky)
plot(pred3w.travniky, main="Travniky (No. diagnostic species ~ E)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred3w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred3w.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred3w.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred3w.travniky.mask.tif", sep=""))
crs(pred3w.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod3w.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

#historical factors
set.seed(1234)
mod4w.travniky <- gbm.step(data=travniky, gbm.x = c(env.vars[-1], hist.vars[-6], "Releve_area"), 
                          gbm.y = "Pocet_vybranych_druhu", family = "poisson", site.weights = w.travniky,
                          tree.complexity = 5, learning.rate = 0.003, 
                          bag.fraction = 0.5, step.size=100, max.trees = 30000)

#only history
set.seed(1234)
mod42w.travniky <- gbm.step(data=travniky, gbm.x = c(hist.vars[-6], "Releve_area"), 
                           gbm.y = "Pocet_vybranych_druhu", family = "poisson", site.weights = w.travniky,
                           tree.complexity = 5, learning.rate = 0.003, 
                           bag.fraction = 0.5, step.size=100, max.trees = 30000)


mod4w.travniky$gbm.call$best.trees
hist(mod4w.travniky$residuals)
mod4w.travniky$contributions

cor(travniky$Pocet_vybranych_druhu, mod4w.travniky$fitted)^2 #pseudo-R2
mean(mod4w.travniky$residuals * mod4w.travniky$residuals) #MSE
mod4w.travniky$cv.statistics$deviance.mean #min cv deviance
mod4w.travniky$cv.statistics$deviance.se #cv deviance se

pred4w.travniky <- predict(env.hist, mod4w.travniky, const=Releve_area, n.trees=mod4w.travniky$gbm.call$best.trees, type="response")
plot(pred4w.travniky, main="Travniky (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T)
writeRaster(pred4w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred4w.travniky.tif", sep=""), format="GTiff", overwrite=TRUE)

pred4w.travniky <- mask(pred4w.travniky, ranges[["travniky_1km"]], maskvalue=0)
pred4w.travniky <- round(pred4w.travniky)
plot(pred4w.travniky, main="Travniky (No. diagnostic species ~ E + H)", col=my.ramp(255))
plot(country, add=T, lwd=1)
writeRaster(pred4w.travniky, filename = paste(getwd(), "/brt_rasters/travniky/pred4w.travniky.mask.tif", sep=""), format="GTiff", overwrite=TRUE)
pred4w.travniky <- raster(paste(getwd(), "/brt_rasters/travniky/pred4w.travniky.mask.tif", sep=""))
crs(pred4w.travniky) <- crs(Alt)

#Spatial autocorrelation in residuals
mor.cor <- correlog(travniky[, c("POINT_X", "POINT_Y")], mod4w.travniky$residuals)
plot(mor.cor)
round(mor.cor[,3],3)

###Cross-validation progress
my.col <- viridis(3)

x.lim <- c(0, max(c(mod3w.travniky$gbm.call$best.trees, mod4w.travniky$gbm.call$best.trees, mod42w.travniky$gbm.call$best.trees)))
y.lim <- c(min(c(mod3w.travniky$cv.values - mod3w.travniky$cv.loss.ses,
                 mod4w.travniky$cv.values - mod4w.travniky$cv.loss.ses,
                 mod42w.travniky$cv.values - mod42w.travniky$cv.loss.ses)),
           max(c(mod3w.travniky$cv.values + mod3w.travniky$cv.loss.ses,
                 mod4w.travniky$cv.values + mod4w.travniky$cv.loss.ses,
                 mod42w.travniky$cv.values + mod42w.travniky$cv.loss.ses)))

plot(mod3w.travniky$trees.fitted, mod3w.travniky$cv.values - mod3w.travniky$cv.loss.ses, pch=NA,
     ylim= y.lim, xlim= x.lim,
     las=1, xlab="No. trees", ylab="Holdout deviance", main="Specialist species richness of semi-dry and steppe grasslands\nd = 5, lr = 0.003")

polygon(c(mod3w.travniky$trees.fitted, rev(mod3w.travniky$trees.fitted)),
        c(mod3w.travniky$cv.values - mod3w.travniky$cv.loss.ses, rev(mod3w.travniky$cv.values + mod3w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[3])[1,], col2rgb(my.col[3])[2,], col2rgb(my.col[3])[3,], 80, maxColorValue = 255))
lines(mod3w.travniky$trees.fitted, mod3w.travniky$cv.values, lwd=2, col=my.col[3])

polygon(c(mod42w.travniky$trees.fitted, rev(mod42w.travniky$trees.fitted)),
        c(mod42w.travniky$cv.values - mod42w.travniky$cv.loss.ses, rev(mod42w.travniky$cv.values + mod42w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[2])[1,], col2rgb(my.col[2])[2,], col2rgb(my.col[2])[3,], 80, maxColorValue = 255))
lines(mod42w.travniky$trees.fitted, mod42w.travniky$cv.values, lwd=2, col=my.col[2])

polygon(c(mod4w.travniky$trees.fitted, rev(mod4w.travniky$trees.fitted)),
        c(mod4w.travniky$cv.values - mod4w.travniky$cv.loss.ses, rev(mod4w.travniky$cv.values + mod4w.travniky$cv.loss.ses)),
        border = F, col=rgb(col2rgb(my.col[1])[1,], col2rgb(my.col[1])[2,], col2rgb(my.col[1])[3,], 80, maxColorValue = 255))
lines(mod4w.travniky$trees.fitted, mod4w.travniky$cv.values, lwd=2, col=my.col[1])

lines(rep(mod3w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[3])
lines(rep(mod42w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[2])
lines(rep(mod4w.travniky$gbm.call$best.trees, 2), y.lim + c(-1, 1), col=my.col[1])
lines(c(-800, x.lim[2] + 800), rep(min(c(mod3w.travniky$cv.values, mod4w.travniky$cv.values, mod42w.travniky$cv.values)),2), col="red")
legend("top", legend=c("~ E", "~ H", "~ E + H"), lty = 1, lwd = 3,
       col = my.col[c(3,2,1)], horiz=TRUE, bty = "o", bg = "white")

###extract variable contributions from the models--------------------------------
##all species + current factors
contrib1 <- list(mod1w.stijeh$contributions[c(env.vars[-1], "Releve_area"),],
                 mod1w.svetle$contributions[c(env.vars[-1], "Releve_area"),],
                 mod1w.travniky$contributions[c(env.vars[-1], "Releve_area"),]) 
contrib1 <- do.call(cbind.data.frame, contrib1)[,c(2,4,6)]
colnames(contrib1) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib1, 1)

#all species + current & historical factors
contrib2 <- list(mod2w.stijeh$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod2w.svetle$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod2w.travniky$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),]) 
contrib2 <- do.call(cbind.data.frame, contrib2)[,c(2,4,6)]
colnames(contrib2) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib2, 1)

##diagnostic species + current factors
contrib3 <- list(mod3w.stijeh$contributions[c(env.vars[-1], "Releve_area"),],
                 mod3w.svetle$contributions[c(env.vars[-1], "Releve_area"),],
                 mod3w.travniky$contributions[c(env.vars[-1], "Releve_area"),]) 
contrib3 <- do.call(cbind.data.frame, contrib3)[,c(2,4,6)]
colnames(contrib3) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib3, 1)

#diagnostic species + current & historical factors
contrib4 <- list(mod4w.stijeh$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod4w.svetle$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),],
                 mod4w.travniky$contributions[c(env.vars[-1], hist.vars[-6], "Releve_area"),]) 
contrib4 <- do.call(cbind.data.frame, contrib4)[,c(2,4,6)]
colnames(contrib4) <- c("Dark coniferous forests", "Light forests", "Semi-dry and steppe grasslands")
round(contrib4, 1)


###PARTIAL DEPENDENCE PLOTS---------------------------------------------------------------------

###Dark coniferous forests-----------------------------------------------------
mod1w.stijeh$contributions

pdpw.stijeh <- list()
pdpw.stijeh$mod1w <- list()
pdpw.stijeh$mod2w <- list()
pdpw.stijeh$mod3w <- list()
pdpw.stijeh$mod4w <- list()

brt.models <- list(mod1w.stijeh, mod2w.stijeh, mod3w.stijeh, mod4w.stijeh)

for(q in c(1:4))#1:4
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdpw.stijeh[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")

  }
  names(pdpw.stijeh[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)


par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdpw.stijeh$mod4w[[i]], type="l", lwd=1.3, main=names(pdpw.stijeh$mod4w)[i],
       ylab="No. species", col="blue")
  points(pdpw.stijeh$mod4w[[i]], pch=16)
}

plot(pdpw.stijeh$mod4w[[i]], type="p", lwd=1.3, main=names(pdpw.stijeh$mod4w)[i],
     ylab="No. species", col="gray90")



###Light forests-------------------------------------------------------------
mod1w.svetle$contributions

pdpw.svetle <- list()
pdpw.svetle$mod1w <- list()
pdpw.svetle$mod2w <- list()
pdpw.svetle$mod3w <- list()
pdpw.svetle$mod4w <- list()

brt.models <- list(mod1w.svetle, mod2w.svetle, mod3w.svetle, mod4w.svetle)

for(q in c(1:4))
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdpw.svetle[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")
    
  }
  names(pdpw.svetle[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)


par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdpw.svetle$mod4w[[i]], type="l", lwd=1.3, main=names(pdpw.svetle$mod4w)[i],
       ylab="No. species", col="blue")
  points(pdpw.svetle$mod4w[[i]], pch=16)
}

###Semi-dry and steppe grasslands-------------------------------------------------------------
mod1w.travniky$contributions

pdpw.travniky <- list()
pdpw.travniky$mod1w <- list()
pdpw.travniky$mod2w <- list()
pdpw.travniky$mod3w <- list()
pdpw.travniky$mod4w <- list()

brt.models <- list(mod1w.travniky, mod2w.travniky, mod3w.travniky, mod4w.travniky)

for(q in c(1:4))
{
  for(i in 1:length(brt.models[[q]]$contributions$var))
  {
    pdpw.travniky[[q]][[i]] <- plotq.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var)[i], 
                                      continuous.resolution = 51, return.grid = TRUE, type = "response")
    
  }
  names(pdpw.travniky[[q]]) <- as.character(brt.models[[q]]$contributions$var)#[1:6]
}
rm(brt.models)

par(mfrow=c(2,3), mar=c(5,5,3,1))
for(i in 1:6)
{
  plot(pdpw.travniky$mod4w[[i]], type="l", lwd=1.3, main=names(pdpw.travniky$mod4w)[i],
       ylab="No. species", col="blue")
  points(pdpw.travniky$mod4w[[i]], pch=16)
}



