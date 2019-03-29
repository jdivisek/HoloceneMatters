#########################################################################
#			PARTIAL DEPENDENCE PLOTS			#
#########################################################################

library(dismo)
library(gbm)

##All species-------------------------------------------------------------------------
windows()

tiff(width = 17.5, height = 26, units="cm", res=400, compression = "lzw")
layout(matrix(1:15, ncol=3, nrow=5, byrow = F))

par(mar=c(5,2.5,1,0.5), oma=c(0,2,2,0.5))

plot.data <- list(pdp.svetle$mod2, pdp.travniky$mod2, pdp.stijeh$mod2)
names(plot.data) <- c("Light forests", "Semi-dry and steppe grasslands", "Dark coniferous forests")

plot.data[[1]]$bio1[ ,1] <- plot.data[[1]]$bio1[ ,1]/10
plot.data[[1]]$bio4[ ,1] <- plot.data[[1]]$bio4[ ,1]/1000
plot.data[[1]]$pH[ ,1] <- plot.data[[1]]$pH[ ,1]/10
plot.data[[2]]$bio1[ ,1] <- plot.data[[2]]$bio1[ ,1]/10
plot.data[[2]]$bio4[ ,1] <- plot.data[[2]]$bio4[ ,1]/1000
plot.data[[2]]$pH[ ,1] <- plot.data[[2]]$pH[ ,1]/10
plot.data[[3]]$bio1[ ,1] <- plot.data[[3]]$bio1[ ,1]/10
plot.data[[3]]$bio4[ ,1] <- plot.data[[3]]$bio4[ ,1]/1000
plot.data[[3]]$pH[ ,1] <- plot.data[[3]]$pH[ ,1]/10

brt.models <- list(mod2.svetle, mod2.travniky, mod2.stijeh)

var.names <- c("Terrain ruggedness (index)", "Heat load (index)", "Topographic wetness (index)", "Soil pH", "Mean annual temperature (°C)",
               "Temperature seasonality (Std. dev.)", "Annual precipitation (mm/year)", "Precipitation seasonality\n(Coef. of variation)", "Area of forests (%)",
               "Area of grasslands (%)", "Area of arable land (%)", "Build-up area (%)", "Landscape diversity (index)", "Temperateness in Late Glacial\n(PCoA 1)",
               "Temperateness at Holocene onset\n(PCoA 1)", "Landscape openness in Pre-Neolithic\n(PCoA 1)", "Landscape openness in Late Neolithic\n(PCoA 1)", "Representation of taiga in Late Prehistory\n(PCoA 1)",
               "Distance from Neolithic settlement (km)", NA)
names(var.names) <- c(env.vars[-1], hist.vars, "Releve_area")
head(var.names)


for(q in 1:3)
{
  # plot.means <- lapply(plot.data[[q]], FUN=function(x){aggregate(. ~ cut, data=x, FUN=mean)})
  ylim <- c(min(unlist(lapply(plot.data[[q]], FUN=function(x){min(x$y)}))),
              max(unlist(lapply(plot.data[[q]], FUN=function(x){max(x$y)}))))
  
  for(i in 16:20)
  {
    if(i <= length(plot.data[[q]]))
    {
      plot(plot.data[[q]][[i]][ ,1], plot.data[[q]][[i]]$y, type="b", pch=16, col="gray75", ylim=ylim, 
           las=1, xlab=ifelse(is.na(var.names[names(plot.data[[q]])[i]]), expression(Plot~size~(m^2)), var.names[names(plot.data[[q]])[i]]), 
           ylab="", cex.lab =0.9, cex.axis=0.9, cex=0.6)
      
      lo.data <- gbm::plot.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var[i]), continuous.resolution = 100, return.grid = TRUE, type = "response")
      
      if(as.character(brt.models[[q]]$contributions$var[i]) %in% c("pH", "bio1")) {lo.data[,1] <- lo.data[,1]/10}
      if(as.character(brt.models[[q]]$contributions$var[i]) == "bio4") {lo.data[,1] <- lo.data[,1]/1000}
      
      temp.lo <- loess(lo.data$y ~ lo.data[,1], span = 0.5)
      lines(lo.data[,1], fitted(temp.lo), lty = 2, col = "red")
      
      text(min(plot.data[[q]][[i]][ ,1]), max(ylim)*0.96, 
           paste(round(brt.models[[q]]$contributions$rel.inf[i],1), "%", sep = ""), font=2, cex=1, pos=4)
      
      if(i %in% c(1,6,11,16))
      {
        mtext(names(plot.data)[q], side=3, line = 0.5, outer = FALSE, font=2, cex=0.8)
      }
      
      if(q==1)
      {
        mtext("No. species", side=2, line = 3, outer = FALSE, font=1, cex=0.7)
        
      }
    }
    else
    {
      plot(1, type='n', axes=F, xlab = "")
    }
  }
  
}

rm(brt.models, plot.data, temp.lo, lo.data)

dev.off()

##Diagnostic species-------------------------------------------------------------------------
windows()

tiff(width = 17.5, height = 26, units="cm", res=400, compression = "lzw")
layout(matrix(1:15, ncol=3, nrow=5, byrow = F))

par(mar=c(5,2.5,1,0.5), oma=c(0,2,2,0.5))

plot.data <- list(pdp.svetle$mod4, pdp.travniky$mod4, pdp.stijeh$mod4)
names(plot.data) <- c("Light forests", "Semi-dry and steppe grasslands", "Dark coniferous forests")

plot.data[[1]]$bio1[ ,1] <- plot.data[[1]]$bio1[ ,1]/10
plot.data[[1]]$bio4[ ,1] <- plot.data[[1]]$bio4[ ,1]/1000
plot.data[[1]]$pH[ ,1] <- plot.data[[1]]$pH[ ,1]/10
plot.data[[2]]$bio1[ ,1] <- plot.data[[2]]$bio1[ ,1]/10
plot.data[[2]]$bio4[ ,1] <- plot.data[[2]]$bio4[ ,1]/1000
plot.data[[2]]$pH[ ,1] <- plot.data[[2]]$pH[ ,1]/10
plot.data[[3]]$bio1[ ,1] <- plot.data[[3]]$bio1[ ,1]/10
plot.data[[3]]$bio4[ ,1] <- plot.data[[3]]$bio4[ ,1]/1000
plot.data[[3]]$pH[ ,1] <- plot.data[[3]]$pH[ ,1]/10

brt.models <- list(mod4.svetle, mod4.travniky, mod4.stijeh)

var.names <- c("Terrain ruggedness (index)", "Heat load (index)", "Topographic wetness (index)", "Soil pH", "Mean annual temperature (°C)",
               "Temperature seasonality (Std. dev.)", "Annual precipitation (mm/year)", "Precipitation seasonality\n(Coef. of variation)", "Area of forests (%)",
               "Area of grasslands (%)", "Area of arable land (%)", "Build-up area (%)", "Landscape diversity (index)", "Temperateness in Late Glacial\n(PCoA 1)",
               "Temperateness at Holocene onset\n(PCoA 1)", "Landscape openness in Pre-Neolithic\n(PCoA 1)", "Landscape openness in Late Neolithic\n(PCoA 1)", "Representation of taiga in Late Prehistory\n(PCoA 1)",
               "Distance from Neolithic settlement (km)", NA)
names(var.names) <- c(env.vars[-1], hist.vars, "Releve_area")
head(var.names)


for(q in 1:3)
{
  # plot.means <- lapply(plot.data[[q]], FUN=function(x){aggregate(. ~ cut, data=x, FUN=mean)})
  ylim <- c(min(unlist(lapply(plot.data[[q]], FUN=function(x){min(x$y)}))),
            max(unlist(lapply(plot.data[[q]], FUN=function(x){max(x$y)}))))
  
  for(i in 16:20)
  {
    if(i <= length(plot.data[[q]]))
    {
      plot(plot.data[[q]][[i]][ ,1], plot.data[[q]][[i]]$y, type="b", pch=16, col="gray75", ylim=ylim, 
           las=1, xlab=ifelse(is.na(var.names[names(plot.data[[q]])[i]]), expression(Plot~size~(m^2)), var.names[names(plot.data[[q]])[i]]), 
           ylab="", cex.lab =0.9, cex.axis=0.9, cex=0.6)
      
      lo.data <- gbm::plot.gbm(brt.models[[q]], as.character(brt.models[[q]]$contributions$var[i]), continuous.resolution = 100, return.grid = TRUE, type = "response")
      
      if(as.character(brt.models[[q]]$contributions$var[i]) %in% c("pH", "bio1")) {lo.data[,1] <- lo.data[,1]/10}
      if(as.character(brt.models[[q]]$contributions$var[i]) == "bio4") {lo.data[,1] <- lo.data[,1]/1000}
      
      temp.lo <- loess(lo.data$y ~ lo.data[,1], span = 0.5)
      lines(lo.data[,1], fitted(temp.lo), lty = 2, col = "red")
      
      text(min(plot.data[[q]][[i]][ ,1]), max(ylim)*0.96, 
           paste(round(brt.models[[q]]$contributions$rel.inf[i],1), "%", sep = ""), font=2, cex=1, pos=4)
      
      if(i %in% c(1,6,11,16))
      {
        mtext(names(plot.data)[q], side=3, line = 0.5, outer = FALSE, font=2, cex=0.8)
      }
      
      if(q==1)
      {
        mtext("No. specialists", side=2, line = 3, outer = FALSE, font=1, cex=0.7)
        
      }
    }
    else
    {
      plot(1, type='n', axes=F, xlab = "")
    }
  }
  
}

rm(brt.models, plot.data, temp.lo, lo.data)

dev.off()