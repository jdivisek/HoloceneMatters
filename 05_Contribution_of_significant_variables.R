#################################################################################
#		PLOT CONTRIBUTIONS OF SIGNIFICANT VARIABLES			#
#################################################################################

###IMPORT P-VALUES-----------------------------------------------------------------------------------------------
paths <- list.files(path=paste(getwd(), "/BRT - var. contribution tests", sep=""), pattern='txt', full.names=TRUE ) 

###total species richness----------------------------------------------------------------------------------------
mod2.p <- list()

for(i in 1:3)
{
  mod2.p[[i]] <- read.delim(paths[4:6][i], header=T, row.names=1)
}
names(mod2.p) <- c("p2.stijeh", "p2.svetle", "p2.travniky")
mod2.p


set.seed(1234)
mod2sig.stijeh <-  gbm.step(data=stijeh, gbm.x = rownames(mod2.p$p2.stijeh)[mod2.p$p2.stijeh$p2.stijeh < 0.05], 
                         gbm.y = "tot_rich", family = "poisson",
                         tree.complexity = 5, learning.rate = 0.001, 
                         bag.fraction = 0.5, step.size=100, max.trees = 30000)

set.seed(1234)
mod2sig.svetle <-  gbm.step(data=svetle, gbm.x = rownames(mod2.p$p2.svetle)[mod2.p$p2.svetle$p2.svetle < 0.05], 
                            gbm.y = "tot_rich", family = "poisson",
                            tree.complexity = 5, learning.rate = 0.001, 
                            bag.fraction = 0.5, step.size=100, max.trees = 30000)

set.seed(1234)
mod2sig.travniky <-  gbm.step(data=travniky, gbm.x = rownames(mod2.p$p2.travniky)[mod2.p$p2.travniky$p2.travniky < 0.05], 
                            gbm.y = "Pocet_druhu_celkem", family = "poisson",
                            tree.complexity = 5, learning.rate = 0.003, 
                            bag.fraction = 0.5, step.size=100, max.trees = 30000)

plot.data <- matrix(data=NA, ncol=3, nrow=3, byrow = F)

plot.data[1,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[1:5])
plot.data[2,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[6:9])
plot.data[3,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[10])

plot.data[1,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[1:7])
plot.data[2,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[8:11])
plot.data[3,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[12])

plot.data[1,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[1:5])
plot.data[2,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[6:10])
plot.data[3,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[11])
plot.data

head(plot.data)
barplot(plot.data[, c(3,2,1)], beside=F, horiz=T, col=c("gold1", "dodgerblue3", "gray40"),
        border=F, main="Total species richness", names.arg = c("Dark coniferous\nforests", "Semi-dry and steppe\ngrasslands", "Light forests"),
        las=1)


###diagnostic species richness--------------------------------------------------------------
mod4.p <- list()

for(i in 1:3)
{
  mod4.p[[i]] <- read.delim(paths[10:12][i], header=T, row.names=1)
}
names(mod4.p) <- c("p4.stijeh", "p4.svetle", "p4.travniky")
mod4.p

set.seed(1234)
mod4sig.stijeh <-  gbm.step(data=stijeh, gbm.x = rownames(mod4.p$p4.stijeh)[mod4.p$p4.stijeh$p4.stijeh < 0.05], 
                            gbm.y = "dg_rich", family = "poisson",
                            tree.complexity = 5, learning.rate = 0.001, 
                            bag.fraction = 0.5, step.size=100, max.trees = 30000)

set.seed(1234)
mod4sig.svetle <-  gbm.step(data=svetle, gbm.x = rownames(mod4.p$p4.svetle)[mod4.p$p4.svetle$p4.svetle < 0.05], 
                            gbm.y = "dg_rich", family = "poisson",
                            tree.complexity = 5, learning.rate = 0.001, 
                            bag.fraction = 0.5, step.size=100, max.trees = 30000)

set.seed(1234)
mod4sig.travniky <-  gbm.step(data=travniky, gbm.x = rownames(mod4.p$p4.travniky)[mod4.p$p4.travniky$p4.travniky < 0.05], 
                              gbm.y = "Pocet_vybranych_druhu", family = "poisson",
                              tree.complexity = 5, learning.rate = 0.003, 
                              bag.fraction = 0.5, step.size=100, max.trees = 30000)

plot.data <- matrix(data=NA, ncol=3, nrow=3, byrow = F)

plot.data[1,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[1:5])
plot.data[2,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[6:7])
plot.data[3,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[8])

plot.data[1,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[1:2])
plot.data[2,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[3:5])
plot.data[3,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[6])

plot.data[1,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[1:2])
plot.data[2,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[3:5])
plot.data[3,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[6])
plot.data

head(plot.data)
barplot(plot.data[, c(3,2,1)], beside=F, horiz=T, col=c("gold1", "dodgerblue3", "gray40"),
        border=F, main="Specialist species richness", names.arg = c("Dark coniferous\nforests", "Semi-dry and steppe\ngrasslands", "Light forests"),
        las=1)

###BARPLOTS-------------------------------------------------------------------------------------

par(mfrow=c(1,2), mar=c(5,0.5,2,1), oma=c(0,9,0,6), xpd=TRUE)

plot.data <- matrix(data=NA, ncol=3, nrow=3, byrow = F)

plot.data[1,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[1:5])
plot.data[2,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[6:9])
plot.data[3,1] <- sum(mod2sig.svetle$contributions[mod2sig.svetle$var.names,]$rel.inf[10])

plot.data[1,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[1:7])
plot.data[2,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[8:11])
plot.data[3,2] <- sum(mod2sig.travniky$contributions[mod2sig.travniky$var.names,]$rel.inf[12])

plot.data[1,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[1:5])
plot.data[2,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[6:10])
plot.data[3,3] <- sum(mod2sig.stijeh$contributions[mod2sig.stijeh$var.names,]$rel.inf[11])
plot.data

head(plot.data)
barplot(plot.data[, c(3,2,1)], beside=F, horiz=T, col=c("gray70", "gray40", "black"),
        border=F, main="Total species richness", names.arg = c("Dark coniferous\nforests", "Semi-dry and steppe\ngrasslands", "Light forests"),
        las=1)

plot.data <- matrix(data=NA, ncol=3, nrow=3, byrow = F)

plot.data[1,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[1:5])
plot.data[2,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[6:7])
plot.data[3,1] <- sum(mod4sig.svetle$contributions[mod4sig.svetle$var.names,]$rel.inf[8])

plot.data[1,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[1:2])
plot.data[2,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[3:5])
plot.data[3,2] <- sum(mod4sig.travniky$contributions[mod4sig.travniky$var.names,]$rel.inf[6])

plot.data[1,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[1:2])
plot.data[2,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[3:5])
plot.data[3,3] <- sum(mod4sig.stijeh$contributions[mod4sig.stijeh$var.names,]$rel.inf[6])
plot.data

head(plot.data)
barplot(plot.data[, c(3,2,1)], beside=F, horiz=T, col=c("gray70", "gray40", "black"),
        border=F, main="Specialist species richness", las=1)
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(3, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend("right", legend=c("Environment", "History", "Plot size"),
       pch=15, horiz = FALSE, bty = "n", pt.cex = 2, col=c("gray70", "gray40", "black"), inset=c(0,0), xpd=TRUE)
text(0.06, -1.25, "Sum of variable contributions (%)")
