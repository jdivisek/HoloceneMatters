paleokriging <- function(dat, coord, y, x, new.coord, x.new, radius=NULL, k=8, do.R2=FALSE, v.fit=NULL)
{
  #Based on methods described by Martin Striz in 2008
  #These methods were used also to prepare climatic maps in Tolasz et al. 2007
  
  #dat = data frame with paleo and elevation data for each point
  #coord = data frame or matrix with coordinates of points (in metres)
  #y = character name of a variable to be interpolated
  #x = character name of a variable describing altitude (in metres) of each climatic station
  #new.coord = coordinates of new localities for which prediction will be calculated (column names must correspond to those in coord)
  #x.new = numberic vector describing altitude (in metres) of new localities
  #radius = radius (i metres) for local regression; if NULL (default), the k nearest neighbours is used
  #k = number of neighbouring points considered in regression; if radius is NOT NULL, k is not used 
  #do.R2 = weigth predictions with relationship with altitude?
  #v.fit = fitted semivariogram model
  
  require(spdep)
  require(raster)
  require(gstat)
  require(vegan)
  
  if(x %in% colnames(dat) == FALSE)
  {
    stop("Explanatory variable is missing in your data")
  }
  
  if(is.null(v.fit) == TRUE)
  {
    #calculate semivariogram
    f.loc <- paste("~", colnames(coord)[1], "+", colnames(coord)[2], sep="")
    f <- paste(y, "~", 1, sep="")
    v <- variogram(as.formula(f), locations = as.formula(f.loc), data=dat)
    v.fit <- fit.variogram(v, vgm("Exp"))
  }
  
  #calculate neighbourhood for each point
  if(is.null(radius) == TRUE)
  {
    neighs <- knearneigh(as.matrix(coord), k=k)
  }
  else
  {
    neighs <- dnearneigh(as.matrix(coord), d1=0, d2=radius)
  }
  
  a <- vector("numeric")
  b <- vector("numeric")
  res <- vector("numeric")
  R2 <- vector("numeric")
  
  f <- paste(y, "~", x, sep="")
  
  #calculte linear model within the neighbourhood and extract model parametres
  for(i in seq(1, nrow(dat)))
  {
    
    if(is.null(radius) == TRUE)
    {
      sel <- c(i, neighs$nn[i,])
    }
    else
    {
      sel <- c(i, neighs[[i]])
    }
    
    sub <- dat[sel,]
    mod <- lm(as.formula(f), data=sub, na.action="na.exclude")
    a[i] <- coef(mod)[1]
    b[i] <- coef(mod)[2]
    res[i] <- residuals(mod)[1]
    R2[i] <- RsquareAdj(mod)$r.squared
  }
  
  params <- cbind(coord, a, b, res, R2)
  params <- params[complete.cases(params),]
  
  #interpolate model parametres using krige
  f.loc <- paste("~", colnames(coord)[1], "+", colnames(coord)[2], sep="")
  
  a.krige <- krige(a~1, locations = as.formula(f.loc), data=params, newdata=new.coord, model=v.fit)
  b.krige <- krige(b~1, locations = as.formula(f.loc), data=params, newdata=new.coord, model=v.fit)
  res.krige <- krige(res~1, locations = as.formula(f.loc), data=params, newdata=new.coord, model=v.fit)
  
  if(do.R2)
  {
    R2.krige <- krige(R2~1, locations = as.formula(f.loc), data=params, newdata=new.coord, model=v.fit)
    dat2 <- cbind(dat,coord)
    y.krige <- krige(as.formula(paste(y, "~1", sep="")), locations = as.formula(f.loc), data=dat2[is.na(dat2[, y]) == FALSE,], newdata=new.coord, model=v.fit)
  }
  
  #calculate predicted values
  result <- a.krige$var1.pred + (x.new * b.krige$var1.pred) + res.krige$var1.pred
  
  if(do.R2)
  {
    result.R2 <- (R2.krige$var1.pred * result) + (1 - R2.krige$var1.pred) * y.krige$var1.pred
    return(result.R2)
  }
  else
  {
    return(result)
  }
}
