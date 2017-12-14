#install.packages(c("gstat", "sp))


#####################################################################
###################### KRIGING FUNCTION #############################
#####################################################################
KRIG <- function(data, x, y, z, grid=NULL){
 library(sp)
 library(gstat)
 

 df <- data
 colnames(df) <- gsub(x, "x", colnames(df))
 colnames(df) <- gsub(y, "y", colnames(df))
 colnames(df) <- gsub(z, "z", colnames(df))

 coordinates(df) <- ~ x + y

 z.vgm <- variogram(z~1, df) # calculates sample variogram values 
 z.fit <- fit.variogram(z.vgm, model=vgm(1, "Sph", 900, 1)) # fit model


 #plot(z.vgm, z.fit)
 
 if(is.null(grid)==T){
  df.grid <- data.frame(x=seq(min(df$x),max(df$x),length.out=100),
				y=seq(min(df$y),max(df$y),length.out=100))
  df.grid <- expand.grid(df.grid)
 } else { df.grid = grid }


 coordinates(df.grid) <- ~ x + y
 z.kriged <- krige(z ~ 1, df, df.grid, model=z.fit)

 z.kriged <- as.data.frame(z.kriged)


 par(mfrow=c(1,4))

 plot(y~x, data=df, cex=(df$z/max(df$z))+1,
	main="Values of 'z', represented by point size")
 plot(y~x, data=df, pch=19, main="Points w/ Measurements")
 plot(y~x, data=df.grid, main="Points at which to estimate",
	pch=19, cex=.7)
 plot(y~x, data=z.kriged, pch=19,
	col=rgb((min(z.kriged$var1.pred)/z.kriged$var1.pred),
		0,0,(z.kriged$var1.pred/max(z.kriged$var1.pred))),
	main="Kriged Points")
}

######################################################################
######################################################################
######################################################################

######################################################################
## The "meuse" dataset comes from the package 'sp'
## Meuse contains chemical concentrations at various coordinate
## locations. Useful for practicing kriging

KRIG(data=meuse, x="x", y="y", z="zinc", grid=meuse.grid)
KRIG(data=meuse, x="x", y="y", z="zinc")









#########################################################################
############### PRACTICING JUST WITH THE MEUSE DATASET ##################
############### without the function "KRIG"            ##################
#########################################################################
df <- meuse
head(df)
str(df)

par(mfrow=c(1,1))
plot(y~x, data=df, cex=(df$zinc/max(df$zinc))+1)

coordinates(df) <- ~ x + y
 str(df)

lzn.vgm <- variogram(log(zinc)~1, df) # calculates sample variogram values 
lzn.fit <- fit.variogram(lzn.vgm, model=vgm(1, "Sph", 900, 1)) # fit model

plot(lzn.vgm, lzn.fit)
df.grid <- meuse.grid
head(df.grid)
head(df)

par(mfrow=c(1,2))
plot(y~x, data=df, pch=19, main="Points w/ Measurements")
plot(y~x, data=df.grid, main="Points at which to estimate",
	pch=19, cex=.7)

coordinates(df.grid) <- ~ x + y
lzn.kriged <- krige(log(zinc) ~ 1, df, df.grid, model=lzn.fit)

lzn.kriged <- as.data.frame(lzn.kriged)
str(lzn.kriged)
head(lzn.kriged)

range(lzn.kriged$var1.pred)


par(mfrow=c(1,1))
plot(y~x, data=lzn.kriged, pch=19,
	col=rgb((min(lzn.kriged$var1.pred)/lzn.kriged$var1.pred),
		0,0,(lzn.kriged$var1.pred/max(lzn.kriged$var1.pred))))














