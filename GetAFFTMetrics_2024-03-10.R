library(terra);library(AFFT)

rawfiles<-list.files("0_raw",full.names = F)
indexfiles<-list.files("1_intermediate/ndvi")#,full.names = F))
rawfiles<-rawfiles[!rawfiles %in% indexfiles]
rawfiles<-paste("0_raw",rawfiles,sep = "/")
terraOptions(memfrac = .9,datatype = "INT1U")

#GetBandIndices(rawfiles)


indexfiles<-list.files("1_intermediate/ndvi")#,full.names = T)
aggfiles<-list.files("2_aggregated/ndvi")
indexfiles<-indexfiles[!indexfiles %in% aggfiles]
ind<-0

#GetAFFT(indexfiles)
# 
# fl<-list.files("2_aggregated/ndvi",full.names = T)
# r1<-rast(fl[3])
# plot(r1,col = rainbow(255))

band<-"ndvi"#readline("Band name:")
cols<-rainbow(255,alpha = ((1:255)+200)/455)

lastind<-length(aggfiles)
aggfiles<-list.files("2_aggregated/ndvi",
                     pattern = ".tif",full.names = F)

junk<-lapply(aggfiles[lastind:length(aggfiles)],function(x){
  r0<-rast(paste("0_raw",x,sep = "/"))
  r2<-rast(paste("1_intermediate",band,x,sep = "/"))
  r3<-rast(paste("2_aggregated",band,x,sep = "/"))
  #plot(subset(r3,1:6), range =c(0,255),col = cols)
  par(mfrow =c(2,3))
  #r4<-app(r3,function(x){c(x[1],x[2],sum(x[3:4]))/sum(x[1:4])*255})
  r4<-subset(r3,1:3+2)
  plotRGB(r4)
  plotRGB(r4,stretch = "lin")
  plotRGB(r4,stretch = "hist")
  plotRGB(r0)#,stretch = 'lin')
  plotRGB(c(r2,r2,r2),stretch = "lin")
  plotRGB(c(r3$mean,r3$mean,r3$mean),stretch = "lin")
  readline(paste(which(aggfiles==x),"out of",length(aggfiles)))
  })


