library(terra);library(AFFT)
rm(list = ls())
gc()
localpath<-"D:/LocalNaip"
rawpath<-"N:/mpsg_naip"
indpath<-paste(localpath,"1_intermediate",sep = "/")
ndvipath<-paste(indpath,"ndvi",sep = "/")
aggpath<-paste(localpath,"2_aggregated",sep = "/")
ndvipath_a<-paste(aggpath,"ndvi",sep = "/")

rfl<-rawfiles<-list.files(rawpath,pattern = ".tif")
indexfiles<-list.files(ndvipath)

rawfiles<-list.files(rawpath,full.names = T,pattern = ".tif")
names(rawfiles)<-rfl

rawfiles<-rawfiles[!rfl %in% indexfiles]

terraOptions(memfrac = .9,datatype = "INT1U")

setwd(localpath)
GetBandIndices(rawfiles[3],ncpu = 4)

indexfiles<-list.files(ndvipath,full.names = T)
aggfiles<-paste(ndvipath,list.files(ndvipath_a,full.names = F),sep = "/")
indexfiles<-indexfiles[!indexfiles %in% aggfiles]
ind<-0

GetAFFT(indexfiles,
        rawpath = rawpath,
        indpath = paste(localpath,"1_intermediate",sep = "/"),
        aggpath = paste(localpath,"2_aggregated",sep = "/"),
        ncpu = 7,overwrite = T
        )
# 
# fl<-list.files("2_aggregated/ndvi",full.names = T)
# r1<-rast(fl[3])
# plot(r1,col = rainbow(255))

band<-"ndvi"#readline("Band name:")
cols<-rainbow(255,alpha = ((1:255)+200)/455)

lastind<-length(aggfiles)
aggfiles<-list.files(paste(localpath,"2_aggregated/ndvi",sep ="/"),
                     pattern = ".tif",full.names = F)
aggfiles<-"m_3710261_ne_13_060_20210729.tif"
junk<-lapply(aggfiles[lastind:length(aggfiles)],function(x,p = 2,band = "ndvi"){
  r0<-rast(paste(rawpath,x,sep = "/"))
  r2<-rast(paste(localpath,"1_intermediate",band,x,sep = "/"))
  r3<-rast(paste(localpath,"2_aggregated",band,x,sep = "/"))
  if(p==2){
        par(mfrow =c(2,2))
   # r4<-app(r3,function(x){c(x[1],x[2],sum(x[3:4]))/sum(x[1:4])*255})
    r4<-subset(r3,c(7:9))
  
    e1<-ext(r4)#/5
    plotRGB(r4,ext = e1,smooth = F,stretch = "lin")
    plotRGB(c(r2,r2,r2),stretch = "lin",ext = e1,smooth = F)
    #plotRGB(r4,ext = e1)
    #plotRGB(r4,stretch = "hist",ext = e1)
    plotRGB(c(r3$mean,r3$mean,r3$mean),stretch = "lin",ext = e1,smooth = F)
    plotRGB(r0,ext = e1,smooth = F)#,stretch = 'lin')
  }else if (p==1){
    plot(subset(r3,1:12), col = cols,smooth = F)
  }

  readline(paste(which(aggfiles==x),"out of",length(aggfiles)))
  })



## Things to think about ----------------
### 1) Workflow. Should GetBandIndices be rolled in to GetAFFT?
### 2) Parallel handling. GetBandIndices currently uses 4 cpus (one per index). Could re-jigger the parallelization to send each image out to a separate cpu, but how would memory work?
### 3) File structure -- improvements needed?
### 4) Tracking -- how to keep track of what is and isn't done when working on multiple platforms?

