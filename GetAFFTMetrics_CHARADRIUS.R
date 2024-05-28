library(terra);library(AFFT)
rm(list = ls())
gc()
localpath<-"D:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath<-"N:/mpsg_naip"
localrawpath<-"C:/TEMP/0_raw"
aggpath<-"D:/LocalNaip/2_aggregated"

indpath<-paste(localcpath,"1_intermediate",sep = "/")
indpaths<-paste(indpath,c("ndvi","bri","ndgr","ndng"),sep = "/")

aggpaths<-paste(aggpath,c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")



batches<-read.csv("N:/mpsg_naip_batch_files/BatchLog_v2.csv")
print(batches)

curbatches<-14
batches$Status[curbatches]<-"Running"
batches$Location[curbatches]<-"Charadrius"
write.csv(batches,"N:/mpsg_naip_batch_files/BatchLog_v2.csv",row.names = F)

rawfiles<-do.call(c,lapply(curbatches,function(curbatch){read.csv(paste("N:/mpsg_naip_batch_files/naip_batch_",curbatch,".txt",sep = ""))[,1]}))
fileinfo<-file.info(rawfiles)

rfl<-split(rawfiles,round(cumsum(fileinfo$size)/(10^9)/75))

rawfiles<-rfl[[1]]

#pbapply::pblapply(indexfiles,function(x){file.copy(paste(rawpath,x,sep = "/"),paste(localrawpath,x,sep = "/"))})
rawfiles<-batchlist<-sapply(rawfiles,function(x){x<-strsplit(x,"/")[[1]][3];x})
#indexfiles<-list.files(ndvipath)
indfiles<-table(list.files(indpaths));indfiles<-names(indfiles[indfiles>3])
aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])

rfl<-rawfiles<-rawfiles[!rawfiles %in% indfiles]

pbapply::pblapply(rawfiles,function(x){file.copy(paste(rawpath,x,sep = "/"),paste(localrawpath,x,sep = "/"))})

rawfiles<-paste(localrawpath,rawfiles,sep = "/");names(rawfiles)<-rfl

terraOptions(memfrac = .9,datatype = "INT1U")


setwd(localpath)

GetBandIndices(rawfiles,ncpu = 4,outpath = localcpath)
# gc()
# library(snowfall)
# time2<-system.time({
#   sfInit(parallel = T, cpus = 7)
#   sfLibrary(AFFT)
#   sfLibrary(terra)
#   sfExport(list =c("rawpath","aggpath","localpath","localcpath","indpath","rawfiles"))
#   stime<-Sys.time()
#   sfExport("stime")
#   write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
#   sfClusterApplyLB(rawfiles,function(x){
#     terraOptions(datatype = "FLT4S")
#     GetBandIndices(x,ncpu = 1,outpath = localcpath)
#     write(x,paste(localpath,"donefiles.txt",sep = "/"),append = T)
#     donefiles<-readLines(paste(localpath,"donefiles.txt",sep = "/"))
#     donefiles<-donefiles[2:length(donefiles)]
#     te<-difftime(Sys.time(), stime,units = "secs")
#     prop.done<-length(donefiles)/length(rawfiles)
#     tr<-te*(1-prop.done)/(prop.done)
#     hms<-SDMap::SecsToHMS(tr)
#     write(paste(round(prop.done * 100,3),"% done"),paste(localpath,"progress.txt",sep = "/"),append = F)
#     write("",paste(localpath,"progress.txt",sep = "/"),append = T)
#     hms<-SDMap::SecsToHMS(te)
#     
#     write("Time Elapsed:",paste(localpath,"progress.txt",sep = "/"),append = T)
#     write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
#     
#     hms<-SDMap::SecsToHMS(tr)
#     
#     write("Time Remaining:",paste(localpath,"progress.txt",sep = "/"),append = T)
#     write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
#   })
#   sfStop()
# })



indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3])
aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])

indexfiles<-indexfiles[indexfiles %in% batchlist]
aggfiles<-list.files(paste(aggpath,list.files(aggpath),sep = "/"),full.names = F)
aggfiles<-names(table(aggfiles)[table(aggfiles)==7])
indexfiles<-indexfiles[!indexfiles %in% aggfiles]
ind<-0

# time1<-system.time({GetAFFT(indexfiles[1:10],
#         rawpath = rawpath,
#         indpath = paste(localpath,"1_intermediate",sep = "/"),
#         aggpath = paste(localpath,"2_aggregated",sep = "/"),
#         ncpu = 7,overwrite = T)})
library(snowfall)

localrawpath<-"C:/TEMP/0_raw"
localcpath<-"C:/TEMP"
time2<-system.time({
  sfInit(parallel = T, cpus = 1)
  sfLibrary(AFFT)
  sfLibrary(terra)
  sfExport(list =c("rawpath","aggpath","localpath","localrawpath","localcpath","indexfiles"))
  stime<-Sys.time()
  sfExport("stime")
  write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
  sfLapply(indexfiles,function(x){
    browser()
    GetAFFT(x,
              rawpath = localrawpath,
              indpath = paste(localcpath,"1_intermediate",sep = "/"),
              aggpath = paste(localpath,"2_aggregated",sep = "/"),
              ncpu = 1,overwrite = T
              )
    write(x,paste(localpath,"donefiles.txt",sep = "/"),append = T)
    donefiles<-readLines(paste(localpath,"donefiles.txt",sep = "/"))
    donefiles<-donefiles[2:length(donefiles)]
    te<-difftime(Sys.time(), stime,units = "secs")
    prop.done<-length(donefiles)/length(indexfiles)
    tr<-te*(1-prop.done)/(prop.done)
    hms<-SDMap::SecsToHMS(tr)
    write(paste(round(prop.done * 100,3),"% done"),paste(localpath,"progress.txt",sep = "/"),append = F)
    write("",paste(localpath,"progress.txt",sep = "/"),append = T)
    hms<-SDMap::SecsToHMS(te)
    
    write("Time Elapsed:",paste(localpath,"progress.txt",sep = "/"),append = T)
    write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
    
    hms<-SDMap::SecsToHMS(tr)
    
    write("Time Remaining:",paste(localpath,"progress.txt",sep = "/"),append = T)
    write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
  })
  sfStop()
})


## Spot-check aggregated files -------------
# 
# fl<-list.files("2_aggregated/ndvi",full.names = T)
# r1<-rast(fl[3])
# plot(r1,col = rainbow(255))

band<-"ndvi"#readline("Band name:")
cols<-rainbow(255,alpha = ((1:255)+200)/455)

lastind<-length(aggfiles)
aggfiles<-table(list.files(aggpaths,recursive = T, pattern = ".tif",full.names = F))
aggfiles<-names(aggfiles[aggfiles == 7])


junk<-lapply(aggfiles,function(x,p = 2,band = "ndvi"){
  r0<-rast(paste(rawpath,x,sep = "/"))
  r2<-rast(paste(localpath,"1_intermediate",band,x,sep = "/"))
  r3<-rast(paste(localpath,"2_aggregated",band,x,sep = "/"))
  if(p==2){
        par(mfrow =c(2,2))
   # r4<-app(r3,function(x){c(x[1],x[2],sum(x[3:4]))/sum(x[1:4])*255})
    r4<-subset(r3,c(7:9))
  
    #e1<-ext(r4)#/5
    plotRGB(r4,stretch = "hist",smooth = F)
    plotRGB(c(r2,r2,r2),stretch = "lin",smooth = F)
    #plotRGB(r4,ext = e1)
    #plotRGB(r4,stretch = "hist",ext = e1)
    plotRGB(c(r3$mean,r3$mean,r3$mean),stretch = "lin",smooth = F)
    plotRGB(r0,smooth = F)#,stretch = 'lin')
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

