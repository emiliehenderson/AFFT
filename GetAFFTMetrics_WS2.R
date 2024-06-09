library(terra);library(AFFT);library(SDMap)
rm(list = ls())
gc()
localpath<-"D:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath<-"N:/mpsg_naip"

localrawpath<-MakeDirIfNeeded("0_raw",localcpath,goto = F)
indpath<-MakeDirIfNeeded("1_intermediate",localcpath,goto = T)
indpaths<-sapply(c("ndvi","ndgr","ndng","bri"),function(x){MakeDirIfNeeded(x,goto = F)})

aggpath<-MakeDirIfNeeded("2_aggregated",localpath,goto = T)
aggpaths<-sapply(c("r","g","n","ndvi","ndgr","ndng","bri"),function(x){MakeDirIfNeeded(x,goto = F)})
aggpaths2<-paste("N:/mpsg_naip_afft/2_aggregated",c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")



batches<-read.csv("N:/mpsg_naip_batch_files/BatchLog_v2.csv")
View(batches)

curbatches<-c(14)
rawfiles<-do.call(c,lapply(curbatches,function(curbatch){read.csv(paste("N:/mpsg_naip_batch_files/naip_batch_",curbatch,".txt",sep = ""))[,1]}))
alldone<-list.files(aggpaths2)
rawfiles<-rawfiles[!rawfiles %in% paste("N:/mpsg_naip/",alldone,sep = "")]

fileinfo<-file.info(rawfiles)

rfs<-split(rawfiles,round(cumsum(fileinfo$size)/(10^9)/75))

for(i in 1:length(rfs)){
  cat("#############################\n",i,"out of", length(rfs),"chunks in batch", curbatches,"\n")
  
  rawfiles<-rfs[[i]]
  rawfiles<-batchlist<-sapply(rawfiles,function(x){x<-strsplit(x,"/")[[1]][3];x})
  indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3])
  aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])
  rfl<-rawfiles<-rawfiles[!rawfiles %in% c(indexfiles,aggfiles,alldone)]
  
  if(length(rawfiles)>0){
  pbapply::pblapply(rawfiles,function(x){file.copy(paste(rawpath,x,sep = "/"),paste(localrawpath,x,sep = "/"),overwrite = F)})
  
  rawfiles<-paste(localrawpath,rawfiles,sep = "/");names(rawfiles)<-rfl
  
  terraOptions(memfrac = .9,datatype = "INT1U")
  gc()
  library(snowfall)

  rawfiles<-gsub("N:/mpsg_naip","C:/TEMP/0_raw",rawfiles)
  time2<-system.time({
  sfInit(parallel = T, cpus = 10)
  sfLibrary(AFFT)
  sfLibrary(terra)
  sfExport(list =c("rawpath","aggpath","localpath","localcpath","indpath","rawfiles"))
  stime<-Sys.time()
  sfExport("stime")
  write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
  sfClusterApplyLB(rawfiles,function(x){
    terraOptions(memfrac = .9/5,datatype = "FLT4S")
    GetBandIndices(x,ncpu = 1,outpath = localcpath)
    write(x,paste(localpath,"donefiles.txt",sep = "/"),append = T)
    donefiles<-readLines(paste(localpath,"donefiles.txt",sep = "/"))
    donefiles<-donefiles[2:length(donefiles)]
    te<-difftime(Sys.time(), stime,units = "secs")
    prop.done<-length(donefiles)/length(rawfiles)
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
}


  indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3])
  aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])
  
  indexfiles<-indexfiles[indexfiles %in% batchlist]
  aggfiles<-list.files(paste(aggpath,list.files(aggpath),sep = "/"),full.names = F)
  aggfiles<-names(table(aggfiles)[table(aggfiles)==7])
  indexfiles<-indexfiles[!indexfiles %in% aggfiles]
  
  library(snowfall)
  # GetAFFT(indexfiles[1],
  #         rawpath = localrawpath,
  #         indpath = paste(localcpath,"1_intermediate",sep = "/"),
  #         aggpath = paste(localpath,"2_aggregated",sep = "/"),
  #         ncpu = 1,overwrite = T)
  time2<-system.time({
    sfInit(parallel = T, cpus = 15)
    sfLibrary(AFFT)
    sfLibrary(terra)
    sfExport(list =c("rawpath","aggpath","localpath","indexfiles","localrawpath","localcpath"))
    stime<-Sys.time()
    sfExport("stime")
    write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
    sfClusterApplyLB(indexfiles,function(x){
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
  
  setwd(localrawpath);fl<-list.files();file.remove(fl);setwd("..")
  fl<-list.files(indpaths,full.names = T);file.remove(fl)
  setwd(localpath)
}

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

