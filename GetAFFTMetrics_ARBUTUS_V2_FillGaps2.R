library(terra);library(AFFT)
rm(list = ls())
gc()
localpath<-"D:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath2<-"N:/mpsg_naip_query2"
rawpath1<-"N:/mpsg_naip"
rawpath3<-"N:/mpsg_naip_missing"

localrawpath<-"C:/TEMP/0_raw"

donepaths<-list.files("N:/mpsg_naip_afft/2_aggregated",full.names = T)
names(donepaths)<-list.files("N:/mpsg_naip_afft/2_aggregated")
fll<-list.files(donepaths,recursive = T, pattern = "tif")
donefiles<-lapply(fll,list.files,pattern = ".tif")
donefiles<-table(do.call(c,donefiles))
missing.files<-list.files(rawpath3,full.names = F);names(missingfiles)<-gsub(".tif","",list.files(rawpath3))#names(donefiles[donefiles <7])



indpath<-paste(localcpath,"1_intermediate",sep = "/")
indpaths<-paste(indpath,c("ndvi","bri","ndgr","ndng"),sep = "/")

aggpath<-paste(localpath,"2_aggregated",sep = "/")
aggpaths<-paste(aggpath,c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")


rawfiles<-missing.files
  
    rl<-list.files(localrawpath,pattern = ".tif")
    copy.me<-rawfiles[!rawfiles %in% list.files(localrawpath)]
    tmp<-pbapply::pblapply(copy.me,function(x){file.copy(paste(rawpath3,x,sep = "/"),paste(localrawpath,x,sep = "/"),overwrite = F)})
    curfilelist<-rawfiles<-list.files(localrawpath,pattern = ".tif")
    rawfiles<-paste(localrawpath,rawfiles,sep = "/");names(rawfiles)<-curfilelist
    rawfiles<-rawfiles[missing.files]
    terraOptions(memfrac = .9,datatype = "INT1U")
    gc()
     indexfiles<-table(list.files(indpaths,full.names = F))
    indexfiles<-names(indexfiles[indexfiles %in% 4])
    need.index<-rawfiles[!names(rawfiles) %in% indexfiles]

    if(length(need.index) > 0){
      
      GetBandIndices(need.index,ncpu = 4,outpath = localcpath)
    #   
    #   library(snowfall)
    # rawpath<-rawpath3
    #   sfInit(parallel = T, cpus = 15)
    #   sfLibrary(AFFT)
    #   sfLibrary(terra)
    #   sfExport(list =c("rawpath","aggpath","localpath","localcpath","indpath","indpaths","rawfiles"))
    #   stime<-Sys.time()
    #   sfExport("stime")
    #   write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
    #  
    #   sfClusterApplyLB(need.index,function(x){
    #     terraOptions(memfrac = .9/5,datatype = "FLT4S")
    #     GetBandIndices(x,ncpu = 1,outpath = localcpath)
    #     write(x,paste(localpath,"donefiles.txt",sep = "/"),append = T)
    #     donefiles<-readLines(paste(localpath,"donefiles.txt",sep = "/"))
    #     donefiles<-donefiles[2:length(donefiles)]
    #     te<-difftime(Sys.time(), stime,units = "secs")
    #     prop.done<-length(donefiles)/length(rawfiles)
    #     tr<-te*(1-prop.done)/(prop.done)
    #     hms<-SDMap::SecsToHMS(tr)
    #     write("Progress on Raw Indices",paste(localpath,"progress.txt",sep = "/"),append = F)
    #     write(paste(round(prop.done * 100,3),"% done"),paste(localpath,"progress.txt",sep = "/"),append = T)
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
    }
  
  indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3]);names(indexfiles)<-indexfiles
  indexfiles<-indexfiles[missing.files]
  
  aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])
  aggfiles<-list.files(paste(aggpath,list.files(aggpath),sep = "/"),full.names = F)
  aggfiles<-names(table(aggfiles)[table(aggfiles)==7])
  indexfiles<-indexfiles[!indexfiles %in% aggfiles]
  
  if(length(indexfiles)>0){
    library(snowfall)
    sfInit(parallel = T, cpus = 17)
    sfLibrary(AFFT)
    sfLibrary(terra)
    rawpath<-rawpath1
    sfExport(list =c("rawpath","aggpath","localpath","indexfiles","localrawpath","localcpath","indpaths","rawfiles"))
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
      prop.done<-length(donefiles)/length(rawfiles)
      tr<-te*(1-prop.done)/(prop.done)
      hms<-SDMap::SecsToHMS(tr)
      
      write("Progress on AFFTs",paste(localpath,"progress.txt",sep = "/"),append = F)
      write(paste(round(prop.done * 100,3),"% done"),paste(localpath,"progress.txt",sep = "/"),append = T)
      write("",paste(localpath,"progress.txt",sep = "/"),append = T)
      hms<-SDMap::SecsToHMS(te)
      
      write("Time Elapsed:",paste(localpath,"progress.txt",sep = "/"),append = T)
      write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
        
      hms<-SDMap::SecsToHMS(tr)
      
      write("Time Remaining:",paste(localpath,"progress.txt",sep = "/"),append = T)
      write(paste(paste(hms,names(hms)),collapse = "   "),
            paste(localpath,"progress.txt",sep = "/"),append = T)
  
    })
    sfStop()
  }
  
  outfilelist<-list.files("D:/LocalNaip/2_aggregated",full.names = T, recursive = T)
  copytolist<-gsub("D:/LocalNaip","N:/mpsg_naip_afft",outfilelist)
  pbapply::pblapply(1:length(outfilelist),function(x){file.copy(outfilelist[x],copytolist[x],overwrite = T)})
  if(all(file.exists(copytolist))){
    file.remove(outfilelist)
    file.remove(list.files(localrawpath,full.names = T))
    file.remove(list.files(indpath,full.names = T,recursive = T))
  }else{print("MISSING FILES! CLEANUP STOPPED");browser()}
  
    
  rawfilebatch<-read.csv("N:/mpsg_naip_batch_files/FileLog.csv")
  rawfilebatch$Status[rawfilebatch$rawfiles%in%curfilelist]<-"Complete"
  rawfilebatch$Location[rawfilebatch$rawfiles%in%curfilelist]<-"Copied To Sequoia"
  write.csv(rawfilebatch,"N:/mpsg_naip_batch_files/FileLog.csv",row.names = F)
  
}

## Spot-check aggregated files -------------
# 
# fl<-list.files("2_aggregated/ndvi",full.names = T)
# r1<-rast(fl[3])
# plot(r1,col = rainbow(255))
review.images<-F
if(review.images){
   
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

}


## Things to think about ----------------
### 1) Workflow. Should GetBandIndices be rolled in to GetAFFT?
### 2) Parallel handling. GetBandIndices currently uses 4 cpus (one per index). Could re-jigger the parallelization to send each image out to a separate cpu, but how would memory work?
### 3) File structure -- improvements needed?
### 4) Tracking -- how to keep track of what is and isn't done when working on multiple platforms?

