library(terra);library(AFFT)
rm(list = ls())
gc()
localpath<-"D:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath<-"J:/WildCarrot/2018_1ft_NAIP_raw"#N:/mpsg_naip_query2"
localrawpath<-"C:/TEMP/0_raw"
indpath<-paste(localcpath,"1_intermediate",sep = "/")
indpaths<-paste(indpath,c("ndvi","bri","ndgr","ndng"),sep = "/")

aggpath<-paste(localpath,"2_aggregated",sep = "/")
aggpaths<-paste(aggpath,c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")
batchfile<-read.csv(paste(localpath,"/OregonTestBatchFiles.csv",sep = ""),row.names = 1)


#write.csv(batchfile,paste(localpath,"/OregonTestBatchFiles.csv",sep = ""), row.names = T)
batchfile<-batchfile[!batchfile$done,]
while(nrow(batchfile)>0){
  cat("start while loop")
  batchlist<-unique(batchfile$batches)
  for(batch in batchlist){
    cat("\n  ####################\n    ## Batch:",batch," ##\n####################\n")
    curbatch<-batchfile[batchfile$batches == batch,]
    rawfiles<-curbatch$tifnames
    rl<-list.files(localrawpath)
    copy.me<-rawfiles[!rawfiles %in% list.files(localrawpath)]
    tmp<-file.copy(curbatch[copy.me,"rawfiles"],localrawpath)
    rawfiles<-list.files(localrawpath,full.names = T);names(rawfiles)<-list.files(localrawpath,full.names = F)
#    terraOptions(memfrac = .9,datatype = "INT1U")
    gc()
    ## Get indices for raw naip files -------------
    indexfiles<-table(list.files(indpaths,full.names = F))
    indexfiles<-names(indexfiles[indexfiles %in% 4])
    need.index<-rawfiles[!names(rawfiles) %in% indexfiles]
    if(length(need.index) > 0){
      cat("getting index files                          \n")
      library(snowfall)
      sfInit(parallel = T, cpus = 15)
      sfLibrary(AFFT)
      sfLibrary(terra)
      sfExport(list =c("rawpath","aggpath","localpath","localcpath","indpath","indpaths","rawfiles"))
      stime<-Sys.time()
      sfExport("stime")
      write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
      
      sfClusterApplyLB(need.index,function(x){
        terraOptions(memfrac = .9/5,datatype = "FLT4S")
        GetBandIndices(x,ncpu = 1,outpath = localcpath)
        write(x,paste(localpath,"donefiles.txt",sep = "/"),append = T)
        donefiles<-readLines(paste(localpath,"donefiles.txt",sep = "/"))
        donefiles<-donefiles[2:length(donefiles)]
        te<-difftime(Sys.time(), stime,units = "secs")
        prop.done<-length(donefiles)/length(rawfiles)
        tr<-te*(1-prop.done)/(prop.done)
        hms<-SDMap::SecsToHMS(tr)
        write("Progress on Raw Indices",paste(localpath,"progress.txt",sep = "/"),append = F)
        write(paste(round(prop.done * 100,3),"% done"),paste(localpath,"progress.txt",sep = "/"),append = T)
        write("",paste(localpath,"progress.txt",sep = "/"),append = T)
        hms<-SDMap::SecsToHMS(te)
        
        write("Time Elapsed:",paste(localpath,"progress.txt",sep = "/"),append = T)
        write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
        
        hms<-SDMap::SecsToHMS(tr)
        
        write("Time Remaining:",paste(localpath,"progress.txt",sep = "/"),append = T)
        write(paste(paste(hms,names(hms)),collapse = "   "),paste(localpath,"progress.txt",sep = "/"),append = T)
      })
      sfStop()
    }
    cat("got index files                          \n")
  
  
    indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3])
    aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])
    aggfiles<-list.files(paste(aggpath,list.files(aggpath),sep = "/"),full.names = F)
    aggfiles<-names(table(aggfiles)[table(aggfiles)==7])
    indexfiles<-indexfiles[!indexfiles %in% aggfiles]

    if(length(indexfiles)>0){
      
      cat("getting agg files                          \n")
      library(snowfall)
      sfInit(parallel = T, cpus = 22)
      sfLibrary(AFFT)
      sfLibrary(terra)
      sfExport(list =c("rawpath","aggpath","localpath","indexfiles","localrawpath","localcpath","indpaths","rawfiles"))
      stime<-Sys.time()
      sfExport("stime")
      write("",paste(localpath,"donefiles.txt",sep = "/"),append = F)
      sfClusterApplyLB(indexfiles,function(x){
        GetAFFT(x,
                zradii = c(0.4,1,5,25,60),
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
    
    cat("got agg files                          \n")
    ## Batch Tracking, and C-drive cleanup ---------------------
    
    cat("recording and cleanup                          \n")
    batchfile<-read.csv(paste(localpath,"/OregonTestBatchFiles.csv",sep = ""),row.names = 1)
    batchfile[batchfile$tifnames %in% c(aggfiles,indexfiles),"done"]<-T
    write.csv(batchfile,paste(localpath,"/OregonTestBatchFiles.csv",sep = ""), row.names = T)
    batchfile<-batchfile[!batchfile$done,]
    rmlist<-unique(c(indexfiles,aggfiles))
    tmp<-file.remove(paste(localrawpath,rmlist,sep = "/"))
    tmp<-lapply(indpaths,function(x){file.remove(paste(x,rmlist,sep = "/"))})
  } 
} 
## Things to think about ----------------
### 1) Workflow. Should GetBandIndices be rolled in to GetAFFT?
### 2) Parallel handling. GetBandIndices currently uses 4 cpus (one per index). Could re-jigger the parallelization to send each image out to a separate cpu, but how would memory work?
### 3) File structure -- improvements needed?
### 4) Tracking -- how to keep track of what is and isn't done when working on multiple platforms?

