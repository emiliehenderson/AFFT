library(terra);library(AFFT)
rm(list = ls())
gc()
localpath<-"E:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath<-"N:/mpsg_naip_query2"
localrawpath<-"C:/TEMP/0_raw"
indpath<-paste(localcpath,"1_intermediate",sep = "/")
indpaths<-paste(indpath,c("ndvi","bri","ndgr","ndng"),sep = "/")

aggpath<-paste(localpath,"2_aggregated",sep = "/")
aggpaths<-paste(aggpath,c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")


rawfilebatch<-read.csv("N:/mpsg_naip_batch_files/FileLog.csv")

 rawfilebatch[rawfilebatch$Location %in% "Castilleja",c(2:3)]<-NA
 write.csv(rawfilebatch,"N:/mpsg_naip_batch_files/FileLog.csv",row.names = F)

while(any(rawfilebatch$Status %in% NA)){
  rawfilebatch<-read.csv("N:/mpsg_naip_batch_files/FileLog.csv")
  rawfiles<-rawfilebatch[is.na(rawfilebatch$Status),1]
  curfilelist<-rawfiles<-rawfiles[1:50]
  rawfilebatch$Status[rawfilebatch$rawfiles%in%rawfiles]<-"Running"
  rawfilebatch$Location[rawfilebatch$rawfiles%in%rawfiles]<-"Castilleja"
  write.csv(rawfilebatch,"N:/mpsg_naip_batch_files/FileLog.csv",row.names = F)
  cat("Files Completed:",sum(rawfilebatch$Status %in% "Complete"),"\n")
  cat("Files Remaining:",sum(rawfilebatch$Status %in% NA),"\n")
  cat("Current Batch:", length(rawfiles),"\n")
  
  if(length(rawfiles)>0){
    rl<-list.files(localrawpath)
    copy.me<-rawfiles[!rawfiles %in% list.files(localrawpath)]
    tmp<-pbapply::pblapply(copy.me,function(x){file.copy(paste(rawpath,x,sep = "/"),paste(localrawpath,x,sep = "/"),overwrite = F)})
    rawfiles<-paste(localrawpath,rawfiles,sep = "/");names(rawfiles)<-curfilelist
    
    terraOptions(memfrac = .9,datatype = "INT1U")
    gc()
     indexfiles<-table(list.files(indpaths,full.names = F))
    indexfiles<-names(indexfiles[indexfiles %in% 4])
    need.index<-rawfiles[!names(rawfiles) %in% indexfiles]

    if(length(need.index) > 0){
        terraOptions(datatype = "FLT4S")
        GetBandIndices(need.index,ncpu = 4,outpath = localcpath)
    }
  }
  indexfiles<-table(list.files(indpaths));indexfiles<-names(indexfiles[indexfiles>3])
  aggfiles<-table(list.files(aggpaths));aggfiles<-names(aggfiles[aggfiles>6])
  aggfiles<-list.files(paste(aggpath,list.files(aggpath),sep = "/"),full.names = F)
  aggfiles<-names(table(aggfiles)[table(aggfiles)==7])
  need.aggfiles<-indexfiles[!indexfiles %in% aggfiles]
  
  if(length(indexfiles)>0){
    
      GetAFFT(need.aggfiles,
              rawpath = localrawpath,
              indpath = paste(localcpath,"1_intermediate",sep = "/"),
              aggpath = paste(localpath,"2_aggregated",sep = "/"),
              ncpu = 7,overwrite = T
      )
      
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


