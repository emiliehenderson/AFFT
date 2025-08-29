## Setup -------------------------------------
  setwd("C:/TEMP")
  library(terra)
  library(AFFT)
  library(snowfall)
localpath<-"D:/LocalNaip"
localcpath<-"C:/TEMP"
rawpath<-"J:/WildCarrot/2018_1ft_NAIP_raw"#N:/mpsg_naip_query2"
localrawpath<-"C:/TEMP/0_raw"
indpath<-paste(localcpath,"1_intermediate",sep = "/")
indpaths<-paste(indpath,c("ndvi","bri","ndgr","ndng"),sep = "/")

aggpath<-paste(localpath,"2_aggregated",sep = "/")
aggpaths<-paste(aggpath,c("r","g","n","ndvi","bri","ndgr","ndng"),sep = "/")
batchfile<-read.csv(paste(localpath,"/OregonTestBatchFiles.csv",sep = ""),row.names = 1)

vl<-lapply(list.files("D:/LocalNaip/2_aggregated",full.names = T),
       function(x){
         flx<-list.files(x,full.names = T)
         y<-strsplit(x,"/")[[1]];y<-y[length(y)]
         v<-vrt(flx)
         nm<-paste(names(rast(flx[1])),y,sep = "_")
        names(v)<-nm
        v
       })
  r1<-rast(vl)

  samp1<-spatSample(r1,100000,"random",xy = F,na.rm = T)
  
  
## identify problematic bands ----------------------


 ### Artifacts ---------------
  if(!file.exists("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")){
    windows(10,10);plot(v2,4);plot(bnd0,add = T);e1<-draw();dev.off()
    artifacts<-data.frame(do.call(rbind,pbapply::pblapply(bb,function(b){

      fl2<-fl[grepl(paste("/",b," ",sep = ""),fl)]
      v1<-rast(fl2)
      par(mfrow =c(2,2))
      plot(v1,main = b); plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin")
      plot(v1,ext = e1);plot(bnd0,add = T)
      plotRGB(c(v1,v1,v1),stretch = "lin",ext = e1)
      plot(bnd0,border = "red",add = T)
      cat("\n################\n ",b,"\n###############\n\n")
      return(c(Flightline = readline("Flightline score: "),Phenology = readline("Phenology Score: "), image = b))
    })))
    write.csv(artifacts,"C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv",row.names = F)
  }else{artifacts<-read.csv("C:/TEMP/ComancheCimarrone/Corrected/artifactscores.csv")}
  

  patchy<-artifacts$image[artifacts$Phenology %in% c(4,5)]
  stripey<-artifacts$image[artifacts$Flightline %in% c(4,5)]
  drop.me<-unique(c(patchy,stripey))
  drop.me<-unique(drop.me)
  
   # nl<-lapply(drop.me,function(n){par(mfrow =c(2,1),mar =c(0,0,0,0));y<-subset(v2,n)
   #    plot(y,main = n,legend = F);plot(vect(e1),add = T)
   #    plot(y,ext = e1,legend = F)
   #    #plot(y,ext = e2,legend = F);plot(vect(e4),add = T)
   #    #plot(y,ext = e4,legend = F)
   #    ;readline(n)})
   
## Scroll through all vrt layers, build rast object vrt containing all layers --------------
  
## Build princomp -------------
  
  
  ### Build PCA -------------------
    drop.me<-list()
    stripey<-patchy<-c()
    drop.me[[1]]<-c("x","y","jdate","sl","imgyr","stateid",unique(c(patchy,stripey)))

    samp1<-samp1[,!colnames(samp1) %in% do.call(c,drop.me)]
    pca1<-princomp(samp1)
     
    
  pred1<-function(x,y,keepvec = which(cumsum(x$sdev/sum(x$sdev))<.991)){
    z<-round(predict(x,y)[,keepvec] * 100,0)}
  
  
  imglist<-list.files(aggpaths[1])
  tl<-pbapply::pblapply(imglist,function(x){
    r<-rast(paste(aggpaths,x,sep = "/"))
    names(r)<-names(r1)
    p<-predict(r,pca1,cores = 1)
    p
  })
  pca<-predict(r1,pca1,fun = pred1,cores = 30)

  names(tl)<-imglist
  
  lapply(imglist,function(x){
    r1<-rast(paste(localrawpath,x,sep = "/"))
    par(mfrow =c(1,3))
    plotRGB(r1,4,1,2, stretch = "lin")
    plotRGB(tl[[x]],1,2,3,stretch = "lin")
    plotRGB(tl[[x]],4,5,6,stretch = "hist")
    readline(x)
  })

  
    