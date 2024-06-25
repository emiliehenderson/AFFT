
## Setup -------------------------------------
setwd("Y:/MPSG_VegMapping/Data/Raster/Source/AFFT")
library(readxl)
library(terra)
seamlines <- data.frame(read_excel("Y:/MPSG_VegMapping/Data/Raster/Source/afft/naip_collection_year_V2.xlsx"))
rownames(seamlines)<-nm<-seamlines$State;names(nm)<-nm

batch<-vect("N:/mpsg_naip_batch_files/batch_spatial/batch_spatial.shp")
## Function -----------------------
CorrectLayers<-function(y,x,maxsample= 100000 ,sampdens = 50,cs = curstate,xy = NULL){
  cat("correcting layers\n")
  if(is.null(xy)){
    cat("\n   Sample step 1        ")
    pts<-spatSample(c(x$JDate),size = max(sampdens,min(maxsample,ncell(y)/sampdens)),xy = T,method = "stratified")
    cat("\r  freq               ")
    f1<-freq(x$JDate)
    f1$prop<-f1$count/sum(f1$count)
    nm<-as.numeric(unique(f1$value));names(nm)<-nm
    cat("\r   sample step 2               \n")
    p<-lapply(nm,function(y){
      z<-pts[pts$JDate %in% y,]
      z<-z[sample(1:nrow(z),maxsample*f1$prop[f1$value %in% y]),]
      z})
    pts<-do.call(rbind,p)
  }else{pts<-xy}
  cat("\r  extracting               \n")
  #y<-crop(y,x)
  mat1<- extract(c(x,y),pts[,1:2])
  mat1<-mat1[apply(mat1,1,function(x){!any(is.na(x))}),]
  nm1<-names(y);names(nm1)<-nm1
  
  outrastl<-pbapply::pblapply(nm1,function(z){
    mat2<-mat1[,c(z,names(x))]
    colnames(mat2)[1]<-"yvar"
    lm1<-lm(yvar~.,data = mat2)
    terraOptions(memfrac = .8)
    outr1<-predict(x,lm1,cores = 2)
    outr2<-subset(y,z)-outr1
    fn<-paste("Y:/MPSG_VegMapping/Data/Raster/Source/afft/Corrected/",cs,"_",z,".tif",sep = "")
    writeRaster(outr2,file = fn,datatype = "INT4S",overwrite = T)
    outr2
  })
  outrastl
}


## Alter this to scroll through states, bands. Consider running bands in parallel on 7 cpus?  
## States in parallel on 3 cpus?
## Workflow --------------------
## Note: Currently wrapping up SD.  Eventually, add 'r' back in, and loop through the remaining states.

states<-c( "CO",   "ND", "NE", "SD")# "WY", "MT", "NM","TX","UT","KS","OK" "seamlines$State[!seamlines$State %in% c("OK")]
#bands<-c("g","n","ndvi","ndgr","ndng","bri")#"r",

for(curstate in states){
  if(curstate %in% c("MT","SD","WY")){bands<-c("g","n","ndvi","ndgr","ndng","bri")
  }else{bands<-c("r","g","n","ndvi","ndgr","ndng","bri")}
  
  cat("\n#####",curstate,"#####\n")
  vrt0<-rast("r/r_vrt.vrt")
  cat("\r making seamline/image date raster layers  ")
  sl<-vect(seamlines[curstate,4])
  sl<-project(sl,vrt0)
  plot(sl)
  e1<-vect(ext(vrt0))
  e1<-buffer(e1,20000)
  sl<-crop(sl,e1)
  plot(sl)
  cat("\r cropping vrt0                            \n ")
  tmp<-crop(vrt0,sl)
  tmp<-trim(tmp)
  sl<-crop(sl,tmp)
  
  sl$JDate<-julian(as.Date(sl$IDATE),origin = as.Date("2019-01-01"))
  #sl2<-crop(sl,buffer(vect(ext(tmp)),20))
  sl2<-as.lines(sl)
  
  cat("\r rasterizing jdate                                   ")
  jdate<-rasterize(sl,tmp,"JDate")
  cat("\r cropping jdate                                   ")
  jdate<-crop(jdate,tmp)
  jdate<-as.factor(jdate)
  tmp<-mask(tmp,jdate)
  cat("\r masking jdate                                   ")
  jdate<-mask(jdate,tmp$r_vrt_1)
  
  cat("\r making seamline distance                                   ")
  sl2<-distance(tmp,sl2,rasterize = T)
  cat("\r masking seamline distance                                   ")
  sl2<-mask(sl2,subset(tmp,1))
  sl2<-trim(sl2)
  jdate<-crop(jdate,sl2)
  tmp<-crop(tmp,sl2)
  
  cat("\n   Sample step 1        ")
  sampdens<-50;maxsample<-100000
  pts<-spatSample(jdate,size = max(sampdens,min(maxsample,ncell(jdate)/sampdens)),xy = T,method = "stratified")
  cat("\r  freq               ")
  f1<-freq(jdate)
  f1$prop<-f1$count/sum(f1$count)
  nm<-as.numeric(as.character(unique(f1$value)));names(nm)<-nm
  cat("\r   sample step 2               \n")
  p<-lapply(nm,function(y){
    z<-pts[pts$JDate %in% y,]
    sampsize<-maxsample*f1$prop[f1$value %in% y]
    if(sampsize > nrow(z)){return(z)
    }else{z<-z[sample(1:nrow(z),sampsize),];return(z)}})
  pts<-do.call(rbind,p)
  
  cat("\r processing band corrections                                   ")
  
  for(band in bands){
    cat("\n#####",curstate,":",band,"#####\n")
    if(which(bands == band)>1){
      vrt1<-rast(paste(band,"/",band,"_vrt.vrt",sep = ""))
      cat("      cropping band                    ")
      tmp<-crop(vrt1,sl2)
      cat("\r      masking band                     ")
      tmp<-mask(tmp,jdate)
    }
    
    cat("\r      correcting band \n")
    cl<-CorrectLayers(tmp,c(jdate,sl2),cs = curstate,xy = pts)
  }
  
}
 

