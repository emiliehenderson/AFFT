
library(terra)
library(SDMap)
library(snowfall)

setwd("C:/TEMP/RockyMountains/Predictors")
fl0<-list.files(,pattern = "AFFT_PCA_")
setwd("..")
outdir<-MakeDirIfNeeded("PCA_PATCH",goto = F)
lapply(fl0,function(x){MakeDirIfNeeded(x,outdir,goto = F)})
setwd("Predictors")
fll1<-list.files(getwd(),recursive = T, full.names = T,pattern = "tif")
fll1<-fll1[!grepl("pca1_",fll1)]

outfiles<-gsub("/Predictors/","/PCA_PATCH/",fll1);names(outfiles)<-outfiles
fll1<-fll1[!file.exists(outfiles)]

  sfInit(T,24)
  sfLibrary(terra)
tmp<-sfClusterApplyLB(fll1,function(x){
    outfile<-gsub("/Predictors/","/PCA_PATCH/",x)
    if(!file.exists(outfile)){
      r2<-r1<-rast(x)
      for(i in 1:30)r2<-focal(r2,w = 3,fun = median,na.policy = "only",na.rm = T)
      names(r2)<-names(r1)
      writeRaster(r2,filename = outfile,overwrite = T)
    }
    gc()
    outfile
})

fl<-list.files("Y:/MPSG_VegMapping/Data/Raster/Predictors",pattern = "_RM_",recursive = T)
