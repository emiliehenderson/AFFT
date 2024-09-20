
library(terra)
library(SDMap)
setwd("C:/TEMP/ComancheCimarrone/PCA")
fl0<-list.files(,pattern = "AFFT_PCA_CC")
setwd("..")
outdir<-MakeDirIfNeeded("PCA_PATCH",goto = F)
lapply(fl0,function(x){MakeDirIfNeeded(x,outdir,goto = F)})
setwd("PCA")
fll1<-list.files(getwd(),recursive = T, full.names = T,pattern = "tif")
fll1<-fll1[!grepl("pca1_",fll1)]

pbapply::pblapply(fll1,function(x){
  r2<-r1<-rast(x)
  for(i in 1:20)r2<-focal(r2,w = 3,fun = mean,na.policy = "only")
  names(r2)<-names(r1)
  outfile<-gsub("/PCA/","/PCA_PATCH/",sources(r1))
  writeRaster(r2,filename = outfile,overwrite = T)
  return(T)
})

