


MakeTiles<-function(fp1,fl){
  #src1<-sources(rast)
  # sfExport("src1")
  #write("","C:/TEMP/progress.txt",append = F)
      sfClusterApplyLB(names(rast),function(p){
        progfile<-paste("C:/TEMP/",p,"_progress.txt",sep = "")
        write(paste("Starting",p,"at:",Sys.time()),progfile,append = F)
       # z<-subset(rast(src1),p,filename = paste("C:/TEMP/",p,".tif",sep = ""),overwrite = T,memfrac = .9)
        #write(paste("got temporary tif: ",Sys.time()),progfile,append = F)
        z<-rast(paste("C:/TEMP/",p,".tif",sep = ""))
      #   g <- st_read("Y:/MPSG_VegMapping/Data/Spatial/spatial_database.gpkg", layer = "non-overlapping tiles")
      #   pbapply::pblapply(fp1$id,function(tid){
      #     
      #     fp2 <- vect(g[g$id %in% tid,])
      #     cat(tid)
      #     z<-rast("C:/TEMP/temp.tif")
      # 
      #     cat(" cropping");z<-crop(z,fp2,memfrac = .9)
      #     cat(" extending");z<-extend(z,fp2)
      #     SDMap::MakeDirIfNeeded(p,"C:/TEMP/Dakotas/Predictors/")
      #   
      #     cat(" writing \n")
      #     writeRaster(z,filename = paste("C:/TEMP/Dakotas/Predictors/",p,"/",p,"_",tid,".tif",sep = ""),overwrite = T,memfrac = .9)
      #     #file.remove(paste("C:/TEMP/",p,".tif",sep = ""))
      #     write(paste("   ",tid,Sys.time()),progfile,append = T)
      # })
      })
}
library(SDMap)
library(sf)
library(snowfall)


#fl<-list.files("C:/TEMP",pattern = ".tif")
#fl<-fl[!substr(fl,1,1)=="t"]

library(terra);library(AFFT);library(sf)
setwd("C:/TEMP/RockyMountains/PCA")
fl<-list.files(,pattern = "tif")
fl<-fl[!grepl("_",fl)]

#r1<-rast(fl[1])
#names(r1)<-paste("AFFT_PCA_RM_",1:dim(r1)[3],sep = "")
#r2<-writeRaster(r1,"C:/TEMP/temptiles/tmp1.tif",overwrite = T)
#r1<-writeRaster(r2,fl[1],overwrite = T)
v1<-vrt(fl,filename = "PCA.vrt",overwrite = T,set_names = T)

sfInit(T,10)
sfLibrary(terra)

terraOptions(memfrac = .9)
progfile<-"C:/TEMP/Progress.txt"
write("",progfile,append = F)
sfExport("progfile")
sfClusterApplyLB(names(v1),function(x){
  terraOptions(memfrac = .08)
  v1<-rast("PCA.vrt")
  fn<-paste("C:/TEMP/",x,".tif",sep = "")
  cat(x, ": subsetting")
  y<-subset(v1,x)
  cat(",  writing")
  writeRaster(y,filename = fn,overwrite = T)
  write(x,progfile,append = T)
  fn
})


bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
bnd0<- bnd<-project(bnd,rast(list.files("Y:/MPSG_VegMapping/Data/Raster/Predictors/elevation",full.names = T,pattern = "tif")[1]))
bnd<-bnd[bnd$FORESTNUMB %in% c("06","10","15"),]#
bnd<-buffer(bnd,100000)
bnd<-aggregate(bnd)
fl<-list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = T)
names(fl)<-substr(list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = F),8,10)

fp<-GetFootprints(fl,return.as.list = F)
fp$TileID<-names(fl)
crs(fp)<-crs(bnd)
bnd1<-crop(fp,bnd)


g <- st_read("Y:/MPSG_VegMapping/Data/Spatial/spatial_database.gpkg", layer = "non-overlapping tiles")
g2<-vect(g);g2<-crop(g2,bnd)
tids<-g2$id
g <- g[g$id %in% tids,]

nm<-names(v1)#gsub(".tif","",fl)
sfLibrary(sf)
sfExport("tids")
sfClusterApplyLB(nm,function(p,tids1 = tids){
  progfile<-paste("C:/TEMP/",p,"_progress.txt",sep = "")
  write(paste("Starting",p,"at:",Sys.time()),progfile,append = F)
  # z<-subset(rast(src1),p,filename = paste("C:/TEMP/",p,".tif",sep = ""),overwrite = T,memfrac = .9)
  #write(paste("got temporary tif: ",Sys.time()),progfile,append = F)
  z1<-rast(paste("C:/TEMP/",p,".tif",sep = ""))
  g <- st_read("Y:/MPSG_VegMapping/Data/Spatial/spatial_database.gpkg", layer = "non-overlapping tiles")
  fp1<-g[g$id %in% tids1,]
  pbapply::pblapply(fp1$id,function(tid){
    fp2 <- vect(g[g$id %in% tid,])
    cat(tid)
    # z<-rast("C:/TEMP/temp.tif")
    
    cat(" cropping");z<-crop(z1,fp2,memfrac = .9)
    cat(" extending");z<-extend(z,fp2)
    SDMap::MakeDirIfNeeded(p,"C:/TEMP/RockyMountains/Predictors/")
    
    cat(" writing \n")
    writeRaster(z,filename = paste("C:/TEMP/RockyMountains/Predictors/",p,"/",p,"_",tid,".tif",sep = ""),overwrite = T,memfrac = .9)
    #file.remove(paste("C:/TEMP/",p,".tif",sep = ""))
    write(paste("   ",tid,Sys.time()),progfile,append = T)
  })
})







library(sf)


library(snowfall)
sfInit(T,5)
sfLibrary(AFFT)
sfLibrary(sf)
sfLibrary(terra)
MakeTiles(vect(g),v1)
# 
# 
# src<-sources(v1)
# sfExport("src")
# sfExport("tids")
# source("D:/RPackages/AFFT/TileRasterFunctions.R")
# sfExport("create_raster_tiles")
# sfExport("create_tiles_ext")
# 
# 
#   write("starting tile functions","progress.txt",append = F)
#   
# sfClusterApplyLB(names(v1),function(var_name){
#   v1<-rast(src)
#   r1<-subset(v1,var_name)
#   outdir <- paste("C:/TEMP/Dakotas/Predictors", var_name, sep = "/")
#   g <- st_read("Y:/MPSG_VegMapping/Data/Spatial/spatial_database.gpkg", layer = "non-overlapping tiles")
# 
#   g <- g[g$id %in% tids,]
#   create_raster_tiles(r1, g, outdir, variable_name = var_name, prefix = NULL, dtype = "FLT4S",overwrite = T)
#   write(var_name,"progress.txt",append = T)
# })
# 
# sfStop()
# 
# 
# 
# 
#   
#   
# 
# MakeTiles(fp[fp$TileID %in% tids,],v1)
# 
# ## Make a rastersource table.
