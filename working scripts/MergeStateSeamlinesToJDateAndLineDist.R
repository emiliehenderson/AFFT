setwd("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full")
library(terra)
library(readxl)
library(AFFT)

tilemap<-vect("r/r_process_files/r_tile_extents.shp")

seamlines <- data.frame(read_excel("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full/naip_collection_year_V2.xlsx"))

rownames(seamlines)<-nm<-seamlines$State;names(nm)<-nm

setwd("..");setwd("..")

vrt0<-subset(rast("r/r_vrt.vrt"),1)

fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)

fp<-vect(fp)
fp$TileID<- substr(fl,18,20)

cat("\r making seamline/image date raster layers  ")

maplist<-pbapply::pblapply(seamlines$State[seamlines$State !="CO"],function(curstate){
  cat("\n",curstate)
  cat("\n making seamlines vect               ")
  sl<-vect(seamlines[curstate,4])
  sl<-project(sl,vrt0)
  curtiles<-intersect(fp,sl)
  tl<-paste("r/r_tiles/r_tile_",unique(curtiles$TileID),".tif",sep = "")
  names(tl)<-unique(curtiles$TileID)
  tl2<-paste("D:/LocalNaip/TEMP/tile_",names(tl),".tif",sep = "")
  
  vrt0<-vrt(tl2)
  
  sl<-crop(sl,ext(vrt0)*1.01)
  cat("\r cropping vrt                      ")
  tmp<-crop(vrt0,sl)
  slx<-crop(sl,tmp)
  tmp<-mask(tmp,slx)
  tmp<-trim(tmp)
  slx<-crop(sl,tmp)
  
  cat("\r making jdate                      ")
  
  odate<-paste(substr(sl$IDATE[1],1,4),"01","01",sep = "-")
  sl$JDATE<-julian(as.Date(sl$IDATE),origin = as.Date(odate))
  jdate<-rasterize(sl,tmp,"JDATE")

  cat("\r making seamline dist               ")

  e2<-ext(sl) *.999
  sl<-crop(as.lines(sl),e2)
  sl<-distance(tmp,sl,rasterize = T)
  slx$imgyr<-as.numeric(substr(odate,3,4))
  slx$stateid<-which(seamlines$State == curstate)
  imgyr<-rasterize(slx,tmp,"imgyr")
  stateid<-rasterize(slx,tmp,"stateid")
  cat("\r combining and masking output maps   ")
  twomaps<-rast(list(jdate = jdate,sl=sl,imgyr = imgyr,stateid = stateid))
  twomaps<-mask(twomaps,slx,filename = paste("D:/LocalNaip/Suppl/",curstate,"_Supplemental.tif",sep = ""),overwrite = T)
  twomaps
})
setwd("D:/LocalNaip/Suppl")
v1<-vrt(fl<-list.files(,".tif"),filename = "suppl.vrt",overwrite = T,set_names = T)
