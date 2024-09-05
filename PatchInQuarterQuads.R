## Get Tile Templates, and footprints ----------
  template<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles/r_tile_044.tif")
  fp<-vect(GetFootprints(list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = T,pattern = "tif")))
  tl<-list.files("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/r/r_tiles",full.names = F)
  tl<-substr(tl,8,10)
  fp$TileIDs<-tl
## Read in files for r only ---------
  rf<-c('C:/TEMP/0_raw/m_3610001_ne_14_060_20210605.tif','C:/TEMP/0_raw/m_3610002_ne_14_060_20210429.tif','C:/TEMP/0_raw/m_3710250_sw_13_060_20210708.tif','C:/TEMP/0_raw/m_3810351_ne_13_060_20210730.tif','C:/TEMP/0_raw/m_3810823_sw_12_060_20210905.tif','C:/TEMP/0_raw/m_3910717_se_13_060_20210812.tif','C:/TEMP/0_raw/m_4010363_nw_13_060_20210717.tif','C:/TEMP/0_raw/m_4010754_sw_13_060_20210828.tif','C:/TEMP/0_raw/m_4110207_sw_13_060_20200707.tif','C:/TEMP/0_raw/m_4110623_sw_13_060_20190825.tif','C:/TEMP/0_raw/m_4210101_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4210456_se_13_060_20200711.tif','C:/TEMP/0_raw/m_4310101_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310103_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310109_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310109_sw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310123_se_14_060_20210917.tif','C:/TEMP/0_raw/m_4310123_sw_14_060_20210829.tif','C:/TEMP/0_raw/m_4310124_nw_14_060_20210829.tif','C:/TEMP/0_raw/m_4310124_se_14_060_20210829.tif','C:/TEMP/0_raw/m_4310124_sw_14_060_20210829.tif','C:/TEMP/0_raw/m_4310125_ne_14_060_20210824.tif','C:/TEMP/0_raw/m_4310125_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310125_se_14_060_20210824.tif','C:/TEMP/0_raw/m_4310125_sw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310126_ne_14_060_20210824.tif','C:/TEMP/0_raw/m_4310126_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310126_se_14_060_20210824.tif','C:/TEMP/0_raw/m_4310126_sw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310127_ne_14_060_20210824.tif','C:/TEMP/0_raw/m_4310127_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310128_ne_14_060_20210824.tif','C:/TEMP/0_raw/m_4310133_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310133_se_14_060_20210824.tif','C:/TEMP/0_raw/m_4310133_sw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310141_ne_14_060_20210824.tif','C:/TEMP/0_raw/m_4310141_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4310141_se_14_060_20210824.tif','C:/TEMP/0_raw/m_4310149_se_14_060_20210824.tif','C:/TEMP/0_raw/m_4310515_se_13_060_20190718.tif','C:/TEMP/0_raw/m_4410057_ne_14_060_20210829.tif','C:/TEMP/0_raw/m_4410149_nw_14_060_20210824.tif','C:/TEMP/0_raw/m_4410458_se_13_060_20190901.tif')
  rf<-gsub("C:/TEMP/0_raw","N:/mpsg_naip_afft/2_aggregated",rf)
  rf1<-gsub("2_aggregated","2_aggregated/r",rf)

  nl1<-list.files("N:/mpsg_naip_afft/2_aggregated/r",full.names = T)
  fpall<-vect(GetFootprints(nl1))  
  fpall$filenames<-gsub("N:/mpsg_naip_afft/2_aggregated/r/","",nl1)
  rf2<-pbapply::pbsapply(rf1,function(x){y<-gsub("N:/mpsg_naip_afft/2_aggregated/r","C:/TEMP/p1",x);r<-rast(x);project(r,template,method = "bilinear",align_only = T,filename = y,overwrite = T);y})
  rf3<-vrt(rf2)
  fpn<-vect(GetFootprints(rf2))
  nm<-gsub("C:/TEMP/p1/","",rf2)
  fpn$File<-nm
  tll<-unique(do.call(c,lapply(1:nrow(fpn),function(x){y<-fpn[x,];y<-crop(fp,y);return(y$TileIDs)})))

  fp<-fp[fp$TileIDs %in% tll,]

## Tile by Tile, replace values.
  q<-pbapply::pblapply(2:nrow(fp),function(x){
    y<-crop(fpn,fp[x,])
    cat("\n\n###",fp$TileIDs[x],"###\n")
    for(b in c("r","g","n","ndvi","ndgr","ndng","bri")){
      cat("---",b,"---\r")
      par(mfrow =c(1,2))
      ct<-rast(paste("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/",b,"/",b,"_tiles/",b,"_tile_",fp$TileIDs[x],".tif",sep = ""))
      fl1<-project(vrt(paste("N:/mpsg_naip_afft/2_aggregated/",b,"/",y$File,sep = "")),ct,method = "near",align_only = T)
      pt1<-plot(ct,"f-20",main = "f-20, unfixed")
      names(fl1)<-names(ct)
      ct1<-merge(fl1,ct)
      ct1<-crop(ct1,ct)
      plot(ct1,"f-20",main = "f-20, fixed")
      
      writeRaster(ct1,sources(ct),overwrite = T)
    }
  })
  
  
  GetAFFT<-function("N:/mpsg_naip_afft/2_aggregated/n/m_4210254_ne_13_060_20200708.tif",
                    zradii =c(1.25,2.084,4.164,8.34,16.4,30,60),outres = 30,#c(0.75, 1.25,2.5, 5, 10, 60)
                    overwrite = T,ncpu = 7,
                    rawpath = "0_raw",
                    indpath = "1_intermediate",
                    aggpath = "2_aggregated")