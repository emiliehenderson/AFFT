
## Setup -------------------------------------
  setwd("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3")
  library(readxl)
  library(terra)
  library(AFFT)
  library(snowfall)
  seamlines <- data.frame(read_excel("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full/naip_collection_year_V2.xlsx"))
  rownames(seamlines)<-nm<-seamlines$State;names(nm)<-nm
  fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  fp<-vect(fp)
  fp$TileID<- substr(fl,18,20)
 
## Buffer Model Region and Select Tiles ------------------------------

  bnd0<-bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  fp<-GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  fp<-vect(fp)
  bnd0<-bnd<-project(bnd,fp)
  bnd<-bnd0[bnd0$FORESTNUMB %in% c("07","18"),]#c("18")
  plot(bnd0);plot(bnd,add = T, col = "skyblue2",border = "skyblue2")
  bnd<-buffer(bnd,100000)
  
  fp$TileID<- substr(fl,18,20)
  curtiles<-intersect(fp,bnd)  ## Need to drop tile 181
  curtiles<-curtiles[!curtiles$TileID %in% c(181,213),]
  
  cat("Getting Supplemental layers             \r")
  
  suppl<-rast("D:/LocalNaip/Suppl/suppl.vrt")

  
  suppl<-crop(suppl,curtiles)
  suppl<-mask(suppl,curtiles)
  nlcdmask<-rast("Y:/MPSG_VegMapping/Data/Raster/Source/nlcd/nlcd_2021_crop.tif")
  nlcdmask<-crop(nlcdmask,suppl)
  nlcdmask<-mask(nlcdmask,suppl)
  nlcdmask<-subset(nlcdmask,1)
  nlcdmask2<-nlcdmask %in% c("Unclassified","Barren Land","Deciduous Forest","Evergreen Forest","Mixed Forest","Shrub/Scrub","Herbaceous","Hay/Pasture","Woody Wetlands","Emergent Herbaceous Wetlands")
  suppl2<-mask(suppl,nlcdmask2,maskvalue = F)
  suppl$imgyr<-as.factor(suppl$imgyr)
  suppl$stateid<-as.factor(suppl$stateid)
  suppl<-writeRaster(suppl,"C:/TEMP/tmp.tif",overwrite = T)
  
## Sample Area of Interest ----------
  cat("Sampling AOI               \r")
  samp1<-spatSample(suppl2$jdate,650,method = "stratified",xy = T,na.rm = T)
  samp2<-spatSample(suppl2$jdate,75000,method = "random",xy = T,na.rm = T)
  samp3<-spatSample(suppl2$imgyr,500,method = "stratified",xy = T,na.rm = T)
  tmp<-mask(suppl2$imgyr,suppl2$imgyr == 20,maskvalue = 0)
  tmp<-trim(tmp)
  
  samp4<-spatSample(tmp,250,method = "stratified",xy = T,na.rm = T)
  
  
  xy<-do.call(rbind,list(samp1[,1:2],samp2[,1:2],samp3[,1:2],samp4[,1:2]))
  print(system.time(mat1<-extract(suppl,xy)))
  
## Band by band -----------------
  bands<-c("bri")#"ndvi" "r","g","n","ndgr","ndng",
  for(b in bands){
    cat("\n\n#################\n#####  ",b,"  #####\n#################\n\n")
    ### Build vrt for only selected tiles. ---------------
    tl<-unique(curtiles$TileID)
    browser()
    fll<-paste("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/",b,"/",b,"_tiles/",b,"_tile_",tl,".tif",sep = "")
    fllt<-paste("C:/TEMP/PreCorrection/",b,"_tile_",tl,".tif",sep = "")
    cat(" Copying files          \n")
    
    # fix.fp<-paste(gsub("0","",substr(fll,1,nchar(fll)-4)[1:23]),".tif",sep = "")
    # fix.fp[14]<-"Y:/MPSG_VegMapping/Data/Raster/Source/afft_full3/bri/bri_tiles/bri_tile_80.tif" 
    # fll[1:23]<-fix.fp
    # 
    pbapply::pblapply(1:length(fll),function(i){file.copy(fll[i],fllt[i])})
    
    cat("building vrt                \n")
    v1<-vrt(fllt);nm<-names(v1)<-names(rast(fllt[1]));names(nm)<-nm
    ### Extract values ---------------------
    cat("extracting                \n")
    #cat(system.time(mat2<-extract(v1,xy[,1:2]))[3]/60,"min.        \r")
    
    sfInit(T,10)
    sfExport("xy")
    sfExport("fllt")
    sfLibrary(terra)
    sfExport("nm")
    write("extracting","D:/LocalNaip/progress.txt",append = F)
    mat3<-do.call(cbind,sfLapply(nm,function(x){
      v<-vrt(fllt);names(v)<-names(rast(fllt[1]))
      v<-subset(v,x)
      e1<-extract(v,xy[,1:2])
      write(paste("extracted",x,":",which(nm == x),"out of",length(nm)),"D:/LocalNaip/progress.txt",append = T)
      e1[,x]
    }))
    sfStop()
    mat3<-cbind(mat1,mat3)
    
    ## Correct Layers -----
    cat("\r Correcting layers                 \n\n")
    terraOptions(memfrac = .8)
    sfInit(T,19)
    sfLibrary(terra)
    s1<-sources(suppl)
    sfExport(list = c("mat3","b","s1"))
    
    write("Correcting","D:/LocalNaip/progress.txt",append = F)
   # if(b == "ndvi"){tl<-tl[19]}
    tmp<-pbapply::pblapply(tl,function(tl1){
      cat(tl1)
      nb<-1:(dim(v1)[3])
      sfExport(list = c("tl1","nb")) 
      write(paste("Tile",tl1,":",which(tl==tl1),"out of",length(tl)),"D:/LocalNaip/Progress.txt",append = T)
      q<-sfLapply(nb,function(b2){
        gc()
        terraOptions(memfrac = .8/20)
        mat4<-mat3[,c(2:5,b2+5)]
        colnames(mat4)[5]<-"y"
        
        mat4<-mat4[apply(mat4,1,function(x){!any(is.na(x))}),]
        suppl<-rast(s1)
        
        
        lm1<-lm(y~jdate+sl+imgyr+stateid+jdate*imgyr+ jdate*sl,data = mat4)
        
        fll2<-list.files("C:/TEMP/PreCorrection",pattern = ".tif",full.names = T)
        fll2<-fll2[grepl(paste(b,"_tile",sep = ""),fll2)]
        fll2<-fll2[grepl(tl1,fll2)]
        fll2<-rast(fll2[grepl(tl1,fll2)])
        fll2<-subset(fll2,b2)
        suppl3<-crop(suppl,fll2)
        #suppl3<-trim(suppl3)
        fll2<-crop(fll2,suppl3)
        fll2<-rast(list(y = fll2,suppl3))
        
        outname<-paste(b,"_",colnames(mat3)[b2+5],sep = "")
        outfile<-paste("C:/TEMP/Dakotas/Corrected/",outname,"_tile_",tl1,".tif",sep = "")
        
        names(fll2)[1]<-"y"
        pred1<-predict(fll2,lm1,na.rm = T)
        pred2<-fll2$y - pred1
        names(pred2)<-outname
        
        pred2<-writeRaster(pred2,filename = outfile,overwrite = T)#
        gc()
        sources(pred2)
      })
    })
    rm("v1")
    file.remove(fllt)
    sfStop()
  }  
    
#### Unify tiles to a vrt -----------------------

