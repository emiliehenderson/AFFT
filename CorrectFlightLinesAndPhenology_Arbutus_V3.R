
## Setup -------------------------------------
  setwd("D:/LocalNaip/3_projected/")
  library(readxl)
  library(terra)
  library(AFFT)
  library(snowfall)
  seamlines <- data.frame(read_excel("Y:/MPSG_VegMapping/Data/Raster/Source/afft_full/naip_collection_year_V2.xlsx"))
  rownames(seamlines)<-nm<-seamlines$State;names(nm)<-nm
  fp<-vect("D:/LocalNaip/3_projected/footprints.shp")#GetFootprints(fl<-list.files("r/r_tiles",full.names = T),return.as.list = T,quiet = T)
  
## Buffer Model Region and Select Tiles ------------------------------

  bnd0<-bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
  bnd0<-bnd<-project(bnd,fp)
  bnd<-bnd0[bnd0$FORESTNAME %in% c("Comanche and Cimarron National Grasslands"),]#c("18")
  plot(bnd0);plot(bnd,add = T, col = "skyblue2",border = "skyblue2")
  bnd<-buffer(bnd,100000)
  
  curtiles<-intersect(fp,bnd)  ## Need to drop tile 181
  
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
  samp1<-spatSample(suppl2$jdate,500,method = "stratified",xy = T,na.rm = T)
  samp2<-spatSample(suppl2$jdate,75000,method = "random",xy = T,na.rm = T)
  samp3<-spatSample(suppl2$imgyr,500,method = "stratified",xy = T,na.rm = T)
  samp4<-spatSample(suppl2$stateid,500,method = "stratified",xy = T,na.rm = T)

  
  
  xy<-do.call(rbind,list(samp1[,1:2],samp2[,1:2],samp3[,1:2],samp4[,1:2]))
  print(system.time(mat1<-extract(suppl,xy)))
  
## Band by band -----------------
  bands<-c("r","g","n","ndgr","ndng","bri","ndvi" )
  for(b in bands){
    cat("\n\n#################\n#####  ",b,"  #####\n#################\n\n")
    ### Build vrt for only selected tiles. ---------------
    tl<-unique(curtiles$ImageName)
 
    fll<-paste("D:/LocalNaip/3_projected/",b,"/",tl,sep = "")
    fllt<-paste("C:/TEMP/PreCorrection/",b,"__",tl,".tif",sep = "")
    cat(" Copying files          \n")
    
    pbapply::pblapply(1:length(fll),function(i){file.copy(fll[i],fllt[i])})
    
    cat("building vrt                \n")
    v1<-vrt(fllt);nm<-names(v1)<-names(rast(fllt[1]));names(nm)<-nm
    ### Extract values ---------------------
    cat("extracting                \n")
    #cat(system.time(mat2<-extract(v1,xy[,1:2]))[3]/60,"min.        \r")
    
  
    sfInit(T,30)
    sfExport("xy")
    sfExport("fllt")
    sfLibrary(terra)
    tmp<-sfClusterApplyLB(fllt,function(x){y<-rast(x);z<-extract(y,xy[,1:2]);z<-z[!is.na(z[,2]),];z})
    mat3<-unique(do.call(rbind,tmp))
    mat3<-unique(merge(mat1,mat3,by = "ID"))
    
    mat4<-split(mat3,mat3$ID)
    mat4<-pbapply::pblapply(mat4,function(x){if(nrow(x)>1){rbind(sapply(x,function(x){if(is.factor(x)){return(x[1])}else{return(mean(x))}}))}else{return(x)}})
    mat4<-do.call(rbind,mat4)
    ## Correct Layers -----
    cat("\r Correcting layers                 \n\n")
    terraOptions(memfrac = .8)
    sfInit(T,30)
    sfLibrary(terra)
    s1<-sources(suppl)
    sfExport(list = c("mat4","b","s1"))
    
    write("Correcting","D:/LocalNaip/progress.txt",append = F)
   # if(b == "ndvi"){tl<-tl[19]}
      nb<-1:(dim(v1)[3])
      sfExport(list = c("nb","fllt","tl")) 
    tmp<-sfClusterApplyLB(tl,function(tl1){
      write(paste("Tile",tl1,":",which(tl==tl1),"out of",length(tl)),"D:/LocalNaip/Progress.txt",append = T)
      q<-lapply(nb,function(b2){
        gc()
        terraOptions(memfrac = .8/20)
        mat5<-mat4[,c(2:5,b2+5)]
        colnames(mat5)[5]<-"y"
        
        mat5<-mat5[apply(mat5,1,function(x){!any(is.na(x))}),]
        suppl<-rast(s1)
        
        
        lm1<-lm(y~jdate+sl+ jdate*sl,data = mat5)#jdate*imgyr+imgyr+stateid+
        
        fll2<-list.files("C:/TEMP/PreCorrection",full.names = T)
        fll2<-rast(fll2[grepl(tl1,fll2)])
        fll2<-subset(fll2,b2)
        suppl3<-crop(suppl,fll2)
        #suppl3<-trim(suppl3)
        fll2<-crop(fll2,suppl3)
        fll2<-rast(list(y = fll2,suppl3))
        
        outname<-paste(b,"_",colnames(mat4)[b2+5],sep = "")
        outfile<-paste("C:/TEMP/ComancheCimarrone/Corrected/",outname,"_",tl1,sep = "")
        
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

