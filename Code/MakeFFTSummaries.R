source("Code/Functions.R")
  
## RunCode ------------------
  library(terra)
  ci<-NULL
  outpath<-"D:/R4_TEMP/AFFT"
  airpath<-"D:/R4_TEMP/AFFT/raw"

## Get tile footprints ---------
    fl1<-paste(airpath,list.files(airpath),sep = "/")
    nm1<-sapply(fl1,function(x){gsub(".tif","",strsplit(x,"/")[[1]][5])})
    names(fl1)<-nm1
    if(is.null(ci)){ci<-1:length(fl1)}else{ci<-1:ci}
    fp1<-GetFootprints(fl1[ci])
    fp2<-GetFootprints(fl1[ci],crs(rast(fl1[1])),T)
    names(fp2)<-nm1[ci]
    
    statebnd <- vect("J:/R4VegMapping_Development/R4VegMapping/Data/Spatial/cb_2019_us_state_500k.shp")
    ID<-statebnd[statebnd$NAME=="Idaho",]
    ID<-project(ID,fp1)
    plot(ID)
    
    bnd<-vect("J:/support/r4_mappingsections_final.shp")
    bnd<-project(bnd,fp1)
    
    plot(bnd)
    text(bnd,bnd$ColumnName)
    plot(fp1, border = "gray",add = T)
    plot(ID, lwd = 2,add =T)

    plot(fp1)
    plot(bnd,add = T)
    text(bnd,bnd$ColumnName,col = "red",font = 2,cex = 3)
    bnd<-bnd[bnd$ColumnName %in% "M332F",]
    bnd1<-buffer(bnd,20000)
  ## Select tiles to process ----------------
    
    {
      windows(10,10)
      plot(bnd1,border = "transparent")
      plot(fp1, border = "gray",add = T)
      plot(ID, lwd = 2,add =T)
      fpz<-list.files(outpath)[grepl("m_",outpath)]
      fpz<-GetFootprints(paste(outpath,fpz,sep = "/"))
      plot(fpz,border = "royalblue",lwd = 3.5, ext = ext(fpz)*2)
      plot(fp1,border = "gray",add = T)
      e1<-draw()
      plot(ID)
      plot(fp1,border = "lightgray",lty = 2,add = T)
      plot(e1,border = "red",add = T)
      readline("pause")
      dev.off()
    }
    plot(bnd, border = "darkgray");text(bnd,bnd$ColumnName)
    plot(fp1, border = "gray",add = T)
    plot(ID, lwd = 2,add =T)
    plot(fpz,border = "royalblue",add = T,lwd = 3.5, ext = ext(fpz)*2)
    
    fp3<-do.call(c,pbapply::pblapply(nm1,function(x,ee = vect(e1)){
      ey<-fp2[[x]]
      plot(ey,add = T,border = "skyblue2",lwd = 1)
      if(is.related(ey, ee, "intersects")){
        plot(ey,add = T, border = "green",lwd = 2)
        return(x)
      }
    }))

  ## Calculate summaries and metrics over tiles @30m scale ------------------
    done.files<-gsub(".tif","",list.files(outpath))
    cat("\r tiles to fully process:",sum(!(fp3 %in% done.files)))
    fp3<-fp3[!fp3 %in% done.files]
    
    imglist<-pbapply::pblapply(fl1[fp3],function(x){GetMetrics(x,outpath,ncpus = 15)})
    cfl<-paste(outpath,list.files(outpath),sep = "/")
    cfl<-cfl[grepl(".tif",cfl)]
    done.files<-lapply(lapply(cfl,rast),function(x,y = ID){
      if(crs(x,T)!= crs(y,T)){x<-project(x,y)}
      x
    })

    r1<-merge(sprc(done.files))
    
  ## Correct flight line and phenology artifacts -----------------
    r1<-list.files(outpath)
    r1<-r1[grepl("afft",r1)]
    r1<-r1[order(as.numeric(substr(r1,5,6)))]
    r1<-rast(paste(outpath,r1,sep = "/"))
    sl<-vect("J:/ID_NAIP/suppl/Seamlines_id16_2019.shp")
    sl<-project(sl,ID)
    sl$JDate<-julian(as.Date(sl$IDATE),origin = as.Date("2019-01-01"))
    sl2<-crop(sl,buffer(vect(ext(r1)),20))
    sl2<-as.lines(sl2)

    jdate<-as.int(rasterize(sl,r1,"JDate"))
    sl2<-distance(r1,sl2,rasterize = T)


    cl<-CorrectLayers(r1,c(jdate,sl2));writeRaster(cl,paste(outpath,"CorrectedAFFT_JdateAsInt.txt",sep = "/"))
  ## Make a PCA from the corrected raster and check results ---------------
    cl<-list.files(outpath);cl<-cl[grepl("CorrectedAFFT",cl)];cl<-rast(paste(outpath,cl,sep = "/"))
    r1<-list.files(outpath);r1<-r1[grepl("afft",r1)];r1<-rast(paste(outpath,r1,sep = "/"))
    
    pca1<- MakePCA(cl)
    pca2<- MakePCA(r1)
    writeRaster(pca2,paste(outpath,"AFFT_PCA_Uncorrected.tif",sep = "/"))
    pca3<-CorrectLayers(subset(pca2,1:18),c(jdate,sl2));writeRaster(pca3,paste(outpath,"AFFT_PCA_CorrectedAfterPCA.tif",sep = "/"))
    jdate2<-!is.na(jdate)
    jdate3<-as.factor(jdate)
    pca4<-CorrectLayers(subset(pca2,1:18),c(jdate2,sl2));writeRaster(pca4,paste(outpath,"AFFT_PCA_CorrectedAfterPCA_FlightLineCorrectionOnly.tif",sep = "/"))
    jdate3<-as.int(jdate3)
    
    pca5<-CorrectLayers(subset(pca2,1:18),c(jdate3,as.int(sl2)));writeRaster(pca5,paste(outpath,"AFFT_PCA_CorrectedAfterPCA_FlightLineCorrectionOnly_bands10-18.tif",sep = "/"))


    par(mfrow =c(1,3))
    
    for(i in -1:floor((dim(pca3)[3]/3)-2)){
      plot.me<-1:3+3*(i<-i+1)
      cls<-lapply(list(pca1,pca2,pca3),function(x){plotRGB(subset(x,plot.me),stretch = "lin")})
      cat("\r",toString(plot.me))
    }
  ## Screen layers ---------------------
    nm<-names(pca1);names(nm)<-nm

  ## Drop layers where the R^2 for flightline/phenology is high.----------------------------------
    #ascores<-ScoreArtifacts(pca1,c(jdate,sl2))
   # lapply(ascores,function(x){
   #   if(readline(min(x$coefficients[,4])) %in% "y"){browser()}
   # }) 
    
  ## Drop layers are entirely noise---------------
    
    snowfall::sfInit(parallel = T, cpus = 15,type = "SOCK")
    snowfall::sfExport("ScoreNoise")
    snowfall::sfExport("MakeDonut")
    
    outfile<-sources(pca1)
    snowfall::sfExport("outfile")
    nm<-names(pca1);names(nm)<-nm
    pcanoise<-do.call(rbind,snowfall::sfLapply(nm,function(v,ofile = outfile){
      nl1<-terra::rast(ofile)
      cm<-terra::subset(nl1,v)
      ScoreNoise(cm,q=T, return.rast = F)
    }))
   
    snowfall::sfStop()
    
    goodvars<-subset(pca1,pcanoise[,4]<325)
    badvars<-subset(pca1,pcanoise[,4]>=730)
    checkvars<-subset(pca1,!names(pca1) %in% c(names(goodvars),names(badvars)))
    
    check2<-ViewStackRGB(checkvars)
    check2<-lapply(split(check2,check2),names)
    check2<-lapply(check2,function(x){ out<-do.call(c,reshape2::colsplit(x,"_",c("one","two","three")));names(out)<-NULL;out})
    
    goodvars<-c(goodvars,subset(pca1,check2[[1]]))
    par(mfrow =c(1,3))
    check3<-sapply(check2[[2]],function(x){y<-subset(pca1,x);plot(y);plotRGB(c(y,y,y),stretch = "lin");plotRGB(c(y,y,y),stretch = "hist");readline(x)})
    goodvars<-c(goodvars,subset(pca1,names(check3)[check3 %in% 1]))
    
    names(goodvars)<-gsub("Comp.","",names(goodvars))
    goodvars<-subset(goodvars,names(goodvars)[order(as.numeric(names(goodvars)))])
    writeRaster(goodvars,paste(outpath,"AFFT_PCA_FINAL_SELECTION.tif",sep ="/"))
    