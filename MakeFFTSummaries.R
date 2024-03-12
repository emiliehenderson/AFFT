## Setup ------------------
  library(terra);library(AFFT)
  ci<-NULL
  airpath0<-"J:/airphoto/WY_NAIP/WY_NAIP2019_4band/naip2019_wy/geotiff"
  airpath_local<-"0_raw"
  outpath1<-"1_intermediate"
  outpath2<-"2_aggregated"
  outpath3<-"3_corrected"
  outpath4<-"4_PCA"
## Get tile footprints ---------
    fl1<-paste(airpath0,list.files(airpath0),sep = "/")
    if(!all(substr(fl1,nchar(fl1)-2,nchar(fl1))=="tif")){
      #if(readline("explore subfolers?")%in% "y"){
        fl2<-do.call(c,lapply(fl1,function(x){y<-paste(x,list.files(x),sep = "/")
          if(!all(grepl("tif",y))){
            if(length(y)==1){z<-paste(y,list.files(y),sep = "/")
            }else{z<-lapply(y,function(p){paste(p,list.files(p),sep = "/")
              
              })
            }
            return(z)
          }else{return(y)}
        }))
        fl3<-do.call(c,pbapply::pblapply(fl2,function(x){paste(x,list.files(x),sep = "/")}))
        
      #}
    }
    df<-reshape2::colsplit(fl3,"/",c("a","b","c","d","e","f","g","h","i"))
    df<-data.frame(df,filepath = fl3)
    df$Name<-gsub(".tif","",df$i)
    fl4<-df$filepath;names(fl4)<-df$Name
    fl4<-fl4[substr(fl4,nchar(fl4)-2,nchar(fl4)) =="tif"]
    ## note: There's a hidden subfolder or two in here. perhaps it's not one I need?
## Select model region to work on -------------------
    
    if(is.null(ci)){ci<-1:length(fl4)}else{ci<-1:ci}
    fp1<-GetFootprints(fl4[ci],quiet = T)
    save(fp1,file = "1_intermediate/fp1.RData")
    #fp2<-GetFootprints(fl4[ci],crs(rast(fl1[1])),T)
    #names(fp1)<-names(fl4)[ci]
    
    #names(fp2)<-nm1[ci]
    
    statebnd <- vect("J:/R4VegMapping_Development/R4VegMapping/Data/Spatial/cb_2019_us_state_500k.shp")
    WY<-statebnd[statebnd$NAME=="Wyoming",]
    WY<-project(WY,fp1[[1]])
    plot(WY)
    
    bnd<-vect("J:/support/r4_mappingsections_final.shp")
    bnd<-project(bnd,fp1[[1]])
    
    plot(bnd,add = T)
    text(bnd,bnd$ColumnName)
    plot(vect(fp1), border = "gray",add = T)
    plot(WY, lwd = 2,add =T)

    plot(vect(fp1),border = "lightgray")
    plot(bnd,add = T,border = "royalblue",lwd = 2)
    text(bnd,bnd$ColumnName,col = "red",font = 2,cex = 1)
    bnd1<-bnd[bnd$ColumnName %in% "X342G",]
    bnd1<-buffer(bnd1,20000)

## Select tiles to process ----------------
    
    plot(bnd1, border = "darkgray");text(bnd1,bnd1$ColumnName)
    
    plot(vect(fp1), border = "gray",add = T)
    plot(WY, lwd = 2,add =T)
    if(length(fpz)>0)plot(fpz,border = "royalblue",add = T,lwd = 3.5, ext = ext(fpz)*2)
    nm1<-names(fp1);names(nm1)<-nm1
    fp3<-do.call(c,pbapply::pblapply(nm1,function(x,ee = bnd1){
      ey<-fp1[[x]]
      plot(ey,add = T,border = "skyblue2",lwd = 1)
      if(is.related(ey, ee, "intersects")){
        plot(ey,add = T, border = "green",lwd = 2)
        return(x)
      }
    }))

  ## Calculate summaries and metrics over tiles @30m scale ------------------
    done.files1<-gsub(".tif","",list.files(outpath1))
    
    cat("\r tiles to fully process:",sum(!(fp3 %in% done.files1)))
    fp3<-fp3[!fp3 %in% done.files1]
    #save(fp3,file = "1_intermediate/fp3.RData")
    
    tmp<-pbapply::pblapply(fl4[names(fp3)],function(x){CopyImage(x)})
    
    rawfiles<-list.files("0_raw",full.names = F)
    indexfiles<-gsub("_Indices","",list.files("1_intermediate",full.names = F))
    rawfiles<-rawfiles[!rawfiles %in% indexfiles]
    rawfiles<-paste("0_raw",rawfiles,sep = "/")
    indices<-pbapply::pblapply(rawfiles,function(x){GetMetrics1(x)})
    
    
    cfl<-paste(outpath,list.files(outpath),sep = "/")
    cfl<-cfl[grepl(".tif",cfl)]
    done.files<-lapply(lapply(cfl,rast),function(x,y = ID){
      if(crs(x,T)!= crs(y,T)){x<-project(x,y)}
      x
    })

      
    
  ## Correct flight line and phenology artifacts -----------------
    r1<-list.files(outpath2,full.names = T);names(r1)<-list.files(outpath2)
    r1<-r1[grepl("afft",r1)]
    

    afft_init<-MergeAggregatedAirphotos(r1,outpath2 = outpath2)
    afft_init<-rast(list.files(outpath2,full.names = T)[grepl("AFFT_Initial",list.files(outpath2))])
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
    