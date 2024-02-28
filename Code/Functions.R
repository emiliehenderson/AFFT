
## Functions -------------------

GetFootprints<-function(imglist,outcrs = NULL,return.as.list = F, quiet = T){
  if(is.character(imglist)){
    if(is.null(outcrs)){outcrs<-crs(rast(imglist[1]))}
    el<-pbapply::pblapply(imglist,function(x){if(!quiet){cat("\r",x)};r<-terra::rast(x);e<-terra::vect(terra::ext(r),crs = crs(r));project(e,outcrs)})
  }else if ("SpatRaster" %in% class(imglist[[1]])){
    if(is.null(outcrs)){outcrs<-crs(imglist[1])}
    el<-pbapply::pblapply(imglist,function(r){e<-vect(ext(r),crs = crs(r));e;project(e,outcrs)})
  }
  if(!return.as.list)el<-vect(el)
  
  el
}

GetMetrics<-function(rasterfile,outpath,outres = 30,ncpus = 15,zradii = c(3,6,10,20,40,55)){
  fn<-strsplit(rasterfile,split = "/")[[1]];fn<-fn[length(fn)]
  outfile<-paste(outpath,paste("afft_",fn,sep = ""),sep = "/")
  if(file.exists(outfile)){
    outrast<-rast(outfile)
    res1<-res(rast(rasterfile))[1]
    fact1<-outres/res1
    donut<-round(c(MakeDonut(zradii,res1,fact1)),2)
  }else{
    cat("  \n  Elapsed Time 1:",system.time({
      cat("\n",rasterfile)
      file.remove(paste(outpath,"/temp.tif",sep = ""))
      file.copy(rasterfile,paste(outpath,"/temp.tif",sep = ""))
      r1<-rast(paste(outpath,"/temp.tif",sep = ""));names(r1)<-c("r","g","b","n")
      cat("  making NDVI")
        num<-subset(r1,4)- subset(r1,1)
        denom<-subset(r1,4)+ subset(r1,1)
        ndvi<-num/denom
      cat("  making NDGR")
        num<-subset(r1,2)- subset(r1,1)
        denom<-subset(r1,2)+ subset(r1,1)
        ndgr<-num/denom
      cat("  making NDNG")
        num<-subset(r1,4)- subset(r1,2)
        denom<-subset(r1,4)+ subset(r1,2)
        ndng<-num/denom
      
      cat("  stretching Indices")
        nd<-c(ndvi,ndgr,ndng)
        nd<-stretch(nd,smin = -1, smax = 1)
        
      cat("  making brightness")
       br<-sum(r1)#app(r1,sum)
       br<-stretch(br, smin = 0, smax = 1020)
      #r1<-c(r1,ndvi= stretch(ndvi,smin = -1, smax = 1),ndgr=stretch(ndgr,smin = -1, smax = 1),ndng = stretch(ndng,smin = -1, smax = 1),br<-stretch(br, smin = 0, smax = 1020))
      r1<-c(r1,nd,br)
      names(r1)<-c("r","g","n","ndvi","ndgr","ndng","bri")# "b",
      r1
      #        r1<-writeRaster(r1,filename = paste(outpath,"/temp.tif",sep = ""),overwrite = T)
    })[[3]]/60)
    fact1<-outres/res(r1)[1]
    donut<-round(c(MakeDonut(zradii,res1 = res(r1)[1],fact1)),2)
    cat("  \n  Elapsed Time 2:",system.time({
      outrast<-aggregate(r1,fact = fact1,cores = ncpus,filename = outfile,zm1 = donut,
                         fun = function(x,zm1 = donut){
                           if(!any(is.na(x))){
                             x<-matrix(x,nrow = ceiling(sqrt(length(x))), byrow = T)
                             ff1<-c(abs(gsignal::fftshift(fft(x - mean(x)),MARGIN =c(1,2))^2))
                             y<-tapply(ff1,zm1,sum)/sum(ff1)
                             y2<-c(mean(x),
                                   quantile(x,c(0,.1,.5,.9,1)),
                                   sd(x)*10,
                                   (moments::skewness(c(x))+10)*10,
                                   log((moments::kurtosis(c(x))))*100)
                             outvec<-round(c(y*100,y2),0)
                           }else{
                             outvec<-rep(NA,(length(unique(zm1))+9))
                           }
                           #print(length(outvec))
                           return(outvec)
                         })
      
    })[[3]]/60,"\n")
    beepr::beep(10)
  }
  names1<-paste("f",xfun::numbers_to_words(ceiling(unique(donut))),letters[1:length(unique(donut))],sep = "_")
  names2<-c("mean","min","Q10","med","Q90","max","sd10","skew_p10t10","kurt_lt100")
  
  basenames<- c(names1,names2)
  nms<-   paste(c("r","g","b","n","ndvi","bri"),rep(basenames,each = 6),sep = "_")
  names(outrast)<-nms
  file.remove(paste(outpath,"/temp.tif",sep = ""))
  return(outrast)
}  

MakeDonut<-function(zr,res1,fact1){
  d<-do.call(c,lapply(zr,function(r,res2 = res1,f2= fact1){
    size <- f2*res2
    e = c(-size/2,size/2,-size/2,size/2)
    pt1<-terra::vect(rbind(c(0,0)))
    ya<-terra::rast(nrows=size/res2, ncols=size/res2, nlyrs=1,  extent = terra::ext(e), 
                    resolution = res2)
    y<-terra::rasterize(terra::buffer(pt1, r),ya)
    y<-y * (1/r) * size/2
  }))
  d<-terra::app(d,function(x){y<-max(x,na.rm = T)})
  d<-matrix(c(d[]),nrow = fact1)
  d
}

MakePCA<-function(r,sampdens = 50,maxsample = 100000){
  samp<-spatSample(r,size = min(maxsample,ncell(r)/sampdens))
  samp<-samp[apply(samp,1,function(x){!any(is.na(x))}),]
  cat(paste("sampsize:",nrow(samp)))
  pca1<-princomp(samp)
  pcamap<-predict(r,pca1)
  pcamap
}

CorrectLayers<-function(y,x,maxsample= 100000 ,sampdens = 50){
  names(x)<-c("x1","x2")
  cat("masking\r")
  x<-mask(x,subset(y,1))
  cat("getting sample\r")
  mat1<-spatSample(c(x,y),size = min(maxsample,ncell(y)/sampdens))
  mat1<-mat1[apply(mat1,1,function(x){!any(is.na(x))}),]
  nm1<-names(y);names(nm1)<-nm1
  cat("correcting",length(nm1),"layers\n")
  outrastl<-pbapply::pblapply(nm1,function(z){cat(" ",z,"   \r")
    mat2<-mat1[,c(z,names(x))]
    colnames(mat2)[1]<-"yvar"
    lm1<-lm(yvar~x1+x2+x1:x2,data = mat2)
    outr1<-predict(x,lm1)
    outr2<-subset(y,z)-outr1
    outr2
  })
  rast(outrastl)
}

ScoreNoise<-function(y,f1 = 256,f2 = 64, q = T,
                     scorefun = function(x, q = c(.05,.1,.5,.9,.95)){quantile(x,q,na.rm = T)},
                     return.rast = F){
  ## power of two block size may run faster.
  ## odd-sized images are 'neat'?
  ## Considered: use fft on histogram-equalized stretch? - doesn't make a huge difference
  ## Note: edge effects -- impacts on high frequency detections (function wraps image). Consider adjusting later?
  ## Note: sharp boundaries in original image will contribute to high frequencies detected.
  ## Note: feeding this to a cluster layer-by-layer is faster than using a multicore call to aggregate function.#, cores = 1 

  if(is.character(y)){y<-terra::rast(y)}
  d1<-terra::rast(MakeDonut(zr =c(f1/10,f1)*terra::res(y)[1],res1 = terra::res(y)[1],fact1 = f1))## 2-pixel radius.
      
  testfun<-function(x,fact1 = f1,d=d1, quiet = q){
    x<-matrix(x,nrow = fact1, byrow = T)
    if(!quiet)par(mfrow =c(1,1))
    if(!any(is.na(x))){
      #yp <- x - mean(x) ## subtract off mean values (to normalize ff)
      #ff3<-fft(yp) ## fft of mean-subtracted image
      #ya<-abs(ff3^2)# ## make the complex numbers into real numbers.
      
      ya<-abs(gsignal::fftshift(fft(x - mean(x)),MARGIN =c(1,2))^2)
      z1<-terra::zonal(terra::rast(ya),d,fun = "sum")
      noisescore<-(z1[,2]/sum(z1[,2]))[1]
      if(!quiet){
        par(mfrow =c(3,1))
        terra::plot(terra::rast(x),main = noisescore * 1000)
        terra::plot(terra::rast(ya))
        terra::plot(d)
        st<-Sys.time();while(abs(difftime(st,Sys.time()))<1){}
      }
      return(as.integer(noisescore * 1000))
    }else{return(NA)}
  }
  #if(q){core1<-6}else{core1<-1}
  z<-terra::aggregate(y,fact = f1, testfun,cores = 1)
  if(!q){
    par(mfrow =c(2,1))
    plot(y,main = "original image")
    plot(z)
  }
  tmp<-terra::extract(z,1:terra::ncell(z))
  if(dim(y)[3]==1){
    if(return.rast)return(list(scores = scorefun(tmp),zones = z))
    else return(scorefun(tmp))
  }else{
     if(return.rast)return(list(scores = apply(tmp,2,scorefun),zones = z))
     else return(apply(tmp,2,scorefun))
  }
}

ScoreArtifacts<-function(y,x,maxsample= 100000 ,sampdens = 50){
  mat1<-spatSample(c(x,y),size = min(maxsample,ncell(y)/sampdens))
  mat1<-mat1[apply(mat1,1,function(x){!any(is.na(x))}),]
  nm1<-names(y);names(nm1)<-nm1
  lapply(nm1,function(z){
    mat2<-mat1[,c(z,names(x))]
    colnames(mat2)[1]<-"yvar"
    lm1<-summary(lm(yvar~.,data = mat2))
    plot(y,z,main = format(lm1$r.squared,digits = 8,scientific = F))
    print(lm1)
    readline(z)
    lm1
    # outr1<-predict(x,lm1)
    #  outr2<-subset(y,z)-outr1
    # outr2
  })
  
}

ViewStackRGB<-function(x){
  vec<-sort(rep(1:(ceiling(dim(x)[3]/3)),length.out = dim(x)[3]))
  vec<-split(names(x),vec);names(vec)<-sapply(vec,function(z){paste(z,collapse = "_")})
  
  sapply(vec,function(y){
    tsy<-paste(y,collapse = "_")
    if(length(y)==3){
      par(mfrow =c(1,2))
      plotRGB(subset(x,y),stretch = "lin", main = "Linear Stretch")
      plotRGB(subset(x,y),stretch = "hist",main = "Histogram-equalize\nStretch")
    }else{
      par(mfrow =c(1,length(y)))
      plot(subset(x,y),main = tsy)
    }
    readline(tsy)
  })
}
