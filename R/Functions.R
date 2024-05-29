## GetFootprints ------------
#' Grabs extents from all airphotos imglist.
#'
#' @description This function is useful for imputation mapping from models built with the yaImpute package.
#'
#' @export
#'
#' @param imglist list of airphotos (text, full paths, with names)
#' @param outcrs projection object if reprojection is needed
#' @param return.as.list sets return to one airphoto extent (terra vect format) per object in a list (as opposed to one vect object)
#' @param quiet set to F if more intermediate messaging is wanted.
#' @return Terra vect object containing airphoto extents, either as a single object, or as a list of extents.
#'
GetFootprints<-function(imglist,outcrs = NULL,return.as.list = T, quiet = F){
  if(is.character(imglist)){
    if(is.null(outcrs)){outcrs<-crs(rast(imglist[1]))}
    el<-pbapply::pblapply(imglist,function(x){if(!quiet){cat("\r",x,"           ")};r<-terra::rast(x);e<-terra::vect(terra::ext(r),crs = crs(r));project(e,outcrs)})
  }else if ("SpatRaster" %in% class(imglist[[1]])){
    if(is.null(outcrs)){outcrs<-crs(imglist[1])}
    el<-pbapply::pblapply(imglist,function(r){cat(r,"     \r");e<-vect(ext(r),crs = crs(r));e;project(e,outcrs)})
  }
  if(!return.as.list)vect(el)
  el
}
## CopyImages ---------------
#' Copies raw airphotos to 0_raw local directory for further processing.
#' @description Copies raw airphotos to 0_raw local directory for further processing.
#' @export
#' @param x full path to images (full path).
#' @param outdir local directory to store image copy
#' @return NULL

CopyImage<-function(x,outdir = "0_raw"){
  fn<-strsplit(x,split = "/")[[1]];fn<-fn[length(fn)]
  outfile<-paste(outdir,fn,sep = "/")
  if(!exists(outfile))file.copy(x,outfile)
  return(NULL)
}
## GetBandIndices ------------
#' Copies source airphotos to intermediate directory, calculates extra indices to add to image. Resolution unchanged.
#'
#' @description Copies source airphotos to intermediate directory, calculates extra indices to add to image. Resolution unchanged.
#'
#' @export
#'
#' @param filelist character vector with tif image names. If full path not included, GetMetrics will assume that files are in 0_raw folder.
#' @param indfuns file path pointer where R will store intermediate (full resolution) image tiles, with calculated indices as well as original bands (although not blue currently)
#' @return filepaths pointing to single-band rasters nested in subfolders inside 1_intermediate. 
#'
GetBandIndices<-function(filelist,indfuns = indexFuns,outpath = "1_intermediate",ncpu = 4){
  if(ncpu >1){
    require(snowfall)
    sfInit(parallel = T, cpus = ncpu)
     sfLibrary(terra)
     sfLibrary(AFFT)
     sfExport("indfuns")
     sfExport("outpath")
    # sfExport("filelist")
    # stime<-Sys.time()
    # sfExport("stime")
    # write(" ","donefiles.txt",append = F)
    indices<-pbapply::pblapply(filelist,function(y){
      GetMetrics1(y,outpath,indfuns,parallel = T)
      # write(y,"donefiles.txt",append = T)
      # donefiles<-readLines("donefiles.txt")
      # donefiles<-donefiles[2:length(donefiles)]
      # te<-difftime(Sys.time(), stime,units = "secs")
      # prop.done<-length(donefiles)/length(filelist)
      # tr<-te*(1-prop.done)/(prop.done)
      # hms<-SDMap::SecsToHMS(tr)
      # write(paste(round(prop.done * 100,3),"% done"),"progress.txt",append = F)
      # write("","progress.txt",append = T)
      # hms<-SDMap::SecsToHMS(te)
      # 
      # write("Time Elapsed:","progress.txt",append = T)
      # write(paste(paste(hms,names(hms)),collapse = "   "),"progress.txt",append = T)
      # hms<-SDMap::SecsToHMS(tr)
      # write("Time Remaining:","progress.txt",append = T)
      # write(paste(paste(hms,names(hms)),collapse = "   "),"progress.txt",append = T)
      # return(y)
      })#,parallel = T
    return(indices)
  }else{
    indices<-pbapply::pblapply(filelist,function(y){GetMetrics1(y,outpath,indfuns,parallel = F)})
    return(indices)
  }
 
}
## GetMetrics1 ------------
#' Calculates indices from RBGN. Resolution unchanged. internal for GetBandIndices
#'
#' @description Copies source airphotos to intermediate directory, calculates extra indices to add to image. Resolution unchanged.
#'
#'
#' @param rasterfile character vector with tif image names. If full path not included, GetMetrics will assume that files are in 0_raw folder.
#' @param outpath file path pointer where R will store intermediate (full resolution) image tiles, with calculated indices as well as original bands (although not blue currently)
#' @return multiband raster with ndvi, ndgr, ndng, brightness.
#'
GetMetrics1<-function(rasterfile,
                      outpath = "1_intermediate",
                      fl = indexFuns,parallel = F){
  require(terra)
  cat("\n",rasterfile,"\n")
  
  nm<-names(fl);names(nm)<-nm
  if(parallel){
     sfExport("rasterfile")
        gc()
        indlist<-sfLapply(nm,function(myfun,rf = rasterfile){fl[[myfun]](rf,outpath)})
        gc()
  }else{
    gc()
    indlist<-lapply(nm,function(myfun,rf = rasterfile){fl[[myfun]](rf,outpath)})
    gc()
  }                               
      indlist
  
  
}  



## GetFullResImage ------------
#' internal.
#'
#' @description internal.
#'
#'
#'
GetFullResImageList<-function(fn,subset1 = c("r","g","n","ndvi","ndng","ndgr","bri"),rawpath = "0_raw",indpath = "1_intermediate",indfolders =c("ndvi","ndng","ndgr","bri")){
  fn<-strsplit(fn,"/")[[1]];fn<-fn[length(fn)]
  raw<-data.frame(Path = paste(rawpath,fn,sep = "/"),Band =c(1:4));rownames(raw)<-c("r","g","b","n")
  ind<-data.frame(Path = paste(indpath,indfolders,fn,sep = "/"),Band = 1);rownames(ind)<-c("ndvi","ndng","ndgr","bri")
  df<-rbind(raw,ind)
  df[subset1,c("Path","Band")]
}
## GetAFFT ------------
#' Generate metrics aggregating information from full-res images to coarser resolution.
#'
#' @description Generate metrics aggregating information from full-res images to coarser resolution.
#'
#' @export
#'
#' @param filelist character vector  containing file names for quarter quad images (no path included here. that gets added later)
#' @param zradii radii for summarizing fft-related statistics.
#' @param outres resolution of output raster
#' @param overwrite logical flag -- overwrite existing images?
#' @param ncpu if > 1, parallel processing will be used internally. Most efficient if set to the number of layers that need processing (e.g., for current config., it's 7: r,g,n,ndvi,ndng,ndgr, and bri).
#' @param rawpath path to raw, native resolution naip images
#' @param indpath path to native resolution indices calculated from raw naip images
#' @param aggpath path to save output files aggregated to 'outres' resolution
#' @return raster containing multiple bands, describing texture at scales described over scales indicated by zradii, as well as an array of non-textural summary statistics.
#'
GetAFFT<-function(filelist,
                  zradii =c(0.75, 1.25,2.5, 5, 10, 60),outres = 30,
                  overwrite = T,ncpu = 7,
                  rawpath = "0_raw",
                  indpath = "1_intermediate",
                  aggpath = "2_aggregated"){
  fn<-strsplit(filelist[1],"/")[[1]];fn<-fn[length(fn)]
  r0<-rast(paste(rawpath,fn,sep = "/"))
  res1<-res(r0)[1]
  fact1<-outres/res1
 # donut<-round(c(MakeDonut(zradii,res1,fact1)),2)
   donut<-round(MakeDonut(zradii,res1,fact1,return.rast = F),2)
  
  aggfun1<-function(x,zm1 = donut,na.rm = T){
    if(!any(is.na(x)) & length(x)== length(donut)){
      x<-matrix(x,nrow = ceiling(sqrt(length(x))), byrow = T)
      ff1<-c(abs(gsignal::fftshift(fft(x - mean(x)),MARGIN =c(1,2))^2))
      #y<-scales::rescale(tapply(ff1,zm1,sum)/sum(ff1),from =c(0,1),to =c(0,254))
      y<-tapply(ff1,zm1,mean)
      yb<-scales::rescale(tapply(ff1,zm1,sum)/sum(ff1),from =c(0,1),to =c(0,254))
      y2<-c(mean(x),quantile(x,c(.025,.5,.95)))
      y3<-c(sd(x))
      y4<-log(moments::skewness(c(x))+10) * 50
      y5<-log(moments::kurtosis(c(x)))*100
      outvec<-round(c(log(y)*100,yb,y2,y3,y4,y5),0)
    }else if( length(x) != length(donut) & !any(is.na(x))){browser()
    }else{
      outvec<-rep(NA,(length(unique(zm1))*2+7))
    }
    outvec
  }
  if(ncpu > 1){
    require(snowfall)
    sfInit(ncpu,parallel = T)
    sfLibrary(gsignal)
    sfLibrary(terra)
    sfLibrary(AFFT)
    sfExport("fact1",local = T)
    sfExport("donut",local = T)
    sfExport("rawpath",local = T)
    sfExport("indpath",local = T)
    sfExport("aggpath",local = T)
    sfExport("aggfun1",local = T)
  }
  affts<-pbapply::pblapply(filelist,function(y,fl){GetAFFT1(y,zradii,outres,fact1,donut,overwrite,ncpu,rawpath,indpath,aggpath,aggfun1)})
  
  if(ncpu > 1){sfStop()}
  affts
}

GetAFFT1<-function(r,zradii =c(2,6,56),outres = 30,fact1,donut ,overwrite = T,ncpu,
                   rp = rawpath,ip = indpath,ap = aggpath,aggfun = aggfun1){
  rl<-GetFullResImageList(r,rawpath = rp, indpath = ip)
  nm<-rownames(rl);names(nm)<-nm
  if(ncpu > 1){
    sfExport("rl")
    ol<-sfLapply(nm,function(x,f1=fact1){
      terraOptions(memfrac = 1/7.01,datatype = "INT4S")
      r1<-rast(rl[x,1])
      r1<-terra::subset(r1,rl[x,2])
      outrast<-aggregate(r1,fact = f1,
                    fun = aggfun,zm1=donut)
      nm<-c(paste("f-",unique(round(donut,2)),sep = ""),paste("fp-",unique(round(donut,2)),sep = ""),c("mean","Q025","med","Q95","sd","skew","kurt"))
      names(outrast)<-nm#c(paste("f-",unique(round(donut,2)),sep = ""),paste("fp-",unique(round(donut,2)),sep = ""),c("mean","Q025","Med","Q95","sd","skew","kurt"))
      rn<-strsplit(r,"/")[[1]];rn<-rn[length(rn)]
      writeRaster(outrast,filename = paste(ap,x,rn,sep = "/"),datatype = "INT4S",overwrite = overwrite)
      rm(list =c("outrast","r1"))
      gc()
      return(x)
  })
  
  }else{
    
    terraOptions(memfrac = .9,datatype = "INT4S")
    ol<-lapply(nm,function(x,f1=fact1){
      r1<-rast(rl[x,1],lyrs = as.numeric(rl[x,2]))
      outrast<-aggregate(r1,fact = f1,fun = aggfun,zm1=donut)
      nm<-c(paste("f-",unique(round(donut,2)),sep = ""),
            paste("fp-",unique(round(donut,2)),sep = ""),
            c("mean","Q025","med","Q95","sd","skew","kurt"))
      names(outrast)<-nm
      rn<-strsplit(r,"/")[[1]];rn<-rn[length(rn)]
      writeRaster(outrast,filename = paste(ap,x,rn,sep = "/"),datatype = "INT4S",overwrite = overwrite)
      rm(list =c("outrast","r1"))
      gc()
      return(x)
    })
    
  }

  ol
}

## GetMetrics ------------
#' Calculates textural summaries of input images at an aggregated scale.
#'
#' @description This function is useful for imputation mapping from models built with the yaImpute package.
#'
#' @export
#'
#' @param rasterfile character vector with tif image names. If full path not included, GetMetrics will assume that these files are in 0_raw folder.
#' @param outpath1 file path pointer where R has already stored intermediate (full resolution) indices 
#' @param outpath2 file path pointer where R will store aggregated image tiles (reduced resolution), texture summaries across all bands and indices.
#' @param outres resolution of output texture summaries
#' @param ncpus specify number of cpus to use for parallel processing. 
#' @param zradii numeric vector specifying frequency-length bins for summarizing fft spectrum. (maximum value of this should be larger than diagonal of output pixel size)
#' @seealso \code{\link{MakeDonut}}
#' @return multiband raster
#'
GetMetrics<-function(rasterfile,outpath1 = "1_intermediate", outpath2 = "2_aggregated", outres = 30,ncpus = 3,zradii = c(3,6,10,20,40,55)){
  require(terra)
  if(!grepl("/",rasterfile)){rasterfile<-paste("0_raw/",rasterfile,sep = "")}
  fn<-strsplit(rasterfile,split = "/")[[1]];fn<-fn[length(fn)]
  intfile<-paste(outpath1,fn,sep = "/")
  outfile<-paste(outpath2,paste("afft_",fn,sep = ""),sep = "/")
  if(file.exists(outfile)){
    outrast<-rast(outfile)
    res1<-res(rast(rasterfile))[1]
    fact1<-outres/res1
    donut<-round(c(MakeDonut(zradii,res1,fact1)),2)
  }else{
    cat("  \n  Elapsed Time 1:",system.time({
      cat("\n",rasterfile)
      file.remove(paste(outpath1,"/temp.tif",sep = ""))
      file.copy(rasterfile,paste(outpath1,"/temp.tif",sep = ""))
      r1<-rast(paste(outpath1,"/temp.tif",sep = ""));names(r1)<-c("r","g","b","n")
      cat("  making NDVI")
        num<-(r1$n - r1$r)
        denom<-(r1$n + r1$r)
        ndvi<-num/denom
      cat("  making NDGR")
        num<-(r1$g - r1$r)
        denom<-(r1$g + r1$r)
        ndgr<-num/denom
      cat("  making NDNG")
        num<-(r1$n - r1$g)
        denom<-(r1$n + r1$g)
        ndng<-num/denom
      
      cat("  stretching Indices")
        nd<-c(ndvi,ndgr,ndng)
        nd<-stretch(nd,smin = -1, smax = 1)
        names(nd)<-c("ndvi","ndgr","ndng")
        
      cat("  making brightness")
       br<-sum(r1)#app(r1,sum)
       br<-stretch(br, smin = 0, smax = 1020)
      r1<-c(r1,nd,br)  
      names(r1)<-c("r","g","b","n","ndvi","ndgr","ndng","bri")
      r1<-terra::subset(r1,c("r","g","n","ndvi","ndgr","ndng","bri"))#,"b" removed
      
      r1<-writeRaster(r1,filename = intfile,overwrite = T)
      r1
    })[[3]]/60)
    fact1<-outres/res(r1)[1]
    donut<-round(c(MakeDonut(zradii,res1 = res(r1)[1],fact1)),2)
    cat("  \n  Elapsed Time 2:",system.time({
      names1<-paste("f",xfun::numbers_to_words(ceiling(unique(donut))),letters[1:length(unique(donut))],sep = "_")
      names2<-c("mean","Q05","Q10","med","Q90","Q95","sdX10","skewX10","kurtX100")
      bandnames<-c("r","g","n","ndvi","ndgr","ndng","bri")
      basenames<- c(names1,names2)
      nms<-do.call(c,lapply(bandnames,function(x){paste(x,basenames,sep = "_")}))
      
      outrast<-aggregate(r1,fact = fact1,cores = ncpus,zm1 = donut,
                         fun = function(x,zm1 = donut){
                           if(!any(is.na(x))){
                             x<-matrix(x,nrow = ceiling(sqrt(length(x))), byrow = T)
                             ff1<-c(abs(gsignal::fftshift(fft(x - mean(x)),MARGIN =c(1,2))^2))
                             y<-tapply(ff1,zm1,sum)/sum(ff1)
                             y2<-c(mean(x),
                                   quantile(x,c(.05,.1,.5,.9,.95)),
                                   sd(x)*10,
                                   (moments::skewness(c(x))+10)*10,
                                   log((moments::kurtosis(c(x))))*100)
                             outvec<-round(c(y*1000,y2*100),0)
                           }else{
                             outvec<-rep(NA,(length(unique(zm1))+9))
                           }
                           #print(length(outvec))
                           return(outvec)
                         })
      
    })[[3]]/60,"\n")
    names(outrast)<-nms
    writeRaster(outrast,filename = outfile,datatype = "INT2U")
    beepr::beep(10)
  }
  names(outrast)<-nms
  file.remove(paste(outpath1,"/temp.tif",sep = ""))
  return(outrast)
}  
## MakeDonut ------------
#' Used within GetMetrics.
#'
#' @description Creates a matrix whose outer size matches the dimension of the image to be summarized (in this case, an airphoto image covering exactly one landsat pixel)
#'
#' @export
#'
#' @param zr Radii for nested circles
#' @param res1 output resolution (full size indicated by the output matrix)
#' @param fact1 resolution of imagery to be summarized
#' @seealso \code{\link{GetMetrics}}
#' @return matrix with integers that is used for extracting zonal summaries of fft spectrum. Donut values indicate scale of variation in meters.
#'
MakeDonut<-function(zr,res1,fact1,return.rast = T){
  d<-do.call(c,lapply(zr,function(r,res2 = res1[1],f2= fact1[1]){
    size <- f2*res2[1]
    e = c(-size/2,size/2,-size/2,size/2)
    pt1<-terra::vect(rbind(c(0,0)))
    ya<-terra::rast(nrows=size/res2, ncols=size/res2, nlyrs=1,  extent = terra::ext(e), 
                    resolution = res2)
    y<-terra::rasterize(terra::buffer(pt1, r),ya)
    y<-y * (1/r) * size/2
  }))
  d<-terra::app(d,function(x){y<-max(x,na.rm = T)})
  if(return.rast)return(d)
  else{return(d[])}
}

## MergeAggregatedAirphotos ------------
#' Used after GetMetrics.
#'
#' @description merges tiles in tiflist, one tif per band in the tiles.
#'
#' @export
#'
#' @param tiflist list of files to merge
#' @param outpath2 file path to folder where output of this function should be stored (usually same directory as tiflist).
#' @return list of images. Mostly called to generate images that are easy to read in and work with later.
#'
MergeAggregatedAirphotos<-function(tiflist,outpath2 = outpath2){
  ## Consider allowing terra to use more RAM?  How to do this?
  ## Handling of names in GetMetrics may require revisions in this function.
  mergefun<-function(x,tfn,tmplist = tmp,outpath.cur = outpath2){
    tmpfile<-paste(outpath.cur,tfn,x,".tif",sep = "")
    if(file.exists(tmpfile)){
      return(tmpfile)
    }else{
      y<-vrt(tmp[[x]],paste(outpath.cur,"/tmp",sep = ""),return_filename = T,overwrite = T)
      y<-rast(y)
      y<-writeRaster(y,tmpfile)
      return(sources(y))
    }
    
  }
  tmp<-reshape2::colsplit(names(tiflist),"_",c(letters[1:8]))
  tmp<-split(tiflist,tmp$c)
  nm<-names(tmp);names(nm)<-nm
  cat("merge1\n")
  m1<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m1_"))
  
  tmp<-split(m1,substr(names(m1),1,nchar(names(m1))-1))
  nm<-names(tmp);names(nm)<-nm
  cat("merge2\n")
  m2<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m2_"))
  
  tmp<-split(m2,substr(names(m2),1,nchar(names(m2))-1))
  nm<-names(tmp);names(nm)<-nm
  cat("merge3\n")
  m3<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m3_"))
  
  tmp<-split(m3,substr(names(m3),1,nchar(names(m3))-2))
  nm<-names(tmp);names(nm)<-nm
  cat("merge4\n")
  m4<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m4_"))
  
  ### at a certain stage, consider: merging individual layers may well be a more sensible 
  ## operation than merging the whole dang thing. 
  
  
  # y<-vrt(m4,paste(outpath2,"/tmp",sep = ""),return_filename = T,overwrite = T)
  # y<-rast(y)
  bandnames<-c("r", "g","n","ndgr","ndng","ndvi","bri" )
  statnames<-c("f-one-a","f-one-b","f-two-c","f-three-d","f-five-e","mean","min","Q10","med","Q90","max", "sdx10","skewx10","kurtx100")
  nms4<-do.call(c,lapply(statnames,function(x){paste(x,bandnames,sep = "_")}))
  names(y)<-nms4
  names(nms4)<-nms4
  if(length(nms4)!=dim(rast(m4[1]))[3]){print("STOP - names and number of bands do not match!");browser()}
  q<-pbapply::pblapply(nms4,function(x){
    cl<-sprc(lapply(m4,function(z){subset(rast(z),which(nms4==x))}))
    outlayer<-merge(cl,filename = paste(outpath2,"/AFFT_Initial_",x,".tif",sep = ""),overwrite = T)
  })
  
  
  
  ## cleanup interrim files
  interrim.files<-c(m1,m2,m3,m4)
  rm1<-sapply(interrim.files,file.remove)
  
  return(y)  
}


## MergeAggregatedAirphotos2 ------------
#' Used after GetMetrics.
#'
#' @description merges tiles in tiflist, one tif per band in the tiles.
#'
#' @export
#'
#' @param tiflist list of files to merge
#' @param filename name of output file
#' @param outpath2 file path to folder where output of this function should be stored (usually same directory as tiflist).
#' @return output image path. Mostly called to generate images that are easy to read in and work with later.
#'
MergeAggregatedAirphotos2<-function(tiflist,filename, outpath = "2_aggregated"){
  ## Consider allowing terra to use more RAM?  How to do this?
  ## Handling of names in GetMetrics may require revisions in this function.
  mergefun<-function(x,tfn,tp = tmp,op = outpath){
    tmpfile<-paste(op,tfn,x,".tif",sep = "")
    if(file.exists(tmpfile)){
      return(tmpfile)
    }else{
      browser()
      vrtfile<-paste(op,"/tmp",sep = "")
      sprc1<-sprc(tp[[x]])
      y<-terra::merge(sprc1,filename = tmpfile)
      y<-rast(tp[[x]][1])
      for(i in 2:length(tp)){y<-merge(y,rast(tp[[x]]))}
      #y<-terra::vrt(tp[[x]],vrtfile,overwrite = T)
      y<-terra::rast(vrtfile)
      y<-terra::writeRaster(y,tmpfile)
      return(terra::sources(y))
    }
  }
  tmp<-reshape2::colsplit(names(tiflist),"_",c(letters[1:8]))
  tmp<-split(tiflist,tmp$b)
  nm<-names(tmp);names(nm)<-nm
  cat("merge1\n")
  m1<-do.call(c,pbapply::pblapply(nm,mergefun,"/m1_"))
  
  tmp<-split(m1,substr(names(m1),1,nchar(names(m1))-1))
  nm<-names(tmp);names(nm)<-nm
  cat("merge2\n")
  m2<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m2_"))
  
  tmp<-split(m2,substr(names(m2),1,nchar(names(m2))-1))
  nm<-names(tmp);names(nm)<-nm
  cat("merge3\n")
  m3<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m3_"))
  
  tmp<-split(m3,substr(names(m3),1,nchar(names(m3))-2))
  nm<-names(tmp);names(nm)<-nm
  cat("merge4\n")
  m4<-do.call(c,pbapply::pblapply(nm,mergefun,tfn = "/m4_"))
  
  ### at a certain stage, consider: merging individual layers may well be a more sensible 
  ## operation than merging the whole dang thing. 
  outfile<-paste(outpath,"/AFFT_Initial_",filename,".tif",sep = "")
  
  
  vrtfile<-paste(outpath,"/tmp",sep = "")
  y<-terra::vrt(m4,vrtfile,overwrite = T)
  y<-terra::rast(vrtfile)
  y<-terra::writeRaster(y,outfile,overwrite = T)
  
  
  
  ## cleanup interrim files
  interrim.files<-c(m1,m2,m3,m4)
  rm1<-sapply(interrim.files,file.remove)
  
  return(y)  
}

## MakePCA ------------
#' Perform principal components analysis to help reduce dimensionality of a multiband image.
#'
#' @description Creates a matrix whose outer size matches the dimension of the image to be summarized (in this case, an airphoto image covering exactly one landsat pixel)
#'
#' @export
#'
#' @param r multiband rast image.
#' @param sampdens sampling density for building PCA
#' @param maxsample maximum number of samples to draw from image.
#' @return multiband raster containing pca axis scores.
#'
MakePCA<-function(r,sampdens = 50,maxsample = 100000){
  samp<-spatSample(r,size = min(maxsample,ncell(r)/sampdens))
  samp<-samp[apply(samp,1,function(x){!any(is.na(x))}),]
  cat(paste("sampsize:",nrow(samp)))
  pca1<-princomp(samp)
  pcamap<-predict(r,pca1)
  pcamap
}

## CorrectLayers ------------
#' Adjusts imagery to reduce artifacts associated with flightlines, and phenology inconsistencies.
#'
#' @description Corrects image artifacts based on seamlines, image dates in airphoto metadata polygons.  Alternate sources of imagery corrections could be used.
#'
#' @export
#'
#' @param y raster object created from seamlines, containing two layers, one describing distance from seamline, the other indicating image date.
#' @param x raster to be corrected
#' @param maxsample maximum pixels to be sampled for corrective model
#' @param sampdens sampling density to be sampled for corrective model.
#' @return multiband raster containing layers adjusted according to information associated with image acquisition date, and also distance to seamline (note: flightline distance would be better, but not currently available)
#'
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
## ScoreNoise ------------
#' Extracts fourier summary of image, returns estimate of how much image variability is at a pixel-to-pixel scale (a.k.a. noise)
#'
#' @description Extracts fourier summary of image, returns estimate of how much image variability is at a pixel-to-pixel scale (a.k.a. noise)
#'
#' @export
#'
#' @param y rast object
#' @param f1 resolution to extract sub-images for fourier analysis.  256 currently used for efficiency
#' @param f2 secondary resolution for summarizing fourier spectrum. Set to 64 for computational efficiency.
#' @param q if set to T, R will print out images as it scores noise for the large image tiles.
#' @param scorefun specifies how variability in noise scores for tiled fft summaries will be stored
#' @param return.rast if set to T, function returns the image of noise scores with a large pixel size (aggregate of original image to to f1)
#' @return statistics describing the proportion of image variation at very small scales.
#'
ScoreNoise<-function(y,f1 = 256,f2 = 64, q = T,
                     scorefun = function(x, q = c(.05,.1,.5,.9,.95)){quantile(x,q,na.rm = T)},
                     return.rast = F){
  ## power of two block size may run faster.
  ## odd-sized images are 'neat'?
  ## Note: edge effects -- impacts on high frequency detections (function wraps image). Consider adjusting later?
  ## Note: sharp boundaries in original image will contribute to high frequencies detected.
  ## Note: feeding this function to a cluster within a wrapper, with one image per core is faster than using a multicore call to aggregate function.#, cores = 1 

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
## ScoreArtifacts ------------
#' Extracts fourier summary of image, returns estimate of how much image variability is at a pixel-to-pixel scale (a.k.a. noise)
#'
#' @description Extracts fourier summary of image, returns estimate of how much image variability is at a pixel-to-pixel scale (a.k.a. noise)
#'
#' @export
#'
#' @param x rast object
#' @param y rast object containing layers that may indicate artifacts in teh imagery (e.g., flightlines, or phenology blocks)
#' @param maxsample maximum pixels to sample
#' @param sampdens density to sample pixels (up to maxsample)
#' @return list of linear models describing correlations between x and y.
#'
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
## ViewStackRGB ------------
#' Useful for viewing multi-band images in rgb space.  Groups bands 1:3, 4:6, ... through the number of layers in the input raster.
#'
#' @description Extracts fourier summary of image, returns estimate of how much image variability is at a pixel-to-pixel scale (a.k.a. noise)
#'
#' @export
#'
#' @param x rast object
#' @return returns NULL.  Used for plotting
#'
ViewStackRGB<-function(x){
  vec<-sort(rep(1:(ceiling(dim(x)[3]/3)),length.out = dim(x)[3]))
  #vec<-split(names(x),vec);names(vec)<-sapply(vec,function(z){paste(z,collapse = "_")})
  vec<-split(1:length(names(x)),vec)
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
  return(NULL)
}

## pause -----------------
#'
#' function for introducing a brief pause into a script
#' @description pauses R's operation 
#' 
#' @export
#' 
#' @param secs seconds for pause
#' @return NULL
#'

pause<-function(secs){
  st<-Sys.time()
  while(abs(difftime(st,Sys.time(),units = "secs"))<secs){}
  return(NULL)
}

## indexFuns ------------
#' functions for tallying indices on raw airphotos. pass this object (or revised indices) to cluster with sfExport.
#'
#' @description List of functions for tallying indices on raw airphotos. pass this object (or revised indices) to cluster with sfExport.
#'


indexFuns<-list(
  ndvi = function(r,outpath_i = getwd()){
    fn<-strsplit(r,"/")[[1]];fn<-fn[length(fn)]
    fn2<-paste(outpath_i,"/1_intermediate/ndvi/",fn,sep = "")
    r<-terra::rast(r,lyrs =c(1,4))
    names(r)<-c("one","two")
    y<-as.int(  ((  (r$one-r$two)/(r$one+r$two)  )+1) *126  )
    writeRaster(y,filename = fn2, datatype = "INT1U",overwrite = T)
    rm(list = c("r","y"))
    
    gc()
    fn2},
  ndgr = function(r,outpath_i = getwd()){
    fn<-strsplit(r,"/")[[1]];fn<-fn[length(fn)]
    fn2<-paste(outpath_i,"/1_intermediate/ndgr/",fn,sep = "")
    r<-terra::rast(r,lyrs =c(2,1))
    names(r)<-c("one","two")
    y<-as.int(  ((  (r$one-r$two)/(r$one+r$two)  )+1) *126  )
    writeRaster(y,filename = fn2, datatype = "INT1U",overwrite = T)
    rm(list = c("r","y"))
    
    gc()
    fn2
    
  },
  ndng = function(r,outpath_i= getwd()){
    fn<-strsplit(r,"/")[[1]];fn<-fn[length(fn)]
    fn2<-paste(outpath_i,"/1_intermediate/ndng/",fn,sep = "")
    r<-terra::rast(r,lyrs =c(4,2))
    names(r)<-c("one","two")
    y<-as.int(  ((  (r$one-r$two)/(r$one+r$two)  )+1) *126  )
    writeRaster(y,filename = fn2, datatype = "INT1U",overwrite = T)
    rm(list = c("r","y"))
    gc()
    fn2
  },
  bri = function(r,outpath_i= getwd()){
    fn<-strsplit(r,"/")[[1]];fn<-fn[length(fn)]
    fn2<-paste(outpath_i,"/1_intermediate/bri/",fn,sep = "")
    r<-terra::rast(r)
    y<-as.int(sum(r)/4)
    writeRaster(y,filename = fn2, datatype = "INT1U",overwrite = T)
    rm(list = c("r","y"))
    gc()
    
    fn2}
)
