
library(terra);library(AFFT)

rawfiles<-list.files("0_raw",full.names = F,pattern = ".tif")
indexfiles<-list.files("1_intermediate/ndvi")#,full.names = F))
rawfiles<-rawfiles[!rawfiles %in% indexfiles]
rawfiles<-paste("0_raw",rawfiles,sep = "/")
terraOptions(memfrac = .9,datatype = "INT4U")




MAA<-function(tiflist,filename, outpath = "2_aggregated"){
  ## Consider allowing terra to use more RAM?  How to do this?
  ## Handling of names in GetMetrics may require revisions in this function.
  mergefun<-function(x,tfn,tp = tmp,op = outpath){
    tmpfile<-paste(op,tfn,x,".tif",sep = "")
    if(file.exists(tmpfile)){return(tmpfile)
    }else{
       if(length(tp[[x]])>1){
        v1<-vrt(tp[[x]])
        names(v1)<-names(rast(tp[[x]][1]))
        y<-terra::writeRaster(v1,tmpfile,overwrite = T)
        return(terra::sources(y))
      }else{return(tp[[x]])}
    }
   
  }
  tmp<-reshape2::colsplit(names(tiflist),"_",c(letters[1:8]))
  tmp<-split(tiflist,tmp$b)
  nm<-names(tmp);names(nm)<-nm
  mergeind<-0
  nfiles<-length(tiflist)
  while(nfiles>4){
    cat("merge",mergeind<-mergeind+1,"\n\n")
    m1<-do.call(c,pbapply::pblapply(nm,mergefun,paste("/m",mergeind,"_",sep = "")))
    nfiles<-length(m1)
    tmp<-split(m1,substr(nm,1,nchar(nm)-1))
    nm<-names(tmp);names(nm)<-nm
    cat(nm)
    
  }
  
    v1<-vrt(m1)
    names(v1)<-names(rast(m1[1]))
  
    out<-writeRaster(v1,filename = paste(outpath,"/AFFT_Initial_",filename,".tif",sep = ""),datatype = "FLT4S",overwrite = T)
  m1<-"X"
  v1<-"X"
  ## cleanup interrim files
    sapply(paste("m",1:mergeind,sep = ""),function(x){sapply(list.files(outpath,x,,T),file.remove)})

  return(sources(out))  
}


fixproj<-function(tiflist){
  t1<-sort(table(pl<-sapply(tiflist,function(x){crs(rast(x),proj = T)})),decreasing = T)
  if(length(t1)>1){
    toproj<-names(t1)[1]
    fromgroups<-names(t1)[2:length(t1)]
    t2<-split(tiflist,pl)
    for(p in fromgroups){
      z<-pbapply::pblapply(t2[[p]],function(d){
        sf<-d
        d<-rast(d)
        q1<-project(d,rast(t2[[toproj]][1]),align = T)
        names(q1)<-names(d)
        writeRaster(q1,sf,overwrite = T,datatype = "INT4S")
     })
    }
  }
  return(NULL)
}
fixborders<-function(tiflist){
  mskfun<-function(x){
    y<-rast(x)
    msk<-as.numeric(subset(y,1) ==0)
    y<-mask(y,msk,maskvalue = 1, updatevalue = NA)
    writeRaster(y,x,overwrite = T,datatype = "INT4S")
  }
  tmp<-pbapply::pbsapply(tiflist,mskfun)
  return(NULL)
}
lapply(c("bri","g","n","ndgr","ndng","r"),function(x){
  fl<-list.files(paste("2_aggregated/",x,sep = ""),".tif",,T)
  nm<-list.files(paste("2_aggregated/",x,sep = ""),".tif",)
  names(fl)<-nm
  length(fl)
  fixproj(fl)
  MAA(fl,x)
})
