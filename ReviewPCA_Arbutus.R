library(terra)
library(AFFT)
setwd("C:/TEMP/RockyMountains")

bnd<-vect("Y:/MPSG_VegMapping/Data/Spatial/Mapping boundaries final/MPSG_mapping_boundaries_final.shp")
bnd0<- bnd<-project(bnd,rast(list.files("Y:/MPSG_VegMapping/Data/Raster/Predictors/elevation",full.names = T,pattern = "tif")[1]))
bnd<-bnd[bnd$FORESTORGC %in% c("0118","0207"),]#
bnd<-buffer(bnd,100000)
bnd<-aggregate(bnd)

fl<-list.files("PCA",full.names = T,pattern = ".tif")
fl<-fl[!grepl("_",fl)]
v1<-vrt(fl)
print(dim(v1))
names(v1)<-names(rast(fl[1]))
windows(10,10);plot(v1,1);e1<-draw();e3<-draw();plot(v1,1,ext = e1);e2<-draw();plot(v1,1,ext = e3);e4<-draw();dev.off()

par(mfrow =c(2,3));sn<-1;ce<-e2
for(i in c(-1+1:floor(dim(v1)[3]/3))+(sn-1)/3){
  j<-i*3;plotRGB(subset(v1,1:3 + j),stretch = "lin",
          main = paste("Components",j+1,"-",j+3)
          ,ext = ce)
  plot(bnd0,add = T,border = "black")
}

  tmp<-cbind(pca1$loadings)

layout.show(lo<-layout(matrix(c(1,2,3,1,4,5),byrow = T,nrow = 2),widths =c(1.25,1,1)))

par(mar =c(0,0,0,0))
mc<-100000#1:dim(v1)[3]
badcomps<-sapply(1:dim(v1)[3],function(i){
  cp<-hcl.colors(255,palette = "viridis")
  v0<-subset(v1,i)
  plot(v0,main = paste("Comp.",i,sep = ""),legend = F,
       mar =c(0,0,0,0),axes = F,maxcell = mc*3,col = cp)
  plot(vect(e1),border = "red",add = T,lwd = 3)
  plot(vect(e3),border = "orange",add = T,lwd = 3)
  plot(v0,ext = e1,legend = F,maxcell = mc*2,axes = F,
       mar =c(0,0,0,0),col = cp)
  plot(vect(e1),border = "red",add = T,lwd = 3)
  plot(vect(e2),border = "yellow",add = T,lwd = 3)
  plot(v0,ext = e2,legend = F,maxcell = mc*2,axes = F,
       mar =c(0,0,0,0),col = cp)
  plot(vect(e2),border = "yellow",add = T,lwd = 3)
  
  plot(v0,ext = e3,legend = F,maxcell = mc*2,axes = F,
       mar =c(0,0,0,0),col = cp)
  plot(vect(e3),border = "orange",add = T,lwd = 3)
  plot(vect(e4),border = "magenta",add = T,lwd = 3)
  plot(v0,ext = e4,legend = F,maxcell = mc*2,axes = F,
       mar =c(0,0,0,0),col = cp)
  plot(vect(e4),border = "magenta",add = T,lwd = 3)
  readline(paste("Comp.",i,": ",sep = ""))
})

  

v2<-rast(list.files("C:/TEMP/RockyMountains/Corrected",pattern = "vrt",full.names = T))
{
topaxis<-list()
for(comp in c(1,3,4,6,10)){#which(badcomps == "x")
  tmp<-tmp[order(abs(tmp[,comp]),decreasing = T),]
  nm<-rownames(tmp);names(nm)<-nm
  topaxis[[as.character(comp)]]<-lapply(nm[1:10],function(x){
    lo<-layout(matrix(c(5,5,1,1,1,5,5,2,3,4),byrow = T,nrow = 2),heights =c(1.5,4), widths =c(.3,1,1,1,1))
    par(mar =c(3,2,3,0))
    barplot(tmp[x,1:dim(v1)[3],drop = F], las = 2, cex.names = .5,main = x)
    r1<-subset(v2,x)
    par(mar =c(0,0,0,0))
    ea<-e3
    eb<-e4
    plot(r1, legend = F,mar =c(0,0,0,0),axes = F);plot(vect(ea),border = "orange",add = T);plot(vect(eb),add = T, border = 'red')
    plot(r1,ext = ea, legend = F,mar =c(0,0,0,0),axes = F);plot(vect(ea),border = "orange",add = T);plot(vect(eb),add = T, border = "red")
    plot(r1,ext = eb, legend = F,mar =c(0,0,0,0),axes = F);plot(vect(eb),add = T, border = "red")
    plot(v1,comp,main = paste("Comp.",comp,sep = ""),legend = F,maxcell = 300000,mar =c(0,0,0,0),axes = F);plot(vect(ea),border = "orange",add = T);plot(vect(eb),add = T, border = 'red')
    y<-readline(x)
    return(y)})
}

cat(paste("'",paste(unique(do.call(c,lapply(topaxis,function(x){y<-do.call(c,x)=="x";names(y)[y]}))),collapse = "','"),"'",sep = ""))
}  


rs<-data.frame(varname = paste("AFFT_PCA_RM_",1:dim(v1)[3],sep = ""))
rs$filepath<-paste("Y:/MPSG_VegMapping/Data/Raster/Predictors",rs$varname,sep = "/")
rs$screening<-1
rs$screening[badcomps %in% "x"]<-0
rs$Dakotas<-0
rs$RockyMountains<-1
rs$ComancheCimarrone<-0
rs$datatype <- "AFFT_PCA"

write.csv(rs,"PCA/Rastersource.csv")
  setwd("PCA")
  fl<-list.files(,pattern = ".tif")
for(i in 1:dim(v1)[3]){
  cat("\n####",i,"####\n")
  d1<-SDMap::MakeDirIfNeeded(rs$varname[i],goto = F)
  
  pbapply::pblapply(fl,function(x){y<-subset(rast(x),i)
    fn1<-paste(d1,rs$varname[i],sep = "/")
    fn2<-gsub("pca1",fn1,x)
    writeRaster(y,filename = fn2,overwrite = T)
  })
}
  