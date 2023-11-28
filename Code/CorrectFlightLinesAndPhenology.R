sl<-vect("Spatial/Imagery/NAIP/co19/shp/Seamlines_co08_2019.shp")## -- not current path.
sl<-project(sl,CO)
sl$JDate<-julian(as.Date(sl$IDATE),origin = as.Date("2019-01-01"))
sl2<-crop(sl,buffer(vect(ext(r1)),20))
sl2<-as.lines(sl2)

jdate<-rasterize(sl,r1,"JDate")
jdate<-crop(jdate,r1)
jdate<-as.factor(jdate)
sl2<-distance(r1,sl2,rasterize = T)


CorrectLayers<-function(y,x,maxsample= 75000 ,sampdens = 50){
  mat1<-spatSample(c(x,y),size = min(maxsample,ncell(y)/sampdens))
  mat1<-mat1[apply(mat1,1,function(x){!any(is.na(x))}),]
  nm1<-names(y);names(nm1)<-nm1
  outrastl<-lapply(nm1,function(z){
    mat2<-mat1[,c(z,names(x))]
    colnames(mat2)[1]<-"yvar"
    lm1<-lm(yvar~.,data = mat2)
    outr1<-predict(x,lm1)
    outr2<-subset(y,z)-outr1
    outr2
  })
  rast(outrastl)
}

cl<-CorrectLayers(r1,c(jdate,sl2))

