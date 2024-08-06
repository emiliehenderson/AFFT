library(terra)
library(AFFT)
setwd("C:/TEMP/ComancheCimarrone")


fl<-list.files("C:/TEMP/ComancheCimarrone/PCA",full.names = T,pattern = "f.tif");v1<-vrt(fl)
print(dim(v1))
# for(i in 0:(dim(v1)[3]/3 - 1)){par(mfrow =c(1,3))
#   v2<-subset(v1,1:3+i*3)
#   plotRGB(v2,stretch = "lin")
#   plot(vect(e1),add = T)
#   plot(vect(e2),add = T)
#   plotRGB(v2,stretch = "lin",ext = e1)
#   plot(vect(e2),add = T)
#   plotRGB(v2,stretch = "lin",ext = e2)
#           readline(i*3+1)}
# 



par(mfrow =c(1,2))
sn<-1
while(sn < dim(v1)[3] - 3){
  for(i in c(0:1)+(sn-1)/3){
    j<-i*3
    plotRGB(subset(v1,1:3 + j),stretch = "lin",main = paste("Components",j+1,"-",j+3))
    plot(bnd0,add = T,border = "red")
  };sn<-sn+6
}




tmp<-cbind(pca1$loadings)

comp<-2
tmp<-tmp[order(abs(tmp[,comp]),decreasing = T),]
nm<-rownames(tmp);names(nm)<-nm
topaxis<-lapply(nm[1:10],function(x){
  lo<-layout(matrix(c(5,2,1,3,1,4),byrow = F,nrow = 2),heights =c(1,1), widths =c(1.5,1,1))
  par(mar =c(3,2,3,0))
  barplot(tmp[x,1:dim(v1)[3],drop = F], las = 2, cex.names = .5,main = x)
  r1<-rast(paste("C:/TEMP/ComancheCimarrone/Corrected/",x,".vrt",sep = ""))
  par(mar =c(0,0,0,0))
  plot(r1, legend = F);plot(vect(e1),border = "orange",add = T);plot(vect(e3),add = T, border = 'red')
  plot(r1,ext = e1, legend = F);plot(vect(e1),border = "orange",add = T)
  plot(r1,ext = e3, legend = F);plot(vect(e3),add = T, border = "red")
  plot(v1,comp,main = paste("Comp.",comp,sep = ""),legend = F);plot(vect(e1),border = "orange",add = T);plot(vect(e3),add = T, border = 'red')
  y<-readline(x)
  return(y)})

  #return(colnames(tmp[which.max(tmp[x,])]))})


dm<-c()
for(i in 1:15){ ## plots band, top 16 band components
  par(mfrow =c(1,2))
  cv<-subset(v1,i)
  plotRGB(c(cv,cv,cv),stretch = "hist")#
  plotRGB(c(cv,cv,cv),stretch = "hist",ext = e1)#
  mtext(paste("Component:",i),col = "red",font = 2, cex = 1.2)
  readline(i)
  fn<-rownames(tmp)[order(abs(tmp[,i]),decreasing = T)[1:5]]
  for(j in fn){
    r1<-rast(paste("C:/TEMP/ComancheCimarrone/Corrected/",j,".vrt",sep = ""))
    plotRGB(c(r1,r1,r1),stretch = "hist")
    plotRGB(c(r1,r1,r1),ext = e1,stretch = "hist")
    mtext(j,line = 0, font = 2)
    if(readline(j) %in% "x"){dm<-c(dm,j)}
  }
}


# 
# Correction options:
#   1) Add suppl to pca -- doesn't work well. Handles phenology poorly.
#   2) No pre-pca corrections -- artifacts grow dominant somewhere round about band 15-18.
#   3) Corrections pre-pca
#   4) Weed out input corrected images based on artifacts, and noise.  
#     a) Noise is consistently high for fp-20 
#     b) Filter out flightline, or phenology scores = 5. (double-check)
#     c) filter out average score >=4? - doublecheck.
#   5) dropped flightline scores >= 4, phenology scores = 5, NoiseScores 10% >=600.  Then, filtered additional variables by hand.
      # 72 Variables used: 
#     [1] "bri_f-12"    "bri_f-20"    "bri_fp-1.5"  "bri_fp-3"    "bri_fp-6"    "bri_kurt"    "bri_skew"   
#[8] "g_f-12"      "g_f-20"      "g_f-6"       "g_fp-12"     "g_fp-3"      "g_fp-6"      "g_kurt"     
#[15] "g_sd"        "g_skew"      "n_fp-0.25"   "n_fp-1.5"    "n_fp-12"     "n_fp-3"      "n_fp-6"     
#[22] "n_kurt"      "n_mean"      "n_med"       "n_Q95"       "n_skew"      "ndgr_f-12"   "ndgr_f-20"  
#[29] "ndgr_f-6"    "ndgr_fp-1.5" "ndgr_fp-12"  "ndgr_fp-3"   "ndgr_fp-6"   "ndgr_kurt"   "ndgr_mean"  
#[36] "ndgr_med"    "ndgr_Q025"   "ndgr_sd"     "ndgr_skew"   "ndng_f-12"   "ndng_f-20"   "ndng_f-3"   
#[43] "ndng_f-6"    "ndng_fp-3"   "ndng_fp-6"   "ndng_kurt"   "ndng_Q025"   "ndng_sd"     "ndng_skew"  
#[50] "ndvi_f-1.5"  "ndvi_f-12"   "ndvi_f-20"   "ndvi_f-3"    "ndvi_f-6"    "ndvi_fp-3"   "ndvi_fp-6"  
#[57] "ndvi_kurt"   "ndvi_mean"   "ndvi_med"    "ndvi_Q025"   "ndvi_Q95"    "ndvi_sd"     "ndvi_skew"  
#[64] "r_f-12"      "r_f-20"      "r_f-6"       "r_fp-12"     "r_fp-3"      "r_fp-6"      "r_kurt"     
#[71] "r_sd"        "r_skew"     
  