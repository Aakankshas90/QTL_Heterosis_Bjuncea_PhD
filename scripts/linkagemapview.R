library(LinkageMapView)
library(qtl) 
library(ASMap)
arraylmv<-read.csv(file = "arraymarkers_3_lmv.csv", header = TRUE, sep = ",")
data(arraymarkers_3_lmv.csv)  
class(arraylmv)
summary(arraylmv)
lmv.linkage.plot(mapthis = arraylmv, outfile = arraylmv, denmap = TRUE, lgw = 0.5, pdf.width = 20, pdf.height = 15)


                

------------uniquemarkers------------
  setwd("D:/R/")
library(LinkageMapView)
library(qtl) 
library(ASMap)
uniquemarkerslmv<-read.csv(file = "uniquemarkers.csv", header = TRUE, sep = ",")
class(uniquemarkerslmv)
summary(uniquemarkerslmv)
lmv.linkage.plot(mapthis = uniquemarkerslmv, outfile = uniquemarkerslmv, denmap = FALSE)


#-----------------ej8--------------------------------------------
setwd("D:/lmv/")
library(LinkageMapView)
library(qtl) 
library(ASMap)
EJZ_lmv<-read.csv(file = "EJZlmv.csv", header = TRUE, sep = ",")
class(EJZ_lmv)
summary(EJZ_lmv)
lmv.linkage.plot(mapthis = EJZ_lmv, outfile = EJZ_lmv, denmap = TRUE, lgw = .35, pdf.width = 15, pdf.height = 10)



#--------------------------------------map7609------------------------------------

setwd("D:/lmv/")
library(LinkageMapView)
library(qtl) 
library(ASMap)
map7609_lmv<-read.csv(file = "map7609.csv", header = TRUE, sep = ",")
class(map7609_lmv)
summary(map7609_lmv)
lmv.linkage.plot(mapthis = map7609_lmv, outfile = map7609_lmv, denmap = TRUE,  pdf.width = 20, pdf.height = 15)


#--------------------------------------MAP4463------------------------------------

setwd("D:/lmv/")
library(LinkageMapView)
library(qtl) 
library(ASMap)
MAP4463_lmv<-read.csv(file = "MAP4463.csv", header = TRUE, sep = ",")
class(MAP4463_lmv)
summary(MAP4463_lmv)
lmv.linkage.plot(mapthis = MAP4463_lmv, outfile = MAP4463_lmv, denmap = TRUE,  pdf.width = 20, pdf.height = 15)
