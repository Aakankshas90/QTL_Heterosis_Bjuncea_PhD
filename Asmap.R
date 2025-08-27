describe(pheno)
library(psych)
setwd("D:/R/")
getwd()
data15<-read.csv(file = "Data14-15.csv", header = TRUE, sep = ",")
data15
View(data15)
options(max.print = 999999)
data15
library(psych)
corPlot(data15)
describe(data15)
corr.test(data15)
print.psych(x, digits = 2, all = FALSE, cut = NULL, sort = FALSE, short = TRUE, signif = NULL)


install.packages("devtools")
library(devtools)
install_github("mckaylab/TSPmap")
install.packages("qtl2", repos="http://rqtl.org/qtl2cran")
mapdata<-read.csv(file = "Akanksha_sampleGBS.csv", header = TRUE)
mapdata
View(mapdata)
data(mapdata)
mapdata<-read.cross("csv","D:/R","Akanksha_sampleGBS.csv","estimate.map=FALSE",genotypes=c("a","b"))
class(mapdata)<- "dh"
summary(mapdata)


mapdata1 <- mstmap(mapdata, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-9, detectBadData = TRUE, miss.thresh = 0.7)

summary.map(VEHMap483)

plot.map(vhoptical2)

vhoptical3 <- pull.map(vhoptical2)
tab <- map2table(vhoptical3)
write.csv(tab, "D:/gbsdata/bjoptical/data/vhopticalmap1234.csv")


library(ASMap)
library(qtl) 




data(newmarkers)

newmarkers_1<-read.cross("csv","D:/R","newmarkers.csv","estimate.map=FALSE",genotypes=c("a","b"))
class(newmarkers_1)[1] <- "dh" 
summary(newmarkers_1)
newmarkers_2 <- mstmap(newmarkers_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-10, detectBadData = TRUE, miss.thresh = 0.7)
newmarkers_summary<-summary.map(newmarkers_2)
plot.map(newmarkers_2)
newmarkers_3 <- pull.map(newmarkers_2)
tab <- map2table(newmarkers_3)
write.csv(tab, "D:/R/newmarkers_3.csv")
write.csv(newmarkers_summary, "D:/R/newmarkers_summary.csv")
plotMissing(newmarkers_1)
heatMap(newmarkers_2)

data(arraymarkers)

arraymarkers_1<-read.cross("csvr","D:/R","arraymarkers.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(arraymarkers_1)[1] <- "dh" 
summary(arraymarkers_1)
arraymarkers_2 <- mstmap(arraymarkers_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-14, detectBadData = TRUE, miss.thresh = 0.7)
arraymarkers_summary<-summary.map(arraymarkers_2)
arraymarkers_3 <- pull.map(arraymarkers_2)
tab <- map2table(arraymarkers_3)
write.csv(tab, "D:/R/arraymarkers_3.csv")
write.csv(arraymarkers_summary, "D:/R/arraymarkers_summary.csv")
plot.map(arraymarkers_2)
plotMissing(arraymarkers_1)
cg <- comparegeno(arraymarkers_1)
heatMap(arraymarkers_2)

  
 

  setwd("D:/R/")  
allmarkersdellines_1<-read.cross("csvr","D:/R/","allmarkersdellines.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(allmarkersdellines_1)[1] <- "dh" 
summary(allmarkersdellines_1)
allmarkersdellines_2 <- mstmap(allmarkersdellines_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-14, detectBadData = TRUE, miss.thresh = 0.7)
allmarkersdellines_summary<-summary.map(allmarkersdellines_2)
plot.map(allmarkersdellines_2)
allmarkersdellines_3 <- pull.map(allmarkersdellines_2)
tab <- map2table(allmarkersdellines_3)
write.csv(tab, "D:/R/allmarkersdellines_3.csv")
write.csv(allmarkersdellines_summary, "D:/R/allmarkersdellines_summary.csv")



----------------totalmarkers187---------------------

  
setwd("D:/R/")
library(ASMap)
library(qtl)
totalmarkers187_1<-read.cross("csvr","D:/R/","totalmarkers187.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(totalmarkers187_1)[1] <- "dh" 
summary(totalmarkers187_1)
totalmarkers187_2 <- mstmap(totalmarkers187_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-20, detectBadData = FALSE, miss.thresh = 0.7)
totalmarkers187_summary<-summary.map(totalmarkers187_2)
plot.map(totalmarkers187_2)
totalmarkers187_3 <- pull.map(totalmarkers187_2)
tab <- map2table(totalmarkers187_3)
write.csv(tab, "D:/R/totalmarkers187_3.csv")
write.csv(totalmarkers187_summary, "D:/R/totalmarkers187_summary.csv")


-------------------rice2019new------------------------

  setwd("D:/RR/")
library(ASMap)
library(qtl)
data(rice2019new)


rice2019nodup_1<-read.cross("csvr","D:/RR/data/", "rice2019nodup.csv","estimate.map=FALSE", genotypes=c("AA","BB"))
rice2019nodup_1 <- convert2riself(rice2019nodup_1)
summary(rice2019nodup_1)

dup <- findDupMarkers(rice2019new1, exact.only=TRUE)
snpmap.nodup <- drop.markers(rice2019new1, unlist(dup))
write.cross(snpmap.nodup, "csv", "Data/rice2019nodup")

data(rice2019nodup)
rice2019nodup_1<-read.cross("csv","D:/RR/data/", "rice2019nodup.csv","estimate.map=FALSE", genotypes=c("AA","BB"))
rice2019nodup_1 <- convert2riself(rice2019nodup_1)
summary(rice2019nodup_1)
rice2019nodup_2 <- mstmap(rice2019nodup_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-8, detectBadData = TRUE, miss.thresh = 0.7)
rice2019nodup_summary<-summary.map(rice2019nodup_2)
plot.map(rice2019nodup_2)
rice2019nodup_3 <- pull.map(rice2019nodup_2)
tab <- map2table(rice2019nodup_3)


-----------------------------totalmarkers187---------------------------------

setwd("D:/R/")
library(ASMap)
library(qtl)
data(totalmarkers187)
totalmarkers187_1<-read.cross("csvr","D:/R/data/", "totalmarkers187.csv","estimate.map=FALSE", genotypes=c("A","B"))
class(totalmarkers187_1)[1] <- "dh"
summary(totalmarkers187_1)
dup <- findDupMarkers(totalmarkers187_1, exact.only=TRUE)
nodup <- drop.markers(totalmarkers187_1, unlist(dup))
write.cross(nodup, "csvr", "data/totalmarkersnodup")


data(totalmarkersnodup)
totalmarkersnodup_1<-read.cross("csvr","D:/R/data/", "totalmarkersnodup.csv","estimate.map=FALSE", genotypes=c("AA","BB"))
class(totalmarkersnodup_1)[1] <- "dh"



setwd("D:/R/")
library(ASMap)
library(qtl)
data(Uniqueplusend)
Uniqueplusend_1<-read.cross("csvr","D:/R/data/", "Uniqueplusend.csv","estimate.map=FALSE", genotypes=c("A","B"))
class(Uniqueplusend_1)[1] <- "dh"
summary(Uniqueplusend_1)
dup <- findDupMarkers(Uniqueplusend_1, exact.only=TRUE)
nodup <- drop.markers(Uniqueplusend_1, unlist(dup))
write.cross(nodup, "csvr", "data/Uniqueplusnodup")

data(Uniqueplusnodup)
Uniqueplusnodup_1<-read.cross("csvr","D:/R/data/", "Uniqueplusnodup.csv","estimate.map=FALSE", genotypes=c("AA","BB"))
class(Uniqueplusnodup_1)[1] <- "dh"
summary(Uniqueplusnodup_1)
Uniqueplusnodup_2 <- mstmap(Uniqueplusnodup_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-14, detectBadData = TRUE, miss.thresh = 0.7)
Uniqueplusnodup_summary<-summary.map(Uniqueplusnodup_2)
plot.map(Uniqueplusnodup_2)
Uniqueplusnodup_3 <- pull.map(Uniqueplusnodup_2)
tab <- map2table(Uniqueplusnodup_3)
write.csv(tab, "D:/R/Uniqueplusnodup_3.csv")
write.csv(Uniqueplusnodup_summary, "D:/R/Uniqueplusnodup_summary.csv")




library(qtl)
library(ASMap)
library(xlsx)
library(qtlTools)


setwd("D:/R/")
data(allmarkers)
allmarkers1<-read.cross("csvr","F:/akansha/data", "allmarkers.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(allmarkers1)[1] <- "dh"
summary(allmarkers1)

allmarkers2 <- mstmap(allmarkers1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-15, detectBadData = TRUE, miss.thresh = 0.3)
summary.map(allmarkers2)



allmarkers3 <- pull.map(allmarkers2)
tab <- map2table(allmarkers3)
write.xlsx(tab, "F:/akansha/allmarkersmap_2.xlsx")

write.cross(allmarkers2, "csv", "F:/akansha/veh14012019")


#------------------unique only----14-01-2019-------------------------------------

data(unique)
unique1<-read.cross("csvr","F:/akansha/data", "unique.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(unique1)[1] <- "dh"
summary(unique1)


dup <- findDupMarkers(unique1, exact.only=TRUE)
unique2 <- drop.markers(unique1, unlist(dup))
write.cross(unique2, "csvr", "F:\\akansha\\data/uniquest")


data(uniquest)
uniquest1<-read.cross("csvr","F:/akansha/data", "uniquest.csv","estimate.map=FALSE",genotypes=c("AA","BB"))
class(uniquest1)[1] <- "dh"
summary(uniquest1)



uniquest22 <- mstmap(uniquest1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-15, detectBadData = TRUE, miss.thresh = 0.3)
summary.map(uniquest22)


uniquest33 <- pull.map(uniquest22)
tab <- map2table(uniquest33)
write.xlsx(tab, "F:/akansha/uniquest_map_4.xlsx")

write.cross(uniquest22, "csv", "F:/akansha/vehqtlcross")



#-------------------------LESS THAN Zero RF-----------MARKERS------------------------------------------------------


cross2<-dropSimilarMarkers(uniquest22, rf.threshold = 0.00000001)

summary(cross2)

write.cross(cross2, "csv", "F:/akansha/reduced")



#--------------------totalmarkersnodup------------------------------------

library(qtl)
library(ASMap)
library(writexl)
library(qtlTools)


setwd("D:/R/")
data(totalmarkersnodup)
totalmarkersnodup1<-read.cross("csvr","D:/R/data", "totalmarkersnodup.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(totalmarkersnodup1)[1] <- "dh"
summary(totalmarkersnodup1)
totalmarkersnodup2 <- mstmap(totalmarkersnodup1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-40, detectBadData = TRUE, miss.thresh = 0.3)
totalmarkersnodupsummary<-summary.map(totalmarkersnodup2)
plot.map(totalmarkersnodup2)
totalmarkersnodup3 <- pull.map(totalmarkersnodup2)
tab <- map2table(totalmarkersnodup3)
write.csv(tab, "D:/R/totalmarkersnodup3.csv")
write.csv(totalmarkersnodupsummary, "D:/R/totalmarkersnodupsummary.csv")


#--------------------------------totalmarkersnodup---------------------------------------

library(qtl)
library(ASMap)
library(writexl)
library(qtlTools)


setwd("D:/R/")
data(totalmarkersnodup)
totalmarkersnodup1<-read.cross("csvr","D:/R/data/","totalmarkersnodup.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(totalmarkersnodup1)[1] <- "dh"
summary(totalmarkersnodup1)
tmnclones<-statMark(totalmarkersnodup1, stat.type =c("marker"), map.function = "kosambi")
tmnclones
write.csv(tmnclones, "D:/R/data/tmnclones.csv")
totalmarkersnodup_2 <- mstmap(totalmarkersnodup_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-25, detectBadData = TRUE, miss.thresh = 0.3)
totalmarkersnodup_summary<-summary.map(totalmarkersnodup_2)
plot.map(totalmarkersnodup_2)
totalmarkersnodup_3 <- pull.map(totalmarkersnodup_2)
tab <- map2table(totalmarkersnodup_3)
write.csv(tab, "D:/R/totalmarkersnodup_3.csv")
write.csv(totalmarkersnodup_summary, "D:/R/totalmarkersnodup_summary.csv")



cross2<-dropSimilarMarkers(Uniqueplusnodup_2, rf.threshold = 0.00000001)
jittermap(cross2,amount = 1e-06)
summary(cross2)

write.csv(cross2, "D:/R/reduced.csv")


#--------------------totalmarkers------------------------------------

library(qtl)
library(ASMap)
library(writexl)
library(qtlTools)

setwd("D:/R/")
data(totalmarkers)
totalmarkers1<-read.cross("csvr","D:/R/data/","totalmarkers.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(totalmarkers1)[1] <- "dh"
summary(totalmarkers1)
plotMissing(totalmarkers1)
totalmiss<-statGen(totalmarkers1, bychr = FALSE, stat.type = "miss")
maptotal<-subset(totalmarkers1, ind = totalmiss$miss<3000)
sd<-profileMark(totalmarkers1, stat.type = c("seg.dist"), layout=c(1,1))
mm <- statMark(totalmarkers1, stat.type = c("marker"))
dm <- subsetCross(totalmarkers1)[(mm > 0.98) | (mm < 0.2)]

-------------------------------Rashmi data------------------------------------
setwd("D:/R/")
library(ASMap)
library(qtl)
tdmap_1<-read.cross("csvr","D:/R/","tdmap.csv","estimate.map=FALSE",genotypes=c("a","b"))
class(tdmap_1)[1] <- "dh" 
summary(tdmap_1)
tdmap_2 <- mstmap(tdmap_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-4, detectBadData = FALSE, miss.thresh = 0.7)
tdmap_summary<-summary.map(tdmap_2)
plot.map(tdmap_2)
tdmap_3 <- pull.map(tdmap_2)
tab <- map2table(tdmap_3)
write.csv(tab, "D:/R/tdmap_3.csv")
write.csv(tdmap_summary, "D:/R/tdmap_summary.csv")

#-----------------------------shangmap-------------------------------
setwd("D:/R/data/")
library(ASMap)
library(qtl)
Shangmap_1<-read.cross("csvr","D:/R/data/","Shangmap.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(Shangmap_1)[1] <- "dh" 
summary(Shangmap_1)
Shangmap_2 <- mstmap(Shangmap_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-4, detectBadData = FALSE, miss.thresh = 0.7)
Shangmap_summary<-summary.map(Shangmap_2)
plot.map(Shangmap_2)
Shangmap_3 <- pull.map(Shangmap_2)
tab <- map2table(Shangmap_3)
write.csv(tab, "D:/R/Shangmap_3.csv")
write.csv(Shangmap_summary, "D:/R/Shangmap_summary.csv")



#----------------------------A1-----------------------------------------

setwd("D:/R/data/")
library(ASMap)
library(qtl)
a1_1<-read.cross("csvr","D:/R/data/","a1.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(a1_1)[1] <- "dh" 
summary(a1_1)
a1_2 <- mstmap(a1_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-12, detectBadData = TRUE, miss.thresh = 0.7)
a1_summary<-summary.map(a1_2)
plot.map(a1_2)
a1_3 <- pull.map(a1_2)
tab <- map2table(a1_3)
write.csv(tab, "D:/R/a1_3.csv")
write.csv(a1_summary, "D:/R/a1_summary.csv")


#-----------------Finding Crossovers and Clones----------------------------------

setwd("D:/R/data/")
library(ASMap)
library(qtl)

a1<-read.cross("csvr","D:/R/data/","arraymarker.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(a1)[1] <- "dh" 
summary(a1)
arraymarkers_xo<-countXO(a1, bychr = FALSE)
plot(countXO(a1, bychr = FALSE))
View(arraymarkers_xo)
mean(arraymarkers_xo)
arraymarkers_cl<-genClones(a1)
View(arraymarkers_cl)
a2<-comparegeno(a1, what = c("proportion"))
hist(a2[lower.tri(a2)], breaks = seq(0,1, length(101)), xlab = "proportion of matching genotypes")
a3<-which(a2>0.99)
a4<-pull.geno(a1)
write.csv(arraymarkers_cl[["cgd"]], "D:/R/arrayclones.csv")
write.csv(a2, "D:/R/arrayclones3.csv")

#----------------------------------map100---------------------------------------------

setwd("D:/R/")
library(ASMap)
library(qtl)
map100_1<-read.cross("csvr","D:/R/data/","map100.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(map100_1)[1] <- "dh" 
summary(map100_1)
map100_2 <- mstmap(map100_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-9, detectBadData = TRUE, miss.thresh = 0.7)
map100_summary<-summary.map(map100_2)
plot.map(map100_2)
map100_3 <- pull.map(map100_2)
tab <- map2table(map100_3)
write.csv(tab, "D:/R/map100_3.csv")
write.csv(map100_summary, "D:/R/map100_summary.csv")

#-----------------------------------jittermap-----------------------------------------

setwd("D:/R/")
library(ASMap)
library(qtl)
data("map")
aa<-read.cross("csvr","D:/R/data/","map.csv","estimate.map=TRUE",genotypes=c("A","B"))
class(aa)[1] <- "dh" 
map_1<-jittermap(aa)
map_2<-pull.map(map_1)
tab<-map2table(map_2)
write.csv(tab, "D:/R/jittmap.csv")


#----------------segregation distortion-----------------------

setwd("D:/R/")
library(ASMap)
library(qtl)
data("map")
aa<-read.cross("csvr","D:/R/data/","map.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(aa)[1] <- "dh"
sg<-geno.table(aa)
sg[sg$P.value<0.05/totmar(aa)]
write.csv(sg, "D:/R/sg.csv")

x<-pull.geno(aa)
y<-apply(x, 1, function(a) table(a))
View(y)

#-----------------------------final map-----------------------------------------

setwd("D:/R/")
library(qtl)
library(ASMap)
library(xlsx)

data(vehjuly19)
vehsummary1<-read.cross("csvr","D:/R/data", "vehjuly19.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(vehsummary1)[1] <- "dh"
summary(vehsummary1)
vehsummary2 <- mstmap(vehsummary1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-14, detectBadData = TRUE, miss.thresh = 0.7)
vehsummarysummary<-summary.map(vehsummary2)
plot.map(vehsummary2)
vehsummary3 <- pull.map(vehsummary2)
tab <- map2table(vehsummary3)
write.xlsx(tab, "D:/R/vehmap.xlsx")
write.xlsx(vehsummarysummary, "D:/R/vehmapsummary.xlsx")

vehXO<-countXO(vehsummary1, bychr=TRUE)
write.xlsx(vehXO, "F:/akansha/vehXO.xlsx")

nxo <- countXO(vehsummary1)

plot(nxo, ylab="No. crossovers")
mean(nxo)
nxo[nxo>45]
write.xlsx(nxo, "F:/akansha/vehN-XO.xlsx")

gc1 <- genClones(vehsummary1, tol = 0.99)
gc1$cgd
write.xlsx(gc1$cgd, "F:/akansha/VEH7609_clones.xlsx")



geno.table(vehsummary1)

gt <- geno.table(vehsummary1)
write.xlsx(gt, "F:/akansha/VEH7609_GT.xlsx")

g <- pull.geno(vehsummary1)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))

#---

data(vehnew)
vehnew1<-read.cross("csvr","D:/R/data", "vehnew.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(vehnew1)[1] <- "dh"
summary(vehnew1)

dupveh <- findDupMarkers(vehnew1, exact.only=TRUE)
vehjulunique <- drop.markers(vehnew1, unlist(dupveh))
write.cross(vehjulunique, "csvr", "D:/R/vehjul19uniq")
View(dupveh)

write.xlsx(gfreq, "F:/akansha/VEH7609_Ind_GT.xlsx")


#----------------------vehuniquemap---------------------------

data(vehjul19uniq)
vehuniq1<-read.cross("csvr","D:/R/data", "vehjul19uniq.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(vehuniq1)[1] <- "dh"
summary(vehuniq1)
vehuniq2 <- mstmap(vehuniq1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-12, detectBadData = TRUE, miss.thresh = 0.7)
vehuniqsummary<-summary.map(vehuniq2)
plot.map(vehuniq2)
vehuniq3 <- pull.map(vehuniq2)
tab <- map2table(vehuniq3)
write.xlsx(tab, "D:/R/vehmap2.xlsx")
write.xlsx(vehuniqsummary, "D:/R/vehmapsummary2.xlsx")


setwd("D:/R/")
library(qtl)
library(ASMap)
library(xlsx)
data(vehjul19uniq)
vehuniq_1<-read.cross("csvr","D:/R/data", "vehjul19uniq.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(vehuniq_1)[1] <- "dh"
summary(vehuniq_1)
vehuniq_2 <- mstmap(vehuniq_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-14, detectBadData = TRUE, miss.thresh = 0.7)
vehuniq_summary<-summary.map(vehuniq_2)
plot.map(vehuniq_2)
vehuniq_3 <- pull.map(vehuniq_2)
tab <- map2table(vehuniq_3)
write.xlsx(tab, "D:/R/vehmap3.xlsx")
write.xlsx(vehuniq_summary, "D:/R/vehmapsummary3.xlsx")


#-------------------------------vehjittermapunique----------------------------------

setwd("D:/R/")
library(ASMap)
library(qtl)
data("vehmapunique")
bb<-read.cross("csvr","D:/R/data/","vehmapunique.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(bb)[1] <- "dh" 
map_1<-jittermap(bb, amount = 1e-4)
map_2<-pull.map(map_1)
tab<-map2table(map_2)
write.csv(tab, "D:/R/vehjittmapunique.csv")


#-----------------------------vehuniquest-----------------------------------------

setwd("D:/R/")
library(qtl)
library(ASMap)
library(xlsx)
data(vehuniquest)
vehuniq_1<-read.cross("csvr","D:/R/data", "vehuniquest.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(vehuniq_1)[1] <- "dh"
summary(vehuniq_1)
vehuniq_2 <- mstmap(vehuniq_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-15, detectBadData = TRUE, miss.thresh = 0.7)
vehuniq_summary<-summary.map(vehuniq_2)
plot.map(vehuniq_2)
vehuniq_3 <- pull.map(vehuniq_2)
tab <- map2table(vehuniq_3)
write.xlsx(tab, "D:/R/vehuniquest3.xlsx")
write.xlsx(vehuniq_summary, "D:/R/vehuniquestsummary3.xlsx")


setwd("D:/R/")
library(qtl)
library(ASMap)
library(xlsx)
data(vehunique)
vehuniq_1<-read.cross("csvr","D:/R/data", "vehunique.csv","estimate.map=FALSE",genotypes=c("A","B"))
class(veh_1)[1] <- "dh"
summary(veh_1)
veh_2 <- mstmap(vehuniq_1, bychr = FALSE, dist.fun ="kosambi", p.value = 1e-11, detectBadData = FALSE, miss.thresh = 0.7)
veh_summary<-summary.map(veh_2)
plot.map(veh_2)
vehuniq_3 <- pull.map(veh_2)
tab <- map2table(vehuniq_3)
write.xlsx(tab, "D:/R/marker1461_bd.xlsx")
write.xlsx(veh_summary, "D:/R/marker1461summary_bd.xlsx")

