setwd("D:/R/data/")
library(corrplot)
library(xlsx)
heterozyg1<-read.csv(file="genomehetero.csv", header = TRUE, row.names = 1)
M <- cor(heterozyg1, use = "pairwise.complete.obs")

corrplot(M, method = "number")

res111 <- cor.mtest(heterozyg1, conf.level = .95)
corrplot(M, method = "number", p.mat = res111$p, insig = "label_sig",sig.level = .01, pch.cex = .9, pch.col = "white")

write.xlsx(M, "D:/R/M.xlsx")


#------------------------------alltraitscorr------------------------------

setwd("D:/R/data/")
library(corrplot)
library(xlsx)
H<-read.csv(file="alltraitscorr.csv", header = TRUE, row.names = 1)
T <- cor(H, use = "pairwise.complete.obs")

corrplot(T, method = "number")

res111 <- cor.mtest(H, conf.level = .95)
corrplot(T, method = "color", p.mat = res111$p, insig = "label_sig",sig.level = .01, pch.cex = .9, pch.col = "white")

write.xlsx(T, "D:/R/T.xlsx")


#------------------------------plhtcorrallpop---------------------------------


setwd("D:/R/data/")
library(corrplot)
library(xlsx)
P<-read.csv(file="plhtcorrallpop.csv", header = TRUE, row.names = 1)
N <- cor(P, use = "pairwise.complete.obs")

corrplot(N, method = "number")

res111 <- cor.mtest(H, conf.level = .95)
corrplot(N, method = "color", p.mat = res111$p, insig = "label_sig",sig.level = .01, pch.cex = .9, pch.col = "white")

write.xlsx(N, "D:/R/N.xlsx")
