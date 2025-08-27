
writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")

setwd("D:/METAN/")
library(metan)
model <- gamem_met(data_ge,
                   env = ENV,
                   gen = GEN,
                   rep = REP,
                   resp = everything())

# Genotype BLUE
gmd(model, "blueg")

# Genotype BLUP
gmd(model, "blupg")

# Genotype-environment BLUE
gmd(model, "bluege")

# Genotype-environment BLUP
gmd(model, "blupge")


#You can also compute the BLUPs or BLUEs for for each environment with


model2 <- gamem(data_ge,
                gen = GEN,
                rep = REP,
                resp = everything(),
                by = ENV)

# BLUE for genotypes in each environment
gmd(model2, "blueg")

# BLUP for genotypes in each environment
gmd(model2, "blupg")



#-------------------DF------------------------
setwd("D:/METAN/")
library(metan)
DF <- read.csv("DF.csv")
model <- gamem_met(DF,
                   env = ENV,
                   gen = GEN,
                   rep = REP,
                   resp = everything())

# Genotype BLUE
gmd(model, "blueg")

# Genotype BLUP
gmd(model, "blupg")

# Genotype-environment BLUE
gmd(model, "bluege")

# Genotype-environment BLUP
gmd(model, "blupge")