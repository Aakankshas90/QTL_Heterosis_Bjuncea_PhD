install.packages("shiny")  
install.packages("circlize")  
install.packages("RColorBrewer")
install.packages("data.table")
install.packages("RLumShiny")  
## try http:// if https:// URLs are not supported  
source("https://bioconductor.org/biocLite.R")  
biocLite("GenomicRanges")
BiocManager::install(c("GenomicFeatures", "AnnotationDbi"))
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicRanges", version = "3.8")
a
BiocInstaller::biocVersion() # recent BiocInstaller
}, error=function(...) {         # no / older BiocInstaller
  BioC_version_associated_with_R_version <-
    get(".BioC_version_associated_with_R_version",
        envir=asNamespace("tools"), inherits=FALSE)
  if (is.function(BioC_version_associated_with_R_version))
    BioC_version_associated_with_R_version()
  else                            # numeric_version
    BioC_version_associated_with_R_version


library(shiny)
runApp("D:/apps/shinyCircos-master", launch.browser = TRUE)