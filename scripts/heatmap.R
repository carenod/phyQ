#limpio la memoria
rm( list=ls() )  #remove all objects
gc()             #garbage collection
Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(edgeR)
library(pheatmap)

# Scripts for ASPli analysis per SRP
setwd('/data4/projects/Metanalisis_light/Dani/phyQ/')

# plot TMM from edgeR