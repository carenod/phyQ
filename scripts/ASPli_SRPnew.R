#limpio la memoria
rm( list=ls() )  #remove all objects
gc()             #garbage collection
Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(ASpli)
library(GenomicFeatures)

# Scripts for ASPli analysis per SRP
setwd('/data4/projects/Metanalisis_light/Dani/phyQ/')
source("./R_scripts/phyQ/scripts/ASpli_pipeline_contrast_strict.R")

# Load genome annotation
ATxDb <- makeTxDbFromGFF("/data/BioData/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/2019/TAIR10_GFF3_genes.gff")
features <- binGenome(ATxDb)
saveDb(ATxDb,file="gene.sqlite")

# Search bam files
bamFiles <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets <- data.frame(bam = bamFiles,
                      genotype = c('WT','WT','WT',
                                   'phyQ','phyQ','phyQ'),
                      light = rep(c('dark', 'red'), each = 6),
                      stringsAsFactors = FALSE)

MergebamFileNames <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', 
                                pattern = ".bam$" ,
                                full.names = TRUE, 
                                recursive = TRUE)


mBAMs <- data.frame(bam = MergebamFileNames[c(1,3,2,4)],
                    condition = getConditions(targets))


##### WT dark vs WT red #####
WTdark_WTred <- ASpli_pipeline_contrast_strict(targets = targets[c(1,2,3,7,8,9),], 
                                               minReadLength = 150, 
                                               libType = 'PE', 
                                               sra = 'WTdark_WTred', 
                                               features = features,
                                               mBAMs = mBAMs[c(1,3),],
                                               strandMode = 0,
                                               contrast = c(1,-1))

# Guardo DEGs
genesDE(WTdark_WTred)


##### phyQ dark vs phyQ red #####
phyQdark_phyQred <- ASpli_pipeline_contrast_strict(targets = targets[c(4,5,6,10,11,12),], 
                                               minReadLength = 150, 
                                               libType = 'PE', 
                                               sra = 'phyQdark_phyQred', 
                                               features = features,
                                               mBAMs = mBAMs[c(2,4),],
                                               strandMode = 0,
                                               contrast = c(1,-1))

##### WT dark vs phyQ dark #####
WTdark_phyQdark <- ASpli_pipeline_contrast_strict(targets = targets[1:6,], 
                                                  minReadLength = 150, 
                                                  libType = 'PE', 
                                                  sra = 'WTdark_phyQdark', 
                                                  features = features,
                                                  mBAMs = mBAMs[c(1,2),],
                                                  strandMode = 0,
                                                  contrast = c(1,-1))

##### WT red vs phyQ red #####
WTred_phyQred <- ASpli_pipeline_contrast_strict(targets = targets[7:12,], 
                                                minReadLength = 150, 
                                                libType = 'PE', 
                                                sra = 'WTred_phyQred', 
                                                features = features,
                                                mBAMs = mBAMs[c(3,4),],
                                                strandMode = 0,
                                                contrast = c(1,-1))

