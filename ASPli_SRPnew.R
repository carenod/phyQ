Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(ASpli)
library(GenomicFeatures)

# Scripts for ASPli analysis per SRP
setwd('/data4/projects/Metanalisis_light/Dani/phyQ/')

# Load genome annotation
ATxDb <- makeTxDbFromGFF("/data/BioData/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/2019/TAIR10_GFF3_genes.gff")
features <- binGenome(ATxDb)
saveDb(ATxDb,file="gene.sqlite")

########### Photoperiod and light pulse #####
##### SRPnew contrast WT vs phyQ Dark #####
# New experiment, performed by Connie
#	DC: A, B, C
# LP: G, H, I

bamFiles_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets_SRPnew <- data.frame(bam = c(bamFiles_SRPnew[c(1:6)]),
                             condition = c('WT','WT','WT',
                                           'phyQ','phyQ','phyQ'),
                             stringsAsFactors = FALSE)

MergebamFileNames_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', pattern = ".bam$" ,
                                       full.names = TRUE, recursive = TRUE)


mBAMs_SRPnew <- data.frame(bam = c(MergebamFileNames_SRPnew[c(1,3)]), 
                           condition = getConditions(targets_SRPnew))

signals_SRPnew_WTvsphyQ_Dark <- ASpli_pipeline_contrast_strict(targets = targets_SRPnew, 
                                                          minReadLength = 150, 
                                                          libType = 'PE', 
                                                          sra = 'SRPnew_WT', 
                                                          features = features,
                                                          mBAMs = mBAMs_SRPnew,
                                                          strandMode = 0,
                                                          form = NULL,
                                                          contrast = c(1,-1))





##### SRPnew contrast phyQ #####
##### SRPnew contrast WT vs phyQ Light #####
#	DC: D, E, F
# LP: J, K, L

bamFiles_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets_SRPnew <- data.frame(bam = c(bamFiles_SRPnew[c(7:12)]),
                             condition = c('WT','WT','WT',
                                           'phyQ','phyQ','phyQ'),
                             stringsAsFactors = FALSE)

MergebamFileNames_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', pattern = ".bam$" ,
                                       full.names = TRUE, recursive = TRUE)


mBAMs_SRPnew <- data.frame(bam = c(MergebamFileNames_SRPnew[c(2,4)]), 
                           condition = getConditions(targets_SRPnew))

signals_SRPnew_DvsLP_phyQ <- ASpli_pipeline_contrast_strict(targets = targets_SRPnew, 
                                                            minReadLength = 150, 
                                                            libType = 'PE', 
                                                            sra = 'SRPnew_PhyQ', 
                                                            features = features,
                                                            mBAMs = mBAMs_SRPnew,
                                                            strandMode = 0,
                                                            form = NULL,
                                                            contrast = c(1,-1))


########### Photoperiod and light pulse #####
##### SRPnew contrast WT light vs WT Dark #####
# New experiment, performed by Connie
#	DC: A, B, C
# LP: G, H, I

bamFiles_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets_SRPnew <- data.frame(bam = c(bamFiles_SRPnew[c(1:3,7:9)]),
                             condition = c('dark','dark','dark',
                                           'light','light','light'),
                             stringsAsFactors = FALSE)

MergebamFileNames_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', pattern = ".bam$" ,
                                       full.names = TRUE, recursive = TRUE)


mBAMs_SRPnew <- data.frame(bam = c(MergebamFileNames_SRPnew[c(1,2)]), 
                           condition = getConditions(targets_SRPnew))

signals_SRPnew_WTvsphyQ_Dark <- ASpli_pipeline_contrast_strict(targets = targets_SRPnew, 
                                                               minReadLength = 150, 
                                                               libType = 'PE', 
                                                               sra = 'SRPnew_WT', 
                                                               features = features,
                                                               mBAMs = mBAMs_SRPnew,
                                                               strandMode = 0,
                                                               form = NULL,
                                                               contrast = c(1,-1))





##### SRPnew contrast phyQ #####
##### SRPnew contrast phyQ light vs phyQ Dark #####
#	DC: D, E, F
# LP: J, K, L

bamFiles_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets_SRPnew <- data.frame(bam = c(bamFiles_SRPnew[c(4:6,10:12)]),
                             condition = c('dark','dark','dark', 
                                           'light','light','light'),
                             stringsAsFactors = FALSE)

MergebamFileNames_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', pattern = ".bam$" ,
                                       full.names = TRUE, recursive = TRUE)


mBAMs_SRPnew <- data.frame(bam = c(MergebamFileNames_SRPnew[c(3,4)]), 
                           condition = getConditions(targets_SRPnew))

signals_SRPnew_DvsLP_phyQ <- ASpli_pipeline_contrast_strict(targets = targets_SRPnew, 
                                                            minReadLength = 150, 
                                                            libType = 'PE', 
                                                            sra = 'SRPnew_PhyQ', 
                                                            features = features,
                                                            mBAMs = mBAMs_SRPnew,
                                                            strandMode = 0,
                                                            form = NULL,
                                                            contrast = c(1,-1))

##### SRPnew formula phyQ #####
#	DC: D, E, F
# LP: J, K, L

bamFiles_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets_SRPnew <- data.frame(bam = bamFiles_SRPnew,
                             genotype = c('dark','dark','dark', 'dark','dark','dark',
                                           'light','light','light', 'light','light','light'),
                             light = c('WT', 'WT', 'WT', 'phyQ', 'phyQ', 'phyQ',
                                       'WT', 'WT', 'WT', 'phyQ', 'phyQ', 'phyQ'),
                             stringsAsFactors = FALSE)

MergebamFileNames_SRPnew <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', pattern = ".bam$" ,
                                       full.names = TRUE, recursive = TRUE)


mBAMs_SRPnew <- data.frame(bam = MergebamFileNames_SRPnew, 
                           condition = c('dark_WT', 'light_WT', 'dark_phyQ', 'light_phyQ'))

formula <- formula(~genotype+light+genotype:light)

signals_formula_2 <- ASpli_pipeline_formula_strict(targets = targets_SRPnew, 
                                                            minReadLength = 150, 
                                                            libType = 'PE', 
                                                            sra = 'SRPnew_PhyQ', 
                                                            features = features,
                                                            mBAMs = mBAMs_SRPnew,
                                                            strandMode = 0,
                                                            form = formula)

signals_formula_3 <- ASpli_pipeline_formula_strict(targets = targets_SRPnew, 
                                                           minReadLength = 150, 
                                                           libType = 'PE', 
                                                           sra = 'SRPnew_PhyQ', 
                                                           features = features,
                                                           mBAMs = mBAMs_SRPnew,
                                                           strandMode = 0,
                                                           form = formula,
                                                           coef = 3)

signals_formula_4 <- ASpli_pipeline_formula_strict(targets = targets_SRPnew, 
                                                           minReadLength = 150, 
                                                           libType = 'PE', 
                                                           sra = 'SRPnew_PhyQ', 
                                                           features = features,
                                                           mBAMs = mBAMs_SRPnew,
                                                           strandMode = 0,
                                                           form = formula,
                                                           coef = 4)

ASpli_pipeline_contrast_Dani(targets_SRPnew, 
                             minReadLength = 150, 
                             libType = 'PE', 
                             contrast = ,
                                                                      form, sra, features, mBAMs, strandMode)
