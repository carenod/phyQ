#limpio la memoria
rm( list=ls() )  #remove all objects
gc()             #garbage collection
Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(ASpli)
library(GenomicFeatures)
library(data.table)
library(edgeR)

# Scripts for ASPli analysis per SRP
setwd('/data4/projects/Metanalisis_light/Dani/phyQ/')
source("./R_scripts/phyQ/scripts/ASpli_pipeline_contrast_strict.R")
source("./R_scripts/phyQ/scripts/ASpli_contrast.R")

# Load genome annotation
# ATxDb <- makeTxDbFromGFF("/data/BioData/Arabidopsis_thaliana/Ensembl/TAIR10/Annotation/Archives/2019/TAIR10_GFF3_genes.gff")
# saveDb(ATxDb,file="gene.sqlite")
ATxDb <- loadDb("gene.sqlite")
features <- binGenome(ATxDb)

# Set library parameters
libType <- "PE"
stranMode <- 0
minReadLength <- 150L

# Search bam files
bamFiles <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/02_BAM', pattern = "out.bam$" ,
                              full.names = TRUE, recursive = TRUE)

targets <- data.frame(bam = bamFiles,
                      genotype = c('WT','WT','WT',
                                   'phyQ','phyQ','phyQ'),
                      light = rep(c('dark', 'red'), each = 6),
                      stringsAsFactors = FALSE)

getConditions(targets)
# "WT_dark"   "phyQ_dark" "WT_red"    "phyQ_red" 

MergebamFileNames <- list.files(path = '/data4/projects/ASpliDB/Data/SRPnew/03_mergedBAMs', 
                                pattern = ".bam$" ,
                                full.names = TRUE, 
                                recursive = TRUE)


mBAMs <- data.frame(bam = MergebamFileNames[c(1,3,2,4)],
                    condition = getConditions(targets))


##### Unnormalized read count #####

gb_all <- gbCounts(features, targets, 
                   minReadLength = minReadLength,
                   maxISize = 50000, 
                   minAnchor = 5, 
                   libType= libType,
                   strandMode= stranMode, 
                   alignFastq = FALSE, 
                   dropBAM = FALSE)

asd_all <- jCounts(counts = gb_all, 
                   features = features, 
                   minReadLength = minReadLength, 
                   libType= libType,
                   threshold = 5,                 
                   minAnchor = 10,
                   strandMode = strandMode)

gene_counts <- countsg(gb_all)
bin_counts <- countsb(gb_all)
fwrite(gene_counts, "./counts/raw/gene_counts")
fwrite(bin_counts, "./counts/raw/bin_counts")

##### plot MDS ####

# filter by expression
group <- c(1,1,1,2,2,2,3,3,3,4,4,4)
x <- gene_counts[8:19]
x <- bin_counts[10:21]
y <- DGEList(counts = x, group = group)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]

# normalize
y <- calcNormFactors(y)
y$samples

# plot
plotMDS(y, top = 2000)

##### WT dark vs WT red #####
# WTdark_WTred <- ASpli_pipeline_contrast_strict(targets = targets[c(1,2,3,7,8,9),], 
#                                                minReadLength = 150, 
#                                                libType = 'PE', 
#                                                sra = 'WTdark_WTred', 
#                                                features = features,
#                                                mBAMs = mBAMs[c(1,3),],
#                                                strandMode = 0,
#                                                contrast = c(1,-1))

WTdark_WTred <- ASpli_pipeline_contrast_strict(contrast = c(1, 0, -1, 0),
                                               coef = NULL,
                                               formula = NULL,
                                               sra = 'WTdark_WTred', 
                                               features = features, 
                                               mBAMs = mBAMs[c(1,3),], 
                                               gbcounts = gb_all, 
                                               asd = as_all)

# Guardo DEGs
WTdark_WTred_DEG <- genesDE(WTdark_WTred$gbDUreport)
fwrite(WTdark_WTred_DEG, "./R_scripts/phyQ/results/DEG/WTdark_WTred_DEG.txt")

##### phyQ dark vs phyQ red #####
# phyQdark_phyQred <- ASpli_pipeline_contrast_strict(targets = targets[c(4,5,6,10,11,12),], 
#                                                minReadLength = 150, 
#                                                libType = 'PE', 
#                                                sra = 'phyQdark_phyQred', 
#                                                features = features,
#                                                mBAMs = mBAMs[c(2,4),],
#                                                strandMode = 0,
#                                                contrast = c(1,-1))

phyQdark_phyQred <- ASpli_pipeline_contrast_strict(contrast = c(0, 1, 0, -1),
                                                   coef = NULL,
                                                   formula = NULL,
                                                   sra = 'phyQdark_phyQred', 
                                                   features = features, 
                                                   mBAMs = mBAMs[c(2,4),], 
                                                   gbcounts = gb_all, 
                                                   asd = as_all)

# Guardo DEGs
phyQdark_phyQred_DEG <- genesDE(phyQdark_phyQred$gbDUreport)
fwrite(phyQdark_phyQred_DEG, "./R_scripts/phyQ/results/DEG/phyQdark_phyQred_DEG.txt")

##### WT dark vs phyQ dark #####
# WTdark_phyQdark <- ASpli_pipeline_contrast_strict(targets = targets[1:6,], 
#                                                   minReadLength = 150, 
#                                                   libType = 'PE', 
#                                                   sra = 'WTdark_phyQdark', 
#                                                   features = features,
#                                                   mBAMs = mBAMs[c(1,2),],
#                                                   strandMode = 0,
#                                                   contrast = c(1,-1))

WTdark_phyQdark <- ASpli_pipeline_contrast_strict(contrast = c(1, -1, 0, 0),
                                                  coef = NULL,
                                                  formula = NULL,
                                                  sra = 'WTdark_phyQdark', 
                                                  features = features, 
                                                  mBAMs = mBAMs[c(1,2),], 
                                                  gbcounts = gb_all, 
                                                  asd = as_all)

# Guardo DEGs
WTdark_phyQdark_DEG <- genesDE(WTdark_phyQdark$gbDUreport)
fwrite(WTdark_phyQdark_DEG, "./R_scripts/phyQ/results/DEG/WTdark_phyQdark_DEG.txt")

##### WT red vs phyQ red #####
# WTred_phyQred <- ASpli_pipeline_contrast_strict(targets = targets[7:12,], 
#                                                 minReadLength = 150, 
#                                                 libType = 'PE', 
#                                                 sra = 'WTred_phyQred', 
#                                                 features = features,
#                                                 mBAMs = mBAMs[c(3,4),],
#                                                 strandMode = 0,
#                                                 contrast = c(1,-1))

WTred_phyQred <- ASpli_pipeline_contrast_strict(contrast = c(0, 0, 1, -1),
                                                coef = NULL,
                                                formula = NULL,
                                                sra = 'WTred_phyQred', 
                                                features = features, 
                                                mBAMs = mBAMs[c(3,4),], 
                                                gbcounts = gb_all, 
                                                asd = as_all)

# Guardo DEGs
WTred_phyQred_DEG <- genesDE(WTred_phyQred$gbDUreport)
fwrite(WTred_phyQred_DEG, "./R_scripts/phyQ/results/DEG/WTred_phyQred_DEG.txt")


