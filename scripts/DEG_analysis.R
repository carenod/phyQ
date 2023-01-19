#limpio la memoria
rm( list=ls() )  #remove all objects
gc()             #garbage collection
Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(data.table)
library(ggplot2)
library(topGO)

# Scripts for ASPli analysis per SRP
setwd('/data4/projects/Metanalisis_light/Dani/phyQ/')

# Cargo tablas
wt_original <- fread("./R_scripts/phyQ/results/DEG/WTdark_WTred_DEG.txt") # 18114 genes
phyq_original <- fread("./R_scripts/phyQ/results/DEG/phyQdark_phyQred_DEG.txt") # 18506 genes

# filto por fdr
wt <- wt_original[gen.fdr < 0.05, ] # quedan 10507
phyq <- phyq_original[gen.fdr < 0.05, ] # quedan 9517

# hago histograma de logfc
hist(wt$logFC, xlim = c(-2, 2), breaks = 1000)
hist(phyq$logFC, xlim = c(-2, 2), breaks = 1000)  

wt[, .SD[which.min(abs(logFC))]] # 0.2449089
phyq[, .SD[which.min(abs(logFC))]] # -0.2737039
# no hay mucho logfc chicos, asi que no filtro por logfc

# veo el overlap
length(unique(intersect(wt$symbol, phyq$symbol))) # 12296
length(unique(setdiff(wt$symbol, phyq$symbol))) # 2779
length(unique(setdiff(phyq$symbol, wt$symbol))) # 1789
# la mayoria del cambio en expresión de los genes no depende de los PHY

# que pasa si filto por logfc?
wt_logfc <- wt[abs(logFC) > 0.58, ] 
phyq_logfc <- phyq[abs(logFC) > 0.58, ]

# veo el overlap
length(unique(intersect(wt_logfc$symbol, phyq_logfc$symbol))) # 5527 0.5934078
length(unique(setdiff(wt_logfc$symbol, phyq_logfc$symbol))) # 2299 0.2468327
length(unique(setdiff(phyq_logfc$symbol, wt_logfc$symbol))) # 1488 0.1597595
# Si filto por logFC hay menos genes cuya expresion es independiente de PHY

# que pasa si filto por logfc aun mas?
wt_logfc_strict <- wt[abs(logFC) > 1, ] 
phyq_logfc_strict <- phyq[abs(logFC) > 1, ]

# veo el overlap
length(unique(intersect(wt_logfc_strict$symbol, phyq_logfc_strict$symbol))) # 2547 0.5385824
length(unique(setdiff(wt_logfc_strict$symbol, phyq_logfc_strict$symbol))) # 1743 0.3422344
length(unique(setdiff(phyq_logfc_strict$symbol, wt_logfc_strict$symbol))) # 803 0.1576674
# Si filto con un logfc estricto aumento el porcentaje de genes

# que pasa si filtro por fdr 0.01 y logfc estrico
wt_strict <- wt[gen.fdr < 0.01, ] 
phyq_strict <- phyq[gen.fdr < 0.01, ] #
wt_strict <- wt_strict[abs(logFC) > 1, ] 
phyq_strict <- phyq_strict[abs(logFC) > 1, ]

# veo el overlap
length(unique(intersect(wt_strict$symbol, phyq_strict$symbol))) # 2473 0.502336
length(unique(setdiff(wt_strict$symbol, phyq_strict$symbol))) # 1689 0.3430835
length(unique(setdiff(phyq_strict$symbol, wt_strict$symbol))) # 761 0.1545805
# Si filto por fdr y logfc estricto es parecido a solo por logfc estricto

##### dependientes PHY vs logFC threshold #####
# grafico #degs que dependen de PHYs vs umbral logfc
# Para genes con fdr < 0.05

phyq_dependent <- data.table("threshold" = double(),    # Create an empty data.table
                             "independent" = character(),
                             "dependent" = double(),
                             "total" = double())

for (th in seq(0, 10, 0.01)) {
  wt_tmp <- wt[abs(logFC) > th, ] 
  phyq_tmp <- phyq[abs(logFC) > th, ]
  independent <- length(unique(intersect(wt_tmp$symbol, phyq_tmp$symbol)))
  dependent_in_wt <- length(unique(setdiff(wt_tmp$symbol, phyq_tmp$symbol))) 
  dependent_in_phyq <- length(unique(setdiff(phyq_tmp$symbol, wt_tmp$symbol))) 
  dependent <- dependent_in_wt + dependent_in_phyq
  total <- independent + dependent
  independent <- independent / total * 100
  dependent <- dependent / total * 100
  new_row <- list(th, independent, dependent, total)
  phyq_dependent <- rbind(phyq_dependent, new_row)
  rm(new_row, independent, dependent, th)
}

ggplot(data = phyq_dependent, aes(x= threshold, y = dependent)) + 
  geom_line()

##### Hago GO terms de los dependientes y los independientes #####
# lo hago sobre los genes con fdr < 0.05 y logfc > 0.58

independent_genes <- intersect(wt_logfc$symbol, phyq_logfc$symbol)
dependent_genes <- setdiff(wt_logfc$symbol, phyq_logfc$symbol)
dependent_genes <- c(dependent_genes, setdiff(phyq_logfc$symbol, wt_logfc$symbol))

all_genes <- unique(intersect(wt_original$symbol, phyq_original$symbol))
unique_genes <- unique(data$gene)
# Binary table with spliced genes

geneList_dependent <- ifelse(all_genes %in% dependent_genes, 1, 0)
names(geneList_dependent) <- all_genes

geneList_independent <- ifelse(all_genes %in% independent_genes, 1, 0)
names(geneList_independent) <- all_genes

AtMeta_GO <- new("topGOdata",
                 ontology = "BP",
                 allGenes = geneList_independent,
                 geneSelectionFun = function(x)(x == 1),
                 annot = annFUN.org, mapping = "org.At.tair.db")

# Fisher test to determine significance
resultFisher <- runTest(AtMeta_GO, algorithm = "elim", statistic = "fisher")
topGO <- GenTable(AtMeta_GO, raw.p.value = resultFisher, 
                  topNodes = length(resultFisher@score),
                  numChar = 120)
# Enrichment = ( significant genes / all (spliced?) genes ) 
#                ________________________________________________
#              ( annotated genes / protein-coding genes in genome ) 

topGO$enrichment <- (topGO$Significant / length(unique_genes) ) / (topGO$Annotated / 27416) 
topGO_results <- transform(topGO, raw.p.value = as.numeric(raw.p.value))

