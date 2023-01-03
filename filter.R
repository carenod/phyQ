#limpio la memoria
rm( list=ls() )
gc()

Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")

require(data.table)
require(biomaRt)
require(seqinr)
require(stringr)
require(made4)
require(seqTools)

setwd("/data4/projects/Metanalisis_light/ASpli_reports_STRICT_eline/SRP_new/")

##### Anchor based #####
abr <- fread("./sr_SRP058008/anchorbased.txt")

colnames(abr)

# para evitar que confunda un corte endonucleolitico con splicing
# me quedo con junturas donde J1 y J3 tienen al menos el 60% de J3
# abr1 <- abr[ countsJ1.F_treated > 0.6 * countsJ3.F_control & countsJ2.F_treated >  0.6 * countsJ3.F_control, ]
# abr2 <- abr[ countsJ1.F_control > 0.6 * countsJ3.F_treated & countsJ2.F_control >  0.6 * countsJ3.F_treated, ]
# abr3 <- abr[ countsJ3.F_treated > 0.6 * countsJ1.F_control & countsJ3.F_treated >  0.6 * countsJ2.F_control, ]
# abr4 <- abr[ countsJ3.F_control > 0.6 * countsJ1.F_treated & countsJ3.F_control >  0.6 * countsJ2.F_treated, ]

# abr_concat <- rbindlist(list(abr1, abr2, abr3, abr4))  

# Filtro para quedarme con junturas no anotadas significativas
#abr <- abr[junction.annotated == "No", ]
abr <- abr[junction.fdr < 0.05 & abs(junction.dPIR) > 0.05, ]
abr <- unique(abr)

fwrite(abr, "./sr_SRP058008/anchorbased_filtered.txt", dec = ",", sep = "\t")

##### Bin based #####
bbr <- fread("./sr_SRP058008/binbased.txt")
#colnames(bbr)

# Filtro para quedarme con significativas a nivel coverage y
# dPIR/PSI
bbr <- bbr[bin.fdr < 0.1 & cluster.fdr < 0.1 & abs(bin.logFC) > 0.58 & 
             (abs(cluster.dPIR) > 0.05 | abs(cluster.dPSI) > 0.05),]

fwrite(bbr, "./sr_SRP058008/binbased_filtered.txt", dec = ",", sep = "\t")

##### Locale based #####

lbr <- fread("./sr_SRP058008/localebased.txt") # 1394 rows
is <- fread("./SRP_new_contrastWT/is_SRPnew_WT/D-LP/D - LP.csv")

# arreglo nombre
colnames(is) <- gsub(" ", "_", colnames(is))

is <- is[(Bin_Evidence == 1 & Bin_SJ_Evidence == 1) |
           #(Bin_Evidence == 1 & logFC > 1) |
           Anchor_Evidence == 1 | Locale_Evidence == 1,]

colnames(lbr)
colnames(is)

lbr <- lbr[cluster.fdr < 0.05, ] # 389 rows



# en lbr no hay dPSI
# dPSI aparece solo en is
isl <- is[Locale_Evidence == 1,] 
isl <- isl[abs(Participation - dParticipation) > 0.04, ]

# agrego dPSI
lbr <- lbr[cluster.locus %in% isl$Locus, ] # 191 rows
lbr[isl, dPSI := dParticipation, on= c(cluster.locus = "Locus")]

fwrite(lbr, "./sr_SRP058008/localebased_filtered.txt", dec = ",", sep = "\t")

##### Lista de todos los genes ####
#blue_as <- c(abr$cluster.locus, bbr$locus, lbr$cluster.locus)
blue_as <- unique(is$Locus)


# todos los genes que se expresan
geneX <- fread("/data4/projects/Metanalisis_light/Dani/phyQ/.DEG_light_vs_dark_WT.txt")

genes_blue <- data.table(AGI = geneX$symbol)
genes_blue[, AS := fifelse(AGI %in% blue_as, 1, 0)]

setdiff(blue_as, genes_blue$AGI)
unique(blue_as)

##### uso binary matrix de eline #####
binary <- fread("/data4/projects/Metanalisis_light/Dani/Clustering/binary_matrix.csv")
colnames(binary)[1] <- "AGI"

##### Bajo secuencias de genes #####

# genero mart
phytozome_v13 <- useMart(biomart = "phytozome_mart", 
                         dataset = "phytozome", 
                         host = "https://phytozome-next.jgi.doe.gov") 

# busco las secuencias
genes_seq <- getBM(attributes = c("gene_name1", "gene_exon_intron"), 
                  filters = c("organism_id"), #"gene_name_filter"), 
                  values = list(organism_id = 167), #, gene_name_filter = binary$V1), 
                  mart = phytozome_v13)

# chequeo que esten todas las secuencias

length(unique(genes_seq$gene_name1)) # hay 17603 secuencias cuando en pedi 18114

diferencias <- comparelists(genes_blue$AGI, genes_seq$gene_name1, Set.Diff)
diferencias$Set.Diff

dif_as <- comparelists(diferencias$Set.Diff, blue_as)
dif_as$intersect # solo me pìerdo un gen


##### genero kmers ####

 
for (i in 1:nrow(genes_seq)) {

  s <- genes_seq$gene_exon_intron[i]
  n <- genes_seq$gene_name1[i]

  # si hay un caracter que no es atcg lo elimino
  s <- gsub("[^ATGC]", "", s)
  tmp <- countDnaKmers(s, k = 6, 1, nchar(s)-5)
  
  if (i == 1) {
    # para la primera iteracion genero un dt nuevo
    gene_kmer <- data.table(kmer = rownames(tmp),
                            counts = tmp)
    # pongo el gen como nombre de la columna
    setnames(gene_kmer, "counts.1", n)


  } else {
    # para las siguientes filas
    tmp <- data.table(counts = tmp, 
                      kmer = rownames(tmp))
    
    gene_kmer[tmp, counts := counts.1, on="kmer"]
    
    # pongo el gen como nombre de la columna
    setnames(gene_kmer, "counts", n)

  }

}

gene_kmer <- transpose(gene_kmer, keep.names = "AGI", make.names = "kmer")

# agrego si se splicea o no
#gene_kmer[binary, AS := AS, on="AGI"]

genes_binary <- binary$AGI
gene_kmer[, AS := ifelse(AGI %in% genes_binary, 1, 0)]


# guardo tabla final
fwrite(gene_kmer, "from_binary.txt")
