library(stringr)

# Este script analisa los valores de un gen en todas las tablas de ASpli

# Selecciono gen y experimento
experiment <- WTred_phyQred
goi <- 'AT1G66800'
boi <- 'AT1G66800:I003' # bin of interest

# Gene unnormalized read counts
subset(countsg(experiment$gbcounts), symbol == goi)
# symbol locus_overlap       gene_coordinates    start      end length effective_length   1   2   3   7   8   9
# AT1G66800 AT1G66800             - Chr1:24924813-24926365 24924813 24926365   1553             1128 306 296 153 325 392 309

# Bin unnormalized read counts
subset(countsb(experiment$gbcounts), symbol == goi)
# AT1G66800:I003       I        - AT1G66800             - AT1G66800 Chr1:24924813-24926365 24925567 24925664     98  36   7
# AT1G66800:I003  0  34  33  33

# Junture unnormalized read counts
subset(countsj(experiment$gbcounts), symbol == goi)
# Chr1.24925566.24925665 AT1G66800:E003;AT1G66800:E004            - 64  89 35  50  76  53

subset(countse1i(experiment$gbcounts), symbol == goi)
# 11  5  0 22 24 28

subset(countsie2(experiment$gbcounts), symbol == goi)
#  9 17 0 31 42 30

tmp <- altPSI(experiment$asd)
subset(tmp, boi %in% rownames(tmp))
# nada

# Annotation based: reads and PIR
tmp <- irPIR(experiment$asd)
tmp[rownames(tmp) == boi, ]
# muestra J1, 2 y 3

# Annotation free: reads and PIR
tmp <- junctionsPIR(experiment$asd)
View(subset(tmp, hitIntron == boi))

# Annotation based: bin usage
subset(binsDU(experiment$gbDUreport), symbol == goi)

# Annotated juncture test
tmp <- jir(experiment$jDUreport)
tmp[rownames(tmp) == boi, ]

# Annotation free: Locale
tmp <- localej(experiment$jDUreport)
tmp[str_detect(rownames(tmp), "Chr1.2492"), ]

tmp <- localec(experiment$jDUreport)
subset(tmp, str_detect(range, "Chr1.2492"))

# Annotation free: Anchorage
tmp <- anchorj(experiment$jDUreport)
tmp[str_detect(rownames(tmp), "Chr1.2492"), ]

tmp <- anchorc(experiment$jDUreport)
tmp[str_detect(rownames(tmp), "Chr1.2492"), ]

