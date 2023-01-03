##### ASpli pipeline ######
# ASpli_pipeline function to execute gbcounts, asd, and return 
ASpli_pipeline_formula_strict <- function(targets, minReadLength, libType,  
                                          form, sra, features, mBAMs, strandMode) {
  
  print('Starting read counting')
  
  gbcounts <- gbCounts(features = features, 
                       targets = targets,
                       minReadLength = minReadLength, 
                       maxISize = 50000, 
                       libType = libType,
                       strandMode = strandMode)    
  
  print('Calculating PSI, PIR and PJU')
  
  # Junction-based de-novo counting and splicing signal estimation:
  #  asd (alternative splicing d?) counts junctions (jCounts function), which are reads overlapping bin borders 
  #  (e.g. 1 bp in an exon and 99 in the following intron)
  asd <- jCounts(counts = gbcounts, 
                 features = features, 
                 minReadLength = minReadLength, 
                 libType= libType,
                 threshold = 10, # paso threshold de 5 a 10 (STRICT)                 
                 minAnchor = 10,
                 strandMode = strandMode)
  
  results <- c()
  
  for (i in 2:4) {
    
    print(paste0('Doing gbDUreport ', i))
    # Differential gene expression and bin usage signal estimation
    gbDUreport <- gbDUreport(gbcounts, 
                             minGenReads = 100, # paso de 10 a 100 (STRICT)
                             minBinReads = 50,  # paso de 5 a 50 (STRICT)
                             minRds = 0.5, # paso de .05 a 0.5 (STRICT)
                             contrast = NULL, 
                             formula = form, 
                             coef = i)
    
    print(paste0('Doing jDUreport ', i))
    # Differential junction usage analysis
    jDUreport <- jDUreport(asd, 
                           minAvgCounts = 50, # paso de 5 a 50 (STRICT)
                           filterWithContrasted = TRUE,
                           runUniformityTest = TRUE,
                           mergedBams = mBAMs,
                           maxPValForUniformityCheck = 0.2, 
                           strongFilter = TRUE,
                           maxConditionsForDispersionEstimate = 24,
                           contrast = NULL,
                           formula = form,
                           coef = i,
                           maxFDRForParticipation = 0.05,
                           useSubset = FALSE)
    
    print(paste0('splicing report ', i))
    # Bin and junction signal integration:
    sr <- splicingReport(gbDUreport, jDUreport, gbcounts) 
    
    print(paste0('integrate signals ', i))
    # Summary of integration of splicing signals along genomic regions:
    is <- integrateSignals(sr, 
                           asd, 
                           bin.FC = 2, 
                           bin.fdr = 0.05, 
                           bin.inclussion = 0.1,
                           nonunif = 1, 
                           usenonunif = FALSE,
                           bjs.inclussion = 0.1, 
                           bjs.fdr = 0.05,
                           a.inclussion = 0.1, 
                           a.fdr = 0.05,
                           l.inclussion = 0.1, 
                           l.fdr = 0.05,
                           otherSources = NULL, 
                           overlapType = "any")
    
    
    # Region-based information of splicing signals:
    signals <- signals(is)
    
    print(paste0('export Integrated Signals ', i))
    # Export region-based integrated signal results into interactive HTML pages:
    exportIntegratedSignals(is, 
                            output.dir= paste0("is_", sra, "_", i), 
                            sr, 
                            gbcounts, 
                            features, 
                            asd,
                            useLog = FALSE, 
                            ntop = NULL, 
                            mergedBams = mBAMs,
                            makeGraphs=TRUE)
    
    print(paste0('saving splicing report ', i))
    # Export region-based integrated signal results in tab delimited table:
    writeSplicingReport(sr, output.dir = paste0("sr_", sra, "_", i))
    # Make list with results for this coef:
    name <- paste("coef:", i, sep='_')
    results[[name]] <- list(gbDUreport, jDUreport, sr, is, signals)
  }
  
  print(paste0(sra, ' analysis ended'))
  
  return(list(asd = asd,
              gbcounts = gbcounts,
              results = results))
}
