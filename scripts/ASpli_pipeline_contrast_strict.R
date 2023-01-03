##### ASpli pipeline ######
# ASpli_pipeline function to execute gbcounts, asd, and return 
# It has a strict filter for minBinReads
ASpli_pipeline_contrast_strict <- function(targets, minReadLength, libType, contrast, 
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
                 threshold = 5,                 
                 minAnchor = 10,
                 strandMode = strandMode)
  
  print('Doing gbDUreport')
  
  # Differential gene expression and bin usage signal estimation
  gbDUreport <- gbDUreport(gbcounts, 
                           contrast = contrast, 
                           formula = form, 
                           coef = coef,
                           minBinReads = 50)
  
  print('Doing jDUreport')
  # Differential junction usage analysis
  jDUreport <- jDUreport(asd, minAvgCounts = 50,
                         filterWithContrasted = TRUE,
                         runUniformityTest = TRUE,
                         mergedBams = mBAMs,
                         maxPValForUniformityCheck = 0.2, 
                         strongFilter = TRUE,
                         maxConditionsForDispersionEstimate = 24,
                         contrast = contrast,
                         formula = form,
                         coef = coef,
                         maxFDRForParticipation = 0.05,
                         useSubset = FALSE)

  print('splicing report')
  # Bin and junction signal integration:
  sr <- splicingReport(gbDUreport, jDUreport, gbcounts) 
  
  print('integrate signals')
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
  
  print('export Integrated Signals')
  exportIntegratedSignals(is, 
                          output.dir= paste0("is_", sra), 
                          sr, 
                          gbcounts, 
                          features, 
                          asd,
                          useLog = FALSE, 
                          ntop = NULL, 
                          mergedBams = mBAMs,
                          makeGraphs=TRUE)
  
  print('saving splicing report')
  writeSplicingReport(sr, output.dir = paste0("sr_", sra))


    print(paste0(sra, ' analysis ended'))
    
    return(list(gbcounts = gbcounts, 
                asd = asd, 
                gbDUreport = gbDUreport, 
                jDUreport = jDUreport, 
                sr = sr, 
                is = is,  
                signals = signals(is)))
}

