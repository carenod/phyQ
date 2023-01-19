##### ASpli pipeline ######
# ASpli_pipeline function to perform comparisons

require(ASpli)
require(GenomicFeatures)

ASpli_contrast <- function(contrast = NULL,
                           coef = NULL,
                           formula = NULL,
                           sra = NULL, 
                           features = NULL, 
                           mBAMs = NULL, 
                           gbcounts = NULL, 
                           asd = NULL) {

  print('Doing gbDUreport')
  
  # Differential gene expression and bin usage signal estimation
  gbDUreport <- gbDUreport(counts = gbcounts, 
                           contrast = contrast, 
                           coef = coef,
                           formula = formula,
                           minBinReads = 5)
  
  print('Doing jDUreport')
  # Differential junction usage analysis
  jDUreport <- jDUreport(asd = asd, 
                         minAvgCounts = 5,
                         filterWithContrasted = TRUE,
                         runUniformityTest = TRUE,
                         mergedBams = mBAMs,
                         maxPValForUniformityCheck = 0.2, 
                         strongFilter = TRUE,
                         maxConditionsForDispersionEstimate = 24,
                         contrast = contrast,
                         coef = coef,
                         formula = formula,
                         maxFDRForParticipation = 0.05,
                         useSubset = FALSE)
  
  print('splicing report')
  # Bin and junction signal integration:
  sr <- splicingReport(gbDUreport, jDUreport, gbcounts) 
  
  print('integrate signals')
  # Summary of integration of splicing signals along genomic regions:
  is <- integrateSignals(sr, 
                         asd, 
                         bin.FC = 0.58, 
                         bin.fdr = 0.1, 
                         bin.inclussion = 0.05,
                         nonunif = 1, 
                         usenonunif = FALSE,
                         bjs.inclussion = 0.05, 
                         bjs.fdr = 0.1,
                         a.inclussion = 0.05, 
                         a.fdr = 0.1,
                         l.inclussion = 0.05, 
                         l.fdr = 0.1,
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

