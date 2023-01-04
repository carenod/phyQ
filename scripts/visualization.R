#limpio la memoria
rm( list=ls() )  #remove all objects
gc()             #garbage collection
Sys.setenv(http_proxy="172.16.254.254:3128")
Sys.setenv(https_proxy="172.16.254.254:3128")
library(gviz)

# helpful tutorials
# https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/
# https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/presentations/singlepage/Viz_part_1.html#TxDb_to_GeneRegionTrack
# https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/presentations/singlepage/Viz_part_2.html#Recap2
# https://rockefelleruniversity.github.io/RU_VisualizingGenomicsData/presentations/singlepage/Viz_part_3.html
