################################################################################
################################  SEMTree ######################################
################################################################################
library(devtools)
devtools::install_github("brandmaier/semtree",force=TRUE, build_opts = c()) #aktuellstes Package von github runterladen!


source("/dss/dsshome1/0B/ra35tik2/paper2/lvforest/lvforest_semtree.R")
start_time <- Sys.time()
univ_forest <- lvforest(
  data = simu,
  split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5',
  pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5
                Eta1 ~ 1
                simuvar1|0*t1',
  input = c("num1","ord1","cat1"),  
  idvar = "id",
  ntrees= 10,
  bagging=NULL,
  cutoff_rmsea = .05,
  ordered = TRUE,
  minn = 200,
  split=2,   #sqrt(length(input))
  estmethod = "EBM",
  stab_bonferroni=FALSE
)
end_time <- Sys.time()
runtime_forest <- end_time - start_time

#write result in simu
simu$forest_scores = univ_forest[["Scores"]][["Eta1"]][match(simu$id,univ_forest[["Scores"]][["Eta1"]]$id),"lvscore"]
simu = simu[,c(ncol(simu),1:(ncol(simu)-1))]

#write results
save(simu, univ_forest, runtime_forest,file="/dss/dsshome1/0B/ra35tik2/paper2/simulation/unidimdensional/lvforest_simulation/240126_single_dataset_result.RData")
