################################################################################
################################  SEMTree ######################################
################################################################################
source("/dss/dsshome1/0B/ra35tik2/paper2/lvforest/lvforest_semtree.R")
start_time <- Sys.time()
univ_forest <- lvforest(
  data = simu,
  split_model ='Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8',
  pred_model = 'Eta1 =~ simuvar1 + simuvar2 + simuvar3 + simuvar4 + simuvar5 + simuvar6 + simuvar7 + simuvar8
                Eta1 ~ 1
                simuvar1|0*t1',
  input = input_sampled,  
  idvar = "id",
  ntrees= 20,
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
simu$forest_scores_part = univ_forest[["Scores"]][["Eta1"]][match(simu$id,univ_forest[["Scores"]][["Eta1"]]$id),"lvscore"]
simu = simu[,c(ncol(simu),1:(ncol(simu)-1))]