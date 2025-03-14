## ---- include = FALSE---------------------------------------------------------
dpi = 125

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  dpi=dpi, 
  fig.retina=1, 
  fig.width=1440/dpi, #4:3 FHD
  fig.height=1080/dpi, 
  out.width="100%",
  crop = NULL,
  warning = T, 
  error = T
)

rm(dpi)

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("Coxmos")

## ---- eval = FALSE------------------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("BiostatOmics/Coxmos", build_vignettes = TRUE)

## ----setup, results = "hide"--------------------------------------------------
# load Coxmos
library(Coxmos)

## ---- eval=FALSE--------------------------------------------------------------
#  # install.packages("RColorConesa")
#  library(RColorConesa)

## -----------------------------------------------------------------------------
# load dataset
data("X_multiomic", package = "Coxmos")
data("Y_multiomic", package = "Coxmos")

X <- X_multiomic
Y <- Y_multiomic

rm(X_multiomic, Y_multiomic)

## ---- echo = FALSE------------------------------------------------------------
knitr::kable(X$mirna[1:5,1:5]);knitr::kable(X$proteomic[1:5,1:5])

knitr::kable(Y[1:5,])

## -----------------------------------------------------------------------------
ggp_density.event <- plot_events(Y = Y, 
                                 categories = c("Censored","Death"), #name for FALSE/0 (Censored) and TRUE/1 (Event)
                                 y.text = "Number of observations", 
                                 roundTo = 0.5, 
                                 max.breaks = 15) 

## ----fig.small = T------------------------------------------------------------
ggp_density.event$plot

## -----------------------------------------------------------------------------
set.seed(123)
index_train <- caret::createDataPartition(Y$event,
                                          p = .7, #70 %
                                          list = FALSE,
                                          times = 1)

X_train <- list()
X_test <- list()
for(omic in names(X)){
  X_train[[omic]] <- X[[omic]][index_train,,drop=F]
  X_test[[omic]] <- X[[omic]][-index_train,,drop=F]
}

Y_train <- Y[index_train,]
Y_test <- Y[-index_train,]

## -----------------------------------------------------------------------------
EPV <- getEPV.mb(X_train, Y_train)
for(b in names(X_train)){
  message(paste0("EPV = ", round(EPV[[b]], 4), ", for block ", b))
}


## ---- message=F---------------------------------------------------------------
x.center = c(mirna = T, proteomic = T) #if vector, must be named
x.scale = c(mirna = F, proteomic = F) #if vector, must be named

## ----warning=T, eval=F--------------------------------------------------------
#  cv.sb.splsicox_res <- cv.sb.splsicox(X = X_train, Y = Y_train,
#                                       max.ncomp = 2, penalty.list = c(0.5,0.9),
#                                       n_run = 2, k_folds = 5,
#                                       x.center = x.center, x.scale = x.scale,
#                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                       remove_variance_at_fold_level = F,
#                                       remove_non_significant_models = F, alpha = 0.05,
#                                       w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                       PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.sb.splsicox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.sb.splsicox_res$plot_AUC

## -----------------------------------------------------------------------------
sb.splsicox_model <- sb.splsicox(X = X_train, Y = Y_train,
                                 n.comp = 1, #cv.sb.splsicox_res$opt.comp, 
                                 penalty = 0.9, #cv.sb.splsicox_res$opt.penalty,
                                 x.center = x.center, x.scale = x.scale,
                                 remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                 remove_non_significant = F, 
                                 alpha = 0.05, MIN_EPV = 5, 
                                 returnData = T, verbose = F)

sb.splsicox_model

## -----------------------------------------------------------------------------
sb.splsicox_model <- sb.splsicox(X = X_train, Y = Y_train,
                                 n.comp = 1, #cv.sb.splsicox_res$opt.comp,
                                 penalty = 0.9, #cv.sb.splsicox_res$opt.penalty,
                                 x.center = x.center, x.scale = x.scale,
                                 remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
                                 remove_non_significant = T,
                                 alpha = 0.05, MIN_EPV = 5, 
                                 returnData = T, verbose = F)

sb.splsicox_model

## ----warning=T, eval=F--------------------------------------------------------
#  cv.isb.splsicox_res <- cv.isb.splsicox(X = X_train, Y = Y_train,
#                                       max.ncomp = 2, penalty.list = c(0.5,0.9),
#                                       n_run = 2, k_folds = 5,
#                                       x.center = x.center, x.scale = x.scale,
#                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                       remove_variance_at_fold_level = F,
#                                       remove_non_significant_models = F, alpha = 0.05,
#                                       w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                       PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.isb.splsicox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.isb.splsicox_res$list_cv_spls_models$mirna$plot_AUC
#  cv.isb.splsicox_res$list_cv_spls_models$proteomic$plot_AUC

## ---- warning=F, eval=F-------------------------------------------------------
#  isb.splsicox_model <- isb.splsicox(X = X_train, Y = Y_train, cv.isb = cv.isb.splsicox_res,
#                                        x.center = x.center, x.scale = x.scale,
#                                        remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
#                                     remove_non_significant = FALSE, alpha = 0.05,
#                                     MIN_EPV = 5, returnData = TRUE, verbose = FALSE)
#  
#  isb.splsicox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.sb.splsdrcox_res <- cv.sb.splsdrcox(X = X_train, Y = Y_train,
#                                         max.ncomp = 2, vector = NULL,
#                                         n_run = 2, k_folds = 10,
#                                         x.center = x.center, x.scale = x.scale,
#                                         #y.center = FALSE, y.scale = FALSE,
#                                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                         remove_variance_at_fold_level = F,
#                                         remove_non_significant_models = F, alpha = 0.05,
#                                         w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                         MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                         pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                         MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                         PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.sb.splsdrcox_res

## -----------------------------------------------------------------------------
sb.splsdrcox_model <- sb.splsdrcox(X = X_train, 
                                   Y = Y_train, 
                                   n.comp = 2, #cv.sb.splsdrcox_res$opt.comp, 
                                   vector = list("mirna" = 484, "proteomic" = 369), #cv.sb.splsdrcox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale,
                                   remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL, 
                                   remove_non_significant = T, alpha = 0.05, MIN_EPV = 5,
                                   returnData = T, verbose = F)

sb.splsdrcox_model

## ----warning=T, eval=F--------------------------------------------------------
#  cv.isb.splsdrcox_res <- cv.isb.splsdrcox(X = X_train, Y = Y_train,
#                                       max.ncomp = 2, vector = NULL,
#                                       n_run = 2, k_folds = 5,
#                                       x.center = x.center, x.scale = x.scale,
#                                       remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                       remove_variance_at_fold_level = F,
#                                       remove_non_significant_models = F, alpha = 0.05,
#                                       w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                       MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                       pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                       MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                       PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.isb.splsdrcox_res

## ---- fig.small=T, eval=F-----------------------------------------------------
#  cv.isb.splsdrcox_res$list_cv_spls_models$mirna$plot_AUC
#  cv.isb.splsdrcox_res$list_cv_spls_models$proteomic$plot_AUC

## ---- warning=F, eval=F-------------------------------------------------------
#  isb.splsdrcox_model <- isb.splsdrcox(X = X_train, Y = Y_train, cv.isb = cv.isb.splsdrcox_res,
#                                     x.center = x.center, x.scale = x.scale,
#                                     remove_near_zero_variance = TRUE, remove_zero_variance = TRUE, toKeep.zv = NULL,
#                                     remove_non_significant = FALSE, alpha = 0.05,
#                                     MIN_EPV = 5, returnData = TRUE, verbose = FALSE)
#  
#  isb.splsdrcox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  cv.mb.splsdrcox_res <- cv.mb.splsdrcox(X = X_train, Y = Y_train,
#                                         max.ncomp = 2, vector = NULL, #NULL - autodetection
#                                         MIN_NVAR = 10, MAX_NVAR = NULL, n.cut_points = 10, EVAL_METHOD = "AUC",
#                                         n_run = 2, k_folds = 4,
#                                         x.center = x.center, x.scale = x.scale,
#                                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                         remove_variance_at_fold_level = F,
#                                         remove_non_significant_models = F, alpha = 0.05,
#                                         w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                         MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                         pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                         MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                         PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.mb.splsdrcox_res

## -----------------------------------------------------------------------------
mb.splsdrcox_model <- mb.splsdrcox(X = X_train, Y = Y_train, 
                                   n.comp = 2, #cv.mb.splsdrcox_res$opt.comp,
                                   vector = list("mirna" = 326, "proteomic" = 369), #cv.mb.splsdrcox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale, 
                                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                                   remove_non_significant = T, alpha = 0.05,
                                   MIN_AUC_INCREASE = 0.01,
                                   pred.method = "cenROC", max.iter = 200,
                                   times = NULL, max_time_points = 15,
                                   MIN_EPV = 5, returnData = T, verbose = F)

mb.splsdrcox_model

## ---- warning=F, eval=F-------------------------------------------------------
#  # run cv.splsdrcox
#  cv.mb.splsdacox_res <- cv.mb.splsdacox(X = X_train, Y = Y_train,
#                                         max.ncomp = 2, vector = NULL, #NULL - autodetection
#                                         n_run = 2, k_folds = 4,
#                                         x.center = x.center, x.scale = x.scale,
#                                         remove_near_zero_variance = T, remove_zero_variance = F, toKeep.zv = NULL,
#                                         remove_variance_at_fold_level = F,
#                                         remove_non_significant_models = F, alpha = 0.05,
#                                         w_AIC = 0, w_C.Index = 0, w_AUC = 1, w_I.BRIER = 0, times = NULL, max_time_points = 15,
#                                         MIN_AUC_INCREASE = 0.01, MIN_AUC = 0.8, MIN_COMP_TO_CHECK = 3,
#                                         pred.attr = "mean", pred.method = "cenROC", fast_mode = F,
#                                         MIN_EPV = 5, return_models = F, remove_non_significant = F, returnData = F,
#                                         PARALLEL = F, verbose = F, seed = 123)
#  
#  cv.mb.splsdacox_res

## -----------------------------------------------------------------------------
mb.splsdacox_model <- mb.splsdacox(X = X_train, Y = Y_train, 
                                   n.comp = 2, #cv.mb.splsdacox_res$opt.comp,
                                   vector = list("mirna" = 326, "proteomic" = 10), #cv.mb.splsdacox_res$opt.nvar,
                                   x.center = x.center, x.scale = x.scale, 
                                   remove_near_zero_variance = T, remove_zero_variance = T, toKeep.zv = NULL,
                                   remove_non_significant = T, alpha = 0.05,
                                   MIN_AUC_INCREASE = 0.01,
                                   pred.method = "cenROC", max.iter = 200,
                                   times = NULL, max_time_points = 15,
                                   MIN_EPV = 5, returnData = T, verbose = F)

mb.splsdacox_model

## -----------------------------------------------------------------------------
lst_models <- list("SB.sPLS-ICOX" = sb.splsicox_model,
                   #"iSB.sPLS-ICOX" = isb.splsicox_model,
                   "SB.sPLS-DRCOX-Dynamic" = sb.splsdrcox_model,
                   #"iSB.sPLS-DRCOX-Dynamic" = isb.splsdrcox_model,
                   #"SB.sPLS-DRCOX-Penalty" = sb.splsdrcox_penalty_model,
                   #"iSB.sPLS-DRCOX-Penalty" = isb.splsdrcox_penalty_model,
                   #"SB.sPLS-DACOX-Dynamic" = sb.splsdacox_model,
                   #"iSB.sPLS-DACOX-Dynamic" = isb.splsdacox_model,
                   "MB.sPLS-DRCOX" = mb.splsdrcox_model,
                   "MB.sPLS-DACOX" = mb.splsdacox_model)

eval_results <- eval_Coxmos_models(lst_models = lst_models,
                                  X_test = X_test, Y_test = Y_test, 
                                  pred.method = "cenROC",
                                  pred.attr = "mean",
                                  times = NULL, max_time_points = 15, 
                                  PARALLEL = F)

## ---- eval=FALSE--------------------------------------------------------------
#  lst_evaluators <- c(cenROC = "cenROC", risksetROC = "risksetROC")
#  
#  eval_results <- purrr::map(lst_evaluators, ~eval_Coxmos_models(lst_models = lst_models,
#                                                                X_test = X_test, Y_test = Y_test,
#                                                                pred.method = .,
#                                                                pred.attr = "mean",
#                                                                times = NULL,
#                                                                max_time_points = 15,
#                                                                PARALLEL = F))

## -----------------------------------------------------------------------------
eval_results

## ---- warning=F---------------------------------------------------------------
lst_eval_results <- plot_evaluation(eval_results, evaluation = "AUC")
lst_eval_results_brier <- plot_evaluation(eval_results, evaluation = "IBS")

## ---- fig.small=T, warning=F--------------------------------------------------
lst_eval_results$lst_plots$lineplot.mean
lst_eval_results$lst_plot_comparisons$t.test

# lst_eval_results$cenROC$lst_plots$lineplot.mean
# lst_eval_results$cenROC$lst_plot_comparisons$t.test

## -----------------------------------------------------------------------------
lst_models_time <- list(#cv.sb.splsicox_res,
                        sb.splsicox_model,
                        #isb.splsicox_model,
                        #cv.sb.splsdrcox_res,
                        sb.splsdrcox_model,
                        #isb.splsdrcox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdrcox_model,
                        #cv.mb.splsdrcox_res,
                        mb.splsdacox_model,
                        eval_results)

## -----------------------------------------------------------------------------
ggp_time <- plot_time.list(lst_models_time, txt.x.angle = 90)

## ---- fig.small=T-------------------------------------------------------------
ggp_time

## -----------------------------------------------------------------------------
#lst_forest_plot <- plot_forest.list(lst_models)
lst_forest_plot <- plot_forest(lst_models$`SB.sPLS-DRCOX`)

## ---- fig.small=T-------------------------------------------------------------
#lst_forest_plot$`SB.sPLS-DRCOX`
lst_forest_plot

## -----------------------------------------------------------------------------
#lst_ph_ggplot <- plot_proportionalHazard.list(lst_models)
lst_ph_ggplot <- plot_proportionalHazard(lst_models$`SB.sPLS-DRCOX`)

## ---- fig.small=T-------------------------------------------------------------
#lst_ph_ggplot$`SB.sPLS-DRCOX`
lst_ph_ggplot

## -----------------------------------------------------------------------------
#density.plots.lp <- plot_cox.event.list(lst_models, type = "lp")
density.plots.lp <- plot_cox.event(lst_models$`SB.sPLS-DRCOX`, type = "lp")

## ---- fig.small=T-------------------------------------------------------------
density.plots.lp$plot.density
density.plots.lp$plot.histogram

## -----------------------------------------------------------------------------
ggp_scores <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`,
                             comp = c(1,2), mode = "scores")

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_scores$plot_block

## -----------------------------------------------------------------------------
ggp_loadings <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`, 
                               comp = c(1,2), mode = "loadings",
                               top = 10)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_loadings$plot_block

## -----------------------------------------------------------------------------
ggp_biplot <- plot_PLS_Coxmos(model = lst_models$`SB.sPLS-DRCOX`, 
                             comp = c(1,2), mode = "biplot",
                             top = 15,
                             only_top = T)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp_biplot$plot_block

## ---- warning=F---------------------------------------------------------------
variable_auc_results <- eval_Coxmos_model_per_variable(model = lst_models$`SB.sPLS-DRCOX`, 
                                                      X_test = lst_models$`SB.sPLS-DRCOX`$X_input, 
                                                      Y_test = lst_models$`SB.sPLS-DRCOX`$Y_input,
                                                      pred.method = "cenROC", pred.attr = "mean",
                                                      times = NULL, max_time_points = 15,
                                                      PARALLEL = FALSE)

variable_auc_plot_train <- plot_evaluation(variable_auc_results, evaluation = "AUC")

## ---- fig.small=T, warning=F--------------------------------------------------
variable_auc_plot_train$lst_plots$lineplot.mean

## ---- warning=FALSE-----------------------------------------------------------
# ggp.simulated_beta <- plot_pseudobpenalty.list(lst_models = lst_models, 
#                                            error.bar = T, onlySig = T, alpha = 0.05, 
#                                            zero.rm = T, auto.limits = T, top = 20,
#                                            show_percentage = T, size_percentage = 2, verbose = F)

ggp.simulated_beta <- plot_pseudobeta(model = lst_models$`SB.sPLS-DRCOX`, 
                                      error.bar = T, onlySig = T, alpha = 0.05, 
                                      zero.rm = T, auto.limits = T, top = 20,
                                      show_percentage = T, size_percentage = 2)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp.simulated_beta$plot

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp.simulated_beta$mb_plot$plot

## ---- warning=F---------------------------------------------------------------
# LST_KM_RES_LP <- getAutoKM.list(type = "LP",
#                                 lst_models = lst_models,
#                                 comp = 1:4,
#                                 top = 10,
#                                 ori_data = T,
#                                 BREAKTIME = NULL,
#                                 only_sig = T, alpha = 0.05)

LST_KM_RES_LP <- getAutoKM(type = "LP",
                           model = lst_models$`SB.sPLS-DRCOX`,
                           comp = 1:4,
                           top = 10,
                           ori_data = T,
                           BREAKTIME = NULL,
                           only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_LP$LST_PLOTS$LP

## ---- warning=FALSE-----------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_LP)
# LST_KM_TEST_LP <- getTestKM.list(lst_models = lst_models, 
#                                  X_test = X_test, Y_test = Y_test, 
#                                  type = "LP",
#                                  BREAKTIME = NULL, n.breaks = 20,
#                                  lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_LP)
LST_KM_TEST_LP <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                            X_test = X_test, Y_test = Y_test, 
                            type = "LP",
                            BREAKTIME = NULL, n.breaks = 20,
                            cutoff = lst_cutoff)

## ---- warning=FALSE-----------------------------------------------------------
LST_KM_TEST_LP

## ---- warning=F---------------------------------------------------------------
# LST_KM_RES_COMP <- getAutoKM.list(type = "COMP",
#                                   lst_models = lst_models,
#                                   comp = 1:4,
#                                   top = 10,
#                                   ori_data = T,
#                                   BREAKTIME = NULL,
#                                   only_sig = T, alpha = 0.05)

LST_KM_RES_COMP <- getAutoKM(type = "COMP",
                             model = lst_models$`SB.sPLS-DRCOX`,
                             comp = 1:4,
                             top = 10,
                             ori_data = T,
                             BREAKTIME = NULL,
                             only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_COMP$LST_PLOTS$mirna$comp_2
LST_KM_RES_COMP$LST_PLOTS$proteomic$comp_1
LST_KM_RES_COMP$LST_PLOTS$proteomic$comp_2

## ---- warning=F---------------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_COMP)
# LST_KM_TEST_COMP <- getTestKM.list(lst_models = lst_models, 
#                                    X_test = X_test, Y_test = Y_test, 
#                                    type = "COMP",
#                                    BREAKTIME = NULL, n.breaks = 20,
#                                    lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_COMP)
LST_KM_TEST_COMP <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                              X_test = X_test, Y_test = Y_test, 
                              type = "COMP",
                              BREAKTIME = NULL, n.breaks = 20,
                              cutoff = lst_cutoff)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_TEST_COMP$comp_2_mirna
LST_KM_TEST_COMP$comp_1_proteomic
LST_KM_TEST_COMP$comp_2_proteomic

## ---- warning=F---------------------------------------------------------------
# LST_KM_RES_VAR <- getAutoKM.list(type = "VAR",
#                                  lst_models = lst_models,
#                                  comp = 1:4,
#                                  top = 10,
#                                  ori_data = T,
#                                  BREAKTIME = NULL,
#                                  only_sig = T, alpha = 0.05)

LST_KM_RES_VAR <- getAutoKM(type = "VAR",
                            model = lst_models$`SB.sPLS-DRCOX`,
                            comp = 1:4,
                            top = 10,
                            ori_data = T,
                            BREAKTIME = NULL,
                            only_sig = T, alpha = 0.05)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_RES_VAR$LST_PLOTS$mirna$hsa.minus.miR.minus.21.minus.5p
LST_KM_RES_VAR$LST_PLOTS$proteomic$var_840
LST_KM_RES_VAR$LST_PLOTS$proteomic$var_7535

## ---- warning=FALSE-----------------------------------------------------------
# lst_cutoff <- getCutoffAutoKM.list(LST_KM_RES_VAR)
# LST_KM_TEST_VAR <- getTestKM.list(lst_models = lst_models, 
#                                   X_test = X_test, Y_test = Y_test, 
#                                   type = "VAR", ori_data = T,
#                                   BREAKTIME = NULL, n.breaks = 20,
#                                   lst_cutoff = lst_cutoff)

lst_cutoff <- getCutoffAutoKM(LST_KM_RES_VAR)
LST_KM_TEST_VAR <- getTestKM(model = lst_models$`SB.sPLS-DRCOX`, 
                             X_test = X_test, Y_test = Y_test, 
                             type = "VAR", ori_data = T,
                             BREAKTIME = NULL, n.breaks = 20,
                             cutoff = lst_cutoff)

## ---- fig.small=T, warning=FALSE----------------------------------------------
LST_KM_TEST_VAR$mirna$hsa.minus.miR.minus.21.minus.5p
LST_KM_TEST_VAR$proteomic$var_840
LST_KM_TEST_VAR$proteomic$var_7535

## -----------------------------------------------------------------------------
new_pat <- list()
for(b in names(X_test)){
  new_pat[[b]] <- X_test[[b]][1,,drop=F]
}


## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(new_pat$mirna),])

## ---- warning=FALSE-----------------------------------------------------------
# ggp.simulated_beta_newPat <- plot_observation.pseudobeta.list(lst_models = lst_models, 
#                                                              observation = new_pat,
#                                                              error.bar = T, onlySig = T, alpha = 0.05,
#                                                              zero.rm = T, auto.limits = T, show.betas = T, top = 20)

ggp.simulated_beta_newPat <- plot_observation.pseudobeta(model = lst_models$`SB.sPLS-DRCOX`,
                                                        observation = new_pat,
                                                        error.bar = T, onlySig = T, alpha = 0.05,
                                                        zero.rm = T, auto.limits = T, show.betas = T, 
                                                        top = 20, txt.x.angle = 90)

## ---- fig.small=T, warning=FALSE----------------------------------------------
ggp.simulated_beta_newPat$plot$mirna
ggp.simulated_beta_newPat$plot$proteomic

## ---- warning=FALSE-----------------------------------------------------------
pat_density <- plot_observation.eventDensity(observation = new_pat, 
                                             model = lst_models$`SB.sPLS-DRCOX`, 
                                             time = NULL, 
                                             type = "lp")

## ---- fig.small=T, warning=FALSE----------------------------------------------
pat_density

## -----------------------------------------------------------------------------
pat_histogram <- plot_observation.eventHistogram(observation = new_pat, 
                                                 model = lst_models$`SB.sPLS-DRCOX`, 
                                                 time = NULL, 
                                                 type = "lp")

## ---- fig.small=T, warning=FALSE----------------------------------------------
pat_histogram

## -----------------------------------------------------------------------------
sub_X_test <- list()
for(b in names(X_test)){
  sub_X_test[[b]] <- X_test[[b]][1:5,]
}

## -----------------------------------------------------------------------------
knitr::kable(Y_test[rownames(sub_X_test$proteomic),])

## ---- warning=FALSE-----------------------------------------------------------
# lst_cox.comparison <- plot_multipleObservations.LP.list(lst_models = lst_models,
#                                                     observations = sub_X_test,
#                                                     error.bar = T, zero.rm = T, onlySig = T,
#                                                     alpha = 0.05, top = 5)

lst_cox.comparison <- plot_multipleObservations.LP(model = lst_models$`SB.sPLS-DRCOX`,
                                                   observations = sub_X_test,
                                                   error.bar = T, zero.rm = T, onlySig = T,
                                                   alpha = 0.05, top = 5)

## ---- fig.small=T, warning=FALSE----------------------------------------------
lst_cox.comparison$plot

