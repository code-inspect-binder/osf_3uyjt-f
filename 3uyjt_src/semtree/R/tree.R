# script to grow a SEM tree with the quadratic LGCM
# probably takes half an hour on a 2023-ish laptop

# load semtree library
library(semtree)

# make sure to run preprocessing first (path is relative to project root)
#source("R/preprocessing_hrs.R")
rndhrs_subset <- read.csv("data/rndhrs_subset.csv")

df <- rndhrs_subset

# filter only those needed for the model
modelData <- rndhrs_subset %>%  select(
  R4AGEY_B:R11AGEY_B, R4TR20, R5TR20, R6TR20, R7TR20, 
  R8TR20, R9TR20, R10TR20,R11TR20,CAGE_T4:CAGE_T11, CAGE2_T4:CAGE2_T11)

# OpenMx quadraticLGCM with definition variables
manifests<-c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20")
latents<-c("icept","slope","quad")
model_quadratic <- mxModel("LGCMl", 
                           type="RAM",
                           manifestVars = manifests,
                           latentVars = latents,
                           mxPath(from="icept",to=c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20"),
                                  free=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE), 
                                  value=c(1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0) , arrows=1, 
                                  label=c("icept__R4TR20","icept__R5TR20","icept__R6TR20","icept__R7TR20","icept__R8TR20","icept__R9TR20","icept__R10TR20","icept__R11TR20") ),
                           mxPath(from="slope",to=c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20"), 
                                  free=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),  arrows=1, 
                                  label=c("data.CAGE_T4","data.CAGE_T5","data.CAGE_T6","data.CAGE_T7","data.CAGE_T8","data.CAGE_T9","data.CAGE_T10","data.CAGE_T11") ),
                           mxPath(from="quad",to=c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20"), 
                                  free=c(FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),  arrows=1, 
                                  label=c("data.CAGE2_T4","data.CAGE2_T5","data.CAGE2_T6","data.CAGE2_T7","data.CAGE2_T8","data.CAGE2_T9","data.CAGE2_T10","data.CAGE2_T11") ),                 
                           mxPath(from="one",to=c("icept","slope","quad"), free=TRUE, value=c(1.0,1.0,1.0) , arrows=1, 
                                  label=c("icept_mu","slope_mu","quad_mu") ),
                           mxPath(from="icept",to=c("slope","icept","quad"), free=TRUE, value=c(0.0,1.0,0.0) , arrows=2, 
                                  label=c("icept_slope_cov","icept_var","icept_quad_cov") ),
                           mxPath(from="slope",to=c("slope","quad"), free=c(TRUE), value=c(1.0,0.0) ,
                                  arrows=2, label=c("slope_var","slope_quad_cov") ),
                           mxPath(from="quad",to=c("quad"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("quad_var") ),
                           mxPath(from="R4TR20",to=c("R4TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R5TR20",to=c("R5TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R6TR20",to=c("R6TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R7TR20",to=c("R7TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R8TR20",to=c("R8TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R9TR20",to=c("R9TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R10TR20",to=c("R10TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R11TR20",to=c("R11TR20"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("residual_var") ),
                           mxPath(from="one",to=c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20"), free=F, value=0, arrows=1),
                           mxData(as.data.frame(modelData), type = "raw")
                           
)

# fit model and display summary
fitted_quadratic_lgcm <- mxRun(model_quadratic)
summary(fitted_quadratic_lgcm)


# fit tree with quadratic model
hrstree <- semtree(fitted_quadratic_lgcm, 
                   as.data.frame(rndhrs_subset),
                   control=semtree.control(verbose=TRUE,
                                           method="naive",
                                           missing="party",
                                           min.bucket = 500,
                                           min.N = 250,
                                           exclude.heywood = FALSE
                                           #                                           heywood
                   ),
                   covariates=vars_AB 
)

# save results
save(hrstree, file="results/hrs-lgcm-tree-quad.sav")

# plot first two levels of the tree
plot(prune(hrstree, 2))
