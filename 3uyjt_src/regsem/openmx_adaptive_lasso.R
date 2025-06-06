library(regsem);library(OpenMx)
# run preprocessing_hrs.R first

#df <- rndhrs_subset

# filter only those needed for the model
modelData <- rndhrs_subset %>%  select(
  R4AGEY_B:R11AGEY_B, R4TR20, R5TR20, R6TR20, R7TR20, 
  R8TR20, R9TR20, R10TR20,R11TR20,CAGE_T4:CAGE_T11, CAGE2_T4:CAGE2_T11,
  all_of(vars_AB))

modelData <- as.data.frame(modelData)

modelData$RAEDYRS = modelData$RAEDYRS - 12
modelData$RAGENDER = modelData$RAGENDER - 1
#modelData$R2WEIGHT = scale(modelData$R2WEIGHT)
#modelData$R2BMI = scale(modelData$R2BMI)

#sub_dat = modelData[,66:ncol(modelData)]
#recode_id1 = sub_dat == 9 & is.na(sub_dat)==F
#sub_dat[recode_id1]  <- NA
#recode_id2 = sub_dat > 1 & is.na(sub_dat)==F
#sub_dat[recode_id2]  <- 1
modelData$R2PUSH[modelData$R2PUSH > 1 & is.na(modelData$R2PUSH)==F] <- 1

#outt = readRDS("modelData.rds")
#modelData = outt$modelData
#cov_est = outt$cov_est

cov.mat = cov(modelData[,vars_AB],use="pairwise.complete.obs")
cor.mat = cor(modelData[,vars_AB],use="pairwise.complete.obs")
cov_est = cov.mat[lower.tri(cov.mat,diag=T)]

vars_AB22 <-c("RAEDYRS","RAGENDER","R2SHLT","R2HOSP",  "RAHISPAN","race_black","race_other",
              "R2SMOKEV","R3CHOLST","R3FLUSHT","R2DRINK","R2CHAIRA","R2BACK","R2JOGA")

# include covariates

indicators<-c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20")
covariates = c(vars_AB22)
latents<-c("icept","slope","quad")
baseModel <- mxModel("baseModel", 
                           type="RAM",
                           manifestVars = c(indicators,covariates),
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
                           mxPath(from="one",to=c("icept","slope","quad"), free=TRUE, value=c(10.5,-.13,-.31) , arrows=1, 
                                  label=c("icept_mu","slope_mu","quad_mu") ),
                           mxPath(from="icept",to=c("slope","icept","quad"), free=TRUE, value=c(-.026,5.6,-.62) , arrows=2, 
                                  label=c("icept_slope_cov","icept_var","icept_quad_cov") ),
                           mxPath(from="slope",to=c("slope","quad"), free=c(TRUE), value=c(.01,.01) ,
                                  arrows=2, label=c("slope_var","slope_quad_cov") ),
                           mxPath(from="quad",to=c("quad"), free=c(TRUE), value=c(.043) , arrows=2, label=c("quad_var") ),
                           mxPath(from="R4TR20",to=c("R4TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R5TR20",to=c("R5TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R6TR20",to=c("R6TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R7TR20",to=c("R7TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R8TR20",to=c("R8TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R9TR20",to=c("R9TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R10TR20",to=c("R10TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from="R11TR20",to=c("R11TR20"), free=c(TRUE), value=c(5.1) , arrows=2, label=c("residual_var") ),
                           mxPath(from=covariates,to="slope",arrows=1,free=T,values=0,label=paste("reg2",1:length(covariates),sep="")),
                           mxPath(from=covariates,to="icept",arrows=1,free=T,values=0,label=paste("reg1",1:length(covariates),sep="")),
                           mxPath(from=covariates,to="quad",arrows=1,free=T,values=0,label=paste("reg3",1:length(covariates),sep="")),
                           mxPath(from=covariates,arrows=2,connect="unique.pairs",free=F,values=cov_est),
                           mxPath(from="one",to=c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20"), free=F, value=0, arrows=1),
                           mxData(as.data.frame(modelData), type = "raw")
                           
)


#fitted_quadratic_lgcm2 <- mxRun(model_quadratic2)


pen.seq = 0.0001 * 2^(0:30)

pars.mat = matrix(NA,length(pen.seq),length(covariates)*3)
bic = rep(NA,length(pen.seq))

#free <- which(baseModel$A$free == T,arr.ind=T)

base.out <- outt$base.out #readRDS("mimic_base_AB.rds")

#saveRDS(list(modelData=modelData,cov_est=cov_est),"C:/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/modelData.rds")

test_fun <- function(seq){

    loc = as.numeric(base.out$A$free == T) # penalize all parameters
    # penalty value:
    penalty_value <- mxMatrix(type = "Full", nrow = 1,ncol = 1,free = F,values = pen.seq[seq],name = "penalty_value")
    ml_mat_nan = 1/base.out$A$values
    ml_mat_nan[is.nan(ml_mat_nan)] <- 0
    ml_mat_nan[is.infinite(ml_mat_nan)] <- 0
    
    ml_mat <- mxMatrix(type="Full",nrow = nrow(baseModel$A$values),ncol=ncol(baseModel$A$values),free=F,values=ml_mat_nan,name="ml_mat")
    selectedA <- mxMatrix(type = "Full", nrow = nrow(baseModel$A$values)*ncol(baseModel$A$values),ncol = 1,free = F,values = loc,name = "selectedA")
  
      penalty_combined <- mxAlgebra(1*penalty_value*(t(abs(cvectorize(ml_mat * baseModel.A)))%*%selectedA)+
                                      0*penalty_value*(t((cvectorize(baseModel.A))**2)%*%selectedA)
                                      , name = "penalty_combined") # elastic net
    
      t(abs(cvectorize(ml_mat$values * base.out$A$values))) %*% selectedA$values
    
    
    regfit_algebra <- mxAlgebra(baseModel.fitfunction + penalty_combined, name = "regfit_algebra")
    
    regfit_fun <- mxFitFunctionAlgebra("regfit_algebra")
    
    reg_Model <- mxModel(name = "reg_Model", submodel = baseModel,selectedA,
                         penalty_value, selectedA, penalty_combined,ml_mat,
                         regfit_algebra, # data_cov,
                         regfit_fun)
    reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No")
    reg_Model <- mxOption(reg_Model, "Standard Errors", "No")
    
    fit_reg_Model <- mxRun(reg_Model,silent=T)
    
    statistics = summary(fit_reg_Model)
    # have to calculate own bic
    #pars = round(fit_reg_Model$output$matrices$baseModel.A[loc],3)
    pars = round(fit_reg_Model$output$estimate,3)
    bic = statistics$Minus2LogLikelihood + log(nrow(modelData))*(statistics$ChiDoF-sum(pars==0))
    #pars.mat[j,] = pars
    
    
 # id = which(bic==min(bic))
  #par.ret[i,] = pars.mat[id,]
    return(list(pars,bic))
}


library(parallel)
no_cores <- detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, c("base.out","pen.seq","baseModel","modelData","cov_est")) # make sure to name baseModel in openmx

ret <- parLapply(cl,1:30,test_fun) # have to use parLapply if not passing grid, just vector

saveRDS(ret,"openmx_alasso30_ab2.rds")
