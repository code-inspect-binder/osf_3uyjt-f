library(regsem);library(OpenMx)
# run preprocessing_hrs.R first

df <- rndhrs_subset

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

#modelData[,66:ncol(modelData)] = sub_dat

cov.mat = cov(modelData[,vars_AB],use="pairwise.complete.obs")
cor.mat = cor(modelData[,vars_AB],use="pairwise.complete.obs")
cov_est = cov.mat[lower.tri(cov.mat,diag=T)]



# code gets MLE estimates from OpenMx
# include covariates

indicators<-c("R4TR20","R5TR20","R6TR20","R7TR20","R8TR20","R9TR20","R10TR20","R11TR20")
covariates = c(vars_AB)
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


fitted_quadratic_lgcm2 <- mxRun(baseModel)


saveRDS(fitted_quadratic_lgcm2,"/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/mimic_base_AB.rds")


base.out <- readRDS("/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/mimic_base_AB.rds")

(1/base.out$A$values) * base.out$A$values

summary(fitted_quadratic_lgcm2)


pen.seq = c(0,(seq(1,sqrt(nrow(modelData)),length.out=30)**2))

pars.mat = matrix(NA,length(pen.seq),length(covariates)*3)
bic = rep(NA,length(pen.seq))

#free <- which(baseModel$A$free == T,arr.ind=T)

#base.out <- readRDS("C:/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/mimic_base.rds")
saveRDS(list(modelData=modelData,cov_est=cov_est),"C:/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/modelData_all.rds")

#base.out$A$values[loc]

openmx_lasso30_all <- readRDS("~/GitHub/handbook_of_sem_chapter/R/openmx_lasso30_all.rds")


ml_coefs = abs(openmx_lasso30_all[[1]][[1]][1:150])

test_fun <- function(seq){


    loc = as.numeric(baseModel$A$free == T) # penalize all parameters
    # penalty value:
    penalty_value <- mxMatrix(type = "Full", nrow = 1,ncol = 1,free = F,values = pen.seq[seq],name = "penalty_value")
    selectedA <- mxMatrix(type = "Full", nrow = nrow(baseModel$A$values)*ncol(baseModel$A$values),ncol = 1,free = F,values = loc,name = "selectedA")
  
    ml_mat <- mxMatrix(type="Full",nrow = nrow(baseModel$A$values),ncol=ncol(baseModel$A$values),free=F,values=1/base.out$A$values,name="ml_mat")
    #sum((t(abs(cvectorize(baseModel.A[loc]))))*(1/ml_coefs))
    

    
      penalty_combined <- mxAlgebra(penalty_value*sum((t(abs(cvectorize(baseModel.A[selectedA]))))*(1/ml_coefs)), name = "penalty_combined") # elastic net
    
    
    
    regfit_algebra <- mxAlgebra(baseModel.fitfunction + penalty_combined, name = "regfit_algebra")
    
    regfit_fun <- mxFitFunctionAlgebra("regfit_algebra")
    
    reg_Model <- mxModel(name = "reg_Model", submodel = baseModel,selectedA,
                         penalty_value, selectedA, penalty_combined,
                         regfit_algebra, # data_cov,
                         regfit_fun)
    reg_Model <- mxOption(reg_Model, "Calculate Hessian", "No")
    reg_Model <- mxOption(reg_Model, "Standard Errors", "No")
    
    fit_reg_Model <- mxRun(reg_Model,silent=T)
    
    statistics = summary(fit_reg_Model)
    # have to calculate own bic
    pars = round(fit_reg_Model$output$matrices$baseModel.A[loc],3)
    bic = statistics$Minus2LogLikelihood + log(nrow(modelData))*(statistics$ChiDoF-sum(pars==0))
    #pars.mat[j,] = pars
    
    
 # id = which(bic==min(bic))
  #par.ret[i,] = pars.mat[id,]
    return(list(pars,bic))
}


library(parallel)
no_cores <- 1#detectCores()
cl <- makeCluster(no_cores)
clusterExport(cl, c("ml_coefs","pen.seq","baseModel")) # make sure to name baseModel in openmx

ret <- parLapply(cl,1:1,test_fun) # have to use parLapply if not passing grid, just vector


str(openmx_alasso30_ab2)
openmx_alasso30[[30]]

bic = rep(NA,30)

for(i in 1:30){bic[i] = openmx_alasso30_ab2[[i]][[2]]}

which(bic == min(bic))
length(openmx_alasso30_ab2[[17]][[1]])

mat_alasso_pars = matrix(openmx_alasso30_ab2[[17]][[1]][1:42],14,3)
rownames(mat_alasso_pars) <- vars_AB22
colnames(mat_alasso_pars) <- c("icept","slope","quad")
saveRDS(mat_alasso_pars,"/Users/rjacobuc/Documents/Github/handbook_of_sem_chapter/R/alasso_pars_mat_ab2.rds")

bic[bic > 600000] <- NA


tiff("bic_ab2.tiff",res=300,width=6,height=4,units="in")
plot(1:30,bic,xlab="Lambda",xaxt="n",yaxt="n",ylab="BIC")
axis(1,at=c(1,4,9,14,19,24,29),labels=round(0.0001 * 2^(0:29),3)[c(1,4,9,14,19,24,29)],cex.axis=.8)
abline(v=17,lty=3)
axis(2,cex.axis=.6)
dev.off()
pars.mat <- matrix(NA,30,42)
for(i in 1:30){
  pars.mat[i,] = openmx_alasso30_ab2[[i]][[1]][1:42]
}

pars.mat[27,] <- 0

tiff("alasso_ab2_pars.tiff",res=300,width=6,height=4,units="in")
plot(1:30,pars.mat[,1],type="l",ylim = c(min(pars.mat) * 0.95,
                                         max(pars.mat) * 1.05),ylab="Coefficient",xlab="Lambda",xaxt="n")
abline(a=0,b=0)
for(j in 2:ncol(pars.mat)){
  lines(1:30,pars.mat[,j])
}
axis(1,at=c(1,4,9,14,19,24,29),labels=round(0.0001 * 2^(0:29),3)[c(1,4,9,14,19,24,29)],cex.axis=.8)
abline(v=17,lty=3)

