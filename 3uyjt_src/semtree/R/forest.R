num_forests <- 7

# need to run tree.R before forest.R to
# establish template model and load data

# selection of potential predictors for tree and forest
vars_AB <- c("RAEDYRS",
             
             "RAGENDER",
             "R2SHLT","R2HOSP","R2NRSHOM","R2DOCTOR","R2HOMCAR","R2WALKRA", # "R2HLTHLM"
             "R2DRESSA","R2BATHA","R2EATA","R2BEDA","R2PHONEA", # ,"R2TOILTA""R2MAPA","R2CALCA"
             "R2MONEYA","R2MEDSA","R2SMOKEV",
             "R2PUSH")
#
# initiate parallel execution
#

#future::plan("sequential")
future::plan("multisession")

#
# run forest with SEM trees
# this will take ages - check out the faster score-based trees!
#
hrsforest <- semtree::semforest(model=fitted_quadratic_lgcm, 
                                data=as.data.frame(rndhrs_subset), 
                                covariates=vars_AB,
                                control=semforest.control(num.trees=num_forests)
                                )
# compute variable importance
hrsvim <- varimp(hrsforest)

# save variable importance
save(hrsvim, file="results/hrs-vim.sav")


# plot vim
plot(hrsvim)