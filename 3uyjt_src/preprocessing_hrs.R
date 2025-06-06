library(tidyverse)
library(haven)

# load the SPSS data file (provided by Ross)

if (!exists("rndhrs"))
  if (Sys.info()["login"]=="brandmaier") {
    # it's Andy!
    cat("Welcome, Andy!\n")
    rndhrs <- haven::read_spss("rndhrs_o.sav")
  } else {
    # it's Ross!
    cat("Welcome, Ross!\n")
    rndhrs <- read_spss("/Volumes/GoogleDrive/Shared drives/hrs_paper/rndhrs_o.sav")
  }

# filter cohorts to include
#  0.Hrs/Ahead ovrlap
#             1.Ahead
#              2.Coda
#               3.Hrs
#         4.WarBabies
# AND
# to include people at least 50 at Wave 4 and younger than 88 at Wave 11 
# TODO: be more inclusive?, also this includes people outside this age range because people
# can be older than 88 at Wave 4 and not have data at Wave 11
rndhrs_subset <- rndhrs %>% 
  filter(HACOHORT %in% c(0,1,2,3,4))  %>% 
  # select age range
  filter(R4AGEY_B > 50) %>% 
  filter(R11AGEY_B < 88) 

# create a centered age variable, centered at mean age at Wave 4
cage <-  round(mean(rndhrs$R4AGEY_B,na.rm=TRUE))
rndhrs_subset <- rndhrs_subset %>%
  mutate(CAGE_T4 = R4AGEY_B-cage,
         CAGE_T5 = R5AGEY_B-cage,
         CAGE_T6 = R6AGEY_B-cage,
         CAGE_T7 = R7AGEY_B-cage,
         CAGE_T8 = R8AGEY_B-cage,
         CAGE_T9 = R9AGEY_B-cage,
         CAGE_T10 = R10AGEY_B-cage,
         CAGE_T11 = R11AGEY_B-cage
  )

# create centered age^2 variable for quadratic LGCM
qscale <- 100
rndhrs_subset <- rndhrs_subset %>%
  mutate(CAGE2_T4 = CAGE_T4^2 / qscale,
         CAGE2_T5 = CAGE_T5^2 / qscale,
         CAGE2_T6 = CAGE_T6^2 / qscale,
         CAGE2_T7 = CAGE_T7^2 / qscale,
         CAGE2_T8 = CAGE_T8^2 / qscale,
         CAGE2_T9 = CAGE_T9^2 / qscale,
         CAGE2_T10 = CAGE_T10^2 / qscale,
         CAGE2_T11 = CAGE_T11^2 / qscale
  )

# only for plotting:
# select only time points and TR20 outcomes
# take a random subset and convert from wide to long
rndhrs_subset_plot <- rndhrs_subset %>%  
  select(HHIDPN, R4AGEY_B:R11AGEY_B, R4TR20, R5TR20, R6TR20, R7TR20, R8TR20, R9TR20, R10TR20,R11TR20) %>% 
  sample_n(40) %>% 
  pivot_longer(c(R4AGEY_B:R11AGEY_B,R4TR20:R11TR20),
               names_to=c("Wave",".value"),
               names_pattern = "R(\\d*)(.*)")


# plot the random subset
rndhrs_subset_plot %>%
  ggplot(aes(x=AGEY_B,y=TR20,group=factor(HHIDPN)))+
  geom_line()+
  theme_light()

# because definitions variable cannot be NA, recode those to a fake value
# and set the corresponding observed variable to NA
# this way, we can estimate LGCMs with definition variables and FIML
rndhrs_subset[is.na(rndhrs_subset$CAGE_T4), "R4TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T4), "CAGE_T4"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T4), "CAGE2_T4"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T5), "R5TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T5), "CAGE_T5"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T5), "CAGE2_T5"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T6), "R6TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T6), "CAGE_T6"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T6), "CAGE2_T6"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T7), "R7TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T7), "CAGE_T7"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T7), "CAGE2_T7"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T8), "R8TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T8), "CAGE_T8"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T8), "CAGE2_T8"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T9), "R9TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T9), "CAGE_T9"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T9), "CAGE2_T9"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T10), "R10TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T10), "CAGE_T10"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T10), "CAGE2_T10"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE_T11), "R11TR20"] <- NA
rndhrs_subset[is.na(rndhrs_subset$CAGE_T11), "CAGE_T11"] <- -99999
rndhrs_subset[is.na(rndhrs_subset$CAGE2_T11), "CAGE2_T11"] <- -99999


vars_AB2 <-c("RAEDYRS","RAGENDER","R2SHLT","R2HOSP", "RARACEM",  "RAHISPAN",
             "R2SMOKEV","R3CHOLST","R3FLUSHT","R2DRINK","R2CHAIRA","R2BACK","R2JOGA")


save(vars_AB2, file="results/varlist.Rda")


# some conversions
rndhrs_subset$RAGENDER <- factor(ifelse(rndhrs_subset$RAGENDER==1,"male","female"))
rndhrs_subset$RAEDYRS <- factor(rndhrs_subset$RAEDYRS, ordered = TRUE) 
rndhrs_subset$RARACE <- factor(rndhrs_subset$RARACE, ordered=FALSE)
rndhrs_subset$RAHISPAN <- factor(rndhrs_subset$RAHISPAN, ordered=FALSE)


# and once again with only a subset of covariates for testing
rndhrs_subset_tree_ab2 <- rndhrs_subset %>%  select(
  R4AGEY_B:R11AGEY_B, R4TR20, R5TR20, R6TR20, R7TR20, 
  R8TR20, R9TR20, R10TR20,R11TR20,CAGE_T4:CAGE_T11, CAGE2_T4:CAGE2_T11,
  !!!vars_AB2)  # !!! is the awesome unmasking operator

save(rndhrs_subset_tree_ab2, file="rndhrs_subset_tree_ab2.sav")
