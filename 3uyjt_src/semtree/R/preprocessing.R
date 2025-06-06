library(tidyverse)
library(haven)

# load the SPSS data file (provided by Ross Jacobucci)
# this has a size of about 500 Megabyte and takes a while
# to load
rndhrs <- haven::read_spss("data/rndhrs_o.sav")


# filter cohorts to include
#  0.Hrs/Ahead ovrlap
#             1.Ahead
#              2.Coda
#               3.Hrs
#         4.WarBabies
# AND
# to include people at least 50 at Wave 4 and younger than 88 at Wave 11 

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
# rescale by qscale to avoid huge numbers
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



# selection of potential predictors for tree and forest
vars_AB <- c("RAEDYRS",
             
             "RAGENDER",
             "R2SHLT","R2HOSP","R2NRSHOM","R2DOCTOR","R2HOMCAR","R2WALKRA", # "R2HLTHLM"
             "R2DRESSA","R2BATHA","R2EATA","R2BEDA","R2PHONEA", # ,"R2TOILTA""R2MAPA","R2CALCA"
             "R2MONEYA","R2MEDSA","R2SMOKEV",
             "R2PUSH")

#R4AGEY_B:R11AGEY_B, R4TR20, R5TR20, R6TR20, R7TR20, 
#R8TR20, R9TR20, R10TR20,R11TR20,CAGE_T4:CAGE_T11, CAGE2_T4:CAGE2_T11

vars_select <- c(vars_AB, "R4TR20","R5TR20","R6TR20",
                 "R7TR20","R8TR20","R9TR20","R10TR20","R11TR20",
                 paste0("CAGE_T",4:11),paste0("CAGE2_T",4:11),
                 paste0("R",4:11,"AGEY_B")
)


# make subset selection
rndhrs_subset <- rndhrs_subset %>% select(all_of(vars_select))

# save subset
write.csv(rndhrs_subset, file="data/rndhrs_subset.csv", 
          quote=FALSE, row.names=FALSE)


