library(stats)

setwd("~/Dropbox/Ami-work/Cajal_school/")

# MRiShare sample data
master_morph <- read.table("sample_mrishare_morphometry.csv", header=TRUE, sep=",", dec=".")
head(master_morph)
summary(master_morph)

# MRiShare subject info
master_subinfo <- read.table("sample_mrishare_subinfo.csv", header=TRUE, sep=",", dec=".")
head(master_subinfo)
summary(master_subinfo)

#### Analysis plan ###########################################################
# First analyze male and female data separately, and in each group
# find the simplest model that explain the 10 variables of interest.
#  -> If the selected models are the same for each sex, test the 
#     Sex effect and interactions with the selected variables in combined group
#
# I will have a sample script doing this for the first metric (FS6_DKT_ct),
# but you can try it out with other metrics of your choice.
#
# Pay particular attention to whether your metric is influenced by eTIV, and
# how the methods to take it into account can modify your results.
##############################################################################
#### Data ####

# First combine two dataframes
master_dat <- merge(master_subinfo, master_morph,by="mrishare_id")

## Female data
female_dat <- master_dat[ which(master_dat$Sex=='F'), ]
female_dat$Age_c = as.numeric(scale(female_dat$Age, scale=F))
female_dat$eTIV_c = as.numeric(scale(female_dat$FS6_eTIV, scale=F))
female_dat$Sq_Age_c = as.numeric(female_dat$Age_c ^2)
summary(female_dat)

## Male data
male_dat <- master_dat[ which(master_dat$Sex=='M'), ]
male_dat$Age_c = as.numeric(scale(male_dat$Age, scale=F))
male_dat$eTIV_c = as.numeric(scale(male_dat$FS6_eTIV, scale=F))
male_dat$Sq_Age_c = as.numeric(male_dat$Age_c ^2)
summary(male_dat)

## Combined data
both_dat <-data.frame(master_dat)
both_dat$Age_c = as.numeric(scale(both_dat$Age, scale=F))
both_dat$eTIV_c = as.numeric(scale(both_dat$FS6_eTIV, scale=F))
both_dat$Sq_Age_c = as.numeric(both_dat$Age_c ^2)


#######################################################
###  1) Cortical thickness
#######################################################
dir.create("lm1_CT")

###############################
# STEP 1: Group specific tests
###############################
# Use centered Age/eTIV
#################
# 1-1) Female
#################

# Test with most complex model you are willing to test
lm1_CT_f <- lm( FS6_DKT_T_ct ~ Age_c + Sq_Age_c + eTIV_c + Age_c:eTIV_c + 
                    Sq_Age_c:eTIV_c, data=female_dat)
lm1_CT_f_step <- step(lm1_CT_f, scope=list(lower=FS6_DKT_T_ct ~ Age_c, upper=lm1_CT_f))
summary(lm1_CT_f_step)

# This means that in females, effect of CT on age may be better explained by 
# quadratic model and that there may be differences between those who have high or low
# eTIV. We will visualize this effect below, but try removing interaction effect to
# see how much model fit (sq.R value changes).

lm1_CT_f_noInt <- lm(FS6_DKT_T_ct ~ Age_c + Sq_Age_c + eTIV_c, data=female_dat)
summary(lm1_CT_f_noInt)

# step test for this model
lm1_CT_f_noInt_step <- step(lm1_CT_f_noInt, scope=list(lower=FS6_DKT_T_ct ~ Age_c, upper=lm1_CT_f_noInt)) 
summary(lm1_CT_f_noInt_step)

# W/o allowing for the interaction, the models do not seem to explain data much.

### Save the model results
sink("lm1_CT/lm1_CT_f_step.txt")
print(step(lm1_CT_f, scope=list(lower=FS6_DKT_T_ct ~ Age_c, upper=lm1_CT_f)))
sink()

sink("lm1_CT/lm1_CT_f_step_model.txt")
print(summary(lm1_CT_f_step))
sink()

sink("lm1_CT/lm1_CT_f_noInt_model.txt")
print(summary(lm1_CT_f_noInt))
sink()

sink("lm1_CT/lm1_CT_f_noInt_step_model.txt")
print(summary(lm1_CT_f_noInt_step))
sink()

### Plot predictions by the first model
# 1) Age effect on CT

# To visualize the effect, you can predict CT values for any 
# theoretical sex/age/eTIV. We will create a df to be used for 
# predicting CT at different age at fixed eTIV, or vice versa.

# let's say we will predict CT at the range of age in this dataset,
# bet 18~28
agerange <- seq(18, 28, 0.2)

# df to use for predicting age effect at mean eTIV for females
f_pred_dat_mean_eTIV <- data.frame(Age=agerange,
                         Age_c=agerange - mean(female_dat$Age),
                         eTIV=mean(female_dat$eTIV),
                         eTIV_c=mean(female_dat$eTIV_c))
f_pred_dat_mean_eTIV$Sq_Age_c <- (f_pred_dat_mean_eTIV$Age_c)^2

# To look at interaction between age and eTIV, let's also create 
# another df for hypothetical female data with low eTIV
f_pred_dat_low_eTIV <- data.frame(Age=agerange,
                                   Age_c=agerange - mean(female_dat$Age),
                                   eTIV=mean(female_dat$eTIV) - 2*sd(female_dat$eTIV))
f_pred_dat_low_eTIV$Sq_Age_c <- (f_pred_dat_low_eTIV$Age_c)^2
f_pred_dat_low_eTIV$eTIV_c <- f_pred_dat_low_eTIV$eTIV - mean(female_dat$eTIV)

# another with high eTIV
f_pred_dat_high_eTIV <- data.frame(Age=agerange,
                                  Age_c=agerange - mean(female_dat$Age),
                                  eTIV=mean(female_dat$eTIV) + 2*sd(female_dat$eTIV))
f_pred_dat_high_eTIV$Sq_Age_c <- (f_pred_dat_high_eTIV$Age_c)^2
f_pred_dat_high_eTIV$eTIV_c <- f_pred_dat_high_eTIV$eTIV - mean(female_dat$eTIV)

# Now predict the CT using the 3 DF above
f_CT_pred_mean_eTIV <- predict(lm1_CT_f_step, f_pred_dat_mean_eTIV)
f_CT_pred_low_eTIV <- predict(lm1_CT_f_step, f_pred_dat_low_eTIV)
f_CT_pred_high_eTIV <- predict(lm1_CT_f_step, f_pred_dat_high_eTIV)

# Now that we have predicted values for 3 set of hypothetical data,
# plot them on top of scatter plot of CT vs Age

par(mar=c(4, 4, 4, 4) + 0.1)
plot(female_dat$Age, female_dat$FS6_DKT_T_ct,
     col=adjustcolor("coral1", alpha.f = 0.5), pch=16 , cex=1.3,
     main="Cortical thickness change with age in females",
     xlab="Age (years)", ylab="Mean cortical thickness (mm)")

# Lines of prediction
lines(agerange, f_CT_pred_mean_eTIV, col="firebrick4", lwd=2)
lines(agerange, f_CT_pred_low_eTIV, col="firebrick4", lty=3, lwd=2)
lines(agerange, f_CT_pred_high_eTIV, col="firebrick4", lty=2, lwd=2)

text(20, 2.95 , pos=4, col="firebrick4", cex=0.8,
     paste("Model : ", round(lm1_CT_f_step$coefficients[1],3) , " + " ,
           round(lm1_CT_f_step$coefficients[2], 3) , "*Age"  , "+" ,
           round(lm1_CT_f_step$coefficients[3], 3) , "*Age^2", "+" ,
           round(lm1_CT_f_step$coefficients[4]/1000, 3) , "*eTIV", "+" , # per eTIV in cc
           round(lm1_CT_f_step$coefficients[5]/1000, 3) , "*eTIV X Age^2", "+" ,
           "\n" , "squared R adjusted = ",round(summary(lm1_CT_f_step)$adj.r.squared, 3)))


# 1) The relationship between CT and eTIV
# Can you make more predictions, this time at fixed age (mean, low, high), at 
# a hypothetical (but realistic) range of eTIV?

# Plot this prediction lines on top of the scatter plot of CT vs eTIV (hint:
# You just have to change the x axis...)




################
# 1-2) Male
################

# Can you repeat the analysis for male data? (It's OK to copy and modify...!)


### Save the model results


### Plot predictions by the first model
# 1) Age effect on CT

# 2) eTIV and CT


###############################
# STEP 2: Combined group
###############################
# Now, run the analysis with
# entire data, including Sex as a categorical variable.


### Save the model results

### Plot prediction
# First create df for prediction


# Now plot

