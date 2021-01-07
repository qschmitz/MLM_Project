###############################################################################
#-------------------------Mixed Linear Models Project-------------------------#
#--------------------------------------by-------------------------------------#
#--------------------------------------------<Please put your name here>------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------Hervégil Voegeli-----------------#
###############################################################################
#Here is the ordered code, how it will be presented in our report
#to not miss all our trials&errors, we can do it in the Fun.R file
# Research questions:
#   1)if the in-vitro sensitivity score can be used to predict the treatment effect on the patient at the hospital
#   2)if the in-vitro sensitivity score decreases with time and is stronger for line 1 than for lines 2 and 3+ 

# Initialization  ------------------------------------------------------------
#Ctrl+Shift+R to create a new section
#Ctrl+Alt+T to run only the current section
#Ctrl+Shift+Enter to run all file
#Ctrl-Shift-C to comment/uncomment quickly
library(scales)#This library allows to control the nuance of colors with the function alpha()
library(nlme) #Groupdata
library(lmerTest)
library(MuMIn)
library(pbkrtest)
library(car)
library(mvtnorm)
library(dplyr)
library(circlize)
library(lmtest)
Chemo = read.table(file="Chemotherapy-version1.csv",header=TRUE,sep=",")
attach(Chemo)

Chemo[,"line"] <- as.factor(Chemo[,"line"]);
Chemo[,"month"]  <- as.factor(Chemo[,"month"]);
Chemo[,"month"]  <- as.integer(Chemo[,"month"]);
sapply(Chemo,class)
Chemo[,"line"] <- as.factor(Chemo[,"line"])

#Maybe we do this only to answer the 2nd question?
#merge levels 3,4 and 5 of line in one level (+3) as asked in the 2nd question
levels(Chemo[,"line"])<-c("1","2","+3","4","5")
levels(Chemo[,"line"])[levels(Chemo[,"line"])=="4"] <-"+3"
levels(Chemo[,"line"])[levels(Chemo[,"line"])=="5"] <-"+3"

#re-indexing months
chemo=Chemo
for (i in 1:19){
  chemo[patient==LETTERS[i],"month"]=seq(1,length(chemo[patient==LETTERS[i],"month"]))
}

# 4.1 First look of the data ----------------------------------------------

# • Define:
#– what are the independent experiment units (subjects, households, lots, ...)
  #We have only within-subject factors
  #Patients are the independent experiment units
#– the role of the different variables in your analysis:
#   ∗ which variable is the response ?
      #tumour
#   ∗ which variables are covariates ?
      # all the others
#   ∗ which variables are factors ?
      #patient, line   
#   ∗ Are the covariates and factors fixed or random effects ? (justify)
    #R.E. : patients (we do not exhaust the population, they do not brings explaination,levels chosen randomly)
    #F.E. : line (regressor), month (nested in line), sensitivity
# ∗ Are the factors crossed, nested, ... ?
    # months is nested in line
# • Check if the class of your variables correspond to your needs:
#   – class numeric for the response:
        #tumour
# – class numeric or integer covariates like age, ... : 
        # month,sensitivity
# – class factor or ordered for factors like subject, gender, ...:
        # patient, line
# • Write a first version of your model that include potential interactions
#patient: i, line: l, month: m, sensitivity: s
  #Y_{ilms} = \mu + \alpha_l \+ \beta_s + (\alpha_l \cdot \beta_s) + \delta_m + \gamma_s + zeta_i +epsilon_ilms


#month produces a zero variance matrix-> we should avoid it
#and such result is logical->month is not indexed correctly, it is set as a repeating loop
#but line does not produce such singular matrix and is a better measure of time I believe

lmer_Chemo = lmer(tumour~line*sensitivity+month*line+(month|patient),data=Chemo)
  summary(lmer_Chemo)
lme_Chemo = lme(fixed=tumour~line*sensitivity+month*line, random=~1|patient,data=Chemo, method="ML")
  summary(lme_Chemo)

lmer_Chemo = lmer(tumour~sensitivity*line+month*line+(1|patient),data=Chemo)
  summary(lmer_Chemo)
lme_Chemo = lme(fixed=tumour~sensitivity*line+month*line, random=~1|patient,data=Chemo, method="ML")
  summary(lme_Chemo)
# the variance of the random slope is almost null...is it really a random slope?

# 4.2 Exploratory data analysis ----------------------------------------------

# • check if the model’s assumptions hold:
#   – linearity of the relation between the covariates and the response
source('sections/4.2.linearity.r')

# normality of random effects and residuals
  #almost normal here
qqnorm(ranef(lme_Chemo)[,1],pch=16);qqline(ranef(lme_Chemo)[,1],col=2,lwd=2,lty=2)
  #we need fatter tails for the residuals
qqnorm(lme_Chemo$residuals,pch=16);qqline(lme_Chemo$residuals,col=2,lwd=2,lty=2)

#sensitivity and tumour should be used as response variable
#clear interaction between month and line to determine sensitivity or tumour
interaction.plot(line,month,tumour)
interaction.plot(month,line,tumour)
#interaction between line and sensitivity? the interaction plot are useless here
# interaction.plot(line,sensitivity,tumour)
# interaction.plot(sensitivity,line,tumour)

#but what interactions should we keep?
lme_Chemo_no_interact = lme(fixed=tumour~sensitivity+month+line, random=~1|patient,data=Chemo, method="ML")
lme_Chemo_interact_month = lme(fixed=tumour~sensitivity+line*month, random=~1|patient,data=Chemo, method="ML")
lme_Chemo_interact_sensitivity = lme(fixed=tumour~line*sensitivity+month, random=~1|patient,data=Chemo, method="ML")

anova(lme_Chemo_interact_month,lme_Chemo)#we do not get a better model with the interaction between sensitivity and line
anova(lme_Chemo_interact_sensitivity,lme_Chemo)#but a better one by having month*line interaction


# and having interaction month*line is better than no interaction at all
anova(lme_Chemo_no_interact,lme_Chemo_interact_month)

#we should keep interaction between month and line and retrieve the other.

# – equal variance of the responses around the different factor’s levels

plot.design(Chemo)

    # to equal variance is to use the ANOVA tests
anova(lme_Chemo,type="sequential")#effect of sequentially adding each new element
    # sensitivity, its interaction and the interaction of line and month not significant
Anova(lme_Chemo,type="II")#effect of each element to the whole
    # all become more significant, line:sensitivity is statistically significant
Anova(lme_Chemo,type="III")#effect of each element to the whole (with interactions)
    #they are all significant apart from interaction between line and month.
#CONCLUSION: get interaction of line and month out of the model

#1) Overall data no clear heteroskedasticity, but kinda different variances per line
plot(lme_Chemo,resid(.,type="p")~fitted(.),abline=0,id=0.05)
plot(lme_Chemo,resid(.,type="p")~fitted(.)|line,abline=0,id=0.05)
# robust needed for patients A,P,C and maybe K

#2) Breusch-Pagan test: we cannot reject H_0: Homoskedasticity
bptest(lme_Chemo)

#3) We see the same residual patterns for all month
Chemo.gpdata = groupedData(tumour~sensitivity|month,data=Chemo,order.groups=F,
                            label=list(x="Sensitivity",y="Tumour response"),
                            units=list(y="(percent by weight)"))
plot(Chemo.gpdata)
#Patients do not have the same sensitivity evolution at all!
#Not the same reaction pattern from sensitivity to tumour for each patient
Chemo.gpdata = groupedData(tumour~sensitivity|patient,data=Chemo,order.groups=F,
                           label=list(x="Sensitivity",y="Tumour response"),
                           units=list(y="(percent by weight)"))
plot(Chemo.gpdata)
#Tumour cannot serve as random slope for sensitivity
Chemo.gpdata = groupedData(sensitivity~tumour|patient,data=Chemo,order.groups=F,
                           label=list(x="Tumour",y="Sensitivity"),
                           units=list(y="(percent by weight)"))
plot(Chemo.gpdata)
#Month looks like a candidate for a random slope but with a lot of disturbance
Chemo.gpdata = groupedData(tumour~month|patient,data=Chemo,order.groups=F,
                           label=list(x="Tumour",y="Sensitivity"),
                           units=list(y="(percent by weight)"))
plot(Chemo.gpdata)
#No patients seems strange
#Line can be used as random slope:
Chemo.gpdata = groupedData(tumour~line|patient,data=Chemo,order.groups=F,
                           label=list(x="Tumour",y="Sensitivity"),
                           units=list(y="(percent by weight)"))
plot(Chemo.gpdata,outer =~line)
# 4) Boxplots shows line, month and patient heteroskedasticity

# – presence or absence of interactions between the different “explanatory variables”
boxplot(tumour~patient,main="Tumour responses per patient")#is 19 enough to diagnose patient heteroskedasticity?
abline(h=mean(tumour),col="red")
boxplot(tumour~month,data = Chemo,main="Tumour responses per month of treatment")#that's more accurate!
abline(h=mean(tumour),col="red")
boxplot(tumour~month,data = chemo,main="Tumour responses per month of treatment")#that's more accurate!
abline(h=mean(tumour),col="red")
boxplot(tumour~line,main="Tumour response per line")#which is consistent with less effect of lines through time, less variance
abline(h=mean(tumour),col="red")
boxplot(sensitivity~patient)#clear heteroskedasticity of sensitivity distribution through patients, is it correlated to survival or anything else?
abline(h=mean(sensitivity),col="red")

# – correlated or uncorrelated random effects ?
#   • get an impression of the variables potentially explanatory
#we have only one random effect, hence no possible correlations
# pairs(lme_Chemo,~ranef(.)|patient,id=0.05,grid=T)




# 4.3 Model building and hypotheses ----------------------------------------------
# • Rewrite your model to take into account the results of the exploratory data analysis. If
# you have a doubt about the presence of a parameter (interaction between fixed factors,
#                                                     correlation between random effects), include it in the full model. It’s pertinence in
# the model will be tested.
# • Define precisely the hypotheses you want to test depending on the questions you wish
# to answer (Are you interested on the main effect of a factor or on its marginal effect
#            ?).

# We need to add:
  # heteroskedasticity for line
  # robust estimates
  # we should keep interaction between month and line and retrieve the other.


# 4.4 Model estimation ----------------------------------------------

# • Estimate the postulated model.

# 4.5 Diagnostic analyses of the full model ----------------------------------------------

# Before analyzing the results of your estimations, you have to determine the model is valid
# • check the residuals : are they iid N(0, σ2) ?
#   • check the random effects : are they iid N(0, σ2) ?
#   If not, try
# • a scale’s transformation of some variable (log of the response in the case of heteroscedasticity for example)
# • a quadratic form if non linearity of the influence of a covariate is suspected
# • to relax the assumption of equal variance of the observations around the levels of a
# factor
# • ...
# 4.6 Model selection ----------------------------------------------

# Once the model is estimated and its validity checked, you may wish to determine with
# likelihood ratio tests if some of its parameters are significant or not with the aim of selecting
# a more parsimonious model. Suggestions:
#   • if the variance of random effect is small compared to the residual variance, this random
# effect may not be useful.
# • if the means of the responses for the levels of a factor seem the same, a model without
# that fixed factor may be more adequate.
# • ...
# The selected model should be re-checked as in 2.5.
# 4.7 Hypotheses
# Once you have selected the appropriate model, test the hypotheses you defined earlier.




