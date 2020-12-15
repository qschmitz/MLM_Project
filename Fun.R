#Simply having fun with the data, but keeping it here
library(scales)#This library allows to control the nuance of colors with the function alpha()
library(nlme) #Groupdata
library(lmerTest)
library(MuMIn)
library(pbkrtest)
library(car)
library(mvtnorm)
library(dplyr)
setwd("E:/Universite de Genève/Faculté d'économie/Mon Master en statistique/A1S1/Mixed Linear Models/DATA exam/MLM_Project")


#What are we dealing with
Chemo = read.table(file="Chemotherapy-version1.csv",header=TRUE,sep=",")
head(Chemo);View(Chemo);#19 patients
attach(Chemo)

#check data type
sapply(Chemo, class)#line should be factor
Chemo[,"line"] <- as.factor(Chemo[,"line"]);
#also month?
# Chemo[,"month"] <- as.factor(Chemo[,"month"])

#study data with plots
  #Boxplots: use it for heteroskedasticity, distribution
boxplot(tumour~patient)
boxplot(tumour~month)
boxplot(tumour~line)
boxplot(tumour~sensitivity)

boxplot(patient~month)
boxplot(patient~line)
boxplot(patient~sensitivity)

boxplot(month~line)
boxplot(month~sensitivity)

boxplot(line~sensitivity)

  #Plot design: shows variability of your factors, i.e. heteroskedasticity
plot.design(Chemo)

  #Interaction plot: does one marginal effect depend on its relative position on the characteristics space?
interaction.plot(sensitivity,month,tumour)#is sensitivity affecting the month effect on tumour?
interaction.plot(month,sensitivity,tumour)#month has no interaction with sensitivity
interaction.plot(month,patient,tumour)#neither has it on any patient
interaction.plot(month,line,tumour)#HA-HAAAAA MONTH HAS AN INTERACTION WITH LINE 5

interaction.plot(line,patient,tumour)#line has an increasing effect on the patient effect on tumour

#setting contrasts sum for line, month: there is not a reference line or month, so sum contrasts
contrasts(Chemo$line)="contr.sum"
contrasts(Chemo$month)="contr.sum"

#setting model
#with month as r.v. I get 0 variance...so it is fixed? but in Orth data, age was a r.v. ...
lmer_Chemo = lmer(tumour~line*sensitivity+month+(1|patient),data=Chemo)
lme_Chemo = lme(fixed=tumour~line*sensitivity+month, random=~1|patient,data=Chemo, method="ML")


