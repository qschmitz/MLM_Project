###############################################################################
#-------------------------Mixed Linear Models Project-------------------------#
#--------------------------------------by-------------------------------------#
#--------------------------------------------<You can put your name here>-----#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------Hervégil Voegeli-----------------#
###############################################################################

# Initialization  ------------------------------------------------------------
#Ctrl+Shift+R to create a new section
#Ctrl+Alt+T to run only the current section
#Ctrl+Shift+Enter to run all file

#Simply having fun with the data, but keeping it here
library(scales)#This library allows to control the nuance of colors with the function alpha()
library(nlme) #Groupdata
library(lmerTest)
library(MuMIn)
library(pbkrtest)
library(car)
library(mvtnorm)
library(dplyr)
library(circlize)
setwd("E:/Universite de Genève/Faculté d'économie/Mon Master en statistique/A1S1/Mixed Linear Models/DATA exam/MLM_Project")


# Overview data -----------------------------------------------------------

#What are we dealing with
Chemo = read.table(file="Chemotherapy-version1.csv",header=TRUE,sep=",")
head(Chemo);#View(Chemo);#19 patients

attach(Chemo)


# Studying variables -------------------------------------------------------

#check data type
sapply(Chemo, class)#line should be factor
Chemo[,"line"] <- as.factor(Chemo[,"line"]);
#also month?->ofc not
# Chemo[,"month"] <- as.factor(Chemo[,"month"])


# Overview tumour/sensitivity -------------------------------------------------------------------

#I want to see the evolution of the tumour through time:
#month should be re-indexed and restart to 0 for each new patient
for (i in 1:19){
  Chemo[patient==LETTERS[i],"month"]=seq(1,length(Chemo[patient==LETTERS[i],"month"]))
}
#we can plot the tumour growth
plot(Chemo[patient==LETTERS[1],"month"],Chemo[patient==LETTERS[1],"tumour"],
     type="l",main="Tumour growth through time",
     xlab = "Time in months",ylab = "Tumour evolution",
     xlim=c(1,14),ylim=c(-3,1))
for (i in 2:19){
  lines(Chemo[patient==LETTERS[i],"month"],Chemo[patient==LETTERS[i],"tumour"],type="l",col=rand_color(1))
}
#and display the mean for the sake of clarity
moyenne=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    moyenne[j]=moyenne[j]+Chemo[patient==LETTERS[i],"tumour"][j]/19
  }
}
lines(seq(14),moyenne,type="l",col="red",lwd=5)
text(5,-0.25,labels="Mean",col="red")

#J-S mean
#first, we need its variance
variance=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    variance[j]=variance[j]+(Chemo[patient==LETTERS[i],"tumour"][j]-moyenne[j])^2/19
  }
}
#then we compute it
moyenne_JS=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    #I directly sum the p_i_JS and take their mean
    moyenne_JS[j]=moyenne_JS[j]+ (moyenne[j] + (1-16*(moyenne[j]*(1-moyenne[j])/19)/variance[j])*(Chemo[patient==LETTERS[i],"tumour"][j]-moyenne[j]))/19
  }
}
lines(seq(14),moyenne_JS,type="l",col="blue",lwd=5)
text(5,-1.5,labels="James Stein Mean",col="blue")



#Now the same with the evolution of sensitivity
plot(Chemo[patient==LETTERS[1],"month"],Chemo[patient==LETTERS[1],"sensitivity"],
     type="l",main="Sensitivity to the drug combination of interest",
     xlab = "Time in months",ylab = "Sensitivity",
     xlim=c(0,14),ylim =c(-0.1,1))
for (i in 2:19){
  lines(Chemo[patient==LETTERS[i],"month"],Chemo[patient==LETTERS[i],"sensitivity"],type="l",col=rand_color(1))
}
#and display the mean for the sake of clarity
moyenne_sens=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    moyenne_sens[j]=moyenne_sens[j]+Chemo[patient==LETTERS[i],"sensitivity"][j]/19
  }
}
lines(seq(14),moyenne_sens,type="l",col="red",lwd=5)
text(0.25,0.65,labels="Mean",col="red")

#J-S mean #acts weird, I surely made a mistake
#first, we need its variance
variance_sens=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    variance_sens[j]=variance_sens[j]+(Chemo[patient==LETTERS[i],"sensitivity"][j]-moyenne_sens[j])^2/19
  }
}
#then we compute it
moyenne_JS_sens=rep(0,14)
for (i in 1:19){
  for (j in Chemo[patient==LETTERS[i],"month"]){
    #I directly sum the p_i_JS and take their mean
    moyenne_JS_sens[j]=moyenne_JS_sens[j]+ (moyenne_sens[j] + (1-16*(moyenne_sens[j]*(1-moyenne_sens[j])/19)/variance_sens[j])*(Chemo[patient==LETTERS[i],"sensitivity"][j]-moyenne_sens[j]))/19
  }
}
lines(seq(14),moyenne_JS_sens,type="l",col="blue",lwd=5)
text(1.5,0.25,labels="James Stein Mean",col="blue")



# PLOTS -------------------------------------------------------------------
#I transformed data, I have to re-set it to its original values
Chemo = read.table(file="Chemotherapy-version1.csv",header=TRUE,sep=",")
Chemo[,"line"] <- as.factor(Chemo[,"line"]);


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


