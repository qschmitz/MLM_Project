###############################################################################
#-------------------------Mixed Linear Models Project-------------------------#
#--------------------------------------by-------------------------------------#
#--------------------------------------------<You can put your name here>-----#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------<Or here>------------------------#
#--------------------------------------------Hervégil Voegeli-----------------#
###############################################################################
# Research questions:
#   1)if the in-vitro sensitivity score can be used to predict the treatment effect on the patient at the hospital
#   2)if the in-vitro sensitivity score decreases with time and is stronger for line 1 than for lines 2 and 3+ 


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
chemo = Chemo#I'll need this as I will transform somehow Chemo
head(Chemo);#View(Chemo);
#data size :  137x5
#We have 19 patients from A to S
#month are numeroted only for treatments line and reset to 0 after each treatment
#we do not have the same number of treatment per patient nor same number of month per treatment:
# patient D: only 1 line of 1 month (change of patient after that->death)
# patient E, line 2: lasted only 2 month (then start line 3->chimio side effects (anemia, vomit, etc.) intolerable)
#data is therefore unbalanced (->anova type III != type II)

attach(Chemo)


# Studying variables -------------------------------------------------------

#check data type
sapply(Chemo, class)#line should be factor
Chemo[,"line"] <- as.factor(Chemo[,"line"]);
chemo[,"line"] <- as.factor(chemo[,"line"]);
#also month?->ofc not
# Chemo[,"month"] <- as.factor(Chemo[,"month"])


# Overview tumour/sensitivity -------------------------------------------------------------------

#I want to see the evolution of the tumour through time:
#month should be re-indexed and restart to 0 for each new patient
for (i in 1:19){
  chemo[patient==LETTERS[i],"month"]=seq(1,length(chemo[patient==LETTERS[i],"month"]))
}
#we can plot the tumour growth
plot(chemo[patient==LETTERS[1],"month"],chemo[patient==LETTERS[1],"tumour"],
     type="l",main="Tumour growth through time",
     xlab = "Time in months",ylab = "Tumour evolution per month of treatment",
     xlim=c(1,14),ylim=c(-3,1))
for (i in 2:19){
  lines(chemo[patient==LETTERS[i],"month"],chemo[patient==LETTERS[i],"tumour"],type="l",col=rand_color(1))
}
#and display the mean for the sake of clarity
moyenne=rep(0,14)
for (i in 1:19){
  for (j in chemo[patient==LETTERS[i],"month"]){
    moyenne[j]=moyenne[j]+chemo[patient==LETTERS[i],"tumour"][j]/19
  }
}
lines(seq(14),moyenne,type="l",col="red",lwd=5)

#let's add the 95% CI for this line
#unknown variance-> Student
#unbalanced data-> df changes over time
df=rep(0,14)
for (i in 1:19){
  for (j in chemo[patient==LETTERS[i],"month"]){
    df[j]=df[j]+1
  }
}
Up_moyenne=moyenne+qt(0.975,df-1)*sqrt(variance)/sqrt(19)
Lo_moyenne=moyenne+qt(0.025,df-1)*sqrt(variance)/sqrt(19)
lines(seq(14),Up_moyenne,type="l",col="coral",lwd=2)
lines(seq(14),Lo_moyenne,type="l",col="coral",lwd=2)
legend(12,-1,legend=c("Mean","95% CI"),col=c("red","coral"),fill=c("red","coral"))
# Conclusions based on graph:
# Drug treatments have a smaller effect through time
# we cannot conclude the effect of a change of drug from this graph as it appears at different point of time depending of the patient


# #J-S mean
# #first, we need its variance
# variance=rep(0,14)
# for (i in 1:19){
#   for (j in chemo[patient==LETTERS[i],"month"]){
#     variance[j]=variance[j]+(chemo[patient==LETTERS[i],"tumour"][j]-moyenne[j])^2/19
#   }
# }
# #then we compute it
# moyenne_JS=rep(0,14)
# for (i in 1:19){
#   for (j in chemo[patient==LETTERS[i],"month"]){
#     #I directly sum the p_i_JS and take their mean
#     moyenne_JS[j]=moyenne_JS[j]+ (moyenne[j] + (1-16*(moyenne[j]*(1-moyenne[j])/19)/variance[j])*(chemo[patient==LETTERS[i],"tumour"][j]-moyenne[j]))/19
#   }
# }
# lines(seq(14),moyenne_JS,type="l",col="blue",lwd=5)
# text(5,-1.5,labels="James Stein Mean",col="blue")


#Now the same with the evolution of sensitivity
plot(chemo[patient==LETTERS[1],"month"],chemo[patient==LETTERS[1],"sensitivity"],
     type="l",main="Sensitivity to the drug combination of interest during treatment",
     xlab = "Time in months",ylab = "Sensitivity",
     xlim=c(0,14),ylim =c(-0.1,1))
for (i in 2:19){
  lines(chemo[patient==LETTERS[i],"month"],chemo[patient==LETTERS[i],"sensitivity"],type="l",col=rand_color(1))
}
#and display the mean for the sake of clarity
moyenne_sens=rep(0,14)
for (i in 1:19){
  for (j in chemo[patient==LETTERS[i],"month"]){
    moyenne_sens[j]=moyenne_sens[j]+chemo[patient==LETTERS[i],"sensitivity"][j]/19
  }
}
lines(seq(14),moyenne_sens,type="l",col="red",lwd=5)

#95% CI
#df has been calculated in the previous graph
Up_moyenne_sens=moyenne_sens+qt(0.975,df-1)*sqrt(variance_sens)/sqrt(19)
Lo_moyenne_sens=moyenne_sens+qt(0.025,df-1)*sqrt(variance_sens)/sqrt(19)
lines(seq(14),Up_moyenne_sens,type="l",col="coral",lwd=2)
lines(seq(14),Lo_moyenne_sens,type="l",col="coral",lwd=2)
legend(12,1,legend=c("Mean","95% CI"),col=c("red","coral"),fill=c("red","coral"))

#Conclusions:
#Overall sensitivity is significantly reducing through time.

# #J-S mean #acts weird, I surely made a mistake
# #first, we need its variance
# variance_sens=rep(0,14)
# for (i in 1:19){
#   for (j in chemo[patient==LETTERS[i],"month"]){
#     variance_sens[j]=variance_sens[j]+(chemo[patient==LETTERS[i],"sensitivity"][j]-moyenne_sens[j])^2/19
#   }
# }
# #then we compute it
# moyenne_JS_sens=rep(0,14)
# for (i in 1:19){
#   for (j in chemo[patient==LETTERS[i],"month"]){
#     #I directly sum the p_i_JS and take their mean
#     moyenne_JS_sens[j]=moyenne_JS_sens[j]+ (moyenne_sens[j] + (1-16*(moyenne_sens[j]*(1-moyenne_sens[j])/19)/variance_sens[j])*(chemo[patient==LETTERS[i],"sensitivity"][j]-moyenne_sens[j]))/19
#   }
# }
# lines(seq(14),moyenne_JS_sens,type="l",col="blue",lwd=5)
# text(1.5,0.25,labels="James Stein Mean",col="blue")

# PLOTS -------------------------------------------------------------------

#study data with plots
  #Boxplots: use it for heteroskedasticity, distribution
boxplot(tumour~patient)#is 19 enough to diagnose patient heteroskedasticity?

#this next boxplot is kinda biased since some "month 3" appears befort some "month 1"
boxplot(tumour~month)#left-skewed distribution per month, the third month shows the smallest decline in tumour regression
boxplot(tumour~month,data = chemo,main="Tumour evolution per month of treatment")#that's more accurate!

boxplot(tumour~line)#which is consistent with less effect of lines through time, less variance
boxplot(sensitivity~patient)#clear heteroskedasticity of sensitivity distribution through patients


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



