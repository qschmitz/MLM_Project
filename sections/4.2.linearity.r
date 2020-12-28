# 4.2 Exploratory data analysis ----------------------------------------------

# • check if the model’s assumptions hold:
#   – linearity of the relation between the covariates and the response
# -> only for integer or continuous covariates, by definition useless for factors
plot(tumour~sensitivity, data=Chemo) # linearity seems ok, no quadratic (or else) pattern emerges
plot(tumour~month, data=chemo) #used indexed months, linearity also seems ok