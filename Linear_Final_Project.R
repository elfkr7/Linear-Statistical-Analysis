##load the data
data_tumor_raw=read.table("tumor_size.txt")
data_tumor=data_tumor_raw
colnames(data_tumor)=c("y","x1","x2","x3","x4","x5")
colnames(data_tumor_raw)=c("y","x1","x2","x3","x4","x5")
head(data_tumor)


##Q1
library(leaps)
library(MASS)
library(olsrr)
library(MPV)

##Strategy for Vaiable Selection and Model Building

##Step 1: Fit the largest model possible to the data and perform thorough analysis

data_tumor$x5=as.factor(data_tumor$x5) ##since x5 is an indicator variable
model=lm(y~.,data=data_tumor)
summary(model)
anova(model)##MSRES 187.

##Intercept and the regressors of x1,x3,x4,x52,x53 are statistically significant and adjusted R2 is 86.7%.
##Pvalue is almost 0 and this says regression is signifcant.

##Step 2: Residual Analysis

#Normal Probability Plot of Residuals
Rstudent_residual = rstudent(model)
sorted_Rstudent_residual=sort(Rstudent_residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_Rstudent_residual))
for(i in 1:length(sorted_Rstudent_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_Rstudent_residual)
}

plot(x=sorted_Rstudent_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_Rstudent_residual))
qqnorm(sorted_Rstudent_residual,xlim = c(-4,4),ylim = c(-4,4),main='Normal Q-Q Plot of the Residuals')
qqline(sorted_Rstudent_residual)

## Deviations exist in the upper and the lower tails. There are problems with the normality of the residual.

#Plot of Residuals against the Fitted Values

yhat = model$fit 
plot(x=yhat,y=Rstudent_residual)
abline(0,0)

##First there is nonlinear pattern in the graph and resiudals for some observations are out of the band of -2,2.
##These observations having studentized residuals are greater than 2 and smaller than -2 are potential outliers
##that need further analyses. Also, nonlinear pattern in the graph indicates some other regressors or transformation
##on the regressors and/or response variable may be needed to deal with this nonlinear pattern in the graph.

library(car)
avPlots(model = model,marginal.scale = T,id=T)

##In Added-Variable plots, there is a problem with linearity for x1 and x2 variables since all the points are located at almost the same point
##This may indicate a sort of collinearity. X4 also needs further analyses.

##Pair plots of x1 and x2

plot(data_tumor$x1,data_tumor$x2)

##Almost a perfect linearity between x1 and x2.

plot(data_tumor$x1,data_tumor$x4)
##There is a linear pattern at some level between x1 and x4. For seriousness of multicollinearity, further tests are needed.

plot(data_tumor$x2,data_tumor$x4)
##There is a linear pattern at some level between x2 and x4. For seriousness of multicollinearity, further tests are needed.

##Transformation

##Since we have problem with the normally assumption of the residual, we should transform on the response
##variable to handle this problem with the normally we see in normality plot of the 
##residual and nonlinear pattern in the graph of residual vs fitted values. 

bc = boxcox(object = model)
lambda = bc$x[which.max(bc$y)]
lambda

##apply lambda to the response variable

data_tumor$y=(data_tumor$y)^lambda
model2=lm(y~.,data=data_tumor)
summary(model2) ##Adjusted R square is 87.3%.
anova(model2)##MSRES 16.3.
#Normal Probability Plot of Residuals
Rstudent_residual = rstudent(model2)
sorted_Rstudent_residual=sort(Rstudent_residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_Rstudent_residual))
for(i in 1:length(sorted_Rstudent_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_Rstudent_residual)
}

plot(x=sorted_Rstudent_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_Rstudent_residual))

## It seems like deviations in the tails are healed. Except for some points located around tails and having deviations, 
##the normality seems okay. But for those points we need to do outlier detection analyses.

#Plot of Residuals against the Fitted Values

yhat = model2$fit 
plot(x=yhat,y=Rstudent_residual)
abline(0,0)

##We see potential outlier points having R student residuals more than 2 and smaller than -2. Also, we can observe nonlinear pattern existed in the previous 
##graph has been handled.Most of the points are located around r student=0 line without any pattern between the band of -2 and 2.

##To check relations of response variables and regressors and transformation on regressor x4.

plot(data_tumor$x1,data_tumor$y)
plot(data_tumor$x2,data_tumor$y)
plot(data_tumor$x3,data_tumor$y)
plot(data_tumor$x4,data_tumor$y)## non linear relationship between x4 and y. Taking square of x4 may help.

data_tumor$x4=data_tumor$x4^2
model3=lm(y~.,data=data_tumor)
summary(model3) ##Adjusted R square is 93.7%.
anova(model3)##MSRES 8.1.


##Outlier Analysis
Rstudent_residual = rstudent(model3)
sorted_Rstudent_residual=sort(Rstudent_residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_Rstudent_residual))
for(i in 1:length(sorted_Rstudent_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_Rstudent_residual)
}

plot(x=sorted_Rstudent_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_Rstudent_residual))

#Plot of Residuals against the Fitted Values

yhat = model3$fit 
plot(x=yhat,y=Rstudent_residual)
abline(0,0)
obs_potential_outliers=Rstudent_residual[which(abs(Rstudent_residual) > 2.25)]
obs_potential_outliers

##Based on the graph, the cut off value to distinguih outliers can be set as band of 2.25. Observations out of the band of abs(2.30)
##could be named as outlier. To get to know them better, move farward to leverage and influential points analyses.

##Leverage and influential points analysis
summary(influence.measures(model3))

##It looks like the measures of cook.d ,dfbetas,and hat are pretty normal but the measure that could catch the problem are diffits and covratios.
##Cutoff for the dffits is 2*sqrt(p/n) and based on the cutoff value, observations of 52,59,83,119,and 127 are potential influential points.
##All of them have also high residuals based on the outlier detection analysis.
##cutoff value for covratio is COVRATIO_i < 1 - 3*p / n or COVRATIO_i > 1 +3*p / n   ( we can trust this ratio since in our case n=130>3*p=7). Observations are out of the range are displayed in the list with star.
##OVerlaping observation for both measures(dffits and covraito) are 52,59,83,and 119, which are potential influential points. Alll of them have high residuals.
cutoff_dffits=2*sqrt(7/nrow(data_tumor))
cutoff_dffits
cutoff_covr_lower_boundary=1-3*7/130
cutoff_covr_lower_boundary
cutoff_covr_upper_boundary=1+3*7/130
cutoff_covr_upper_boundary


##Alternative Model Analysis without outliers
obs_num_outl=as.numeric(names(obs_potential_outliers))
obs_num_outl

data_tumor_without_outlier=data_tumor[-obs_num_outl,]

model4=lm(y~.,data_tumor_without_outlier)
summary(model4)##ADjusted R2 is 95.9% and coefficients of variables are pretty similar to ones in the model before outliers were removed.
anova(model4) ##MRES 4.6

#Normal Probability Plot of Residuals
Rstudent_residual = rstudent(model4)
sorted_Rstudent_residual=sort(Rstudent_residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_Rstudent_residual))
for(i in 1:length(sorted_Rstudent_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_Rstudent_residual)
}

plot(x=sorted_Rstudent_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_Rstudent_residual))

## Except for one point having Rstudent residual around4, the resiudals fit normality pretty well.

#Plot of Residuals against the Fitted Values

yhat = model4$fit 
plot(x=yhat,y=Rstudent_residual)
abline(0,0)

##Since p=7<30,  all possible regressions is feasible.
##The Best Regression Models

all_poss_reg=ols_step_all_possible(model3)
all_poss_reg
all_poss_reg[which.min(all_poss_reg$cp),]
all_poss_reg[which.max(all_poss_reg$adjr),]
PRESS(model3)
PRESS(lm(y~x1+x3+x4+x5,data=data_tumor))

##The best model which has the smallest MAllow's Cp and the highest adj R2 is y=B0+B1*x1+B2*x2+B3*x3+B4*x4+B5*x52+B6*x53+E
##The best model's Cp is 5, adj R2 is 93.7% and PRESS stat is 1143.369. 
##The second best model with regressors x1,x3,x4 and x5 has PRESS SCORE 1212.582. So, the best model has the smallest PRESS score. 

##stepwise algorithm

stepwise_algorithm=ols_step_both_p(model3,pent = 0.05,prem = 0.05)
stepwise_algorithm

##Stepwise algorithm uses F stat as threshold and F stat are calculated using MSRES statistics. Stepwise algorithm also
##chose the same best model with all regressors.

##Q2
cov(data_tumor_raw[-1])

##Based on the covariance matrix, we see that x1 and x2 has an important correlation.

model_with_interactions=lm(y~.^5,data=data_tumor)
summary(model_with_interactions)
ols_coll_diag(model_with_interactions)

##One condition index exceeds 1000 so there is at least one serious multicollinearity in the data. There are other 4 condition indices 
##greater than 100. For vifs, especially for x1 and x2 are super high. We can expect interaction term may have high vif
##but the smallest vif is 13.48 which is greater than 5 and 10. We can conclude the data has serious mulitcollinearity issue.

#Q3

model_ridge=lm.ridge(model_with_interactions,lambda=seq(0, 0.1, 0.0001))
plot(model_ridge)
MASS::select(model_ridge)

#Refit the model with the best lambda = 0.0242
model_ridge_best=lm.ridge(model_with_interactions,lambda=0.0113)
coef(model_ridge_best)

data_tumor_full = data.frame(model.matrix(~(.)^5,data_tumor[-1]))
data_tumor_full=data.frame(data_tumor$y,data_tumor_full)[-2] ##drop intercept

##find y_hat for each observation with rigde model
predicted = model_ridge_best$ym +
  scale(data_tumor_full[,-1], center = model_ridge_best$xm,
        scale = model_ridge_best$scales) %*%
  model_ridge_best$coef

##metrics for ridge regression
SSRES_ridge = sum((predicted - data_tumor_full$data_tumor.y)^2)
SSRES_ridge

MSRES_ridge = SSRES_ridge/(nrow(data_tumor_full)-ncol(data_tumor_full)-1)
MSRES_ridge

SST_ridge = sum((data_tumor$y - mean(data_tumor$y))^2)
SST_ridge

SSR_ridge=SST_ridge-SSRES_ridge
SSR_ridge

R2_ridge= (1-SSRES_ridge/SST_ridge)
R2_ridge

adjR2_ridge = 1 -  ((SSRES_ridge/(nrow(data_tumor_full)-ncol(data_tumor_full)-1))/(SST_ridge/(nrow(data_tumor_full)-1)))
adjR2_ridge

  MSE_ridge= sqrt(mean((predicted - data_tumor_full$data_tumor.y)^2))
MSE_ridge

##Q4

library(pls)
set.seed (1000)

model_pcr <- pcr(data_tumor.y~., data = data_tumor_full, scale = TRUE)
summary(model_pcr)

for (i in 1:47){
predicted_pcr=predict(model_pcr,data_tumor_full[-1],ncomp = i)

##metrics for PCR regression
SSRES_pcr = sum((predicted_pcr - data_tumor_full$data_tumor.y)^2)
print(SSRES_pcr)

MSRES_pcr = SSRES_pcr/(nrow(data_tumor_full)-i)
MSRES_pcr

SST_pcr = sum((data_tumor_full$data_tumor.y - mean(data_tumor_full$data_tumor.y))^2)
SST_pcr

SSR_pcr=SST_pcr-SSRES_pcr
SSR_pcr

R2_pcr= (1-SSRES_pcr/SST_pcr)
R2_pcr

adjR2_pcr = 1 -  ((SSRES_pcr/(nrow(data_tumor_full)-i))/(SST_pcr/(nrow(data_tumor_full)-1)))
print(adjR2_pcr)

MSE_pcr= sqrt(mean((predicted_pcr - data_tumor_full$data_tumor.y)^2))
MSE_pcr

}

##select ncomp=7 based on x variance approach.

predicted_pcr=predict(model_pcr,data_tumor_full[-1],ncomp = 7)

##metrics for PCR regression
SSRES_pcr = sum((predicted_pcr - data_tumor_full$data_tumor.y)^2)
SSRES_pcr

MSRES_pcr = SSRES_pcr/(nrow(data_tumor_full)-7)
MSRES_pcr

SST_pcr = sum((data_tumor_full$data_tumor.y - mean(data_tumor_full$data_tumor.y))^2)
SST_pcr

SSR_pcr=SST_pcr-SSRES_pcr
SSR_pcr

R2_pcr= (1-SSRES_pcr/SST_pcr)
R2_pcr

adjR2_pcr = 1 -  ((SSRES_pcr/(nrow(data_tumor_full)-7))/(SST_pcr/(nrow(data_tumor_full)-1)))
adjR2_pcr

MSE_pcr= sqrt(mean((predicted_pcr - data_tumor_full$data_tumor.y)^2))
MSE_pcr
##Q5

##pcr residual plot analysis

residual = (data_tumor_full$data_tumor.y-predicted_pcr)
sorted_residual=sort(residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_residual))
for(i in 1:length(sorted_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_residual)
}

plot(x=sorted_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_residual))


#Plot of Residuals against the Fitted Values

plot(x=predicted_pcr,y=residual)
abline(0,0)

##ridge residual plot analysis

residual = (data_tumor_full$data_tumor.y-predicted)
sorted_residual=sort(residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_residual))
for(i in 1:length(sorted_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_residual)
}

plot(x=sorted_residual,y=p_i[1,],xlim =c(-1,1),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_residual))


#Plot of Residuals against the Fitted Values

plot(x=predicted,y=residual)
abline(0,0)

##best model in part a residual plot analysis
Rstudent_residual = rstudent(model3)
sorted_Rstudent_residual=sort(Rstudent_residual)

p_i=matrix(NA,nrow=1,ncol=length(sorted_Rstudent_residual))
for(i in 1:length(sorted_Rstudent_residual)){
  p_i[[i]]=(i-0.5)/length(sorted_Rstudent_residual)
}

plot(x=sorted_Rstudent_residual,y=p_i[1,],xlim =c(-4,4),ylim=c(0,1),main='Normal Probability Plot of the Residuals')
abline(lm(p_i[1,]~sorted_Rstudent_residual))

#Plot of Residuals against the Fitted Values

yhat = model3$fit 
plot(x=yhat,y=Rstudent_residual)
abline(0,0)

##Q6

set.seed(1923)
train_=sample(c(1:dim(data_tumor_full)[1]), dim(data_tumor_full)[1]*0.7)

estimation=data_tumor_full[train_, ]
prediction=data_tumor_full[-train_, ]

model_pcr_train = pcr(data_tumor.y~., data = estimation, scale = TRUE)
summary(model_pcr_train)
#95% of variability is explained by the first 6 component but the difference between 7th and 6th is so small.
##To interpret the result better, ncom=7 will be kept same.

test_predicted=predict(model_pcr_train,prediction[-1],ncomp = 7)

##metrics for PCR regression

MSE_pcr= sqrt(mean((test_predicted - prediction$data_tumor.y)^2))
MSE_pcr

a=as.data.frame(model_pcr_train$coefficients)
view(a)
round(prcomp(estimation[-1])$rotation,3)

pca1=prcomp(estimation[-1],scale. = T)
scores_pca1=pca1$x[,1:7]
reg_pca=lm(estimation$data_tumor.y~scores_pca1)
summary(reg_pca)
