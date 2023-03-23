#################################################
#  Modeling Malaria  Data  using NS-MBG         #
#################################################
rm(list=ls())
# Set working directory
setwd("C:\\Users\\Bedilu\\Desktop\\Bedilu\\ResearchProjects\\WorkingPapers\\PI\\Malaria\\Data-Code\\Malaria")

# Loading libraries and our function to fit the non-stationary model
library(geoR)
library(PrevMap)
library(ncf)
source("NSgeo_MLE.R") # Calling our non-stationary R-module

#Load the data
d=read.csv("d.csv")
#Adding emperical logit transformation of malaria prevalence in the data frame
d$logit <- log((d$positive + 0.5)/(d$examined-d$positive + 0.5))

# Explore how residuals from the classical LM are location dependent
# Which implies that spatial models are needed.
# futher implies  the outcome variable is influenced by location specific covariates 
ols=lm(d$logit ~alt+temp+hum+prec+dist_aqua, data=d)
plot(d$longitude, d$latitude, col=c("blue", "red")[sign(resid(ols))/2+1.5], pch=19,
     cex=abs(resid(ols))/max(resid(ols))*2, xlab="geographical X-coordinates", ylab="geographical Y-coordinates")

# Split the data into different groups to see how covaraites affect the variability
d$prec_quartile<-cut(d$prec,quantile(d$prec),include.lowest=TRUE,labels=FALSE)
d$water_quartile<-cut(d$dist_aqua,quantile(d$dist_aqua),include.lowest=TRUE,labels=FALSE)
boxplot(logit ~prec_quartile, data= d, xlab = "Precipitation quartile",
        ylab = "Emperical logit transformation of malaria prevalence", main = "")
boxplot(logit ~water_quartile, data= d, xlab = "Distance from water quartile",
        ylab = "Emperical logit transformation of malaria prevalence", main = "")

############################################
# Variogram and Correlogram Exploration    #
############################################
# a) variogam
# Case 1: Based on geographical distance
geodata.vg=as.geodata(d,coords=1:2, data.col=11)
dist.vg=variog(geodata.vg);plot(dist.vg)
varfit.dist=variofit(dist.vg, ini.cov.pars = c(3,5), cov.model ="matern")
summary(varfit.dist)
# Case 2: Based on covariates (Temperature, precipitation, distance from water)
# after exploration prec semivariogram resembles ecludian distance variogram
d$PREC=d$prec
geodata.prec=as.geodata(d,coords=c(8,12), data.col =11)
prec.vg=variog(geodata.prec);plot(prec.vg)

# Variogram and correlogram Plot
par(mfrow = c(2,2))
plot(dist.vg,xlab="Euclidean distance (u)", main="(a)")
plot(prec.vg,xlab="Precipitation differennce", main="(b)")
plot(dist.cor,xlab="Euclidean distance (u)", main="(c)")
plot(prec.cor,xlab="Precipitation differennce", main="(d)")


##############################################
#      Geostatistical Modeling               #
##############################################
# a) Standard MBG
ID.coords0 <- create.ID.coords(d,~longitude + latitude)
# Excluding precipitation from the mean structure
MBG.1 <- linear.model.MLE(formula=logit~temp+alt+hum+dist_aqua, coords = ~longitude + latitude, data = d,
                          start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0, kappa = 0.5, method="nlminb")
summary(MBG.1, log.cov.pars = FALSE)
# 95% CI for parameter estimates
beta.hat.geo1<- summary(MBG.1)$coefficients[,1:2]
cbind(beta.hat.geo1,
      beta.hat.geo1[,1]-qnorm(0.975)*beta.hat.geo1[,2],
      beta.hat.geo1[,1]+qnorm(0.975)*beta.hat.geo1[,2])

# Including precipitation from the mean structure
MBG.2 <- linear.model.MLE(formula=logit~temp+alt+prec+hum+dist_aqua, coords = ~longitude + latitude, data = d,
                              start.cov.pars = c(2, 1),ID.coords=ID.coords0,fixed.rel.nugget=0, kappa = 0.5, method="nlminb")
summary(MBG.2, log.cov.pars = FALSE)

# 95% CI for parameter estimates
beta.hat.geo2 <- summary(MBG.2)$coefficients[,1:2]
cbind(beta.hat.geo2,
      beta.hat.geo2[,1]-qnorm(0.975)*beta.hat.geo2[,2],
      beta.hat.geo2[,1]+qnorm(0.975)*beta.hat.geo2[,2])

# Model validation for linear geostatistical model using Monte Carlo methods based on the variogram
variog.diag <-variog.diagnostic.lm (MBG.2, n.sim = 1000,
                                    uvec = NULL,plot.results = TRUE,
                                    range.fact = 1, which.test = "both",
                                    param.uncertainty = FALSE)

variog.diag <-variog.diagnostic.lm (MBG.2, n.sim = 1000,
                                    uvec = NULL,plot.results = TRUE,
                                    range.fact = 1, which.test = "variogram",
                                    param.uncertainty = FALSE)

# Model Validation using cross- validation 
#Creating training data as 80% of the dataset
library(caret)
set.seed(123)
rs.d <-createDataPartition(d$logit, p = 0.8, list = FALSE)
#Generating training dataset (td) from the random sample(rs)
td.d  <- d[rs.d, ]
# Generating testing dataset from rows which are not included in random_sample
testing_dataset <- d[-rs.d, ]
# Fit the model using the training dataset
ID.coords_training <- create.ID.coords(td.d,~longitude + latitude)
fit.training <- linear.model.MLE(formula=logit~temp+alt+prec+hum+dist_aqua, coords = ~longitude + latitude, data = td.d,
                                 start.cov.pars = c(2, 1),ID.coords=ID.coords_training ,fixed.rel.nugget=0, kappa = 0.5, method="nlminb")
summary(fit.training, log.cov.pars = FALSE)
# Predict test dataset
test.Pred<-spatial.pred.linear.MLE(fit.training, grid.pred=testing_dataset[,c(1,2)],
                                   predictors = testing_dataset[,c(7,6,8,9,10)],
                                   scale.predictions ="logit",
                                   standard.errors = TRUE,
                                   thresholds = -2,
                                   scale.thresholds = "logit")
# Error
pred=test.Pred$logit[[1]]
e= pred-testing_dataset$logit;
hist(e)
Error=cbind(pred,testing_dataset$logit, pred-testing_dataset$logit )
plot(pred,testing_dataset$logit)
plot(e,pred)
cat(paste ("Mean Square error: ", round(sum(e*e)/length(e),3), "\n"))

#b) Nonstationary MBG (NS-MBG)

#Case 1: Including precipitation  only in the covariance structure
ID.coords_prec <- create.ID.coords(d,~longitude + latitude+prec)
NS_MBG.1 <- NSgeo.MLE(logit ~temp+alt+hum+dist_aqua, 
                         coords = ~longitude + latitude,
                         rand.effect.domain = ~ prec,
                         start.cov.pars = c(2,1),fixed.rel.nugget=0,
                         data=d,ID.coords = ID.coords_prec,method="nlminb",
                         messages = TRUE,return_se=TRUE)
NS_MBG.1 
#To get parameter estimates with their respective 95% CI
beta1=NS_MBG.1$par; se.beta1=NS_MBG.1$se
lcl1=beta1-1.96*se.beta1; ucl1=beta1+1.96*se.beta1
CI1=cbind(beta1,se.beta1,lcl1,ucl1);CI1

#Case 2: Including precipitation  both in the mean and covariance structure
NS_MBG.2 <- NSgeo.MLE(logit ~temp+alt+prec+hum+dist_aqua, 
                         coords = ~longitude + latitude,
                         rand.effect.domain = ~prec,
                         start.cov.pars = c(2,1),fixed.rel.nugget=0,
                         data=d,ID.coords = ID.coords_prec,method="nlminb",
                         messages = TRUE,return_se=TRUE)
NS_MBG.2
# Parameter estimates with 95%CI
beta2=NS_MBG.2$par; se.beta2=NS_MBG.2$se
lcl2=beta2-1.96*se.beta2; ucl2=beta2+1.96*se.beta2
CI2=cbind(beta2,se.beta2,lcl2,ucl2);CI2

# Model Comparison using AIC
cat(paste ("Stationary MBG.1 AIC: ",  round(-2*MBG.1$log.lik + 2*length(MBG.1$par),3), "\n"))
cat(paste ("Stationary MBG.2 AIC: ",  round(-2*MBG.2$log.lik + 2*length(MBG.2$par),3), "\n"))
cat(paste ("Nonstationary Model NS_MBG.1 : ", round(2*length(NS_MBG.1$par)+2*NS_MBG.1$objective,3), "\n"))
cat(paste ("Nonstationary Model NS_MBG.2: ", round(2*length(NS_MBG.2$par)+2*NS_MBG.2$objective,3), "\n"))

#######################################################
#      Malaria Prediction  to unsampled location      #
#######################################################

dp=read.csv("dp.csv")
dim(dp); names(dp);head(dp)
dim(unique(dp[, c("longitude", "latitude")])) # To know unique no of locations
# Getting Mozambique boundary
library(prevR)
mz.bdry=create.boundary(countries="Mozambique", multiple = F, proj = "+proj=longlat +datum=WGS84" )

# Case 1: Using the stationary MBG
stand.Pred<-spatial.pred.linear.MLE(MBG.2, grid.pred=dp[,c(1,2)],
                                    predictors = dp[,c(4,3,5,6,7)],
                                    scale.predictions ="logit",
                                    standard.errors = TRUE,
                                    thresholds = -2,
                                    scale.thresholds = "logit")

stand.Pred$exceedance.prob <- 1-stand.Pred$exceedance.prob
plot(stand.Pred,summary="exceedance.prob",zlim=c(0,1),main="Malaria risk")
lines(mz.bdry)


# Case 2: Using NS-MBG approach
NS_MBG.Pred<-spatial.pred.linear.MLE(NS_MBG.2, grid.pred=dp[,c(1,2)],
                          predictors = dp[,c(4,3,5,6,7)],
                          scale.predictions ="logit",
                          standard.errors = TRUE,
                          thresholds = -2,
                          scale.thresholds = "logit")
NS_MBG.Pred$exceedance.prob <- 1-NS_MBG.Pred$exceedance.prob
plot(NS_MBG.Pred,summary="exceedance.prob",zlim=c(0,1),main="Malaria risk")
lines(mz.bdry)


