#model diagnosis 
# generates traceplots and density plots calculates as well as summary statistics 

#Packages####
#use this of you haven't installed many of these packages
#list.of.packages <- c("vegan","plotrix","MASS","data.table", "gdata","lattice","plyr","dplyr", "lme4", "arm", "gridExtra", "ggplot2", "eeptools", "taxize", "lmtest", "Biostrings", "tidyr", "googlesheets", "wesanderson", "drc", "zoo","unmarked","jagsUI", "data.table","vegan","digest")
#new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
#if(length(new.packages)) install.packages(new.packages)

libs=c("vegan","plotrix","MASS","data.table", "gdata","lattice","plyr","dplyr", "lme4", "arm", "gridExtra", "ggplot2", "eeptools", "taxize", "lmtest", "Biostrings", "tidyr", "googlesheets", "wesanderson", "drc", "zoo","unmarked","jagsUI", "data.table","vegan","digest","tidyverse")
lapply(libs, require, character.only = TRUE)  #"gsubfn"


#Functions####
# This function is later used for printing simulation results on the screen
printsummres<-function(thetahat,thename="estimated parameter"){
  cat(paste0("\n",thename,": mean = ",round(mean(thetahat),2),
             ", 2.5% = ",round(quantile(thetahat,prob=0.025),2),
             ", 97.5% = ",round(quantile(thetahat,prob=0.975),2)))
}

ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
cumProb=function(x) 1-prod(1-x)  #define cumulative prob function
sitesFun<-function(x, Nreplicates=3){rep(1:x, each=Nreplicates)}

#functions

printsummres<-function(thetahat,thename="estimated parameter"){
  cat(paste0("\n",thename,": mean = ",round(mean(thetahat),2),
             ", 2.5% = ",round(quantile(thetahat,prob=0.025),2),
             ", 97.5% = ",round(quantile(thetahat,prob=0.975),2)))
}

ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
cumProb=function(x) 1-prod(1-x)  #define cumulative prob function
sitesFun<-function(x, Nreplicates=3){rep(1:x, each=Nreplicates)}

# the model ####
p10_max=0.1                             
##############################################
##############################################
#The Actual Model
##############################################
##############################################
##############################################

####Bayesian version
### Run first this part once to create the file with the JAGS model
sink("RoyleLink_prior.txt")
cat("model {
                                    # Priors
                                    psi ~ dunif(0,1)
                                    p11 ~ dunif(0.01,1)
                                    p10 ~ dunif(0.001, p10_max)
                                    
                                    # Likelihood 
                                    for (i in 1:S){
                                    z[i] ~ dbern(psi)
                                    p[i] <- z[i]*p11 + (1-z[i])*p10
                                    for (j in 1:K){
                                    Y[i,j] ~ dbern(p[i])
                                    }
                                    }
                                    } ",fill=TRUE)
sink()


model_Bayesian <- function(datalist,nOTUs=length(datalist), S=46, K=3, doprint=FALSE,p10_max=0.1,
                           ni=5000,nt=2,nc=10,nb=10,myparallel=TRUE) {   
  psihat<-p11hat<-p10hat<-rep(nOTUs)
  modelSummaries<-list()
  hh<-datalist[[1]]
  # fit the model    
  jags.inits <-function()(list(psi=runif(1,0.05,0.95),p11=runif(1, 0.01,1),p10=runif(1,0.001,p10_max)))
  jags.data  <-list(Y=hh,S=S[[1]],K=K,p10_max=p10_max)
  jags.params<-c("psi","p11","p10")
  model<-jags(data = jags.data, inits = jags.inits, parameters.to.save= jags.params, 
              model.file= "RoyleLink_prior.txt", n.thin= nt, n.chains= nc, 
              n.iter= ni, n.burnin = nb, parallel=myparallel)
  #jpeg(paste0(format(Sys.time(), "%H_%M_%S"),"_ModelParamsPlotDEMO.jpg"))
  plot(model)
  dev.off()
  plot(model)
  
  psihat[1] <- model$summary["psi","50%"]
  p11hat[1] <- model$summary["p11","50%"]
  p10hat[1] <- model$summary["p10","50%"]    
  modelSummaries[[1]]<-model$summary
  
  
  if (doprint){
    printsummres(psihat,thename="estimated psi")
    printsummres(p11hat,thename="estimated p11")
    printsummres(p10hat,thename="estimated p10")
  }
  #saveRDS(modelSummaries, paste0(format(Sys.time(), "%H_%M_%S"),"_ModelSummaries_DEMO.rds"))
  BayesResults<-list(psihat=psihat,p11hat=p11hat,p10hat=p10hat,modelSummaries=modelSummaries)
  return(BayesResults)
  return(plot(model))
}
#############################
ASVlist <- readRDS(file="Data/2022_10_31/derived_data/scratch/ASVlist.rds")

#let's compare a few different detetion histories
#a few examples of histories with differing numbers of PCR detections (n = 168)
p1 <- as.data.frame(ASVlist[1]) # found in 243 survey reps
sum(p1)
p2 <- as.data.frame(ASVlist[5]) #found in 193 of survey reps
sum(p2)
p3 <- as.data.frame(ASVlist[10]) #found in 122 survey reps
sum(p3)
p4 <- as.data.frame(ASVlist[20]) #found in 74 survey reps
sum(p4)
p5 <- as.data.frame(ASVlist[50]) #found in 31 PCR reps
sum(p5)
p6 <- as.data.frame(ASVlist[100]) #found in 9 PCR reps
sum(p6)
p7 <- as.data.frame(ASVlist[150]) #found in 2 PCR reps
sum(p7)
p8 <- as.data.frame(ASVlist[200]) #found in 2 PCR reps
sum(p8)

#for ASVlist[1] ####
d1 <- as.data.frame(ASVlist[1])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 243 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")


#for ASVlist[5] ####
d1 <- as.data.frame(ASVlist[5])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 193 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[10] ####
d1 <- as.data.frame(ASVlist[10])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 122 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[20] ####
d1 <- as.data.frame(ASVlist[20])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 74 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[50] ####
d1 <- as.data.frame(ASVlist[50])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 31 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[100] ####
d1 <- as.data.frame(ASVlist[100])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 9 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[150] ####
d1 <- as.data.frame(ASVlist[150])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 2 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

#for ASVlist[200] ####
d1 <- as.data.frame(ASVlist[200])
datalist=list(d1)
modelsOut<-model_Bayesian(datalist)
modelsOut$modelSummaries
ProbOcc=function(x, psi=psi, p11=p11, p10=p10, K=K){
  (psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
}
K=3
psi=modelsOut$modelSummaries[[1]][1]
p11=modelsOut$modelSummaries[[1]][2]
p10=modelsOut$modelSummaries[[1]][3]
psi_975<-modelsOut$modelSummaries[[1]][25]
psi_025<-modelsOut$modelSummaries[[1]][9]  
p11_975<-modelsOut$modelSummaries[[1]][26]
p11_025<-modelsOut$modelSummaries[[1]][10]  
p01_975<-modelsOut$modelSummaries[[1]][27]
p01_025<-modelsOut$modelSummaries[[1]][11]  
plot(x=0:5, y=ProbOcc(x=c(0:5), psi=psi, p11=p11, p10=p10, K=K), pch=19,
     xlab="N detections out of 5 replicates", ylab="Probability of Occupancy", 
     type="b", main="Probability of Occupancy with constructed 95% Credibility Interval - 1 detections",
     ylim=range(1e-9:1))
#low CI has low psi, low p11, high p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_025, p11=p11_025, p10=p01_975, K=K), col="grey75", type="b")
#high CI has high psi, high p11, low p10
points(x=0:5, y=ProbOcc(x=c(0:5), psi=psi_975, p11=p11_975, p10=p01_025, K=K), col="grey75", type="b")

