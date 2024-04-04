library(data.table)
library(hal9001)
source("Data_generation-Masters_thesis")
set.seed(210399)

#HAL-TMLE estimator
HAL_TMLE <- function(d){
  
  d <- as.data.frame(d)
  
  #Data for estimating the propensity score
  treatment_data <- d[, !names(d) %in% "Y"]
  
  #HAL fit of the conditional mean
  Q_fit <- fit_hal(X = as.matrix(d[, !names(d) %in% "Y"]),
                   Y = d$Y,
                   family = "binomial",
                   lambda = 0.001,
                   fit_control = list(cv_select = FALSE,
                                      lambda.min.ratio = 1e-04,
                                      prediction_bounds = "default"))
  
  #HAL fit of the propensity score
  G_fit <- fit_hal(X = as.matrix(treatment_data[, !(names(treatment_data)) %in% "A"]),
                   Y = treatment_data$A,
                   family = "binomial",
                   lambda = 0.001,
                   fit_control = list(cv_select = FALSE,
                                      lambda.min.ratio = 1e-04,
                                      prediction_bounds = "default"))
  
  d <- data.table(d)
  
  #HAL estimate of the conditional mean under A=1
  d[, pred.Q0:=predict(Q_fit, type = "response",new_data=copy(d)[, A:=1])]
  
  #HAL estimate of the intervention specific mean
  HAL.est <- d[,mean(pred.Q0)]
  
  #HAL estimate of the conditional mean
  d[, pred.Q:=predict(Q_fit, type = "response", new_data=d)]
  
  #HAL estimate of the propensity score
  d[, pred.G:=predict(G_fit, type = "response" ,new_data=d)]
  
  #Clever covariate
  d[, H:=A/pred.G]
  
  #In case we get estimates that are 1 or 0
  d$pred.Q[d$pred.Q==0] <- 0.00000001
  d$pred.Q[d$pred.Q==1] <- 1-0.00000001
  
  #Targeting step
  fit_tmle <- glm(Y ~ offset(qlogis(pred.Q)) + H - 1, family=binomial, data=d) 
  
  eps.hat <- fit_tmle$coef
  
  d[, pred.tmle1:=plogis(qlogis(pred.Q0) + eps.hat/(pred.G))]
  d[, pred.tmle:=plogis(qlogis(pred.Q) + eps.hat*H)]
  
  #HAL-TMLE estimate of the intervention specific mean
  tmle.est <- d[, mean(pred.tmle1)]
  
  #Estimated standard deviation
  tmle.se <- d[, sqrt(mean((H*(Y-pred.tmle) + pred.tmle1 - tmle.est)^2)/nrow(d))]
  return(c(tmle.est=tmle.est, tmle.se=tmle.se, HAL.est=HAL.est))
}

#GLM-TMLE estimator
GLM_TMLE <- function(d){
  
  d <- as.data.frame(d)
  
  treatment_data <- d[, !names(d) %in% "Y"]
  
  form1 <- Y~.
  form2 <- A~.
  
  #GLM fit of the conditional mean
  Q_fit <- glm(form1, family = "binomial", d)
  
  #GLM fit of the propensity score
  G_fit <- glm(form2, family = "binomial", treatment_data)
  
  d <- data.table(d)
  
  #HAL estimate of the conditional mean under A=1
  d[, pred.Q0:=predict(Q_fit, type = "response",newdata=copy(d)[, A:=1])]
  
  #GLM estimate of the intervention specific mean
  glm.est <- d[, mean(pred.Q0)]
  
  #GLM estimate of the conditional mean
  d[, pred.Q:=predict(Q_fit, type = "response", newdata=d)]
  
  #GLm estimate of the propensity score
  d[, pred.G:=predict(G_fit, type = "response" ,newdata=d)]
  
  #GLM estimate of the conditional mean
  d[, H:=A/pred.G]
  
  #Targeting step
  fit_tmle <- glm(Y ~ offset(qlogis(pred.Q)) + H - 1, family=binomial, data=d) 
  
  eps.hat <- fit_tmle$coef
  
  d[, pred.tmle1:=plogis(qlogis(pred.Q0) + eps.hat/(pred.G))]
  d[, pred.tmle:=plogis(qlogis(pred.Q) + eps.hat*H)]
  
  #GLM-TMLE estimate of intervention specific mean
  tmle.est <- d[, mean(pred.tmle1)]
  
  #Estimated standard deviation
  tmle.se <- d[, sqrt(mean((H*(Y-pred.tmle) + pred.tmle1 - tmle.est)^2)/nrow(d))]
  return(c(tmle.est=tmle.est, tmle.se=tmle.se, glm.est=glm.est))
  
}



#Simulation study funtion
sim_study <- function(DGP, estimator, num_repetitions,n){
  
  sim <- matrix(ncol = 2, nrow = 2*num_repetitions)
  
  for (i in 1:num_repetitions){
    d <- DGP(n)
    
    est <- estimator(d)
    
    sim[i,1] <- est[1]
    
    sim[i,2] <- est[2]
    
    sim[i+num_repetitions,1] <- est[3]
    
    sim[i+num_repetitions,2] <- NA
    
  }
  
  return(sim)
  
}

#Number of repetitions in our simulation study
num_repetitions <- 200

#True value of the intervention specific mean under the linear model
ATE_Lin <- sim_func1(1e6,1)

#GLM-TMLE simulation under the linear model
sim22 <- sim_study(DGP = sim_func1,
                   estimator = GLM_TMLE, 
                   num_repetitions = num_repetitions,
                   n = 200)

#HAL-TMLE simulation under the linear model
sim21 <- sim_study(DGP = sim_func1,
                   estimator = HAL_TMLE, 
                   num_repetitions = num_repetitions,
                   n = 200)


#True value of the intervention specific mean under the linear model
ATE_jumps <- sim_func2(1e6,1)

#GLM-TMLE simulation under the non-linear model
sim12 <- sim_study(DGP = sim_func2,
                   estimator = GLM_TMLE, 
                   num_repetitions = num_repetitions,
                   n = 200)

#HAL-TMLE simulation under the non-linear model
sim11 <- sim_study(DGP = sim_func2,
                   estimator = HAL_TMLE, 
                   num_repetitions = num_repetitions,
                   n = 200)

simHALTMLEJumpLAST <- sim11
write.csv(simHALTMLEJumpLAST, file = "simHALTMLEJumpLAST", row.names = FALSE)

simGLMTMLEJumpLAST <- sim12
write.csv(simGLMTMLEJumpLAST, file = "simGLMTMLEJumpLAST", row.names = FALSE)

simHALTMLELinLAST <- sim21
write.csv(simHALTMLELinLAST, file = "simHALTMLELinLAST", row.names = FALSE)

simGLMTMLELinLAST <- sim22
write.csv(simGLMTMLELinLAST, file = "simGLMTMLELinLAST", row.names = FALSE)

#Estimate coverage

#Compute summary statistics and estimated densities

GLM.se <- mean(simGLMTMLEJumpLAST$V2[!is.na(simGLMTMLEJumpLAST$V2)])
HAL.se <- mean(simHALTMLEJumpLAST$V2[!is.na(simHALTMLEJumpLAST$V2)])

GLM.CI <- ATE_jumps + c(-1.96,1.96)*GLM.se 
HAL.CI <- ATE_jumps + c(-1.96,1.96)*HAL.se

fun1 <- function(x){dnorm(x,mean = ATE_jumps, sd = GLM.se)}
fun2 <- function(x){dnorm(x,mean = ATE_jumps, sd = HAL.se)}


GLM.se1 <- mean(simGLMTMLELinLAST$V2[!is.na(simGLMTMLELinLAST$V2)])
HAL.se1 <- mean(simHALTMLELinLAST$V2[!is.na(simHALTMLELinLAST$V2)])

GLM.CI1 <- ATE_Lin + c(-1.96,1.96)*GLM.se1 
HAL.CI1 <- ATE_Lin + c(-1.96,1.96)*HAL.se1

fun11 <- function(x){dnorm(x,mean = ATE_Lin, sd = GLM.se1)}
fun21 <- function(x){dnorm(x,mean = ATE_Lin, sd = HAL.se1)}


#Computing bias and discrepancy in median
biasLinGLM <- mean(simGLMTMLELinLAST$V1[is.na(simGLMTMLELinLAST$V2)]) - ATE_Lin
biasLinGLMTMLE <- mean(simGLMTMLELinLAST$V1[!is.na(simGLMTMLELinLAST$V2)]) - ATE_Lin
biasLinHAL <- mean(simHALTMLELinLAST$V1[is.na(simHALTMLELinLAST$V2)]) - ATE_Lin
biasLinHALTMLE <- mean(simHALTMLELinLAST$V1[!is.na(simHALTMLELinLAST$V2)]) - ATE_Lin

biasJumpGLM <- mean(simGLMTMLEJumpLAST$V1[is.na(simGLMTMLEJumpLAST$V2)]) - ATE_jumps
biasJumpGLMTMLE <- mean(simGLMTMLEJumpLAST$V1[!is.na(simGLMTMLEJumpLAST$V2)]) - ATE_jumps
biasJumpHAL <- mean(simHALTMLEJumpLAST$V1[is.na(simHALTMLEJumpLAST$V2)]) - ATE_jumps
biasJumpHALTMLE <- mean(simHALTMLEJumpLAST$V1[!is.na(simHALTMLEJumpLAST$V2)]) - ATE_jumps

medLinGLM <- median(simGLMTMLELinLAST$V1[is.na(simGLMTMLELinLAST$V2)]) - ATE_Lin
medLinGLMTMLE <- median(simGLMTMLELinLAST$V1[!is.na(simGLMTMLELinLAST$V2)]) - ATE_Lin
medLinHAL <- median(simHALTMLELinLAST$V1[is.na(simHALTMLELinLAST$V2)]) - ATE_Lin
medLinHALTMLE <- median(simHALTMLELinLAST$V1[!is.na(simHALTMLELinLAST$V2)]) - ATE_Lin

medJumpGLM <- median(simGLMTMLEJumpLAST$V1[is.na(simGLMTMLEJumpLAST$V2)]) - ATE_jumps
medJumpGLMTMLE <- median(simGLMTMLEJumpLAST$V1[!is.na(simGLMTMLEJumpLAST$V2)]) - ATE_jumps
medJumpHAL <- median(simHALTMLEJumpLAST$V1[is.na(simHALTMLEJumpLAST$V2)]) - ATE_jumps
medJumpHALTMLE <- median(simHALTMLEJumpLAST$V1[!is.na(simHALTMLEJumpLAST$V2)]) - ATE_jumps


#Compute simple models

simsimple1 <- sim_study(DGP = sim_func2simple,
                        estimator = HAL_TMLE,
                        num_repetitions = num_repetitions,
                        n=500)
simsimple2 <- sim_study(DGP = sim_func2simple,
                        estimator = HAL_TMLE,
                        num_repetitions = num_repetitions,
                        n=1000)
simsimple3 <- sim_study(DGP = sim_func2simple,
                        estimator = HAL_TMLE,
                        num_repetitions = num_repetitions,
                        n=2000)

simsimple1GLM <- sim_study(DGP = sim_func2simple,
                        estimator = GLM_TMLE,
                        num_repetitions = num_repetitions,
                        n=500)
simsimple2GLM <- sim_study(DGP = sim_func2simple,
                        estimator = GLM_TMLE,
                        num_repetitions = num_repetitions,
                        n=1000)
simsimple3GLM <- sim_study(DGP = sim_func2simple,
                        estimator = GLM_TMLE,
                        num_repetitions = num_repetitions,
                        n=2000)

write.csv(simsimple1, file = "simsimple1", row.names = FALSE)
write.csv(simsimple2, file = "simsimple2", row.names = FALSE)
write.csv(simsimple3, file = "simsimple3", row.names = FALSE)
write.csv(simsimple1GLM, file = "simsimple1GLM", row.names = FALSE)
write.csv(simsimple2GLM, file = "simsimple2GLM", row.names = FALSE)
write.csv(simsimple3GLM, file = "simsimple3GLM", row.names = FALSE)

#Compute coverage

cover1 <- mean(simsimple1[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple1[1:200,2]) 
               & simsimple1[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple1[1:200,2]))

cover2 <- mean(simsimple2[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple2[1:200,2]) 
               & simsimple2[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple2[1:200,2]))

cover3 <- mean(simsimple3[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple3[1:200,2]) 
               & simsimple3[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple3[1:200,2]))

cover1glm <- mean(simsimple1GLM[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple1GLM[1:200,2]) 
               & simsimple1GLM[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple1GLM[1:200,2]))

cover2glm <- mean(simsimple2GLM[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple2GLM[1:200,2]) 
               & simsimple2GLM[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple2GLM[1:200,2]))

cover3glm <- mean(simsimple3GLM[1:200,1]>=ATE_simplejumps - 1.96*mean(simsimple3GLM[1:200,2]) 
               & simsimple3GLM[1:200,1] <=ATE_simplejumps + 1.96*mean(simsimple3GLM[1:200,2]))

num_obs_cover <- c(500,500,1000,1000,2000,2000)
est_name_cover <- c("HAL-TMLE","GLM-TMLE","HAL-TMLE","GLM-TMLE","HAL-TMLE","GLM-TMLE")

coverage <- rbind(cover1,cover1glm,cover2,cover2glm,cover3,cover3glm)

dfcover <- data.frame(num_obs_cover = as.numeric(num_obs_cover),
                      est_name_cover = as.factor(est_name_cover),
                      coverage)
