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
                                      use_min = TRUE,
                                      lambda.min.ratio = 1e-04,
                                      prediction_bounds = "default"))
  
  #HAL fit of the propensity score
  G_fit <- fit_hal(X = as.matrix(treatment_data[, !(names(treatment_data)) %in% "A"]),
                   Y = treatment_data$A,
                   family = "binomial",
                   lambda = 0.001,
                   fit_control = list(cv_select = FALSE,
                                      use_min = TRUE,
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
  
  #GLM fit of the conditional mean
  Q_fit <- glm(Y~A+X1+X2+X3, family = "binomial", d)
  
  #GLM fit of the propensity score
  G_fit <- glm(A~X1+X2+X3, family = "binomial", d)
  
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

#Alternative GLM-TMLE for non-linear model
GLM_TMLE <- function(d){
  
  d <- as.data.frame(d)
  
  #GLM fit of the conditional mean
  Q_fit <- glm(Y~A+X1+X2+X3, family = "binomial", d)
  
  #GLM fit of the propensity score
  G_fit <- glm(A~X1+X2+X3, family = "binomial", d)
  
  d <- data.table(d)
  
  #HAL estimate of the conditional mean under A=1
  d[, pred.Q0:=predict(Q_fit, type = "response",new_data=copy(d)[, A:=1])]
  
  #GLM estimate of the intervention specific mean
  glm.est <- d[, mean(pred.Q0)]
  
  #GLM estimate of the conditional mean
  d[, pred.Q:=predict(Q_fit, type = "response", new_data=d)]
  
  #GLm estimate of the propensity score
  d[, pred.G:=predict(G_fit, type = "response" ,new_data=d)]
  
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

