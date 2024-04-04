#load necessary packages
library(rstan)
library(data.table)
library(hal9001)
library(dplyr)
library(readr)
library(haldensify)
library(tidyverse)

#set seed for reproducability
set.seed(210399)

#prepare data
pph_subset <- read_csv("Downloads/pph-subset.csv")

n <- 4000

data <- pph_subset
data <- na.omit(data)
data <- data[sample(nrow(data),n),]

data <- data %>%
  mutate(PPHbin = ifelse(PPHbin == "Yes", 1, 0), CS = ifelse(CS == "Yes", 1, 0))

data <- data %>%
  mutate_if(is.character, as.factor)

data <- data %>%
  mutate_all(as.numeric)


covars_name <- c("CS","Instrumental",
                 #"MAlder",
                 "Macrosomia4500",
                 "flerfold",
                 "Abruptio",
                 "PraevBin",
                 "Argumented",
                 #"Praeterm",
                 "Praeecl",
                 "Induced",
                 "Episiotomi",
                 "PrevTotal",
                 #"PrevRetained",
                 "Retained",
                 "PrevAbruptio",
                 "PrevCS")

data_names <- c("PPHbin", "CS","Instrumental",
                #"MAlder",
                "Macrosomia4500",
                "flerfold",
                "Abruptio",
                "PraevBin",
                "Argumented",
                #"Praeterm",
                "Praeecl",
                "Induced",
                "Episiotomi",
                "PrevTotal",
                "PrevAbruptio",
                "PrevCS",
                #"PrevRetained",
                "Retained")

data <- data[, names(data) %in% data_names]
covars_all <- data[, names(data) %in% covars_name]
covars_propscore <- covars_all[, !names(covars_all) == "CS"]


HAL_TMLE_practical <- function(d,covars_all,covars_propscore,a){
  
  fit_trmt <- fit_hal(X = as.matrix(covars_all), Y = d$PPHbin,
                      family = "binomial",
                      lambda = NULL,
                      fit_control = list(cv_select = TRUE,
                                         use_min = FALSE,
                                         prediction_bounds = "default"))
  
  fit_propscore <- fit_hal(X = as.matrix(covars_propscore), Y = d$CS,
                           family = "binomial",
                           lambda = NULL,
                           fit_control = list(cv_select = TRUE,
                                              use_min = FALSE,
                                              prediction_bounds = "default"))
  
  d <- data.table(data)
  
  pred.Q0 <- predict(fit_trmt, type = "response", new_data=copy(d)[, CS:=a])
  pred.Q <- predict(fit_trmt, type = "response", new_data=d)
  pred.G <- predict(fit_propscore, type = "response", new_data=d)
  
  H <- (d$CS==a)/(pred.G^a*(1-pred.G)^(1-a))
  
  fit_tmle <- glm(PPHbin ~ offset(qlogis(pred.Q)) + H - 1, family=binomial, data=d)
  
  eps.hat <- fit_tmle$coef
  
  pred.tmle1 <- plogis(qlogis(pred.Q0) + eps.hat/(pred.G^a*(1-pred.G)^(1-a)))
  pred.tmle <- plogis(qlogis(pred.Q) + eps.hat*H)
  
  tmle.est <- mean(pred.tmle1)
  
  tmle.se <- sqrt(mean((H*(d$PPHbin-pred.tmle) + pred.tmle1 - tmle.est)^2)/nrow(d))
  
  return(c(tmle.est = tmle.est, tmle.se = tmle.se))
}

TSM1 <- HAL_TMLE_practical(d=data,covars_all=covars_all,
                           covars_propscore=covars_propscore,a=1)

TSM0 <- HAL_TMLE_practical(d=data,covars_all=covars_all,
                           covars_propscore=covars_propscore,a=0)

TSM1

TSM1[1] + c(-1.96,1.96)*TSM1[2]

TSM0

TSM0[1] + c(-1.96,1.96)*TSM0[2]



#motivating example
firsttry <- glm(PPHbin~. , family = binomial, data = data)
firsttry$coefficients["CS"]
confint(firsttry)
