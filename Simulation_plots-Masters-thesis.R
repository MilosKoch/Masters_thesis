library(ggplot2)
library(tidyverse)

est_name <- c(rep("HAL-TMLE",200),rep("HAL",200),rep("GLM-TMLE",200),rep("GLM",200))
est_name <- as.factor(est_name)

dfJumps <- rbind(simHALTMLEJumpLAST,simGLMTMLEJumpLAST)

dfJumps <- cbind(dfJumps,est_name)
dfJumps <- dfJumps %>%
  mutate(est_name = factor(est_name))

ggplot(data = dfJumps, aes(x = est_name, y = V1)) +
  geom_boxplot(aes(fill=est_name)) +
  geom_hline(yintercept = ATE_jumps, color = "red") +
  labs(x = "Estimator", y = "Values", fill = "Estimator")

GLM.se <- mean(dfJumps$V2[dfJumps$est_name=="GLM-TMLE"])
HAL.se <- mean(dfJumps$V2[dfJumps$est_name=="HAL-TMLE"])

GLM.CI <- ATE_jumps + c(-1.96,1.96)*GLM.se 
HAL.CI <- ATE_jumps + c(-1.96,1.96)*HAL.se

fun1 <- function(x){dnorm(x,mean = ATE_jumps, sd = GLM.se)}
fun2 <- function(x){dnorm(x,mean = ATE_jumps, sd = HAL.se)}

ggplot(data = dfJumps, aes(x = V1, fill = est_name)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.025, position = "dodge") +
  facet_wrap(~est_name) +
  geom_vline(aes(xintercept = ATE_jumps), color = "red") +
  geom_function(data = data.frame(V1 = 0, est_name = "GLM-TMLE"),
                fun = fun1) +
  geom_function(data = data.frame(V1 = 0, est_name = "HAL-TMLE"),
                fun = fun2) + 
  geom_vline(data = dfJumps %>% filter(est_name == "GLM-TMLE"),
             aes(xintercept=GLM.CI[1]), linetype = "dashed", color = "red") +
  geom_vline(data = dfJumps %>% filter(est_name == "GLM-TMLE"),
             aes(xintercept=GLM.CI[2]), linetype = "dashed", color = "red") +
  geom_vline(data = dfJumps %>% filter(est_name == "HAL-TMLE"),
             aes(xintercept=HAL.CI[1]), linetype = "dashed", color = "red") +
  geom_vline(data = dfJumps %>% filter(est_name == "HAL-TMLE"),
             aes(xintercept=HAL.CI[2]), linetype = "dashed", color = "red") +
  labs(x = "Values", y = "Density", fill = "Estimator")

simGLMTMLELinLAST <- read_csv("simGLMTMLELinLAST")
simHALTMLELinLAST <- read_csv("simHALTMLELinLAST")

est_name1 <- c(rep("HAL-TMLE",200),rep("HAL",200),rep("GLM-TMLE",200),rep("GLM",200))


dfLin <- rbind(simHALTMLELinLAST,simGLMTMLELinLAST)

dfLin <- cbind(dfLin,est_name1)
dfLin <- dfLin %>%
  mutate(est_name1 = factor(est_name1))

ggplot(data = dfLin, aes(x = est_name1, y = V1)) +
  geom_boxplot(aes(fill=est_name1)) +
  geom_hline(yintercept = ATE_Lin, color = "red") +
  labs(x = "Estimator", y = "Values", fill = "Estimator")

GLM.se1 <- mean(dfLin$V2[dfLin$est_name1=="GLM-TMLE"])
HAL.se1 <- mean(dfLin$V2[dfLin$est_name1=="HAL-TMLE"])

GLM.CI1 <- ATE_Lin + c(-1.96,1.96)*GLM.se1 
HAL.CI1 <- ATE_Lin + c(-1.96,1.96)*HAL.se1

fun11 <- function(x){dnorm(x,mean = ATE_Lin, sd = GLM.se1)}
fun21 <- function(x){dnorm(x,mean = ATE_Lin, sd = HAL.se1)}

ggplot(data = dfLin, aes(x = V1, fill = est_name1)) +
  geom_histogram(aes(y=after_stat(density)), binwidth = 0.025, position = "dodge") +
  facet_wrap(~est_name1) +
  geom_vline(aes(xintercept = ATE_Lin), color = "red") +
  geom_function(data = data.frame(V1 = 0, est_name1 = "GLM-TMLE"),
                fun = fun11) +
  geom_function(data = data.frame(V1 = 0, est_name1 = "HAL-TMLE"),
                fun = fun21) + 
  geom_vline(data = dfLin %>% filter(est_name1 == "GLM-TMLE"),
             aes(xintercept=GLM.CI1[1]), linetype = "dashed", color = "red") +
  geom_vline(data = dfLin %>% filter(est_name1 == "GLM-TMLE"),
             aes(xintercept=GLM.CI1[2]), linetype = "dashed", color = "red") +
  geom_vline(data = dfLin %>% filter(est_name1 == "HAL-TMLE"),
             aes(xintercept=HAL.CI1[1]), linetype = "dashed", color = "red") +
  geom_vline(data = dfLin %>% filter(est_name1 == "HAL-TMLE"),
             aes(xintercept=HAL.CI1[2]), linetype = "dashed", color = "red") +
  labs(x = "Values", y = "Density", fill = "Estimator")


?geom_ribbon
CI <- list(NA,c(ATE_jumps - 1.96*HAL.se,ATE_jumps + 1.96*HAL.se),
           NA,c(ATE_jumps - 1.96*GLM.se,ATE_jumps + 1.96*GLM.se))
est_name <- c("HAL","HAL-TMLE","GLM","GLM-TMLE")
est_name <- as.factor(est_name)

dfCI <- as.data.frame(cbind(CI,est_name))

dfCI <- dfCI %>%
  mutate(est_name = factor(est_name))

df <- left_join(dfJumps,dfCI,copy=TRUE)

HAL_width <- 0.7372438 - 0.5738202

ggplot(data = dfJumps, aes(x = V1, fill = est_name)) +
  geom_histogram(aes(y = after_stat(density)), binwidth = 0.025, position = "dodge") +
  facet_wrap(~ est_name) +
  geom_vline(aes(xintercept = ATE_jumps), color = "red") +
  geom_function(data = data.frame(V1 = 0, est_name = "GLM-TMLE"), fun = fun1) +
  geom_function(data = data.frame(V1 = 0, est_name = "HAL-TMLE"), fun = fun2) +
  geom_ribbon(data = dfJumps[dfJumps$est_name == "GLM-TMLE",], 
              aes(ymin = 0, ymax = Inf, xmin = GLM.CI[1], xmax = GLM.CI[2]), 
              fill = "gray", alpha = 0.5) +
  geom_ribbon(data = dfJumps[dfJumps$est_name == "HAL-TMLE",], 
              aes(ymin = 0, ymax = Inf, xmin = 0.5738202, xmax = 0.5738202 + HAL_width), 
              fill = "gray", alpha = 0.5) +
  labs(x = "Values", y = "Density", fill = "Estimator")
