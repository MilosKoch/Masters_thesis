library(ggplot2)
library(tidyverse)
library(patchwork)
source("Data_generation-Masters_thesis")
source("Estimation_Simulation-Masters_thesis")

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



est_name_simple <- c(rep("HAL-TMLE",200),rep("GLM-TMLE",200))

dfsimple1 <- rbind(drop_na(simsimple1),drop_na(simsimple1GLM))
dfsimple2 <- rbind(drop_na(simsimple2),drop_na(simsimple2GLM))
dfsimple3 <- rbind(drop_na(simsimple3),drop_na(simsimple3GLM))

dfsimple1 <- cbind(dfsimple1,est_name_simple)
dfsimple2 <- cbind(dfsimple2,est_name_simple)
dfsimple3 <- cbind(dfsimple3,est_name_simple)

p1 <- ggplot(data = dfsimple1, aes(x = V1, color= est_name_simple)) +
  geom_density() +
  stat_function(fun = fun3, color = "black") +
  labs(x = "Values", y = "Density", color = "Estimator") +
  theme(legend.position="none") + 
  ggtitle("500 observations")

p2 <- ggplot(data = dfsimple2, aes(x = V1, color= est_name_simple)) +
  geom_density() +
  stat_function(fun = fun4, color = "black") +
  labs(x = "Values", y = "Density", color = "Estimator") +
  theme(legend.position="none") + 
  ggtitle("1000 observations")

p3 <- ggplot(data = dfsimple3, aes(x = V1, color= est_name_simple)) +
  geom_density() +
  stat_function(fun = fun5, color = "black") +
  labs(x = "Values", y = "Density", color = "Estimator") + 
  ggtitle("2000 observations")



p4 <- ggplot(data = dfcover, aes(x=num_obs_cover, y = coverage, color = est_name_cover)) +
  geom_line() +
  geom_point(aes(shape = factor(est_name_cover)), size = 3) +
  geom_hline(yintercept = 0.95) +
  labs(x = "Number of observatiosn", y = "Coverage", color = "Estimator") +
  guides(shape = "none")

p1 + p3 + p2 + p4 +
  plot_layout(ncol = 2)


p4
dfsimple <- dfsimple %>%
  mutate(est_name_simple = factor(est_name_simple))

fun3 <- function(x){dnorm(x,mean=ATE_simplejumps,sd=mean(simsimple1$V2[1:200]))}
fun4 <- function(x){dnorm(x,mean=ATE_simplejumps,sd=mean(simsimple2$V2[1:200]))}
fun5 <- function(x){dnorm(x,mean=ATE_simplejumps,sd=mean(simsimple3$V2[1:200]))}



ggplot(data = dfsimple, aes(x = V1, color= est_name_simple)) +
  facet_wrap(~num_obs) +
  geom_density() +
  stat_function(fun = fun5, color = "black")
  
  
  
  geom_function(data = data.frame(V1 = 0, num_obs = "500"),
                fun = fun3) +
  geom_function(data = data.frame(V1 = 0, num_obs = "1000"),
                fun = fun4) +
  geom_function(data = data.frame(V1 = 0, num_obs = "2000"),
                fun = fun5)

ggplot(data = dfsimple, aes(x = V1, colour = est_name_simple)) +
  facet_wrap(~num_obs) +
  geom_density() +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("500")),
                fun = fun3) +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("1000")),
                fun = fun4) +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("2000")),
                fun = fun5)

ggplot(data = dfsimple, aes(x = V1, color = est_name_simple)) +
  facet_grid(~factor(num_obs, levels = c("500", "1000", "2000"))) +
  geom_density() +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("500")),
                fun = fun3, aes(fill = num_obs)) +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("1000")),
                fun = fun4, aes(fill = num_obs)) +
  geom_function(data = data.frame(V1 = 0, num_obs = factor("2000")),
                fun = fun5, aes(fill = num_obs))
