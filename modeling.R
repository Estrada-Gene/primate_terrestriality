####### statistical modeling for Estrada & Marshall (2024) #######
library(tidyverse)
library(brms)


#primate data 
primates <- read.csv(file = "data/primates_final.csv", header = TRUE)