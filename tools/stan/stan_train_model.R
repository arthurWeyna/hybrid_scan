library(rstan)
args = commandArgs(trailingOnly=TRUE)

stan_model(file = args[1], auto_write=T)
