library(rstan)

#Initial values for sampling
init_fun <- function(...) list(theta=0.001, gamma=0.000001)

#This function operates bayesian estimation with stan, computes simple statistics about the input file, and outputs results in list format.
stan_div_estim <- function(model, k, l, chains = 4, warmup = 1000, iter = 2000, cores = 4, refresh = 0){
	difdata <- list(N=length(k), k=k, l=l)
    nb_loc = length(k)
	He_mean <- mean(k/l, na.rm=T)
	He_var <- var(k/l, na.rm=T)
	fit <- sampling(
	     object = model,  # Stan program
	     data = difdata,  # named list of data
         init=init_fun,   # initial values
	     chains = chains, # number of Markov chains
	     warmup = warmup, # number of warmup iterations per chain
	     iter = iter,     # total number of iterations per chain
	     cores = cores,   # number of cores (could use one per chain)
	     refresh = refresh# no progress shown
	)
    postdist <- do.call("rbind", lapply(fit@sim$samples[1:chains], function(x) as.data.frame(x[c("theta", "gamma", "lp__")])[(warmup+1):iter,]))
    summ <- summary(fit)$summary
	rez <- c(nb_loc=nb_loc, He_mean=He_mean, He_var=He_var) 
	for(i in 1:nrow(summ)){
		lin <- summ[i,]
		names(lin) <- paste("stan", sub("__$", "", rownames(summ)[i]), names(lin), sep="_")
		rez <- c(rez, lin)
	}
	rez <- c(rez, chains)
	names(rez)[length(rez)] <- "stan_chains"
	return(list(rez, postdist))
}	


#1: absolute path to input file (required format described in README)
#2: output file with estimations 
#3: output file with full posteriors
#4: absolute path to rds file produced by stan_train_model.R
#5: number of Markov chains to run
#6: number of warm-up iterations for each chain
#7: number of iterations for each chain
#8: cpus 
args <- commandArgs(trailingOnly=TRUE)

#read input files
dataf <- read.table(args[1], sep=" ", header=F)
model <- readRDS(args[4])

if(nrow(dataf)>1) {
dataf[,3] <- dataf[,3]+0.00001
#run estimation
rez <- stan_div_estim(model=model, k=dataf[,ncol(dataf)], l=dataf[,ncol(dataf)-1], chains = as.numeric(args[5]), warmup = as.numeric(args[6]), iter = as.numeric(args[7]), cores = as.numeric(args[8]), refresh = 0)

#output estimates
estim <- c(args[1], rez[[1]])
names(estim)[1] <- "file"
estim <- t(as.matrix(estim))
print(estim)
write.table(x=estim, file=args[2], quote=F, sep=" ", dec=".", row.names=F, col.names=T) 

#output full posteriors
post <- rez[[2]]
names(post)[3] <- "lp"
write.table(x=post, file=args[3], quote=F, sep=" ", dec=".", row.names=F, col.names=T) 

} else {

estim <- rbind(c("file", "nb_loc", "He_mean", "He_var", "stan_theta_mean", "stan_theta_se_mean", "stan_theta_sd", "stan_theta_2.5%", "stan_theta_25%", "stan_theta_50%", "stan_theta_75%", "stan_theta_97.5%", "stan_theta_n_eff", "stan_theta_Rhat", "stan_gamma_mean", "stan_gamma_se_mean", "stan_gamma_sd", "stan_gamma_2.5%", "stan_gamma_25%", "stan_gamma_50%", "stan_gamma_75%", "stan_gamma_97.5%", "stan_gamma_n_eff", "stan_gamma_Rhat", "stan_ratio_mean", "stan_ratio_se_mean", "stan_ratio_sd", "stan_ratio_2.5%", "stan_ratio_25%", "stan_ratio_50%", "stan_ratio_75%", "stan_ratio_97.5%", "stan_ratio_n_eff", "stan_ratio_Rhat", "stan_lp_mean", "stan_lp_se_mean", "stan_lp_sd", "stan_lp_2.5%", "stan_lp_25%", "stan_lp_50%", "stan_lp_75%", "stan_lp_97.5%", "stan_lp_n_eff", "stan_lp_Rhat", "stan_chains"), "NA")
estim[2, 1] <- args[1]
estim[2, 2] <- as.character(nrow(dataf))
write.table(x=estim, file=args[2], quote=F, sep=" ", dec=".", row.names=F, col.names=F) 

post <- rbind(c("theta", "gamma", "lp"), "NA")
write.table(x=post, file=args[3], quote=F, sep=" ", dec=".", row.names=F, col.names=F) 

}
