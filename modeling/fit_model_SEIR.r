## BJS Jan 2025
## Fitting an SEIR model to prevalence time series data, with quantified uncertainty in the parameters
library(odin2)
library(dust2)
library(monty)

# download data
data <- read.csv("sim_200.csv")
colnames(data)[1] <- "time"

# Stochastic discrete SEIR model
seir <- odin({
  initial(S) <- N - E0 - R0
  initial(E) <- E0
  initial(I) <- 0
  initial(R) <- R0
  initial(prevalence, zero_every = 1) <- 0
  update(S) <- S - n_SE
  update(E) <- E + n_SE - n_EI
  update(I) <- I + n_EI - n_IR
  update(R) <- R + n_IR
  update(prevalence) <- (I+E)/N
  n_SE <- Binomial(S, p_SE)
  n_EI <- Binomial(E, p_EI)
  n_IR <- Binomial(I, p_IR)
  p_SE <- 1 - exp(-beta * (1 + sigma*cos(2*pi*(time-phi)/365)) * I / N * dt)
  p_EI <- 1 - exp(-alpha * dt)
  p_IR <- 1 - exp(-gamma * dt)
  beta <- parameter()
  alpha <- parameter()
  gamma <- parameter()
  sigma <- parameter()
  phi <- parameter()
  E0 <- parameter()
  R0 <- parameter()
  N <- 10000

  # comparison to data
  prev <- data()
  active_cases <- prev * N
  active_cases ~ Poisson(I)
}, quiet = TRUE)

pars <- list(beta = 3.7, gamma = 2, alpha=1, E0 = 7, R0 = 6000, sigma=0.45, phi=3)

prior <- monty_dsl({
  beta ~ Exponential(mean=0.15)
  gamma ~ Exponential(mean=0.1)
  E0 ~ Exponential(mean=3)
  R0 ~ Exponential(mean=1500)
  sigma ~ Uniform(0, 0.5)
  phi ~ Uniform(0, 365)
})

seir_packer <- monty_packer(c("beta", "gamma", "E0", "R0", "sigma", "phi"), fixed = list(alpha=1))
print(seir_packer$pack(pars))

one_run <- data[1:36,]
unfilter <- dust_unfilter_create(seir(), 0, one_run)
likelihood <- dust_likelihood_monty(unfilter, seir_packer)

posterior <- prior + likelihood

# sampler with step of 0.01 for beta and gamma, and 1 for E0 and R0
# vcv <- matrix(c(0.004, 0, 0, 0,
#                 0, 0.005, 0, 0,
#                 0, 0, 0.1, 0,
#                 0, 0, 0, 0.1), 4, 4)
vcv <- matrix(c(0.0005, 0, 0, 0, 0.0002, 0,
                0, 0.0005, 0, 0, 0, 0,
                0, 0, 0.1, 0, 0, 0,
                0, 0, 0, 0.1, 0, 0,
                0.0002, 0, 0, 0, 0.0002, 0,
                0, 0, 0, 0, 0, 0.005), 6, 6)
sampler <- monty_sampler_random_walk(vcv)

n_samples = 5e6
# samples <- monty_sample(posterior, sampler, n_samples, initial = seir_packer$pack(pars))
# plot(samples$density, type="l")
# print(samples$pars[((n_samples-1)*nrow(vcv)+1):(n_samples*nrow(vcv))])


# acceptance_rate <- 1 - sum(apply(samples$pars,1,duplicated))/(n_samples*nrow(vcv))
# print(paste("Acceptance rate: ", acceptance_rate))
# # plot(t(drop(samples$pars)))

# # save results
# write.csv(samples$pars, "samples_SEIR.csv")

# # calculate medians and 95% credible intervals of the parameters
# pars <- samples$pars
# med <- apply(pars, 1, median)
# ci <- apply(pars, 1, function(x) quantile(x, c(0.025, 0.975)))
# print(med)
# print(ci)


# # for each 36 weeks of data, save the samples in a separate file
# n_sims <- nrow(data)/36
# for (i in 1:n_sims) {
#   one_run <- data[((i-1)*36+1):(i*36),]
#   unfilter <- dust_unfilter_create(seir(), 0, one_run)
#   pars <- list(beta = 3.7, gamma = 2, alpha=1, E0 = 7, R0 = 6000, sigma=0.45, phi=3)
#   likelihood <- dust_likelihood_monty(unfilter, seir_packer)
#   posterior <- prior + likelihood
#   samples <- monty_sample(posterior, sampler, n_samples, initial = seir_packer$pack(pars))
#   write.csv(samples$pars, paste("samples_SEIR_", i, ".csv"))
# }

# parallel version
library(parallel)
n_cores <- detectCores()
# cl <- makeCluster(n_cores)
# clusterExport(cl, c("data","seir", "unfilter", "prior", "seir_packer", "sampler", "n_samples", "pars"))
# clusterEvalQ(cl, library(odin2))
# clusterEvalQ(cl, library(dust2))
# clusterEvalQ(cl, library(monty))
# clusterEvalQ(cl, library(parallel))
# clusterEvalQ(cl, library(stats))
# clusterEvalQ(cl, library(base))
# # define function to run the model
run_model <- function(i, n_samples) {
  one_run <- data[((i-1)*36+1):(i*36),]
  unfilter <- dust_unfilter_create(seir(), 0, one_run)
  pars <- list(beta = 3.7, gamma = 2, alpha=1, E0 = 7, R0 = 6000, sigma=0.45, phi=3)
  likelihood <- dust_likelihood_monty(unfilter, seir_packer)
  posterior <- prior + likelihood
  samples <- monty_sample(posterior, sampler, n_samples, initial = seir_packer$pack(pars))
  pars <- samples$pars
  med <- apply(pars, 1, median)
  ci <- apply(pars, 1, function(x) quantile(x, c(0.025, 0.975)))
  # write to median and CI to same file
  write.csv(rbind(med, ci), paste0("MCMC_runs_SEIR/med_ci", i, ".csv"))
  # write.csv(samples$pars, paste0("MCMC_runs_SEIR/pars", i, ".csv"))
  # write.csv(samples$density, paste0("MCMC_runs_SEIR/density", i, ".csv"))
}
# # run the model in parallel
# parLapply(cl, 1:300, run_model)
# stopCluster(cl)

## use mclapply instead of parLapply
model_partial <- function(i) {
  run_model(i, 2e6)
}
mclapply(1:500, model_partial, mc.cores = n_cores)
# model_partial <- function(i) {
#   run_model(i, 1e6)
# }
# mclapply(5:300, model_partial, mc.cores = n_cores)

# # histogram of gamma and beta distribution, excluding first 2000 samples
# hist(samples$pars[1,4001:10000,], breaks=50, col="blue", xlab="gamma")
# hist(samples$pars[2,4001:10000,], breaks=50, col="red")

# # multipanel plot of histograms of all four parameters
# par(mfrow=c(2,2))
# for (i in 1:4) {
#   hist(samples$pars[i,4001:10000,], breaks=50, col="blue", xlab=names(pars)[i])
# }

# # Plot a bunch of runs of the model against the data
# sys <- dust_system_create(seir(), pars, n_particles = 20, dt = 0.25)
# dust_system_set_state_initial(sys)
# x <- dust_likelihood_run(unfilter, pars,save_trajectories = TRUE)
# h <- dust_likelihood_last_trajectories(unfilter)
# prevalence <- dust_unpack_state(sys,h)$prevalence
# matplot(t(prevalence), type = "l", col = "grey", lty = 1, lwd = 1, xlab = "Time", ylab = "Prevalence", main = "Prevalence trajectories")
# points(one_run$time, one_run$prev, col = "black", pch = 19, cex = 0.5)