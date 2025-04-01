## BJS Jan 2025
## Fitting an SEIR model to prevalence time series data, with quantified uncertainty in the parameters
library(odin2)
library(dust2)
library(monty)
library(parallel)

args <- commandArgs(trailingOnly = TRUE)
simulation_file <- args[1] # e.g. "sim_200"
n_simulations <- as.numeric(args[2]) # e.g. 500
MCMC_length <- as.numeric(args[3]) # recommended value 2e6

# download data
data <- read.csv(paste0(simulation_file,".csv"))
colnames(data)[1] <- "time"
sim_length <- nrow(data)%/%n_simulations

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

prior <- monty_dsl({
  beta ~ Exponential(mean=0.15)
  gamma ~ Exponential(mean=0.1)
  E0 ~ Exponential(mean=3)
  R0 ~ Exponential(mean=1500)
  sigma ~ Uniform(0, 0.5)
  phi ~ Uniform(0, 365)
})

seir_packer <- monty_packer(c("beta", "gamma", "E0", "R0", "sigma", "phi"), fixed = list(alpha=1))

# step sizes chosen to give acceptance rate of around 0.23
vcv <- matrix(c(0.0008, 0, 0, 0, 0.0004, 0,
                0, 0.0008, 0, 0, 0, 0,
                0, 0, 0.1, 0, 0, 0,
                0, 0, 0, 1, 0, 0,
                0.0004, 0, 0, 0, 0.0004, 0,
                0, 0, 0, 0, 0, 0.004), 6, 6)
sampler <- monty_sampler_random_walk(vcv)

### run MCMCs in parallel for all simulations
# define function to run the model
run_model <- function(i, n_samples) {
  one_run <- data[((i-1)*sim_length+1):(i*sim_length),]
  unfilter <- dust_unfilter_create(seir(), 0, one_run)
  pars <- list(beta = 3.7, gamma = 2, alpha=1, E0 = 7, R0 = 1500, sigma=0.45, phi=3)
  likelihood <- dust_likelihood_monty(unfilter, seir_packer)
  posterior <- prior + likelihood
  samples <- monty_sample(posterior, sampler, n_samples, initial = seir_packer$pack(pars))
  pars_samples <- samples$pars[,(n_samples%/%4):n_samples,1]
  med <- apply(pars_samples, 1, median)
  ci <- apply(pars_samples, 1, function(x) quantile(x, c(0.025, 0.975)))
  attack_rates <- numeric((n_samples%/%4)*3)
  for (par_n in 1:((n_samples%/%4)*3)) {
    sys_pars <- list(beta = pars_samples[1,par_n], gamma = pars_samples[2,par_n], alpha=1, E0 = pars_samples[3,par_n], R0 = pars_samples[4,par_n], sigma=pars_samples[5,par_n], phi=pars_samples[6,par_n]) 
    sys <- dust_system_create(seir(), sys_pars, dt = 0.25,deterministic = TRUE)
    dust_system_set_state_initial(sys)
    time <- 0:150/4
    y <- dust_system_simulate(sys, time)
    attack_rates[par_n] <- dust_unpack_state(sys,y)$R[151]-dust_unpack_state(sys,y)$R[1]
  }
  attack_med <- median(attack_rates)
  attack_ci <- quantile(attack_rates, c(0.025, 0.975))
  # write to median and CI to the same file
  med_ci <- rbind(i,med, ci)
  attack_med_ci <- t(c(i,attack_med, attack_ci))
  write.table(med_ci, file = paste0("parameter_med_ci_",simulation_file,".csv"), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
  write.table(attack_med_ci, file = paste0("attack_med_ci",simulation_file,".csv"), append = TRUE, sep = ",", col.names = FALSE, row.names = FALSE)
}

# set MCMC length
model_partial <- function(i) {
  run_model(i, MCMC_length)
}
# write header to the file
header <- t(c("n_simulation", "beta_med", "beta_lower", "beta_upper", "gamma_med", "gamma_lower", "gamma_upper", "E0_med", "E0_lower", "E0_upper", "R0_med", "R0_lower", "R0_upper", "sigma_med", "sigma_lower", "sigma_upper", "phi_med", "phi_lower", "phi_upper"))
write.table(header, file = paste0("parameter_med_ci_",simulation_file,".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
header_attack <- t(c("n_simulation", "median", "lower", "upper"))
write.table(header_attack, file = paste0("attack_med_ci",simulation_file,".csv"), sep = ",", col.names = FALSE, row.names = FALSE)
# run the model
mclapply(1:n_simulations, model_partial, mc.cores = detectCores())

# calculate median and ci of attack rate medians
attack_med_ci <- read.csv(paste0("attack_med_ci",simulation_file,".csv"), header = TRUE)
attack_med <- attack_med_ci$median
print(attack_med)
print(quantile(attack_med, c(0.025, 0.5, 0.975)))