## BJS Jan 2025
## Fitting an SEIR model to prevalence time series data, with quantified uncertainty in the parameters
library(odin2)
library(dust2)
library(monty)

# download data
data <- read.csv("sim_200.csv")
one_run <- data[1:36,]
# rename "week" column to "time"
colnames(one_run)[1] <- "time"

# Stochastic discrete SIR model
sir <- odin({
  initial(S) <- N - I0 - R0
  initial(I) <- I0
  initial(R) <- R0
  initial(prevalence, zero_every = 1) <- 0
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  update(prevalence) <- I/N
  n_SI <- Binomial(S, p_SI)
  n_IR <- Binomial(I, p_IR)
  p_SI <- 1 - exp((-beta * (1 + sigma*cos(2*pi*(time-1-phi)/52)) ^6) * dt)
  p_IR <- 1 - exp(-gamma * dt)
  beta <- parameter()
  gamma <- parameter()
  sigma <- parameter()
  phi <- parameter()
  I0 <- parameter()
  R0 <- parameter()
  N <- 10000

  # comparison to data
  prev <- data()
  active_cases <- prev * N
  active_cases ~ Poisson(I)
}, quiet = TRUE)

filter <- dust_filter_create(sir(), 0, one_run, n_particles = 1000)

pars <- list(beta = 4.24e-4, gamma = 1, I0 = 0, R0 = 1.5e3, sigma = 1, phi = 18)

prior <- monty_dsl({
  beta ~ Exponential(mean=0.0007)
  # R0 ~ Exponential(mean=1500)
  # sigma ~ Exponential(mean=1)
  # phi ~ Exponential(mean=15)
})

sir_packer <- monty_packer(c("beta"), fixed = list(gamma=1, I0 = 0,sigma=1, phi=18, R0=1500))

likelihood <- dust_likelihood_monty(filter, sir_packer)

posterior <- prior + likelihood

# sampler with step of 0.01 for beta and gamma, and 1 for I0 and R0
acceptances <- c()
step_sizes <- c(1,1e-1,1e-2,1e-3,1e-13,1e-14,1e-15,1e-16)
for (step_size in step_sizes) {
  vcv <- diag(c(step_size),1,1)
  sampler <- monty_sampler_random_walk(vcv)
  n_samples <- 1e4
  samples <- monty_sample(posterior, sampler, n_samples, initial = sir_packer$pack(pars))
  acceptance_rate <- 1 - sum(apply(samples$pars,1,duplicated))/n_samples
  print(paste("Step size: ", step_size, "Acceptance rate: ", acceptance_rate))
  acceptances <- c(acceptances, acceptance_rate)
  # print(apply(samples$pars,1,mean))
}
# log x axis
plot(step_sizes, acceptances, type="l", log="x", xlab="Step size", ylab="Acceptance rate")
# vcv <- diag(c(1e-11),1,1)
# # vcv <- diag(c(1e-10,1e-7),2,2)
# # vcv <- matrix(c(0.0005, 0.0003, 0, 0.0005,
# #                 0.0003, 0.0003, 0.0003, 0,
# #                 0, 0.0003, 0.01, 0,
# #                 0.0005, 0, 0, 0.1), 4, 4)
# sampler <- monty_sampler_random_walk(vcv)

# n_samples = 1e5
# samples <- monty_sample(posterior, sampler, n_samples, initial = sir_packer$pack(pars))
# plot(samples$density, type="l")
# # save results
# # write.csv(samples$pars, "samples_smaller.csv")
# # calculate acceptance rate - i.e. 1 - the proportion of values in pars that are identical to the previous value
# acceptance_rate <- 1 - sum(apply(samples$pars,1,duplicated))/n_samples
# print(acceptance_rate)
# # print(samples$pars)
# print(apply(samples$pars,1,mean))
# plot(t(drop(samples$pars)))

# # Plot a bunch of runs of the model against the data
# sys <- dust_system_create(sir(), pars, n_particles = 30, dt = 0.25)
# dust_system_set_state_initial(sys)
# # x <- dust_likelihood_run(filter, pars,save_trajectories = TRUE)
# # h <- dust_likelihood_last_trajectories(filter)
# y <- dust_system_simulate(sys, 1:36)
# prevalence <- dust_unpack_state(sys,y)$prevalence
# # filtered_prevalence <- dust_unpack_state(sys,h)$prevalence
# print(prevalence)
# # print(t(prevalence))
# # matplot(t(filtered_prevalence), type = "l", col = "red", lty = 1, lwd = 1, xlab = "Time", ylab = "Prevalence", main = "Prevalence trajectories")
# plot(one_run$time, one_run$prev, col = "black", pch = 19, cex = 0.5)
# for (i in 1:30) {
#   matlines(1:36,prevalence[i,], type = "l", col = "grey", lty = 1, lwd = 1)
# }
