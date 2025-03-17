## BJS Jan 2025
## Fitting an SIR model to prevalence time series data, with quantified uncertainty in the parameters
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
  p_SI <- 1 - exp(-beta * I / N * dt)
  p_IR <- 1 - exp(-gamma * dt)
  beta <- parameter()
  gamma <- parameter()
  I0 <- parameter()
  R0 <- parameter()
  N <- 10000

  # comparison to data
  prev <- data()
  active_cases <- prev * N
  active_cases ~ Poisson(I)
}, quiet = TRUE)

filter <- dust_filter_create(sir(), 0, one_run, n_particles = 1000)

pars <- list(beta = 0.6, gamma = 0.5, I0 = 3, R0 = 100)

prior <- monty_dsl({
  beta ~ Exponential(mean=0.15)
  gamma ~ Exponential(mean=0.1)
  I0 ~ Exponential(mean=3)
  R0 ~ Exponential(mean=100)
})

sir_packer <- monty_packer(c("beta", "gamma", "I0", "R0"))

likelihood <- dust_likelihood_monty(filter, sir_packer)

posterior <- prior + likelihood

# sampler with step of 0.01 for beta and gamma, and 1 for I0 and R0
vcv <- matrix(c(0.0005, 0.0003, 0, 0.0005,
                0.0003, 0.0003, 0.0003, 0,
                0, 0.0003, 0.01, 0,
                0.0005, 0, 0, 0.1), 4, 4)
sampler <- monty_sampler_random_walk(vcv)

# samples <- monty_sample(posterior, sampler, 10000, initial = sir_packer$pack(pars))
# plot(samples$density, type="l")
# plot(t(drop(samples$pars)))


# # histogram of gamma and beta distribution, excluding first 2000 samples
# hist(samples$pars[1,4001:10000,], breaks=50, col="blue", xlab="gamma")
# hist(samples$pars[2,4001:10000,], breaks=50, col="red")

# # multipanel plot of histograms of all four parameters
# par(mfrow=c(2,2))
# for (i in 1:4) {
#   hist(samples$pars[i,4001:10000,], breaks=50, col="blue", xlab=names(pars)[i])
# }

# # Plot a bunch of runs of the model against the data
# x <- dust_likelihood_run(filter, list(beta = 0.14, gamma = 0.1, I0 = 1, R0 = 0),save_trajectories = TRUE)
# h <- dust_likelihood_last_trajectories(filter)
# prevalence <- dust_unpack_state(sys,h)$prevalence
# matplot(t(prevalence), type = "l", col = "grey", lty = 1, lwd = 1, xlab = "Time", ylab = "Prevalence", main = "Prevalence trajectories")
# points(one_run$time, one_run$prev, col = "black", pch = 19, cex = 0.5)


sys <- dust_system_create(sir(), pars, n_particles = 30, dt = 0.25)
dust_system_set_state_initial(sys)
# x <- dust_likelihood_run(filter, pars,save_trajectories = TRUE)
# h <- dust_likelihood_last_trajectories(filter)
y <- dust_system_simulate(sys, 1:36)
prevalence <- dust_unpack_state(sys,y)$prevalence
# filtered_prevalence <- dust_unpack_state(sys,h)$prevalence
print(prevalence)
# print(t(prevalence))
# matplot(t(filtered_prevalence), type = "l", col = "red", lty = 1, lwd = 1, xlab = "Time", ylab = "Prevalence", main = "Prevalence trajectories")
plot(one_run$time, one_run$prev, col = "black", pch = 19, cex = 0.5)
for (i in 1:30) {
  matlines(1:36,prevalence[i,], type = "l", col = "grey", lty = 1, lwd = 1)
}
