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

filter <- dust_filter_create(seir(), 0, one_run, n_particles = 1000)

pars <- list(beta = 0.4, gamma = 0.5, alpha=0.5, E0 = 100, R0 = 100, sigma=0.1, phi=10)

prior <- monty_dsl({
  beta ~ Exponential(mean=0.15)
  alpha ~ Exponential(mean=0.1)
  gamma ~ Exponential(mean=0.1)
  E0 ~ Exponential(mean=3)
  R0 ~ Exponential(mean=100)
})

seir_packer <- monty_packer(c("beta", "gamma", "alpha", "E0", "R0"), fixed = list(sigma=1, phi=26))
print(seir_packer$pack(pars))

likelihood <- dust_likelihood_monty(filter, seir_packer)

posterior <- prior + likelihood

# sampler with step of 0.01 for beta and gamma, and 1 for E0 and R0
vcv <- matrix(c(0.0005, 0, 0.0003, 0, 0,
                0, 0.0003, 0, 0, 0,
                0.0003, 0, 0.0003, 0, 0,
                0, 0, 0, 0.01, 0,
                0, 0, 0, 0, 0.01), 5, 5)
sampler <- monty_sampler_random_walk(vcv)

samples <- monty_sample(posterior, sampler, 10000, initial = seir_packer$pack(pars))
plot(samples$density, type="l")
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
# sys <- dust_system_create(seir(), pars, n_particles = 20, dt = 0.25)
# dust_system_set_state_initial(sys)
# x <- dust_likelihood_run(filter, pars,save_trajectories = TRUE)
# h <- dust_likelihood_last_trajectories(filter)
# prevalence <- dust_unpack_state(sys,h)$prevalence
# matplot(t(prevalence), type = "l", col = "grey", lty = 1, lwd = 1, xlab = "Time", ylab = "Prevalence", main = "Prevalence trajectories")
# points(one_run$time, one_run$prev, col = "black", pch = 19, cex = 0.5)