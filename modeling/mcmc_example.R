## BJS March 2025
## Fitting an seasonal model to prevalence time series data
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

unfilter <- dust_unfilter_create(sir(), 0, one_run)

pars <- list(beta = 4.24e-4, gamma = 1, I0 = 0, R0 = 1.5e3, sigma = 1, phi = 18)

prior <- monty_dsl({
  beta ~ Exponential(mean=0.0007)
  R0 ~ Exponential(mean=1500)
})

sir_packer <- monty_packer(c("beta","R0"), fixed = list(gamma=1, I0 = 0, sigma=1, phi=18))

likelihood <- dust_likelihood_monty(unfilter, sir_packer)

posterior <- prior + likelihood

vcv <- diag(c(3e-9,3e-3),2,2)
sampler <- monty_sampler_random_walk(vcv)

n_samples = 1e7
samples <- monty_sample(posterior, sampler, n_samples, initial = sir_packer$pack(pars))
plot(samples$density, type="l")
# save results
#adjust max print
# options(max.print=1000000)
# reset max print to default
options(max.print=99999)
# print value of beta
plot(t(drop(samples$pars)), pch = 19, col = "#00000033")
write.csv(samples$pars, "samples_deterministic_2d.csv")
# # calculate acceptance rate - i.e. 1 - the proportion of values in pars that are identical to the previous value
acceptance_rate <- 1 - sum(apply(samples$pars,1,duplicated))/(n_samples*2)
print(paste("Acceptance rate: ", acceptance_rate))