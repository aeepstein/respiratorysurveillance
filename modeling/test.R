# plot density trajectories for four runs
# "MCMC_runs_SEIR/density", i, ".csv"

# for (i in 1:4) {
#   density <- read.csv(paste0("MCMC_runs_SEIR/density", i, ".csv"))
#   plot(density, type="l", col="blue", xlab="Time", ylab="Density")
# }

# four panel version
par(mfrow=c(2,2))
for (i in 1:4) {
  density <- read.csv(paste0("MCMC_runs_SEIR/density", i, ".csv"))
  plot(density, type="l", col="blue", xlab="Index", ylab="Density")
}
# save plot
dev.copy(png, "density_trajectories.png")
dev.off()
