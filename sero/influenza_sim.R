rm(list=ls())


library(serosim)
library(tidyverse)
library(patchwork)
library(viridis)

## Specify the number of time periods to simulate 
times <- seq(1,4.5*10*12,by=1) 
simulation_settings <- list("t_start"=1,"t_end"=max(times))

## Generate the population demography tibble
N<-5000
## See help file(?generate_pop_demography) for more information on function arguments.
## age_min is set to 0 which allows births to occur until the last time step
## Let's assume that no individuals are removed from the population and set prob_removal to 0
demography <- generate_pop_demography(N, times, age_min=0, prob_removal=0.2)

n_strains <- 10

## Create a matrix giving cross-reactivity between each strain pairing
x_coord <- numeric(n_strains)
y_coord <- numeric(n_strains)

## Simulate antigenic drift as a random walk in two dimensions
x_coord[1] <- 0
y_coord[1] <- 0
## Each strain is antigenically drifted from the previous one with log-normally distributed 
## moves through 2-dimensional antigenic space
for(i in 2:n_strains){
  x_coord[i] <- x_coord[i-1] + rlnorm(1, mean=log(1), sd=0.5)
  y_coord[i] <- y_coord[i-1] + rlnorm(1, mean=log(1), sd=0.5)
}

## Cross-reactivity is a function of antigenic distance of the form e^-b*x, where x is the 
## euclidean distance between strains and b is a parameter
antigenic_map <- data.frame(x_coord=x_coord, y_coord=y_coord)
antigenic_map <- as.matrix(dist(antigenic_map, diag=TRUE, upper=TRUE))
##############
antigenic_map <- exp(-0.25*antigenic_map)
##############

antigenic_map <- antigenic_map %>% as_tibble() %>% 
  mutate(exposure_id=1:n()) %>% 
  pivot_longer(-exposure_id) %>% rename(biomarker_id=name) %>% 
  mutate(biomarker_id = as.numeric(biomarker_id))


## Seasonal forcing with a peak force of infection of 0.2 and minimum of 0
## If we assume an infectious period of 5 days, the simulation gives an annually oscillating R0 between 0 and 2 
beta <- 0.2
sigma <- 1
seasonality <- beta*(1 + sigma*cos(2*pi*(times-1)/12))

plot(1:540, seasonality, type = "l")

## Create an empty array to store the force of exposure for all exposure types
foe_pars <- array(0, dim=c(1,max(times),n_strains))
t_start <- 1

## Set the circulation duration for each strain
for(strain in seq_len(n_strains)){
  t_end <- t_start + 4.5*12 - 1
  foe_pars[,t_start:t_end,strain] <- seasonality[t_start:t_end]
  t_start <- t_end + 1
}

## Specify a simple exposure model which calculates the probability of exposure directly from the force of exposure at that time step. In this selected model, the probability of exposure is 1-exp(-FOE) where FOE is the force of exposure at that time.
exposure_model<-exposure_model_simple_FOE

plot_exposure_model(exposure_model=exposure_model_simple_FOE, times=times,
                    n_groups = 1, n_exposures = 10, foe_pars=foe_pars) +
  facet_wrap(~exposure_id)


## Specify immunity model within the runserosim function below 
immunity_model<-immunity_model_ifxn_biomarker_prot

## Define antibody-mediated protection parameters
model_pars_immunity <- tibble(name=c("biomarker_prot_midpoint","biomarker_prot_width"),
                              mean=c(2.8,1.3),sd=c(NA,NA),distribution=c(NA,NA))
model_pars_immunity <- expand_grid(biomarker_id=1:n_strains, model_pars_immunity) %>% 
  mutate(exposure_id=paste0("Strain_",biomarker_id), biomarker_id=paste0("Biomarker_",biomarker_id))

## Parameters chosen from Coudeville et al. 2010
plot_biomarker_mediated_protection(seq(0,8,by=0.1), biomarker_prot_midpoint=2.8, biomarker_prot_width=1.3) 

## Specify the antibody model 
antibody_model<-antibody_model_monophasic_cross_reactivity

## Antibody model parameters
model_pars_kinetics <- tibble(name=c("boost","wane"),mean=c(4,0.01),sd=c(1,0.01),distribution=c("log-normal","log-normal"))
model_pars_kinetics <- expand_grid(exposure_id =seq_len(n_strains),model_pars_kinetics)  %>%
  mutate(biomarker_id=paste0("Biomarker_", exposure_id), exposure_id = paste0("Strain_",exposure_id))
model_pars_observation <- expand_grid(biomarker_id=paste0("Biomarker_",seq_len(n_strains)),name="obs_sd",
                                      mean=NA,sd=1,distribution="normal") %>% mutate(exposure_id=NA)
model_pars_original <- bind_rows(model_pars_kinetics,model_pars_immunity, model_pars_observation) %>% arrange(exposure_id, biomarker_id)

## Bring in the antibody parameters needed for the antibody model
## Note that the observation error parameter needed for the observation model (Section 7) is defined here too.

## Reformat model_pars for runserosim
model_pars<-reformat_biomarker_map(model_pars_original)
head(model_pars)

## Specify the draw_parameters function
draw_parameters<-draw_parameters_random_fx

# plot_antibody_model(antibody_model_monophasic_cross_reactivity, N=10,
#                     model_pars=model_pars,draw_parameters_fn = draw_parameters_random_fx,
#                     biomarker_map=antigenic_map)

## Specify the observation model 
observation_model<-observation_model_discrete_noise

## Specify the cutoffs for the discrete assay. This will depend on the dilutions of the haemagglutination assay.
## This is a matrix with each row containing all of the cutoffs for that biomarker. Here, we have set the same cutoffs for all biomarkers 
breaks <- seq(0,8,by=1)
cutoffs <- matrix(breaks,nrow=n_strains,ncol=length(breaks), byrow=TRUE)

## Specify assay sensitivity and specificity needed for the observation model
sensitivity<-0.95
specificity<-0.99

sample_month_1 <- 488
sample_month_2 <- 488+6
sample_month_3 <- 500
sample_month_4 <- 500+6
sample_month_5 <- 530
sample_month_6 <- 530+6

## Specify observation_times (serological survey sampling design) to observe all biomarkers across all individuals at the end of the simulation
obs1 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_1, 3))), b=1:n_strains)
obs2 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_2, 3))), b=1:n_strains)
obs3 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_3, 3))), b=1:n_strains)
obs4 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_4, 3))), b=1:n_strains)
obs5 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_5, 3))), b=1:n_strains)
obs6 <- expand_grid(tibble(i=1:N, t = floor(rnorm(N, sample_month_6, 3))), b=1:n_strains)

obs1$t <- sample_month_1
obs2$t <- sample_month_2
obs3$t <- sample_month_3
obs4$t <- sample_month_4
obs5$t <- sample_month_5
obs6$t <- sample_month_6

observation_times <- bind_rows(obs1, obs2, obs3, obs4, obs5, obs6)

## Run the core simulation
res<- runserosim(
  simulation_settings=simulation_settings,
  demography=demography,
  observation_times=observation_times,
  foe_pars=foe_pars,
  biomarker_map=antigenic_map,
  model_pars=model_pars,
  exposure_model=exposure_model_simple_FOE,
  immunity_model=immunity_model_ifxn_biomarker_prot,
  antibody_model=antibody_model_monophasic_cross_reactivity,
  observation_model=observation_model_discrete_noise,
  draw_parameters=draw_parameters_random_fx,
  
  ## Other arguments needed
  cutoffs=cutoffs,
  sensitivity=sensitivity,
  specificity=specificity,
  VERBOSE=NULL,
  attempt_precomputation = TRUE,
  parallel=TRUE,
  n_cores=60
)






save.image(here('/simulation_10k.RData'))
