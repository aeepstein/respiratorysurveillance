library(here)
library(dplyr)
library(ggplot2)
library(caret)

load(here('simulation_10k.RData'))

times_to_retain <- observation_times %>% filter(b==10) #keep only the final strain
ids <- unique(times_to_retain$i)  #individuals around while final strain was circulating
true_inf <-  res$immune_histories_long %>% filter(x==10) #true infection status

true_inf_high <- true_inf %>% filter(t>= sample_month_1 & t<= sample_month_2) #true infection status during high season
true_inf_med <- true_inf %>% filter(t>= sample_month_3 & t<= sample_month_4) #true infection status during med season
true_inf_low <- true_inf %>% filter(t>= sample_month_5 & t<= sample_month_6) #true infection status during low season

obs_states <- res$observed_biomarker_states 
obs_states <- obs_states %>% filter(b==10) #observed test results of 10th strain
table(obs_states$i)

time <- rep(1:6, N)
obs_states$time <- time

obs_wide <- reshape(obs_states, idvar = "i", timevar = "time", direction = "wide")

obs_wide_high <- obs_wide %>% select(1, ends_with(".1"), ends_with(".2")) #high season
obs_wide_med <- obs_wide %>% select(1, ends_with(".3"), ends_with(".4")) #med season
obs_wide_low <- obs_wide %>% select(1, ends_with(".5"), ends_with(".6")) #low season

obs_wide_high$obs_increase <- obs_wide$observed.2 / obs_wide$observed.1 #observed fold-increase (log2)
obs_wide_med$obs_increase <- obs_wide$observed.4 / obs_wide$observed.3
obs_wide_low$obs_increase <- obs_wide$observed.6 / obs_wide$observed.5

obs_wide_high$sero_convert <- ifelse(obs_wide_high$obs_increase >= 2, 1, 0) #measured seroconversion
obs_wide_med$sero_convert <- ifelse(obs_wide_med$obs_increase >= 2, 1, 0)
obs_wide_low$sero_convert <- ifelse(obs_wide_low$obs_increase >= 2, 1, 0)

obs_wide_high$sero_convert <- ifelse(is.nan(obs_wide_high$obs_increase), 0, obs_wide_high$sero_convert)
obs_wide_med$sero_convert <- ifelse(is.nan(obs_wide_med$obs_increase), 0, obs_wide_med$sero_convert)
obs_wide_low$sero_convert <- ifelse(is.nan(obs_wide_low$obs_increase), 0, obs_wide_low$sero_convert)

obs_wide_high <- obs_wide_high %>% filter(!is.na(obs_increase))
obs_wide_med <- obs_wide_med %>% filter(!is.na(obs_increase))
obs_wide_low <- obs_wide_low %>% filter(!is.na(obs_increase))

#Measured seroconversion rates among entire population
prop.table(table(obs_wide_high$sero_convert))
prop.table(table(obs_wide_med$sero_convert))
prop.table(table(obs_wide_low$sero_convert))

true_inf_high <- true_inf_high %>% filter(value==1 & x==10) #individuals truly infected during each season
true_inf_med <- true_inf_med %>% filter(value==1 & x==10)
true_inf_low <- true_inf_low %>% filter(value==1 & x==10)

denominator_high <- nrow(obs_wide_high) 
denominator_med <- nrow(obs_wide_med) 
denominator_low <- nrow(obs_wide_low) 

true_attackrate_high <- length(unique(true_inf_high$i))/denominator_high
true_attackrate_med <- length(unique(true_inf_med$i))/denominator_med
true_attackrate_low <- length(unique(true_inf_low$i))/denominator_low


#############SAMPLING##############
########high AR###########
obs_wide_high$true_inf <- ifelse(obs_wide_high$i %in% true_inf_high$i, 1, 0) 
table(obs_wide_high$sero_convert, obs_wide_high$true_inf )

sensitivity <- table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,2] /
  (table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,2] +
     table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,2])

specificity <- table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,1] /
  (table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,1] +
     table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,1])

ppv <- table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,2] /
  (table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,2] +
     table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[2,1])

npv <- table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,1] /
  (table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,1] +
     table(obs_wide_high$sero_convert, obs_wide_high$true_inf)[1,2])

sample_loop <- function(db, true_attackrate) {
  df <- data.frame(matrix(ncol = 2, nrow = 2))
  colnames(df) <- c('cat', 'sero_conv')
  df$cat <- c("truth", "measured_seroconvert")
  df$sero_conv <- c(true_attackrate,
                    prop.table(table(db$sero_convert ))[[2]])

  df <- data.frame(matrix(ncol = 2))
  colnames(df) <- c('cat', 'sero_conv')
  
  sample_sizes <- seq(1000, 7000, by = 1000)
  sample_sizes <- c(250, 500, sample_sizes)
  for(i in 1:length(sample_sizes)){
    ss <- c()
    for(j in 1:1000) {
      sampled_indices <- sample(nrow(db), size = sample_sizes[i], replace = F)
      sampled_rows <- db[sampled_indices, ]
      sero <- prop.table(table(sampled_rows$sero_convert))[[2]]
      ss <- append(ss, sero)
    }
    df_tobind <- data.frame(matrix(ncol = 2, nrow = 1000))
    colnames(df_tobind) <- c('cat', 'sero_conv')
    df_tobind$cat <- sample_sizes[i]
    df_tobind$sero_conv <- ss
    df <- bind_rows(df, df_tobind)
  }
  df <- df %>% filter(!is.na(sero_conv))
  return(df)
}

df <- sample_loop(db = obs_wide_high, true_attackrate = true_attackrate_high)

plot_estimates <- function(df, true_attackrate) {
  colors <- c( "True attack rate" = "goldenrod")
  ggplot() + geom_violin(data = df, aes(x = factor(cat), y = sero_conv), fill = "blue", alpha = 0.4, color = "blue") +
    ylim(0.0, 0.5) +
    labs(y = "Estimated seroconversion rate",
         x = "Sample size",
         color = "") +
    geom_hline(aes(yintercept=true_attackrate, color = "True attack rate"),linetype = "dashed",linewidth=1) +
    scale_color_manual(values = colors)
}

plot_estimates(df = df, true_attackrate = true_attackrate_high)

df$truth <- true_attackrate_high
df$diff_abs <- abs(df$sero_conv - df$truth)
df$diff_rel <- df$sero_conv / df$truth

df$within_margin <- df$diff_abs <= 0.01
power_1 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.025
power_25 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.05
power_5 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

sample_sizes <- seq(1000, 9000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)

colors <- c("Estimate within +/- 1pp of truth" = "darkgreen",
            "Estimate within +/- 2.5pp of truth" = "red",
            "Estimate within +/- 5pp of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_point(data =power_1, aes(x=cat, y=power), color = "darkgreen" ) +
  geom_line(data =power_1, aes(x=cat, y=power, color = "Estimate within +/- 1pp of truth" )) +
  geom_line(data =power_25, aes(x=cat, y=power, color = "Estimate within +/- 2.5pp of truth" )) +
  geom_point(data =power_25, aes(x=cat, y=power), color = "red") +
  geom_line(data =power_5, aes(x=cat, y=power, color = "Estimate within +/- 5pp of truth" )) +
  geom_point(data =power_5, aes(x=cat, y=power), color = "blue") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

df$within_margin <- ifelse(df$diff_rel >= 0.9 & df$diff_rel <= 1.1, 1, 0)
power_10 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.80 & df$diff_rel <= 1.2, 1, 0)
power_20 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.5 & df$diff_rel <= 1.5, 1, 0)
power_50 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=cat, y=power, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_20, aes(x=cat, y=power, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_50, aes(x=cat, y=power, color = "Estimate within +/- 50% of truth" )) +
  geom_point(data =power_20, aes(x=cat, y=power), color = "red") +
  geom_point(data =power_50, aes(x=cat, y=power), color = "blue") +
  geom_point(data =power_10, aes(x=cat, y=power), color = "darkgreen" ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 





#########med AR##########
obs_wide_med$true_inf <- ifelse(obs_wide_med$i %in% true_inf_med$i, 1, 0) 
table(obs_wide_med$sero_convert, obs_wide_med$true_inf )

sensitivity <- table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,2] /
  (table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,2] +
     table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,2])

specificity <- table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,1] /
  (table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,1] +
     table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,1])

ppv <- table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,2] /
  (table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,2] +
     table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[2,1])

npv <- table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,1] /
  (table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,1] +
     table(obs_wide_med$sero_convert, obs_wide_med$true_inf)[1,2])

df <- sample_loop(db = obs_wide_med, true_attackrate = true_attackrate_med)

plot_estimates(df = df, true_attackrate = true_attackrate_med)
df$truth <- true_attackrate_med
df$diff_abs <- abs(df$sero_conv - df$truth)
df$diff_rel <- df$sero_conv / df$truth

df$within_margin <- df$diff_abs <= 0.01
power_1 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.025
power_25 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.05
power_5 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

sample_sizes <- seq(1000, 9000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)

colors <- c("Estimate within +/- 1pp of truth" = "darkgreen",
            "Estimate within +/- 2.5pp of truth" = "red",
            "Estimate within +/- 5pp of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_point(data =power_1, aes(x=cat, y=power), color = "darkgreen" ) +
  geom_line(data =power_1, aes(x=cat, y=power, color = "Estimate within +/- 1pp of truth" )) +
  geom_line(data =power_25, aes(x=cat, y=power, color = "Estimate within +/- 2.5pp of truth" )) +
  geom_point(data =power_25, aes(x=cat, y=power), color = "red") +
  geom_line(data =power_5, aes(x=cat, y=power, color = "Estimate within +/- 5pp of truth" )) +
  geom_point(data =power_5, aes(x=cat, y=power), color = "blue") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

df$within_margin <- ifelse(df$diff_rel >= 0.9 & df$diff_rel <= 1.1, 1, 0)
power_10 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.80 & df$diff_rel <= 1.2, 1, 0)
power_20 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.5 & df$diff_rel <= 1.5, 1, 0)
power_50 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=cat, y=power, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_20, aes(x=cat, y=power, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_50, aes(x=cat, y=power, color = "Estimate within +/- 50% of truth" )) +
  geom_point(data =power_20, aes(x=cat, y=power), color = "red") +
  geom_point(data =power_50, aes(x=cat, y=power), color = "blue") +
  geom_point(data =power_10, aes(x=cat, y=power), color = "darkgreen" ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 




#########low AR##########
obs_wide_low$true_inf <- ifelse(obs_wide_low$i %in% true_inf_low$i, 1, 0) 
table(obs_wide_low$sero_convert, obs_wide_low$true_inf )

sensitivity <- table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,2] /
  (table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,2] +
     table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,2])

specificity <- table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,1] /
  (table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,1] +
     table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,1])

ppv <- table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,2] /
  (table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,2] +
     table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[2,1])

npv <- table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,1] /
  (table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,1] +
     table(obs_wide_low$sero_convert, obs_wide_low$true_inf)[1,2])

df <- sample_loop(db = obs_wide_low, true_attackrate = true_attackrate_low)

plot_estimates(df = df, true_attackrate = true_attackrate_low)
df$truth <- true_attackrate_low
df$diff_abs <- abs(df$sero_conv - df$truth)
df$diff_rel <- df$sero_conv / df$truth

df$within_margin <- df$diff_abs <= 0.01
power_1 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.025
power_25 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- df$diff_abs <= 0.05
power_5 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

sample_sizes <- seq(1000, 9000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)

colors <- c("Estimate within +/- 1pp of truth" = "darkgreen",
            "Estimate within +/- 2.5pp of truth" = "red",
            "Estimate within +/- 5pp of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_point(data =power_1, aes(x=cat, y=power), color = "darkgreen" ) +
  geom_line(data =power_1, aes(x=cat, y=power, color = "Estimate within +/- 1pp of truth" )) +
  geom_line(data =power_25, aes(x=cat, y=power, color = "Estimate within +/- 2.5pp of truth" )) +
  geom_point(data =power_25, aes(x=cat, y=power), color = "red") +
  geom_line(data =power_5, aes(x=cat, y=power, color = "Estimate within +/- 5pp of truth" )) +
  geom_point(data =power_5, aes(x=cat, y=power), color = "blue") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

df$within_margin <- ifelse(df$diff_rel >= 0.9 & df$diff_rel <= 1.1, 1, 0)
power_10 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.80 & df$diff_rel <= 1.2, 1, 0)
power_20 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))
df$within_margin <- ifelse(df$diff_rel >= 0.5 & df$diff_rel <= 1.5, 1, 0)
power_50 <- df %>% group_by(cat) %>% summarise(power = mean(within_margin))

colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")

ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=cat, y=power, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_20, aes(x=cat, y=power, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_50, aes(x=cat, y=power, color = "Estimate within +/- 50% of truth" )) +
  geom_point(data =power_20, aes(x=cat, y=power), color = "red") +
  geom_point(data =power_50, aes(x=cat, y=power), color = "blue") +
  geom_point(data =power_10, aes(x=cat, y=power), color = "darkgreen" ) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

