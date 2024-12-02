
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

obs_wide_high$obs_increase <- obs_wide$observed.2 / obs_wide$observed.1
obs_wide_med$obs_increase <- obs_wide$observed.4 / obs_wide$observed.3
obs_wide_low$obs_increase <- obs_wide$observed.6 / obs_wide$observed.5

obs_wide_high$sero_convert <- ifelse(obs_wide_high$obs_increase >= log2(2.5), 1, 0) #lower threshold here than for longitudinals
obs_wide_med$sero_convert <- ifelse(obs_wide_med$obs_increase >= log2(2.5), 1, 0)
obs_wide_low$sero_convert <- ifelse(obs_wide_low$obs_increase >= log2(2.5), 1, 0)

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

true_inf_high <- true_inf_high %>% filter(value==1 & x==10)
true_inf_med <- true_inf_med %>% filter(value==1 & x==10)
true_inf_low <- true_inf_low %>% filter(value==1 & x==10)

denominator_high <- nrow(obs_wide_high) 
denominator_med <- nrow(obs_wide_med) 
denominator_low <- nrow(obs_wide_low) 

true_attackrate_high <- length(unique(true_inf_high$i))/denominator_high
true_attackrate_med <- length(unique(true_inf_med$i))/denominator_med
true_attackrate_low <- length(unique(true_inf_low$i))/denominator_low


###Sampled cross-sections, with a subset of longitudinals####

#######high AR##########
infected <- true_inf_high %>% select(i)
infected$infected <- 1
obs_wide_high <- merge(obs_wide_high, infected, by = "i", all.x = T)
obs_wide_high$infected <- ifelse(is.na(obs_wide_high$infected), 0 , obs_wide_high$infected)
obs_wide_high$infected <- ifelse(is.na(obs_wide_high$obs_increase), NA , obs_wide_high$infected)
obs_wide_high$obs_increase_abs <- obs_wide_high$observed.2 - obs_wide_high$observed.1

sero_conv <- obs_wide_high %>% filter(sero_convert==1) 
infected <- obs_wide_high %>% filter(infected==1) 
summary(sero_conv$obs_increase_abs) # observed increase in titers among people who are measured as seroconverters
summary(infected$obs_increase_abs)  # observed increase in titers among true infecteds

true_positives <- obs_wide_high %>% filter(sero_convert==1 & infected==1)
summary(true_positives$obs_increase_abs)
false_positives <- obs_wide_high %>% filter(sero_convert==1 & infected==0)
summary(false_positives$obs_increase_abs)
false_negatives <- obs_wide_high %>% filter(sero_convert==0 & infected==1)
summary(false_negatives$obs_increase_abs)

colors <- c("Observed titer increase among observed sero-converters" = "red",
            "Observed titer increase among true infecteds" = "blue")
ggplot() +
  geom_histogram(data = sero_conv, aes(x = obs_increase_abs, fill = "Observed titer increase among observed sero-converters"),
                 alpha = 0.5) +
  geom_histogram(data = infected, aes(x = obs_increase_abs, fill = "Observed titer increase among true infecteds"),
                 alpha = 0.5) +
  scale_fill_manual(values = c("Observed titer increase among observed sero-converters" = "red",
                               "Observed titer increase among true infecteds" = "blue")) +
  theme(legend.title=element_blank())


sample_loop_longitudinals <- function(db) {
  #n_long <- 500
  n_repeats <- 500
  
  long_ss <- c(100, 500, 1000) #sample size of longitudinals
  ss <- c()
  mean_diff <- c()
  for(b in 1:length(long_ss)){
    for(i in 1:n_repeats) {
      sampled_indices <- sample(nrow(db), size = long_ss[b], replace = F)
      sampled_rows <- db[sampled_indices, ]
      serconverted_longitudinal <- sampled_rows %>% filter(sero_convert == 1)
      mean_obs <- mean(serconverted_longitudinal$obs_increase_abs , na.rm = T)
      print(long_ss[b])
      print(mean_obs)
      ss <- append(ss, long_ss[b])
      mean_diff <- append(mean_diff, mean_obs)
    }
  }
  df_long_ss <- data.frame(matrix(ncol = 2, nrow = n_repeats*3))
  colnames(df_long_ss) <- c('ss', 'mean_diff')
  df_long_ss$ss <- ss
  df_long_ss$mean_diff <- mean_diff
  return(df_long_ss)
}

df_long_ss <- sample_loop_longitudinals(obs_wide_high)

df_long_100 <- df_long_ss %>% filter(ss==100)
df_long_500 <- df_long_ss %>% filter(ss==500)
df_long_1000 <- df_long_ss %>% filter(ss==1000)

obs_t1 <- obs_states %>% filter(time==1) %>% filter(!is.na(observed))
obs_t2 <- obs_states %>% filter(time==2)  %>% filter(!is.na(observed))

sample_loop_crosssectionals <- function(db_t1, db_t2, n_repeats, n_long) {
  df <- data.frame(matrix(ncol = 2))
  colnames(df) <- c('sample_size', 'estimated_diff')
  sample_sizes <- seq(1000, 7000, by = 1000)
  sample_sizes <- c(250, 500, sample_sizes)
  
  for(i in 1:length(sample_sizes)){
    ss <- c()
    for(j in 1:n_repeats) {
      sampled_indices <- sample(nrow(db_t1), size = sample_sizes[i], replace = F)
      sampled_rows_t1 <- db_t1[sampled_indices, ]
      sampled_indices <- sample(nrow(db_t2), size = sample_sizes[i], replace = F)
      sampled_rows_t2 <- db_t2[sampled_indices, ]
      t1_value <- mean(sampled_rows_t1$observed)
      t2_value <- mean(sampled_rows_t2$observed)
      sero <- t2_value-t1_value
      ss <- append(ss, sero)
    }
    df_tobind <- data.frame(matrix(ncol = 2, nrow = n_repeats))
    colnames(df_tobind) <- c('sample_size', 'estimated_diff')
    df_tobind$sample_size <- sample_sizes[i]
    df_tobind$estimated_diff <- ss
    df <- bind_rows(df, df_tobind)
    #df$scaling_factor <- df_long_100$mean_diff[b]
    #df$estimated_conversion <- df$estimated_diff/df_long_100$mean_diff[b]
  }
  
  df <- df %>% filter(!is.na(estimated_diff))
  df <- df %>% slice(rep(1:n(), each = n_long))
  df_long_100 <- df_long_100[rep(row.names(df_long_100), times = n_repeats), ]
  df_long_500 <- df_long_500[rep(row.names(df_long_500), times = n_repeats), ]
  df_long_1000 <- df_long_1000[rep(row.names(df_long_1000), times = n_repeats), ]
  
  df$scaling_factor_100 <- df_long_100$mean_diff
  df$scaling_factor_500 <- df_long_500$mean_diff
  df$scaling_factor_1000 <- df_long_1000$mean_diff
  
  df$estimated_seroconv_100 <- df$estimated_diff/df$scaling_factor_100
  df$estimated_seroconv_500 <-  df$estimated_diff/df$scaling_factor_500
  df$estimated_seroconv_1000 <-  df$estimated_diff/df$scaling_factor_1000
  return(df)
}
  
df <- sample_loop_crosssectionals(db_t1 = obs_t1, db_t2 = obs_t2, n_repeats = 1000, n_long = 500)

create_plot <- function(db, true_attackrate) {
  colors <- c("True attack rate" = "darkgreen",
              "100 longitudinals" = "blue",
              "500 longitudinals" = "red",
              "1000 longitudinals" = "purple")
  
  ggplot() +
    geom_violin(data = db, aes(x = factor(sample_size), y = estimated_seroconv_100, color = "100 longitudinals"), alpha = 0.2, fill = "blue") +
    geom_violin(data = db, aes(x = factor(sample_size), y = estimated_seroconv_500, color = "500 longitudinals"), alpha = 0.2, fill = "red") +
    geom_violin(data = db, aes(x = factor(sample_size), y = estimated_seroconv_1000, color = "1000 longitudinals"), alpha = 0.2, fill = "purple") +
    ylim(0.0, 0.7) +
    labs(y = "Estimated seroconversion rate",
         x = "Sample size",
         color = "") +
    geom_hline(aes(yintercept=true_attackrate, color = "True attack rate"),linetype = "dashed",size=1) +
    scale_color_manual(values = colors)
}

create_plot(db = df, true_attackrate = true_attackrate_high)

estimate_power_absolute <- function(db, margin, attack_rate) {
  db$truth <- attack_rate
  db$diff_100 <- abs(db$estimated_seroconv_100 - db$truth)
  db$diff_500 <- abs(db$estimated_seroconv_500 - db$truth)
  db$diff_1000 <- abs(db$estimated_seroconv_1000 - db$truth)
  
  db$within_margin_100 <- db$diff_100 <= margin
  db$within_margin_500 <- db$diff_500 <= margin
  db$within_margin_1000 <- db$diff_1000 <= margin
  power <- db %>% group_by(sample_size) %>% summarise(power_100 = mean(within_margin_100),
                                                      power_500 = mean(within_margin_500),
                                                      power_1000 = mean(within_margin_1000))
  return(power)
  
}

power_5 <- estimate_power_absolute(db = df, margin = 0.05, attack_rate = true_attackrate_high)
power_1 <- estimate_power_absolute(db = df, margin = 0.1, attack_rate = true_attackrate_high)

colors <- c("Estimate within +/- 5pp of truth" = "red",
            "Estimate within +/- 10pp of truth" = "blue")
sample_sizes <- seq(1000, 7000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_5, aes(x=sample_size, y=power_100, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_500, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 5pp of truth" )) +
  
  geom_line(data =power_1, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10pp of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

estimate_power_relative <- function(db, margin, attack_rate) {
  db$truth <- attack_rate
  db$diff_rel_100 <- db$estimated_seroconv_100 / db$truth
  db$diff_rel_500 <- db$estimated_seroconv_500 / db$truth
  db$diff_rel_1000 <- db$estimated_seroconv_1000 / db$truth
  db$within_margin_100 <- ifelse(db$diff_rel_100 >= (1-margin) & db$diff_rel_100 <= (1+margin), 1, 0)
  db$within_margin_500 <- ifelse(db$diff_rel_500 >= (1-margin) & db$diff_rel_500 <= (1+margin), 1, 0)
  db$within_margin_1000 <- ifelse(db$diff_rel_1000 >= (1-margin) & db$diff_rel_1000 <= (1+margin), 1, 0)
  power <- db %>% group_by(sample_size) %>% summarise(power_100 = mean(within_margin_100),
                                                      power_500 = mean(within_margin_500),
                                                      power_1000 = mean(within_margin_1000))
}


power_10 <- estimate_power_relative(db = df, margin = 0.1, attack_rate = true_attackrate_high)
power_20 <- estimate_power_relative(db = df, margin = 0.2, attack_rate = true_attackrate_high)
power_50 <- estimate_power_relative(db = df, margin = 0.5, attack_rate = true_attackrate_high)



colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10% of truth" )) +
  
  
  geom_line(data =power_20, aes(x=sample_size, y=power_100, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_500, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 20% of truth" )) +
  
  
  geom_line(data =power_50, aes(x=sample_size, y=power_100, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_500, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 50% of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 


#######medium AR##########
infected <- true_inf_med %>% select(i)
infected$infected <- 1
obs_wide_med <- merge(obs_wide_med, infected, by = "i", all.x = T)
obs_wide_med$infected <- ifelse(is.na(obs_wide_med$infected), 0 , obs_wide_med$infected)
obs_wide_med$infected <- ifelse(is.na(obs_wide_med$obs_increase), NA , obs_wide_med$infected)
obs_wide_med$obs_increase_abs <- obs_wide_med$observed.4 - obs_wide_med$observed.3

sero_conv <- obs_wide_med %>% filter(sero_convert==1) 
infected <- obs_wide_med %>% filter(infected==1) 
summary(sero_conv$obs_increase_abs) # observed increase in titers among people who are measured as seroconverters
summary(infected$obs_increase_abs)  # observed increase in titers among true infecteds

true_positives <- obs_wide_med %>% filter(sero_convert==1 & infected==1)
summary(true_positives$obs_increase_abs)
false_positives <- obs_wide_med %>% filter(sero_convert==1 & infected==0)
summary(false_positives$obs_increase_abs)
false_negatives <- obs_wide_med %>% filter(sero_convert==0 & infected==1)
summary(false_negatives$obs_increase_abs)

colors <- c("Observed titer increase among observed sero-converters" = "red",
            "Observed titer increase among true infecteds" = "blue")
ggplot() +
  geom_histogram(data = sero_conv, aes(x = obs_increase_abs, fill = "Observed titer increase among observed sero-converters"),
                 alpha = 0.5) +
  geom_histogram(data = infected, aes(x = obs_increase_abs, fill = "Observed titer increase among true infecteds"),
                 alpha = 0.5) +
  scale_fill_manual(values = c("Observed titer increase among observed sero-converters" = "red",
                               "Observed titer increase among true infecteds" = "blue")) +
  theme(legend.title=element_blank())

df_long_ss <- sample_loop_longitudinals(obs_wide_med)

df_long_100 <- df_long_ss %>% filter(ss==100)
df_long_500 <- df_long_ss %>% filter(ss==500)
df_long_1000 <- df_long_ss %>% filter(ss==1000)

obs_t1 <- obs_states %>% filter(time==3) %>% filter(!is.na(observed))
obs_t2 <- obs_states %>% filter(time==4)  %>% filter(!is.na(observed))

df <- sample_loop_crosssectionals(db_t1 = obs_t1, db_t2 = obs_t2, n_repeats = 1000, n_long = 500)

create_plot(db = df, true_attackrate = true_attackrate_med)

power_5 <- estimate_power_absolute(db = df, margin = 0.05, attack_rate = true_attackrate_med)
power_1 <- estimate_power_absolute(db = df, margin = 0.1, attack_rate = true_attackrate_med)

colors <- c("Estimate within +/- 5pp of truth" = "red",
            "Estimate within +/- 10pp of truth" = "blue")
sample_sizes <- seq(1000, 7000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_5, aes(x=sample_size, y=power_100, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_500, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 5pp of truth" )) +
  
  geom_line(data =power_1, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10pp of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 


power_10 <- estimate_power_relative(db = df, margin = 0.1, attack_rate = true_attackrate_med)
power_20 <- estimate_power_relative(db = df, margin = 0.2, attack_rate = true_attackrate_med)
power_50 <- estimate_power_relative(db = df, margin = 0.5, attack_rate = true_attackrate_med)


colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10% of truth" )) +
  
  
  geom_line(data =power_20, aes(x=sample_size, y=power_100, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_500, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 20% of truth" )) +
  
  
  geom_line(data =power_50, aes(x=sample_size, y=power_100, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_500, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 50% of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 



#######low AR##########
infected <- true_inf_low %>% select(i)
infected$infected <- 1
obs_wide_low <- merge(obs_wide_low, infected, by = "i", all.x = T)
obs_wide_low$infected <- ifelse(is.na(obs_wide_low$infected), 0 , obs_wide_low$infected)
obs_wide_low$infected <- ifelse(is.na(obs_wide_low$obs_increase), NA , obs_wide_low$infected)
obs_wide_low$obs_increase_abs <- obs_wide_low$observed.6 - obs_wide_low$observed.5

sero_conv <- obs_wide_low %>% filter(sero_convert==1) 
infected <- obs_wide_low %>% filter(infected==1) 
summary(sero_conv$obs_increase_abs) # observed increase in titers among people who are measured as seroconverters
summary(infected$obs_increase_abs)  # observed increase in titers among true infecteds

true_positives <- obs_wide_low %>% filter(sero_convert==1 & infected==1)
summary(true_positives$obs_increase_abs)
false_positives <- obs_wide_low %>% filter(sero_convert==1 & infected==0)
summary(false_positives$obs_increase_abs)
false_negatives <- obs_wide_low %>% filter(sero_convert==0 & infected==1)
summary(false_negatives$obs_increase_abs)

colors <- c("Observed titer increase among observed sero-converters" = "red",
            "Observed titer increase among true infecteds" = "blue")
ggplot() +
  geom_histogram(data = sero_conv, aes(x = obs_increase_abs, fill = "Observed titer increase among observed sero-converters"),
                 alpha = 0.5) +
  geom_histogram(data = infected, aes(x = obs_increase_abs, fill = "Observed titer increase among true infecteds"),
                 alpha = 0.5) +
  scale_fill_manual(values = c("Observed titer increase among observed sero-converters" = "red",
                               "Observed titer increase among true infecteds" = "blue")) +
  theme(legend.title=element_blank())

df_long_ss <- sample_loop_longitudinals(obs_wide_low)

df_long_100 <- df_long_ss %>% filter(ss==100)
df_long_500 <- df_long_ss %>% filter(ss==500)
df_long_1000 <- df_long_ss %>% filter(ss==1000)

obs_t1 <- obs_states %>% filter(time==5) %>% filter(!is.na(observed))
obs_t2 <- obs_states %>% filter(time==6)  %>% filter(!is.na(observed))

df <- sample_loop_crosssectionals(db_t1 = obs_t1, db_t2 = obs_t2, n_repeats = 1000, n_long = 500)

create_plot(db = df, true_attackrate = true_attackrate_low)

power_5 <- estimate_power_absolute(db = df, margin = 0.05, attack_rate = true_attackrate_low)
power_1 <- estimate_power_absolute(db = df, margin = 0.1, attack_rate = true_attackrate_low)

colors <- c("Estimate within +/- 5pp of truth" = "red",
            "Estimate within +/- 10pp of truth" = "blue")
sample_sizes <- seq(1000, 7000, by = 1000)
sample_sizes <- c(250, 500, sample_sizes)
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_5, aes(x=sample_size, y=power_100, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_500, color = "Estimate within +/- 5pp of truth" )) +
  geom_line(data =power_5, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 5pp of truth" )) +
  
  geom_line(data =power_1, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10pp of truth" )) +
  geom_line(data =power_1, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10pp of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 


power_10 <- estimate_power_relative(db = df, margin = 0.1, attack_rate = true_attackrate_low)
power_20 <- estimate_power_relative(db = df, margin = 0.2, attack_rate = true_attackrate_low)
power_50 <- estimate_power_relative(db = df, margin = 0.5, attack_rate = true_attackrate_low)


colors <- c("Estimate within +/- 10% of truth" = "darkgreen",
            "Estimate within +/- 20% of truth" = "red",
            "Estimate within +/- 50% of truth" = "blue")
ggplot() +
  geom_hline(aes(yintercept=0.8),
             color = "black", linetype = "dashed") +
  geom_line(data =power_10, aes(x=sample_size, y=power_100, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_500, color = "Estimate within +/- 10% of truth" )) +
  geom_line(data =power_10, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 10% of truth" )) +
  
  
  geom_line(data =power_20, aes(x=sample_size, y=power_100, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_500, color = "Estimate within +/- 20% of truth" )) +
  geom_line(data =power_20, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 20% of truth" )) +
  
  
  geom_line(data =power_50, aes(x=sample_size, y=power_100, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_500, color = "Estimate within +/- 50% of truth" )) +
  geom_line(data =power_50, aes(x=sample_size, y=power_1000, color = "Estimate within +/- 50% of truth" )) +
  
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.1)) +
  scale_x_continuous(limits = c(250, 7000), breaks = sample_sizes) +
  theme_minimal()+
  scale_color_manual(values = colors) +
  labs(y = "Power",
       x = "Sample size",
       color = "") 

