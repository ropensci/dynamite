library(TraMineR)
library(dplyr)
library(tidyr)
library(btvcm)
data(mvad, package = "TraMineR")
d <- pivot_longer(mvad, 15:86, "time") %>% select(id, time, value)
d$time <- rep(1:72, length = nrow(d))
#d <-  d %>% filter(id < 100)
fit <- btvcm:::btvcmfit(
    obs(value ~ 1, family = categorical()) + lags() +
        splines(df = 10),
    d, id, time, chains = 1, refresh = 10)

cf <- coef(fit) %>% group_by(time, variable) %>%
    summarise(mean = mean(value),
        lwr = quantile(value, 0.025), upr = quantile(value, 0.975))
cf %>%
    ggplot(aes(time, mean)) + theme_bw() +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.7) +
    geom_line() +
    facet_wrap(~ variable, scales = "free_y")

fits <- fitted(fit)

sumr <- fits %>% filter(id == 1) %>% group_by(time) %>%
    summarise(across(starts_with("value"),
        list(mean = mean,
            ymin = ~ quantile(., 0.025, na.rm = TRUE),
            ymax = ~ quantile(., 0.975, na.rm = TRUE))))

sumr %>% ggplot(aes(time, value_training_mean)) +
    geom_ribbon(aes(ymin = value_training_ymin, ymax = value_training_ymax), alpha = 0.3) +
    geom_line() +
    theme_bw()

sumr %>% ggplot(aes(time, value_joblessness_mean)) +
    geom_ribbon(aes(ymin = value_joblessness_ymin, ymax = value_joblessness_ymax), alpha = 0.3) +
    geom_line() +
    theme_bw()

sumr %>% ggplot(aes(time, value_employment_mean)) +
    geom_ribbon(aes(ymin = value_employment_ymin, ymax = value_employment_ymax), alpha = 0.3) +
    geom_line() +
    theme_bw()

seqdplot(seqdef(mvad, 15:86)[1,])



newdata <- d
newdata$value[newdata$time > 1] <- NA
system.time(pred <- predict(fit, newdata = newdata))

sumr <- pred %>% group_by(time) %>%
    summarise(data.frame(proportions(table(value))))

# plot
seqdata <- seqdef(mvad, 15:86)
seqdplot(seqdata)

sumr$value <- factor(sumr$value, levels = rev(alphabet(seqdata)), ordered = TRUE)
cpal <- attributes(seqdata)$cpal
ggplot(sumr, aes(fill = value, y = Freq, x = time)) +
    geom_bar(stat = "identity",
        width = 1, colour = "black") +
    scale_fill_manual(values = rev(cpal),
        labels = rev(alphabet(seqdata))) +
    labs(x = "", y = "probability") +
    theme_bw()

library(seqHMM)
# simple time-homogenous markov model
model <- build_mm(seqdata) # automatically estimates transition probabilities
sim <- simulate_hmm(n_sequences = 1000, model$initial_probs, model$transition_probs, model$emission_probs, 72)
seqdplot(sim$obs)


newdata <- d
newdata$value[newdata$time > 1] <- NA
newdata$value[newdata$time == 1] <- "school" # fix everyone to school first
system.time(pred <- predict(fit, newdata = newdata))

sumr <- pred %>% group_by(time) %>%
    summarise(data.frame(proportions(table(value))))
sumr$value <- factor(sumr$value, levels = rev(alphabet(seqdata)), ordered = TRUE)
cpal <- attributes(seqdata)$cpal
ggplot(sumr, aes(fill = value, y = Freq, x = time)) +
    geom_bar(stat = "identity",
        width = 1, colour = "black") +
    scale_fill_manual(values = rev(cpal),
        labels = rev(alphabet(seqdata))) +
    labs(x = "", y = "probability") +
    theme_bw()

sim <- simulate_hmm(n_sequences = 1000, c(0,0,0,0,1,0), model$transition_probs, model$emission_probs, 72)
seqdplot(sim$obs)

