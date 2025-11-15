install.packages(c("tidyverse", "lme4","emmeans"))
# Load required packages
library(tidyverse)
library(lme4)
library(emmeans)

df <- read.csv("Data_T2.csv", header = TRUE)
str(df)

# Preprocess data
df <- df %>%
  mutate(
    ID = as.factor(ID),
    Treatment = as.factor(Treatment),
    Model = as.factor(Model),
    # create friendly labels
    Treatment_lab = recode(Treatment, `0` = "Control", `1` = "ATXC10X"),
    Model_lab = recode(Model, `0` = "Induced", `1` = "Xenograft")
  )

summary(df$DV)
table(df$Treatment, df$Model)

# Scatter and lines
ggplot(df, aes(x = Time, y = DV, group = ID, color = Treatment_lab)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 1) + 
  facet_wrap(~ Model_lab) +
  scale_y_continuous(trans = "log1p") +
labs(title = "Tumour volume trajectories by subject",
     y = "Tumour volume (mm^3)",
     color = "Treatment") +
  theme_minimal()

# Calculation for SE
df_summary <- df %>%
  group_by(Model_lab, Treatment_lab, Time) %>%
  summarise(
    n = n(),
    mean_DV = mean(DV, na.rm = TRUE),
    sd_DV = sd(DV, na.rm = TRUE),
    se_DV = sd_DV / sqrt(n)
  ) %>% ungroup()

ggplot(df_summary, aes(x = Time, y = mean_DV, color = Treatment_lab)) +
  geom_line() +
  geom_ribbon(aes(ymin = mean_DV - se_DV, ymax = mean_DV + se_DV, fill = Treatment_lab), alpha = 0.2, color = NA) +
  facet_wrap(~ Model_lab, scales = "free_y") +
  labs(title = "Mean tumour volume ± SE by time, treatment and model",
       y = "Mean tumour volume (mm^3)") +
  theme_minimal()

# Histogram(beginning+end)
ggplot(df %>% filter(Time == min(Time)), aes(x = DV, fill = Treatment_lab)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 20) +
  facet_wrap(~ Model_lab, scales = "free_y") +
  labs(title = "Baseline tumour volumes (first timepoint)",
       x = "Tumour volume (mm^3)") +
  theme_minimal()
ggplot(df %>% filter(Time == max(Time)), aes(x = DV, fill = Treatment_lab)) +
  geom_histogram(position = "identity", alpha = 0.6, bins = 20) +
  facet_wrap(~ Model_lab, scales = "free_y") +
  labs(title = "Baseline tumour volumes (last timepoint)",
       x = "Tumour volume (mm^3)") +
  theme_minimal()
#Transform DV
df <- df %>% mutate(logDV = log(DV + 1))
#LME Model
LME_model<- lmer(logDV ~ Time * Treatment * Model + (1 + Time | ID),
                 data = df, REML = FALSE,
                 control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
summary(LME_model)
VarCorr(LME_model)

LME_no3 <- lmer(logDV ~ Time * Treatment + Time * Model + Treatment * Model + (1 + Time | ID),
                data = df, REML = FALSE,
                control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
anova(LME_no3, LME_model)
summary(LME_no3)
VarCorr(LME_no3)
fixef(LME_no3)
confint(LME_no3, method = "Wald")
# Residual plot
par(mfrow = c(2,2))
plot(LME_no3) 
qqnorm(resid(LME_no3)); qqline(resid(LME_no3))
# Random Effects
ranef_df <- as.data.frame(ranef(LME_no3)$ID) %>% tibble::rownames_to_column("ID")
head(ranef_df)


# Predict
time_grid <- seq(min(df$Time, na.rm = TRUE),
                 max(df$Time, na.rm = TRUE),
                 length.out = 50)
emm_grid <- expand.grid(Time = time_grid,
                        Treatment = levels(df$Treatment),
                        Model = levels(df$Model))
emm_grid$Treatment <- factor(emm_grid$Treatment, levels = levels(df$Treatment))
emm_grid$Model <- factor(emm_grid$Model, levels = levels(df$Model))

emm <- emmeans(LME_no3, ~ Time | Treatment * Model, 
               at = list(Time = time_grid), cov.reduce = FALSE)
emm_df <- as.data.frame(emm)
head(emm_df)

emm_df <- as.data.frame(emm) %>%
  mutate(
    pred_logDV = emmean,
    pred_DV = exp(pred_logDV) - 1,
    lower_DV = exp(pred_logDV - 1.96*SE) - 1,  # Wald CI
    upper_DV = exp(pred_logDV + 1.96*SE) - 1,
    Treatment = factor(Treatment),
    Model = factor(Model),
    Treatment_lab = recode(Treatment, `0` = "Control", `1` = "ATXC10X"),
    Model_lab = recode(Model, `0` = "Induced", `1` = "Xenograft")
  )



p_pred <- ggplot() +
  geom_point(data = df, aes(x = Time, y = DV, color = Treatment_lab), alpha = 0.3, size = 0.8) +
  geom_line(data = emm_df, aes(x = Time, y = pred_DV, color = Treatment_lab), size = 1.1) +
  geom_ribbon(data = emm_df, aes(x = Time, ymin = lower_DV, ymax = upper_DV, fill = Treatment_lab), alpha = 0.15) +
  facet_wrap(~ Model_lab, scales = "free_y") +
  labs(title = "Population-level predicted tumour trajectories (back-transformed)",
       y = "Tumour volume (mm^3)",
       x = "Time (days)",
       color = "Treatment",
       fill = "Treatment") +
  theme_minimal()
print(p_pred)



# predict original data and fitted for that subject:
some_id <- levels(df$ID)[2]
df_sub <- df %>% filter(ID == some_id)
df_sub$pred_log <- predict(LME_no3, newdata = df_sub, re.form = NULL)  # includes random effects by default
df_sub$pred_DV <- exp(df_sub$pred_log) - 1

ggplot(df_sub, aes(x = Time)) +
  geom_point(aes(y = DV), color = "black") +
  geom_line(aes(y = pred_DV), color = "red") +
  labs(title = paste("Subject-level fit for ID:", some_id),
       y = "Tumour volume (mm^3)") +
  theme_minimal()





