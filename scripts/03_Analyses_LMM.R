# ------------------------------------------------------------
# Script: 03_Analyses_LMM.R
#
# Cardiac feedback processing – statistical analyses
#
# This script contains the final preprocessing, model, and figures
# reported in the thesis.
#
# Raw and processed data are not included due to data
# protection regulations.
# ------------------------------------------------------------


## ============================================================
## 0. Setup
## ============================================================

rm(list = ls())

library(dplyr)
library(lme4)
library(lmerTest)
library(emmeans)
library(performance)
library(ggplot2)
library(effects)
library(showtext)


data_dir <- "data"

## ============================================================
## 1. Data Import
## ============================================================

# Load newest final_table_clean automatically
files <- list.files(
  "Outputs",
  pattern = "final_table_clean_.*\\.csv$", 
  full.names = TRUE
)
latest_file <- files[which.max(file.info(files)$mtime)]

df <- read.csv(latest_file)
show (df)

cat("Loaded file:", latest_file, "\n")

## ============================================================
## 2. Preprocessing
## ============================================================

df <- df %>%
  rename(
    ID = CleanID,
    Task = task,
    IBI = ibi_name,
    IBI_change = ibi_change,
    Valence = feedback
  ) %>%
  mutate(
    Group = factor(Group, levels = c("HC", "Distress", "Fear")),
    Valence = factor(Valence, levels = c("neg", "pos")),
    Task = factor(Task, levels = c("Doors", "Reversal")),
    Choice = factor(Choice, levels = c("neg_stay", "neg_switch", "pos_stay", "pos_switch")),
    Gender = case_when(
      Gender %in% c(0, "0") ~ "male",
      Gender %in% c(1, "1") ~ "female",
      Gender %in% c(2, "2") ~ "diverse",
      TRUE ~ as.character(Gender)
    ),
    Gender = factor(Gender, levels = c("male", "female", "diverse")),
    ID = factor(ID),
    IBI = factor(IBI, levels = c("-2", "-1", "0", "1", "2"))
  ) %>%
  mutate(
    IBI_change = as.numeric(IBI_change),
    Age_c = as.numeric(Age) - mean(as.numeric(Age), na.rm = TRUE),
    BMI_c = as.numeric(BMI) - mean(as.numeric(BMI), na.rm = TRUE),
  )

# Sum contrasts 
contrasts(df$Group)   <- contr.sum(length(levels(df$Group)))
contrasts(df$Valence) <- contr.sum(length(levels(df$Valence)))
contrasts(df$Task)    <- contr.sum(length(levels(df$Task)))
contrasts(df$IBI)     <- contr.sum(length(levels(df$IBI)))

contrasts(df$Choice)  <- contr.sum(length(levels(df$Choice)))
contrasts(df$Gender)  <- contr.sum(length(levels(df$Gender)))



# ------------------------------------------------------------
# Covariate evaluation
# ------------------------------------------------------------

covariates <- c("Age_c", "BMI_c", "Gender")

for (cov in covariates) {
  
  formula <- as.formula(
    paste0("IBI_change ~ Group * Valence * Task * IBI + ", cov, " + (1 | ID)")
  )
  
  model <- lmer(formula, data = df)
  
  cat("\n\n========================\n")
  cat("Model with covariate:", cov, "\n")
  cat("========================\n\n")
  
  print(summary(model))
}


## ============================================================
## 4. FINAL MODEL (reported)
## ============================================================

final_model <- lmer (
  IBI_change ~ Group * Valence * Task * IBI + Age_c +
    (1 + Valence*Task | ID),
  control = lmerControl (optimizer = "bobyqa", 
                         optCtrl = list (maxfun = 1e6)), 
  data = df
)

didLmerConverge(final_model)
isSingular(final_model)
summary(final_model)
anova (final_model, type = 3)
r2 (final_model)
icc (final_model)

## ============================================================
## 5. Model Assumptions
## ============================================================
# Normally distributed residuals
qqnorm(resid(final_model))
qqline(resid(final_model))

#Homoscedasticity of residuals
plot(final_model, resid(.) ~ fitted(.))

#Influential observations / outliers
library(performance)
check_outliers(final_model)

## ============================================================
## 6. Estimated Marginal Means and Contrasts (reported results)
## ============================================================
emm_group_ibi        <- emmeans(final_model, ~ Group | IBI)
emm_valence_ibi  <- emmeans(final_model, ~ Valence | IBI)
emm_task_ibi     <- emmeans(final_model, ~ Task | IBI)
emm_group_task_ibi <- emmeans(final_model, ~ Group | Task * IBI)
emm_group_val_ibi         <- emmeans(final_model, ~ Group | Valence * IBI)
emm_gvi         <- emmeans(final_model, ~ Group * Valence * IBI)
emm_group_val_task_ibi <- emmeans (final_model, ~ Group | Valence * Task * IBI)
emm_gvti         <- emmeans(final_model, ~ Group * Valence * Task * IBI)



## ------------------------------------------------------------
## Valence effects across time (Valence × IBI)
## ------------------------------------------------------------

summary(emm_valence_ibi)

contrast(
  emm_valence_ibi,
  method = "pairwise",
  adjust = "fdr"
)

## ------------------------------------------------------------
## Task effects across time (Task × IBI)
## ------------------------------------------------------------

summary(emm_task_ibi)

contrast(
  emm_task_ibi,
  method = "pairwise",
  adjust = "fdr"
)

## ------------------------------------------------------------
## Group × Task × IBI (reported interaction)
## ------------------------------------------------------------

summary(emm_group_task_ibi)

contrast(
  emm_group_task_ibi,
  method = "pairwise",
  adjust = "fdr"
)

## ------------------------------------------------------------
## H1 / RQ1: Group differences across time (Group × IBI)
## ------------------------------------------------------------
summary(emm_group_ibi)

contrast(
  emm_group_ibi,
  method = "pairwise",
  adjust = "fdr"
)

## ------------------------------------------------------------
## H2a/H2b: Valence effects across time in Controls 
## ------------------------------------------------------------
summary (emm_group_val_ibi)

contrast(
  emm_group_val_ibi,
  method = "pairwise",
  by = c("Group", "IBI"),
)

## ------------------------------------------------------------
## RQ2: Group differences in valence-specific responses
##       (Group × Valence × IBI)
## ------------------------------------------------------------

summary(emm_group_val_ibi)

contrast(
  emm_group_val_ibi,
  method = "pairwise",
  adjust = "fdr"
)

## ------------------------------------------------------------
## RQ3: Group × Valence × Task × IBI
## ------------------------------------------------------------
summary(emm_group_val_task_ibi )

contrast(
  emm_group_val_task_ibi ,
  method = "pairwise",
  adjust = "fdr"
)

## ============================================================
## 7. Figures
## ============================================================
prep_emm_df <- function(emm_obj) {
  df <- as.data.frame(emm_obj)
  
  if ("Group" %in% names(df)) {
    df <- df %>%
      mutate(Group = recode(Group, "HC" = "Controls"))
  }
  
  if ("Valence" %in% names(df)) {
    df <- df %>%
      mutate(Valence = factor(Valence, levels = c("pos", "neg")))
  }
  
  if ("Task" %in% names(df)) {
    df <- df %>%
      mutate(Task = factor(Task, levels = c("Doors", "Reversal")))
  }
  
  df
}


plot_groups_over_ibi <- function(df, facet_var = NULL) {
  
  p <- ggplot(
    df,
    aes(
      x = factor(IBI),
      y = emmean,
      color = Group,
      group = Group,
      shape = Group
    )
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(
      aes(ymin = asymp.LCL, ymax = asymp.UCL),
      width = 0.1
    ) +
    scale_color_manual(
      values = c(
        "Controls" = "#0072B2",  # blau
        "Distress" = "#E69F00",  # orange
        "Fear"     = "#B50D0D"   # rot
      )
    ) +
    scale_shape_manual(
      values = c(
        "Controls" = 16,
        "Distress" = 17,
        "Fear"     = 15
      )
    ) +
    scale_y_reverse(
      breaks = scales::breaks_width(5)
    ) +
    labs(
      x = "IBI",
      y = "Estimated IBI change (ms)",
      color = "Group"
    ) +
    guides(
      shape = "none",
      color = guide_legend(
        override.aes = list(shape = c(16, 17, 15))
      )
    ) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "bottom",
      axis.title = element_text(size = 16),
      axis.text  = element_text(size = 16),
      legend.title = element_text(size = 16),
      legend.text  = element_text(size = 16)
    )
  
  # Facet NUR wenn explizit gewünscht
  if (!is.null(facet_var)) {
    p <- p +
      facet_wrap(vars(.data[[facet_var]])) +
      theme(
        strip.text = element_text(size = 16, face = "bold"),
        strip.background = element_blank()
      )
  }
  
  return(p)
}


# Figure 1: Valence × IBI
df_fig1 <- prep_emm_df(emm_valence_ibi)

ggplot(
  df_fig1,
  aes(x = factor(IBI), y = emmean, color = Valence, group = Valence)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1) +
  scale_color_manual(
    values = c("pos" = "#16CC9A", "neg" = "#e7298a"),
    labels = c("Positive feedback", "Negative feedback")
  ) +
  scale_y_reverse(
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(x = "IBI", y = "Estimated IBI change (ms)", color = "Feedback valence") +
  theme_classic(base_size = 12)+
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    panel.grid = element_blank(),
  )



# Figure 2: Task × IBI
df_fig2 <- prep_emm_df(emm_task_ibi)

ggplot(
  df_fig2,
  aes(x = factor(IBI), y = emmean, color = Task, group = Task)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1) +
  scale_color_manual(
    values = c("Doors" = "#56B4E9", "Reversal" = "#D55E00")
  ) +
  scale_y_reverse( 
    breaks = scales::pretty_breaks(n = 6)
  )+
  labs(x = "IBI", y = "Estimated IBI change (ms)", color = "Task") +
  theme_classic(base_size = 12)+
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    panel.grid = element_blank(),
  )
# Figure 3: Group × Task × IBI
df_fig3 <- prep_emm_df(emm_group_task_ibi)

plot_groups_over_ibi(df_fig3, facet_var = "Task")

# Figure 4: Group × IBI
df_fig4 <- prep_emm_df(emm_group_ibi)

plot_groups_over_ibi(df_fig4)

# Figure 5: Group × Valence × IBI
df_fig5 <- prep_emm_df(emm_group_val_ibi)

plot_groups_over_ibi(df_fig5, facet_var = "Valence")

# Figure 6: Group × Valence × Task × IBI
df_fig6 <- prep_emm_df(emm_gvti) %>%
  filter(Task %in% c("Doors", "Reversal"))

y_min <- floor(min(df_fig6$asymp.LCL, na.rm = TRUE) / 5) * 5
y_max <- ceiling(max(df_fig6$asymp.UCL, na.rm = TRUE) / 5) * 5

y_min <- floor(min(df_fig6$asymp.LCL, na.rm = TRUE))
y_max <- ceiling(max(df_fig6$asymp.UCL, na.rm = TRUE))

y_scale_fig6 <- scale_y_reverse(
  limits = c(y_max, y_min),
  breaks = pretty(c(y_min, y_max), n = 6)
)


plot_groups_over_ibi(
  subset(df_fig6, Task == "Doors"),
  facet_var = "Valence"
) +
  labs(title = "Doors Task") +
  y_scale_fig6 +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

plot_groups_over_ibi(
  subset(df_fig6, Task == "Reversal"),
  facet_var = "Valence"
) +
  labs(title = "Reversal Learning Task") +
  y_scale_fig6 +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
  )

# Figure 7: Valence × Task × IBI
emm_vti         <- emmeans(final_model, ~  Valence * Task * IBI)
df_fig7 <- prep_emm_df(emm_vti)

ggplot(
  df_fig7,
  aes(
    x = factor(IBI),
    y = emmean,
    color = Valence,
    group = Valence
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = asymp.LCL, ymax = asymp.UCL),
    width = 0.1
  ) +
  facet_wrap(~ Task) +
  scale_color_manual(
    values = c("pos" = "#16CC9A", "neg" = "#e7298a"),
    labels = c("Positive feedback", "Negative feedback")
  ) +
  scale_y_reverse(
    breaks = scales::pretty_breaks(n = 4)
  ) +
  labs(
    x = "IBI",
    y = "Estimated IBI change (ms)",
    color = "Feedback valence"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    legend.position = "bottom",
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 16),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 16),
    panel.grid   = element_blank()
  )

