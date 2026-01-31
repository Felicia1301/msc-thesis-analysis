# ------------------------------------------------------------
# 01_Demographics.R
#
# Descriptive statistics and group comparisons for demographic
# and clinical variables.
#
# Groups:
#   HC        = Healthy Controls
#   Distress  = HiTOP Distress disorders
#   Fear      = HiTOP Fear disorders
#
# This script produces all demographic results reported in the
# Methods and Results sections of the thesis.
# ------------------------------------------------------------

# Load required packages
library(here)
library(readxl)
library(dplyr)
library(tableone)
library(gtsummary)
library(car)


# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------

# Note: The data file is not included in this repository due
# to data protection regulations.
Demographics <- read_excel(
  here("data", "Demographics.xlsx"),
  sheet = "All"
)

# ------------------------------------------------------------
# Variable preparation
# ------------------------------------------------------------

# Recode group variable according to HiTOP subdomains
Demographics$Group <- factor(
  Demographics$Group,
  levels = c(0, 1, 2),
  labels = c("HC", "Distress", "Fear")
)

# # Variables included in the demographic table 
vars <- c("Age", "Gender", "BMI", "BDI_2", "STAI_T", "OCI_R")

# Categorical variables
catVars <- c("Gender")



# Table 1 (reported in the manuscript)
Demographics %>%
  select(-ID) %>% 
  tbl_summary(
    by = Group,
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    digits = all_continuous() ~ 2
  ) %>%
  add_overall() %>% 
  as_gt()


# ------------------------------------------------------------
# Assumption checks
# ------------------------------------------------------------

# Convert variables to numeric where required
Demographics$Age <- as.numeric(Demographics$Age)
Demographics$BMI <- as.numeric(Demographics$BMI)
Demographics$STAI_T <- as.numeric(Demographics$STAI_T)
Demographics$BDI_2 <- as.numeric(Demographics$BDI_2)
Demographics$OCI_R <- as.numeric(Demographics$OCI_R)


# Check homogeneity of variances using Levene's test
# For clinical symptom scales with skewed distributions,
# the median-centered (Brown–Forsythe) Levene test was applied (BDI-II, OCI-R).
leveneTest(Age ~ Group, data = Demographics, center = mean) # marginally unequal variances, p = 0.04591
leveneTest(BMI ~ Group, data = Demographics, center = mean) # unequal variances, p = 0.0223
leveneTest(STAI_T ~ Group, data = Demographics) # equal variances
leveneTest(BDI_2 ~ Group, data = Demographics, center = median) # unequal variances, p < 0.001
leveneTest(OCI_R ~ Group, data = Demographics, center = median) # unequal variances, p < 0.05

# ------------------------------------------------------------
# Group comparisons
# ------------------------------------------------------------

# One-way ANOVAs for variables with homogeneous variances

anova_STAI <- aov(STAI_T ~ Group, data = Demographics)
summary (anova_STAI) # F = 132.2; p = <2e-16 ***


# Welch ANOVAs for variables with unequal variances
oneway.test(Age   ~ Group, data = Demographics, var.equal = FALSE) # F = 0.89887, num df = 2.00, denom df = 135.62, p-value = 0.4094
oneway.test(BMI   ~ Group, data = Demographics, var.equal = FALSE) # F = 1.4506, num df = 2.00, denom df = 172.22, p-value = 0.2373
oneway.test(BDI_2 ~ Group, data = Demographics, var.equal = FALSE) # F = 340.74, num df = 2.00, denom df = 217.15, p-value < 2.2e-16
oneway.test(OCI_R ~ Group, data = Demographics, var.equal = FALSE) # F = 30.078, num df = 2.00, denom df = 150.96, p-value = 1.013e-11



# Fisher's exact test used due to small expected cell counts
tbl_gender <- table(Demographics$Gender, Demographics$Group)
fisher.test(tbl_gender)



# Post-hoc pairwise comparisons with FDR correction
t_test_STAI <- pairwise.t.test(
  Demographics$STAI_T, 
  Demographics$Group, p.adjust.method = "fdr",
  pool.sd = TRUE
  )

t_test_BDI <- pairwise.t.test(
  Demographics$BDI_2, 
  Demographics$Group, 
  p.adjust.method = "fdr",
  pool.sd = FALSE
  )

t_test_OCI <- pairwise.t.test(
  Demographics$OCI_R, 
  Demographics$Group, 
  p.adjust.method = "fdr",
  pool.sd = FALSE)

# show results
t_test_STAI
t_test_BDI
t_test_OCI


# Effect sizes:
# Cohen's d (pooled SD) for STAI_T
# d_av (variance-robust) for BDI_2 and OCI_R

compute_cohen_d_pooled <- function(group1, group2) {
  # Remove missing values
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  
  # Return NA if one of the groups is empty
  if (length(group1) == 0 || length(group2) == 0) {
    return(NA)
  }
  
  # Group means
  mean1 <- mean(group1)
  mean2 <- mean(group2)
  
  # Group standard deviations
  sd1 <- sd(group1)
  sd2 <- sd(group2)
  
  # Sample sizes
  n1 <- length(group1)
  n2 <- length(group2)
  
  # Pooled standard deviation
  pooled_sd <- sqrt(((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2))
  
  # Cohen's d
  cohen_d <- (mean1 - mean2) / pooled_sd
  return(cohen_d)
}

# d_av calculated 
compute_d_av <- function(group1, group2) {
  group1 <- na.omit(group1)
  group2 <- na.omit(group2)
  if (length(group1) == 0 || length(group2) == 0) return(NA)
  
  m1 <- mean(group1); 
  m2 <- mean(group2)
  sd1 <- sd(group1); 
  sd2 <- sd(group2)
  
  denom <- sqrt((sd1^2 + sd2^2) / 2)
  return ((m1 - m2) / denom)
}


# Pairwise group contrasts for effect size estimation

variables <- c("STAI_T", "BDI_2", "OCI_R")

comparisons <- list(
  c("HC", "Fear"),
  c("HC", "Distress"),
  c("Fear", "Distress")
)


# Initialize empty data frame for all effect size results
all_results <- data.frame(
  Variable = character(),
  Comparison = character(),
  EffectSizeType = character(),
  EffectSizeValue = numeric(),
  stringsAsFactors = FALSE
)
# Loop over variables and group contrasts
for (var in variables) {
  for (comparison in comparisons) {
    group1 <- Demographics[[var]][Demographics$Group == comparison[1]]
    group2 <- Demographics[[var]][Demographics$Group == comparison[2]]
    
    if (var == "STAI_T") {
      d_val <- compute_cohen_d_pooled(group1, group2)
      es_type <- "Cohen's d (pooled SD)"
    } else {
      d_val <- compute_d_av(group1, group2)
      es_type <- "d_av (variance-robust)"
    }
    
    all_results <- rbind(
      all_results,
      data.frame(
        Variable = var,
        Comparison = paste(comparison[1], "vs.", comparison[2]),
        EffectSizeType = es_type,
        EffectSizeValue = d_val,
        stringsAsFactors = FALSE
      )
    )
  }
}


#Inspect results
print(all_results)
