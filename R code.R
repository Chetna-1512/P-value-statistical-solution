# Load required libraries
library(readxl)  
library(dplyr)  
library(tidyr)  
library(ggplot2)  
library(lme4)  
library(lmerTest) 
library(emmeans) 
library(shiny)  
library(DT) 

# Load the dataset
file_path <- "C:/Users/DELL/Downloads/conc_data.xlsx"
data <- read_excel(file_path)

# Inspect the dataset
str(data)
head(data)

# Reshape to long format
colnames(data)
data_long <- data %>%
  pivot_longer(
    cols = matches("^\\d+|\\d+\\.\\d+$"), # Matches numerical column names like 0, 0.33
    names_to = "time",
    values_to = "concentration"
  ) %>%
  mutate(time = as.numeric(time))

# Inspect reshaped data
print(data_long, n=50)

################################ Calculations of parameters
### Cmax
cmax <- data_long %>%
  group_by(sub, TRT) %>% # Group by subject and treatment
  summarise(Cmax = max(concentration, na.rm = TRUE), .groups = "drop") # Calculate Cmax

# View results
print(cmax)

### AUCt
# Function to calculate AUC using Linear Trapezoidal Rule
calculate_auc <- function(time, concentration) {
  # Sort by time (just in case it's not sorted)
  sorted_data <- data.frame(time = time, concentration = concentration) %>%
    arrange(time)
  
  # Apply Linear Trapezoidal Rule
  auc <- sum((diff(sorted_data$time) * 
                (head(sorted_data$concentration, -1) + tail(sorted_data$concentration, -1)) / 2))
  return(auc)
}

# Calculate AUCt for each subject and treatment
auc_data <- data_long %>%
  group_by(sub, TRT) %>% # Group by subject and treatment
  summarise(AUCt = calculate_auc(time, concentration), .groups = "drop") # Calculate AUC

# View results
print(auc_data)

### kel
# Small constant to avoid log(0)
epsilon <- 1e-6

# Adjust concentration
data_long <- data_long %>%
  mutate(adjusted_conc = ifelse(concentration == 0, epsilon, concentration))

# Log transformation of adjusted concentration
data_long <- data_long %>%
  mutate(log_conc = log(adjusted_conc))

# Fit linear model for each subject and treatment group
lambda_z_results <- data_long %>%
  group_by(sub, TRT) %>%  # Group by subject and treatment
  summarise(
    k = abs(coef(lm(log_conc ~ time))[2]),  # k is the absolute value of the slope
    .groups = "drop"
  )

# View results for all 20 subjects with treatments T & R
print(lambda_z_results)

# k values
k = as.numeric(lambda_z_results$k)
print(k)

### AUCi
#C_last
C_last <- data_long %>%
  filter(time == 24) %>%
  select(TRT, concentration)

# Calculate AUC_i
AUC_i <- auc_data$AUCt + (C_last$concentration / k)

# Print the result
print(AUC_i)

### t_half
t_half = log(2)/k

### Tmax 
tmax_results <- data_long %>%
  group_by(sub, TRT) %>%
  summarise(
    Tmax = time[which.max(concentration)], # Time corresponding to the maximum concentration
    Cmax = max(concentration, na.rm = TRUE), # Maximum concentration
    .groups = "drop"
  )

# View Tmax results
print(tmax_results)

result = data.frame(Sub=data$sub, Trmt=data$TRT, Seq = data$seq, Per = data$per, kel = k,Tmax = tmax_results$Tmax, 
                    Cmax = tmax_results$Cmax, AUCt = auc_data$AUCt, AUCi = AUC_i, t_half = t_half)

############################ Visualization of Test vs. Reference Ratios
# Separate Test and Reference data
test_data <- result %>% filter(Trmt == "T") %>% rename_with(~ paste0(., "_T"), -Sub)
ref_data <- result %>% filter(Trmt == "R") %>% rename_with(~ paste0(., "_R"), -Sub)

# Merge Test and Reference data by subject
merged_data <- test_data %>%
  inner_join(ref_data, by = c("Sub" = "Sub"))

# Calculate Ratios
ratios <- merged_data %>%
  mutate(
    Cmax_ratio = Cmax_T / Cmax_R,
    AUCt_ratio = AUCt_T / AUCt_R,
    AUCi_ratio = AUCi_T / AUCi_R
  ) %>%
  select(Sub, Cmax_ratio, AUCt_ratio, AUCi_ratio)  # Keep relevant columns

# View Results
print(ratios)

# Convert ratios data to long format for plotting
ratios_long <- ratios %>%
  pivot_longer(
    cols = c(Cmax_ratio, AUCt_ratio, AUCi_ratio),
    names_to = "Parameter",
    values_to = "Ratio"
  )

# Create the plot
ggplot(ratios_long, aes(x = Sub, y = Ratio, color = Parameter, group = Parameter)) +
  geom_line(size = 1) +                           # Line connecting ratios
  geom_point(size = 2) +                          # Points for individual ratios
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", size = 0.7) +  # Reference line at 0.8
  geom_hline(yintercept = 1.25, linetype = "dashed", color = "red", size = 0.7) + # Reference line at 1.25
  labs(
    title = "Subject vs. Test/Reference Ratios for Cmax, AUCt, and AUCi",
    x = "Subject",
    y = "Test/Reference Ratio",
    color = "Parameter"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels if needed
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

################################ Statistical analysis
# Load necessary libraries
library(lme4)       # For mixed-effects models
library(lmerTest)   # For p-values in mixed-effects models
library(emmeans)    # For LS means
library(dplyr)      # For data manipulation

# Log-transform the PK parameters
result <- result %>%
  mutate(
    log_Cmax = log(Cmax),
    log_AUCt = log(AUCt),
    log_AUCi = log(AUC_i)
  )

# Build mixed-effects models for each PK parameter
# Cmax model
model_Cmax <- lmer(log_Cmax ~ Trmt + Seq + Per + (1|Seq:Sub), data = result)

# AUCt model
model_AUCt <- lmer(log_AUCt ~ Trmt + Seq + Per + (1|Seq:Sub), data = result)

# AUCi model
model_AUCi <- lmer(log_AUCi ~ Trmt + Seq + Per + (1|Seq:Sub), data = result)

# Calculate LS means for treatments
lsmeans_Cmax <- emmeans(model_Cmax, ~ Trmt)
lsmeans_AUCt <- emmeans(model_AUCt, ~ Trmt)
lsmeans_AUCi <- emmeans(model_AUCi, ~ Trmt)

# Perform ANOVA to get p-values
anova_Cmax <- anova(model_Cmax)
anova_AUCt <- anova(model_AUCt)
anova_AUCi <- anova(model_AUCi)

# Display results
# LS Means
print("LS Means for Cmax:")
print(lsmeans_Cmax)

print("LS Means for AUCt:")
print(lsmeans_AUCt)

print("LS Means for AUCi:")
print(lsmeans_AUCi)

# ANOVA Results
print("ANOVA for Cmax:")
print(anova_Cmax)

print("ANOVA for AUCt:")
print(anova_AUCt)

print("ANOVA for AUCi:")
print(anova_AUCi)

