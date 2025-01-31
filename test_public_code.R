# Test the public pipeline

library(dplyr)
library(purrr)
library(foreach)
# source all the files
# assuming working directory is the one with all this code
other_files <- setdiff(list.files("."), "test_public_code.R")
walk(other_files, source)
# Simulate assignments data

indications <- data.frame(
  indication = c("Breast", "CLL", "Colorectal", "DLBCL", "Gastric", "Head and Neck", "Hepatocellular Carcinoma", "Melanoma",
                 "Multiple Myeloma", "Non-Small Cell", "Ovarian", "Pancreatic", "Prostate", "Renal", "Small Cell", "Urothelial"),
  proportion = c(0.170, 0.004, 0.181, 0.005, 0.061, 0.014, 0.007, 0.026, 0.012, 0.254, 0.084, 0.068, 0.044, 0.022, 0.016, 0.032)
  )

n_patients <- 30000

set.seed(1)

assignments = data.frame(
  patientid = 1:n_patients,
  indication = sample(indications$indication, size = n_patients, replace = TRUE, prob = indications$proportion),
  assignment = sample(c("Training","Testing"), size = n_patients, replace = TRUE, prob = c(0.8,0.2))
)

age <- sample(18:85, n_patients, replace = TRUE)
gender <- sample(c("Male", "Female"), n_patients, replace = TRUE)
baseline_l1_ecog <- sample(0:4, n_patients, replace = TRUE)
race <- sample(c("White","Nonwhite"), n_patients, replace = TRUE, prob = c(0.73, 0.27))
groupstage <- sample(1:4, n_patients, replace = TRUE, prob = c(0.10, 0.12, 0.24, 0.54))
ml_smoking_status <- sample(c("History of smoking","No history of smoking"), n_patients, replace = TRUE, prob = c(0.58, 0.42))
time_dx_to_l1 <- rnorm(n_patients, mean = 740, sd = 800) # basic assumption
time_dx_to_l1[time_dx_to_l1 < 0] <- 14 # truncate neg values at 14 days: start treatment 2 weeks after diagnosis

# For entry_days_l1_predictor, it's a long right-tailed distribution with lots of 0s. 

# Generate a large number of zeros
zeros <- rep(0, n_patients * 0.8) # 80% zeros

# Generate data for the long tail
set.seed(0)
tail <- rgamma(n_patients * 0.2, shape = 0.5, scale = 30) 

# Combine
entry_days_l1_predictor <- c(zeros, tail)

# Shuffle
set.seed(1)
entry_days_l1_predictor <- sample(entry_days_l1_predictor)

public_data <- data.frame(
  patientid = assignments$patientid,
  indication = assignments$indication,
  mice.imp = 1,
  age = age,
  gender = gender,
  race = race,
  groupstage = groupstage,
  baseline_l1_ecog = baseline_l1_ecog,
  time_dx_to_l1 = time_dx_to_l1,
  entry_days_l1 = entry_days_l1_predictor,
  os_days_l1 = entry_days_l1_predictor + sample.int(365*5, n_patients, replace = TRUE),
  died = sample(c(0,1), n_patients,prob = c(0.3, 0.7), replace = TRUE)
)


# bootstrap_data ----------------------------------------------------------

# Will perform just one for illustration

public_data_boot <- bootstrap_data(1, public_data, assignments)

# make_model_data ---------------------------------------------------------

model_data_sim <- make_model_data(public_data, assignments)
model_outputs <- run_pipeline(config = tibble(algorithm = "cv.glmnet",
                                              features = "benchmark_1",
                                              univariate_screen = FALSE,
                                              indication_method = "factor"),
                              m = 1,
                              row = 1,
                              model_data = model_data_sim,
                              return_fits = TRUE,
                              verbose = 4)
model_data_sim_boot <- make_model_data(public_data_boot, assignments)
model_outputs_boot <- run_pipeline(config = tibble(algorithm = "cv.glmnet",
                                              features = "benchmark_1",
                                              univariate_screen = FALSE,
                                              indication_method = "factor"),
                              m = 1,
                              row = 1,
                              model_data = model_data_sim_boot,
                              return_fits = TRUE,
                              verbose = 4)

