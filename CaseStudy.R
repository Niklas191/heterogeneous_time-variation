# ==============================================================================
# Capturing heterogeneous time-variation in covariate effects in non-proportional hazard regression models
# Niklas Hagemann, Thomas Kneib and Kathrin Moellenhoff
#
# Code for the case study 
# ==============================================================================

# Packages =====================================================================

library(dplyr) # package for data manipulation
library(pammtools) # package for PEM data transformation
library(mgcv) # package for fitting GAMs
library(survival) # package for handling survival datasets
library(pec) # package for IBS

# Dataset 1 ====================================================================

data1 <- read.delim("difg_tcga_gdc_clinical_data.tsv", stringsAsFactors = TRUE)

data1 <- data1 %>% dplyr::filter(! is.na(Overall.Survival..Months.), ! is.na(Overall.Survival.Status))
data1$event <- as.numeric(data1$Overall.Survival.Status =="1:DECEASED")
data1$time <- data1$Overall.Survival..Months.
names(data1)[c(4)] <- c("Age")
data1 <- data1 %>% filter(Sample.Type == "Primary solid Tumor")
data1 <- data1 %>% select(Patient.ID, Age, Fraction.Genome.Altered, Primary.Diagnosis, Sex, event, time)

# Dataset 2 ====================================================================

data2 <- read.delim("gbm_tcga_gdc_clinical_data.tsv", stringsAsFactors = TRUE)

data2 <- data2 %>% dplyr::filter(! is.na(Overall.Survival..Months.), ! is.na(Overall.Survival.Status))
data2$event <- as.numeric(data2$Overall.Survival.Status =="1:DECEASED")
data2$time <- data2$Overall.Survival..Months.
names(data2)[c(4)] <- c("Age")
data2 <- data2 %>% filter(Sample.Type == "Primary solid Tumor")
data2 <- data2 %>% select(Patient.ID, Age, Fraction.Genome.Altered, Primary.Diagnosis, Sex, event, time)


data <- bind_rows(data1, data2)

# Remove patients with missing values and survival times <0 ====================
data <- data %>% filter(time >= 0)
data <- data %>% na.omit()

# Set survival time to half a day for patients who died on their admission date
data$time[data$time == 0] <- 1/60

# Add end of study after 8 years ===============================================

data$event[data$time >96] <- 0
data$time[data$time >96] <- 96

# Transform data for PAMM ======================================================

data_ped <- data %>% as_ped(Surv(time = time, event = event)~., cut = unique(data$time), id = "id")

# Fit models ===================================================================

set.seed(123456)
model_fs <- pamm(formula = ped_status ~ 
                   s(tend,bs = "ps", k = 12, m = c(3, 1)) +  
                   s(tend, Primary.Diagnosis, bs = "fs", xt  = list(bs = "ps"), k = 12, m = c(3, 1), by = Fraction.Genome.Altered) +
                   Age +
                   Primary.Diagnosis + 
                   Sex, 
                 data = data_ped,
                 family = poisson(), 
                 offset = offset, 
                 method = "REML")

set.seed(123456)
model_add <- pamm(formula = ped_status ~ 
                    s(tend, bs = "ps", k = 12, m = c(3, 1)) +  
                    s(tend, bs = "ps", k = 12, m = c(3, 1), by = Fraction.Genome.Altered) +
                    s(Fraction.Genome.Altered, Primary.Diagnosis, bs = "re") +
                    Age +
                    Primary.Diagnosis + 
                    Sex, 
                  data = data_ped,
                  family = poisson(), 
                  offset = offset, 
                  method = "REML")


set.seed(123456)
model_re <- pamm(formula = ped_status ~ 
                   s(tend,bs = "ps", k = 12, m = c(3, 1)) +  
                   s(Fraction.Genome.Altered, Primary.Diagnosis, bs = "re") +
                   Age +
                   Primary.Diagnosis + 
                   Sex, 
                 data = data_ped,
                 family = poisson(), 
                 offset = offset, 
                 method = "REML")

set.seed(123456)
model_ps <- pamm(formula = ped_status ~ 
                   s(tend, bs = "ps", k = 12, m = c(3, 1)) +  
                   s(tend, bs = "ps", k = 12, m = c(3, 1), by = Fraction.Genome.Altered) +
                   Age +
                   Primary.Diagnosis + 
                   Sex, 
                 data = data_ped,
                 family = poisson(), 
                 offset = offset, 
                 method = "REML")

set.seed(123456)
model_linear <- pamm(formula = ped_status ~ 
                       s(tend, bs = "ps", k = 12, m = c(3, 1)) +  
                       Fraction.Genome.Altered +
                       Age +
                       Primary.Diagnosis + 
                       Sex, 
                     data = data_ped,
                     family = poisson(), 
                     offset = offset, 
                     method = "REML")

# Calculate fit measures =======================================================

pec_fs = pec(model_fs, 
             Surv(time, event) ~ 1, 
             data = data, 
             times = sort(unique(data$time)), 
             start = min(unique(data$time)), 
             exact = FALSE,
             verbose = FALSE)

pec_add = pec(model_add, 
              Surv(time, event) ~ 1, 
              data = data, 
              times = sort(unique(data$time)),
              start = min(unique(data$time)), 
              exact = FALSE,
              verbose = FALSE)

pec_ps = pec(model_ps, 
             Surv(time, event) ~ 1, 
             data = data, 
             times = sort(unique(data$time)), 
             start = min(unique(data$time)), 
             exact = FALSE,
             verbose = FALSE)

pec_re = pec(model_re, 
             Surv(time, event) ~ 1, 
             data = data, 
             times = sort(unique(data$time)), 
             start = min(unique(data$time)), 
             exact = FALSE,
             verbose = FALSE)

pec_lin = pec(model_linear, 
              Surv(time, event) ~ 1, 
              data = data, 
              times = sort(unique(data$time)), 
              start = min(unique(data$time)), 
              exact = FALSE,
              verbose = FALSE)

fit_measures <- cbind(c("FS", "ADD", "RE", "PS", "LIN"), 
                      round(c(logLik(model_fs), logLik(model_add), logLik(model_re), logLik(model_ps), logLik(model_linear)), 2),
                      round(c(
                        pec::crps(pec_fs)[2,1], 
                        pec::crps(pec_add)[2,1], 
                        pec::crps(pec_re)[2,1], 
                        pec::crps(pec_ps)[2,1], 
                        pec::crps(pec_lin)[2,1]), 4),
                      round(c(AIC(model_fs), AIC(model_add), AIC(model_re), AIC(model_ps), AIC(model_linear)), 2)
)

print(fit_measures)
