#### Packages used ####
#install.packages("janitor")
#install.packages("sjPlot")
library(tidyverse)
library(readxl)
library(pROC)
library(modelr)


#### Data availability and import ####
## Data can be downloaded from https://doi.org/10.7910/DVN/9FMMJI

training <- read_excel("data/Oldfield et al. T3cDM biomarker dataset.xlsx",
                       sheet = "Set 2 (training)")
validation <- read_excel("data/Oldfield et al. T3cDM biomarker dataset.xlsx",
                         sheet = "Set 3 (validation)")

## Clean up variable names
### Set 2
training <- training %>%
  rename(
    adiponectin_serum_ug_ml = `Adiponectin Serum  (µg/mL)`,
    il1ra_plasma_pg_ml = `IL-1Ra Plasma (pg/mL)`,
    ca19_9_serum_u_ml_elisa = `CA-19-9   Serum (pg/mL)`,
    aact_serum_ug_ml = `AACT      Sreum  (µg/mL)`
  ) %>%
  janitor::clean_names()

### Set 3
validation <- validation %>%
  rename(
    adiponectin_serum_ug_ml = `Adiponectin    Serum (µg/mL)`,
    il1ra_plasma_pg_ml = `IL-1Ra Plasma (pg/mL)`,
    ca19_9_serum_u_ml_elisa = `CA-19-9   Serum  (U/mL)`
  ) %>%
  janitor::clean_names() %>%
  mutate(simple_group_fig =
           case_when(simple_group == "PDAC" ~ "PDAC",
                     simple_group == "CP" ~ "CP",
                     simple_group == "HC" ~ "HC",
                     main_group == "LSDM" ~ "DM",
                     main_group == "NOD" ~ "NOD"
           ))


#### ROC statistics for all biomarkers in Table 3 ####
## Set 2 (training)
### Create classifier variable based on PDAC versus PDAC-DM, and PDAC-DM vs LSDM

training <- training %>%
  mutate("pdac_pdacdm" = case_when(
    main_group == "PDAC-DM" ~ 0,
    main_group == "PDAC" ~ 1,
    TRUE ~ NA_real_
  ),
  "pdacdm_lsdm" = case_when(
    main_group == "DM" ~ 0,
    main_group == "PDAC-DM" ~ 1,
    TRUE ~ NA_real_
  )
  )

### Long form of dataset for easy nested operations
training_long <- training %>%
  pivot_longer(
    cols = contains(c("serum", "plasma", "vwf")),
    names_to = "biomarker",
    values_to = "concentration"
  )

### Nest by biomarker, then for each comparison calculate AUC, 95% CIs
### and sensitivities at different 0.99, 0.975, 0.95 and 0.9 sensitivity

training_table3_auc <- training_long %>%
  nest(data = -biomarker) %>%
  mutate(
    pdac_pdacdm_auc = map_dbl(data, ~ auc(
      roc(pdac_pdacdm ~ concentration, data = .))),
    pdac_pdacdm_lower_ci = map_dbl(data, ~ ci.auc(
      roc(pdac_pdacdm ~ concentration, data = .))[[1]]
    ),
    pdac_pdacdm_upper_ci = map_dbl(data, ~ ci.auc(
      roc(pdac_pdacdm ~ concentration, data = .))[[3]]
    ),
    pdac_pdacdm_sens_99_spec = map_dbl(data, ~ coords(
      roc(pdac_pdacdm ~ concentration, data = .),
      x = 0.99, input = "specificity")[[3]]
    ),
    pdac_pdacdm_sens_975_spec = map_dbl(data, ~ coords(
      roc(pdac_pdacdm ~ concentration, data = .),
      x = 0.975, input = "specificity")[[3]]
    ),
    pdac_pdacdm_sens_95_spec = map_dbl(data, ~ coords(
      roc(pdac_pdacdm ~ concentration, data = .),
      x = 0.95, input = "specificity")[[3]]
    ),
    pdac_pdacdm_sens_90_spec = map_dbl(data, ~ coords(
      roc(pdac_pdacdm ~ concentration, data = .),
      x = 0.9, input = "specificity")[[3]]
    ),
    pdacdm_lsdm_auc = map_dbl(data, ~ auc(
      roc(pdacdm_lsdm ~ concentration, data = .))),
    pdacdm_lsdm_lower_ci = map_dbl(data, ~ ci.auc(
      roc(pdacdm_lsdm ~ concentration, data = .))[[1]]
    ),
    pdacdm_lsdm_upper_ci = map_dbl(data, ~ ci.auc(
      roc(pdacdm_lsdm ~ concentration, data = .))[[3]]
    ),
    pdacdm_lsdm_sens_99_spec = map_dbl(data, ~ coords(
      roc(pdacdm_lsdm ~ concentration, data = .),
      x = 0.99, input = "specificity")[[3]]
    ),
    pdacdm_lsdm_sens_975_spec = map_dbl(data, ~ coords(
      roc(pdacdm_lsdm ~ concentration, data = .),
      x = 0.975, input = "specificity")[[3]]
    ),
    pdacdm_lsdm_sens_95_spec = map_dbl(data, ~ coords(
      roc(pdacdm_lsdm ~ concentration, data = .),
      x = 0.95, input = "specificity")[[3]]
    ),
    pdacdm_lsdm_sens_90_spec = map_dbl(data, ~ coords(
      roc(pdacdm_lsdm ~ concentration, data = .),
      x = 0.90, input = "specificity")[[3]]
    )
  ) %>%
  select(-data) %>%
  mutate(across(where(is.double), round, digits = 3)
  )

### Export as .csv file
write_csv(training_table3_auc,
          "output/roc statistics/table 3 ROC stats for training set.csv")

## Set 3 (validation)
### Create classifier variable based on PDAC-DM vs NOD

validation <- validation %>%
  mutate("pdacdm_nod" = case_when(
    main_group == "PDAC-DM" ~ 1,
    main_group == "NOD" ~ 0,
    TRUE ~ NA_real_
  )
  )

### Long form of dataset for easy nested operations
validation_long <- validation %>%
  pivot_longer(
    cols = contains(c("serum", "plasma", "vwf")),
    names_to = "biomarker",
    values_to = "concentration"
  )

### Nest by biomarker, then for each comparison calculate AUC, 95% CIs
### and sensitivities at different 0.99, 0.975, 0.95 and 0.9 sensitivity
validation_table3_auc <- validation_long %>%
  nest(data = -biomarker) %>%
  mutate(
    pdacdm_nod_auc = map_dbl(data, ~ auc(
      roc(pdacdm_nod ~ concentration, data = .))),
    pdacdm_nod_auc_lower_ci = map_dbl(data, ~ ci.auc(
      roc(pdacdm_nod ~ concentration, data = .))[[1]]
    ),
    pdacdm_nod_auc_upper_ci = map_dbl(data, ~ ci.auc(
      roc(pdacdm_nod ~ concentration, data = .))[[3]]
    ),
    pdacdm_nod_sens_99_spec = map_dbl(data, ~ coords(
      roc(pdacdm_nod ~ concentration, data = .),
      x = 0.99, input = "specificity")[[3]]
    ),
    pdacdm_nod_sens_975_spec = map_dbl(data, ~ coords(
      roc(pdacdm_nod ~ concentration, data = .),
      x = 0.975, input = "specificity")[[3]]
    ),
    pdacdm_nod_sens_95_spec = map_dbl(data, ~ coords(
      roc(pdacdm_nod ~ concentration, data = .),
      x = 0.95, input = "specificity")[[3]]
    ),
    pdacdm_nod_sens_90_spec = map_dbl(data, ~ coords(
      roc(pdacdm_nod ~ concentration, data = .),
      x = 0.9, input = "specificity")[[3]]
    )
  ) %>%
  select(-data) %>%
  mutate(across(where(is.double), round, digits = 3)
  )

### Export as .csv file
write_csv(validation_table3_auc,
          "output/roc statistics/table 3 ROC stats for validation set.csv")

#### Correlation of adiponectin and age/BMI in Set 3 (validation) within disease group ####

## BMi vs adiponectin Spearman Rank correlation, split by simple group
bmi_vs_adp_cor_simple_group <- validation %>%
  group_by(simple_group_fig) %>%
  summarise(r = cor.test(bmi_kg_m2, adiponectin_serum_ug_ml,
                         method = "spearman")$estimate,
            p_value = cor.test(bmi_kg_m2, adiponectin_serum_ug_ml,
                               method = "spearman")$p.value
  )
#bmi_vs_adp_cor_simple_group

## Age vs adiponectin, split by simple group
age_vs_adp_cor_simple_group <- validation %>%
  group_by(simple_group_fig) %>%
  summarise(r = cor.test(age, adiponectin_serum_ug_ml,
                         method = "spearman")$estimate,
            p_value = cor.test(age, adiponectin_serum_ug_ml,
                               method = "spearman")$p.value
  )
#age_vs_adp_cor_simple_group

#### Assessing the addition of CA19-9, age, sex and BMI to IL1-Ra and adiponectin in Set 3####
#### Type3c vs NOD

## Create classifier variable based on Type3cDM versus NOD

validation <- validation %>%
  mutate("t3c_nod" = case_when(
    main_group == "PDAC-DM" | main_group == "CP-DM" ~ 1,
    main_group == "NOD" ~ 0,
    TRUE ~ NA_real_
  )
  )

## Type3c vs NOD logistic regression - IL1-Ra + adiponectin
t3c_nod_glm_adp_il1ra <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                               il1ra_plasma_pg_ml,
                             data = validation,
                             family = "binomial"
)

# Type3c vs NOD logistic regression - IL1-Ra + adiponectin + CA19-9
t3c_nod_glm_adp_il1ra_ca199 <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                                     il1ra_plasma_pg_ml + ca19_9_serum_u_ml_elisa,
                                   data = validation,
                                   family = "binomial"
)

# Type3c vs NOD logistic regression - IL1-Ra + adiponectin + age
t3c_nod_glm_adp_il1ra_age <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                                     il1ra_plasma_pg_ml + age,
                                   data = validation,
                                   family = "binomial"
)

# Type3c vs NOD logistic regression - IL1-Ra + adiponectin + sex
t3c_nod_glm_adp_il1ra_sex <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                                   il1ra_plasma_pg_ml + gender,
                                 data = validation,
                                 family = "binomial"
)

# Type3c vs NOD logistic regression - IL1-Ra + adiponectin + BMI
t3c_nod_glm_adp_il1ra_bmi <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                                   il1ra_plasma_pg_ml + bmi_kg_m2,
                                 data = validation,
                                 family = "binomial"
)

# Type3c vs NOD logistic regression - IL1-Ra + adiponectin + age + BMI
t3c_nod_glm_adp_il1ra_age_bmi <- glm(t3c_nod ~ adiponectin_serum_ug_ml +
                                   il1ra_plasma_pg_ml + age + bmi_kg_m2,
                                 data = validation,
                                 family = "binomial"
)

## Type3c vs NOD logistic regression - IL1-Ra + adiponectin,
## filtered to include only cases with age and BMI data

t3c_nod_glm_adp_il1ra_filt_bmi_age <- validation %>%
  filter(!is.na(age) & !is.na(bmi_kg_m2)) %>%
  glm(t3c_nod ~ adiponectin_serum_ug_ml + il1ra_plasma_pg_ml,
      data = .,
      family = "binomial"
  )

#summary(t3c_nod_glm_adp_il1ra)
#summary(t3c_nod_glm_adp_il1ra_ca199)
#summary(t3c_nod_glm_adp_il1ra_age)
#summary(t3c_nod_glm_adp_il1ra_sex)
#summary(t3c_nod_glm_adp_il1ra_bmi)
#summary(t3c_nod_glm_adp_il1ra_age_bmi)
#summary(t3c_nod_glm_adp_il1ra_filt_bmi_age)

## Export as formatted tables
sjPlot::tab_model(t3c_nod_glm_adp_il1ra, t3c_nod_glm_adp_il1ra_ca199,
                  file = "output/model estimates/T3c DM vs NOD - IL1-Ra + Adiponectin + CA19-9.html")
sjPlot::tab_model(t3c_nod_glm_adp_il1ra_filt_bmi_age, t3c_nod_glm_adp_il1ra_age_bmi,
                  file = "output/model estimates/Supplementary Table S4, T3c DM vs NOD - IL1-Ra + Adiponectin + age + BMI.html")
sjPlot::tab_model(t3c_nod_glm_adp_il1ra, t3c_nod_glm_adp_il1ra_sex,
                  file = "output/model estimates/T3c DM vs NOD - IL1-Ra + Adiponectin + sex.html")


## Add predicted probabilities for each patient from both models
validation <- validation %>%
  add_predictions(t3c_nod_glm_adp_il1ra, 
                  var = "t3c_nod_adp_il1ra_prob",
                  type = "response") %>%
  add_predictions(t3c_nod_glm_adp_il1ra_ca199, 
                  var = "t3c_nod_adp_il1ra_ca199_prob",
                  type = "response")

## Calculate AUC of ROC curve for T3cDM vs NOD - Adiponectin + IL1-Ra
roc_t3c_nod_adp_il1ra <- validation %>%
  roc(t3c_nod, t3c_nod_adp_il1ra_prob)
auc(roc_t3c_nod_adp_il1ra)
  
## Calculate AUC of ROC curve for T3cDM vs NOD - Adiponectin + IL1-Ra + CA19-9
roc_t3c_nod_adp_il1ra_ca199 <- validation %>%
  roc(t3c_nod, t3c_nod_adp_il1ra_ca199_prob)
auc(roc_t3c_nod_adp_il1ra_ca199)

#### Type3c DM vs all diabetes mellitus

## Create classifier variable based on Type3c DM versus all DM
validation <- validation %>%
  mutate("t3c_dm" = case_when(
    main_group == "PDAC-DM" | main_group == "CP-DM" ~ 1,
    main_group == "NOD" | main_group == "LSDM" ~ 0 ,
    TRUE ~ NA_real_
  )
  )

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin
t3c_dm_glm_adp_il1ra <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                              il1ra_plasma_pg_ml,
                            data = validation,
                            family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin + CA19-9
t3c_dm_glm_adp_il1ra_ca199 <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                                    il1ra_plasma_pg_ml + ca19_9_serum_u_ml_elisa,
                                  data = validation,
                                  family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin + age
t3c_dm_glm_adp_il1ra_age <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                                    il1ra_plasma_pg_ml + age,
                                  data = validation,
                                  family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin + sex
t3c_dm_glm_adp_il1ra_sex <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                                  il1ra_plasma_pg_ml + gender,
                                data = validation,
                                family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin + BMI
t3c_dm_glm_adp_il1ra_bmi <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                                    il1ra_plasma_pg_ml + bmi_kg_m2,
                                  data = validation,
                                  family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin + age + BMI
t3c_dm_glm_adp_il1ra_age_bmi <- glm(t3c_dm ~ adiponectin_serum_ug_ml +
                                  il1ra_plasma_pg_ml + age + bmi_kg_m2,
                                data = validation,
                                family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin,
## filtered to include only cases with age and BMI data

t3c_dm_glm_adp_il1ra_filt_bmi_age <- validation %>%
  filter(!is.na(age) & !is.na(bmi_kg_m2)) %>%
  glm(t3c_dm ~ adiponectin_serum_ug_ml + il1ra_plasma_pg_ml,
      data = .,
      family = "binomial"
      )

#summary(t3c_dm_glm_adp_il1ra)
#summary(t3c_dm_glm_adp_il1ra_ca199)
#summary(t3c_dm_glm_adp_il1ra_age)
#summary(t3c_dm_glm_adp_il1ra_sex)
#summary(t3c_dm_glm_adp_il1ra_bmi)
#summary(t3c_dm_glm_adp_il1ra_age_bmi)
#summary(t3c_dm_glm_adp_il1ra_filt_bmi_age)

## Export as formatted tables
sjPlot::tab_model(t3c_dm_glm_adp_il1ra, t3c_dm_glm_adp_il1ra_ca199,
          file = "output/model estimates/T3c DM vs all DM - IL1-Ra + Adiponectin + CA19-9.html")
          
sjPlot::tab_model(t3c_dm_glm_adp_il1ra_filt_bmi_age, t3c_dm_glm_adp_il1ra_age_bmi,
          file = "output/model estimates/Supplementary Table S3, T3c DM vs all DM - IL1-Ra + Adiponectin + age + BMI.html")

sjPlot::tab_model(t3c_dm_glm_adp_il1ra, t3c_dm_glm_adp_il1ra_sex,
                  file = "output/model estimates/T3c DM vs all DM - IL1-Ra + Adiponectin + sex.html")

## Add predicted probabilities for each patient from IL1-Ra + adiponectin +
## age + bmi models
validation <- validation %>%
  add_predictions(t3c_nod_glm_adp_il1ra_age_bmi, 
                  var = "t3c_nod_adp_il1ra_age_bmi_prob",
                  type = "response") %>%
  add_predictions(t3c_dm_glm_adp_il1ra_age_bmi, 
                  var = "t3c_dm_adp_il1ra_age_bmi_prob",
                  type = "response")

## Calculate AUC of ROC curve for T3cDM vs all DM - 
## adiponectin + IL-1Ra + age + BMI
roc_t3c_dm_adp_il1ra_age_bmi <- validation %>%
  roc(t3c_dm, t3c_dm_adp_il1ra_age_bmi_prob)
#auc(roc_t3c_dm_adp_il1ra_age_bmi)
#ci.auc(roc_t3c_dm_adp_il1ra_age_bmi)

## Calculate AUC of ROC curve for T3cDM vs NOD - 
## adiponectin + IL-1Ra + age + BMI
roc_t3c_nod_adp_il1ra_age_bmi <- validation %>%
  roc(t3c_nod, t3c_nod_adp_il1ra_age_bmi_prob)
#auc(roc_t3c_nod_adp_il1ra_age_bmi)
#ci.auc(roc_t3c_nod_adp_il1ra_age_bmi)

## Add predicted probabilities for each patient from IL1-Ra + adiponectin +
## CA19-9 models
validation <- validation %>%
  add_predictions(t3c_dm_glm_adp_il1ra, 
                  var = "t3c_dm_adp_il1ra_prob",
                  type = "response") %>%
  add_predictions(t3c_dm_glm_adp_il1ra_ca199, 
                  var = "t3c_dm_adp_il1ra_ca199_prob",
                  type = "response")

## Calculate AUC of ROC curve for T3cDM vs all DM - Adiponectin + IL1-Ra
roc_t3c_dm_adp_il1ra <- validation %>%
  roc(t3c_dm, t3c_dm_adp_il1ra_prob)
#auc(roc_t3c_dm_adp_il1ra)

## Calculate AUC of ROC curve for T3cDM vs all DM - Adiponectin + IL1-Ra + CA19-9
roc_t3c_dm_adp_il1ra_ca199 <- validation %>%
  roc(t3c_dm, t3c_dm_adp_il1ra_ca199_prob)
#auc(roc_t3c_dm_adp_il1ra_ca199)

#### Assessment of CA19-9 alone in Type 3c DM vs NOD
## Type3cDM vs NOD
t3c_nod_ca199_roc <- validation %>%
  roc(t3c_nod, ca19_9_serum_u_ml_elisa)

### AUC and 95% CI
#auc(t3c_nod_ca199_roc)
#ci.auc(t3c_nod_ca199_roc)

## Type3c DM vs all DM
t3c_dm_ca199_roc <- validation %>%
  roc(t3c_dm, ca19_9_serum_u_ml_elisa)

### AUC and 95% CI
#auc(t3c_dm_ca199_roc)
#ci.auc(t3c_dm_ca199_roc)

