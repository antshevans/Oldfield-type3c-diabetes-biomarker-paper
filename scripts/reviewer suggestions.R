#### Packages used ####
#install.packages("janitor")
#install.packages("sjPlot")
library(tidyverse)
library(readxl)
library(pROC)
library(modelr)
library(gt)
library(corrplot)

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

### Add detailed labels for each biomarker
training_table3_auc <- training_table3_auc %>%
  mutate(
    biomarker_labs = # Create labels for table headings
      case_when(
        biomarker == "adiponectin_serum_ug_ml" ~ "Adiponectin (\U00B5g/mL)",
        biomarker == "adrenomedullin_plasma_ng_m_l" ~ "Adrenomedullin plasma (ng/mL)",
        biomarker == "adrenomedullin_serum_ng_m_l" ~ "Adrenomedullin serum (ng/mL)",
        biomarker == "aact_serum_ug_ml" ~ "AACT (\U00B5g/mL)",
        biomarker == "apo_a1_serum_mg_d_l" ~ "Apo-A1 (mg/dL)",
        biomarker == "ca19_9_serum_u_ml_elisa" ~ "CA19-9 (U/mL)",
        biomarker == "chemerin_serum_ng_m_l" ~ "Chemerin (ng/mL)",
        biomarker == "clusterin_serum_ng_m_l" ~ "Clusterin (ngmL)",
        biomarker == "c_peptide_serum_pg_m_l" ~ "C-peptide (pg/mL)",
        biomarker == "ghrelin_serum_pg_m_l" ~ "Ghrelin (pg/mL)",
        biomarker == "gip_serum_pg_m_l" ~ "GIP (pg/mL)",
        biomarker == "glp_1_serum_pg_m_l" ~ "GLP-1 (pg/mL)",
        biomarker == "glucagon_serum_pg_m_l" ~ "Glucagon (pg/mL)",
        biomarker == "ifn_g_plasma_pg_m_l" ~ "IFN-G (pg/mL)",
        biomarker == "il_12_plasma_pg_m_l" ~ "IL-12 (pg/mL)",
        biomarker == "il1ra_plasma_pg_ml" ~ "IL-1Ra (pg/mL)",
        biomarker == "il_4_plasma_pg_m_l" ~ "IL-4 (pg/mL)",
        biomarker == "il_6_plasma_pg_m_l" ~ "IL-6 (pg/mL)",
        biomarker == "il_7_plasma_pg_m_l" ~ "IL-7 (pg/mL)",
        biomarker == "il_8_plasma_pg_m_l" ~ "IL-8 (pg/mL)",
        biomarker == "il_9_plasma_pg_m_l" ~ "IL-9 (pg/mL)",
        biomarker == "insulin_serum_pg_m_l" ~ "Insulin (pg/mL)",
        biomarker == "leptin_serum_pg_m_l" ~ "Leptin (pg/mL)",
        biomarker == "mip_1a_plasma_pg_m_l" ~ "MIP-1A (pg/mL)",
        biomarker == "mip_1b_plasma_pg_m_l" ~ "MIP-1B (pg/mL)",
        biomarker == "pai_1_serum_pg_m_l" ~ "PAI-1 (pg/mL)",
        biomarker == "pdgf_bb_plasma_pg_m_l" ~ "PDGF-BB (pg/mL)",
        biomarker == "sparc_serum_ng_m_l" ~ "SPARC (ng/mL)",
        biomarker == "transferin_serum_mg_m_l" ~ "Transferrin (mg/mL)",
        biomarker == "tsp_1_serum_ug_m_l" ~ "TSP-1 (\U00B5g/mL)",
        biomarker == "vwf_m_u_m_l" ~ "vWF (mU/mL)"
      ),
    pdac_pdacdm_auc_ci =
      paste0(round(pdac_pdacdm_auc, 2), "\n(", round(pdac_pdacdm_lower_ci, 2), "-", round(pdac_pdacdm_upper_ci, 2), ")"),
    pdacdm_lsdm_auc_ci =
      paste0(round(pdacdm_lsdm_auc, 2), "\n(", round(pdacdm_lsdm_lower_ci, 2), "-", round(pdacdm_lsdm_upper_ci ,2), ")")
  ) %>%
  select(biomarker, biomarker_labs, everything()) %>%
  arrange(biomarker_labs)

### Create gt html/pdf tables for export to supplementary table
training_table3_auc_gt <- training_table3_auc %>%
  select(biomarker_labs,
         pdacdm_lsdm_auc_ci,
         contains("pdacdm_lsdm_sens")
         ) %>%
gt(rowname_col = "biomarker_labs") %>%
  tab_spanner(
    label = "PDAC-DM vs LSDM",
    columns = contains("pdacdm_lsdm")
  ) %>%
  cols_label(
    pdacdm_lsdm_auc_ci = "AUC (95% CI)",
    pdacdm_lsdm_sens_99_spec = "Sensitivity at 99% specificity",
    pdacdm_lsdm_sens_975_spec = "Sensitivity at 97.5% specificity",
    pdacdm_lsdm_sens_95_spec = "Sensitivity at 95% specificity",
    pdacdm_lsdm_sens_90_spec = "Sensitivity at 90% specificity"
  ) %>%
  fmt_percent(
    columns = contains("_sens"),
    decimals = 1
  ) %>%
  tab_source_note(
    md(
      "**AUC**, Area under the curve; **PDAC-DM**, pancreatic cancer-associated diabetes; **LSDM**, long-standing diabetes (>3yr post-diagnosis of DM)"
    )
  ) %>%
  tab_options(
    table.font.names = c("Helvetica", "Arial"),
    table.border.top.style = "hidden",
    table.border.bottom.style = "hidden",
  ) %>%
  cols_width(
    biomarker_labs ~ px(300),
    contains("AUC") ~ px(150),
    contains("_sens") ~ px(100)
  ) %>%
  tab_style(
    style = "padding-left:10px;padding-right:10px;",
    locations = list(cells_body(), cells_column_labels(), cells_column_spanners())
  )

gtsave(training_table3_auc_gt,
       "output/roc statistics/Supplementary Table S4 - table 3 ROC stats for training set.pdf"
       )


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

### Add detailed labels for each biomarker
validation_table3_auc <-  validation_table3_auc %>%
  mutate(
    biomarker_labs = # Create labels for table headings
      case_when(
        biomarker == "adiponectin_serum_ug_ml" ~ "Adiponectin (\U00B5g/mL)",
        biomarker == "ca19_9_serum_u_ml_elisa" ~ "CA19-9 (U/mL)",
        biomarker == "c_peptide_serum_pg_m_l" ~ "C-peptide (pg/mL)",
        biomarker == "ghrelin_serum_pg_m_l" ~ "Ghrelin (pg/mL)",
        biomarker == "gip_serum_pg_m_l" ~ "GIP (pg/mL)",
        biomarker == "glp_1_serum_pg_m_l" ~ "GLP-1 (pg/mL)",
        biomarker == "glucagon_serum_pg_m_l" ~ "Glucagon (pg/mL)",
        biomarker == "ifn_g_plasma_pg_m_l" ~ "IFN-G (pg/mL)",
        biomarker == "il_12_plasma_pg_m_l" ~ "IL-12 (pg/mL)",
        biomarker == "il1ra_plasma_pg_ml" ~ "IL-1Ra (pg/mL)",
        biomarker == "il_4_plasma_pg_m_l" ~ "IL-4 (pg/mL)",
        biomarker == "il_6_plasma_pg_m_l" ~ "IL-6 (pg/mL)",
        biomarker == "il_8_plasma_pg_m_l" ~ "IL-8 (pg/mL)",
        biomarker == "il_9_plasma_pg_m_l" ~ "IL-9 (pg/mL)",
        biomarker == "insulin_serum_pg_m_l" ~ "Insulin (pg/mL)",
        biomarker == "leptin_serum_pg_m_l" ~ "Leptin (pg/mL)",
        biomarker == "mip_1a_plasma_pg_m_l" ~ "MIP-1A (pg/mL)",
        biomarker == "mip_1b_plasma_pg_m_l" ~ "MIP-1B (pg/mL)",
        biomarker == "pai_1_serum_pg_m_l" ~ "PAI-1 (pg/mL)",
        biomarker == "pdgf_bb_plasma_pg_m_l" ~ "PDGF-BB (pg/mL)",
        biomarker == "rantes_plasma_pg_m_l" ~ "RANTES (pg/mL)",
        biomarker == "tsp_1_serum_mg_m_l" ~ "TSP-1 (\U00B5g/mL)"
      ),
    pdacdm_nod_auc_ci =
      paste0(round(pdacdm_nod_auc, 2), "\n(", round(pdacdm_nod_auc_lower_ci, 2), "-", round(pdacdm_nod_auc_upper_ci, 2), ")")
    ) %>%
  select(biomarker, biomarker_labs, everything()) %>%
  arrange(biomarker_labs)

### Create gt html/pdf tables for export to supplementary table
validation_table3_auc_gt <- validation_table3_auc %>%
  select(
    biomarker_labs, pdacdm_nod_auc_ci, contains("pdacdm_nod_sens")
  ) %>%
  gt(rowname_col = "biomarker_labs") %>%
  tab_spanner(
    label = "PDAC-DM vs NOD",
    columns = contains("pdacdm_nod")
  ) %>%
  cols_label(
    pdacdm_nod_auc_ci = "AUC (95% CI)",
    pdacdm_nod_sens_99_spec = "Sensitivity at 99% specificity",
    pdacdm_nod_sens_975_spec = "Sensitivity at 97.5% specificity",
    pdacdm_nod_sens_95_spec = "Sensitivity at 95% specificity",
    pdacdm_nod_sens_90_spec = "Sensitivity at 90% specificity"
    ) %>%
  fmt_percent(
    columns = contains("_sens"),
    decimals = 1
  ) %>%
  tab_source_note(
    md(
      "**AUC**, Area under the curve; **PDAC-DM**, pancreatic cancer-associated diabetes; **NOD**, new-onset diabetes (<3yr post-diagnosis of DM)"
    )
  ) %>%
  tab_options(
    table.font.names = c("Helvetica", "Arial"),
    table.border.top.style = "hidden",
    table.border.bottom.style = "hidden",
  ) %>%
  cols_width(
    biomarker_labs ~ px(220),
    contains("AUC") ~ px(150),
    contains("_sens") ~ px(120)
  ) %>%
  tab_style(
    style = "padding-left:10px;padding-right:10px;",
    locations = list(cells_body(), cells_column_labels(), cells_column_spanners())
  )

gtsave(validation_table3_auc_gt,
       "output/roc statistics/Supplementary Table S5 - table 3 ROC stats for validation set.pdf"
)

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
                  file = "output/model estimates/Supplementary Table S7, T3c DM vs NOD - IL1-Ra + Adiponectin + age + BMI.html")
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
#auc(roc_t3c_nod_adp_il1ra)
#ci.auc(roc_t3c_nod_adp_il1ra)
  
## Calculate AUC of ROC curve for T3cDM vs NOD - Adiponectin + IL1-Ra + CA19-9
roc_t3c_nod_adp_il1ra_ca199 <- validation %>%
  roc(t3c_nod, t3c_nod_adp_il1ra_ca199_prob)
#auc(roc_t3c_nod_adp_il1ra_ca199)
#ci.auc(roc_t3c_nod_adp_il1ra_ca199)

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
          file = "output/model estimates/Supplementary Table S6, T3c DM vs all DM - IL1-Ra + Adiponectin + age + BMI.html")

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
#ci.auc(roc_t3c_dm_adp_il1ra)

## Calculate AUC of ROC curve for T3cDM vs all DM - Adiponectin + IL1-Ra + CA19-9
roc_t3c_dm_adp_il1ra_ca199 <- validation %>%
  roc(t3c_dm, t3c_dm_adp_il1ra_ca199_prob)
#auc(roc_t3c_dm_adp_il1ra_ca199)
#ci.auc(roc_t3c_dm_adp_il1ra_ca199)


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

#### Correlation matrices of all markers ####

#### Training set
training_spearman <- training %>%
  rename(
    "Adiponectin (\U00B5g/mL)" = "adiponectin_serum_ug_ml",
    "Adrenomedullin plasma (ng/mL)" = "adrenomedullin_plasma_ng_m_l",
    "Adrenomedullin serum (ng/mL)" = "adrenomedullin_serum_ng_m_l",
    "AACT (\U00B5g/mL)" = "aact_serum_ug_ml",
    "Apo-A1 (mg/dL)" = "apo_a1_serum_mg_d_l",
    "CA19-9 (U/mL)" = "ca19_9_serum_u_ml_elisa",
    "Chemerin (ng/mL)" = "chemerin_serum_ng_m_l",
    "Clusterin (ngmL)" = "clusterin_serum_ng_m_l",
    "C-peptide (pg/mL)" = "c_peptide_serum_pg_m_l",
    "Ghrelin (pg/mL)" = "ghrelin_serum_pg_m_l",
    "GIP (pg/mL)" = "gip_serum_pg_m_l",
    "GLP-1 (pg/mL)" = "glp_1_serum_pg_m_l",
    "Glucagon (pg/mL)" = "glucagon_serum_pg_m_l",
    "IFN-G (pg/mL)" = "ifn_g_plasma_pg_m_l",
    "IL-12 (pg/mL)" = "il_12_plasma_pg_m_l",
    "IL-1Ra (pg/mL)" = "il1ra_plasma_pg_ml",
    "IL-4 (pg/mL)" = "il_4_plasma_pg_m_l",
    "IL-6 (pg/mL)" = "il_6_plasma_pg_m_l",
    "IL-7 (pg/mL)" = "il_7_plasma_pg_m_l",
    "IL-8 (pg/mL)" = "il_8_plasma_pg_m_l",
    "IL-9 (pg/mL)" = "il_9_plasma_pg_m_l",
    "Insulin (pg/mL)" = "insulin_serum_pg_m_l",
    "Leptin (pg/mL)" = "leptin_serum_pg_m_l",
    "MIP-1A (pg/mL)" = "mip_1a_plasma_pg_m_l",
    "MIP-1B (pg/mL)" = "mip_1b_plasma_pg_m_l",
    "PAI-1 (pg/mL)" = "pai_1_serum_pg_m_l",
    "PDGF-BB (pg/mL)" = "pdgf_bb_plasma_pg_m_l",
    "SPARC (ng/mL)" = "sparc_serum_ng_m_l",
    "Transferrin (mg/mL)" = "transferin_serum_mg_m_l",
    "vWF (mU/mL)" = "vwf_m_u_m_l"
  ) %>%
  select(contains(c("mL", "dL"))) %>%
  cor(use = "complete.obs",
      method = "spearman")

### Create correlation matrix using Spearman's rank method
pdf(file = "output/figures/supplementary figure 3a - spearman's correlation matrix of all training markers.pdf")
corrplot(training_spearman,
         method = "square",
         type = "lower",
         tl.col = "black",
         tl.cex = 0.7,
         tl.offset = 0.4,
         addCoef.col = "black",
         number.cex = 0.3)
dev.off()

#### Validation set
validation_spearman <- validation %>%
  rename(
    "Adiponectin (\U00B5g/mL)" = "adiponectin_serum_ug_ml",
    "CA19-9 (U/mL)" = "ca19_9_serum_u_ml_elisa",
    "C-peptide (pg/mL)" = "c_peptide_serum_pg_m_l",
    "Ghrelin (pg/mL)" = "ghrelin_serum_pg_m_l",
    "GIP (pg/mL)" = "gip_serum_pg_m_l",
    "GLP-1 (pg/mL)" = "glp_1_serum_pg_m_l",
    "Glucagon (pg/mL)" = "glucagon_serum_pg_m_l",
    "IFN-G (pg/mL)" = "ifn_g_plasma_pg_m_l",
    "IL-12 (pg/mL)" = "il_12_plasma_pg_m_l",
    "IL-1Ra (pg/mL)" = "il1ra_plasma_pg_ml",
    "IL-4 (pg/mL)" = "il_4_plasma_pg_m_l",
    "IL-6 (pg/mL)" = "il_6_plasma_pg_m_l",
    "IL-8 (pg/mL)" = "il_8_plasma_pg_m_l",
    "IL-9 (pg/mL)" = "il_9_plasma_pg_m_l",
    "Insulin (pg/mL)" = "insulin_serum_pg_m_l",
    "Leptin (pg/mL)" = "leptin_serum_pg_m_l",
    "MIP-1A (pg/mL)" = "mip_1a_plasma_pg_m_l",
    "MIP-1B (pg/mL)" = "mip_1b_plasma_pg_m_l",
    "PAI-1 (pg/mL)" = "pai_1_serum_pg_m_l",
    "PDGF-BB (pg/mL)" = "pdgf_bb_plasma_pg_m_l",
    "RANTES (pg/mL)" = "rantes_plasma_pg_m_l"
    ) %>%
  select(contains("mL")) %>%
  cor(use = "complete.obs",
      method = "spearman")

### Create correlation matrix using Spearman's rank method
pdf(file = "output/figures/supplementary figure 3b - spearman's correlation matrix of all validation markers.pdf")
corrplot(validation_spearman,
         method = "square",
         type = "lower",
         tl.col = "black",
         tl.cex = 0.7,
         tl.offset = 0.4,
         addCoef.col = "black",
         number.cex = 0.5)
dev.off()
