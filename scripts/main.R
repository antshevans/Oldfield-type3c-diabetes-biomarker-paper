#### Packages used ####
#install.packages("janitor")
library(data.table)
setDTthreads(0) # Uses all logical CPUs available
library(readxl)
library(ggplot2)
library(ggbeeswarm)
library(cowplot)
library(purrr)
library(patchwork)
library(modelr)
library(pROC)
library(plotROC)

#### Data import and availability ####
## Data can be downloaded from https://doi.org/10.7910/DVN/9FMMJI

## Import data
set2 <- janitor::clean_names(
  setDT(
    read_excel("data/Oldfield et al. T3cDM biomarker dataset.xlsx",
               sheet = "Set 2 (training)")
  )
)

set3 <- janitor::clean_names(
  setDT(
    read_excel("data/Oldfield et al. T3cDM biomarker dataset.xlsx",
               sheet = "Set 3 (validation)")
  )
)

## Combine both data frames
d <- rbindlist(list(set2, set3), fill = TRUE)

## Correct units within cleaned biomarker names,
## from "m" for "micro" to "u" for "mu"
setnames(d,
         old = c(
           "adiponectin_serum_mg_m_l",
           "transferin_serum_mg_m_l",
           "tsp_1_serum_mg_m_l",
           "aact_sreum_mg_m_l"
         ),
         new = c(
           "adiponectin_serum_ug_m_l",
           "transferin_serum_ug_m_l",
           "tsp_1_serum_ug_m_l",
           "aact_serum_ug_m_l"
         )
)

## Change biomarker labels for easier fitting on plot later
d[,
  plot_group := fcase(
    main_group == "PDAC-DM", "PDAC\nDM",
    main_group == "PDAC", "PDAC",
    main_group == "CP DM", "CP\nDM",
    main_group == "CP", "CP",
    main_group == "DM", "LSDM",
    main_group == "HC", "HC",
    main_group == "CP-DM", "CP\nDM",
    main_group == "LSDM", "LSDM",
    main_group == "NOD", "NOD"
  )]

## Create disease groups for Figure 2D
d[,
     plot_group_scatter := factor(
       fcase(
       main_group == "PDAC-DM", "PDAC",
       main_group == "PDAC", "PDAC",
       main_group == "CP DM", "CP",
       main_group == "CP", "CP",
       main_group == "DM", "DM",
       main_group == "HC", "HC",
       main_group == "CP-DM", "CP",
       main_group == "LSDM", "DM",
       main_group == "NOD", "NOD"
     ),
     levels = c(
       "PDAC", "CP", "NOD", "DM", "HC"
     )
     )
     ]

## Reorder disease groups for supplementary figures
d[,
  plot_group_supp := factor(plot_group,
                            levels = c(
                              "PDAC", "PDAC\nDM",
                              "LSDM", "NOD",
                              "CP", "CP\nDM", 
                              "HC"
                            )
  )
]

## Convert concentrations to reasonable units
d[,
  `:=`(ghrelin_serum_ng_m_l = ghrelin_serum_pg_m_l / 1000,
       clusterin_serum_ug_m_l = clusterin_serum_ng_m_l / 1000,
       c_peptide_serum_ng_m_l = c_peptide_serum_pg_m_l / 1000,
       gip_serum_ng_m_l = gip_serum_pg_m_l / 1000,
       glp_1_serum_ng_m_l = glp_1_serum_pg_m_l / 1000,
       glucagon_serum_ng_m_l = glucagon_serum_pg_m_l / 1000,
       il_1ra_plasma_ng_m_l = il_1ra_plasma_pg_m_l / 1000,
       insulin_serum_ng_m_l = insulin_serum_pg_m_l / 1000,
       leptin_serum_ng_m_l = leptin_serum_pg_m_l / 1000,
       pai_1_serum_ng_m_l = pai_1_serum_pg_m_l / 1000,
       pdgf_bb_plasma_ng_m_l = pdgf_bb_plasma_pg_m_l / 1000,
       sparc_serum_ug_m_l = sparc_serum_ng_m_l / 1000,
       transferin_serum_mg_m_l = transferin_serum_ug_m_l / 1000,
       vwf_u_m_l = vwf_m_u_m_l / 1000,
       rantes_plasma_ng_m_l = rantes_plasma_pg_m_l / 1000
  )
]

#### Figure 2D ####
## Spearman's rank correlation of BMI vs adiponectin in Set 3 (validation)
## Note that BMI data were unavailable for Set 2
bmi_adipo_all_cor <- round(cor(d$bmi_kg_m2,
                               d$adiponectin_serum_ug_m_l,
                               method = "spearman",
                               use = "complete.obs"
), 2)

## Colour palette for disease groups
colour_palette <- c("#E69F00", "#F0E442", "#56B4E9", "#0072B2",  "#009E73")

## Scatterplot of BMI versus ADP in Set 3 (validation)
bmi_adipo_simple_scatter <- d[sample_set == "Validation", ] |>
  ggplot(
    aes(bmi_kg_m2,
        adiponectin_serum_ug_m_l,
        colour = plot_group_scatter,
        shape = plot_group_scatter)
  ) +
  geom_point(size = 1.5) +
  scale_colour_manual(values = colour_palette) +
  scale_shape_manual(values=c(15, 19, 18, 17, 20)) +
  theme_cowplot() +
  theme(legend.position = c(0.7, 0.75)) +
  coord_cartesian(xlim = c(15, 50), y = c(0, 50), expand = FALSE) +
  labs(colour = NULL,
       shape = NULL,
       x = bquote("BMI (kg/m"^2*")"),
       y = "Serum adiponectin (µg/mL)"
  ) +
  annotate("text",
           label = paste("r[s]", "==", bmi_adipo_all_cor),
           parse = TRUE, x = 40, y = 20, hjust = 0)

#bmi_adipo_simple_scatter

save_plot("output/figures/Figure 2D - BMI vs ADP scatterplot.pdf",
          bmi_adipo_simple_scatter,
          base_aspect_ratio = 1.1
)

#### Figures 4C and 4D ####

## Create classifier variable based on Type3cDM versus NOD
## and Type3cDM vs all DM

d[,
  `:=`(
    "t3c_nod" = fcase(
      main_group == "PDAC-DM" | main_group == "CP-DM", 1,
      main_group == "NOD", 0
    ),
    "t3c_dm" = fcase(
      main_group == "PDAC-DM" | main_group == "CP-DM", 1,
      main_group == "NOD" | main_group == "LSDM", 0
    )
  )
  ]

## Type3c vs NOD logistic regression - IL1-Ra + adiponectin
t3c_nod_glm_adp_il1ra <- glm(t3c_nod ~ adiponectin_serum_ug_m_l +
                               il_1ra_plasma_pg_m_l,
                             data = d[sample_set == "Validation"],
                             family = "binomial"
)

## Type3c vs all DM logistic regression - IL1-Ra + adiponectin
t3c_dm_glm_adp_il1ra <- glm(t3c_dm ~ adiponectin_serum_ug_m_l +
                              il_1ra_plasma_pg_m_l,
                            data = d[sample_set == "Validation"],
                            family = "binomial"
)

## Add predicted probabilities for each patient from both models
 d <- d |> add_predictions(t3c_nod_glm_adp_il1ra,
                           var = "t3c_nod_adp_il1ra_prob",
                           type = "response"
                           ) |> 
   add_predictions(t3c_dm_glm_adp_il1ra,
                   var = "t3c_dm_adp_il1ra_prob",
                   type = "response"
                   )

 #### Figure 4C ####
 ## ROC curves of T3Cdm vs all DM with adiponectin + IL-1Ra,
 ## in combination and individually
 t3c_dm_melt <- melt_roc(d[sample_set == "Validation", ],
                          d = "t3c_dm",
                          m = c("adiponectin_serum_ug_m_l",
                                "il_1ra_plasma_pg_m_l",
                                "t3c_dm_adp_il1ra_prob"))
 
 roc_t3c_dm <- ggplot(t3c_dm_melt,
                       aes(d = D, m = M, colour = name, linetype = name)) +
   geom_roc(labels = FALSE,
            n.cuts = 0) +
   scale_y_continuous(expand = c(0,0), limits = c(0, 1.)) +
   scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
   theme_cowplot(font_family = "Helvetica") +
   theme(legend.position = c(0.3, 0.3),
         legend.text = element_text(size = 12),
         legend.key.width=unit(1.5,"cm")
   ) +
   labs(linetype = NULL,
        colour = NULL,
        x = "1 - Specificity",
        y = "Sensitivity") +
   scale_colour_manual(labels = c("Adiponectin (0.76)", "IL-1Ra (0.79)", "Adiponectin + IL-1Ra (0.90)"),
                       values = c("#F8766D", "#00BA38", "Blue")) +
   scale_linetype_manual(labels = c("Adiponectin (0.76)", "IL-1Ra (0.79)", "Adiponectin + IL-1Ra (0.90)"),
                         values = c("solid", "longdash", "solid")
   )
 
## Export figure 
 save_plot(
   "output/figures/Figure 4C - IL-1Ra + Adiponectin ROC curves, typ3c vs all DM.pdf",
   roc_t3c_dm,
   base_asp = 1.1
 )
 
## View ROC curves
#roc_t3c_dm    
 
## Calculate AUCs plus 95% confidence intervals by DeLong method 
#calc_auc(roc_t3c_dm) 

### Using pROC package to calculate 95% CIs
proc_t3c_dm_adp <- roc(t3c_dm ~ adiponectin_serum_ug_m_l,
                              data = d[sample_set == "Validation", ]) 
proc_t3c_dm_il1ra <- roc(t3c_dm ~ il_1ra_plasma_pg_m_l,
                              data = d[sample_set == "Validation", ]) 
proc_t3c_dm_adp_il1ra <- roc(t3c_dm ~ t3c_dm_adp_il1ra_prob,
                             data = d[sample_set == "Validation", ])

#auc(proc_t3c_dm_adp)
#ci.auc(proc_t3c_dm_adp)

#auc(proc_t3c_dm_il1ra)
#ci.auc(proc_t3c_dm_il1ra)

#auc(proc_t3c_dm_adp_il1ra)
#ci.auc(proc_t3c_dm_adp_il1ra)

 
#### Figure 4D ####
## ROC curves of T3Cdm vs NOD with adiponectin + IL-1Ra,
## in combination and individually
t3c_nod_melt <- melt_roc(d[sample_set == "Validation", ],
                                  d = "t3c_nod",
                                  m = c("adiponectin_serum_ug_m_l",
                                        "il_1ra_plasma_pg_m_l",
                                        "t3c_nod_adp_il1ra_prob"))

roc_t3c_nod <- ggplot(t3c_nod_melt,
                          aes(d = D, m = M, colour = name, linetype = name)) +
  geom_roc(labels = FALSE,
           n.cuts = 0) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1.)) +
  scale_x_continuous(expand = c(0,0), limits = c(0, 1)) +
  theme_cowplot(font_family = "Helvetica") +
  theme(legend.position = c(0.3, 0.3),
        legend.text = element_text(size = 12),
        legend.key.width=unit(1.5,"cm")
  ) +
  labs(linetype = NULL,
       colour = NULL,
       x = "1 - Specificity",
       y = "Sensitivity") +
  scale_colour_manual(labels = c("Adiponectin (0.75)", "IL-1Ra (0.78)", "Adiponectin + IL-1Ra (0.91)"),
                      values = c("#F8766D", "#00BA38", "Blue")) +
  scale_linetype_manual(labels = c("Adiponectin (0.75)", "IL-1Ra (0.78)", "Adiponectin + IL-1Ra (0.91)"),
                        values = c("solid", "longdash", "solid")
  )

## Export figure
save_plot(
  "output/figures/Figure 4D - IL-1Ra + Adiponectin ROC curves, typ3c vs NOD.pdf",
  roc_t3c_nod,
          base_asp = 1.1
)

## View ROC curves
#roc_t3c_nod

## Calculate AUCs plus 95% confidence intervals by DeLong method 
#calc_auc(roc_t3c_nod) 

### Using pROC package to calculate 95% CIs
proc_t3c_nod_adp <- roc(t3c_nod ~ adiponectin_serum_ug_m_l,
                       data = d[sample_set == "Validation", ]) 
proc_t3c_nod_il1ra <- roc(t3c_nod ~ il_1ra_plasma_pg_m_l,
                         data = d[sample_set == "Validation", ]) 
proc_t3c_nod_adp_il1ra <- roc(t3c_nod ~ t3c_nod_adp_il1ra_prob,
                              data = d[sample_set == "Validation", ])

#auc(proc_t3c_nod_adp)
#ci.auc(proc_t3c_nod_adp)

#auc(proc_t3c_nod_il1ra)
#ci.auc(proc_t3c_nod_il1ra)

#auc(proc_t3c_nod_adp_il1ra)
#ci.auc(proc_t3c_nod_adp_il1ra)


#### Supplementary figures ####

## Create longer form of data frame for easy faceting by biomarker
d_plot <- melt(
  d,
  id.vars =
    c("sample_id", "sample_set", "main_group", "subgroup",
      "simple_group", "plot_group_supp", "plot_group", "plot_group_scatter",
      "gender", "age", "bmi_kg_m2"),
  variable.name = "biomarker",
  value.name = "concentration"
)

## Create label variable to label y axes
d_plot[,
       biomarker_labs :=
         fcase(
           biomarker == "adiponectin_serum_ug_m_l", "Adiponectin (\U00B5g/mL)",
           biomarker == "adrenomedullin_plasma_ng_m_l", "Adrenomedullin (ng/mL)",
           biomarker == "aact_serum_ug_m_l", "AACT(\U00B5g/mL)",
           biomarker == "apo_a1_serum_mg_d_l", "Apo-A1 (mg/dL)",
           biomarker == "chemerin_serum_ng_m_l", "Chemerin (ng/mL)",
           biomarker == "clusterin_serum_ug_m_l", "Clusterin (\U00B5g/mL)",
           biomarker == "c_peptide_serum_ng_m_l", "C-peptide (ng/mL)",
           biomarker == "ghrelin_serum_ng_m_l", "Ghrelin\ (ng/mL)",
           biomarker == "gip_serum_ng_m_l", "GIP (ng/mL)",
           biomarker == "glp_1_serum_ng_m_l", "GLP-1 (ng/mL)",
           biomarker == "glucagon_serum_ng_m_l", "Glucagon (ng/mL)",
           biomarker == "ifn_g_plasma_pg_m_l", "IFN-G (pg/mL)",
           biomarker == "il_12_plasma_pg_m_l", "IL-12 (pg/mL)",
           biomarker == "il_1ra_plasma_ng_m_l", "IL-1Ra (ng/mL)",
           biomarker == "il_4_plasma_pg_m_l", "IL-4 (pg/mL)",
           biomarker == "il_6_plasma_pg_m_l", "IL-6 (pg/mL)",
           biomarker == "il_7_plasma_pg_m_l", "IL-7 (pg/mL)",
           biomarker == "il_8_plasma_pg_m_l", "IL-8 (pg/mL)",
           biomarker == "il_9_plasma_pg_m_l", "IL-9 (pg/mL)",
           biomarker == "insulin_serum_ng_m_l", "Insulin (ng/mL)",
           biomarker == "leptin_serum_ng_m_l", "Leptin (ng/mL)",
           biomarker == "mip_1a_plasma_pg_m_l", "MIP-1A (pg/mL)",
           biomarker == "mip_1b_plasma_pg_m_l", "MIP-1B (pg/mL)",
           biomarker == "pai_1_serum_ng_m_l", "PAI-1 (ng/mL)",
           biomarker == "pdgf_bb_plasma_ng_m_l", "PDGF-BB (ng/mL)",
           biomarker == "rantes_plasma_ng_m_l", "RANTES (ng/mL)",
           biomarker == "sparc_serum_ug_m_l", "SPARC (\U00B5g/mL)",
           biomarker == "transferin_serum_mg_m_l", "Transferrin (mg/mL)",
           biomarker == "tsp_1_serum_ug_m_l", "TSP-1 (\U00B5g/mL)",
           biomarker == "vwf_u_m_l", "VWF (U/mL)"
         )
]

#### Supplementary Figure 1 ####
## Plot all biomarkers from Set 2

### Function to create plot
plot_training <- function(biomarker, y_label){
  ggplot(
    d[sample_set == "Training", ],
    aes(plot_group_supp, !!sym(biomarker))) +
    stat_summary(geom = "crossbar", size = 0.2, fun = median) +
    stat_summary(geom = "errorbar",
                 fun.min = function(z){quantile(z, 0.25) },
                 fun.max = function(z){ quantile(z, 0.75)},
                 width = 0.5
    ) +
    geom_beeswarm(
      size = 0.5
    ) +
    labs(x = NULL,
         y = y_label) +
    theme_cowplot(
      font_family = "Helvetica"
    ) +
    theme(
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7)
    )
}

## Set 2 (Training)
### Plots are created individually to allow for individual y-axis scaling

## Biomarkers to iterate over
set2_biomarkers <- as.vector(unique(d_plot[sample_set == "Training" &
                                             !is.na(biomarker_labs),
                                           biomarker]))
## Axis labels to iterate over
set2_labels <- as.vector(unique(d_plot[sample_set == "Training" &
                                         !is.na(biomarker_labs),
                                       biomarker_labs]))

## Produce all plots for Supplementary Figure 1
plots2 <- map2(set2_biomarkers, set2_labels, plot_training)

## Extract individual plots for final figure assembly,
## rescaling y-axis if appropriate
set2_aact <- plots2[[3]]
set2_adiponectin <- plots2[[1]]
set2_adrenomedullin <- plots2[[2]]
set2_apoa1 <- plots2[[4]]
set2_chemerin <- plots2[[5]]
set2_ifng <- plots2[[6]]
set2_il4 <- plots2[[7]]
set2_il6 <- plots2[[8]]
set2_il7 <- plots2[[9]] + scale_y_continuous(trans = "log2")
set2_il8 <- plots2[[10]]
set2_il9 <- plots2[[11]] + scale_y_continuous(trans = "log2")
set2_il12 <- plots2[[12]]
set2_mip1a <- plots2[[13]]
set2_mip1b <- plots2[[14]]
set2_ghrelin <- plots2[[16]]
set2_clusterin <- plots2[[17]]
set2_cpeptide <- plots2[[18]]
set2_gip <- plots2[[19]] + scale_y_continuous(trans = "log2")
set2_glp1 <- plots2[[20]]
set2_glucagon <- plots2[[21]]
set2_il1ra <- plots2[[22]]
set2_insulin <- plots2[[23]]
set2_leptin <- plots2[[24]]
set2_pai1 <- plots2[[25]]
set2_pdgfbb <- plots2[[26]]
set2_sparc <- plots2[[27]]
set2_transferrin <- plots2[[28]] + scale_y_continuous(trans = "log2")
set2_vwf <- plots2[[29]]

## Create figure with patchwork package
### Letters to label plots
tags <- list(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m",
               "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z",
               "aa", "ab"))

### Assemble figure
set2_fig <- set2_adrenomedullin + set2_chemerin + set2_aact + set2_apoa1 +
  set2_clusterin + set2_sparc + set2_transferrin + set2_vwf +
  set2_ifng + set2_il1ra + set2_il4 + set2_il6 + set2_il7 +
  set2_il8 + set2_il9 + set2_il12 + 
  set2_mip1a + set2_mip1b + set2_pdgfbb + set2_adiponectin + set2_cpeptide +
  set2_ghrelin + set2_gip + set2_glp1 + set2_glucagon + 
  set2_insulin + set2_leptin + set2_pai1 +
  plot_layout(ncol = 4) +
  plot_annotation(
    tag_levels = tags
  )

## Save figure
ggsave2("output/figures/supplementary figure 1.pdf",
        set2_fig,
        width = 210,
        height = 297,
        units = "mm"
)

#### Supplementary Figure 2 ####
## Plot all biomarkers from Set 3 (validation)

## Function to plot validation data

plot_validation <- function(biomarker, y_label){
  ggplot(
    d[sample_set == "Validation", ],
    aes(plot_group_supp, !!sym(biomarker))) +
    stat_summary(geom = "crossbar", size = 0.2, fun = median) +
    stat_summary(geom = "errorbar",
                 fun.min = function(z){quantile(z, 0.25) },
                 fun.max = function(z){ quantile(z, 0.75)},
                 width = 0.5
    ) +
    geom_beeswarm(
      size = 0.5
    ) +
    labs(x = NULL,
         y = y_label) +
    theme_cowplot(
      font_family = "Helvetica"
    ) +
    theme(
      axis.text.x = element_text(size = 5),
      axis.text.y = element_text(size = 6),
      axis.title.y = element_text(size = 7)
    )
}

## Biomarkers to iterate over
set3_biomarkers <- as.vector(unique(d_plot[sample_set == "Validation" &
                                             !is.na(biomarker_labs),
                                           biomarker]))

## Axis labels to iterate over
set3_labels <- as.vector(unique(d_plot[sample_set == "Validation" &
                                         !is.na(biomarker_labs),
                                       biomarker_labs]))

## Produce all plots for Supplementary Figure 2
plots3 <- map2(set3_biomarkers, set3_labels, plot_validation)

## Extract individual plots for final figure assembly,
## rescaling y-axis if appropriate
set3_adiponectin <- plots3[[1]] + scale_y_continuous(trans = "log2")
set3_ifng <- plots3[[6]] + scale_y_continuous(trans = "log2")
set3_il4 <- plots3[[7]]
set3_il6 <- plots3[[8]]
set3_il8 <- plots3[[10]]
set3_il9 <- plots3[[11]]
set3_il12 <- plots3[[12]]
set3_mip1a <- plots3[[13]]
set3_mip1b <- plots3[[14]]
set3_tsp1 <- plots3[[15]]
set3_ghrelin <- plots3[[16]]
set3_cpeptide <- plots3[[18]]
set3_gip <- plots3[[19]]
set3_glp1 <- plots3[[20]]
set3_glucagon <- plots3[[21]]
set3_il1ra <- plots3[[22]]
set3_insulin <- plots3[[23]]
set3_leptin <- plots3[[24]]
set3_pai1 <- plots3[[25]]
set3_pdgfbb <- plots3[[26]]
set3_rantes <- plots3[[30]]

## Create figure with patchwork
set3_fig <- set3_tsp1 + set3_ifng + set3_il1ra + set3_il4 +
  set3_il6 + set3_il8 + set3_il9 + set3_il12 +
  set3_mip1a + set3_mip1b + set3_pdgfbb +   set3_rantes +
  set3_adiponectin + set3_cpeptide + set3_ghrelin + set3_gip +
  set3_glp1 + set3_glucagon + set3_insulin + set3_leptin +
  set3_pai1 +
  plot_layout(ncol = 3) +
  plot_annotation(
    tag_levels = "a"
  )

## Save plot
ggsave2("output/figures/supplementary figure 2.pdf",
        set3_fig,
        width = 210,
        height = 297,
        units = "mm"
)

#### P value multiple comparison adjustments ####

### Figure 2B
## Filter to Set 2 (training) and disease group of interest
d_filt <- d[sample_set == "Training" &
              main_group %in% c("PDAC", "PDAC-DM", "DM"),
]

## Wilcoxon test, holm adjustment
wilcox_pvals_2b <- pairwise.wilcox.test(d_filt$adiponectin_serum_ug_m_l,
                                                d_filt$main_group
)
#wilcox_pvals_2b


### Figure 2C
#### All group comparisons, no adjustment
## filter to Set 3 (validation) and disease groups of interest
d_filt <- d[sample_set == "Validation" &
              main_group %in% c("LSDM", "NOD", "PDAC", "PDAC-DM"),
]

## Wilcoxon test, holm adjustment
wilcox_pvals_2c <- pairwise.wilcox.test(d_filt$adiponectin_serum_ug_m_l, d_filt$main_group
)
#wilcox_pvals_2c

### Figure 3A
## filter to Set 2 (training) and disease groups of interest
d_filt <- d[sample_set == "Training" &
              main_group %in% c("PDAC-DM", "PDAC", "DM"),
]
## Wilcoxon test, holm adjustment
wilcox_pvals_3a <-  pairwise.wilcox.test(d_filt$il_1ra_plasma_pg_m_l,
                                         d_filt$main_group)
#wilcox_pvals_3a


### Figure 3B
#### Filter to set 3 (validation) and disease groups of interest
d_filt <- d[sample_set == "Validation" &
              main_group %in% c("LSDM", "NOD", "PDAC", "PDAC-DM"),
]


wilcox_pvals_3b <- pairwise.wilcox.test(d_filt$il_1ra_plasma_pg_m_l,
                                             d_filt$main_group)

#wilcox_pvals_3b

