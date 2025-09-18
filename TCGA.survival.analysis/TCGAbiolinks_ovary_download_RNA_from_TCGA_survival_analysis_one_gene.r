# Load libraries

library(TCGAbiolinks)
library(dplyr)
library(purrr)
library(DT)
library(SummarizedExperiment)
library(sesame)
library(sesameData)
library(janitor)
library(maftools)
library(survival)
library(survminer)

# Parallel processing
library(future.apply)
plan(multisession)

# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
# Arguments
# GDCdownload
# GDCquery
# GDCprepare

print("Download the clinical data : for all patients")

# -------------------------------------------
# Step 1: Query and download BCR Biotab clinical data
# -------------------------------------------
query <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

GDCdownload(query)
clinical_tab_all <- GDCprepare(query)

# -------------------------------------------
# Step 2: Explore available clinical tables
# -------------------------------------------
cat("Clinical tables available in TCGA-OV Biotab:\n")
print(names(clinical_tab_all))

# -------------------------------------------
# Step 3: Extract tables of interest
# -------------------------------------------
clinical_patient <- clinical_tab_all[["clinical_patient_ov"]]
clinical_follow_up <- clinical_tab_all[["clinical_follow_up_v1.0_ov"]]

clinical_drug <- clinical_tab_all[["clinical_drug_ov"]]

clinical_radiation <- clinical_tab_all[["clinical_radiation_ov"]]
clinical_omf <- clinical_tab_all[["clinical_omf_v4.0_ov"]]

clinical_nte <- clinical_tab_all[["clinical_nte_ov"]]
clinical_follow_up2 <- clinical_tab_all[["clinical_follow_up_v1.0_ov"]]

# Clinical tables available in TCGA-OV Biotab:
# [1] * "clinical_patient_ov"            * "clinical_nte_ov"               
# [3] * "clinical_drug_ov"               * "clinical_follow_up_v1.0_ov"    
# [5] * "clinical_follow_up_v1.0_nte_ov" * "clinical_radiation_ov"         
# [7] * "clinical_omf_v4.0_ov"    

# list.files()

print("Survival Analysis")

# Step 1: Query and download clinical data
query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Clinical",
  data.type = "Clinical Supplement",
  data.format = "BCR Biotab",
  access = "open"
)
GDCdownload(query)
clinical_list <- GDCprepare(query)

# Step 2: Extract patient-level clinical data
clin.OV <- clinical_list$clinical_patient_OV

# Basic structure
# str(clinical_list)
# class(clinical_list)
# length(clinical_list)

# First few elements
# head(clinical_list)

# If it's a data frame
summary(clinical_list)
names(clinical_list)

clinical_list 

# Step 2: Extract patient-level clinical data
clin.OV <- clinical_list$clinical_patient_ov

head(clin.OV)
print(colnames(clin.OV))
print("dim size clinOV:")
dim(clin.OV)

# If patient died: use days_to_death
# If patient is alive: use days_to_last_followup (censored observation)
# This creates a single time variable combining both scenarios.

# Clean and prepare data using actual column names
clin.OV <- clin.OV %>%
  filter(!is.na(gender), !is.na(vital_status)) %>%
  mutate(
    gender = as.character(gender),
    age = as.numeric(age_at_initial_pathologic_diagnosis),
    days_to_death = as.numeric(death_days_to),
    days_to_last_followup = as.numeric(last_contact_days_to),
    vital_status_bin = ifelse(vital_status == "Dead", 1, 0),
    time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death),  ## to verify the time
    TIME = coalesce(days_to_death, days_to_last_followup),                      ## to verify the time  
    Time = case_when(                                                           ## to verify the time
      !is.na(days_to_death) ~ days_to_death,
      !is.na(days_to_last_followup) ~ days_to_last_followup,        
      TRUE ~ NA_real_)
  ) %>%
  filter(!is.na(time), time > 0)

# If patient died: use days_to_death
# If patient is alive: use days_to_last_followup (censored observation)
# This creates a single time variable combining both scenarios.

# Optional: Create age groups (you can change median to another cutoff)
clin.OV <- clin.OV %>%
            mutate(age_group = ifelse(age > median(age, na.rm = TRUE), "Older", "Younger"))

# Step 4: Kaplan-Meier plots

# 4.1 Survival by Gender
fit_gender <- survfit(Surv(time, vital_status_bin) ~ gender, data = clin.OV)
ggsurvplot(fit_gender, data = clin.OV, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Survival by Gender")

# 4.2 Survival by Age Group
fit_age <- survfit(Surv(time, vital_status_bin) ~ age_group, data = clin.OV)
ggsurvplot(fit_age, data = clin.OV, pval = TRUE, risk.table = TRUE,
           title = "Kaplan-Meier Survival by Age Group")

# 4.3 Survival by Tumor Status (if available)

table(clin.OV$tumor_status)

# 4.3 Survival by Tumor Status (if available)
if ("tumor_status" %in% colnames(clin.OV)) {
  clin.OV <- clin.OV %>%
                     filter(tumor_status %in% c("TUMOR FREE", "WITH TUMOR"))
  
  fit_tumor <- survfit(Surv(time, vital_status_bin) ~ tumor_status, data = clin.OV)
  
  ggsurvplot(fit_tumor, data = clin.OV, pval = TRUE, risk.table = TRUE,
             title = "Kaplan-Meier Survival by Tumor Status (Filtered)")
}

# Count how many rows differ between columns
sum(clin.OV$time != clin.OV$TIME, na.rm = TRUE)
sum(clin.OV$time != clin.OV$Time, na.rm = TRUE)
sum(clin.OV$TIME != clin.OV$Time, na.rm = TRUE)

# Count how many rows differ between columns
sum(is.na(clin.OV$TIME))

head(clin.OV, 3)
tail(clin.OV, 3)



print("Survival analysis")

print("Query and download gene expression data")

# Step 1: Query RNA-seq gene expression (HTSeq - FPKM)
query_exp <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

# Step 2: Download the data
GDCdownload(query_exp)

# Step 3: Prepare the expression matrix (SummarizedExperiment object)
OV_se <- GDCprepare(query_exp)

# Save after first run
saveRDS(OV_se, file = "sep18.OV_expression_TCGA.rds")

# Later: load it directly
# OV_se <- readRDS("OV_se.rds")
# OV_rnaseq <- assay(OV_se)

# str(OV_se)

slotNames(OV_se)       # shows S4 slots (e.g., "assays", "colData", "rowRanges", etc.)
OV_se
cat("\n")
cat("the size of the data : genes * samples")
cat("\n")
dim(assay(OV_se))

head(colData(OV_se), 2)

head(rowData(OV_se), 2)

metadata(OV_se)

head(assay(OV_se), 3)
tail(assay(OV_se), 3)

# Step 4: Extract expression matrix
OV_rnaseq <- assay(OV_se)  # rows = genes, columns = samples

# Optional: log2 transform the matrix (add pseudocount to avoid log(0))
OV_rnaseq <- log2(OV_rnaseq + 1)
print(head(rownames(OV_rnaseq), 3))

# Remove version numbers from ENSG IDs (e.g., ENSG00000000003.15 -> ENSG00000000003)
rownames(OV_rnaseq) <- sub("\\.[0-9]+$", "", rownames(OV_rnaseq))

colnames(OV_rnaseq) <- substr(colnames(OV_rnaseq), 1, 12)
head(OV_rnaseq, 2)
tail(OV_rnaseq, 2)

dataOVcomplete <- OV_rnaseq 

# print(head(rownames(dataOVcomplete), 3))

print(colnames(clin.OV))
clinical_patient_Cancer = clin.OV

# Step 2: Filter and create required survival fields
clinical_patient_Cancer <- clinical_patient_Cancer  %>%
  mutate(
    death_days_to = as.numeric(death_days_to),
    last_contact_days_to = as.numeric(last_contact_days_to),
    vital_status = as.character(vital_status),
    bcr_patient_barcode = as.character(bcr_patient_barcode),
    time = ifelse(is.na(death_days_to), last_contact_days_to, death_days_to),
    vital_status_bin = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time), !is.na(vital_status_bin))

head(clinical_patient_Cancer, 2)
tail(clinical_patient_Cancer, 2)

cat("number of samples and features")
dim(clinical_patient_Cancer)

table(clinical_patient_Cancer$vital_status)
table(clinical_patient_Cancer$vital_status_bin)

# intersecting the barcode ID of the patientssubstr(samples, 1, 12)

print("codes : clinical metadata")
# unique(clinical_patient_Cancer$bcr_patient_barcode)
length(unique(clinical_patient_Cancer$bcr_patient_barcode))

print("codes : RNA seq")
# unique(colnames(OV_rnaseq))
length(unique(colnames(OV_rnaseq)))

# Get unique patient IDs from clinical metadata
clinical_codes <- unique(clinical_patient_Cancer$bcr_patient_barcode)

# Get unique column names from RNA-seq data
rna_codes <- unique(colnames(OV_rnaseq))

# Intersect
common_codes <- intersect(clinical_codes, rna_codes)

# Print results
cat("Number of clinical codes:", length(clinical_codes), "\n")
cat("Number of RNA-seq codes:", length(rna_codes), "\n")
cat("Number of common codes:", length(common_codes), "\n")

# If you want to see the first few common codes
head(common_codes)



clinical_patient_Cancer2 = clinical_patient_Cancer

grep("last", colnames(clinical_patient_Cancer2), value = TRUE)
grep("death", colnames(clinical_patient_Cancer2), value = TRUE)
grep("last|death", colnames(clinical_patient_Cancer2), value = TRUE, ignore.case = TRUE)

# Double check they exist
stopifnot("last_contact_days_to" %in% colnames(clinical_patient_Cancer2))
stopifnot("death_days_to" %in% colnames(clinical_patient_Cancer2))

# Clean up column names
colnames(clinical_patient_Cancer2) <- trimws(colnames(clinical_patient_Cancer2))  # Remove leading/trailing spaces

library(janitor)

clinical_patient_Cancer2 <- janitor::clean_names(clinical_patient_Cancer2)        # Make names consistent (snake_case)

# Rename using base R
colnames(clinical_patient_Cancer2)[colnames(clinical_patient_Cancer2) == "last_contact_days_to"] <- "days_to_last_follow_up"
colnames(clinical_patient_Cancer2)[colnames(clinical_patient_Cancer2) == "death_days_to"] <- "days_to_death"

# Check the result
print(colnames(clinical_patient_Cancer2))

table(clinical_patient_Cancer2$vital_status)
table(clinical_patient_Cancer2$vital_status_bin)

head(dataOVcomplete)
dim(dataOVcomplete)

# check overlap between the patient barcodes

colnames(dataOVcomplete) <- substr(colnames(dataOVcomplete ), 1, 12)

shared_barcodes <- intersect(
  colnames(dataOVcomplete ),
  clinical_patient_Cancer2$bcr_patient_barcode
)

print("shared barcodes")
length(shared_barcodes)  
head(clinical_patient_Cancer2, 2)

dataOVcomplete <- dataOVcomplete [, shared_barcodes]
dim(dataOVcomplete)

head(dataOVcomplete, 3)
dim(dataOVcomplete)

# You can check both duplicate rows and duplicate columns in your dataframe clinical_patient_Cancer2

# --- Duplicate rows ---
dup_rows <- sum(duplicated(clinical_patient_Cancer2))
cat("Number of duplicate rows:", dup_rows, "\n")

# If you want to see the duplicate rows themselves:
clinical_patient_Cancer2[duplicated(clinical_patient_Cancer2), ]

# --- Duplicate columns ---
dup_cols <- sum(duplicated(as.list(clinical_patient_Cancer2)))
cat("Number of duplicate columns:", dup_cols, "\n")

# If you want the names of duplicate columns:
names(clinical_patient_Cancer2)[duplicated(names(clinical_patient_Cancer2))]

# to make the names unique :
# names(clinical_patient_Cancer2) <- make.unique(names(clinical_patient_Cancer2))

dup_cols <- duplicated(as.list(clinical_patient_Cancer2))
which(dup_cols)

clinical_patient_Cancer2 <- clinical_patient_Cancer2[, !duplicated(as.list(clinical_patient_Cancer2))]
dim(clinical_patient_Cancer2)
head((clinical_patient_Cancer2$bcr_patient_barcode), 2)

length(shared_barcodes)
head(shared_barcodes, 2)

#shared_barcodes <- intersect(
#  colnames(dataOVcomplete ),
#  clinical_patient_Cancer2$bcr_patient_barcode
#)

clinical_patient_Cancer3 <- clinical_patient_Cancer2 %>% filter(bcr_patient_barcode %in% shared_barcodes)
dim(clinical_patient_Cancer3)

print(colnames(clinical_patient_Cancer3))

table(clinical_patient_Cancer3$vital_status_bin)
table(clinical_patient_Cancer3$vital_status)

# time = continuous variable (days, months, years)
# vital_status_bin = binary variable (1 = event, 0 = censored)
# vital_status = original categorical variable ("Dead"/"Alive")

# Remove unwanted columns

clinical_patient_Cancer3 <- clinical_patient_Cancer3 %>%
  select(-any_of(c(
    "form_completion_date", "prospective_collection", "retrospective_collection",
    "ethnicity", "jewish_religion_heritage_indicator", "history_other_malignancy",
    "method_initial_path_dx_other", "lymphovascular_invasion_indicator",
    "vascular_invasion_indicator", "karnofsky_score", "ecog_score",
    "performance_status_timing", "radiation_treatment_adjuvant",
    "pharmaceutical_tx_adjuvant", "treatment_outcome_first_course",
    "days_to_tumor_progression", "new_tumor_event_dx_indicator", "clinical_m",
    "clinical_n", "clinical_t", "days_to_patient_progression_free",
    "disease_code", "extranodal_involvement", "icd_10", "icd_o_3_histology",
    "icd_o_3_site", "informed_consent_verified", "pathologic_m", "pathologic_n",
    "pathologic_t", "pathologic_stage", "project_code", "stage_other",
    "system_version", "birth_days_to", "race", "history_neoadjuvant_treatment",
    "initial_pathologic_dx_year", "tumor_tissue_site"
  )))


# ! histological_type

# Check result
cat("Columns after removal:\n")
print(colnames(clinical_patient_Cancer3))

print(dim(clinical_patient_Cancer3))

clinical_patient_Cancer3 <- clinical_patient_Cancer3 %>%
  filter(
    !is.na(vital_status_bin),
    !is.na(time)
  )

head(clinical_patient_Cancer3, 2)
dim(clinical_patient_Cancer3)

if (all(c("vital_status_bin", "time") %in% colnames(clinical_patient_Cancer3))) {
  cat("Both 'vital_status_bin' and 'time' are present in the clinical data.\n")
} else {
  missing_cols <- c("vital_status_bin", "time")[!(c("vital_status_bin", "time") %in% colnames(clinical_patient_Cancer3))]
  cat("Missing column(s):", paste(missing_cols, collapse = ", "), "\n")
}

dim(clinical_patient_Cancer3) 

length(colnames(dataOVcomplete))
length(clinical_patient_Cancer2$bcr_patient_barcode)
length(shared_barcodes)

table(clinical_patient_Cancer2$vital_status)
length(clinical_patient_Cancer3$bcr_patient_barcode)

print(shared_barcodes)

# Check for NA or empty rows in expression
# Remove genes with too many zeros or no variance:

dataOVcomplete <- dataOVcomplete [, shared_barcodes]

keep_genes <- apply(dataOVcomplete , 1, function(x) {
  var_x <- var(x, na.rm = TRUE)
  nonzero <- sum(x > 0)
  !is.na(var_x) && var_x > 0 && nonzero >= 10
})

dim(dataOVcomplete)
dataOVcomplete <- dataOVcomplete [keep_genes, ]
dim(dataOVcomplete )

# Simply check if duplicates exist (TRUE/FALSE)
any(duplicated(colnames(dataOVcomplete)))  # TRUE if duplicates exist
any(duplicated(rownames(dataOVcomplete)))  # TRUE if duplicates exist

head(dataOVcomplete , 2)
tail(dataOVcomplete , 2)
# rownames(dataOVcomplete , 2)

# Remove version numbers from ENSG IDs
rownames(dataOVcomplete) <- sub("\\.[0-9]+$", "", rownames(dataOVcomplete))

head(dataOVcomplete , 2)
tail(dataOVcomplete , 2)

# Simply check if duplicates exist (TRUE/FALSE)
any(duplicated(colnames(dataOVcomplete)))  # TRUE if duplicates exist
any(duplicated(rownames(dataOVcomplete)))  # TRUE if duplicates exist


write.csv(dataOVcomplete, "sep18.OV_se.expression.complete.csv", row.names = TRUE, col.names= TRUE, quote = FALSE)
saveRDS(dataOVcomplete, "sep18.OV_se.expression.complete.rds")

gene_expr =   dataOVcomplete 
gene_expr_df <- as.data.frame(t(gene_expr))
gene_expr_df$barcode <- rownames(gene_expr_df)

dim(gene_expr_df)

# Prepare clinical data with matching barcodes

clinical_data = clinical_patient_Cancer3
c("vital_status_bin", "time") %in% colnames(clinical_patient_Cancer3)

clinical_data <- clinical_data %>%
                  mutate(barcode = substr(bcr_patient_barcode, 1, 12))

print(colnames(clinical_data))
head(clinical_data,2)
tail(clinical_data,2)
dim(clinical_data)

# Simply check if duplicates exist (TRUE/FALSE)
any(duplicated(colnames(clinical_data)))  # TRUE if duplicates exist
any(duplicated(rownames(clinical_data)))  # TRUE if duplicates exist

dim(clinical_data)
colnames(clinical_data)
print(colnames(clinical_data))

head(gene_expr_df)

# Merge clinical and expression data
merged_df <- inner_join(clinical_data, gene_expr_df, by = "barcode")

cat("integrated matrix : gene expression, tumor samples, and metadata")
dim(merged_df)

head(merged_df, 2)
tail(merged_df, 2)

# verifications : 
table(merged_df$vital_status)
table(merged_df$vital_status_bin)

# if NA values in the time column :
sum(is.na(merged_df$time))

write.csv(merged_df, "sep18.OV_se.expression.and.clinical.info.csv", row.names = FALSE, col.names= TRUE, quote = FALSE)
print(colnames(merged_df[1:20]))



set.seed(1234)  

df = merged_df

colnames(df[1:50])
# rownames(df[1:10])

table(df$vital_status)
table(df$vital_status_bin)

print("survival analysis : test one gene")
set.seed(1234)

# Make a copy of the merged dataset to avoid mutating the original
df <- merged_df    # ← This is important!

table(df$vital_status)
table(df$vital_status_bin)

# HR > 1: High expression increases death risk (worse survival)
# HR < 1: High expression decreases death risk (better survival)
# HR = 1: No difference between groups

# Interpretation of the forest plot

# Reference group (Low expression)
# Low (N=123) is set as the reference, hazard ratio = 1.0 by definition.
# High expression group
# High (N=123) has a hazard ratio (HR) of 0.82
# 95% CI = (0.58 – 1.20).

# Interpretation: Patients with High expression have about an 18% lower hazard of death compared to Low, 
# but the confidence interval crosses 1.0 → the effect is not statistically significant.

# p-value
# Shown as 0.267 (Cox regression, Wald test).
# This is >0.05, so the null hypothesis (no survival difference between groups) cannot be rejected.

# Global log-rank test
# At the bottom: Global p-value (Log-Rank): 0.26697.
# This is the same conclusion as above: no significant survival difference between groups.

# Other metrics

# Events: 133 → number of deaths/events observed.
# AIC: 1168.36 → model fit criterion (lower AIC is better, useful for comparing models).
# Concordance index: 0.52 → predictive accuracy (0.5 = random chance, 1.0 = perfect; 
# so here the model is weakly predictive).

# High expression appears protective (HR < 1), but not significantly.
# Both the Cox model and log-rank test agree: no significant survival difference between High vs Low groups.

# * Hazard Ratio (HR) < 1
# The Cox model gave HR = 0.82 for High expression vs Low.
# HR < 1 means the High group had a lower risk of death compared to the reference (Low).
# Specifically, 0.82 means about an 18% reduction in hazard.
# So if this were statistically significant, we’d conclude that high expression is protective (linked to better survival).
# * Confidence Interval (CI) includes 1
# The 95% CI was (0.58 – 1.20).
# Since the interval crosses 1.0, we can’t rule out the possibility that there’s no difference at all.
# * p-values (Cox + log-rank)
# Cox model p-value (Wald test): 0.267
# Global log-rank test p-value: 0.267
# Both > 0.05 → not statistically significant.
# That means the apparent protective trend could be due to chance, not a real biological effect.
# * Putting it together
# The direction of the effect (HR < 1) suggests that high expression might be beneficial.
# But the evidence is weak (not statistically significant).
# Both methods (CoxPH and log-rank test) tell the same story → no clear survival difference between groups.

# To test the proportional hazards (PH) assumption in a Cox model using Schoenfeld residuals. 
# Output interpretation

# The cox.zph() test gives you:
# A p-value for each covariate → if > 0.05 → assumption holds (no strong evidence of non-proportional hazards).
# A global p-value → checks the whole model.
# 👉 If p < 0.05 → evidence that the PH assumption may be violated.

# dfbeta values (influence diagnostics) for your Cox model, not Schoenfeld residuals.

# What this plot shows:

# Y-axis: Change in regression coefficient (dfbeta values)
# X-axis: Individual observations (patients) in your dataset
# Red dashed line: Reference line at 0

# Interpretation: Good news: Your model looks stable!

# Most values are close to 0 - This means most patients don't dramatically influence the model
# No extreme outliers - Values range roughly from -0.04 to +0.02, which is reasonable
# Random scatter pattern - No systematic trends or clusters of influential points

# What to look for (red flags):

# Large spikes (|dfbeta| > 0.1) - Would indicate highly influential patients
# Systematic patterns - Would suggest model misspecification
# Clusters of high values - Could indicate subgroups with different effects

# Rule of thumb:

# |dfbeta| < 2/√n is usually acceptable
# For our data: 2/√n ≈ 2/√250 ≈ 0.13
# Your values are well below this threshold

# What to Do If PH Assumption Is Violated
# Stratify the model: coxph(... + strata(variable))
# Add time-dependent covariates: e.g., coxph(Surv(...) ~ var + tt(var), tt = function(x, t, ...) x * log(t))



cat("gene of interest : PPP4R2 : ENSG00000163605")

# ---- Display options (Jupyter/RMarkdown) ----
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 150)

# ---- Packages ----
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)

# ---- Gene of interest ----
gene_to_test <- "ENSG00000163605"
gene_to_test <- gsub("[^[:alnum:]_.-]", "_", gene_to_test)  # safe for filenames

# ---- Sanity checks ----
stopifnot(gene_to_test %in% colnames(df))
stopifnot(all(c("time","vital_status_bin") %in% colnames(df)))

# ---- Filter to rows actually used in survival analysis ----
df0 <- df %>%
  filter(!is.na(time), !is.na(vital_status_bin), !is.na(.data[[gene_to_test]]))

cat("N with valid survival & gene values:", nrow(df0), "\n")
print(summary(df0[[gene_to_test]]))

# ---- Define cutoffs (tertiles) on the ANALYSIS subset ----
if (!exists("cutoff_high") || !exists("cutoff_low") ||
    is.na(cutoff_high) || is.na(cutoff_low)) {
  qs <- quantile(df0[[gene_to_test]], probs = c(0.33, 0.67), na.rm = TRUE, type = 1)
  cutoff_low  <- qs[1]
  cutoff_high <- qs[2]
}
cat("Cutoff low:", cutoff_low, " | Cutoff high:", cutoff_high, "\n")

# ---- Create expression groups ----
df0 <- df0 %>%
  mutate(
    expr_group = case_when(
      .data[[gene_to_test]] >= cutoff_high ~ "High",
      .data[[gene_to_test]] <= cutoff_low  ~ "Low",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(expr_group)) %>%
  mutate(expr_group = factor(expr_group, levels = c("Low","High")))

cat("N after grouping:", nrow(df0), "\n"); print(table(df0$expr_group))

# Guard against degenerate splits
stopifnot(length(unique(df0$expr_group)) == 2 && all(table(df0$expr_group) > 0))

# ---- Survival objects ----
surv_obj <- Surv(time = df0$time, event = df0$vital_status_bin)

# ---- KM + explicit log-rank ----
fit <- survfit(surv_obj ~ expr_group, data = df0)
logrank_test <- survdiff(surv_obj ~ expr_group, data = df0)
p_logrank <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
cat("\nLog-rank p-value:", signif(p_logrank, 4), "\n")

km_plot <- ggsurvplot(
  fit, data = df0, pval = TRUE, conf.int = TRUE,
  risk.table = TRUE, risk.table.height = 0.22,
  censor.shape = 124, censor.size = 3,
  xlab = "Time (days)", ylab = "Overall survival probability",
  title = paste("Survival by expression of", gene_to_test),
  legend.title = "Expression", legend.labs = c("Low","High"),
  ggtheme = theme_minimal(base_size = 12)
)
print(km_plot$plot); print(km_plot$table)

ggsave(paste0(gene_to_test, "_KM_plot_expression.jpg"),
       plot = km_plot$plot, width = 8, height = 6, dpi = 300)
ggsave(paste0(gene_to_test, "_KM_plot_risktable.jpg"),
       plot = km_plot$table, width = 8, height = 3, dpi = 300)

# ---- Cox PH + forest ----
cox_fit <- coxph(surv_obj ~ expr_group, data = df0)
cox_summary <- summary(cox_fit)
cat("\nCox PH summary (High vs Low):\n"); print(cox_summary)

hr       <- cox_summary$coefficients[1, "exp(coef)"]
p_cox    <- cox_summary$coefficients[1, "Pr(>|z|)"]
ci_lo    <- cox_summary$conf.int[1, "lower .95"]
ci_hi    <- cox_summary$conf.int[1, "upper .95"]

cox_results <- data.frame(
  Variable = rownames(cox_summary$coefficients)[1],
  HR = round(hr, 3),
  CI = sprintf("(%0.3f – %0.3f)", ci_lo, ci_hi),
  p.value = signif(p_cox, 3)
)
print(cox_results)

forest <- ggforest(cox_fit, data = df0)
print(forest)
ggsave(paste0(gene_to_test, "_cox_forest.png"),
       plot = forest, width = 7, height = 5, dpi = 300)

# ---- Median survival by group ----
cat("\nMedian survival by group:\n"); print(surv_median(fit))

# ---- PH diagnostics (Schoenfeld) ----
cat("\nTest proportional hazards assumptions\n")
ph_test <- cox.zph(cox_fit)
print(ph_test)

# Save panel with a proper outer title
png(paste0(gene_to_test, "_schoenfeld_residuals.png"),
    width = 8, height = 6, units = "in", res = 300)
op <- par(no.readonly = TRUE); par(oma = c(0, 0, 2.2, 0))
plot(ph_test)
mtext("Schoenfeld Residuals Test", outer = TRUE, line = 0.8, cex = 1.1)
par(op); dev.off()

cat("\n--- Scaled Schoenfeld Residuals (transform = 'identity') ---\n")
scaled_ph_test <- cox.zph(cox_fit, transform = "identity")
print(scaled_ph_test)

png(paste0(gene_to_test, "_schoenfeld_residuals_scaled.png"),
    width = 8, height = 6, units = "in", res = 300)
op <- par(no.readonly = TRUE); par(oma = c(0, 0, 2.2, 0))
plot(scaled_ph_test)
mtext("Scaled Schoenfeld Residuals (identity transform)", outer = TRUE, line = 0.8, cex = 1.1)
par(op); dev.off()

# ---- Influential observations (dfbeta) ----
dfbeta_vals <- residuals(cox_fit, type = "dfbeta")
matplot(dfbeta_vals, type = "l", lty = 1,
        main = "Influential Observations (dfbeta)",
        xlab = "Observation index", ylab = "dfbeta values")
abline(h = 0, lty = 2)

png(paste0(gene_to_test, "_dfbeta_influential.png"),
    width = 8, height = 6, units = "in", res = 300)
matplot(dfbeta_vals, type = "l", lty = 1,
        main = paste("Influential Observations (dfbeta) -", gene_to_test),
        xlab = "Observation index", ylab = "dfbeta values")
abline(h = 0, lty = 2)
dev.off()

# ---- Model fit (concordance) ----
c_index <- cox_summary$concordance[1]
cat("\nC-index (concordance):", round(c_index, 3), "\n")


dim(merged_df)
dim(df)

# ?Surv
# ?survfit
# ?survdiff
# ?coxph
# ?cox.zph

# cat("the cox object\n")
# cox_fit
# cat("\nsummary of cox\n")
# summary(cox_fit)
# cat("\nthe structure of cox object\n")
# str(cox_fit)

# survdiff performs a log-rank test (which compares survival curves), it doesn't explicitly fit a Cox proportional hazards model. 
# If you wanted to assess hazard ratios and adjust for covariates, you would use the coxph() function from the survival package.

# scaling on multiple genes 

print("Survival Analysis : all the genes in A2780_targets.txt")

library(readr)

# Read, trim, drop empties, deduplicate, strip version suffixes like ".12"

genes <- read_lines("A2780_targets.txt", progress = FALSE)
genes <- unique(trimws(genes))
genes <- genes[nzchar(genes)]
genes <- sub("\\.[0-9]+$", "", genes)  # ENSG00000149948.12 -> ENSG00000149948

# Check which are present in your data frame 'df' (genes are columns)
genes_present <- intersect(genes, colnames(df))
genes_missing <- setdiff(genes, colnames(df))
cat("Present:", length(genes_present), " | Missing:", length(genes_missing), "\n")
if (length(genes_missing)) writeLines(genes_missing, "missing_genes.txt")


gene_list = genes_present
# print(gene_list)
print(paste("Number of genes:", length(gene_list)))

# To be continued in a script that loops over each gene ! 

# Prepare expression data
# plan(multisession)  # use parallel processing

# Other functions :
# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/analysis.html

# TCGAanalyze: Analyze data from TCGA.
# TCGAanalyze_Preprocessing: Preprocessing of Gene Expression data (IlluminaHiSeq_RNASeqV2)

# TCGAanalyze_DEA & TCGAanalyze_LevelTab: Differential expression analysis (DEA)
# TCGAanalyze_EAcomplete & TCGAvisualize_EAbarplot: Enrichment Analysis

# TCGAanalyze_survival: Survival Analysis
# TCGAanalyze_SurvivalKM: Correlating gene expression and Survival Analysis
# TCGAanalyze_DMR: Differentially methylated regions Analysis

# TCGAvisualize_Heatmap: Create heatmaps with cluster bars
# TCGAvisualize_Volcano: Create volcano plot
# TCGAvisualize_PCA: Principal Component Analysis plot for differentially expressed genes
# TCGAvisualize_meanMethylation: Mean DNA Methylation Analysis

# TCGAvisualize_starburst: Integration of gene expression and DNA methylation data


