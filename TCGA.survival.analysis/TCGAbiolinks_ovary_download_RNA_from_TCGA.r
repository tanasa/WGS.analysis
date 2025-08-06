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

list.files()

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

clinical_list 

# Step 2: Extract patient-level clinical data
clin.OV <- clinical_list$clinical_patient_ov

head(clin.OV)
colnames(clin.OV)

# Clean and prepare data using actual column names
clin.OV <- clin.OV %>%
  filter(!is.na(gender), !is.na(vital_status)) %>%
  mutate(
    gender = as.character(gender),
    age = as.numeric(age_at_initial_pathologic_diagnosis),
    days_to_death = as.numeric(death_days_to),
    days_to_last_followup = as.numeric(last_contact_days_to),
    vital_status_bin = ifelse(vital_status == "Dead", 1, 0),
    time = ifelse(is.na(days_to_death), days_to_last_followup, days_to_death)
  ) %>%
  filter(!is.na(time), time > 0)

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

head(clin.OV, 3)



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

# Step 4: Extract expression matrix
OV_rnaseq <- assay(OV_se)  # rows = genes, columns = samples

# Optional: log2 transform the matrix (add pseudocount to avoid log(0))
OV_rnaseq <- log2(OV_rnaseq + 1)

head(OV_rnaseq, 2)

# Save after first run
saveRDS(OV_se, file = "OV_se.rds")

# Later: load it directly
# OV_se <- readRDS("OV_se.rds")
# OV_rnaseq <- assay(OV_se)

str(OV_se)

slotNames(OV_se)       # shows S4 slots (e.g., "assays", "colData", "rowRanges", etc.)

OV_se

colData(OV_se)

rowData(OV_se)



colnames(OV_rnaseq) <- substr(colnames(OV_rnaseq), 1, 12)
dataOVcomplete <- log2(OV_rnaseq + 1)
colnames(clin.OV)
clinical_patient_Cancer = clin.OV

# Step 2: Optional — filter and create required survival fields
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
colnames(clinical_patient_Cancer)


