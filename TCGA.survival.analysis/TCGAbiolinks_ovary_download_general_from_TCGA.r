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

# Parallel processing

library(future.apply)
plan(multisession)

# Harmonized database examples
# DNA methylation data: Recurrent tumor samples

# query <- GDCquery(
#    project = c("TCGA-OV", "TCGA-LGG"),
#    data.category = "DNA Methylation",
#    platform = c("Illumina Human Methylation 450"),
#    sample.type = "Recurrent Tumor"
# )

# datatable(
#    getResults(query), 
#    filter = 'top',
#    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#    rownames = FALSE
# )

# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html
# Arguments
# GDCdownload
# GDCquery
# GDCprepare

# If a SummarizedExperiment object was chosen, the data can be accessed with three different accessors: 
# assay for the data information, rowRanges to gets the range of values in each row and colData to get 
# the sample information (patient, batch, sample type, etc) (Huber et al. 2015; H. J. Morgan M Obenchain V and H., n.d.). 
# An example is shown in listing below.



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


# ============================================================================
# STEP 2: Explore Available Tables
# ============================================================================

print("Available clinical tables in TCGA-OV:")
print(names(clinical_tab_all))

# ============================================================================
# STEP 3: Extract and View Each Table
# ============================================================================

# Extract each table
clinical_patient <- clinical_tab_all[["clinical_patient_ov"]]
clinical_nte <- clinical_tab_all[["clinical_nte_ov"]]
clinical_drug <- clinical_tab_all[["clinical_drug_ov"]]
clinical_follow_up <- clinical_tab_all[["clinical_follow_up_v1.0_ov"]]
clinical_follow_up_nte <- clinical_tab_all[["clinical_follow_up_v1.0_nte_ov"]]
clinical_radiation <- clinical_tab_all[["clinical_radiation_ov"]]
clinical_omf <- clinical_tab_all[["clinical_omf_v4.0_ov"]]

# ============================================================================
# STEP 4: Display Table Information
# ============================================================================

# Function to display table summary
display_table_info <- function(table_name, data) {
  cat("\n", "=", 60, "\n")
  cat("TABLE:", table_name, "\n")
  cat("=", 60, "\n")
  cat("Dimensions:", dim(data), "\n")
  cat("Number of rows:", nrow(data), "\n")
  cat("Number of columns:", ncol(data), "\n")
  cat("\nColumn names:\n")
  print(colnames(data))
  cat("\nFirst few rows:\n")
  print(head(data, 3))
  cat("\nData types:\n")
  print(sapply(data, class))
  cat("\n", "=", 60, "\n")
}

# Display information for each table
display_table_info("clinical_patient_ov", clinical_patient)

# Display information for each table
display_table_info("clinical_nte_ov", clinical_nte)

# Display information for each table
display_table_info("clinical_drug_ov", clinical_drug)

# Display information for each table
display_table_info("clinical_follow_up_v1.0_ov", clinical_follow_up)

# Display information for each table
display_table_info("clinical_follow_up_v1.0_nte_ov", clinical_follow_up_nte)

# Display information for each table
display_table_info("clinical_radiation_ov", clinical_radiation)

# Display information for each table
display_table_info("clinical_omf_v4.0_ov", clinical_omf)



# Check if the column exists
if ('pharmaceutical_therapy_drug_name' %in% colnames(clinical_drug)) {
  cat('Unique drug names in pharmaceutical_therapy_drug_name column:\n')
  unique_drugs <- unique(clinical_drug$pharmaceutical_therapy_drug_name)
  unique_drugs <- unique_drugs[!is.na(unique_drugs) & unique_drugs != '']
  print(sort(unique_drugs))
  cat('\nTotal unique drugs:', length(unique_drugs), '\n')
} else {
  cat('Column pharmaceutical_therapy_drug_name not found.\n')
  cat('Available columns:\n')
  print(colnames(clinical_drug))
}

# Extract Cisplatin Data from TCGA-OV Clinical Drug Data
# Fixed version to handle data structure issues

library(TCGAbiolinks)
library(dplyr)

# ============================================================================
# STEP 1: Download and Prepare Clinical Drug Data
# ============================================================================

print("Downloading TCGA-OV clinical data...")

# Query TCGA-OV clinical data
query <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Clinical",
  data.type = "Clinical Supplement", 
  data.format = "BCR Biotab"
)

# Download the data
GDCdownload(query)

# Prepare the data
clinical_tab_all <- GDCprepare(query)

# Extract drug data and ensure it's a data frame
clinical_drug <- clinical_tab_all[["clinical_drug_ov"]]

# Check data structure
cat("Data structure check:\n")
cat("Class of clinical_drug:", class(clinical_drug), "\n")
cat("Dimensions:", dim(clinical_drug), "\n")

# Convert to data frame if it's a list
if (is.list(clinical_drug) && !is.data.frame(clinical_drug)) {
  print("Converting list to data frame...")
  clinical_drug <- as.data.frame(clinical_drug, stringsAsFactors = FALSE)
}

# Check columns
cat("Available columns:\n")
print(colnames(clinical_drug))

# ============================================================================
# STEP 2: Extract Cisplatin Data from pharmaceutical_therapy_drug_name
# ============================================================================

print("Extracting cisplatin data from pharmaceutical_therapy_drug_name column...")

# Check if the column exists
if (!"pharmaceutical_therapy_drug_name" %in% colnames(clinical_drug)) {
  cat("Column 'pharmaceutical_therapy_drug_name' not found!\n")
  cat("Available columns:\n")
  print(colnames(clinical_drug))
  stop("Please check the correct column name")
}

# Extract rows where pharmaceutical_therapy_drug_name contains "cisplatin"
cisplatin_data <- clinical_drug %>%
  filter(grepl("cisplatin", pharmaceutical_therapy_drug_name, ignore.case = TRUE))

# Ensure cisplatin_data is a data frame
if (is.list(cisplatin_data) && !is.data.frame(cisplatin_data)) {
  cisplatin_data <- as.data.frame(cisplatin_data, stringsAsFactors = FALSE)
}

# ============================================================================
# STEP 3: Display Results
# ============================================================================

cat("\n", "=", 60, "\n")
cat("CISPLATIN TREATMENT DATA\n")
cat("=", 60, "\n")
cat("Total rows found:", nrow(cisplatin_data), "\n")
cat("Column used: pharmaceutical_therapy_drug_name\n")
cat("Data class:", class(cisplatin_data), "\n")
cat("=", 60, "\n")

# Display all columns
print("All columns in cisplatin data:")
print(colnames(cisplatin_data))

cat("\n", "=", 60, "\n")
cat("CISPLATIN DATA PREVIEW\n")
cat("=", 60, "\n")

# Show first few rows
print(head(cisplatin_data, 5))

# Show all unique drug names found (to see variations)
cat("\n", "=", 60, "\n")
cat("UNIQUE DRUG NAMES CONTAINING 'CISPLATIN'\n")
cat("=", 60, "\n")
unique_cisplatin_names <- unique(cisplatin_data$pharmaceutical_therapy_drug_name)
print(sort(unique_cisplatin_names))

# ============================================================================
# STEP 4: Summary Statistics
# ============================================================================

cat("\n", "=", 60, "\n")
cat("SUMMARY OF KEY VARIABLES\n")
cat("=", 60, "\n")

# Patient count
if ("bcr_patient_barcode" %in% colnames(cisplatin_data)) {
  unique_patients <- unique(cisplatin_data$bcr_patient_barcode)
  cat("Unique patients treated with cisplatin:", length(unique_patients), "\n")
  
  # Show unique patients
  cat("\nUnique patient barcodes:\n")
  print(unique_patients)
} else {
  cat("Column 'bcr_patient_barcode' not found in cisplatin data\n")
}

# Show summary of other key variables
if ("pharmaceutical_therapy_type" %in% colnames(cisplatin_data)) {
  cat("\nTherapy types:\n")
  print(table(cisplatin_data$pharmaceutical_therapy_type))
}

if ("pharmaceutical_therapy_dosage" %in% colnames(cisplatin_data)) {
  cat("\nDosage summary:\n")
  print(summary(cisplatin_data$pharmaceutical_therapy_dosage))
}

if ("pharmaceutical_therapy_duration" %in% colnames(cisplatin_data)) {
  cat("\nDuration summary:\n")
  print(summary(cisplatin_data$pharmaceutical_therapy_duration))
}

# ============================================================================
# STEP 5: Treatments per Patient Analysis
# ============================================================================

cat("\n", "=", 60, "\n")
cat("TREATMENTS PER PATIENT\n")
cat("=", 60, "\n")

# Count cisplatin treatments per patient (with error handling)
if ("bcr_patient_barcode" %in% colnames(cisplatin_data)) {
  tryCatch({
    treatments_per_patient <- cisplatin_data %>%
      count(bcr_patient_barcode, sort = TRUE)
    
    cat("Cisplatin treatments per patient:\n")
    print(treatments_per_patient)
    
    # Show patients with multiple cisplatin treatments
    multiple_treatments <- treatments_per_patient %>%
      filter(n > 1)
    
    if (nrow(multiple_treatments) > 0) {
      cat("\nPatients with multiple cisplatin treatments:\n")
      print(multiple_treatments)
      
      cat("\nDetailed data for patients with multiple treatments:\n")
      for (patient in multiple_treatments$bcr_patient_barcode) {
        cat("\nPatient:", patient, "\n")
        patient_data <- cisplatin_data %>%
          filter(bcr_patient_barcode == patient)
        print(patient_data)
      }
    }
  }, error = function(e) {
    cat("Error in treatments per patient analysis:", e$message, "\n")
    cat("Trying alternative approach...\n")
    
    # Alternative approach using base R
    patient_counts <- table(cisplatin_data$bcr_patient_barcode)
    cat("Cisplatin treatments per patient (base R):\n")
    print(patient_counts)
  })
} else {
  cat("Cannot analyze treatments per patient - bcr_patient_barcode column not found\n")
}

# ============================================================================
# STEP 6: Save Results
# ============================================================================

# Save to CSV
tryCatch({
  write.csv(cisplatin_data, "tcga_ov_cisplatin_treatment_data.csv", row.names = FALSE)
  
  cat("\n", "=", 60, "\n")
  cat("DATA SAVED\n")
  cat("=", 60, "\n")
  cat("Cisplatin data saved to: tcga_ov_cisplatin_treatment_data.csv\n")
  cat("=", 60, "\n")
}, error = function(e) {
  cat("Error saving data:", e$message, "\n")
})

# ============================================================================
# STEP 7: Additional Analysis
# ============================================================================

cat("\n", "=", 60, "\n")
cat("ADDITIONAL ANALYSIS\n")
cat("=", 60, "\n")

# Show combination therapies (if cisplatin is combined with other drugs)
cat("Drug combinations with cisplatin:\n")
for (i in 1:min(nrow(cisplatin_data), 10)) {  # Limit to first 10 for display
  drug_name <- cisplatin_data$pharmaceutical_therapy_drug_name[i]
  if ("bcr_patient_barcode" %in% colnames(cisplatin_data)) {
    patient <- cisplatin_data$bcr_patient_barcode[i]
    cat("Patient:", patient, "- Drug:", drug_name, "\n")
  } else {
    cat("Row", i, "- Drug:", drug_name, "\n")
  }
}

# Show therapy timing if available
if ("pharmaceutical_therapy_start_date" %in% colnames(cisplatin_data)) {
  cat("\nTherapy start dates:\n")
  print(table(cisplatin_data$pharmaceutical_therapy_start_date))
}

if ("pharmaceutical_therapy_end_date" %in% colnames(cisplatin_data)) {
  cat("\nTherapy end dates:\n")
  print(table(cisplatin_data$pharmaceutical_therapy_end_date))
}

# ============================================================================
# STEP 8: Data Structure Summary
# ============================================================================

cat("\n", "=", 60, "\n")
cat("DATA STRUCTURE SUMMARY\n")
cat("=", 60, "\n")
cat("Original clinical_drug class:", class(clinical_tab_all[["clinical_drug_ov"]]), "\n")
cat("Processed clinical_drug class:", class(clinical_drug), "\n")
cat("Cisplatin data class:", class(cisplatin_data), "\n")
cat("Cisplatin data dimensions:", dim(cisplatin_data), "\n")

print("Cisplatin data extraction completed!") 



query <- GDCquery(
    project = "TCGA-OV", 
    data.category = "Sequencing Reads",
    data.type = "Aligned Reads", 
    data.format = "bam",
    workflow.type = "STAR 2-Pass Chimeric"
)
# Only first 10 to make render faster
datatable(
    getResults(query, rows = 1:10,cols = c("file_name","cases")), 
    filter = 'top',
    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
    rownames = FALSE
)

# getManifest(query, save = FALSE) 

print("Download ATAC files")

# Query ATAC-seq metadata
atac_query <- TCGAbiolinks:::GDCquery_ATAC_seq()

# Extract results and select desired columns
atac_results <- getResults(atac_query)[, c("file_name", "file_size")]

# Display as an interactive datatable
datatable(
  atac_results,
  filter = 'top',
  options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
  rownames = FALSE
)

# https://rpubs.com/tiagochst/atac_seq_workshop

# query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "rds") 
# GDCdownload(query, method = "client")

# query <- TCGAbiolinks:::GDCquery_ATAC_seq(file.type = "bigWigs") 
# GDCdownload(query, method = "client")

# Samples with DNA methylation and gene expression data
# We will access the harmonized database and search for all patients with DNA methylation (platform HumanMethylation450k) 
# and gene expression data for Colon Adenocarcinoma tumor (TCGA-COAD).



print("Download Methylation and Gene Expression files")

query_met <- GDCquery(
  project = "TCGA-OV",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",        
  platform = "Illumina Human Methylation 450"
)

query_exp <- GDCquery(
    project = "TCGA-OV",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts"
)

# Get all patients that have DNA methylation and gene expression.
common.patients <- intersect(
    substr(getResults(query_met, cols = "cases"), 1, 12),
    substr(getResults(query_exp, cols = "cases"), 1, 12)
)

# print(common.patients)

# Only seelct the first 5 patients
# query_met <- GDCquery(
#    project = "TCGA-COAD",
#    data.category = "DNA Methylation",
#    platform = c("Illumina Human Methylation 450"),
#    barcode = common.patients[1:5]
# )

# query_exp <- GDCquery(
#    project = "TCGA-COAD",
#    data.category = "Transcriptome Profiling",
#    data.type = "Gene Expression Quantification", 
#    workflow.type = "STAR - Counts",
#    barcode = common.patients[1:5]
# )

# Extract barcodes and identify common patients
all_met_barcodes <- getResults(query_met, cols = "cases")
all_exp_barcodes <- getResults(query_exp, cols = "cases")

all_met_patients <- substr(all_met_barcodes, 1, 12)
all_exp_patients <- substr(all_exp_barcodes, 1, 12)

common_patients <- intersect(all_met_patients, all_exp_patients)
print(common_patients)

# ------------------ SELECT FIRST PATIENT ------------------

first_patient <- common_patients[1]

# Get methylation barcodes for first patient
met_patient_barcodes <- all_met_barcodes[all_met_patients == first_patient]

# Get expression barcodes for first patient
exp_patient_barcodes <- all_exp_barcodes[all_exp_patients == first_patient]

print(exp_patient_barcodes)

# datatable(
#    getResults(query_met, cols = c("data_type","cases")),
#    filter = 'top',
#    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#    rownames = FALSE
# )

# datatable(
#    getResults(query_exp, cols = c("data_type","cases")), 
#    filter = 'top',
#    options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#    rownames = FALSE
# )



print("METHYLATION VALUES")

# ------------------ METHYLATION DATA ------------------

query_met_patient <- GDCquery(
  project = "TCGA-OV",
  data.category = "DNA Methylation",
  data.type = "Methylation Beta Value",
  platform = "Illumina Human Methylation 450",
  barcode = met_patient_barcodes
)

GDCdownload(query_met_patient)
met_data_patient <- GDCprepare(query_met_patient)

# Get beta values
beta_values <- assay(met_data_patient)
summary(beta_values[, 1])
mean(beta_values[, 1], na.rm = TRUE)

# Plot density of beta values
plot(density(beta_values[, 1], na.rm = TRUE),
     main = "Methylation Beta Values (First Sample)",
     xlab = "Beta Value", col = "steelblue")

# Top hypo/hypermethylated CpGs
top_hyper <- head(sort(beta_values[, 1], decreasing = TRUE), 10)
top_hypo  <- head(sort(beta_values[, 1], decreasing = FALSE), 10)

cat("Top Hypermethylated CpG")
head(top_hyper, 2)
cat("Top Hypomethylated CpG")
head(top_hypo, 2)



# ------------------ EXPRESSION DATA ------------------

query_exp_patient <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  barcode = exp_patient_barcodes
)

GDCdownload(query_exp_patient)
exp_data_patient <- GDCprepare(query_exp_patient)

# Expression matrix
expr_matrix <- assay(exp_data_patient)
summary(expr_matrix[, 1])
mean(expr_matrix[, 1], na.rm = TRUE)

# Plot density of expression values
plot(density(expr_matrix[, 1], na.rm = TRUE),
     main = "Gene Expression (First Sample)",
     xlab = "Read Counts",
     col = "darkgreen",
     xlim = c(0, 200000))  

# Top expressed genes
cat("Top expressed genes")
head(sort(expr_matrix[, 1], decreasing = TRUE), 10)
cat("Less expressed genes")
head(sort(expr_matrix[, 1], decreasing = FALSE), 10)



# Preparing OV data for downstream analysis: Differential Expression Analysis
# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/extension.html

# TCGAbiolinks is providing a new design for the TCGAanalyze_DEA() function to perform differential expression analysis (DEA) 
# either contrasting two groups or multiple sets by providing a formula for contrast. 
# Limma pipeline was also added on top of the previously implemented one in EdgeR.

# log.trans boolean to perform log cpm transformation. Set to TRUE for log transformation
# voom boolean to voom transform data. Set to true to apply transformation
# trend boolean to perform limma-trend pipeline. Set to TRUE to go through limma-trend

# TCGAbatch_correction: Handle batch correction and lima-voom transformation
# This function calls ComBat from sva package to handle batch correction, display some plots for exploratory data analysis purposes. 
# More specifically, it shows a plot with black as a kernel estimate of the empirical batch effect density and red as the parametric.

# TCGAtumor_purity: Filter TCGA samples according to tumor purity
# Function takes TCGA barcodes as input and filter them according to values specified by user. 
# Default value for all of them is 0. Arguments are estimates using different measures imported from the TCGA consortium, 
# and they are as the following (Refer to [1]):
# D. Aran, M. Sirota, and A. J. Butte. Systematic pan-cancer analysis of tumour purity. Nat Commun, 6:8971, 2015.

# estimate uses gene expression profiles of 141 immune genes and 141 stromal genes
# absolute which uses somatic copy-number data (estimations were available for only 11 cancer types)
# lump (leukocytes unmethylation for purity), which averages 44 non-methylated immune-specific CpG sites

cat("Download GTEx data available through the Recount2 project: 

Recount2 project has made gene count data available for the scientific community to download. 
Recount2 is an online resource consisting of RNA-seq gene and exon counts as well as coverage bigWig files for 2041 different studies. 

It is the second generation of the ReCount project. 
The raw sequencing data were processed with Rail-RNA as described in the recount2 paper.")

cat("Calculate stemness score with TCGAanalyze_Stemness")

cat( "A stemness score reflects the degree to which a tumor exhibits features of stem cells — 
such as self-renewal, dedifferentiation, and high plasticity — 
which are often associated with cancer aggressiveness, poor prognosis, and therapy resistance.

These scores were originally derived using machine learning models trained on:

Expression data from pluripotent stem cells.

Methylation profiles from stem cells vs. differentiated tissues" )

# stemness_scores <- TCGAanalyze_Stemness(
#  data = exp_data,
#  stemSig = "mRNAsi",   # or "mDNAsi" for methylation-based
# )

# stemness_scores <- TCGAanalyze_Stemness(
#  data = exp_data,
#  stemSig = "mDNAsi",   # or "mDNAsi" for methylation-based
# )



# Transcriptome Profiling
# miRNA expression

print("Download miRNA expression")

# Define a TCGA patient barcode (first 3 segments = patient ID)
patient_barcode <- "TCGA-WR-A838"

# Step 1: Query miRNA expression data for one patient
query.mirna <- GDCquery(
  project = "TCGA-OV", 
  experimental.strategy = "miRNA-Seq",
  data.category = "Transcriptome Profiling", 
  data.type = "miRNA Expression Quantification",
  barcode = patient_barcode
)

# Step 2: Download the data
GDCdownload(query.mirna)

# Step 3: Prepare the data
mirna_data <- GDCprepare(query.mirna)

# Step 4: Check class and view content accordingly
cat("Class of returned object:", class(mirna_data), "\n")

if ("SummarizedExperiment" %in% class(mirna_data)) {
  cat("It's a SummarizedExperiment\n")
  cat("Expression matrix (assay):\n")
  print(head(assay(mirna_data)))
  cat("Sample metadata (colData):\n")
  print(head(colData(mirna_data)))
  cat("miRNA annotations (rowData):\n")
  print(head(rowData(mirna_data)))
  
} else if ("data.frame" %in% class(mirna_data) || "tbl_df" %in% class(mirna_data)) {
  cat("It's a data.frame or tibble\n")
  cat("Here are the first few rows:\n")
  print(head(mirna_data))
  
} else {
  cat("Unknown object type. Here’s a structure preview:\n")
  str(mirna_data)
}


print("Download Copy Number Variation")

# Define a TCGA patient barcode (first 3 segments = patient ID)
patient_barcode <- "TCGA-WR-A838"

# Step 1: Query CNV expression data for one patient
query.cv <- GDCquery(
  project = "TCGA-OV", 
  data.category = "Copy Number Variation",
  data.type = "Copy Number Segment",
  barcode = patient_barcode
)

# Step 2: Download the data
GDCdownload(query.cv)

# Step 3: Prepare the data
cv_data <- GDCprepare(query.cv)

# Step 4: Check class and view content accordingly
cat("Class of returned object:", class(cv_data), "\n")

if ("SummarizedExperiment" %in% class(cv_data)) {
  cat("It's a SummarizedExperiment\n")
  cat("Expression matrix (assay):\n")
  print(head(assay(cv_data)))
  cat("Sample metadata (colData):\n")
  print(head(colData(cv_data)))
  cat("miRNA annotations (rowData):\n")
  print(head(rowData(cv_data)))
  
} else if ("data.frame" %in% class(cv_data) || "tbl_df" %in% class(cv_data)) {
  cat("It's a data.frame or tibble\n")
  cat("Here are the first few rows:\n")
  print(head(cv_data))
  
} else {
  cat("Unknown object type. Here’s a structure preview:\n")
  str(cv_data)
}


# https://bioconductor.org/packages/release/bioc/vignettes/TCGAbiolinks/inst/doc/download_prepare.html

library(TCGAbiolinks)

# Define patient barcode
patient_barcode <- "TCGA-WR-A838"

# Step 1: Query gene-level CNV data for the patient
query.cnv <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Gene Level Copy Number",
  access = "open",
  barcode = patient_barcode
)

# Step 2: Get results and view available files
results.cnv <- getResults(query.cnv)
print(results.cnv[, c("file_id", "cases", "analysis_workflow_type")])

# Step 3: Keep only files where 'cases' do NOT contain semicolon (i.e., single sample)
results.cnv <- results.cnv[!grepl(";", results.cnv$cases), ]

# Optional: Keep only tumor samples (e.g., ending in "01A")
results.cnv <- results.cnv[grepl("01A", results.cnv$cases), ]

# Step 4: Remove duplicate case entries
filtered_file_ids <- results.cnv[!duplicated(results.cnv$cases), "file_id"]

# Step 5: Subset the query to only include those files
query.cnv$results[[1]] <- query.cnv$results[[1]][query.cnv$results[[1]]$file_id %in% filtered_file_ids, ]

# Step 6: Download the files
GDCdownload(query.cnv)

# Step 7: Prepare the data
cnv_data <- GDCprepare(query.cnv)

# Step 8: Preview
head(cnv_data)



print("Allele-specific Copy Number")

print("Find the barcodes in the OV data")

# Define the TCGA project
project_id <- "TCGA-OV"

# Step 1: Query Allele-specific Copy Number Segment data
query <- GDCquery(
  project = project_id,
  data.category = "Copy Number Variation",
  data.type = "Allele-specific Copy Number Segment",
  access = "open"
)

# Step 2: Extract query results
results <- getResults(query)

# Step 3: Safety check — did we get results?
if (nrow(results) == 0) {
  stop("No files found for the specified project and data type.")
}

# Step 4: Print summary
cat("Files returned:\n")
# print(results[, c("file_id", "cases", "submitter_id", "analysis_workflow_type")])

# Step 5: Remove duplicated samples (by keeping one per sample barcode)
results$sample_barcode <- substr(results$cases, 1, 16)
results <- results[!duplicated(results$sample_barcode), ]

# Step 6: Filter the original query to just these non-duplicated files
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$file_id %in% results$file_id,
]

# Step 7: Download filtered files
GDCdownload(query)

# Step 8: Prepare the data safely (no duplicate samples now)
data <- GDCprepare(query)

# Step 9: Preview in console
cat("Preview of allele-specific CNV data:\n")
print(head(data, 3))

# Step 10: Display interactive table
# if (inherits(data, "SummarizedExperiment")) {
#  datatable(
#    as.data.frame(assay(data)),
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = TRUE
#  )
# } else if (is.data.frame(data)) {
#  datatable(
#    data,
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = FALSE
#  )
# } else {
#  cat("Data is not a SummarizedExperiment or data.frame. Class is:\n")
#  print(class(data))
# }


print("Extract CNA data for the patient TCGA-WR-A838")

# Define TCGA patient barcode
patient_barcode <- "TCGA-WR-A838"

# Step 1: Query Allele-specific CNV data for the patient
query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Copy Number Variation",
  data.type = "Allele-specific Copy Number Segment",
  access = "open",
  barcode = patient_barcode
)

# Step 2: Extract results
results <- getResults(query)

# Step 3: Check and filter duplicates (some files have tumor/normal pairs)
if (nrow(results) == 0) {
  stop("No files found for the specified patient.")
}

# Print what's available
cat("Files returned for patient:\n")
print(results[, c("file_id", "cases", "submitter_id", "analysis_workflow_type")])

# Step 4: Extract one file per sample (use first 16 characters = sample barcode)
results$sample_barcode <- substr(results$cases, 1, 16)
results <- results[!duplicated(results$sample_barcode), ]

# Step 5: Subset the original query to the selected file IDs
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$file_id %in% results$file_id,
]

# Step 6: Download only the necessary files
GDCdownload(query)

# Step 7: Prepare the data (no duplication error now)
data <- GDCprepare(query)

# Step 8: Preview result
cat("Preview of allele-specific CNV data for", patient_barcode, ":\n")
print(head(data, 3))

# Step 9: Display in an interactive table
# if (inherits(data, "SummarizedExperiment")) {
#  datatable(
#    as.data.frame(assay(data)),
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = TRUE
#  )
# } else if (is.data.frame(data)) {
#  datatable(
#    data,
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = FALSE
#  )
# } else {
#  cat("Unexpected data format:\n")
#  print(class(data))
# }


print("Extract mutation data for all OV patients")

# Define the TCGA project
project_id <- "TCGA-OV"

# Step 1: Query mutation (MAF) data
query <- GDCquery(
  project = project_id,
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

# Step 2: Extract query results
results <- getResults(query)

# Step 3: Check what columns are available
cat("Columns in results:\n")
print(names(results))

# Step 4: Preview available data
cat("Mutation file preview:\n")
print(head(results))

# Step 5: Remove duplicated samples using the correct column
# Usually 'cases' or 'case_submitter_id' or 'submitter_id' contains the barcode
if ("cases" %in% names(results)) {
  results$sample_barcode <- substr(results$cases, 1, 16)
} else if ("case_submitter_id" %in% names(results)) {
  results$sample_barcode <- substr(results$case_submitter_id, 1, 16)
} else {
  stop("Could not find a barcode column to deduplicate samples.")
}

results <- results[!duplicated(results$sample_barcode), ]

# Step 6: Filter the original query object
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$file_id %in% results$file_id,
]

# Step 7: Download the data
GDCdownload(query)

# Step 8: Prepare the data (MAF)
mutation_data <- GDCprepare(query)

# Step 9: Display preview
cat("Preview of mutation data:\n")
print(head(mutation_data, 3))

# Step 10: Display interactive table
# if (is.data.frame(mutation_data)) {
#  datatable(
#    mutation_data,
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = FALSE
#  )
# } else {
#  cat("Mutation data is not a data.frame. Class:\n")
#  print(class(mutation_data))
# }


# Define the TCGA project and patient barcode
patient_barcode <- "TCGA-WR-A838"

# Step 1: Query mutation (MAF) data for a specific patient
query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  barcode = patient_barcode,
  access = "open"
)

# Step 2: Extract query results
results <- getResults(query)

# Step 3: Check what columns are available
cat("Columns in results:\n")
print(names(results))

# Step 4: Preview available data
cat("Mutation file preview:\n")
print(head(results))

# (Optional) Step 5: Print matched barcodes
if ("cases" %in% names(results)) {
  results$sample_barcode <- substr(results$cases, 1, 16)
  print(unique(results$sample_barcode))
}

# Step 6: Filter the query to valid file_ids (just a safety step)
query$results[[1]] <- query$results[[1]][
  query$results[[1]]$file_id %in% results$file_id,
]

# Step 7: Download the data (skips if already downloaded)
GDCdownload(query)

# Step 8: Prepare the mutation (MAF) data
mutation_data <- GDCprepare(query)

# Step 9: Show preview
cat("Preview of mutation data:\n")
print(head(mutation_data, 3))

# Step 10: Display in interactive DT table
# if (is.data.frame(mutation_data)) {
#  datatable(
#    mutation_data,
#    options = list(scrollX = TRUE, pageLength = 10),
#    rownames = FALSE
#  )
# } else {
#  cat("Mutation data is not a data.frame. Class:\n")
#  print(class(mutation_data))
# }

# Only first 50 to make render faster
# datatable(mutation_data[1:10,],
#          filter = 'top',
#          options = list(scrollX = TRUE, keys = TRUE, pageLength = 5), 
#          rownames = FALSE)



# Load required libraries

library(maftools)

# Step 1: Query mutation data using the correct workflow type
query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking",
  access = "open"
)

# Step 2: Download the data
GDCdownload(query)

# Step 3: Prepare the data as a data.frame (likely a tibble)
maf_df <- GDCprepare(query)

# Step 4: Convert to a maftools MAF object
maf_object <- read.maf(maf = as.data.frame(maf_df))

# Step 5: Get a summary
getSampleSummary(maf_object)

# Optional: Visualize MAF stats
plotmafSummary(maf = maf_object)

# Set plot height for Jupyter output (adjust this as needed)
options(repr.plot.width = 12, repr.plot.height = 8)

# Oncoplot for top 10 mutated genes
oncoplot(maf = maf_object, top = 10, removeNonMutated = TRUE)

# Transitions and transversions
titv_result <- titv(maf = maf_object, plot = FALSE, useSyn = TRUE)

# Plot TiTv summary
plotTiTv(res = titv_result)

# plot <- mafSurvival(
#  maf = maf_object,
#  genes = "TP53",
#  time = 'time',
#  Status = 'Overall_Survival_Status',
#  isTCGA = TRUE
# )

#  maf,
#  genes = NULL,
#  samples = NULL,
#  clinicalData = NULL,
#  time = "Time",



# survdiff performs a log-rank test (which compares survival curves), it doesn't explicitly fit a Cox proportional hazards model. 
# If you wanted to assess hazard ratios and adjust for covariates, you would use the coxph() function from the survival package.

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


