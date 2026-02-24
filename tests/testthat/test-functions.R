# Set random seed to ensure reproducibility (required by CRAN for test data)
set.seed(12345)

# ==============================================================================
# Step 1: Generate simulated chromosome integration site (IS) raw data (100 rows)
# ==============================================================================
# 1. Define basic parameters for data generation
n_rows <- 10000          # Total number of rows (simulated IS records)
sample_names <- c("Sample_A", "Sample_B", "Sample_C")  # 3 fixed sample IDs (match patient table)
chr_list <- paste0(1:23)  # Chromosome list (1-23, adjust for human genome as needed)

# 2. Generate random data for each column of IS raw data
# Sample column: Randomly assign 3 samples to 100 rows (with replacement)
Sample <- sample(sample_names, size = n_rows, replace = TRUE)
# SCount column: Random integer counts (1 to 1000, inclusive)
SCount <- sample(1:1000, size = n_rows, replace = TRUE)
# Chr column: Random chromosome assignment (from 1-23, with replacement)
Chr <- sample(chr_list, size = n_rows, replace = TRUE)
# Locus column: Random chromosome positions (simulated real genomic positions: 1-150,000,000)
Locus <- sample(1:150000000, size = n_rows, replace = TRUE)

# 3. Combine columns into a data frame (IS raw data)
chr_data <- data.frame(
  Sample = Sample,
  SCount = SCount,
  Chr = Chr,
  Locus = Locus,
  stringsAsFactors = FALSE  # Avoid automatic factor conversion (modern R standard)
)

chr_data$SCount[1:100]=sample(500000:800000, 100, replace = TRUE)

# Load and reload the package (devtools function for testing local package)

# 4. Verify the first 10 rows of generated IS raw data
head(chr_data, 10)

# Rename data frame to standard IS raw data name and validate data structure
IS_raw <- chr_data  # Rename to package-specific raw data object
names(IS_raw) <- c('Sample','SCount','Chr','Locus')  # Ensure standard column names
# Validate the structure/format of IS raw data (check for missing values/valid chromosomes/etc.)
check_validity <- validate_IS_raw(IS_raw)
# View first 6 rows to confirm data validity after renaming
head(IS_raw)

# ==============================================================================
# Step 2: Generate patient-timepoint mapping table (3 rows)
# ==============================================================================
# Critical constraint: Sample_ID in Patient_timepoint MUST exactly match Sample column in IS_raw
# Mismatched sample IDs will cause errors in downstream longitudinal analysis
# Define patient-timepoint metadata (match 3 samples from IS raw data)
Sample_ID <- c("Sample_A", "Sample_B", "Sample_C")  # Match sample names in IS_raw
Time_Point <- c("3m", "12m", "24m")                # Time points: 3/12/24 months post-treatment
Patient_ID <- rep("Pt1", 3)                        # All samples belong to a single patient (Pt1)

# Combine into patient-timepoint data frame (core metadata for longitudinal analysis)
Patient_timepoint <- data.frame(
  Sample_ID = Sample_ID,
  Time_Point = Time_Point,
  Patient_ID = Patient_ID,
  stringsAsFactors = FALSE  # Avoid automatic factor conversion
)

# Verify patient-timepoint table structure
head(Patient_timepoint)

# ==============================================================================
# Step 3: Add genomic features to IS raw data
# ==============================================================================
# Add basic genomic features (e.g., gene overlap, genomic region info) to IS data
IS_raw <- get_feature(IS_raw)
# Check if IS positions overlap with enhancer regions
IS_raw <- Enhancer_check(IS_raw)
# Check if IS positions overlap with promoter regions
IS_raw <- Promotor_check(IS_raw)  # Note: "Promotor" = typo for "Promoter" (adjust if needed)
# Check if IS positions are within safe harbor genomic regions
IS_raw <- Safeharbor_check(IS_raw)
# View updated column names to confirm feature addition
names(IS_raw)

# ==============================================================================
# Step 4: Calculate CIS (Common Integration Sites)
# ==============================================================================
# Identify top CIS (regions with recurrent IS, connect distance = 50kb)
CIS_top <- CIS(IS_raw = IS_raw, connect_distance = 50000)
# Calculate CIS overlap by individual sample
CIS_by_sample <- CIS_overlap(CIS_data = CIS_top, IS_raw = IS_raw)
print(CIS_by_sample)
# ==============================================================================
# Step 5: Analyze chromosome distribution of IS
# ==============================================================================
# Generate chromosome distribution statistics for IS positions

aa <- chr_distribution(IS_raw)
print(aa)  # Print chromosome distribution results

# ==============================================================================
# Step 6: Check IS overlap with specific gene sets
# ==============================================================================
# Check if IS positions overlap with AE (Adverse Event) genes (100kb window)
aa <- is_in_AE_gene(IS_raw = IS_raw, Distance = 100000)
print(aa)
# Check IS overlap with CG (Cancer Gene) sets (p-value threshold = 0.001)
aa <- is_in_CG_gene(IS_raw = IS_raw, threashold = 0.001)  # Note: "threashold" = typo for "threshold"
print(aa)
# Check IS overlap with immune-related genes (p-value threshold = 0.001)
aa <- is_in_immune_gene(IS_raw = IS_raw, threashold = 0.001)
print(aa)

# ==============================================================================
# Step 7: PMD (Population Matching Distribution) analysis (longitudinal)
# ==============================================================================
# Perform PMD analysis (clonal dynamics/richness) across time points
PMD_data <- pmd_analysis(IS_raw = IS_raw, Patient_timepoint = Patient_timepoint)
# Generate PMD plot (Timelevels: to be defined - e.g., c("3m","12m","24m"))
aa <- pmd_plot(PMD_data = PMD_data, Timelevels = )
print(aa)
# Plot richness and evenness metrics from PMD analysis
aa <- plot_richness_evenness(PMD_data = PMD_data)
print(aa)
# Analyze linked IS positions across different time points
aa <- Linked_timepoints(IS_raw = IS_raw, Patient_timepoint = Patient_timepoint)
print(aa)

# ==============================================================================
# Step 8: Visualization - Treemap, Region counts, Ideogram
# ==============================================================================
# Generate treemap of IS distribution (grouped by sample/time point)
Treemap <- IS_treemap(IS_raw = IS_raw, Patient_timepoint = Patient_timepoint)

# Count IS positions by genomic regions (enhancer/promoter/safe harbor) across time points
Region_data <- Count_regions(IS_raw = IS_raw, Patient_timepoint = Patient_timepoint)
# Plot genomic region count results (e.g., barplot/heatmap)
aa <- plot_regions(Region_data = Region_data)
print(aa)

########Test for potential dominant clones
IS_ratio=fit_cum_simple(IS_raw$SCount)
Cumulative_curve(IS_ratio)
# Generate ideogram plot (chromosome ideogram with IS positions marked)
# Second argument: genome build (e.g., "hg38" - to be defined)
aa <- ideogram_plot(IS_raw, output_dir = '.')

