W# Server Locations
`/home/projects/AD_CERC`

# Local Locations
`/home/reedrich/Wrighton-Lab/Projects/Jorge_AD_CERC`

## Gene database
`/home/projects/AD_CERC/DRAM_2500_contigs/DRAM_combined/DRAM_merged/genes.fna`


# Software versions
Prodigal v2.6.3
cd-hit v4.6
FastQC v0.11.2
sickle v1.33
Bowtie2 v2.4.5
samtools v1.9
BBMap version 38.89

# Building a gene database
Reference: [Josue's Burned analysis](https://github.com/jrr-microbio/forest_soil_burn_treatment/tree/main/gene_resolved)
- Call genes using DRAM on the >2500 scaffolds
- cluster all genes.faa
- map reads back to those clusters
- use those abundance to find patterns


1) Cluster genes from DRAM (`/home/projects/AD_CERC/DRAM_2500_contigs/DRAM_combined/DRAM_merged/genes_2500_DRAM_all_merged.fna)
Script: `1_cdhit_est_cluster_genes_2500_DRAM-out.sh`


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=128gb
#SBATCH --job-name=cd-hit
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reed.woyda@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

cd /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/

#====== CD-HIT version 4.6 (built on Nov  6 2014) ======
cd-hit-est -i /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/genes_2500_DRAM_all_merged.fna -o /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/clustered_genes_2500_DRAM_all_merged.fna -c 0.95 -aS 0.80 -M 128000 -T 20 -d 0 -B 1 -bak 1
```


2) Trim raw reads: `3_trimming_raw_reads.sh`

Notes: 
- Running fastQC on a single sample to see if I need to trim the reads - DONE - QC NEEDED
- Need to trim - DONE


```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reed.woyda@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

# Set the input directory containing your FASTQ files 
input_dir="/home/ORG-Data-2/AD_CERC/metaG_CU_sequencing_16Mar2020" 

# Set the output directory 
output_dir="/home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/processed_reads"


# Iterate through R1 FASTQ files 
for r1_file in $input_dir/*_R1_001.fastq.gz; do 
	# Extract the base name of the R1 file 
	base_name=$(basename "$r1_file" _R1_001.fastq.gz) 
	
	# Construct corresponding R2 file name 
	r2_file="${base_name}_R2_001.fastq.gz" 
	
	
	# Run sickle command 
	sickle pe -f "$r1_file" -r "$input_dir/${base_name}_R2_001.fastq.gz"  -t sanger -o "$output_dir/${base_name}_R1_trimmed.fastq.gz"  -p "$output_dir/${base_name}_R2_trimmed.fastq.gz"  -s "$output_dir/${base_name}_singles.fastq"  
	
	# Run fastqc on processed reads 
	fastqc -t 20 "$output_dir/${base_name}_R1_trimmed.fastq.gz" "$output_dir/${base_name}_R2_trimmed.fastq.gz" -o "$output_dir"
done

```




3) Build bowtie database of the clustered genes and then map each set of reads back to this clustered gene database
Script: `3_read_mapping_2500_genes_DRAM_output_clustered.sh`

Notes: 
- Running fastQC on a single sample to see if I need to trim the reads - DONE - NEEDED
- Trimming first with sickle - DONE - GOOD
- Create database using bowtie on the concatenated genes.fna from DRAM
- Map the trimmed reads to the database using the params
	- -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 20
- Convert SAM to BAM
- Filter using Reformat.sh at 95% ID


```

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=20
#SBATCH --time=14-00:00:00
#SBATCH --mem=300gb
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=reed.woyda@colostate.edu
#SBATCH --partition=wrighton-hi,wrighton-low

# pre-processing:
# concatenate MAGs into single file for mapping reference.
#Select medium-high quality genomes
# fastq_trimmed_list.txt is a list of trimmed read names.

cd /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/

#create bowtie2 index with a concatenated genomes fasta.
bowtie2-build /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/clustered_genes_2500_DRAM_all_merged.fna gene_DB --large-index --threads 20 

input_file="/home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/scripts/fastq_trimmed_list.txt"

# Check if the input file exists and is a regular file
if [[ -f "$input_file" ]]; then
    while read -r r1_file && read -r r2_file; do
        # Ensure that both R1 and R2 files are provided
        if [[ -z "$r1_file" || -z "$r2_file" ]]; then
            echo "Error: Incomplete FASTQ file pair."
            continue
        fi
        
		filename=$(basename "$r1_file")

        # Get the element name from the R1 file (adjust as needed based on your file naming convention)
        element=$(echo "$filename" | sed -E 's/_S[0-9]+.*_R1_trimmed.fastq//')


        echo "begin bowtie2 '$element' mapping"
        mkdir "$element"_mapping
        cd "$element"_mapping

        # run bowtie2 with indexed genes.
        bowtie2 -D 10 -R 2 -N 1 -L 22 -i S,0,2.50 -p 20 -x /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/gene_DB -S "$element".sam -1 "$r1_file" -2 "$r2_file"

        echo "begin samtools"
        # convert sam to bam
        samtools view -@ 20 -bS "$element".sam > "$element".bam

        # this will then give you the actual 95% ID mapping.
        reformat.sh -Xmx100g idfilter=0.95 pairedonly=t primaryonly=t in="$element".bam out="$element"_mapped_95ID.FILTERED.bam 

        # sort the bam file
        samtools sort -T "$element".95ID.sorted -o "$element".95ID.sorted.bam "$element"_mapped_95ID.FILTERED.bam -@ 20 

        mv "$element".95ID.sorted.bam /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/read_mapping_genes_output/

        echo "Finished mapping and samtools for '$element' in list"
        # clean up the files
        # rm -r "$element".sam
        # rm -r "$element".bam
        # rm -r "$element"_mapped_95ID.FILTERED.bam
        cd /home/projects/AD_CERC/gene_resolved_analysis_RW_08312023/DRAM_unclustered_genes/
        # move back up to root and restart the process.
    done < "$input_file"
else
    echo "Input file '$input_file' not found."
fi

```

4) Read Mapping and coverM
script: `4_read_mapping_coverM_gene_resolved.sh`

Notes:
- CoverM reads_per_base at 95%
- CoverM coverage at 75%
- CoverM proper pairs only

```
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=30
#SBATCH --time=14-00:00:00
#SBATCH --mem=200gb
#SBATCH --partition=wrighton-hi,wrighton-low

coverm contig --proper-pairs-only -m reads_per_base --bam-files C1A_DALY7_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file C1A_coverm_reads_per_base.txt &> C1A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files C3A_DALY8_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file C3A_coverm_reads_per_base.txt &> C3A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files C7A_DALY9_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file C7A_coverm_reads_per_base.txt &> C7A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files C8A_DALY10_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file C8A_coverm_reads_per_base.txt &> C8A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files I1A_DALY11_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file I1A_coverm_reads_per_base.txt &> I1A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files I3A_DALY12_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file I3A_coverm_reads_per_base.txt &> I3A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files I7A_DALY13_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file I7A_coverm_reads_per_base.txt &> I7A_reads_per_base_stats.txt
coverm contig --proper-pairs-only -m reads_per_base --bam-files I8A_DALY14_mapped_95ID.FILTERED.sorted.bam --threads 20 --min-read-percent-identity-pair 0.95 --min-covered-fraction 0 --output-format dense --output-file I8A_coverm_reads_per_base.txt &> I8A_reads_per_base_stats.txt

coverm contig -b C1A_DALY7_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file C1A_coverm_min75.txt &> C1A_min75_stats.txt
coverm contig -b C3A_DALY8_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file C3A_coverm_min75.txt &> C3A_min75_stats.txt
coverm contig -b C7A_DALY9_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file C7A_coverm_min75.txt &> C7A_min75_stats.txt
coverm contig -b C8A_DALY10_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file C8A_coverm_min75.txt &> C8A_min75_stats.txt
coverm contig -b I1A_DALY11_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file I1A_coverm_min75.txt &> I1A_min75_stats.txt
coverm contig -b I3A_DALY12_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file I3A_coverm_min75.txt &> I3A_min75_stats.txt
coverm contig -b I7A_DALY13_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file I7A_coverm_min75.txt &> I7A_min75_stats.txt
coverm contig -b I8A_DALY14_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0.75 --threads 20 --output-file I8A_coverm_min75.txt &> I8A_min75_stats.txt

coverm contig --proper-pairs-only -m count --bam-files C1A_DALY7_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file C1A_counts.txt  &> C1A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files C3A_DALY8_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file C3A_counts.txt  &> C3A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files C7A_DALY9_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file C7A_counts.txt  &> C7A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files C8A_DALY10_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file C8A_counts.txt  &> C8A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files I1A_DALY11_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file I1A_counts.txt  &> I1A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files I3A_DALY12_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file I3A_counts.txt  &> I3A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files I7A_DALY13_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file I7A_counts.txt  &> I7A_coverM_counts.txt
coverm contig --proper-pairs-only -m count --bam-files I8A_DALY14_mapped_95ID.FILTERED.sorted.bam --min-covered-fraction 0 -t 20 --output-file I8A_counts.txt  &> I8A_coverM_counts.txt

```

5) Merge coverM outputs in R

script: `5_coverM_merging.R`

Notes:
- I want to keep hits that were >97% ID, >=3x depth, and >=75% coverage.
- 

```
## ---------------------------
## Script merge_coverM_outputs
## 
## Purpose merge MetaG coverM outputs
##
## Author: Reed Woyda
##
## Date Created: September 15, 2023
##
## Copyright (c) Reed Woyda, 2023
## Email: reed.woyda@colostate.edu
## ---------------------------
library(tidyverse)

# Where Reed started and modified until line 43
# /home/reedrich/Wrighton-Lab/Projects/Tiffany_Weir_Blueberry/CoverM_outputs
coverage_new = read.csv("/home/reedrich/Wrighton-Lab/Projects/Jorge_AD_CERC/Scripts/coverm_reads_per_base.txt", sep = "\t", row.names=1, header = TRUE)
depth_new = read.table("/home/reedrich/Wrighton-Lab/Projects/Jorge_AD_CERC/Scripts/coverm_min75.txt", sep = "\t", row.names=1, header = TRUE)
counts_new = read.table("/home/reedrich/Wrighton-Lab/Projects/Jorge_AD_CERC/Scripts/trimmed_mean.txt", sep = "\t", row.names=1, header = TRUE)

# Remove the specified text from column names
colnames(counts_new) <- gsub("\\.95ID\\.sorted\\.Trimmed\\.Mean", "", colnames(counts_new))
colnames(coverage_new) <- gsub("\\.95ID\\.sorted\\.Reads\\.per\\.base", "", colnames(coverage_new))
colnames(depth_new) <- gsub("\\.95ID\\.sorted\\.Mean", "", colnames(depth_new))

# Multiply all values in coverage_new by 151
coverage_new <- coverage_new * 151

# Create a new data frame with the same structure as counts_new, initialized with 0s
trimmed_mean_final <- as.data.frame(matrix(0, nrow = nrow(counts_new), ncol = ncol(counts_new)))

# Iterate through all rows and columns
for (row_index in 1:nrow(counts_new)) {
  for (col_index in 1:ncol(counts_new)) {
    # Check conditions for the current element
    if (coverage_new[row_index, col_index] >= 3 && depth_new[row_index, col_index] > 0) {
      # If conditions are met, assign the corresponding value from counts_new
      trimmed_mean_final[row_index, col_index] <- counts_new[row_index, col_index]
    }
  }
}

# Set row names
rownames(trimmed_mean_final) <- row.names(coverage_new)

# Set column names
colnames(trimmed_mean_final) <- colnames(counts_new)



write.csv(trimmed_mean_final, file="AD_CERC_rel_abunds_95ID_75cov_3xdepth.csv")

```

# Differential gene abundance analysis

6) Merging DRAM functions and MaAslin

script: `6_DRAM_merging_functions_and_MaAslin.R`


```
library(dplyr)
library(tidyr)
library(tibble)
library(janitor)

if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")

# Install 'Maaslin2' from Bioconductor
BiocManager::install("Maaslin2", update = FALSE, ask = FALSE, force = TRUE)


# read in annotation ids
ids= read.delim("DRAM_2500_merged_annotations.tsv", header=TRUE, sep = "\t", fill = T)
# create new variable "primary", selecting one id per gene. priority cazy > merops > kegg > vogdb
ids$cazy_hits_summ=as.character(ids$cazy_hits) #bring in cazy id
ids$peptidase_hits_summ=as.character(ids$peptidase_id) #bring in peptidase id
ids$kegg_hits_summ=as.character(ids$kegg_id) #bring in ko id
ids$vogdb_hits_summ=as.character(ids$vogdb_id) #bring in ko id
ids$gene_ids=as.character(ids$gene) #bring in scaffold id

ids$cazy_hits_summ <- trimws(ids$cazy_hits_summ)
ids$peptidase_hits_summ <- trimws(ids$peptidase_hits_summ)
ids$kegg_hits_summ <- trimws(ids$kegg_hits_summ)
ids$vogdb_hits_summ <- trimws(ids$vogdb_hits_summ)

ids$primary <- character(nrow(ids))  # Create an empty character vector for the primary column

for (i in 1:nrow(ids)) {
  if (!is.na(ids$cazy_hits_summ[i]) && ids$cazy_hits_summ[i] != "") {
    ids$primary[i] <- ids$cazy_hits_summ[i]
  }  else if (!is.na(ids$peptidase_hits_summ[i]) && ids$peptidase_hits_summ[i] != "") {
    ids$primary[i] <- ids$peptidase_hits_summ[i]
  } else if (!is.na(ids$kegg_hits_summ[i]) && ids$kegg_hits_summ[i] != "") {
    ids$primary[i] <- ids$kegg_hits_summ[i]
  } else if (!is.na(ids$vogdb_hits_summ[i]) && ids$vogdb_hits_summ[i] != "") {
    ids$primary[i] <- ids$vogdb_hits_summ[i]
  }
    else {
    ids$primary[i] <- "none"
  }
}


ids_only = ids %>% select(gene_ids,primary) #if else statement to give CAZY, if not MEROPS, if not KEGG, if not vOGDB, if not "none".
rm(ids)#just clear the env.
annotation_ids = ids_only #make this cleaner name
rm(ids_only) #just clear the env.

# pulled DRAM module IDs from the genome summary form in the DRAM Github.
# RW just combined all metab summary sheets into one
module_info = read.delim("metab_summary_all.csv", sep = "\t", header=T, na.strings = c("","NA"))
joined_annots_and_module=left_join(annotation_ids, module_info,by=c("primary"="gene_id"), relationship = "many-to-many") #merge the two files.
sum(duplicated(joined_annots_and_module$gene_ids)) #check no duplicates just in case


# Table 1:
filtered_module_info <- module_info[module_info$gene_id %in% joined_annots_and_module$primary, ]

# Perform an inner join between the filtered "module_info" and "filtered_joined_annots" based on the "gene_id" and "primary" columns
merged_metadata <- merge(filtered_module_info, joined_annots_and_module, by.x = "gene_id", by.y = "primary")

# Create the "gene_metadata" table with the specified columns
gene_metadata <- data.frame(
  contig = merged_metadata$gene_ids,             # Use the "gene_ids" column from "filtered_joined_annots"
  gene_id = merged_metadata$gene_id,            # Use the "gene_id" column from "filtered_module_info"
  gene_description = merged_metadata$gene_description.x,
  module = merged_metadata$module.x,
  header = merged_metadata$header.x,
  subheader = merged_metadata$subheader.x
)

write.csv(gene_metadata, "Gene_metadata.csv")

# Table 2:
# read in abundances
abunds = read.csv("CoverM_abundances.csv")

# Merge the "abunds" and "joined_annots_and_module" data frames by matching the first columns
merged_table <- merge(abunds, joined_annots_and_module, by.x = "X", by.y = "gene_ids", all.x = TRUE)

# Create a new data frame with the desired columns
new_table <- merged_table[, c("primary", names(abunds)[2:9])]

# Remove rows where "gene_id" is "none"
maaslin_input <- subset(new_table, primary != "none")

# Remove rows with duplicate values in the first column
maaslin_input <- maaslin_input[!duplicated(maaslin_input$primary), ]

gene_ids_abunds <- maaslin_input


# Set row names to the first column values
rownames(maaslin_input) <- maaslin_input$primary

# Remove the first column
maaslin_input <- maaslin_input[, -1]

write.csv(maaslin_input, "maaslin_input_abunds.csv")



library(Maaslin2)

# RW bring in metadata
metadata <- read.table("metadata.csv", sep = ",", header=T, row.names=1)

# Remove text between pipes or parentheses in row names
rownames(maaslin_input) <- gsub("\\|.*?\\||\\(.*?\\)", "|", rownames(maaslin_input))

# Find row names longer than 10,000 bytes
long_row_names <- rownames(maaslin_input)[nchar(rownames(maaslin_input)) > 10000]

# Create a mapping table
mapping_table <- data.frame(
  original_row_names = long_row_names,
  keys = paste("key", seq_along(long_row_names), sep = "_")
)

# Replace long row names with keys
rownames(maaslin_input)[nchar(rownames(maaslin_input)) > 10000] <- mapping_table$keys


#And now run the same commands.
fit_data = Maaslin2(input_data = maaslin_input, 
                               input_metadata = metadata, 
                               min_prevalence = 0,
                               normalization  = "TSS",
                               transform = "LOG",
                               output         = "output_sample_time_Maaslin_RW_09142023", 
                                fixed_effects  = "Sample_Time")

### Collapsing to functional level
# use abunds table and gene_metadata table to collapse on the "module" in gene_metadata

# remove extra contig naming
gene_metadata$contig <- sub("^2500-contigs-renamed_", "", gene_metadata$contig)




# Merge the "abunds" and "gene_metadata" tables based on the "X" and "contig" columns
merged_data <- abunds %>%
  left_join(gene_metadata, by = c("X" = "contig"))

# Remove rows with NA values in the module column
merged_data <- merged_data %>%
  filter(!is.na(module))

# Remove rows with duplicate gene_id values
merged_data <- merged_data %>%
  distinct(gene_id, .keep_all = TRUE)

# Group by "module" and summarize the abundance values for all sample columns
result_table <- merged_data %>%
  group_by(module) %>%
  summarize(across(starts_with("C") | starts_with("I"), sum, na.rm = TRUE), .groups = "drop")

# Debugging: Print out the first few rows of the result_table
cat("First few rows of result_table:\n")
print(head(result_table))

# Debugging: Check for any NA values in the result_table
cat("NA values in result_table:\n")
print(sum(is.na(result_table)))

# Debugging: Check for any NA values in the module column
cat("NA values in module column:\n")
print(sum(is.na(result_table$module)))

# Debugging: Check for any NA values in sample columns
sample_cols <- names(result_table)[2:length(names(result_table))]
cat("NA values in sample columns:\n")
print(colSums(is.na(result_table[sample_cols])))

# Debugging: Check the structure of the result_table
cat("Structure of result_table:\n")
str(result_table)



# create input to maaslin2
maaslin_input_2 <- as.data.frame(result_table)

# Set row names to the first column values
row.names(maaslin_input_2) <- maaslin_input_2$module

# Remove the first column
maaslin_input_2 <- maaslin_input_2[, -1]




#And now run the same commands.
fit_data = Maaslin2(input_data = maaslin_input_2, 
                    input_metadata = metadata, 
                    min_prevalence = 0,
                    normalization  = "TSS",
                    transform = "LOG",
                    output         = "Maaslin2_output_Sample_Time_MODULE_Second_Try", 
                    fixed_effects  = "Sample_Time")

```
