library(tidyverse)
library(fs)

# Input parameters
organism1 = 'Past'
organism2 = 'Smic'
organisms = paste(organism1, organism2, sep='_')
cat("ORGANISMS: ", organisms, "\n", sep='')

# derived constants
KALLISTO_ABUNDANCE_PATH = paste("/results_Kallisto_", organisms, "-reefGenomics/abundance.tsv", sep="")
BWA_SALMON_QUANT_SF_PATH = paste("/results_bwa_Salmon_", organisms, "-reefGenomics/salmon_quant/quant.sf", sep="")
BWA_SALMON_REGEXP = paste('*', BWA_SALMON_QUANT_SF_PATH, sep='')

## Function to extract TPM from each file
quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$Name, data$TPM))
}

# Remove the path components from the column names of the tibble
# x: a vector of column names
remove_quant_path = function(x) {
   return (unlist(lapply(x, function(s) {
     if (is.character(s) && endsWith(s, "quant.sf")) {
       return (unlist(strsplit(s, '/'))[1])
     }
     return(s)
   })));
}

## Extract Salmon calculated quants
STAR_SALMON_REGEXP = "*/*salmon_quant/quant.sf"
STAR_SALMON_QUANT_SF_PATH = paste("/results_STAR_Salmon/salmon_quant/quant.sf", sep="")
extract_salmon_quants <- function(analysis_dir, outdir, regexp=STAR_SALMON_REGEXP) {
  message("\nExtract Salmon calculated quants")
  message(paste("dir: [", analysis_dir, "]", sep=''))
  message(paste("regexp: [", regexp, "]", sep=''))
  salmon_files <- fs::dir_ls(analysis_dir, regexp=regexp, recurse=T)
  #message("SALMON FILES FOUND: ")
  #print(salmon_files)
  #message("\n**** END SALMON FILES FOUND ****\n")
  # get the file list and pipe it into our extractor function
  salmon_df <- salmon_files %>%
    set_names(.) %>%
    map_dfr(quant_extractor, .id = "file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    # strip the input directory part out of the column name
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    # remove the last path part from the name
    rename_with(~ remove_quant_path(.)) %>%
    rename(gene_id = `data$Name`)

  salmon_df_org1 <- salmon_df %>% filter(grepl(organism1, gene_id))
  salmon_df_org2 <- salmon_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/STAR_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/STAR_Salmon_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/STAR_Salmon_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  write_csv(salmon_df, file=merged_file)
  write_csv(salmon_df_org1, file=org1_file)
  write_csv(salmon_df_org2, file=org2_file)
}

## Function to extract TPM from each file
rsem_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$gene_id, data$TPM))
}

## Extract Bowtie2/RSEM quants
extract_bowtie2rsem_quants <- function(analysis_dir, outdir) {
  message("Extract Bowtie2/RSEM calculated quants")
  rsem_files <- fs::dir_ls(analysis_dir, glob="*.genes.results", recurse=T)

  # get the file list and pipe it into our extractor function
  rsem_df <- rsem_files %>%
    set_names(.) %>%
    map_dfr(rsem_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from=`file.ID`, values_from=`data$TPM`) %>%
    rename_with(~ basename(.x)) %>%
    rename_with(~ gsub(".genes.results", "", .x, fixed=TRUE)) %>%
    rename(gene_id=`data$gene_id`)

  # filter for Acerv
  rsem_df_org1 <- rsem_df %>%
    filter(grepl(organism1, gene_id))

  # filter for Smic
  rsem_df_org2 <- rsem_df %>%
    filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Bowtie_RSEM_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  write_csv(rsem_df, file=merged_file)
  write_csv(rsem_df_org1, file=org1_file)
  write_csv(rsem_df_org2, file=org2_file)
}

## Function to extract TPM from each file
kallisto_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$target_id, data$tpm))
}

## Extract Kallisto results
extract_kallisto_quants <- function(analysis_dir, outdir) {
  message("Extract Kallisto calculated quants")
  kallisto_files <- fs::dir_ls(analysis_dir, glob="*/abundance.tsv", recurse=T)

  # get the file list and pipe it into our extractor function
  kallisto_df <- kallisto_files %>%
    set_names(.) %>%
    map_dfr(kallisto_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$tpm`) %>%
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    rename_with(~ gsub(KALLISTO_ABUNDANCE_PATH, "", .x, fixed=TRUE)) %>%
    rename(gene_id=`data$target_id`)

  kallisto_df_org1 <- kallisto_df %>% filter(grepl(organism1, gene_id))
  kallisto_df_org2 <- kallisto_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/Kallisto_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')
  write_csv(kallisto_df, file=merged_file)
  write_csv(kallisto_df_org1, file=org1_file)
  write_csv(kallisto_df_org2, file=org2_file)
}

## Function to extract TPM from each file
bwa_quant_extractor = function(fname) {
  data = read_delim(file=fname, delim="\t", show_col_types=FALSE)

  # create and return a data frame like this:
  return(tibble(data$Name, data$TPM))
}

## Extract bwa/Salmon results
extract_bwasalmon_results <- function(analysis_dir, outdir) {
  message("Extract bwa/Salmon results")
  bwa_files <- fs::dir_ls(analysis_dir, regexp=BWA_SALMON_REGEXP, recurse=T)

  # get the file list and pipe it into our extractor function
  bwa_df <- bwa_files %>%
    set_names(.) %>%
    map_dfr(bwa_quant_extractor, .id="file.ID") %>%
    pivot_wider(names_from = `file.ID`, values_from=`data$TPM`) %>%
    rename_with(~ gsub(analysis_dir, "", .x, fixed=TRUE)) %>%
    rename_with(~ gsub(BWA_SALMON_QUANT_SF_PATH, "", .x, fixed=TRUE)) %>%
    rename(gene_id = `data$Name`)

  bwa_df_org1 <- bwa_df %>% filter(grepl(organism1, gene_id))
  bwa_df_org2 <- bwa_df %>% filter(grepl(organism2, gene_id))

  # write into  file
  merged_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_Merged.csv", sep='')
  org1_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism1, ".csv", sep='')
  org2_file = paste(outdir, "/bwa_Salmon_", organisms, "_TPM_matrix_", organism2, ".csv", sep='')

  write_csv(bwa_df, file=BWA_DF_MERGED_FILE)
  write_csv(bwa_df_org1, file=BWA_DF_ORG1_FILE)
  write_csv(bwa_df_org2, file=BWA_DF_ORG2_FILE)
}

RNA_SEQ_ANALYSIS_DIR = "/proj/omics4tb2/wwu/Global_Search/pilot_fail_concat-output/"
OUTDIR = '/proj/omics4tb2/wwu/Global_Search/RNASeq_TPMs'
if (!dir.exists(OUTDIR)) {
  dir.create(OUTDIR)
}

extract_salmon_quants(RNA_SEQ_ANALYSIS_DIR, OUTDIR, STAR_SALMON_REGEXP)

#select <- order(rowMeans(kallisto_df_Acerv)[2:126], decreasing=TRUE)[1:20]
#df <- as.data.frame(colData(dds)[,c("condition","type")])
#pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=FALSE, annotation_col=df)
