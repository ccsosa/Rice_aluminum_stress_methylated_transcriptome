library(data.table)
library(ggpubr)

# Optimized coverage calculation (vectorized)
calculate_coverage_vectorized <- function(dmr_start, dmr_end, gene_start, gene_end) {
  dmr_length <- dmr_end - dmr_start + 1
  overlap_start <- pmax(dmr_start, gene_start)
  overlap_end <- pmin(dmr_end, gene_end)
  overlap_length <- pmax(0, overlap_end - overlap_start + 1)
  coverage_percentage <- (overlap_length / dmr_length) * 100
  return(coverage_percentage)
}

# Read genes file
genes <- fread("D:/JJG_OCT25/Data/Genomic_features_sativa/genes.csv")
setDT(genes)

# Define input parameters
inDir <- "D:/JJG_OCT25/Data/Epigenome/DMRs_annotated"
files <- c("AZU_CG_AC_AT", "AZU_CHG_AC_AT", "AZU_CHH_AC_AT",
           "BGI_CG_BC_BT", "BGI_CHG_BC_BT", "BGI_CHH_BC_BT")

chromosomes <- paste0("chr", sprintf("%02d", 1:12))

# Main processing loop
files_genes <- list()

for(h in seq_along(files)) {
  # h <- 1
  cat("\nProcessing:", files[h], "\n")
  # Read DMR file
  file_DMR <- fread(paste0(inDir, "/", files[h], "_StateCalls-filtered.txt_Annotated"))
  file_DMR <- file_DMR[, c(1:4,11:12)]
  setDT(file_DMR)
  
  # Process all chromosomes at once using data.table join
  all_overlaps <- list()
  
  for(i in seq_along(chromosomes)) {
    cat("  Chromosome:", chromosomes[i], "\n")
    
    # Subset data for current chromosome
    genes_chr <- genes[Chromosome == chromosomes[i]]
    dmr_chr <- file_DMR[seqnames == chromosomes[i]]
    
    if(nrow(genes_chr) == 0 || nrow(dmr_chr) == 0) next
    
    # Create keys for overlap join
    setkey(genes_chr, Start, End)
    setkey(dmr_chr, start, end)
    
    # Perform overlap join using foverlaps
    overlaps <- foverlaps(
      dmr_chr,
      genes_chr,
      by.x = c("start", "end"),
      by.y = c("Start", "End"),
      type = "any",
      nomatch = NULL
    )
    
    if(nrow(overlaps) > 0) {
      # Calculate coverage percentage
      overlaps[, coverage_percentage := calculate_coverage_vectorized(
        start, end, Start, End
      )]
      
      # Filter by coverage threshold
      overlaps <- overlaps[coverage_percentage >= 60]
      
      if(nrow(overlaps) > 0) {
        # Select relevant columns and rename
        overlaps[, Gene := ID_gene]
        overlaps[, overlap := TRUE]
        
        all_overlaps[[i]] <- overlaps
      }
    }
  }
  
  # Combine all chromosomes
  if(length(all_overlaps) > 0) {
    chr_i_list <- rbindlist(all_overlaps, fill = TRUE)
    chr_i_list[, file := files[h]]
    files_genes[[h]] <- chr_i_list
  }
  
  cat("  Found", if(length(all_overlaps) > 0) nrow(chr_i_list) else 0, "overlapping DMRs\n")
}

# Combine all files
final_result <- rbindlist(files_genes, fill = TRUE)

# Save output
# fwrite(final_result, "D:/JJG_OCT25/file_DMR_overlap.csv")

cat("\nProcessing complete! Total overlaps found:", nrow(final_result), "\n")

# Summary statistics
cat("\nSummary by file:\n")
print(final_result[, .N, by = file])

cat("\nSummary by chromosome:\n")
print(final_result[, .N, by = seqnames][order(seqnames)])

replace_fi <- c("Azucena (CG context)", "Azucena (CHG context)", "Azucena (CHH context)",
                "BGI (CG context)", "BGI (CHG context)", "BGI (CHH context)")

for(i in 1:length(files)){
  final_result$file[(final_result$file==files[[i]])] <- replace_fi[[i]]
}

# Save output
fwrite(final_result, "D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_GENES.csv")
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Read upstream file
genes <- fread("D:/JJG_OCT25/Data/Genomic_features_sativa/upstream.csv")
setDT(genes)

# Define input parameters
inDir <- "D:/JJG_OCT25/Data/Epigenome/DMRs_annotated"
files <- c("AZU_CG_AC_AT", "AZU_CHG_AC_AT", "AZU_CHH_AC_AT",
           "BGI_CG_BC_BT", "BGI_CHG_BC_BT", "BGI_CHH_BC_BT")

chromosomes <- paste0("chr", sprintf("%02d", 1:12))

# Main processing loop
files_genes <- list()

for(h in seq_along(files)) {
  cat("\nProcessing:", files[h], "\n")
  # Read DMR file
  file_DMR <- fread(paste0(inDir, "/", files[h], "_StateCalls-filtered.txt_Annotated"))
  file_DMR <- file_DMR[, c(1:4,11:12)]
  setDT(file_DMR)
  
  # Process all chromosomes at once using data.table join
  all_overlaps <- list()
  
  for(i in seq_along(chromosomes)) {
    cat("  Chromosome:", chromosomes[i], "\n")
    
    # Subset data for current chromosome
    genes_chr <- genes[Chromosome == chromosomes[i]]
    dmr_chr <- file_DMR[seqnames == chromosomes[i]]
    
    if(nrow(genes_chr) == 0 || nrow(dmr_chr) == 0) next
    
    # Create keys for overlap join
    setkey(genes_chr, Start, End)
    setkey(dmr_chr, start, end)
    
    # Perform overlap join using foverlaps
    overlaps <- foverlaps(
      dmr_chr,
      genes_chr,
      by.x = c("start", "end"),
      by.y = c("Start", "End"),
      type = "any",
      nomatch = NULL
    )
    
    if(nrow(overlaps) > 0) {
      # Calculate coverage percentage
      overlaps[, coverage_percentage := calculate_coverage_vectorized(
        start, end, Start, End
      )]
      
      # Filter by coverage threshold
      overlaps <- overlaps[coverage_percentage >= 60]
      
      if(nrow(overlaps) > 0) {
        # Select relevant columns and rename
        overlaps[, Gene := ID_gene]
        overlaps[, overlap := TRUE]
        
        all_overlaps[[i]] <- overlaps
      }
    }
  }
  
  # Combine all chromosomes
  if(length(all_overlaps) > 0) {
    chr_i_list <- rbindlist(all_overlaps, fill = TRUE)
    chr_i_list[, file := files[h]]
    files_genes[[h]] <- chr_i_list
  }
  
  cat("  Found", if(length(all_overlaps) > 0) nrow(chr_i_list) else 0, "overlapping DMRs\n")
}

# Combine all files
final_result <- rbindlist(files_genes, fill = TRUE)

# Save output
# fwrite(final_result, "D:/JJG_OCT25/file_DMR_overlap.csv")

cat("\nProcessing complete! Total overlaps found:", nrow(final_result), "\n")

# Summary statistics
cat("\nSummary by file:\n")
print(final_result[, .N, by = file])

cat("\nSummary by chromosome:\n")
print(final_result[, .N, by = seqnames][order(seqnames)])

replace_fi <- c("Azucena (CG context)", "Azucena (CHG context)", "Azucena (CHH context)",
                "BGI (CG context)", "BGI (CHG context)", "BGI (CHH context)")

for(i in 1:length(files)){
  final_result$file[(final_result$file==files[[i]])] <- replace_fi[[i]]
}


fwrite(final_result, "D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_Upstream.csv")

#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# # Read genes file
# genes <- fread("D:/JJG_OCT25/Data/Genomic_features_sativa/downstream.csv")
# setDT(genes)
# 
# # Define input parameters
# inDir <- "D:/JJG_OCT25/Data/Epigenome/DMRs_annotated"
# files <- c("AZU_CG_AC_AT", "AZU_CHG_AC_AT", "AZU_CHH_AC_AT",
#            "BGI_CG_BC_BT", "BGI_CHG_BC_BT", "BGI_CHH_BC_BT")
# 
# chromosomes <- paste0("chr", sprintf("%02d", 1:12))
# 
# # Main processing loop
# files_genes <- list()
# 
# for(h in seq_along(files)) {
#   cat("\nProcessing:", files[h], "\n")
#   # Read DMR file
#   file_DMR <- fread(paste0(inDir, "/", files[h], "_StateCalls-filtered.txt_Annotated"))
#   file_DMR <- file_DMR[, c(1:4,11:12)]
#   setDT(file_DMR)
#   
#   # Process all chromosomes at once using data.table join
#   all_overlaps <- list()
#   
#   for(i in seq_along(chromosomes)) {
#     cat("  Chromosome:", chromosomes[i], "\n")
#     
#     # Subset data for current chromosome
#     genes_chr <- genes[Chromosome == chromosomes[i]]
#     dmr_chr <- file_DMR[seqnames == chromosomes[i]]
#     
#     if(nrow(genes_chr) == 0 || nrow(dmr_chr) == 0) next
#     
#     # Create keys for overlap join
#     setkey(genes_chr, Start, End)
#     setkey(dmr_chr, start, end)
#     
#     # Perform overlap join using foverlaps
#     overlaps <- foverlaps(
#       dmr_chr,
#       genes_chr,
#       by.x = c("start", "end"),
#       by.y = c("Start", "End"),
#       type = "any",
#       nomatch = NULL
#     )
#     
#     if(nrow(overlaps) > 0) {
#       # Calculate coverage percentage
#       overlaps[, coverage_percentage := calculate_coverage_vectorized(
#         start, end, Start, End
#       )]
#       
#       # Filter by coverage threshold
#       overlaps <- overlaps[coverage_percentage >= 60]
#       
#       if(nrow(overlaps) > 0) {
#         # Select relevant columns and rename
#         overlaps[, Gene := ID_gene]
#         overlaps[, overlap := TRUE]
#         
#         all_overlaps[[i]] <- overlaps
#       }
#     }
#   }
#   
#   # Combine all chromosomes
#   if(length(all_overlaps) > 0) {
#     chr_i_list <- rbindlist(all_overlaps, fill = TRUE)
#     chr_i_list[, file := files[h]]
#     files_genes[[h]] <- chr_i_list
#   }
#   
#   cat("  Found", if(length(all_overlaps) > 0) nrow(chr_i_list) else 0, "overlapping DMRs\n")
# }
# 
# # Combine all files
# final_result <- rbindlist(files_genes, fill = TRUE)
# 
# # Save output
# # fwrite(final_result, "D:/JJG_OCT25/file_DMR_overlap.csv")
# 
# cat("\nProcessing complete! Total overlaps found:", nrow(final_result), "\n")
# 
# # Summary statistics
# cat("\nSummary by file:\n")
# print(final_result[, .N, by = file])
# 
# cat("\nSummary by chromosome:\n")
# print(final_result[, .N, by = seqnames][order(seqnames)])
# 
# replace_fi <- c("Azucena (CG context)", "Azucena (CHG context)", "Azucena (CHH context)",
#                 "BGI (CG context)", "BGI (CHG context)", "BGI (CHH context)")
# 
# for(i in 1:length(files)){
#   final_result$file[(final_result$file==files[[i]])] <- replace_fi[[i]]
# }
# 
# 
# fwrite(final_result, "D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_Downstream.csv")
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
#################################################################################
# Read repeats file
genes <- fread("D:/JJG_OCT25/Data/Genomic_features_sativa/repeats.csv")
setDT(genes)

# Define input parameters
inDir <- "D:/JJG_OCT25/Data/Epigenome/DMRs_annotated"
files <- c("AZU_CG_AC_AT", "AZU_CHG_AC_AT", "AZU_CHH_AC_AT",
           "BGI_CG_BC_BT", "BGI_CHG_BC_BT", "BGI_CHH_BC_BT")

chromosomes <- paste0("chr", sprintf("%02d", 1:12))

# Main processing loop
files_genes <- list()

for(h in seq_along(files)) {
  cat("\nProcessing:", files[h], "\n")
  # Read DMR file
  file_DMR <- fread(paste0(inDir, "/", files[h], "_StateCalls-filtered.txt_Annotated"))
  file_DMR <- file_DMR[, c(1:4,11:12)]
  setDT(file_DMR)
  
  # Process all chromosomes at once using data.table join
  all_overlaps <- list()
  
  for(i in seq_along(chromosomes)) {
    cat("  Chromosome:", chromosomes[i], "\n")
    
    # Subset data for current chromosome
    genes_chr <- genes[Chromosome == chromosomes[i]]
    dmr_chr <- file_DMR[seqnames == chromosomes[i]]
    
    if(nrow(genes_chr) == 0 || nrow(dmr_chr) == 0) next
    
    # Create keys for overlap join
    setkey(genes_chr, Start, End)
    setkey(dmr_chr, start, end)
    
    # Perform overlap join using foverlaps
    overlaps <- foverlaps(
      dmr_chr,
      genes_chr,
      by.x = c("start", "end"),
      by.y = c("Start", "End"),
      type = "any",
      nomatch = NULL
    )
    
    if(nrow(overlaps) > 0) {
      # Calculate coverage percentage
      overlaps[, coverage_percentage := calculate_coverage_vectorized(
        start, end, Start, End
      )]
      
      # Filter by coverage threshold
      overlaps <- overlaps[coverage_percentage >= 60]
      
      if(nrow(overlaps) > 0) {
        # Select relevant columns and rename
        overlaps[, Gene := ID_gene]
        overlaps[, overlap := TRUE]
        
        all_overlaps[[i]] <- overlaps
      }
    }
  }
  
  # Combine all chromosomes
  if(length(all_overlaps) > 0) {
    chr_i_list <- rbindlist(all_overlaps, fill = TRUE)
    chr_i_list[, file := files[h]]
    files_genes[[h]] <- chr_i_list
  }
  
  cat("  Found", if(length(all_overlaps) > 0) nrow(chr_i_list) else 0, "overlapping DMRs\n")
}

# Combine all files
final_result <- rbindlist(files_genes, fill = TRUE)

# Save output
# fwrite(final_result, "D:/JJG_OCT25/file_DMR_overlap.csv")

cat("\nProcessing complete! Total overlaps found:", nrow(final_result), "\n")

# Summary statistics
cat("\nSummary by file:\n")
print(final_result[, .N, by = file])

cat("\nSummary by chromosome:\n")
print(final_result[, .N, by = seqnames][order(seqnames)])

replace_fi <- c("Azucena (CG context)", "Azucena (CHG context)", "Azucena (CHH context)",
                "BGI (CG context)", "BGI (CHG context)", "BGI (CHH context)")

for(i in 1:length(files)){
  final_result$file[(final_result$file==files[[i]])] <- replace_fi[[i]]
}

fwrite(final_result, "D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_repeats.csv")

################################################################################
genes <- fread("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_GENES.csv");genes$dataset <- "genes"
up <- fread("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_Upstream.csv");up$dataset <- "Upstream"
repeats <- fread("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_repeats.csv");repeats$dataset <- "repeats"


total_file <- rbindlist(list(genes, up, repeats), fill = TRUE)
total_file$context <- paste0(total_file$Status," - ",total_file$file," - ",total_file$dataset)
fwrite(total_file, "D:/JJG_OCT25/RESULTS_MET/file_DMR_overlaps.csv")
################################################################################
#load transcriptomic data
Azu_DEG <- fread("D:/JJG_OCT25/Data/Transriptome/Azu_Deseq2_p0.05.csv")
Azu_DEG <- Azu_DEG[which(abs(Azu_DEG$log2FoldChange)>=2),]
BGI_DEG <- fread("D:/JJG_OCT25/Data/Transriptome/BGI_Deseq2_p0.05.csv")
BGI_DEG <- BGI_DEG[which(abs(BGI_DEG$log2FoldChange)>=2),]
################################################################################
#getting lists of unique
cats <- unique(total_file$context)
gene_cats <- list()
for(i in 1:length(cats)){
  gene_cats[[i]] <- unique(total_file$Gene[which(total_file$context==cats[[i]])])
}
names(gene_cats) <- cats
#adding transcriptome
gene_cats$AZU_Expressed_genes <- unique(Azu_DEG$Gene_Id)
gene_cats$BGI_Expressed_genes <- unique(BGI_DEG$Gene_Id)
saveRDS(gene_cats,"D:/JJG_OCT25/RESULTS_MET/Genes_list.RDS")
################################################################################


library(dplyr)
library(tidyr)

# 1. Convert the list to a long, tidy data frame
gene_df <- stack(gene_cats) |> 
  rename(Gene_ID = values, Category = ind)

# 2. Count occurrences and list the categories per gene
gene_overlaps <- gene_df |> 
  group_by(Gene_ID) |> 
  summarize(
    Total_Lists_Present_In = n(),
    Categories = paste(Category, collapse = " | ")
  ) |> 
  arrange(desc(Total_Lists_Present_In))

# View genes that appear in more than one list
shared_genes_all <- filter(gene_overlaps, Total_Lists_Present_In > 1)
print(shared_genes_all)
write.csv(shared_genes_all,"D:/JJG_OCT25/RESULTS_MET/SHARED_GENES_REPEATS.csv")

################################################################################
# 1. Convert your list of genes into a binary 1/0 overlap matrix
require(UpSetR)
upset_data <- fromList(gene_cats)
###

# Transformamos la lista en un data.table de dos columnas
df_contextos <- rbindlist(
  lapply(names(gene_cats), function(nm) {
    data.table(Contexto = nm, Elemento = gene_cats[[nm]])
  })
)

# Ver el resultado
print(df_contextos)
write.csv(df_contextos,"D:/JJG_OCT25/RESULTS_MET/GENES_REPEATS.csv")
###
# 2. Plot the intersections
upset(
  upset_data, 
  nsets = length(gene_cats),          # Include ALL of your categories
  nintersects = 30,                   # Show the top 30 largest overlapping groups
  order.by = "freq",                  # Sort intersections from largest to smallest
  decreasing = TRUE,                  
  mb.ratio = c(0.6, 0.4),             # Balance the height of the top bars vs. bottom matrix
  text.scale = c(1.3, 1.3, 1, 1, 1.5, 1) # Scales text size to ensure long labels fit nicely
)

################################################################################
library(UpSetR)

# 1. Copiar y limpiar nombres de las categorías para evitar colapsos visuales
gene_cats_clean <- gene_cats
names(gene_cats_clean) <- names(gene_cats_clean) |> 
  gsub(pattern = " - Azucena", replacement = " Azu") |> 
  gsub(pattern = " - BGI", replacement = " BGI") |> 
  
  gsub(pattern = "context", replacement = "") |> 
  gsub(pattern = "context ", replacement = "") |> 
  gsub(pattern = " \\( ", replacement = " (") |> 
  gsub(pattern = " \\) ", replacement = ")")

# 1. Obtener la lista única de TODOS los genes globales
all_genes <- unique(unlist(gene_cats_clean))

# 2. Construir la matriz binaria a mano MANTENIENDO los row.names
matrix_data <- as.data.frame(lapply(gene_cats_clean, function(x) as.integer(all_genes %in% x)))
rownames(matrix_data) <- all_genes  # <--- ¡Aquí rescatamos los nombres de los genes!

# 3. Filtrar para dejar SOLO los intersectos (filas que suman más de 1)
matrix_intersect_only <- matrix_data[rowSums(matrix_data) > 1, ]
colnames(matrix_intersect_only) <- names(gene_cats_clean)

png(
  filename = "D:/JJG_OCT25/RESULTS_MET/upset_plot_intersectos_1000dpi.png", 
  width = 36,   # aumentar un poco
  height = 12,  # más altura para 38 sets sin solaparse
  units = "in", 
  res = 1000
)

upset(
  matrix_intersect_only, 
  nsets = length(gene_cats_clean),
  nintersects = 15,
  order.by = "freq",
  mb.ratio = c(0.45, 0.55),
  point.size = 4.3,
  line.size = 1.5,
  text.scale = c(1.3, 0.9, 0.9, 1.0, 1.4, 2.0),
  set_size.show = FALSE,   # <--- elimina el barplot izquierdo, libera espacio para labels
  set_size.angles = 0
)
dev.off()
################################################################################
#now enrichment for only genes!
genes_lists <- gene_cats[c(1:12,37,38)]
names(genes_lists)

x_s2 <-  gprofiler2::gost(query = genes_lists,
                          organism = "osativa", ordered_query = F,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "g_SCS",
                          domain_scope = "annotated", custom_bg = unique(unlist(genes_lists)),
                          numeric_ns = "", sources = c("GO:BP","GO:MF","KEGG"), as_short_link = FALSE)

results_GO <- x_s2$result
fwrite(as.data.table(results_GO), "D:/JJG_OCT25/RESULTS_MET/GOST.csv")
