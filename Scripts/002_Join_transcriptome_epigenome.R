#load transcriptomic data
Azu_DEG <- fread("D:/JJG_OCT25/Data/Transriptome/Azu_Deseq2_p0.05.csv")
Azu_DEG <- Azu_DEG[which(abs(Azu_DEG$log2FoldChange)>=2),]
Azu_DEG$STATUS <- Azu_DEG$log2FoldChange>=2
# Azu_DEG$STATUS[which(Azu_DEG$STATUS==T)] <- "UP"
# Azu_DEG$STATUS[which(Azu_DEG$STATUS==F)] <- "DOWN"
BGI_DEG <- fread("D:/JJG_OCT25/Data/Transriptome/BGI_Deseq2_p0.05.csv")
BGI_DEG$STATUS <- BGI_DEG$log2FoldChange>=2
# BGI_DEG$STATUS_LFC <- NA
# BGI_DEG$STATUS[which(BGI_DEG$STATUS==T)] <- "UP"
# BGI_DEG$STATUS[which(BGI_DEG$STATUS==F)] <- "DOWN"

Azu_DEG$STATUS <- ifelse(Azu_DEG$log2FoldChange >= 2, "UP", "DOWN")
BGI_DEG$STATUS <- ifelse(BGI_DEG$log2FoldChange >= 2, "UP", "DOWN")
DMR <- fread("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlaps.csv")


# cats <- unique(DMR$context)
cats <- c("hyper - Azucena (CG context) - genes",     "hypo - Azucena (CG context) - genes"  ,   
           "hypo - Azucena (CHG context) - genes",     "hyper - Azucena (CHG context) - genes"  , 
           "hypo - Azucena (CHH context) - genes",     "hyper - Azucena (CHH context) - genes"  , 
           "hyper - BGI (CG context) - genes",         "hypo - BGI (CG context) - genes"       ,  
           "hypo - BGI (CHG context) - genes",         "hyper - BGI (CHG context) - genes"      , 
           "hypo - BGI (CHH context) - genes",         "hyper - BGI (CHH context) - genes"     , 
           "hyper - Azucena (CG context) - Upstream",  "hypo - Azucena (CG context) - Upstream"  ,
           "hypo - Azucena (CHG context) - Upstream",  "hyper - Azucena (CHG context) - Upstream",
           "hypo - Azucena (CHH context) - Upstream",  "hyper - Azucena (CHH context) - Upstream",
           "hypo - BGI (CG context) - Upstream",       "hyper - BGI (CG context) - Upstream"     ,
           "hypo - BGI (CHG context) - Upstream",      "hyper - BGI (CHG context) - Upstream"   , 
           "hypo - BGI (CHH context) - Upstream",      "hyper - BGI (CHH context) - Upstream"   )
gene_cats <- list()
for(i in 1:length(cats)){
  gene_cats[[i]] <- unique(DMR$Gene[which(DMR$context==cats[[i]])])
}

names(gene_cats) <- cats
gene_cats$AZU_UP <- unique(Azu_DEG$Gene_Id[which(Azu_DEG$STATUS=="UP")])
gene_cats$AZU_DOWN <- unique(Azu_DEG$Gene_Id[which(Azu_DEG$STATUS=="DOWN")])
gene_cats$BGI_UP <- unique(BGI_DEG$Gene_Id[which(Azu_DEG$STATUS=="UP")])
gene_cats$BGI_DOWN <- unique(BGI_DEG$Gene_Id[which(Azu_DEG$STATUS=="DOWN")])

x_s2 <-  gprofiler2::gost(query = gene_cats,
                          organism = "osativa", ordered_query = F,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "g_SCS",
                          domain_scope = "annotated", custom_bg = NULL,
                          numeric_ns = "", sources = c("GO:BP","GO:MF","KEGG"), as_short_link = FALSE)

results_GO <- x_s2$result
fwrite(as.data.table(results_GO), "D:/JJG_OCT25/RESULTS_MET/GOST_SHARED_TRANS_EPI.csv")
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
  filename = "D:/JJG_OCT25/RESULTS_MET/upset_plot_EPI_TRANS.png", 
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
###################################################
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
write.csv(df_contextos,"D:/JJG_OCT25/RESULTS_MET/EPI_TRANS.csv")
###

unique(df_contextos$Contexto)


# Genes en cualquiera de las primeras 24 categorías (methylation)
methyl_cats <- unique(df_contextos$Contexto)[1:24]
genes_methyl <- unique(df_contextos$Elemento[df_contextos$Contexto %in% methyl_cats])

# Genes en cualquiera de las 4 categorías transcriptómicas
trans_cats <- c("AZU_UP", "AZU_DOWN", "BGI_UP", "BGI_DOWN")
genes_trans <- unique(df_contextos$Elemento[df_contextos$Contexto %in% trans_cats])

# Intersección: genes presentes en ambos grupos
epi_trans_overlap <- intersect(genes_methyl, genes_trans)

length(epi_trans_overlap)



# Genes en overlap
epi_trans_overlap <- intersect(genes_methyl, genes_trans)

# Filtrar df_contextos para esos genes — conserva todas las categorías donde aparecen
df_overlap <- df_contextos[Elemento %in% epi_trans_overlap]

# Ver cuántas categorías por gen
df_overlap <- df_overlap[, .(categorias = paste(Contexto, collapse = " | "),
               n_cats = .N), by = Elemento]

write.csv(df_overlap, "D:/JJG_OCT25/RESULTS_MET/EPI_TRANS_overlap.csv", row.names = FALSE)

 
