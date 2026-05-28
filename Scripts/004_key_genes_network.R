# ============================================================
# Generador de Redes por Categoría Genómica (Región)
# ============================================================
library(tidyverse)
library(igraph)

setwd("D:/DESCARGAS")

df <- read.csv("Supplementary tables.xlsx - Supp. Table2.csv",
               stringsAsFactors = FALSE, na.strings = c("", "NA", NA, "#N/D"))

# ------------------------------------------------------------
# FUNCIÓN: generar_grafos_por_categoria
# ------------------------------------------------------------
generar_grafos_por_categoria <- function(data, region_filter, prefix_output, subset) {
  
  cat(paste0("\n=== PROCESANDO CATEGORÍA: ", region_filter, " | SUBSET: ", subset, " ===\n"))
  
  # Filtrar el dataframe original por la región genómica (ej. "Genebody", "Upstream")
  df_sub <- data %>% filter(Location == region_filter)
  df_sub <- df_sub[which(df_sub$Rice.genotype == subset), ]
  
  # Concatenar Contexto y Dirección de Metilación para formar una feature combinada
  df_sub$Methylation_direction <- paste0(df_sub$Context, "-", df_sub$Methylation_direction)
  
  if (nrow(df_sub) == 0) {
    stop(paste("No se encontraron registros en el df para la categoría:", region_filter, "y subset:", subset))
  }
  
  # ---- 1. Crear nodes_df base para los genes del subset ----
  nodes_df <- df_sub %>%
    group_by(Gene) %>%
    summarise(
      annotation = paste(unique(Gene.annotation), collapse = "//"),
      RXN = paste(unique(Reactions), collapse = "//"),
      .groups = "drop"
    ) %>%
    rename(gene = Gene)
  
  if (nrow(nodes_df) < 2) {
    cat("Muy pocos genes en esta categoría para generar combinaciones.\n")
    return(NULL)
  }
  
  # ---- 2. Crear combinaciones de aristas (Proyección No Dirigida) ----
  gene_pairs <- t(utils::combn(nodes_df$gene, 2))
  
  edges_df <- data.frame(
    FROM = gene_pairs[, 1],
    TO = gene_pairs[, 2],
    stringsAsFactors = FALSE
  )
  
  # Atributos únicos por gen dentro de este subset
  gene_attributes <- df_sub %>%
    group_by(Gene) %>%
    summarise(
      Methylation_direction = if(n_distinct(Methylation_direction, na.rm = TRUE) == 1) first(na.omit(Methylation_direction)) else NA_character_,
      Context = if(n_distinct(Context, na.rm = TRUE) == 1) first(na.omit(Context)) else NA_character_,
      # Azucena_expression = if(n_distinct(Azucena_expression, na.rm = TRUE) == 1) first(na.omit(Azucena_expression)) else NA_character_,
      # BGI_expression = if(n_distinct(BGI_expression, na.rm = TRUE) == 1) first(na.omit(BGI_expression)) else NA_character_,
      Expression = if(n_distinct(Expression, na.rm = TRUE) == 1) first(na.omit(Expression)) else NA_character_,
      
      .groups = "drop"
    )
  
  # Cruzar perfiles y aplicar regla de concordancia estricta
  edges_df <- edges_df %>%
    left_join(gene_attributes, by = c("FROM" = "Gene")) %>%
    left_join(gene_attributes, by = c("TO" = "Gene"), suffix = c("_FROM", "_TO")) %>%
    mutate(
      Methylation_direction = ifelse(Methylation_direction_FROM == Methylation_direction_TO, Methylation_direction_FROM, NA_character_),
      Context               = ifelse(Context_FROM == Context_TO, Context_FROM, NA_character_),
      # Azucena_expression    = ifelse(Azucena_expression_FROM == Azucena_expression_TO, Azucena_expression_FROM, NA_character_),
      # BGI_expression        = ifelse(BGI_expression_FROM == BGI_expression_TO, BGI_expression_FROM, NA_character_)
      Expression        = ifelse(Expression_FROM == Expression_TO, Expression_FROM, NA_character_)
    ) %>%
    select(FROM, TO, Methylation_direction, 
           # Azucena_expression, 
           # BGI_expression, 
           Expression)
  
  # Filtrar filas donde todos los atributos sean NA
  target_cols <- c("Methylation_direction", "Expression")#"Azucena_expression", "BGI_expression")
  edges_df <- edges_df %>% filter(if_any(all_of(target_cols), ~ !is.na(.)))
  
  if (nrow(edges_df) == 0) {
    cat("No existen aristas válidas compartidas para esta categoría.\n")
    return(NULL)
  }
  
  # Guardamos una copia temporal con NAs reales para construir la red dirigida
  edges_df_to_long <- edges_df
  
  # Reemplazar NAs por texto vacío para la visualización del GRAFO 1 en Cytoscape
  edges_df <- edges_df %>% mutate(across(all_of(target_cols), ~ replace_na(.x, "")))
  
  # Filtrar nodos huérfanos del subset de genes
  active_genes <- unique(c(edges_df$FROM, edges_df$TO))
  nodes_df <- nodes_df %>% 
    filter(gene %in% active_genes) %>%
    mutate(across(c(annotation, RXN), ~ replace_na(.x, "")))
  
  # --- GRAFO 1: Red No Dirigida de Genes ---
  g <- graph_from_data_frame(edges_df, vertices = nodes_df, directed = FALSE)
  file_g1 <- paste0(prefix_output, "_", subset, ".graphml")
  write_graph(g, file = file_g1, format = "graphml")
  cat(paste0("-> Guardado Grafo 1 (No Dirigido): ", file_g1, "\n"))
  
  # ---- 3. Transformación a Estructura Bipartita Dirigida (EXCLUYENDO LOCATION) ----
  bipartite_features <- c("Methylation_direction","Expression")
#                          "Azucena_expression", "BGI_expression")
  
  x_edge_list <- edges_df_to_long %>%
    pivot_longer(
      cols = all_of(bipartite_features), 
      names_to = "Variable",
      values_to = "Attribute_Value"
    ) %>%
    filter(Attribute_Value != "" & !is.na(Attribute_Value)) %>%
    mutate(Combined_Feature = paste0(Variable, "_", Attribute_Value)) %>%
    reframe(
      FROM = c(FROM, TO),
      TO = c(Combined_Feature, Combined_Feature)
    ) %>%
    distinct()
  
  # Definir tipos de nodos para el Grafo Bipartito
  nodes_df$type <- "gene"
  
  nuevos_nodos <- data.frame(
    gene = unique(x_edge_list$TO),
    annotation = "",
    RXN = "",
    type = "feature",
    stringsAsFactors = FALSE
  )
  
  nodes_df_bipartite <- rbind(nodes_df, nuevos_nodos) %>%
    distinct(gene, .keep_all = TRUE)
  
  active_nodes_bip <- unique(c(x_edge_list$FROM, x_edge_list$TO))
  nodes_df_bipartite <- nodes_df_bipartite %>% filter(gene %in% active_nodes_bip)
  
  # --- GRAFO 2: Red Dirigida Bipartita ---
  g_v <- graph_from_data_frame(x_edge_list, vertices = nodes_df_bipartite, directed = TRUE)
  
  # 1. Calcular el Out-Degree (Grado de salida) en la red bipartita dirigida
  out_deg <- degree(g_v, mode = "out")
  
  # Filtrar de forma estricta para concentrar sólo genes asociados a un mínimo de 2 features
  # Los nodos tipo 'feature' tienen out-degree = 0, por lo que pasan directo el filtro
  nodos_validos <- names(out_deg[out_deg != 1])
  g_v <- subgraph(g_v, v = nodos_validos)
  
  file_g2 <- paste0(prefix_output, "_", subset, "_directed.graphml")
  write_graph(g_v, file = file_g2, format = "graphml")
  cat(paste0("-> Guardado Grafo 2 (Bipartito Dirigido corregido sin sesgo de Location): ", file_g2, "\n"))
  
  return(list(network_undirected = g, network_bipartite = g_v))
}

# ============================================================
# EJECUCIÓN POR CATEGORÍAS Y SUBSETS (VARIEDADES)
# ============================================================

# --- CATEGORÍA: Genebody ---
región_genebody_AZU <- generar_grafos_por_categoria(
  data = df, 
  region_filter = "Genebody", 
  prefix_output = "Loc_Genebody",
  subset = "Azucena"
)

región_genebody_BGI <- generar_grafos_por_categoria(
  data = df, 
  region_filter = "Genebody", 
  prefix_output = "Loc_Genebody",
  subset = "BGI"
)

# --- CATEGORÍA: Upstream ---
región_upstream_AZU <- generar_grafos_por_categoria(
  data = df, 
  region_filter = "Upstream", 
  prefix_output = "Loc_Upstream",
  subset = "Azucena"
)

región_upstream_BGI <- generar_grafos_por_categoria(
  data = df, 
  region_filter = "Upstream", 
  prefix_output = "Loc_Upstream",
  subset = "BGI"
)

################################################################################
set.seed(1)
################################################################################
# 1. Definir la ruta y dimensiones del archivo (Aumentamos el lienzo para dar más aire)
png(filename = "D:/DESCARGAS/Comparacion_Redes_Bipartitas_Genebody_Labels.png", 
    width = 3600,      # Más ancho para separar los dos paneles
    height = 1800,     # Proporción ideal para evitar solapamientos
    res = 300)         

# 2. Configurar la pantalla
par(mfrow = c(1, 2), mar = c(2, 2, 4, 2)) # Un poco más de margen general

# --- PANEL 1: Azucena ---
CAT_AZU <- región_genebody_AZU$network_bipartite
color_azu <- ifelse(V(CAT_AZU)$type == "gene", "palegreen3", "tomato")

# OPTIMIZACIÓN DEL LAYOUT: Forzamos más iteraciones para que los clústeres se repelan más
layout_azu <- layout_with_fr(CAT_AZU, niter = 1000)

plot(CAT_AZU, 
     layout = layout_azu,
     vertex.color = color_azu,
     vertex.label = V(CAT_AZU)$name,
     vertex.label.color = "black",            
     vertex.size = 7,                          # Nodos ligeramente más pequeños para liberar espacio
     vertex.label.font = ifelse(V(CAT_AZU)$type == "feature", 2, 1), 
     vertex.label.cex = 0.35,                  # Fuente un poco más pequeña pero nítida por los 300 DPI
     vertex.label.dist = 0.5,                  # Distancia corta para que el texto orbite pegado al nodo
     edge.color = "gray88",                    # Aristas más tenues para que no compitan con el texto
     edge.arrow.size = 0.18,                   # Flechas más discretas
     edge.arrow.width = 0.4,
     main = "Genebody - Azucena"
)

# --- PANEL 2: BGI ---
CAT_BGI <- región_genebody_BGI$network_bipartite
color_bgi <- ifelse(V(CAT_BGI)$type == "gene", "palegreen3", "tomato")

# Mismo ajuste de repulsión para el panel BGI
layout_bgi <- layout_with_fr(CAT_BGI, niter = 1000)

plot(CAT_BGI, 
     layout = layout_bgi,
     vertex.color = color_bgi,
     vertex.label = V(CAT_BGI)$name,
     vertex.label.color = "black",
     vertex.size = 7,
     vertex.label.font = ifelse(V(CAT_BGI)$type == "feature", 2, 1), 
     vertex.label.cex = 0.35,                 
     vertex.label.dist = 0.5,                           
     edge.color = "gray88",
     edge.arrow.size = 0.18,
     edge.arrow.width = 0.4,
     main = "Genebody - BGI"
)

# Añadir la leyenda en una zona despejada
legend("bottomleft", 
       legend = c("Genes (>= 2 features)", "Features"), 
       col = c("palegreen3", "tomato"), 
       pch = 21, pt.bg = c("palegreen3", "tomato"), 
       cex = 0.8, bty = "n")

# 3. Cerrar dispositivo y guardar imagen
dev.off()

cat("\n[INFO] Gráfico corregido con etiquetas optimizadas guardado en Descargas.\n")


################################################################################
################################################################################
# 1. Definir la ruta y dimensiones del archivo (Aumentamos el lienzo para dar más aire)
png(filename = "D:/DESCARGAS/Comparacion_Redes_Bipartitas_Upstream_Labels.png", 
    width = 3600,      # Más ancho para separar los dos paneles
    height = 1800,     # Proporción ideal para evitar solapamientos
    res = 300)         

# 2. Configurar la pantalla
par(mfrow = c(1, 2), mar = c(2, 2, 4, 2)) # Un poco más de margen general

# --- PANEL 1: Azucena ---
CAT_AZU <- región_upstream_AZU$network_bipartite
color_azu <- ifelse(V(CAT_AZU)$type == "gene", "palegreen3", "tomato")

# OPTIMIZACIÓN DEL LAYOUT: Forzamos más iteraciones para que los clústeres se repelan más
layout_azu <- layout_with_fr(CAT_AZU, niter = 1000)

plot(CAT_AZU, 
     layout = layout_azu,
     vertex.color = color_azu,
     vertex.label = V(CAT_AZU)$name,
     vertex.label.color = "black",            
     vertex.size = 7,                          # Nodos ligeramente más pequeños para liberar espacio
     vertex.label.font = ifelse(V(CAT_AZU)$type == "feature", 2, 1), 
     vertex.label.cex = 0.35,                  # Fuente un poco más pequeña pero nítida por los 300 DPI
     vertex.label.dist = 0.5,                  # Distancia corta para que el texto orbite pegado al nodo
     edge.color = "gray88",                    # Aristas más tenues para que no compitan con el texto
     edge.arrow.size = 0.18,                   # Flechas más discretas
     edge.arrow.width = 0.4,
     main = "Upstream - Azucena"
)

# --- PANEL 2: BGI ---
CAT_BGI <- región_upstream_BGI$network_bipartite
color_bgi <- ifelse(V(CAT_BGI)$type == "gene", "palegreen3", "tomato")

# Mismo ajuste de repulsión para el panel BGI
layout_bgi <- layout_with_fr(CAT_BGI, niter = 1000)

plot(CAT_BGI, 
     layout = layout_bgi,
     vertex.color = color_bgi,
     vertex.label = V(CAT_BGI)$name,
     vertex.label.color = "black",
     vertex.size = 7,
     vertex.label.font = ifelse(V(CAT_BGI)$type == "feature", 2, 1), 
     vertex.label.cex = 0.35,                 
     vertex.label.dist = 0.5,                           
     edge.color = "gray88",
     edge.arrow.size = 0.18,
     edge.arrow.width = 0.4,
     main = "Upstream - BGI"
)

# Añadir la leyenda en una zona despejada
legend("bottomleft", 
       legend = c("Genes (>= 2 features)", "Features"), 
       col = c("palegreen3", "tomato"), 
       pch = 21, pt.bg = c("palegreen3", "tomato"), 
       cex = 1, bty = "n")

# 3. Cerrar dispositivo y guardar imagen
dev.off()

cat("\n[INFO] Gráfico corregido con etiquetas optimizadas guardado en Descargas.\n")
  