require(gprofiler2)
require(data.table)
x <- read.csv("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_GENES.csv")
genes <- unique(x$ID_gene)

unique_context <- unique(x$file)

CG_hyper <- intersect(x$ID_gene[which(x$file=="Azucena (CG context)" & x$Status=="hyper")],
          x$ID_gene[which(x$file=="BGI (CG context)" & x$Status=="hyper")])


CG_hypo <- intersect(x$ID_gene[which(x$file=="Azucena (CG context)" & x$Status=="hypo")],
                    x$ID_gene[which(x$file=="BGI (CG context)" & x$Status=="hypo")])

CHH_hyper <- intersect(x$ID_gene[which(x$file=="Azucena (CHH context)" & x$Status=="hyper")],
                      x$ID_gene[which(x$file=="BGI (CHH context)" & x$Status=="hyper")])

CHH_hypo <- intersect(x$ID_gene[which(x$file=="Azucena (CHH context)" & x$Status=="hypo")],
                     x$ID_gene[which(x$file=="BGI (CHH context)" & x$Status=="hypo")])

CHG_hyper <- intersect(x$ID_gene[which(x$file=="Azucena (CHG context)" & x$Status=="hyper")],
                       x$ID_gene[which(x$file=="BGI (CHG context)" & x$Status=="hyper")])

CHG_hypo <- intersect(x$ID_gene[which(x$file=="Azucena (CHG context)" & x$Status=="hypo")],
                      x$ID_gene[which(x$file=="BGI (CHG context)" & x$Status=="hypo")])


shared_genes <- unique(c(CG_hyper,CG_hypo,
  CHH_hyper,CHH_hypo,
  CHG_hyper,CHG_hypo)
)

################################################################################

BGI_GEN <- unique(x$ID_gene[x$file %in% c("BGI (CG context)","BGI (CHH context)","BGI (CHG context)")])
AZU_GEN <- unique(x$ID_gene[x$file %in% c("Azucena (CG context)","Azucena (CHH context)","Azucena (CHG context)")])

################################################################################
x <- read.csv("D:/JJG_OCT25/RESULTS_MET/file_DMR_overlap_Upstream.csv")
genes_upstream <- unique(x$ID_gene)
unique_context <- unique(x$file)

CG_hyper_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CG context)" & x$Status=="hyper")],
                      x$ID_gene[which(x$file=="BGI (CG context)" & x$Status=="hyper")])


CG_hypo_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CG context)" & x$Status=="hypo")],
                     x$ID_gene[which(x$file=="BGI (CG context)" & x$Status=="hypo")])

CHH_hyper_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CHH context)" & x$Status=="hyper")],
                       x$ID_gene[which(x$file=="BGI (CHH context)" & x$Status=="hyper")])

CHH_hypo_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CHH context)" & x$Status=="hypo")],
                      x$ID_gene[which(x$file=="BGI (CHH context)" & x$Status=="hypo")])

CHG_hyper_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CHG context)" & x$Status=="hyper")],
                       x$ID_gene[which(x$file=="BGI (CHG context)" & x$Status=="hyper")])

CHG_hypo_UP <- intersect(x$ID_gene[which(x$file=="Azucena (CHG context)" & x$Status=="hypo")],
                      x$ID_gene[which(x$file=="BGI (CHG context)" & x$Status=="hypo")])


shared_upstream <- unique(c(CG_hyper_UP,CG_hypo_UP,
                         CHH_hyper_UP,CHH_hypo_UP,
                         CHG_hyper_UP,CHG_hypo_UP)
)

################################################################################
genes_list <- list(
  shared_genes = shared_genes,
  shared_upstream = shared_upstream
)


################################################################################

BGI_UP <- unique(x$ID_gene[x$file %in% c("BGI (CG context)","BGI (CHH context)","BGI (CHG context)")])
AZU_UP <- unique(x$ID_gene[x$file %in% c("Azucena (CG context)","Azucena (CHH context)","Azucena (CHG context)")])

################################################################################
#genes
# genes_to <- unique(c(genes,genes_upstream))
x_s2 <-  gprofiler2::gost(query = genes_list,
                          organism = "osativa", ordered_query = F,
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                          measure_underrepresentation = FALSE, evcodes = FALSE,
                          user_threshold = 0.05, correction_method = "g_SCS",
                          domain_scope = "annotated", custom_bg = NULL,
                          numeric_ns = "", sources = c("GO:BP","GO:MF","KEGG"), as_short_link = FALSE)

results_GO <- x_s2$result
fwrite(as.data.table(results_GO), "D:/JJG_OCT25/RESULTS_MET/GOST_SHARED.csv")
################################################################################
#entre variedades

var_list <- list(
  BGI_UP = BGI_UP,
  AZU_UP = AZU_UP,
  BGI_GEN = BGI_GEN,
  AZU_GEN = AZU_GEN
)


var_s2 <-  gprofiler2::gost(query = var_list,
                         organism = "osativa", ordered_query = F,
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                         measure_underrepresentation = FALSE, evcodes = FALSE,
                         user_threshold = 0.05, correction_method = "g_SCS",
                         domain_scope = "annotated", custom_bg = NULL,
                         numeric_ns = "", sources = c("GO:BP","GO:MF","KEGG"), as_short_link = FALSE)

results_GO_var <- var_s2$result
fwrite(as.data.table(results_GO_var), "D:/JJG_OCT25/RESULTS_MET/GOST_VARIETY.csv")
################################################################################
