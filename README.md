# Rice_aluminum_stress_methylated_transcriptome
This repository consist of R codes to overlap results of Azucena, and BGI rices varieties for epigenomics, and transcriptomics experiments all together to identify main key genes or mechanism implicated in aluminum stress tolerance.

## This approach consist of the following steps:
- Use files in tabulate format with the Differentially Methylated Regions (DMR) for CG, CHG, and CHH contexts for the Azucena and BGI rice varieties and overlap with the genes, and upstream regions files from RAPD annotation database, available at https://rapdb.dna.naro.go.jp/download/irgsp1.html. Only genes with at least 60% coverage are considered
- Addd the genes with absolute LFC >= 2 from a differentially expressed genes (DEG) and put in a nested listed with the purpose of identify shared genes
- Perform gene enrichment analyses (GSEA) for the differentially methylated and upregulated or downregulated gene lists 
- Obtain a bipartite network with genes and methylated and transcripted features (Methylation context,and up or downregulated genes)
- DMR and methylated Cytosines positions plots per chromosome

> [!IMPORTANT]
> This repo consist of five main folders:

> Scripts:
- 001_Genes_overlap.R
- 002_Join_transcriptome_epigenome.R
- 003_SHARED_GO_ontologies.R
- 004_key_genes_network.R
- DMR_CHR_POS_TOTAL.R
- LFC_mC.R
  
> Inputs:
 - Transcriptome: Differentially expressed genes obtained from DeSEQ2
 - Epigenome_DMRs_annotated: DMR
   
> Results:
 - OVERLAPS: overlap of genes, upstream and DMR
 - SHARED_GENES: overlap of DMR of genes, upstream regions and DEG
 - GSEA: GSEA results
 - GRAPHS: bipartite networks in graphml and PNG formats per rice variety
 - diff_meth_plots: Methylation position
 - LFC_mC_plots: Log2 fold change position and main differentially methylated windows for Azucena and BGI rice varieties

> [!CAUTION]
> Please follow the next order to reproduce the results:
- 001_Genes_overlap.R
- 002_Join_transcriptome_epigenome.R
- 003_SHARED_GO_ontologies.R
- 004_key_genes_network.R
- DMR_CHR_POS_TOTAL.R
- LFC_mC.R

## Citation 

> [!IMPORTANT]
> If you use this approach cite as follows:

Gallo-Franco JJ, Sosa CC, Johannes F, Ghneim-Herrera T and Quimbaya M (2026) Integrative methylome and transcriptome analysis reveals genotype and sequence context-specific responses to aluminum stress in rice. Front. Plant Sci. 17:1896320.  [doi: 10.3389/fpls.2026.1896320](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1896320/full?amp;utm_campaign=imp_impartaut-_05-24_fpls_en_n--ww&amp;utm_medium=social)

## Authors
- Bioinformatics: Chrystian Sosa and Jenny Gallo
- Conceptualization: Jenny Gallo and Mauricio Quimbaya
