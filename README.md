# Rice_aluminum_stress_methylated_transcriptome
This repository consist of R codes to overlap results of Azucena, and BGI rices varieties for epigenomics, and transcriptomics experiments all together to identify main key genes or mechanism implicated in aluminum stress tolerance.

## This approach consist of the following steps:
- It uses files in tabulate format with the Differentially Methylated Regions (DMR) for CG, CHG, and CHH contexts for the Azucena and BGI rice varieties and overlap with the genes, and upstream files from RAPD annotation database, available at https://rapdb.dna.naro.go.jp/download/irgsp1.html. Only genes with at least 60% coverage are considered
- It added the genes with absolute LFC >= 2 and put in a nested listed with the purpose of identify shared genes
- Performs gene enrichments analysis for the differentially methylated and upregulated or downregulated gene lists 
- Obtain a bipartite network with genes and methylated and transcripted features (Methylation context,and up or downregulated genes)

 
