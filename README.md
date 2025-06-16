# EasyVectorOmics

## General Information

EasyVectorOmics is a project currently funded by the Carl-Zeiss-Stiftung ("EasyVectorOmics") and aims at Multi Omics, and textual knowledge integration for high performance computing geometric Biomarker learning. Interactive scientific visualization for expert data exploration and hypothesis generation is supported by an graphical user user interface in the web-browser. 

## Abstract

Transcriptomics is the science of change in gene expression. Current methods provide key insights in a vast variety of life-science disciplines by robust identification of genes whose expression changes significantly between experimental conditions, like e.g. healthy versus wounded tissue. However, no tools yet exist to do comparisons between related genes or to efficiently analyze time series data. 
Tensor Omics is a new framework in applied bioinformatics for analyzing transcriptome data by projecting it into semantically meaningful vector spaces—where each axis represents a biological context such as tissue, dietary regimen, disease state, or developmental stage. Gene families are explored through expression shifts from ancestral centroids, enabling geometric insights into divergence, specificity, and adaptation. Through robust empirical outlier detection the method reliably captures biological signals like genes significantly contributing to wound healing or gut microbial response to protein based dietary regimens.

## Description
EasyVectorOmics is a pipeline designed to analyze gene expression data and evolutionary relationships between genes. This pipeline combines tools for phylogenetic analysis, gene classification, and expression calculations to generate processed results ready for interpretation.

The pipeline includes steps for data normalization, protein similarity calculations, phylogenetic classification, synteny analysis, and evolutionary distance and angle calculations.

### Analysis as a pipeline

Find below the summarized steps of the geometric biomarker detection implemented in Tensor Omics.

1. **Execute [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**

1. **Execute [Orthofinder](https://github.com/davidemms/OrthoFinder)**

1. **Gene expression values (TPM)**
   - TSV file in the format:  
   gene_id | Tissue_X | Tissue_Y | Tissue_Z | ...  

1. **Data Normalization**  
   - Script: `methods/normalization.R`  
   - Input: TPM expression file
   - Output: A normalized TPM expression matrix

2. **Protein Similarity Calculation**  
   - Script: `methods/harmonic_mean_calc.py`  
   - Input: The combined BLAST output file from OrthoFinder2 substitude back with original protein ids.
   - Output: A BLAST file with an added column containing the harmonic mean similarity between protein pairs.

3. **Gene Classification via Phylogeny**  
   - Scripts: `orthofinder_extensions/phylogeny/process_tree.py`, `orthofinder_extensions/phylogeny/tree_rest.py`  
   - Inputs for `process_tree.py`:  
        - Orthologues folder  (From orthofinder)
        - `harmonic_mean.pkl`  
        - `Duplications.tsv`  (From orthofinder)
        - Resolved gene trees folder  (From orthofinder)
        - `Orthogroups.tsv`  (From orthofinder)

   - `process_tree.py` analyzes resolved gene trees to classify gene relationships into conserved orthologs, inparalogs, outparalogs, source copy  inparalogs, and source copy outparalogs.  
   - `tree_rest.py` handles orthogroups not represented in the resolved trees and assigns classifications using remaining input data.

4. **Tandem Gene Classification for Synteny Analysis**  
   - Script: `methods/get_tandems.py`  
   - Input:  
        - All-vs-all BLASTn results  
        - Combined species GTF file   
   - Output: A TSV file listing tandem gene clusters, their species, and the best representative gene for each cluster.

5. **Neighborhood-Based Synteny Classification**  
   - Script: `orthofinder_extensions/synteny/synteny.py`  
   - Input:  
        - Species-specific GTF files  
        - Tandems file from `get_tandems.py`  
        - All-vs-all BLASTn results  
   - Output: A TSV file listing only source-to-source gene relationships.

6. **Gene Relationship Completion Based on Synteny**  
   - Script: `orthofinder_extensions/synteny/classify_genes.py`  
   - Input:  
        - `Orthogroups.tsv` from OrthoFinder  
        - Results from the neighborhood synteny classification  
   - Output: TSV files classifying remaining genes into conserved orthologs, inparalogs, outparalogs, special inparalogs, and special outparalogs.

7. **Filtering of Gene Families**  
   - Script: `methods/filter_orthogroups.py`  
   - Input:  
        - All output TSV files from `classify_genes.py`, or `process_tree.py`  
        - Orthogroups (From orthofinder)
   - Output: Filtered TSV files including only gene families with at least four unique conserved orthologs and copy gene relationships.

8. **Centroid Calculation for Expression Vectors**  
   - Script: `methods/get_centroids.R`  
   - Input:  
        - Normalized TPM value file  
        - Conserved orthologs TSV file  (filtered_orthologs.tsv from step 7)
   - Output: A TSV file containing orthogroup IDs and their expression centroids.

9. **Distance Calculation to Centroid**  
   - Script: `methods/calculate_distances.R`  
   - Input:  
        - Orthogroup centroids tsv file
        - Normalized TPM value file  
        - Conserved orthologs TSV file  (filtered_orthologs.tsv from step 7)
        - Conserved paralogs TSV file  (filtered_paralogs.tsv from step 7)
   - Output:  
        - `all_distances.tsv`: Euclidean distance of each gene to its group's centroid  
        - Source-copy pairwise distance TSV  
        - Outlier classification (top 5%)

10. **Tandem Annotation for Expression Analysis**  
    - Script: `methods/distances_add_tandems.py`  
    - Input:  
        - `Tandems.tsv` from `get_tandems.py`  
        - `distances.tsv` from `calculate_distances.R`  
        - `map.tsv` mapping Gene to Protein IDs  
    - Output: Updated `distances.tsv` annotated with tandem gene information.

11. **Evolutionary Angle Calculations**  
    - Script: `methods/calculate_angles_optimized.py`  
    - Input:  
        - Modified `distances.tsv` with tandem annotations  
        - Orthogroup centroids tsv file
        - Pairwise distance TSV  
    - Output: The input TSV files updated with angle calculations

12. **Functional Annotation and Word Cloud Generation (optional)**  
    - Script: `methods/word_cloud.R`  
    - Input:  
        - Outlier gene table (from `calculate_distances.R`)  
        - Functional annotations generated by [prot-scriber](https://github.com/usadellab/prot-scriber)  
    - Output:  
        - A word cloud image summarizing the most frequent terms in outlier gene annotations  
        - A TSV file with word frequencies and statistical enrichment results  

    - Description:  
        - [prot-scriber](https://github.com/usadellab/prot-scriber) assigns short human-readable descriptions (HRDs) to query biological sequences based on sequence similarity search results (e.g., BLAST or DIAMOND).  
        - The script generates a word cloud to visualize the most frequent terms associated with outlier genes, providing insights into their biological relevance.  
        - Statistical enrichment analysis is performed to identify terms significantly overrepresented in outlier genes compared to non-outliers.  
