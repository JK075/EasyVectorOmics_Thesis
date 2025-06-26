# EasyVectorOmics Analysis Pipeline – Documentation

This README documents the complete data processing and analysis performed within the EasyVectorOmics project. The aim was to utilize expression data from RNA-Seq experiments across multiple animal species to systematically analyze orthologous and paralogous relationships, and to compare their regulatory dynamics across different tissues.

All analyses were conducted on the bioinformatics server at TH Bingen. All scripts, raw data, and results are stored on this system. Precise file paths are specified in each section. The classfication algorithms can be found in orthofinder_extension. All other scripts are stored in methods.

---

## Origin and Content of the Datasets

The RNA-Seq raw data used originates from the GEO study [GSE125483](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE125483). This study includes 307 RNA-Seq samples across 13 tissues, collected from four vertebrate species: *Canis lupus familiaris* (dog), *Mus musculus* (mouse), *Rattus norvegicus* (rat), and *Macaca fascicularis* (a monkey species).

The FASTQ files were downloaded directly from NCBI and are stored under:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/Downloads/`.
A custom script `FastQ_downloader.py ` was made for this purpose.

In addition to the transcription data, the following reference data were downloaded from NCBI for each of the four species:

- Genome annotation in GTF format  
- Protein sequences (FASTA)  
- Transcript sequences (FASTA)  
- Coding sequences (CDS) (FASTA)

These datasets are based on the latest genome assemblies and are documented in the respective `README` files of the downloaded NCBI archives.
A comprehensive list of all genome resources used in this thesis is provided in supplementary_data/Links_IDs.txt. Each link points to the exact reference genome assembly utilized in the study. Within the “Download” section of each corresponding NCBI page, users can access the associated reference transcripts, coding sequences (CDS), protein FASTA files, as well as GTF annotation files and complete NCBI dataset packages that were used
throughout the thesis.

---

## Data Preprocessing: Quality Control and Quantification

After downloading the RNA-Seq data, adapter trimming and quality filtering were performed.

**Trimmomatic** was used for this step. The corresponding script is located at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/trimmomatic_kallisto_scripts/trimmomatic_skript.sh`.

Trimmed files were stored in the same directory as the original raw data. To distinguish them, a "T" was appended to the filenames (e.g., `SRR8474429_1T.fastq.gz`).

Subsequently, transcript quantification was performed using **Kallisto**. The related script is found at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/trimmomatic_kallisto_scripts/kallisto_real.sh`.

The reference transcriptomes used were the transcript FASTA files of each species, located at: ` /storage/EasyVectorOmics/FastQ_GSE125483_JK/reference_transcripts`. Kallisto’s output – particularly the `abundance.tsv` files – is located in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/kallisto_outputs/`.

---

## Protein Sequence Processing and Orthology Analysis

To retain only the longest protein isoform per gene, protein FASTA files were filtered based on the NCBI datasets containing protein lengths and gene IDs. This was performed using a custom Python script: `check_protein_gtf.py`.

The filtered protein FASTA files are stored in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/`.

The ncbi datasets are located in:
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/ncbi_datasets/`.


Next, **OrthoFinder2** was used to conduct a phylogenetically informed classification of orthologous genes. This analysis included pairwise BLASTP comparisons of protein sequences from all species and the construction of gene family trees.

The OrthoFinder launch script is located at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/orthofinder_job.sh`.  
It used the filtered protein FASTA files with only the longest protein isoform per gene as input.

The results – including orthologous groups, duplication events, gene trees, and additional outputs – are located in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/proteom/Orthofinder/`.

Additionally, a custom script (`harmonic_mean_calc`) was used to calculate the **harmonic mean** of similarity scores for each gene pair within the gene trees. This provided a more robust metric for defining orthologous and paralogous relationships. The script processed all BLASTP results into an all-vs-all format and computed the harmonic mean between each pair. The results are stored in:  
`/results/harmonic_means_by_prefix.pkl`.

Based on these metrics, a **tree-based classification** was performed using two key scripts:

- `process_tree.py`  
- `tree_rest.py`

`process_tree.py` combines data from the gene trees, OrthoGroups, duplication events, and harmonic mean scores to assign each gene relationship to a category (ortholog, inparalog, outparalog) for every resolved gene tree produced by OrthoFinder. Gene families too small to yield a resolved gene tree were handled by `tree_rest.py`.

The results are located at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tree/majority/`.

---

## Synteny-Based Classification

In addition to the phylogenetic approach, a second classification strategy was implemented based on **synteny** and genomic neighborhood.

First, CDS files were reduced to include only genes with available protein sequences for the longest isoform. This was done using the script `cds_modifier.py`, located at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/`.
The original cds files are located at:
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/genes/`.



The resulting FASTA file was saved at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/modified_cds_only_genid_matchwith_protein.fasta`.

These sequences were used in an all-vs-all **BLASTn** search, executed with the shell script `blast_all_vs_all.sh`. Results were stored in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/blast/`.

Subsequently, **tandem genes** were identified using the script `TandemsPaper.py`, which takes GTF files of all species as input. Tandem gene pairs are documented in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/tandems/Tandems_output.tsv`.

The classification of orthologous relationships was then carried out using `synteny.py`, which integrates data from the BLAST output, GTF files, and tandem results. The remaining paralogs were classified using `classify_genes.py`, which takes OrthoFinder orthogroups and uses the synteny-defined orthologs as anchors to identify linked paralog pairs.

The outputs were saved in:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/results/neighborhood/`.

---

## Data Consolidation and Preparation for EasyVectorOmics

The orthology relationships derived from synteny- and tree-based classification were unified using `result_formater.py`, keeping only gene pairs that are members of the same orthogroup. While tree-based classification inherently uses orthogroups, synteny-based classification applies them only during the paralog identification step.

Biologically meaningful gene families were then selected with `result_filtering.py`, retaining only orthogroups with at least four conserved gene relationships.

Next, the `abundance.tsv` files from Kallisto were mapped to protein information using `transcript_protein_kallisto_output_maper.py`. Expression values were then normalized using `normalizer.py`.

**Centroid calculation** for each gene family cluster in each tissue was performed using `get_centroids.R`, followed by outlier analysis to identify genes or gene pairs with expression patterns deviating strongly from the centroid. This was done using `calculate_distances.R`.

The expression distance tables and centroids are located in:

- Tree-based classification:  
  `/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results/Majority/`

- Synteny-based classification:  
  `/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results/Neighborhood/`

Additional annotations (e.g., whether a gene is part of a tandem) were added with `distances_add_tandems.py`. The final calculation of angles and distances – the core metrics of EasyVectorOmics – was implemented in `calculate_angles_optimized`.

---

## Functional Annotation and Enrichment Analysis

To explore functional patterns within the gene clusters, the web platform **PANNZER2** was used. Protein sequences (longest isoform per gene) were uploaded and annotated in "annotate" mode. The annotated results were downloaded and saved at:  
`/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results2/pannzer2_results.tsv`.

To visualize functional categories and enrichment patterns, the script `W2test.R` was used. The results can be found at:

- Tree-based:  
  `/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results/Majority/wordcloud/`

- Synteny-based:  
  `/storage/EasyVectorOmics/FastQ_GSE125483_JK/Scripts/EVO/results/Neighborhood/wordcloud/`

---

