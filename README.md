# Python scripts and command-line workflows used in research-based master's thesis on *Vibrio natriegens*

This repository contains all Python scripts and command-line workflows used for all computational analyses conducted in the research-based master's thesis regarding the exploration of the *Vibrio natriegens* pangenome.

## Main Objectives

1.	Delineate the pangenome of *V. natriegens*
2.	Identify homologs of proteins encoded by strain-specific *V. natriegens* genes in eukaryotic, bacterial, archaeal and viral domains and determine the taxonomic distribution of the corresponding matches.
3.	Compile existing Gene Ontology (GO) functional annotations for the *V. natriegens* pangenome.
4.	Identify patterns of amino acid sequence composition and predicted protein structural features that potentially distinguish strain-specific genes within the *V. natriegens* pangenome.
5.	Enrich the functional characterisation of strain-specific *V. natriegens* genes through computationally predicted secretory signal peptides and their cleavage sites.

## Repository Structure

- 'complementary-processing/': supplementary scripts aiding in the analysis of the pangenome.
- 'mmseqs2-analysis/': processing of sequences used as input for MMseqs2 and results summary production.
- 'taxonomic-analysis/': taxonomic classification and summarisation of detectable homologous matches.
- 'gene-ontology-analysis/': processing of genes prior to analysis and fetching of gene ontology annotations. 
- 'protein-features/': processing and summarisation of amino acid composition and predicted protein structural features.
- 'clustering-analysis/': statistical significance tests for feature inclusion, determination of optimal k values, clustering implementation and results summarisation.
- 'signal-peptide-analysis/': processing and summarisation of signal peptide and cleavage site predictions.

## Requirements

pandas
tqdm
numpy
biopython
openpyxl
requests
matplotlib
scikit-learn
seaborn
scipy
statsmodels

## Notes

Local file paths and file names used during the analyses were replaced with placeholders.
