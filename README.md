# BalticFishDB

## Abstract

BalticFishDB is designed for the taxonomic classification of mitochondrial 12S rRNA gene sequences from fish and is tailored exclusively to species from the Baltic Sea. It leverages sequences from public resources that were subject to strict curation, including similarity filtering, iterative phylogeny-based anomaly detection in full-length sequences (930-1030 bp) and phylogenetic placement of partial 12S sequences (<600 bp). The resulting database contains 1553 curated sequences of 253 species and enables efficient and reliable species assignments with modest computational resources using BLAST or DADA2. The focus on the Baltic Sea also excludes misassignments of closely related non-Baltic Sea species. BalticFishDB therefore provides a practical and robust reference for species identification in Baltic Sea fish metabarcoding studies. Additional to the database, we provide the bioinformatic workflow that can be used to updated the reference database. The workflow consists of interactive scripts that document the automatic and manual data retrieval and curation steps.


## Contents of the repository

**assets** 
* [balticfishdb_species_list.tsv](https://github.com/chassenr/BalticFishDB/blob/main/assets/balticfishdb_species_list.tsv): List of known Baltic Sea species (status: 2025-05-20). This list was compiled from [FishBase](https://www.fishbase.se/identification/RegionSpeciesList.php?e_code=104) and the [HELCOM Checklist 2.0 of Baltic Sea Macrospecies](https://helcom.fi/publications/helcom-checklist-2-0-of-baltic-sea-macrospecies/).
* Sequence accession blacklists: During the compilation of the database three manually curated blacklists were used: [exclusions.accnos](https://github.com/chassenr/BalticFishDB/blob/main/assets/exclusions.accnos) is based on the list of accessions with presumably erroneous species assignments in [Meta-Fish-Lib](https://onlinelibrary.wiley.com/doi/10.1111/jfb.14852), [anomalies.accnos](https://github.com/chassenr/BalticFishDB/blob/main/assets/anomalies.accnos) was curated during this workflow based on the output of the [PhyloRef anomaly detection module](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.73159) that inspects taxonomic consistency based on tree topology, [placement.accnos](https://github.com/chassenr/BalticFishDB/blob/main/assets/placement.accnos) contains accessions of sequences shorter than 600bp that were discarded based on an anomalous phylogenetic placement in the tree of full-length 12S sequences due to insufficient region coverage or presumably erroneous species assignments.

**scripts**
* [sequence_selection.sh](https://github.com/chassenr/BalticFishDB/blob/main/scripts/sequence_selection.sh): Interactive bash script for the main workflow, detailing sequence retrieval, customization of sequence selection for the Baltic Sea, curation based on alignment and phylogeny, placement of short sequences to increase species coverage, formatting of final sequence selection for use as reference database with blast and dada2.
* [filter_alignment_sequence_similarity.R](https://github.com/chassenr/BalticFishDB/blob/main/scripts/filter_alignment_sequence_similarity.R): Helper script to calculate sequence similarity based on multiple sequence alignments and to filter sequences based on alignment similarity.
* [reference_db_insilico_pcr.sh](https://github.com/chassenr/BalticFishDB/blob/main/scripts/reference_db_insilico_pcr.sh): Interactive bash script to run in-silico PCR on the reference sequence database and assess species coverage for popular 12S primer pairs.
* [validate_species_gap.R](https://github.com/chassenr/BalticFishDB/blob/main/scripts/validate_species_gap.R): Helper script to calculate sequence similarity based on multiple sequence alignments in order to visualize and quantify the barcode gap per species for popular 12S primer pairs.


## Workflow description

Each step is extensively commented in [sequence_selection.sh](https://github.com/chassenr/BalticFishDB/blob/main/scripts/sequence_selection.sh). Here, we summarized the main analysis steps to guide users through the pipeline.

<img src="https://github.com/chassenr/BalticFishDB/blob/main/images/sequence_selection_overview.png" width="100%" height="100%">

### Preparation of working environment

The workflow was developed and run on Ubuntu 22.04.5 LTS.

The following dependencies are needed to run the workflow. The conda evironments, in which these dependecies are available on the system the workflow was developed on, are indicated in parentheses. Install instructions for these environments can be found [here](https://github.com/chassenr/BalticFishDB/blob/main/envs/install_instructions.sh).

* mafft=7.525 (align-env)
* entrez-direct=24.0 (antismash-8.0.4)
* blast=2.16.0 and fasttree=2.1.11 (anvio-8)
* biopython=1.87 and ete3=3.1.3 for anomaly detection (meme-5.5.9)
* raxml-ng=1.2.2 (orthofinder-2.5.5)
* epa-ng=0.3.8 and gappa=0.8.5 (paprica-env)
* seqkit=2.8.2 (seqkit-2.8.2)
* taxonkit=0.20.0 (taxonkit-0.20.0)

For manual curation, the following programs were used:

* JalView 2.11.5.2
* RStudio 2026.05.0 with R 4.6.0 and the packages: seqinr 4.2.44, stringr 1.6.0, reshape 0.8.10, tidyverse 2.2.0, scales 1.4.0 

The [script for sequence selection and curation](https://github.com/chassenr/BalticFishDB/blob/main/scripts/sequence_selection.sh) is using environment variables for the location of the main working directors and the location of the repository. Change those to adapt them to your system.


### Retrieve available sequence from NCBI, MIDORI2, PhyloRef

While both PhyloRef and MIDORI2 are based on sequences available on NCBI, their update frequency is not keeping pace with new, potentially highly relevant sequence submissions. We therefore retrieve sequences from all three sources. However, since the amount of manual curation differs between the three sources, if there are any duplicate accession numbers, the sequences in PhyloRef are preferred over those in MIDORI2 and both outrank NCBI. Furthermore, we consider prior efforts in 12S reference database curation ([Meta-Fish-Lib](https://github.com/genner-lab/meta-fish-lib)) by considering their sequence accession blacklists. We are also restricting sequences to fish species known to occur in the Baltic Sea.


### Alignment of full length (or near to full length) 12S sequences

Based on sequences between 930bp and 1030bp long of the PhyloRef and MIDORI2 selection, sequences outside this range are trimmed to the 12S region using blast. Confirmed 12S sequences >= 600bp are then aligned and the alignment is trimmed to the 12S region and optionally filtered for divergend sequences within species.


### Anomaly detection

The multiple sequence alignment (MSA) from the previous step is used to reconstruct a phylogenetic tree to run the [anomaly detection module](https://github.com/yannnnmai/PhyloRef/blob/main/workflow/scripts/detect_anomaly.py) from [PhyloRef](https://github.com/yannnnmai/PhyloRef). Angelsharks (*Squatina squatina*) are used as outgroup. Placements of putative anomalous sequences are manually inspected in the tree and a blacklist is curated based on multiple iterations of th anomaly detection. After removal of the blacklisted sequences, the alignment is dereplicated to remove redundant sequences per species. From this sequence selection a BLAST database is constructed to further screen putative 12S sequences shorter than 600bp that were previously not processed.


### Phylogenetic placement of short sequences

Sequences shorter than 600bp are dereplicated within species and BLASTed against the curated 12S alignment. Those identical to sequences in the alignment and those diverging more than 90% from sequences in the alignment of the same species are discarded. The remaining short sequences together with the curated sequences from the alignment are re-aligned together and the aligned is trimmed again to the 12S region. Based on the new MSA, the short sequences (<600bp) are placed into a phylogenetic tree built from the long 12S sequences (>=600bp). The placement of short sequences of species not previously represented among the long sequences was manually inspected and anomalous sequences were blacklisted. Furthermore, short sequences less than 95% similar to long sequences of the same species or those with a better alignment to a long sequence of another species were removed. 


### Finalizing the 12S reference DB

The final 12S reference sequence database consists of the non-redundant sequences that were dereplicated per species. This way potential taxonomic ambiguities are retained while reducing the size of the database. The following files are available for downstream applications:

* reference_12s_aligned.fasta: aligned 12S reference sequence database, fasta headers include: accession|source|rc (reverse complement, if applicable)|path (no spaces)
* reference_12s_seqs.fasta: ungapped 12S reference sequence database, same fasta header as used in the aligned sequence file
* reference_12s_dada2.fasta: ungapped 12S reference sequence database, fasta headers formatted for use with dada2:assignTaxonomy (ranks: kingdom, phylum, class, order, family, genus species)
* reference_12s_blastdb.*: BLAST database


