
# Set working directory to download and curate the reference sequences

WDIR="/srv/bio/Analysis_data/IOWseq000085_KOFI/Intermediate_results/ReferenceDB_12S_mifish"
REPO="/srv/bio/Analysis_data/IOWseq000085_KOFI/Intermediate_results/Repos/BalticFishDB"


### Get available blacklist accession (don't run)

# This section of code retrieves and format the blacklist of sequence accessions provided by Collins et al. (2021) as part of their Meta-Fish-Lib workflow.
# The formatted blacklist is already provided as part of this repository: assets/exclusions.accnos

cd $WDIR
mkdir blacklist
cd blacklist
wget https://raw.githubusercontent.com/genner-lab/meta-fish-lib/refs/heads/main/assets/exclusions.csv
conda activate antismash-8.0.4
sed '1d' exclusions.csv | grep -v '^BOLD' | cut -d',' -f2 | while read gi
do
  efetch -db nucleotide -id "$gi" -format acc 2>/dev/null | sed "s/^/$gi\t/"
done | cut -d'.' -f1 > gi2acc.txt
# There is 1 entry in the exclusion list that is not a gi number: HVDBF493-12. This will be skipped.

# Include corresponding genbank accession for refseq entries, to ensure that all blacklist sequences are removed regardless of their source database.
cut -f2 gi2acc.txt | grep "_" | while read acc
do
  efetch -db nucleotide -id "$acc" -format gb 2>/dev/null | \
  awk '
    /^COMMENT/ { capture = 1; print; next }
    capture && /^[^A-Z]/ { print; next }
    capture && /^[A-Z]/ { capture = 0 }
  ' | \
  grep -oE '[A-Z]{2}[0-9]{4,8}' | sed "s/^/$acc\t/"
done > rs2gb.txt
conda deactivate

# Compile list with all exlcusion accessions
awk '$2 != ""' gi2acc.txt | cut -f2 | cat - <(cut -f2 rs2gb.txt) > exclusions.accnos


### Download available sequences for the Baltic Sea from NCBI 

cd $WDIR
mkdir sequences_ncbi
cd sequences_ncbi
mkdir by_species

# For each valid species name, retrieve all ribosomal gene sequences from the mitochondrion between 150 and 20000bp length
conda activate antismash-8.0.4
cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | while read species
do
  esearch -db nuccore -query "\"$species\"[Organism] AND mitochondrion[filter] AND ribosomal[All] AND (\"150\"[SLEN] : \"20000\"[SLEN])" < /dev/null | efetch -format fasta > by_species/"${species// /_}.fasta"
done
conda deactivate
# The above code was contributed by chatGPT
# To avoid an ambiguous redirect, /dev/null is required
# Efetch returned an error: 
#   curl: (56) OpenSSL SSL_read: OpenSSL/3.6.2: error:0A000126:SSL routines::unexpected eof while reading, errno 0
# However, since this error was followed by: HTTP/1.0 200 OK, it was ignored

# Format fasta headers, removing everything after the first space
# Collect sequences from all species into one fasta file
cd by_species
ls -1 *.fasta | sed 's/\.fasta//' | while read line
do
  sed "/^>/s/ .*$/|$line/" $line.fasta
done > ../ncbi_baltic.fasta
cd ..

# Write sequence accessions to a separate file
grep '^>' ncbi_baltic.fasta | sed 's/^>//' | cut -d'|' -f1 | cut -d'.' -f1 > ncbi_baltic.accnos
cd ..


### Retrieve full length 12S from mitochondrial genomes in MIDORI2

cd $WDIR
mkdir sequences_midori
cd sequences_midori
wget https://www.reference-midori.info/download/Databases/GenBank271_2026-04-07/RDP/uniq/MIDORI2_UNIQ_NUC_GB271_srRNA_RDP.fasta.gz
# To update the MIDORI version, visit their download website: https://www.reference-midori.info/download.php, and select the newest genbank release version
# Navigate to the same data type: MIDORI formatted for RDP containing all unique sequences of the small ribosomal subunit

# Write midori accessions for Baltic Sea species
# Remove hybrid species
cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/$/_/' | zgrep -f - MIDORI2_UNIQ_NUC_GB271_srRNA_RDP.fasta.gz | grep -v " x " | cut -f1 | sed 's/^>//' | cut -d'.' -f1 > midori_baltic.accnos
cd ..


### Retrieve full length 12S from PhyloRef

cd $WDIR
mkdir sequences_phyloref
cd sequences_phyloref

# Actinopterygii downloaded from https://zenodo.org/records/17285318 (date accessed 2026-06-29)
wget https://zenodo.org/records/17285318/files/Actinopterygii_data_v2025Feb.zip?download=1 -O Actinopterygii_data_v2025Feb.zip
unzip Actinopterygii_data_v2025Feb.zip

# Remove unnecessary files
rm -rf __MACOSX
rm Actinopterygii_data_v2025Feb.zip

# Format fasta header
# This assumes that every sequence header contains "|12S|" and that species und accession number are listed before
sed '/^>/s/|12S|.*$//' Actinopterygii_data_v2025Feb/Actinopterygii_cleaned.fa_V2502/12S_cleaned.fa > phyloref_12s.fasta

# Get sequence accessions for Baltic Sea species
# Remove hybrid species
cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/ /_/g' | grep -w -f - phyloref_12s.fasta | grep -v " x " | cut -d'|' -f2 | cut -d'.' -f1 > phyloref_baltic.accnos

# Retrieve genbank accession for refseq entries for compatibility with MIDORI2 to avoid sequence duplication
# Code contributed by GLM-fast 5.2
conda activate antismash-8.0.4
grep "_" phyloref_baltic.accnos | while read acc
do
  efetch -db nucleotide -id "$acc" -format gb 2>/dev/null | \
  awk '
    /^COMMENT/ { capture = 1; print; next }
    capture && /^[^A-Z]/ { print; next }
    capture && /^[A-Z]/ { capture = 0 }
  ' | \
  grep -oE '[A-Z]{2}[0-9]{4,8}' | sed "s/^/$acc\t/"
done > refseq2genbank.accnos
conda deactivate
cd ..


### Remove redundancy between the 3 data sources

cd $WDIR

# Ranking of data sources: PhyloRef over MIDORI2, PhyloRef and MIDORI2 over NCBI
# At this stage, blacklist accessions according to Meta-Fish-Lib will be removed

# Select non-blacklist accessions from PhyloRef
grep -v -w -f $REPO/assets/exclusions.accnos sequences_phyloref/phyloref_baltic.accnos > sequences_phyloref/phyloref_baltic_select.accnos

# Select non-PhyloRef and non-blacklist accessions from MIDORI
cut -f2 sequences_phyloref/refseq2genbank.accnos | cat - sequences_phyloref/phyloref_baltic.accnos | grep -v -w -f - sequences_midori/midori_baltic.accnos | grep -v -w -f $REPO/assets/exclusions.accnos > sequences_midori/midori_baltic_select.accnos

# Select non-PhyloRef, non-MIDORI and non-blacklist accessions from NCBI
cat sequences_phyloref/refseq2genbank.accnos | tr '\t' '\n' | cat - sequences_phyloref/phyloref_baltic.accnos sequences_midori/midori_baltic.accnos | grep -v -w -f - sequences_ncbi/ncbi_baltic.accnos | grep -v -w -f $REPO/assets/exclusions.accnos > sequences_ncbi/ncbi_baltic_select.accnos

# Exceptions (manually confirmed):
# AH013020: take sequence from NCBI, not MIDORI2 (this sequences is split in two in MIDORI)
echo "AH013020" >> sequences_ncbi/ncbi_baltic_select.accnos
sed '/AH013020/d' -i sequences_midori/midori_baltic_select.accnos

# Subset fasta files to selected accessions
# Format fasta header consistently: accnos|species|source
# Remove sequences with ambiguous bases 
# Collect short sequences (<600bp) in a separate file 
# Add data source to fasta header

conda activate seqkit-2.8.2

cd sequences_phyloref
seqkit seq -w 0 phyloref_12s.fasta | grep -A1 -w -f phyloref_baltic_select.accnos | sed '/^--$/d' | sed -E '/^>/s/>([^|]*)\|([^|]*)/>\2|\1/' | sed -e '/^>/s/\..*|/|/' -e '/^>/s/$/|phyloref/' | seqkit seq -m 600 | seqkit grep -w 0 -s -v -r -p [^ATCG] > phyloref_baltic_select.fasta
seqkit seq -w 0 phyloref_12s.fasta | grep -A1 -w -f phyloref_baltic_select.accnos | sed '/^--$/d' | sed -E '/^>/s/>([^|]*)\|([^|]*)/>\2|\1/' | sed -e '/^>/s/\..*|/|/' -e '/^>/s/$/|phyloref/' | seqkit seq -M 599 | seqkit grep -w 0 -s -v -r -p [^ATCG] > phyloref_baltic_select_short.fasta

cd ../sequences_midori
seqkit seq -w 0 MIDORI2_UNIQ_NUC_GB271_srRNA_RDP.fasta.gz | grep -A1 -w -f midori_baltic_select.accnos | sed '/^--$/d' | sed -e '/^>/s/\t.*;/|/' -e '/^>/s/_[0-9]*$//' -e '/^>/s/ /_/' -e '/^>/s/\..*|/|/' -e '/^>/s/$/|midori/' | seqkit seq -m 600 | seqkit grep -w 0 -s -v -r -p [^ATCG] > midori_baltic_select.fasta
seqkit seq -w 0 MIDORI2_UNIQ_NUC_GB271_srRNA_RDP.fasta.gz | grep -A1 -w -f midori_baltic_select.accnos | sed '/^--$/d' | sed -e '/^>/s/\t.*;/|/' -e '/^>/s/_[0-9]*$//' -e '/^>/s/ /_/' -e '/^>/s/\..*|/|/' -e '/^>/s/$/|midori/' | seqkit seq -M 599 | seqkit grep -w 0 -s -v -r -p [^ATCG] > midori_baltic_select_short.fasta

cd ../sequences_ncbi
seqkit seq -w 0 ncbi_baltic.fasta | grep -A1 -w -f ncbi_baltic_select.accnos | sed '/^--$/d' | sed -e '/^>/s/\..*|/|/' -e '/^>/s/$/|ncbi/' | seqkit seq -m 600 | seqkit grep -w 0 -s -v -r -p [^ATCG] > ncbi_baltic_select.fasta
seqkit seq -w 0 ncbi_baltic.fasta | grep -A1 -w -f ncbi_baltic_select.accnos | sed '/^--$/d' | sed -e '/^>/s/\..*|/|/' -e '/^>/s/$/|ncbi/' | seqkit seq -M 599 | seqkit grep -w 0 -s -v -r -p [^ATCG] > ncbi_baltic_select_short.fasta
cd ..
conda deactivate

# Sanity check:
# How many species are represented in each data source
grep '^>' sequences_phyloref/phyloref_baltic_select.fasta | cut -d'|' -f2 | sort | uniq -c | wc -l
# phyloref: 139
grep '^>' sequences_midori/midori_baltic_select.fasta | cut -d'|' -f2 | sort | uniq -c | wc -l
cat sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | grep '^>' | cut -d'|' -f2 | sort | uniq -c | wc -l
# midori: 201 (226 in total)
grep '^>' sequences_ncbi/ncbi_baltic_select.fasta | cut -d'|' -f2 | sort | uniq -c | wc -l
cat sequences_ncbi/ncbi_baltic_select.fasta sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | grep '^>' | cut -d'|' -f2 | sort | uniq -c | wc -l
# ncbi: 150 (233 in total)
cat sequences_phyloref/phyloref_baltic_select_short.fasta sequences_midori/midori_baltic_select_short.fasta sequences_ncbi/ncbi_baltic_select_short.fasta | grep '^>' | cut -d'|' -f2 | sort | uniq -c | wc -l
cat sequences_ncbi/ncbi_baltic_select*.fasta sequences_phyloref/phyloref_baltic_select*.fasta sequences_midori/midori_baltic_select*.fasta | grep '^>' | cut -d'|' -f2 | sort | uniq -c | wc -l
# short sequences: 254 (258 in total)

# Are there more species than expected according to the species list: no
cat sequences_ncbi/ncbi_baltic_select*.fasta sequences_phyloref/phyloref_baltic_select*.fasta sequences_midori/midori_baltic_select*.fasta | grep '^>' | cut -d'|' -f2 | sort -u | grep -v -f <(cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/ /_/')

# Which species are missing when only considering long sequences
cat sequences_ncbi/ncbi_baltic_select.fasta sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | grep '^>' | cut -d'|' -f2 | sort -u | grep -v -f - <(cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/ /_/')
# Argentina_sphyraena
# Babka_gymnotrachelus
# Ballerus_ballerus
# Callionymus_lyra
# Centrolophus_niger
# Cheilopogon_heterurus
# Chelidonichthys_lastoviza
# Ciliata_septentrionalis
# Coregonus_pallasii # also missing in short sequence selection
# Dentex_maroccanus
# Halobatrachus_didactylus
# Lebetus_guilleti
# Lebetus_scorpioides # also missing in short sequence selection
# Liparis_barbatus # also missing in short sequence selection
# Liparis_liparis
# Liparis_montagui
# Lycodes_gracilis
# Misgurnus_fossilis
# Orcynopsis_unicolor # also missing in short sequence selection
# Oxynotus_centrina # also missing in short sequence selection
# Phycis_blennoides
# Platichthys_solemdali # also missing in short sequence selection
# Pomatoschistus_norvegicus
# Pomatoschistus_pictus
# Proterorhinus_marmoratus
# Pterycombus_brama
# Sabanejewia_aurata
# Taractichthys_longipinnis
# Tetronarce_nobiliana
# Vimba_vimba
# Zeugopterus_norvegicus


### Alignment of full length sequences

# Identify 12S in sequences from ncbi and those outside of expected range in midori and phyloref
# Midori and phyloref sequences within the range of 930 and 1030 should represent true 12S

cd $WDIR
mkdir identify_12s
conda activate seqkit-2.8.2

# Concatenate long sequences from all three data sources and filter to those of 930 - 1030bp length
cat sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | seqkit seq -w 0 -m 930 -M 1030 > identify_12s/seed_12s.fasta

# Select all sequences >= 600bp (i.e. not from *short.fasta) that are not part of the seed selection as candidate 12S
cat sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | seqkit seq -w 0 -M 929 > identify_12s/candidate_12s.fasta
cat sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | seqkit seq -w 0 -m 1031 >> identify_12s/candidate_12s.fasta
cat sequences_ncbi/ncbi_baltic_select.fasta >> identify_12s/candidate_12s.fasta
conda deactivate

# Build blast database
cd identify_12s
conda activate anvio-8
makeblastdb -in seed_12s.fasta -dbtype nucl -out seed_12s_blastdb -logfile seed_12s_makeblastdb.log

# Use blast to identify overlap with 12S region
# Set max_target_seqs to include all seed sequences (just in case)
# Include qseq to get aligned part of query sequence
blastn -query candidate_12s.fasta -task blastn -db seed_12s_blastdb -out candidate_12s.blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen qseq" -num_threads 32 -max_target_seqs 1000
conda deactivate

# Assuming that there is only one 12S gene per mitochondrial genome, select the best blast hit based on bitscore to obtain gene coordinates
# Useful link: https://mitofish.aori.u-tokyo.ac.jp/about/overview
# Also remove alignments shorter than 150bp, the minimum length for the *short.fasta selection
sort -k1,1 -k12,12gr candidate_12s.blastout | sort -u -k1,1 --merge | awk '$4 >= 150' > candidate_12s.blastout_best

# Extract identified 12S region from blast result (qseq)
# Separate between long and short (<600bp) 12S sequences
conda activate seqkit-2.8.2
cut -f1,16 candidate_12s.blastout_best | sed 's/^/>/' | tr '\t' '\n' | seqkit seq -w 0 -g -m 600 > confirmed_12s.fasta
cut -f1,16 candidate_12s.blastout_best | sed 's/^/>/' | tr '\t' '\n' | seqkit seq -w 0 -g -M 599 > confirmed_12s_short.fasta
cd ..
conda deactivate

# Prepare long (>=600bp) 12S sequences from all 3 data sources for alignment
mkdir align_12s
cat identify_12s/seed_12s.fasta identify_12s/confirmed_12s.fasta > align_12s/unaligned_12s.fasta

# Sanity check
grep '^>' align_12s/unaligned_12s.fasta | cut -d'|' -f2 | sort | uniq -c | wc -l
# 228 sepcies represented
grep '^>' align_12s/unaligned_12s.fasta | cut -d'|' -f2 | sort -u | grep -v -f - <(cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/ /_/')
# missing species:
# Argentina_sphyraena
# Babka_gymnotrachelus
# Ballerus_ballerus
# Brama_brama # lost during 12S validation step
# Callionymus_lyra
# Centrolophus_niger
# Cheilopogon_heterurus
# Chelidonichthys_lastoviza
# Ciliata_septentrionalis
# Coregonus_pallasii
# Coris_julis # lost during 12S validation step
# Dentex_maroccanus
# Halobatrachus_didactylus
# Lebetus_guilleti
# Lebetus_scorpioides
# Liparis_barbatus
# Liparis_liparis
# Liparis_montagui
# Lycodes_gracilis
# Misgurnus_fossilis
# Orcynopsis_unicolor
# Oxynotus_centrina
# Phycis_blennoides # lost during 12S validation step
# Platichthys_solemdali
# Pomatoschistus_norvegicus
# Pomatoschistus_pictus
# Proterorhinus_marmoratus
# Pterycombus_brama
# Sabanejewia_aurata
# Sarpa_salpa # lost during 12S validation step
# Symphodus_melops # lost during 12S validation step
# Taractichthys_longipinnis
# Tetronarce_nobiliana
# Thorogobius_ephippiatus # lost during 12S validation step
# Vimba_vimba
# Zeugopterus_norvegicus

# Alignment with mafft (default settings)
cd align_12s
conda activate align-env
mafft --adjustdirection unaligned_12s.fasta > aligned_12s.fasta
conda deactivate

# The alignment was visually inspected in JalView to determine start/end positions for 12s

# Trimming alignment to fixed start/end position
# Convert sequence to upper case
conda activate seqkit-2.8.2
seqkit subseq -r 53:1213 aligned_12s.fasta | seqkit seq -w 0 -u > aligned_12s_trimmed.fasta
conda deactivate

# No dereplication was performed on purpose to retain the information about the number of sequences per species
# Compare sequence similarity within species: filter_alignment_sequence_similarity.R

# Subset aligned sequences to those surviving similarity filter (optional)
sed '1d' aligned_12s_similarity_check.txt | awk '$6 == "FALSE"' | cut -f3 | seqkit grep -n -f - aligned_12s_trimmed.fasta > aligned_12s_filtered.fasta
# If all sequences should be kept
cp aligned_12s_trimmed.fasta aligned_12s_filtered.fasta
# Only run one of the above 2 lines of code


### Anomaly detection

cd $WDIR/align_12s

# Manual curation of sequence selection based on phylogenetic tree
# Consider source of sequence (trust phyloref over midori oder ncbi)
# Using PhyloRef detect_anomaly.py script
wget https://raw.githubusercontent.com/yannnnmai/PhyloRef/refs/heads/main/workflow/scripts/detect_anomaly.py
chmod +x detect_anomaly.py

# Recreate phyloref header to use the detect_anomaly script: species|accession|source|o__order|f__family|g__genus
conda activate taxonkit-0.20.0
join -t $'\t' -1 1 -2 1 <(grep '^>' aligned_12s_filtered.fasta | sed 's/^>//' | cut -d'|' -f2 | paste - <(grep '^>' aligned_12s_filtered.fasta | sed 's/^>//') | sort -k1,1) <(grep '^>' aligned_12s_filtered.fasta | sed 's/^>//' | cut -d'|' -f2 | sort -u | sed 's/_/ /' | taxonkit name2taxid --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | taxonkit reformat2 -I 2 --format "o__{order}|f__{family}|g__{genus}" --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | cut -f1,3 | sed 's/ /_/' | sort -k1,1) > tmp
conda deactivate
conda activate seqkit-2.8.2
cut -f2 tmp | paste - <(cut -f1 tmp) | paste -d'|' - <(cut -f2 tmp | cut -d'|' -f1) | paste -d'|' - <(cut -f2 tmp | cut -d'|' -f3) | paste -d'|' - <(cut -f3 tmp) | seqkit replace -w 0 -p '^(.+)$' -r '{kv}' -k - aligned_12s_filtered.fasta > aligned_12s_renamed.fasta
conda deactivate
rm tmp

# Calculate tree
conda activate anvio-8
FastTree -log aligned_12s_fasttree.log -gtr -gamma -nt aligned_12s_renamed.fasta > aligned_12s_fasttree.tree
conda deactivate

# Get outgroup accessions to be used in the anomaly detection (-og)
grep "Squatina_squatina" aligned_12s_filtered.fasta | cut -d'|' -f1 | sed 's/^>//' | tr '\n' ',' | sed 's/,$/\n/'

# Run phyloref anomaly detection
conda activate meme-5.5.9
mkdir anomaly_iter1
python detect_anomaly.py -i aligned_12s_fasttree.tree -o anomaly_iter1/out -og NC_035057,MH166797,MG029174,MH166792,MH166793,MH166794,MH166795,MH166796,MH166798
conda deactivate
# Manually curate anomaly detection output and save accession numbers to be discarded in a new blacklist: anomalies.accnos
# This blacklist will be extended with each iteration of the anomaly detection step
# The curated blacklist is provided on the repositoriy under assets

# Subset alignment, repeat tree reconstruction, repeat anomaly detection
conda activate seqkit-2.8.2
sed -e 's/^/\\|/' -e 's/$/\\|/' $REPO/assets/anomalies.accnos | seqkit grep -w 0 -v -r -f - aligned_12s_renamed.fasta > aligned_12s_clean.fasta
conda deactivate
conda activate anvio-8
FastTree -log aligned_12s_fasttree_clean.log -gtr -gamma -nt aligned_12s_clean.fasta > aligned_12s_fasttree_clean.tree
conda deactivate
conda activate meme-5.5.9
mkdir anomaly_iter2 # manually increase counter in each iteration
python detect_anomaly.py -i aligned_12s_fasttree_clean.tree -o anomaly_iter2/out -og NC_035057,MH166797,MG029174,MH166792,MH166793,MH166794,MH166795,MH166796,MH166798
conda deactivate
# Repeat until no further anomalies occur based on manual inspection

# Once the tree is clean, dereplicate sequences within species to reduce the size of the reference sequence database
mkdir derep_logs
conda activate seqkit-2.8.2
grep '^>' aligned_12s_clean.fasta | cut -d'|' -f1 | sed 's/^>//' | sort -u | while read line
do
  grep -A1 -w "$line" aligned_12s_clean.fasta | sed 's/^--$//' | seqkit rmdup -w 0 -s -D derep_logs/${line}_nr.log
done > aligned_12s_nr_spec.fasta

# Build blast db to screen short sequences
seqkit seq -g -w 0 aligned_12s_nr_spec.fasta > seqs_12s_nr_spec.fasta
conda deactivate
conda activate anvio-8
makeblastdb -in seqs_12s_nr_spec.fasta -dbtype nucl -out seqs_12s_nr_spec_blastdb -logfile seqs_12s_nr_spec_makeblastdb.log
conda deactivate
cd ..


### Phylogenetic placement of short sequences

cd $WDIR

# Collect short sequences, incl. sequences originally part of long selection but moved to short after cutting to 12S region
# Reapply filter to sequences longer than 150bp as absolute minimum length
mkdir placement_12s
conda activate seqkit-2.8.2
cat identify_12s/confirmed_12s_short.fasta sequences_midori/midori_baltic_select_short.fasta sequences_ncbi/ncbi_baltic_select_short.fasta sequences_phyloref/phyloref_baltic_select_short.fasta | seqkit seq -m 150 -w 0 > placement_12s/short_12s.fasta
conda deactivate
cd placement_12s

# Recreate phyloref header in case we want to use the detect_anomaly script again and to be consistent with the header formatting of the long sequences
conda activate taxonkit-0.20.0
join -t $'\t' -1 1 -2 1 <(grep '^>' short_12s.fasta | sed 's/^>//' | cut -d'|' -f2 | paste - <(grep '^>' short_12s.fasta | sed 's/^>//') | sort -k1,1) <(grep '^>' short_12s.fasta | sed 's/^>//' | cut -d'|' -f2 | sort -u | sed 's/_/ /' | taxonkit name2taxid --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | taxonkit reformat2 -I 2 --format "o__{order}|f__{family}|g__{genus}" --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | cut -f1,3 | sed 's/ /_/' | sort -k1,1) > tmp
conda deactivate
conda activate seqkit-2.8.2
cut -f2 tmp | paste - <(cut -f1 tmp) | paste -d'|' - <(cut -f2 tmp | cut -d'|' -f1) | paste -d'|' - <(cut -f2 tmp | cut -d'|' -f3) | paste -d'|' - <(cut -f3 tmp) | seqkit replace -w 0 -p '^(.+)$' -r '{kv}' -k - short_12s.fasta > short_12s_renamed.fasta
rm tmp

# Dereplicate short sequences within species
mkdir derep_logs_short
grep '^>' short_12s_renamed.fasta | cut -d'|' -f1 | sed 's/^>//' | sort -u | while read line
do
  grep -A1 -w "$line" short_12s_renamed.fasta | sed 's/^--$//' | seqkit rmdup -w 0 -s -D derep_logs_short/${line}_nr.log
done > short_12s_nr_spec.fasta
conda deactivate

# Blast short sequences against curated selection of non-redundant (within species) long sequences
# Set max_target_seqs to include all long sequences (just in case)
# Include qseq to get aligned part of query sequence (just in case)
conda activate anvio-8
blastn -query short_12s_nr_spec.fasta -task blastn -db ../align_12s/seqs_12s_nr_spec_blastdb -out short_12s_nr_spec.blastout -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qlen slen qseq" -num_threads 38 -max_target_seqs 2000
conda deactivate

# Within species, remove short sequences that are identical across their full length to sequences already present in the long sequence selection: short_identical.accnos
# Within species, remove short sequences that have a high query coverage and low percent identity: short_divergent.accnos
# see: filter_alignment_sequence_similarity.R
conda activate seqkit-2.8.2
cat short_identical.accnos short_divergent.accnos | seqkit grep -w 0 -n -v -f - short_12s_nr_spec.fasta > short_12s_select.fasta
conda deactivate
# 1935 sequences remaining

# We tested three different approaches for aligning the short sequences and running phylogenetic placement

## Option 1: mafft add fragments (don't run)

mkdir option_1

# Add short sequences to the clean 12S alignment
# I needed to use 6merpair since the automatically detected best algorithm kept running into a segmentation fault
# So far, I am not enforcing keep-length since there are species represented in the short sequences that were missing from the long ones
conda activate align-env
mafft --6merpair --adjustdirection --addfragments short_12s_select.fasta --thread 24 ../align_12s/aligned_12s_nr_spec.fasta > option_1/combined_12s_aligned.fasta
conda deactivate
# Manual inspection showed that the alignment does not look good
# There were many spurious gaps


## option 2: re-align nr_spec sequences

mkdir option_2
cd option_2

# Repeat multiple sequence alignment (MSA) with all sequences: curated long sequences and short sequences
cat ../../align_12s/seqs_12s_nr_spec.fasta ../short_12s_select.fasta > combined_12s_seqs.fasta
conda activate align-env
mafft --adjustdirection combined_12s_seqs.fasta > combined_12s_aligned.fasta
conda deactivate

# Trim alignment to previously determined positions of 12S based on long sequence
# Start/end were manually determined after visual inspection of alignment in JalView
conda activate seqkit-2.8.2
seqkit subseq -r 312:1892 combined_12s_aligned.fasta | seqkit seq -w 0 -u > combined_12s_aligned_trimmed.fasta

# Split alignment into long and short sequences for tree reconstruction and placement
# Remove short sequences that do not meet length threshold of 150bp anymore based on trimmed alignment
grep '^>' ../short_12s_select.fasta | sed -e 's/^>//' -e 's/|/\\|/g' | seqkit grep -w 0 -v -r -f - combined_12s_aligned_trimmed.fasta > long_12s_aligned.fasta
grep '^>' ../short_12s_select.fasta | sed -e 's/^>//' -e 's/|/\\|/g' | seqkit grep -r -f - combined_12s_aligned_trimmed.fasta | seqkit seq -g -m 150 | grep '^>' | sed -e 's/^>//' -e 's/|/\\|/g' | seqkit grep -w 0 -r -f - combined_12s_aligned_trimmed.fasta > short_12s_aligned.fasta

# Since phylogenetic placement won't work with duplicate sequences in tree, remove redundancy between species, but carefully inspect identical sequences
# Identical seuqences across species may need to be resolved in further curation steps
seqkit rmdup -w 0 -s -D multispecies_duplicates.txt long_12s_aligned.fasta > long_12s_aligned_nr.fasta
conda deactivate

# Build non-redundant phylogenetic tree
conda activate anvio-8
FastTree -log long_12s_fasttree_nr.log -gtr -gamma -nt long_12s_aligned_nr.fasta > long_12s_fasttree_nr.tree
conda deactivate

# Extract evolutionary model parameters
conda activate orthofinder-2.5.5
raxml-ng --evaluate --msa long_12s_aligned_nr.fasta --tree long_12s_fasttree_nr.tree --prefix info --model GTR+G
conda deactivate

# Remove short sequences with less than 95% similarity to any of the long sequences per species or those with a better alignment to another species
# See: filter_alignment_sequence_similarity.R

# Each phylogenetic placement for species not already represented in the long sequence selection was manually inspected and a blacklist was compiled for accessions that should be removed
# This was an iterative process
# For each iteration, the blacklisted accessions were removed at this analysis step for the next iteration until no more anomalous placements were observed
# The manually compiled blacklist is provided in repository under assets
conda activate seqkit-2.8.2
cut -d'|' -f2 short_dissimilar.accnos | cat - $REPO/assets/placement.accnos | sed -e 's/^/\\|/' -e 's/$/\\|/' | seqkit grep -w 0 -v -r -f - short_12s_aligned.fasta > short_12s_filtered.fasta
conda deactivate
# 449 sequences remaining

# Run phylogenetic placement
conda activate paprica-env
epa-ng --tree long_12s_fasttree_nr.tree --ref-msa long_12s_aligned_nr.fasta --query short_12s_filtered.fasta --model info.raxml.bestModel
# The placement was manually curated and short sequences with anomalous placements were written to the blacklist: placement.accnos
# This manual curation focused on short sequences of species that were not represented in the long sequence selection
# This was an iterative process, repeating the short sequence filter and placement until no further anomalies were detected

# Tree graft with gappa (optional)
gappa examine graft --jplace-path epa_result.jplace
conda deactivate


## Option 3: placement with SEPP (don't run)

mkdir option_3
cd option_3

# Since phylogenetic placement won't work with duplicate sequences in tree, remove redundancy between species, but carefully inspect identical sequences
conda activate seqkit-2.8.2
seqkit rmdup -w 0 -s -D multispecies_duplicates.txt ../../align_12s/aligned_12s_nr_spec.fasta > aligned_12s_nr.fasta
conda deactivate
# Reconstruct phylogenetic tree
conda activate anvio-8
FastTree -log aligned_12s_fasttree_nr.log -gtr -gamma -nt aligned_12s_nr.fasta > aligned_12s_fasttree_nr.tree
conda deactivate
# Extract evolutionary model parameters and format output for SEPP
conda activate busco-5.8.0
export LC_ALL=C
raxml -g aligned_12s_fasttree_nr.tree -s aligned_12s_nr.fasta -m GTRGAMMA -n nr_tree -p 42
sed '/Partition: 0 with name: No Name Provided/d' RAxML_info.nr_tree > RAxML_info_fixed.nr_tree
# SEPP will generage the alignment and placement at the same time
run_sepp.py -t aligned_12s_fasttree_nr.tree -r RAxML_info_fixed.nr_tree -a aligned_12s_nr.fasta -f ../short_12s_select.fasta -o combined_12s_sepp -rt
conda deactivate
# This approach was discarded after manual inspection of the alignment


### Collect clean sequences for reference DB

cd $WDIR
mkdir finalize_12s
cd finalize_12s
cat ../placement_12s/option_2/long_12s_aligned.fasta ../placement_12s/option_2/short_12s_filtered.fasta > tmp.fasta
# 1553 sequences

# Format fasta header: >accession|source|rc(if applicable)|path(no spaces)
# Take species name from existing fasta header to ensure that valid species name is used
conda activate taxonkit-0.20.0
grep '^>' tmp.fasta | sed 's/^>//' | paste - <(grep '^>' tmp.fasta | cut -d'|' -f2 | sed 's/_R_//') | paste -d'|' - <(grep '^>' tmp.fasta | cut -d'|' -f3) | paste -d'|' - <(grep '^>' tmp.fasta | sed -e '/_R_/! s/^.*$//' -e 's/^.*_R_.*$/rc/') | paste -d'|' - <(grep '^>' tmp.fasta | sed 's/^>//' | cut -d'|' -f1 | sed -e 's/_R_//' -e 's/_/ /' | taxonkit name2taxid --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | taxonkit reformat2 -I 2 --format "k__{kingdom};p__{phylum};c__{class};o__{order};f__{family};g__{genus}" --data-dir /srv/bio/Databases/NCBI_taxdump/20260712 | cut -f3) | paste -d';' - <(grep '^>' tmp.fasta | sed 's/^>//' | cut -d'|' -f1 | sed -e 's/_R_//' -e 's/^/s__/') | sed 's/ /_/g' > tmp

# Fix empty ranks (only necessary for oder level at the moment)
sed 's/c__\([^;]*\);o__;/c__\1;o__\1_incertae_sedis;/' -i tmp

# Replace headers with new name
conda activate seqkit-2.8.2
seqkit replace -w 0 -p '^(.+)$' -r '{kv}' -k tmp tmp.fasta > reference_12s_aligned.fasta
rm tmp tmp.fasta

# Remove gaps for a non-aligned version of the reference sequence database
seqkit seq -w 0 -g reference_12s_aligned.fasta > reference_12s_seqs.fasta

# Formatting of final sequence selection for use as reference database with blast and dada2

# dada2: only retain taxonomic path in fasta header
sed -e '/^>/s/^.*|/>/' -e '/^>/s/[kpcofgs]__//g' -e '/^>/s/$/;/' reference_12s_seqs.fasta > reference_12s_dada2.fasta

# blast: only retain accession and species in fasta header, separated by '|'
sed '/^>/s/|.*;s__/|/' reference_12s_seqs.fasta > reference_12s_blast.fasta
conda activate anvio-8
makeblastdb -in reference_12s_blast.fasta -dbtype nucl -out reference_12s_blastdb -logfile reference_12s_makeblastdb.log
conda deactivate


### Useful summaries

cd $WDIR
conda activate seqkit-2.8.2

# Number of sequences in the reference sequence collection at the beginning of the workflow
wc -l sequences_*/*_baltic.accnos
#  2155 sequences_midori/midori_baltic.accnos
#  7900 sequences_ncbi/ncbi_baltic.accnos
#   339 sequences_phyloref/phyloref_baltic.accnos
# 10394 total

# Number of sequences prior to alignment (after removing redundancy based on accession numbers and meta-fish-lib blacklist accessions)
cat sequences_ncbi/ncbi_baltic_select.fasta sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | seqkit seq -m 930 | grep -c '^>'
# 3078 longer than 930bp
cat sequences_ncbi/ncbi_baltic_select.fasta sequences_phyloref/phyloref_baltic_select.fasta sequences_midori/midori_baltic_select.fasta | seqkit seq -m 600 -M 929 | grep -c '^>'
# 908 between 600 and 930bp
cat sequences_phyloref/phyloref_baltic_select_short.fasta sequences_midori/midori_baltic_select_short.fasta sequences_ncbi/ncbi_baltic_select_short.fasta | grep -c '^>'
# 3384 sequences shorter than 600bp

# Number of curated long (>=600bp) sequences
grep -c '^>' align_12s/aligned_12s_clean.fasta
# 3091 sequences

# Number of non-redundant curated long (>=600bp) sequences (dereplication only within species) 
grep -c '^>' align_12s/aligned_12s_nr_spec.fasta
# 1104 sequences --> those went into final DB

# Number of non-redundant short (<600bp) sequences (dereplication only within species) 
grep -c '^>' placement_12s/short_12s_nr_spec.fasta
# 2329 sequences

# Number of curated non-redundant short (<600bp) sequences (dereplication only within species) 
grep -c '^>' placement_12s/option_2/short_12s_filtered.fasta
# 449 sequences --> those went into final DB

# Number of species in the final DB
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort -u | wc -l
# 253 species in total
seqkit seq -m 930 finalize_12s/reference_12s_seqs.fasta | grep '^>' | cut -d';' -f7 | sort -u | wc -l
# 203 species are represented by at least one full length sequences
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ /\t/' | awk '$1 == 1' | wc -l
# 21 species are only represented by one (non-redundant) sequence
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ /\t/' | awk '$1 == 1' | cut -f2 | sed 's/^s__//' | grep -F -f - placement_12s/derep_logs_short/*.log align_12s/derep_logs/*.log | cut -f2 | cut -d'|' -f1 | sort -u | wc -l
# 10 of those species have more than one published sequence, but those are identical

# The following species are missing in the database
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort -u | sed 's/^s__//' | grep -w -v -F -f - <(cut -f4 $REPO/assets/balticfishdb_species_list.tsv | sed '1d' | sort -u | sed 's/ /_/')
# Argentina_sphyraena
# Coregonus_pallasii
# Dentex_maroccanus
# Lebetus_scorpioides
# Liparis_barbatus
# Orcynopsis_unicolor
# Oxynotus_centrina
# Platichthys_solemdali
# Sarpa_salpa
# Taractichthys_longipinnis
# Tetronarce_nobiliana

# The following species do not have a full length representative
seqkit seq -m 930 finalize_12s/reference_12s_seqs.fasta | grep '^>' | cut -d';' -f7 | sort -u | grep -w -v -F -f - <(grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort -u) | sed 's/^s__//'
# Agonus_cataphractus
# Alburnoides_bipunctatus
# Atherina_presbyter
# Babka_gymnotrachelus
# Ballerus_ballerus
# Belone_belone
# Brama_brama
# Buglossidium_luteum
# Callionymus_lyra
# Callionymus_maculatus
# Centrolophus_niger
# Cheilopogon_heterurus
# Chelidonichthys_lastoviza
# Ciliata_mustela
# Ciliata_septentrionalis
# Conger_conger
# Coris_julis
# Coryphaenoides_rupestris
# Gaidropsarus_vulgaris
# Halobatrachus_didactylus
# Hyperoplus_immaculatus
# Lampetra_fluviatilis
# Lampetra_planeri
# Lebetus_guilleti
# Leuciscus_aspius
# Leuciscus_leuciscus
# Liparis_liparis
# Liparis_montagui
# Lycodes_gracilis
# Misgurnus_fossilis
# Molva_molva
# Myxine_glutinosa
# Nerophis_lumbriciformis
# Nerophis_ophidion
# Petromyzon_marinus
# Phycis_blennoides
# Pomatoschistus_norvegicus
# Pomatoschistus_pictus
# Proterorhinus_marmoratus
# Pterycombus_brama
# Raniceps_raninus
# Sabanejewia_aurata
# Scophthalmus_rhombus
# Spondyliosoma_cantharus
# Symphodus_melops
# Thorogobius_ephippiatus
# Trisopterus_esmarkii
# Trisopterus_minutus
# Vimba_vimba
# Zeugopterus_norvegicus

# The following species do not have any intra-species variability in the database (only 1 non redundant sequence)
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ /\t/' | awk '$1 == 1' | cut -f2 | sed 's/^s__//'
# Agonus_cataphractus
# Alosa_alosa
# Argentina_silus
# Callionymus_lyra
# Chirolophis_ascanii
# Coregonus_albula
# Coregonus_nasus
# Coregonus_peled
# Diplodus_sargus
# Dipturus_nidarosiensis
# Halobatrachus_didactylus
# Lampetra_fluviatilis
# Lebetus_guilleti
# Liopsetta_glacialis
# Misgurnus_fossilis
# Molva_dypterygia
# Neogobius_fluviatilis
# Pelecus_cultratus
# Pomatoschistus_norvegicus
# Pterycombus_brama
# Trachipterus_arcticus

# The following species were only represented by 1 sequence to begin with
grep '^>' finalize_12s/reference_12s_seqs.fasta | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ /\t/' | awk '$1 == 1' | cut -f2 | sed 's/^s__//' | grep -v -w -F -f <(cat placement_12s/derep_logs_short/*.log align_12s/derep_logs/*.log | cut -f2 | cut -d'|' -f1 | sort -u)
# Agonus_cataphractus
# Chirolophis_ascanii
# Coregonus_peled
# Diplodus_sargus
# Halobatrachus_didactylus
# Lebetus_guilleti
# Liopsetta_glacialis
# Pelecus_cultratus
# Pomatoschistus_norvegicus
# Pterycombus_brama
# Trachipterus_arcticus
