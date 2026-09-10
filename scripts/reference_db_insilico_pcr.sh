### in silico PCR
# Configure working environment
WDIR="/srv/bio/Analysis_data/IOWseq000085_KOFI/Intermediate_results/ReferenceDB_12S_mifish"
mkdir insilico_pcr
cd insilico_pcr


### mifish-U

# FWD:GTCGGTAAAACTCGTGCCAGC
# REV:CATAGTGGGGTATCTAATCCCAGTTTG (rc: CAAACTGGGATTAGATACCCCACTATG)

conda activate qiime2-amplicon-2024.5
cutadapt -e 0.2 -O 10 -g "GTCGGTAAAACTCGTGCCAGC...CAAACTGGGATTAGATACCCCACTATG" --discard-untrimmed -o mifish.fasta ../finalize_12s/reference_12s_seqs.fasta
conda deactivate

# overview about how many sequences per species survive in PCR
echo -e "species\treference\tpcr" > mifish_pcr_out.txt
join -t $'\t' -a 1 -1 2 -2 2 <(grep '^>' ../finalize_12s/reference_12s_seqs.fasta | cut -d'|' -f4 | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ s__/\t/') <(grep '^>' mifish.fasta | cut -d'|' -f4 | cut -d';' -f7 | sed 's/ .*$//' | sort | uniq -c | sed -E 's/^ +//' | sed 's/ s__/\t/') >> mifish_pcr_out.txt
# 227 species out of 253 caught by primers

# manual check showed that in alignment, mifish is from 306:590 (no primers)
conda activate seqkit-2.8.2
seqkit subseq -r 306:590 ../finalize_12s/reference_12s_aligned.fasta | seqkit seq -w 0 > mifish_aligned.fasta
conda deactivate

### tele02

# FWD: AAACTCGTGCCAGCCACC
# REV: GGGTATCTAATCCCAGTTTG (rc: CAAACTGGGATTAGATACCC)

conda activate qiime2-amplicon-2024.5
cutadapt -e 0.2 -O 10 -g "AAACTCGTGCCAGCCACC...CAAACTGGGATTAGATACCC" --discard-untrimmed -o tele.fasta ../finalize_12s/reference_12s_seqs.fasta
conda deactivate

# overview about how many sequences per species survive in PCR
echo -e "species\treference\tpcr" > tele_pcr_out.txt
join -t $'\t' -a 1 -1 2 -2 2 <(grep '^>' ../finalize_12s/reference_12s_seqs.fasta | cut -d'|' -f4 | cut -d';' -f7 | sort | uniq -c | sed -E 's/^ +//' | sed 's/ s__/\t/') <(grep '^>' tele.fasta | cut -d'|' -f4 | cut -d';' -f7 | sed 's/ .*$//' | sort | uniq -c | sed -E 's/^ +//' | sed 's/ s__/\t/') >> tele_pcr_out.txt
# 228 species out of 253 caught by primers

# based on manually confirmed mifish coordinates, tele is 310:590 (no primers)
conda activate seqkit-2.8.2
seqkit subseq -r 310:590 ../finalize_12s/reference_12s_aligned.fasta | seqkit seq -w 0 > tele_aligned.fasta
conda deactivate

