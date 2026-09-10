# prepare working environment
setwd("C:/Users/chassenrueck/Documents/IOW_home/IOWseq_projects/IOWseq000085_KOFI/ReferenceDB_12S_mifish/")
require(seqinr)
require(stringr)
require(reshape)
require(tidyverse)

# function to subset msa
subset_msa <- function(a, n) {
  nb = length(n)
  if(is.character(n)) {
    nam = n
    seq = a$seq[match(n, a$nam)]
  } else {
    nam = a$nam[n]
    seq = a$seq[n]
  }
  as.alignment(nb, nam, seq)
}


### calculate pairwise sequence similarity from 12S MSA within each species ####

# read 12S alignment
msa_12s <- read.alignment(
  "aligned_12s_trimmed.fasta",
  format = "fasta"
)

# get species list from sequence names
spec_names <- table(sort(str_split_i(msa_12s$nam, "\\|", 2)))

# calculate pairwise distance and distance to mediod for each species
spec_sequence_similarity <- map_dfr(
  names(spec_names),
  function(x) {
    if(spec_names[x] == 1) {
      out <- data.frame(
        species = x,
        nseq = 1,
        sequence = grep(paste0("|", x, "|"), msa_12s$nam, fixed = T, value = T),
        dist2mediod = 0,
        sim2mediod = 1
      )
    } else {
      # subset to sequences per species
      msa_sub <- subset_msa(msa_12s, grep(paste0("|", x, "|"), msa_12s$nam, fixed = T))
      # calculate pairwise similarity
      msa_dist <- dist.alignment(
        msa_sub, 
        matrix = "identity"
      )
      # calculate distance to mediod 
      # mediod defined as sequences with lowest cumulative distance to all other sequences
      msa_dist_mat <- as.matrix(as.dist(msa_dist))
      msa_dist_medoid <- msa_dist_mat[, which.min(rowSums(msa_dist_mat))] 
      # format output
      out <- data.frame(
        species = x,
        nseq = length(msa_dist_medoid),
        sequence = names(msa_dist_medoid),
        dist2mediod = msa_dist_medoid,
        sim2mediod = 1 - msa_dist_medoid^2
      )
    }
    return(out)
  }
)

# filter by sequence similarity
# remove sequences less than 97% similar to mediod
# but only if more than 2 sequences per species
spec_sequence_similarity <- spec_sequence_similarity %>% 
  mutate(
    discard = ifelse(nseq > 2, sim2mediod < 0.97, FALSE)
  )

# loss of sequences per species if threshold is applied
spec_sequence_loss <- spec_sequence_similarity %>% 
  group_by(species) %>% 
  summarise(nseq_in = n(), nseq_out = sum(!discard)) %>% 
  mutate(loss = (1 - nseq_out/nseq_in) * 100)

# write output
write.table(spec_sequence_similarity, "aligned_similarity_check.txt", row.names = F, sep = "\t", quote = F)
write.table(spec_sequence_loss, "aligned_12s_discard_loss.txt", row.names = F, sep = "\t", quote = F)


### calculate sequence similarity between short and long reference sequences within species ####

blastout <- read.table(
  "short_12s_nr_spec.blastout",
  sep = "\t",
  h = F,
  col.names = c(
    "qseqid", 
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
    "qcovs",
    "qlen",
    "slen",
    "qseq"
  )
) %>% 
  mutate(
    sspec = gsub("\\|.*$", "", sseqid),
    qspec = gsub("\\|.*$", "", qseqid),
    qcov_corrected = (qend - qstart + 1) / qlen
  )
# min(blastout$qend - blastout$qstart) > 0: so all seqs are aligned in fwd orientation (required assumption)

# identify short sequences that are identical to long sequences of the same species
short_identical <- blastout %>% 
  filter(
    mismatch == 0 & gapopen == 0 & qstart == 1 & qend == qlen & sspec == qspec
  ) %>% 
  pull(qseqid) %>% 
  unique()
write(short_identical, "short_identical.accnos")

# identify short sequences that are highly divergent from long sequences of the same species
# most of the query needs to be aligned
short_divergent <- blastout %>% 
  filter(
    qcov_corrected >= 0.98 & pident < 90 & sspec == qspec
  ) %>% 
  pull(qseqid) %>% 
  unique()
write(short_divergent, "short_divergent.accnos")


### Species consistency based on MSA for short sequences ####

# read long alignment
msa_12s_long <- read.alignment(
  "option_2/long_12s_aligned.fasta",
  format = "fasta"
)

# read short alignment
msa_12s_short <- read.alignment(
  "option_2/short_12s_aligned.fasta",
  format = "fasta"
)

# extract species names
spec_names_short <- table(sort(gsub("^_R_", "", str_split_i(msa_12s_short$nam, "\\|", 1))))
spec_names_long <- table(sort(gsub("^_R_", "", str_split_i(msa_12s_long$nam, "\\|", 1))))

# for each species represented by short sequences,
# if this species is also represented by long sequences,
# calulate sequence similarity between each short and long sequences,
# record the maximum similarity.
# also check if the short sequence is more similar to a long sequence from another species
short_long_sim <- map_dfr(
  names(spec_names_short),
  function(x) {
    if(x %in% names(spec_names_long)) {
      # subset to sequences per species
      msa_sub_short <- subset_msa(msa_12s_short, grep(paste0(x, "|"), msa_12s_short$nam, fixed = T))
      msa_sub_long <- subset_msa(msa_12s_long, grep(paste0(x, "|"), msa_12s_long$nam, fixed = T))
      msa_sub_other <- subset_msa(msa_12s_long, grep(paste0(x, "|"), msa_12s_long$nam, invert = T, fixed = T))
      
      out_df <- data.frame(matrix(NA, ncol = 6, nrow = msa_sub_short$nb))
      colnames(out_df) <- c("species", "acc", "max_sim_sp", "ali_len_sp", "max_sim_other", "ali_len_other")
      out_df$species <- x
      
      # calculate pairwise similarity
      for(i in 1:msa_sub_short$nb) {
        print(paste(x, i))
        out_df$acc[i] <- msa_sub_short$nam[i]
        
        tmp_dist <- sapply(
          1:msa_sub_long$nb,
          function(j) {
            msa_tmp <- as.alignment(
              nb = 2, 
              nam = c(msa_sub_short$nam[i], msa_sub_long$nam[j]), 
              seq = c(msa_sub_short$seq[i], msa_sub_long$seq[j])
            )
            c(
              melt(dist.alignment(msa_tmp, matrix = "identity"))[1, 1],
              sum(s2c(msa_tmp$seq[[1]]) != "-" & s2c(msa_tmp$seq[[2]]) != "-")
            )
          }
        ) %>% t() %>% as.data.frame()
        tmp_dist$V1[is.na(tmp_dist$V1)] <- 1
        out_df$max_sim_sp[i] <- max(1 - tmp_dist$V1^2)
        out_df$ali_len_sp[i] <- tmp_dist$V2[which.max(1 - tmp_dist$V1^2)]
        
        tmp_dist <- sapply(
          1:msa_sub_other$nb,
          function(j) {
            msa_tmp <- as.alignment(
              nb = 2, 
              nam = c(msa_sub_short$nam[i], msa_sub_other$nam[j]), 
              seq = c(msa_sub_short$seq[i], msa_sub_other$seq[j])
            )
            c(
              melt(dist.alignment(msa_tmp, matrix = "identity"))[1, 1],
              sum(s2c(msa_tmp$seq[[1]]) != "-" & s2c(msa_tmp$seq[[2]]) != "-")
            )
          }
        ) %>% t() %>% as.data.frame()
        tmp_dist$V1[is.na(tmp_dist$V1)] <- 1
        out_df$max_sim_other[i] <- max(1 - tmp_dist$V1^2)
        out_df$ali_len_other[i] <- tmp_dist$V2[which.max(1 - tmp_dist$V1^2)]
      }
      return(out_df)
    }
  }
)

hist(short_long_sim$max_sim_sp, breaks = 100, xlim = c(0.8, 1))
abline(v = 0.95)

# remove short sequences with less than 95% similarity
# as well as those who have a better alignment (both length and similarity to another species)
short_dissimilar <- short_long_sim %>% 
  filter(max_sim_sp < 0.95 | (ali_len_other > ali_len_sp & max_sim_other > max_sim_sp)) %>% 
  pull(acc) %>% 
  unique()
write(short_dissimilar, "short_dissimilar.accnos")

