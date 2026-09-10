# calculate phylogenetic resolution for mifish and tele primers

# setwd("C:/Users/chassenrueck/Documents/IOW_home/IOWseq_projects/IOWseq000085_KOFI/ReferenceDB_12S_mifish/")
require(seqinr)
require(stringr)
require(reshape)
require(tidyverse)
require(scales)
load("barcode_gap.Rdata")

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

# set gene name, start/end position in alignment
target <- "mifish"
ali_start <- 306
ali_end <- 590

target <- "tele"
ali_start <- 310
ali_end <- 590

# read mifish 12S alignment
ali_target <- read.alignment(
  paste0("insilico_pcr/", target, "_aligned.fasta"),
  format = "fasta"
)

# read full alignment
ali_12s <- read.alignment(
  "finalize_12s/reference_12s_aligned.fasta",
  format = "fasta"
)
all.equal(ali_target$nam, ali_12s$nam)

# remove sequences that do not fully cover amplified region
ali_filt <- subset_msa(
  ali_target, 
  which(
    sapply(ali_12s$seq, function(x) any(which(s2c(x) != "-") <= ali_start) & any(which(s2c(x) != "-") >= ali_end)) &
      sapply(ali_target$seq, function(x) sum(s2c(x) != "-")) >= 150
  )
)
ali_lengths <- sapply(ali_filt$seq, function(x) sum(s2c(x) != "-"))
names(ali_lengths) <- ali_filt$nam

# get species list from sequence names
spec_names <- table(sort(gsub("^s__","", str_split_i(ali_filt$nam, ";", 7))))

# calculate within species similarity
within_sp <- map_dfr(
  names(spec_names),
  function(x) {
    msa_sub <- subset_msa(ali_filt, grep(paste0(";s__", x), ali_filt$nam, fixed = T))
    msa_mat <- dist.alignment(
      msa_sub, 
      matrix = "identity",
      gap = T
    ) %>% 
      as.matrix()
    msa_mat[upper.tri(msa_mat)] <- NA
    msa_df <- msa_mat %>% 
      melt() %>% 
      filter(!is.na(value), X1 != X2) %>% 
      mutate(species = x, comparison = "within", sim = 1 - value^2) %>% 
      select(-value)
    return(msa_df)
  }
)

# calculate between species similarity
# this takes quite a while... 
between_sp <- map_dfr(
  names(spec_names),
  function(x) {
    # subset to sequences per species
    msa_sub_sp <- subset_msa(ali_filt, grep(paste0(";s__", x), ali_filt$nam, fixed = T))
    msa_sub_other <- subset_msa(ali_filt, grep(paste0(";s__", x), ali_filt$nam, invert = T, fixed = T))
    
    msa_df <- map_dfr(
      1:msa_sub_sp$nb,
      function(i) {
        map_dfr(
          1:msa_sub_other$nb,
          function(j) {
            msa_tmp <- as.alignment(
              nb = 2, 
              nam = c(msa_sub_sp$nam[i], msa_sub_other$nam[j]), 
              seq = c(msa_sub_sp$seq[i], msa_sub_other$seq[j])
            )
            data.frame(
              X1 = msa_tmp$nam[1],
              X2 = msa_tmp$nam[2],
              species = x,
              comparison = "between",
              sim = 1 - melt(dist.alignment(msa_tmp, matrix = "identity", gap = T))[1, 1]^2
            )
          }
        )
      }
    )
    return(msa_df)
  }
)
write.table(
  bind_rows(within_sp, between_sp), 
  paste0("insilico_pcr/", target, "_barcode_gap_similarity.txt"),
  quote = F, 
  sep = "\t",
  row.names = F
)

# overall plot
pdf(paste0("insilico_pcr/", target, "_barcode_gap_overview.pdf"), height = 7, width = 7)
p1 <- hist(within_sp$sim, breaks = seq(0, 1, 0.01))
p2 <- hist(between_sp$sim, breaks = seq(0, 1, 0.01))
p1$density <- rescale(p1$density, to = c(0, 1))
p2$density <- rescale(p2$density, to = c(0, 1))
plot(p2, col = adjustcolor("red", alpha.f = 0.2), freq = F, axes = F, xlab = "", ylab = "", main = "Barcode gap")
plot(p1, col = adjustcolor("blue", alpha.f = 0.2), freq = F, add = T)
axis(1)
dev.off()

# per species
pdf(paste0("insilico_pcr/", target, "_barcode_gap_species.pdf"), width = 7, height = 7, onefile = T)
for(i in names(spec_names)) {
  if(sum(within_sp$species == i) == 0) {
    p2 <- hist(between_sp$sim[between_sp$species == i], breaks = seq(0, 1, 0.01), plot = F)
    p2$density <- rescale(p2$density, to = c(0, 1))
    plot(p2, col = adjustcolor("red", alpha.f = 0.2), freq = F, axes = F, xlab = "", ylab = "", main = i)
    abline(v = 1, col = "blue")
    axis(1)
  } else {
    p1 <- hist(within_sp$sim[within_sp$species == i], breaks = seq(0, 1, 0.01), plot = F)
    p2 <- hist(between_sp$sim[between_sp$species == i], breaks = seq(0, 1, 0.01), plot = F)
    p1$density <- rescale(p1$density, to = c(0, 1))
    p2$density <- rescale(p2$density, to = c(0, 1))
    plot(p2, col = adjustcolor("red", alpha.f = 0.2), freq = F, axes = F, xlab = "", ylab = "", main = i)
    plot(p1, col = adjustcolor("blue", alpha.f = 0.2), freq = F, add = T)
    axis(1)
  }
}
dev.off()

# calculate ranges between and within species
sim_ranges <- bind_rows(within_sp, between_sp) %>% 
  group_by(species) %>% 
  summarise(
    min_sp = ifelse(sum(comparison == "within") == 0, 1, min(sim[comparison == "within"])),
    max_sp = ifelse(sum(comparison == "within") == 0, 1, max(sim[comparison == "within"])),
    min_other = min(sim[comparison == "between"]),
    max_other = max(sim[comparison == "between"])
  )
write.table(sim_ranges, paste0("insilico_pcr/", target, "_sim_ranges.txt"), row.names = F, sep = "\t", quote = F)

# save.image(paste0("insilico_pcr/", target, "_barcode_gap.Rdata"))
            
