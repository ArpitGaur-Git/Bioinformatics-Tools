library(vcfR)
library(poppr)
library(ape)

VCF_FILE <- "SWSF_Imputed.vcf"

OUT_DIR <- "GBS_NJ_Phylogeny"

N_BOOT <- 1000
SEED <- 2026

# Maximum allowed missing genotype proportion per SNP
MAX_SNP_MISSING <- 0.20

# Maximum allowed missing genotype proportion per sample
MAX_SAMPLE_MISSING <- 0.20

# Minimum minor allele frequency
MIN_MAF <- 0.05

# Minimum number of shared SNPs for a genotype pair
MIN_SHARED_SNP <- 100

# Display bootstrap values >= this threshold
MIN_BOOT_DISPLAY <- 50

set.seed(SEED)

dir.create(
  OUT_DIR,
  showWarnings = FALSE,
  recursive = TRUE
)


# ================================================================
# 3. READ VCF
# ================================================================

cat("\nReading VCF...\n")

vcf <- read.vcfR(
  VCF_FILE,
  verbose = FALSE
)

gt <- extract.gt(
  vcf,
  element = "GT",
  as.numeric = FALSE
)

cat("Original SNPs:", nrow(gt), "\n")
cat("Original genotypes:", ncol(gt), "\n")


# ================================================================
# 4. RETAIN BIALLELIC SNPs
# ================================================================

ref <- vcf@fix[, "REF"]
alt <- vcf@fix[, "ALT"]

valid_snp <- (
  grepl("^[ACGT]$", ref) &
    grepl("^[ACGT]$", alt)
)

gt <- gt[valid_snp, , drop = FALSE]

cat("Biallelic SNPs:", nrow(gt), "\n")


# ================================================================
# 5. CONVERT VCF GENOTYPES TO DOSAGE
# ================================================================

# 0/0 = 0
# 0/1 = 1
# 1/1 = 2
# Missing = NA
#
# Assumes diploid-style, biallelic genotype calls.

gt <- gsub("\\|", "/", gt)

geno <- matrix(
  NA_real_,
  nrow = nrow(gt),
  ncol = ncol(gt),
  dimnames = dimnames(gt)
)

geno[gt == "0/0"] <- 0
geno[gt %in% c("0/1", "1/0")] <- 1
geno[gt == "1/1"] <- 2

# Genotypes x SNPs
geno <- t(geno)

rm(gt, vcf)
gc()

cat("Dosage matrix created.\n")


# ================================================================
# 6. SAMPLE MISSINGNESS FILTER
# ================================================================

sample_missing <- rowMeans(is.na(geno))

keep_samples <- sample_missing <= MAX_SAMPLE_MISSING

geno <- geno[keep_samples, , drop = FALSE]

cat(
  "Samples retained:",
  nrow(geno),
  "\n"
)

if (nrow(geno) < 4) {
  stop("Fewer than four samples retained.")
}


# ================================================================
# 7. SNP MISSINGNESS FILTER
# ================================================================

snp_missing <- colMeans(is.na(geno))

keep_snps <- snp_missing <= MAX_SNP_MISSING

geno <- geno[, keep_snps, drop = FALSE]

cat(
  "SNPs after missingness filtering:",
  ncol(geno),
  "\n"
)


# ================================================================
# 8. MINOR ALLELE FREQUENCY FILTER
# ================================================================

allele_frequency <- colMeans(
  geno,
  na.rm = TRUE
) / 2

maf <- pmin(
  allele_frequency,
  1 - allele_frequency
)

keep_maf <- is.finite(maf) & maf >= MIN_MAF

geno <- geno[, keep_maf, drop = FALSE]

cat(
  "SNPs after MAF filtering:",
  ncol(geno),
  "\n"
)

if (ncol(geno) < MIN_SHARED_SNP) {
  stop("Insufficient SNPs after filtering.")
}


# ================================================================
# 9. GENETIC DISTANCE FUNCTION
# ================================================================

# Mean absolute allelic dosage difference
#
# 0/0 vs 0/0 = 0
# 0/0 vs 0/1 = 0.5
# 0/0 vs 1/1 = 1
# 0/1 vs 1/1 = 0.5
#
# Missing SNPs are excluded pairwise.

genetic_distance <- function(x) {
  
  n <- nrow(x)
  
  d <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(
      rownames(x),
      rownames(x)
    )
  )
  
  for (i in seq_len(n - 1)) {
    
    for (j in seq.int(i + 1, n)) {
      
      valid <- (
        !is.na(x[i, ]) &
          !is.na(x[j, ])
      )
      
      n_shared <- sum(valid)
      
      if (n_shared < MIN_SHARED_SNP) {
        
        stop(
          paste(
            "Insufficient shared SNPs:",
            rownames(x)[i],
            rownames(x)[j],
            n_shared
          )
        )
      }
      
      dij <- mean(
        abs(
          x[i, valid] -
            x[j, valid]
        ) / 2
      )
      
      d[i, j] <- dij
      d[j, i] <- dij
    }
  }
  
  as.dist(d)
}


# ================================================================
# 10. CALCULATE GENETIC DISTANCE
# ================================================================

cat("\nCalculating genetic distances...\n")

dist_matrix <- genetic_distance(geno)

write.csv(
  as.matrix(dist_matrix),
  file.path(
    OUT_DIR,
    "GBS_Genetic_Distance.csv"
  )
)


# ================================================================
# 11. CONSTRUCT NJ TREE
# ================================================================

cat("\nConstructing NJ tree...\n")

nj_tree <- nj(dist_matrix)

# NJ is naturally unrooted
nj_tree <- unroot(nj_tree)

write.tree(
  nj_tree,
  file = file.path(
    OUT_DIR,
    "GBS_NJ_Original.nwk"
  )
)


# ================================================================
# 12. BOOTSTRAP FUNCTION
# ================================================================

bootstrap_nj <- function(x) {
  
  sampled_snps <- sample(
    seq_len(ncol(x)),
    size = ncol(x),
    replace = TRUE
  )
  
  boot_geno <- x[
    ,
    sampled_snps,
    drop = FALSE
  ]
  
  boot_dist <- genetic_distance(boot_geno)
  
  boot_tree <- nj(boot_dist)
  
  unroot(boot_tree)
}


# ================================================================
# 13. RUN 1000 BOOTSTRAP REPLICATES
# ================================================================

cat("\nStarting bootstrap analysis...\n")

boot_trees <- vector(
  "list",
  N_BOOT
)

start_time <- Sys.time()

for (b in seq_len(N_BOOT)) {
  
  boot_trees[[b]] <- bootstrap_nj(geno)
  
  if (b %% 50 == 0) {
    
    elapsed <- difftime(
      Sys.time(),
      start_time,
      units = "mins"
    )
    
    cat(
      "Bootstrap:",
      b,
      "/",
      N_BOOT,
      "| Elapsed:",
      round(as.numeric(elapsed), 2),
      "minutes\n"
    )
  }
}

class(boot_trees) <- "multiPhylo"


# ================================================================
# 14. CALCULATE BOOTSTRAP SUPPORT
# ================================================================

cat("\nCalculating bootstrap support...\n")

# Unrooted split-based bootstrap support
bp <- prop.clades(
  nj_tree,
  boot_trees,
  rooted = FALSE
)

bootstrap_percentage <- round(
  100 * bp / N_BOOT,
  0
)

nj_tree$node.label <- bootstrap_percentage


# ================================================================
# 15. EXPORT NEWICK FOR iTOL
# ================================================================

cat("\nExporting iTOL tree...\n")

# All bootstrap labels
write.tree(
  nj_tree,
  file = file.path(
    OUT_DIR,
    "GBS_NJ_1000Bootstrap_iTOL.nwk"
  )
)

# Bootstrap values below threshold hidden
nj_tree_filtered <- nj_tree

nj_tree_filtered$node.label[
  bootstrap_percentage < MIN_BOOT_DISPLAY
] <- ""

write.tree(
  nj_tree_filtered,
  file = file.path(
    OUT_DIR,
    "GBS_NJ_Bootstrap50_iTOL.nwk"
  )
)


# ================================================================
# 16. SAVE ALL BOOTSTRAP TREES
# ================================================================

write.tree(
  boot_trees,
  file = file.path(
    OUT_DIR,
    "GBS_1000_Bootstrap_Trees.nwk"
  )
)


# ================================================================
# 17. EXPORT PDF PHYLOGRAM
# ================================================================

pdf(
  file.path(
    OUT_DIR,
    "GBS_NJ_Bootstrap_Phylogram.pdf"
  ),
  width = 14,
  height = 12
)

plot(
  nj_tree_filtered,
  type = "phylogram",
  direction = "rightwards",
  cex = 0.65,
  edge.width = 1.2,
  main = "GBS Neighbor-Joining Phylogram"
)

nodelabels(
  text = nj_tree_filtered$node.label,
  frame = "none",
  cex = 0.6,
  col = "blue"
)

add.scale.bar()

dev.off()


# ================================================================
# 18. EXPORT UNROOTED TREE
# ================================================================

pdf(
  file.path(
    OUT_DIR,
    "GBS_NJ_Unrooted.pdf"
  ),
  width = 14,
  height = 14
)

plot(
  nj_tree_filtered,
  type = "unrooted",
  cex = 0.65,
  edge.width = 1.2,
  main = "GBS NJ Tree - 1000 Bootstraps"
)

nodelabels(
  text = nj_tree_filtered$node.label,
  frame = "none",
  cex = 0.6,
  col = "blue"
)

dev.off()


# ================================================================
# 19. SAVE BOOTSTRAP SUPPORT TABLE
# ================================================================

bootstrap_table <- data.frame(
  Node = seq_along(bootstrap_percentage),
  Bootstrap_Percentage = bootstrap_percentage
)

write.csv(
  bootstrap_table,
  file.path(
    OUT_DIR,
    "GBS_Bootstrap_Support.csv"
  ),
  row.names = FALSE
)


# ================================================================
# 20. FINAL SUMMARY
# ================================================================

cat("\n============================================\n")
cat("GBS NJ PHYLOGENETIC ANALYSIS COMPLETED\n")
cat("============================================\n")

cat("Genotypes:", nrow(geno), "\n")
cat("SNPs:", ncol(geno), "\n")
cat("Bootstrap replicates:", N_BOOT, "\n")

cat(
  "Total runtime:",
  round(
    as.numeric(
      difftime(
        Sys.time(),
        start_time,
        units = "mins"
      )
    ),
    2
  ),
  "minutes\n"
)

cat("\nOutput directory:", OUT_DIR, "\n")

cat("\nMAIN iTOL FILE:\n")
cat("GBS_NJ_1000Bootstrap_iTOL.nwk\n")

cat("\nAnalysis completed successfully.\n")