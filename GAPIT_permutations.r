library(GAPIT)
library(parallel)

# ============================================================
# SETTINGS
# ============================================================

N_PERM <- 1000
SEED   <- 2026

traits <- c("E01", "E02", "E03", "Pooled")

# Requested parallel workers
N_CORES <- 64L

# Check SLURM allocation
slurm_cpus <- suppressWarnings(
  as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
)

if (!is.na(slurm_cpus) && slurm_cpus > 0) {
  
  if (N_CORES > slurm_cpus) {
    warning(
      "N_CORES exceeds SLURM allocation. ",
      "Reducing workers to ", slurm_cpus
    )
    N_CORES <- slurm_cpus
  }
  
  cat("SLURM CPUs allocated:", slurm_cpus, "\n")
  
} else {
  
  cat("SLURM_CPUS_PER_TASK not detected.\n")
  cat("Physical cores detected:",
      detectCores(logical = FALSE), "\n")
}

cat("Parallel workers:", N_CORES, "\n")
cat("Permutations per trait:", N_PERM, "\n")
cat("Total CMLM runs:", N_PERM * length(traits), "\n")


# ============================================================
# CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "Matrix",
  "survival",
  "GAPIT"
)

for (pkg in required_packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      "Required package '",
      pkg,
      "' is not available."
    )
  }
}

cat("\nPackage check passed.\n")
cat("GAPIT version:",
    as.character(packageVersion("GAPIT")), "\n")
cat("Matrix version:",
    as.character(packageVersion("Matrix")), "\n")
cat("survival version:",
    as.character(packageVersion("survival")), "\n")


# ============================================================
# GENOTYPE
# ============================================================

cat("\nReading genotype...\n")

geno <- read.table(
  "SWSF_Revised_Geno.hmp.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE,
  comment.char = "",
  quote = "",
  fill = TRUE
)

geno2 <- geno[, -c(176:178)]

# Same genotype matrix intended for GWAS
G_use <- geno2

cat(
  "Genotype dimensions:",
  nrow(G_use), "x", ncol(G_use), "\n"
)


# ============================================================
# PHENOTYPE
# ============================================================

pheno <- as.data.frame(SWSF_BLUEs)

required_columns <- c(
  "Taxa",
  "E01",
  "E02",
  "E03",
  "Pooled"
)

missing_columns <- setdiff(
  required_columns,
  names(pheno)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing phenotype columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

pheno$E01    <- as.numeric(pheno$E01)
pheno$E02    <- as.numeric(pheno$E02)
pheno$E03    <- as.numeric(pheno$E03)
pheno$Pooled <- as.numeric(pheno$Pooled)


# ============================================================
# ONE PERMUTATION
# ============================================================

run_perm <- function(
    i,
    Y,
    G_use,
    trait,
    seed
) {
  
  # ----------------------------------------------------------
  # Reproducible trait-specific seed
  # ----------------------------------------------------------
  
  trait_number <- match(
    trait,
    c("E01", "E02", "E03", "Pooled")
  )
  
  trait_offset <- trait_number * 100000L
  
  set.seed(
    seed +
      trait_offset +
      i
  )
  
  
  # ----------------------------------------------------------
  # Permute phenotype
  # ----------------------------------------------------------
  
  Yperm <- Y
  
  Yperm[[trait]] <- sample(
    Y[[trait]],
    size = nrow(Y),
    replace = FALSE
  )
  
  
  # ----------------------------------------------------------
  # GAPIT CMLM
  # ----------------------------------------------------------
  
  fit <- tryCatch(
    
    GAPIT(
      Y = Yperm,
      G = G_use,
      model = "CMLM",
      PCA.total = 3,
      file.output = FALSE
    ),
    
    error = function(e) {
      
      message(
        "Trait ", trait,
        " | Permutation ", i,
        " FAILED: ",
        conditionMessage(e)
      )
      
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Check GAPIT result
  # ----------------------------------------------------------
  
  if (
    is.null(fit) ||
    is.null(fit$GWAS)
  ) {
    return(NA_real_)
  }
  
  
  # ----------------------------------------------------------
  # Extract P-values
  # ----------------------------------------------------------
  
  if ("P.value" %in% names(fit$GWAS)) {
    
    p <- fit$GWAS$P.value
    
  } else if ("P" %in% names(fit$GWAS)) {
    
    p <- fit$GWAS$P
    
  } else {
    
    return(NA_real_)
  }
  
  p <- as.numeric(p)
  
  p <- p[
    is.finite(p) &
      p > 0 &
      p <= 1
  ]
  
  if (length(p) == 0) {
    return(NA_real_)
  }
  
  
  # ----------------------------------------------------------
  # Return genome-wide minimum P
  # ----------------------------------------------------------
  
  min(p)
}


# ============================================================
# RUN ONE TRAIT
# ============================================================

run_trait <- function(trait) {
  
  cat("\n\n")
  cat("=============================================\n")
  cat(" TRAIT:", trait, "\n")
  cat("=============================================\n")
  
  
  # ----------------------------------------------------------
  # Prepare phenotype
  # ----------------------------------------------------------
  
  Y <- pheno[
    ,
    c("Taxa", trait),
    drop = FALSE
  ]
  
  Y <- Y[
    !is.na(Y[[trait]]) &
      is.finite(Y[[trait]]),
    ,
    drop = FALSE
  ]
  
  cat("Phenotyped lines :", nrow(Y), "\n")
  cat("Permutations     :", N_PERM, "\n")
  cat("Parallel workers :", N_CORES, "\n")
  
  start_time <- Sys.time()
  
  
  # ==========================================================
  # PARALLEL PERMUTATIONS
  #
  # Linux fork-based parallelization.
  # NO PSOCK cluster / NO socket connections.
  # ==========================================================
  
  minP <- mclapply(
    X = seq_len(N_PERM),
    FUN = run_perm,
    Y = Y,
    G_use = G_use,
    trait = trait,
    seed = SEED,
    mc.cores = N_CORES,
    mc.preschedule = FALSE,
    mc.set.seed = FALSE
  )
  
  minP <- unlist(
    minP,
    use.names = FALSE
  )
  
  
  # ==========================================================
  # SAVE ALL RAW PERMUTATION RESULTS IMMEDIATELY
  # ==========================================================
  
  raw_result <- data.frame(
    Trait = trait,
    Permutation = seq_len(N_PERM),
    MinP = minP
  )
  
  raw_file <- paste0(
    "Permutation_MinP_",
    trait,
    ".csv"
  )
  
  write.csv(
    raw_result,
    raw_file,
    row.names = FALSE
  )
  
  cat("\nRaw results saved:", raw_file, "\n")
  
  
  # ==========================================================
  # CHECK SUCCESS
  # ==========================================================
  
  successful <- is.finite(minP)
  
  n_success <- sum(successful)
  n_failed  <- N_PERM - n_success
  
  cat("\nSuccessful permutations:",
      n_success, "/", N_PERM, "\n")
  
  cat("Failed permutations:",
      n_failed, "\n")
  
  
  if (n_success == 0) {
    
    stop(
      "No successful permutations for ",
      trait,
      ". Check GAPIT/dependencies."
    )
  }
  
  
  # Warn if ANY permutation failed
  if (n_failed > 0) {
    
    warning(
      trait,
      " has ",
      n_failed,
      " failed permutations."
    )
  }
  
  
  minP_success <- minP[successful]
  
  
  # ==========================================================
  # EMPIRICAL GENOME-WIDE THRESHOLDS
  # ==========================================================
  
  P_05 <- as.numeric(
    quantile(
      minP_success,
      probs = 0.05,
      names = FALSE,
      type = 7
    )
  )
  
  P_01 <- as.numeric(
    quantile(
      minP_success,
      probs = 0.01,
      names = FALSE,
      type = 7
    )
  )
  
  
  # ==========================================================
  # RUNTIME
  # ==========================================================
  
  end_time <- Sys.time()
  
  elapsed <- difftime(
    end_time,
    start_time,
    units = "mins"
  )
  
  
  # ==========================================================
  # TRAIT SUMMARY
  # ==========================================================
  
  result <- data.frame(
    
    Trait = trait,
    
    Successful_Permutations = n_success,
    
    Failed_Permutations = n_failed,
    
    P_threshold_0.05 = P_05,
    
    MinusLog10P_0.05 =
      -log10(P_05),
    
    P_threshold_0.01 = P_01,
    
    MinusLog10P_0.01 =
      -log10(P_01),
    
    Runtime_minutes =
      as.numeric(elapsed)
  )
  
  
  cat("\n---------------------------------------------\n")
  cat("RESULT FOR", trait, "\n")
  cat("---------------------------------------------\n")
  
  print(
    result,
    digits = 8
  )
  
  
  # ----------------------------------------------------------
  # Save individual trait summary
  # ----------------------------------------------------------
  
  write.csv(
    result,
    paste0(
      "Permutation_Threshold_",
      trait,
      ".csv"
    ),
    row.names = FALSE
  )
  
  
  # ----------------------------------------------------------
  # Memory cleanup between traits
  # ----------------------------------------------------------
  
  gc()
  
  result
}


# ============================================================
# RUN ALL FOUR TRAITS
# ============================================================

cat("\n\n")
cat("=============================================\n")
cat(" STARTING CMLM PERMUTATION ANALYSIS\n")
cat("=============================================\n")

cat(
  "Total CMLM analyses:",
  N_PERM * length(traits),
  "\n"
)

cat(
  "Parallel workers:",
  N_CORES,
  "\n\n"
)

all_start <- Sys.time()


results_list <- lapply(
  traits,
  run_trait
)


results <- do.call(
  rbind,
  results_list
)


# ============================================================
# SAVE FINAL THRESHOLDS
# ============================================================

write.csv(
  results,
  "GAPIT_CMLM_Permutation_Thresholds.csv",
  row.names = FALSE
)


# ============================================================
# TOTAL RUNTIME
# ============================================================

total_time <- difftime(
  Sys.time(),
  all_start,
  units = "hours"
)


# ============================================================
# FINAL OUTPUT
# ============================================================

cat("\n\n")
cat("=============================================\n")
cat(" FINAL CMLM PERMUTATION THRESHOLDS\n")
cat("=============================================\n\n")

print(
  results,
  digits = 8
)

cat(
  "\nTotal runtime:",
  round(
    as.numeric(total_time),
    2
  ),
  "hours\n"
)


cat("\nFiles saved:\n")

cat(
  "  GAPIT_CMLM_Permutation_Thresholds.csv\n"
)

for (trait in traits) {
  
  cat(
    "  Permutation_MinP_",
    trait,
    ".csv\n",
    sep = ""
  )
  
  cat(
    "  Permutation_Threshold_",
    trait,
    ".csv\n",
    sep = ""
  )
}

cat("\nAnalysis complete.\n")