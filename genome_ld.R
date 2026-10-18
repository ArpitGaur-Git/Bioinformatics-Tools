# ============================================================
# TASSEL 5 — WHEAT LD DECAY ANALYSIS
# Genome-wide | Subgenome-wise | Chromosome-wise
#
# OUTPUT:
# 1. Genome + Subgenome LD decay.pdf
# 2. Chromosome-wise LD decay.pdf
# 3. LD decay summary.csv
# 4. Binned LD data.csv
# 5. Fitted LD curves.csv
#
# ============================================================

# ------------------------------------------------------------
# 1. PACKAGES
# ------------------------------------------------------------

library(data.table)
library(ggplot2)

# ------------------------------------------------------------
# 2. PARAMETERS
# ------------------------------------------------------------

threshold <- 0.20
loess_span <- 0.30
bin_size_kb <- 25
max_distance_mb <- 10

output_dir <- "LD_Decay_Results"

dir.create(
  output_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ------------------------------------------------------------
# 3. LOAD TASSEL 5 TEXT FILE
# ------------------------------------------------------------

ld <- fread(
  file.choose(),
  sep = "\t",
  header = TRUE,
  na.strings = c(
    "NA", "NaN", "nan",
    "NULL", "", "-999", "-999.0"
  ),
  showProgress = TRUE
)

cat("Original SNP pairs:", nrow(ld), "\n")

# ------------------------------------------------------------
# 4. VALIDATE COLUMNS
# ------------------------------------------------------------

required <- c(
  "Locus1", "Position1",
  "Locus2", "Position2",
  "R^2"
)

missing <- setdiff(required, names(ld))

if (length(missing) > 0) {
  stop(
    "Missing columns: ",
    paste(missing, collapse = ", ")
  )
}

# ------------------------------------------------------------
# 5. PREPROCESSING
# ------------------------------------------------------------

chr_levels <- as.vector(
  t(outer(1:7, c("A", "B", "D"), paste0))
)

# Standardize chromosome names
ld[, Chr1 := toupper(
  gsub("^chr", "", trimws(as.character(Locus1)),
       ignore.case = TRUE)
)]

ld[, Chr2 := toupper(
  gsub("^chr", "", trimws(as.character(Locus2)),
       ignore.case = TRUE)
)]

# Retain only intrachromosomal comparisons
ld <- ld[
  Chr1 == Chr2 &
    Chr1 %in% chr_levels
]

cat("Intrachromosomal pairs:", nrow(ld), "\n")

# Convert numeric columns
ld[, Position1 := suppressWarnings(as.numeric(Position1))]
ld[, Position2 := suppressWarnings(as.numeric(Position2))]
ld[, R2 := suppressWarnings(as.numeric(`R^2`))]

# Remove missing and invalid values
ld <- ld[
  is.finite(Position1) &
    is.finite(Position2) &
    is.finite(R2) &
    Position1 > 0 &
    Position2 > 0 &
    R2 >= 0 &
    R2 <= 1
]

# Physical distance
ld[, Dist_bp := abs(Position2 - Position1)]
ld[, Dist_mb := Dist_bp / 1e6]
ld[, Dist_kb := Dist_bp / 1000]

# Filter distance
ld <- ld[
  Dist_bp > 0 &
    Dist_mb <= max_distance_mb
]

# Remove duplicate pairs using unique site indices if available
if (all(c("Site1", "Site2") %in% names(ld))) {
  
  ld[, SNP_A := pmin(as.character(Site1),
                     as.character(Site2))]
  
  ld[, SNP_B := pmax(as.character(Site1),
                     as.character(Site2))]
  
  ld <- unique(
    ld,
    by = c("Chr1", "SNP_A", "SNP_B")
  )
  
} else {
  
  ld[, Pos_A := pmin(Position1, Position2)]
  ld[, Pos_B := pmax(Position1, Position2)]
  
  ld <- unique(
    ld,
    by = c("Chr1", "Pos_A", "Pos_B")
  )
}

# Assign chromosome and subgenome
ld[, Chromosome := Chr1]
ld[, Subgenome := sub(".*([ABD])$", "\\1", Chr1)]

cat("Final SNP pairs:", nrow(ld), "\n")

if (nrow(ld) == 0) {
  stop("No valid LD pairs remain after filtering.")
}

# ------------------------------------------------------------
# 6. DISTANCE BINNING
# ------------------------------------------------------------

# Bin index
ld[, Bin := floor(Dist_kb / bin_size_kb)]

# Genome-wide
genome_bin <- ld[
  ,
  .(
    Dist_mb = mean(Dist_mb),
    Mean_R2 = mean(R2),
    N_pairs = .N
  ),
  by = Bin
]

genome_bin[, Group := "Genome-wide"]
genome_bin[, Level := "Genome"]

# Subgenome-wise
sub_bin <- ld[
  ,
  .(
    Dist_mb = mean(Dist_mb),
    Mean_R2 = mean(R2),
    N_pairs = .N
  ),
  by = .(Subgenome, Bin)
]

sub_bin[, Group := Subgenome]
sub_bin[, Level := "Subgenome"]

# Chromosome-wise
chr_bin <- ld[
  ,
  .(
    Dist_mb = mean(Dist_mb),
    Mean_R2 = mean(R2),
    N_pairs = .N
  ),
  by = .(Chromosome, Bin)
]

chr_bin[, Group := Chromosome]
chr_bin[, Level := "Chromosome"]

# Combine
all_bin <- rbindlist(
  list(
    genome_bin[, .(Level, Group, Bin,
                   Dist_mb, Mean_R2, N_pairs)],
    sub_bin[, .(Level, Group, Bin,
                Dist_mb, Mean_R2, N_pairs)],
    chr_bin[, .(Level, Group, Bin,
                Dist_mb, Mean_R2, N_pairs)]
  )
)

fwrite(
  all_bin,
  file.path(output_dir, "Binned_LD_Data.csv")
)

# Free memory
rm(ld, genome_bin, sub_bin, chr_bin)
gc()

# ------------------------------------------------------------
# 7. LOESS AND LD DECAY FUNCTION
# ------------------------------------------------------------

fit_ld <- function(df) {
  
  setorder(df, Dist_mb)
  
  result <- list(
    decay_mb = NA_real_,
    decay_kb = NA_real_,
    predictions = NULL,
    status = "Insufficient data"
  )
  
  if (
    nrow(df) < 10 ||
    uniqueN(df$Dist_mb) < 10
  ) {
    return(result)
  }
  
  fit <- tryCatch(
    loess(
      Mean_R2 ~ Dist_mb,
      data = df,
      span = loess_span,
      degree = 2,
      control = loess.control(surface = "direct")
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    result$status <- "LOESS failed"
    return(result)
  }
  
  pred <- data.table(
    Dist_mb = seq(
      min(df$Dist_mb),
      max(df$Dist_mb),
      length.out = 2000
    )
  )
  
  pred[, R2_fit := as.numeric(
    predict(fit, newdata = pred)
  )]
  
  pred <- pred[is.finite(R2_fit)]
  
  result$predictions <- pred
  
  if (nrow(pred) < 2) {
    result$status <- "No valid predictions"
    return(result)
  }
  
  # First downward crossing
  idx <- which(
    head(pred$R2_fit, -1) >= threshold &
      tail(pred$R2_fit, -1) < threshold
  )
  
  if (length(idx) == 0) {
    result$status <- "No threshold crossing"
    return(result)
  }
  
  i <- idx[1]
  
  x1 <- pred$Dist_mb[i]
  x2 <- pred$Dist_mb[i + 1]
  
  y1 <- pred$R2_fit[i]
  y2 <- pred$R2_fit[i + 1]
  
  decay_mb <- x1 +
    (threshold - y1) *
    (x2 - x1) / (y2 - y1)
  
  result$decay_mb <- decay_mb
  result$decay_kb <- decay_mb * 1000
  result$status <- "Estimated"
  
  return(result)
}

# ------------------------------------------------------------
# 8. FIT ALL GROUPS
# ------------------------------------------------------------

groups <- unique(
  all_bin[, .(Level, Group)]
)

summary_list <- list()
prediction_list <- list()

for (i in seq_len(nrow(groups))) {
  
  level <- groups$Level[i]
  group <- groups$Group[i]
  
  df <- copy(
    all_bin[
      Level == level & Group == group
    ]
  )
  
  result <- fit_ld(df)
  
  summary_list[[i]] <- data.table(
    Level = level,
    Group = group,
    Decay_Mb = result$decay_mb,
    Decay_kb = result$decay_kb,
    N_pairs = sum(df$N_pairs),
    Status = result$status
  )
  
  if (!is.null(result$predictions)) {
    
    prediction_list[[i]] <- copy(
      result$predictions
    )[
      ,
      `:=`(
        Level = level,
        Group = group
      )
    ]
  }
}

ld_summary <- rbindlist(summary_list)
ld_pred <- rbindlist(prediction_list, fill = TRUE)

# Save CSV
fwrite(
  ld_summary,
  file.path(output_dir, "LD_Decay_Summary.csv")
)

fwrite(
  ld_pred,
  file.path(output_dir, "Fitted_LD_Curves.csv")
)

print(ld_summary)

# ------------------------------------------------------------
# 9. PLOT 1:
# COMBINED GENOME-WIDE + SUBGENOMES
# ------------------------------------------------------------

combined_pred <- ld_pred[
  Level %in% c("Genome", "Subgenome")
]

combined_summary <- ld_summary[
  Level %in% c("Genome", "Subgenome")
]

combined_pred[, Group := factor(
  Group,
  levels = c("Genome-wide", "A", "B", "D")
)]

combined_summary[, Group := factor(
  Group,
  levels = c("Genome-wide", "A", "B", "D")
)]

colors <- c(
  "Genome-wide" = "black",
  "A" = "#0072B2",
  "B" = "#D55E00",
  "D" = "#009E73"
)

p1 <- ggplot(
  combined_pred,
  aes(
    x = Dist_mb,
    y = R2_fit,
    colour = Group
  )
) +
  
  geom_line(linewidth = 1.15) +
  
  geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  
  geom_vline(
    data = combined_summary[
      is.finite(Decay_Mb)
    ],
    aes(
      xintercept = Decay_Mb,
      colour = Group
    ),
    linetype = "dotted",
    linewidth = 0.6,
    show.legend = FALSE
  ) +
  
  scale_colour_manual(
    values = colors,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    limits = c(0, max_distance_mb),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = "Physical distance (Mb)",
    y = expression(r^2),
    colour = NULL
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 11)
  )

ggsave(
  file.path(
    output_dir,
    "01_Genome_Subgenome_LD_Decay.pdf"
  ),
  p1,
  width = 7,
  height = 5
)

# ============================================================
# 10. PLOT 2:
# COMBINED LD DECAY OF ALL 21 WHEAT CHROMOSOMES
# ============================================================

library(ggplot2)

# Extract chromosome-wise fitted curves
chr_pred <- copy(
  ld_pred[Level == "Chromosome"]
)

# Extract chromosome-wise decay estimates
chr_summary <- copy(
  ld_summary[Level == "Chromosome"]
)

# Set chromosome order
chr_pred[, Group := factor(
  Group,
  levels = chr_levels
)]

chr_summary[, Group := factor(
  Group,
  levels = chr_levels
)]

# Generate 21 distinguishable colors
chr_colors <- setNames(
  grDevices::hcl.colors(
    21,
    palette = "Dynamic"
  ),
  chr_levels
)

# ------------------------------------------------------------
# COMBINED CHROMOSOME-WISE LD DECAY PLOT
# ------------------------------------------------------------

p2 <- ggplot(
  chr_pred,
  aes(
    x = Dist_mb,
    y = R2_fit,
    colour = Group,
    group = Group
  )
) +
  
  geom_line(
    linewidth = 0.95,
    alpha = 0.9
  ) +
  
  # LD threshold
  geom_hline(
    yintercept = threshold,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.6
  ) +
  
  scale_colour_manual(
    values = chr_colors,
    drop = FALSE
  ) +
  
  scale_x_continuous(
    limits = c(0, max_distance_mb),
    expand = c(0, 0)
  ) +
  
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.02))
  ) +
  
  labs(
    x = "Physical distance (Mb)",
    y = expression(r^2),
    colour = "Chromosome"
  ) +
  
  guides(
    colour = guide_legend(
      ncol = 3,
      override.aes = list(
        linewidth = 1.5,
        alpha = 1
      )
    )
  ) +
  
  theme_classic(base_size = 14) +
  
  theme(
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 12
    ),
    
    legend.text = element_text(
      size = 10
    ),
    
    axis.title = element_text(
      size = 13
    ),
    
    axis.text = element_text(
      size = 11
    ),
    
    panel.border = element_rect(
      colour = "black",
      fill = NA,
      linewidth = 0.4
    )
  )

# ------------------------------------------------------------
# SAVE FIGURE
# ------------------------------------------------------------

ggsave(
  filename = file.path(
    output_dir,
    "02_Combined_21_Chromosomes_LD_Decay.pdf"
  ),
  plot = p2,
  width = 10,
  height = 6.5
)

# Display
print(p2)
# ------------------------------------------------------------
# 11. DISPLAY PLOTS
# ------------------------------------------------------------

print(p1)
print(p2)

# ------------------------------------------------------------
# 12. FINAL SUMMARY
# ------------------------------------------------------------

cat("\n====================================\n")
cat("LD DECAY ANALYSIS COMPLETED\n")
cat("====================================\n")

print(
  ld_summary[
    ,
    .(
      Level,
      Group,
      Decay_kb,
      N_pairs,
      Status
    )
  ]
)

cat(
  "\nResults saved in:\n",
  normalizePath(output_dir),
  "\n"
)