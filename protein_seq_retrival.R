library(Biostrings)

# Files
gene_file <- "gene_list.txt"
pep_file  <- "IWGSC_RefSeq_v2.1_HC_pep.fa"  # change filename if needed

# ------------------------------------------------------------
# Read gene/transcript list
# DO NOT remove .1
# ------------------------------------------------------------

genes <- suppressWarnings(readLines(gene_file))
genes <- trimws(genes)
genes <- genes[genes != ""]

cat("Number of IDs:", length(genes), "\n")

# Check
head(genes)
tail(genes)


# ------------------------------------------------------------
# Read peptide FASTA
# ------------------------------------------------------------

pep <- readAAStringSet(pep_file)

cat("Number of proteins in FASTA:", length(pep), "\n")

head(names(pep))


# ------------------------------------------------------------
# EXACT MATCH
# ------------------------------------------------------------

idx <- match(genes, names(pep))

cat("Matched:", sum(!is.na(idx)), "\n")
cat("Not matched:", sum(is.na(idx)), "\n")


# Show IDs that did not match
missing <- genes[is.na(idx)]

if(length(missing) > 0) {
  print(missing)
  writeLines(missing, "missing_peptide_IDs.txt")
}


# ------------------------------------------------------------
# Extract proteins
# ------------------------------------------------------------

selected <- pep[idx[!is.na(idx)]]

cat("Extracted:", length(selected), "\n")


# ------------------------------------------------------------
# Save all
# ------------------------------------------------------------

writeXStringSet(
  selected,
  "selected_peptides.fa",
  format = "fasta"
)


# ------------------------------------------------------------
# Split into two files
# ------------------------------------------------------------

n <- length(selected)
half <- ceiling(n / 2)

part1 <- selected[1:half]
part2 <- selected[(half + 1):n]

writeXStringSet(
  part1,
  "selected_peptides_part1.fa",
  format = "fasta"
)

writeXStringSet(
  part2,
  "selected_peptides_part2.fa",
  format = "fasta"
)

cat("\nPart 1:", length(part1), "\n")
cat("Part 2:", length(part2), "\n")
