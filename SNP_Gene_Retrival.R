# ============================================================
# COMPLETE SCRIPT
# SNP -> GENES WITHIN ±1.5 kb
# + v2.1 ID
# + v1.1 ID
# + SNP gene feature
# ============================================================

library(readxl)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(writexl)


# ============================================================
# 1. SETTINGS
# ============================================================

SNP_FILE <- "snp_list.xlsx"
GFF_FILE <- "/home/t15q981/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3"

WINDOW <- 10001000


# ============================================================
# 2. READ SNP FILE
# ============================================================

snps <- read_excel(SNP_FILE)



# ============================================================
# 3. READ GFF3
# ============================================================


gff <- import(GFF_FILE)




# ============================================================
# 4. EXTRACT GENES
# ============================================================

genes <- gff[gff$type == "gene"]


# ============================================================
# 6. CREATE SNP POINTS
# ============================================================

snp_point <- GRanges(
  seqnames = snps$chr,
  
  ranges = IRanges(
    start = snps$pos,
    end   = snps$pos
  ),
  
  SNP = snps$snp,
  SNP_position = snps$pos
)


# ============================================================
# 7. CREATE ±1.5 kb WINDOWS
# ============================================================

snp_window <- GRanges(
  seqnames = snps$chr,
  
  ranges = IRanges(
    start = pmax(1, snps$pos - WINDOW),
    end   = snps$pos + WINDOW
  ),
  
  SNP = snps$snp,
  SNP_position = snps$pos
)


# ============================================================
# 8. FIND GENES WITHIN ±1.5 kb
# ============================================================

gene_hits <- findOverlaps(
  snp_window,
  genes
)

cat(
  "\nTotal SNP-gene combinations:",
  length(gene_hits),
  "\n"
)


# ============================================================
# 9. BUILD SNP-GENE TABLE
# ============================================================

results <- data.frame(
  
  SNP =
    mcols(snp_window)$SNP[
      queryHits(gene_hits)
    ],
  
  Chr =
    as.character(
      seqnames(snp_window)[queryHits(gene_hits)]
    ),
  
  SNP_position =
    mcols(snp_window)$SNP_position[
      queryHits(gene_hits)
    ],
  
  # IWGSC v2.1
  Gene_ID_v2.1 =
    as.character(
      mcols(genes)$ID[
        subjectHits(gene_hits)
      ]
    ),
  
  # IWGSC v1.1
  Gene_ID_v1.1 =
    as.character(
      mcols(genes)$previous_id[
        subjectHits(gene_hits)
      ]
    ),
  
  Gene_start =
    start(genes)[
      subjectHits(gene_hits)
    ],
  
  Gene_end =
    end(genes)[
      subjectHits(gene_hits)
    ],
  
  Strand =
    as.character(
      strand(genes)[
        subjectHits(gene_hits)
      ]
    ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 10. DISTANCE BETWEEN SNP AND GENE
# ============================================================

results <- results %>%
  mutate(
    
    Distance_to_gene = case_when(
      
      SNP_position < Gene_start ~
        Gene_start - SNP_position,
      
      SNP_position > Gene_end ~
        SNP_position - Gene_end,
      
      # SNP is inside gene
      TRUE ~ 0
    )
  )


# ============================================================
# 11. GENE LOCATION
#     Strand-aware upstream/downstream
# ============================================================

results <- results %>%
  mutate(
    
    Gene_location = case_when(
      
      # SNP inside gene
      SNP_position >= Gene_start &
        SNP_position <= Gene_end ~
        "Gene body",
      
      # + strand
      Strand == "+" &
        SNP_position < Gene_start ~
        "Upstream",
      
      Strand == "+" &
        SNP_position > Gene_end ~
        "Downstream",
      
      # - strand
      Strand == "-" &
        SNP_position > Gene_end ~
        "Upstream",
      
      Strand == "-" &
        SNP_position < Gene_start ~
        "Downstream",
      
      TRUE ~ "Unknown"
    )
  )


# ============================================================
# 13. EXTRACT CDS, EXON AND UTR FEATURES
# ============================================================

features <- gff[
  gff$type %in% c(
    "CDS",
    "exon",
    "five_prime_UTR",
    "three_prime_UTR",
    "five_prime_utr",
    "three_prime_utr",
    "5UTR",
    "3UTR"
  )
]


# ============================================================
# 14. FIND WHICH FEATURES CONTAIN EACH SNP
# ============================================================

feature_hits <- findOverlaps(
  snp_point,
  features
)


feature_table <- data.frame(
  
  SNP =
    mcols(snp_point)$SNP[
      queryHits(feature_hits)
    ],
  
  Feature =
    as.character(
      features$type[
        subjectHits(feature_hits)
      ]
    ),
  
  stringsAsFactors = FALSE
)


# ============================================================
# 15. STANDARDIZE FEATURE NAMES
# ============================================================

feature_table <- feature_table %>%
  mutate(
    
    Feature = case_when(
      
      Feature %in% c(
        "five_prime_UTR",
        "five_prime_utr",
        "5UTR"
      ) ~ "5'UTR",
      
      Feature %in% c(
        "three_prime_UTR",
        "three_prime_utr",
        "3UTR"
      ) ~ "3'UTR",
      
      TRUE ~ Feature
    )
  )


# ============================================================
# 16. ASSIGN ONE FEATURE PER SNP
#
# Priority:
# 5'UTR
# 3'UTR
# CDS
# Exon
# ============================================================

feature_summary <- feature_table %>%
  group_by(SNP) %>%
  
  summarise(
    
    SNP_feature = case_when(
      
      any(Feature == "5'UTR") ~
        "5'UTR",
      
      any(Feature == "3'UTR") ~
        "3'UTR",
      
      any(Feature == "CDS") ~
        "CDS",
      
      any(Feature == "exon") ~
        "Exon",
      
      TRUE ~ NA_character_
    ),
    
    .groups = "drop"
  )


# ============================================================
# 17. JOIN FEATURE TO SNP-GENE TABLE
# ============================================================

results <- results %>%
  left_join(
    feature_summary,
    by = "SNP"
  )


# ============================================================
# 18. ASSIGN INTRON / UPSTREAM / DOWNSTREAM
# ============================================================

results <- results %>%
  mutate(
    
    SNP_feature = case_when(
      
      Gene_location == "Upstream" ~
        "Upstream",
      
      Gene_location == "Downstream" ~
        "Downstream",
      
      # SNP inside gene but not CDS/exon/UTR
      Gene_location == "Gene body" &
        is.na(SNP_feature) ~
        "Intron",
      
      TRUE ~ SNP_feature
    )
  )


write.csv( results, "SNP_Gene_Annotation_1Mb.csv", row.names=FALSE)
