##### Analysis of enriched motifs at the 
##### targeted methylation sites

#############################################################################################

#############################################################################################

##### Overview:

### 1. Methylation signal around motifs in targeted methylation regions analysis
### 2. Methylation and motif analysis in the Prrca2 locus
### 3. Methylation and motif analysis in the Limd2 locus

#############################################################################################

library(tidyverse)
library(GenomicRanges)
library(plyranges)

load("~/data/adipocyte_methylation/.RData")

motif_names <-
  c(str_c("Motif_down_", c(1:5)), str_c("Motif_up_", c(1, 7)))
motifsOccurence <- read_tsv("sentinel_with_motifs_v1.tsv")
colnames(motifsOccurence)[1] <- "PeakID"

motifOccurenceGR <-
  motifsOccurence %T>% {
    colnames(.)[22:28] <- motif_names
  } %>%
  pivot_longer(cols = 22:28,
               names_to = "motif",
               values_to = "hit") %>%
  filter(!is.na(motif)) %>%
  mutate(
    hit = str_replace_all(hit,  
                          pattern = "(\\d+\\([A-Za-z]+,[*+-],\\d+\\.\\d+\\)),", 
                          replacement = "\\1_")
    ) %>% 
  separate(hit, into = c("occ1", "occ2", "occ3", "occ4"), "_"
           ) %>% 
  pivot_longer(cols = starts_with("occ"), 
               names_to = "ignore", 
               values_to = "hit"
               ) %>%
  filter(!is.na(hit)) %>% tidyr::extract(
    col = "hit",
    into = c("offset", "sequence", "motif_strand"),
    "(\\d+)\\(([A-Za-z]+),([*+-]).+",
    convert = T
  ) %>%
  mutate(motif_start = Start + offset,
         motif_end = Start + offset + 12) %>%
  makeGRangesFromDataFrame(
    start.field = "motif_start",
    end.field = "motif_end",
    strand.field = "motif_strand",
    keep.extra.columns = T
  )

motif4_instances <-
  motifOccurenceGR %>% filter(motif == "Motif_down_4")

library(EnrichedHeatmap)
library(circlize)

targeted_seq_diff <- read_csv("~/data/adipocyte_methylation/final_with_mean.csv") %>%
  separate(col = "cg_site", into = c("chr", "range"), sep = ":") %>%
  separate(col = "range", into = c("pos1", "pos2"), sep = "_")
targeted_seq_GR <- makeGRangesFromDataFrame(targeted_seq_diff, 
                                            start.field = "pos1",
                                            end.field = "pos2", 
                                            keep.extra.columns = T)

motif4_instances_withCpG <- motif4_instances[countOverlaps(motif4_instances + 150,
                                                           targeted_seq_GR) != 0]

meth_col_fun = colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
mat1 = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 10, mean_mode = "absolute", w = 1,
                         background = 0, smooth = F)
rownames(mat1) <- motif4_instances_withCpG$PeakID

mat1_withNA = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 5, mean_mode = "absolute", w = 1,
                         background = NA, smooth = F)
mat1_smooth_withNA = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 150, mean_mode = "absolute", w = 1,
                         background = NA, smooth = T)
mat1_smooth = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 50, mean_mode = "absolute", w = 1,
                         background = 0, smooth = T)
EnrichedHeatmap(mat1_withNA, col = meth_col_fun )

cpg_per_site <- apply(mat1_withNA, 2, function(x) sum(!is.na(x)))
motif4_in_window <- nrow(mat1_withNA)
mean_diff <- apply(mat1_withNA, 2, mean, na.rm = T)

num_corrected <- mean_diff * -(1 / log2( cpg_per_site / motif4_in_window))


The `echo: false` option disables the printing of code (only output is displayed).

Try to repeat it without the problematic region.


motif4_instances_withCpG <- motif4_instances[countOverlaps(motif4_instances + 50,
                                                           targeted_seq_GR) != 0]
motif4_instances_withCpG <- motif4_instances_withCpG[-23]

meth_col_fun = colorRamp2(c(-0.1, 0, 0.1), c("blue", "white", "red"))
mat1 = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 50, mean_mode = "absolute", w = 1,
                         background = 0, smooth = F)
mat1_withNA = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 50, mean_mode = "absolute", w = 1,
                         background = NA, smooth = F)
mat1_smooth_withNA = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 50, mean_mode = "absolute", w = 1,
                         background = NA, smooth = T)
mat1_smooth = normalizeToMatrix(targeted_seq_GR, motif4_instances_withCpG, 
                         value_column = "meth_diff_sc",
                         extend = 50, mean_mode = "absolute", w = 1,
                         background = 0, smooth = T)
EnrichedHeatmap(mat1, col = meth_col_fun )

cpg_per_site <- apply(mat1_withNA, 2, function(x) sum(!is.na(x)))
motif4_in_window <- nrow(mat1_withNA)
mean_diff <- apply(mat1_withNA, 2, mean, na.rm = T)

num_corrected <- mean_diff * -(1 / log2( cpg_per_site / motif4_in_window))


motif4_withCpG <- which((substr(motif4_instances$sequence, 7, 8) == "CG") | (substr(motif4_instances$sequence, 3, 4) == "CG"))

motif4_subset <- motif4_instances[motif4_withCpG]

# get back to 0 later
lost_methylation <- targeted_seq_GR[targeted_seq_GR$meth_diff_sc <= -0.02]

affected_motif4 <- subsetByOverlaps(motif4_subset, lost_methylation)
affected_sentinels <- affected_motif4$PeakID

gene_associations <- read_tsv("inst/extdata/genic_chip_ensamble_association_v3.tsv")

goi <- gene_associations %>% filter(cpg_site %in% affected_sentinels) %>% 
  dplyr::select(GENEID) %>%
  pull()

goi[21] <-  "ENSG00000161939"
goi[26] <- "ENSG00000258315"

library(clusterProfiler)
library(org.Hs.eg.db)

target_genes <- select(org.Hs.eg.db, keys = goi, columns = c("SYMBOL", "GENENAME"), 
                       keytype = "ENSEMBL")

# GO enrichment didn't really work :(
# motif4_universe <- gene_associations %>% filter(cpg_site %in% motif4_instances$PeakID) %>%
#   dplyr::select(GENEID) %>% pull()
# motif4_universe[57] <- "ENSG00000161939"
# motif4_universe[75] <- "ENSG00000258315"
# 
# go <- enrichGO(goi, OrgDb = org.Hs.eg.db, keyType = "ENSEMBL", ont = "ALL", 
#                universe = motif4_universe)

targeted_motif4_diff <- motif4_subset %>% shift(3) %>% 
  resize(width = 7, fix = "start") %>% mergeByOverlaps(targeted_seq_GR) %>% 
  as_tibble() %>% dplyr::select(PeakID, meth_diff_sc)
colnames(targeted_motif4_diff) <- c("CG", "targeted_diff")

array_motif4_diff <- sentAndNot$sentinel %>% as_tibble() %>% 
  dplyr::select(CG, comb_beta) %>% filter(CG %in% motif4_subset$PeakID)

full_join(targeted_motif4_diff, array_motif4_diff) %>% ggplot(aes(comb_beta, targeted_diff)) + geom_point() + geom_smooth(method ="lm") + stat_cor() + labs(y= "Targeted sequencing", x = "Array sentinel", title = "Correlation of methylation difference between array sentinel and targeted methylation in motif") + xlim(c(-0.2, 0.05)) + ylim(c(-0.2, 0.05))

gene_names_closest <- motif4_instances$`Gene Name`
names(gene_names_closest) <- motif4_instances$PeakID

index <- which(rownames(mat1) %in% motif4_subset$PeakID)
labels_closest <- gene_names_closest[rownames(mat1)[index]]

EnrichedHeatmap(mat1, col = meth_col_fun ) +
rowAnnotation(clos_gene = anno_mark(at = index, labels = labels_closest)) +
rowAnnotation(clos_gene = anno_mark(at = own_index, labels = labels_own))

library(JASPAR2022)
library(TFBSTools)

opts <- list(species = "Homo sapiens",
             collection = "CORE",
             tax_group = "vertebrates",
             all_versions = F)
human_mats <- getMatrixSet(JASPAR2022, opts)
prrc2a <- GRanges("chr6:31607410-31607660")
cpg_in_prrca2a <- subsetByOverlaps(targeted_seq_GR, prrc2a)

library(BSgenome.Hsapiens.UCSC.hg19)
prrc2a_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, cpg_in_prrca2a + 6)

seqsInPrrc2aCpG <- searchSeq(toPWM(human_mats), prrc2a_seqs, min.score = "85%")
prrc2a_tf_match <- writeGFF3(seqsInPrrc2aCpG)

prrc2a_tfs <- prrc2a_tf_match$attributes %>% str_replace("TF=", "") %>% 
  str_replace(";.+", "") %>% unique()

tmp1 <- prrc2a_tfs %>% str_split("::", simplify = T)
tmp2 <- tmp1[tmp1[,2] != "" ,2 ]
prrc2a_tfs <- unique(c(tmp1[,1], tmp2))

prrc2a_tf_ensid <- select(org.Hs.eg.db, keys = prrc2a_tfs, 
                          columns = c("ENSEMBL", "GENENAME"),
                          keytype = "SYMBOL")

expr <- read_rds("~/data/transformed_counts_sa_va.rds")

prrc2a_ensid <- "ENSG00000204469"

prrc2a_expr <- expr %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>%  
  filter(str_detect(sample, "SA")) %>%
  dplyr::select(any_of(c(prrc2a_ensid, prrc2a_tf_ensid$ENSEMBL)))

library(ggpubr)

prrc2a_expr %>% rowid_to_column() %>% pivot_longer(cols = -c(1, 2), names_to = "TF", values_to = "expr") %>% ggplot(aes(ENSG00000204469, expr)) + geom_point() + geom_smooth(method = "lm") + stat_cor() + facet_wrap(~ TF, scales = "free_y")

prrc2a_expr %>% rowid_to_column() %>% pivot_longer(cols = -c(1, 2), names_to = "ENSEMBL", values_to = "expr") %>% full_join(prrc2a_tf_ensid) %>% ggplot(aes(ENSG00000204469, expr)) + geom_point() + geom_smooth(method = "lm") + stat_cor() + facet_wrap(~ SYMBOL, scales = "free_y")

elk1 <- "ENSG00000126767"
essra <- "ENSG00000173153"

prrc2a_cor_scatter <- prrc2a_tf_expression %>% filter(SYMBOL %in% c("ELK1", "ESRRA")) %>% 
     ggplot(aes(ENSG00000204469, expr)) + 
     geom_point() + geom_smooth(method = "lm") + stat_cor() + 
     facet_wrap( ~ SYMBOL, scales = "free_y") + 
     xlab("PRRC2A expression") + ylab("Associated TF expression")

prrc2a_tf_expression %>% group_by(SYMBOL) %>% 
     summarise(cor = tidy(cor.test(ENSG00000204469, expr))) %>% unnest(cols = cor) %>% 
     mutate(SYMBOL = fct_reorder(SYMBOL, estimate, .desc = T)) %>% 
     ggplot(aes(SYMBOL, estimate, fill = -log(p.value))) + geom_col() + 
     theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
     annotation_custom(ggplotGrob(prrc2a_cor_scatter), 
                       xmin = 18, xmax = 43, ymin = 0.1, ymax = 0.7) + 
  scale_fill_viridis() + 
  xlab("Transcription factor") + ylab("Correlation with PRRC2A expression")

prrc2a_sentinel_methylation <- replication_betaMethylation_subset["cg15809217", ] %>%
  data.frame %>% rownames_to_column(var = "Array_ID") %>% as_tibble()

elk1_expression <- expr[elk1, ] %>% data.frame %>% 
  rownames_to_column(var = "Sample.ID") %>% as_tibble()

prrc2a_sentinel_methylation <- rename(prrc2a_sentinel_methylation, Methylation = `.`)

elk1_expression <- rename(elk1_expression, ELK1 = `.`)

samples_arrayID <- replication_phe %>% dplyr::select(Sample.ID, Array_ID) %>% mutate(Sample.ID = str_replace(Sample.ID, " ", ".")) %>% full_join(prrc2a_sentinel_methylation) %>% full_join(elk1_expression) %>% drop_na() %>% filter(str_detect(Sample.ID, "SA"))

samples_arrayID %>% mutate(across(c(3, 4), ~ scale(.x))) %>% ggplot(aes(ELK1, Methylation)) + geom_point() + geom_smooth(method = "lm") + stat_cor()


all_motif4_genes <- unique(c(affected_motif4$`Nearest Ensembl`, goi$GENEID))
all_motif4_gene_names <- select(org.Hs.eg.db, 
                                keys = all_motif4_genes, 
                                columns = c("SYMBOL", "GENENAME"), 
                                keytype = "ENSEMBL")

missing_symbol <- is.na(all_motif4_gene_names$SYMBOL)
all_motif4_gene_names$SYMBOL[missing_symbol] <- all_motif4_gene_names$ENSEMBL[missing_symbol]
all_motif4_gene_names <- all_motif4_gene_names[
  all_motif4_gene_names$ENSEMBL %in% rownames(expr), ]


tf_names <- c("ELK1", "ELK3", "ELF4", "ELF1")
motif4_tf_names <- select(org.Hs.eg.db,
                          keys = tf_names, 
                          columns = c("ENSEMBL", "GENENAME"),
                          keytype = "SYMBOL")

motif4_tf_names <- motif4_tf_names[ , c(2, 1, 3)]

all_motif4_gene_names$role <- "TARGET"
motif4_tf_names$role <- "TF"

# motif4_regulatory_network <- rbind(all_motif4_gene_names, motif4_tf_names)
# motif4_regulatory_network
# as.data.frame(expr[motif4_regulatory_network$ENSEMBL , ])
# motif4_regulatory_network <- cbind(motif4_regulatory_network,
# as.data.frame(expr[motif4_regulatory_network$ENSEMBL , ]))
# motif4_regulatory_network
# motif4_regulatory_network %>% pivot_longer(cols = -c(1:4))
# motif4_regulatory_network %>% pivot_longer(cols = -c(1:4), names_to = "sample", values_to = "expression")
# motif4_regulatory_network %>% pivot_longer(cols = -c(1:4), names_to = "sample", values_to = "expression") %>% filter(str_detect(sample, "SA"))
all_motif4_gene_names <- cbind(all_motif4_gene_names,
                               as.data.frame(expr[all_motif4_gene_names$ENSEMBL , ])) %>%
  pivot_longer(cols = -c(1:4), names_to = "sample", values_to = "expression")

motif4_tf_names <- cbind(motif4_tf_names,
                         as.data.frame(expr[motif4_tf_names$ENSEMBL , ])) %>%
  pivot_longer(cols = -c(1:4), names_to = "sample", values_to = "expression")

# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample")
# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") %>% ggplot(aes(expression.x, expression.y)) + geom_point() + geom_smooth(method = "lm") + facet_wrap(SYMBOL.x ~ SYMBOL.y)
# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") %>% ggplot(aes(expression.x, expression.y)) + geom_point() + geom_smooth(method = "lm") + facet_grid(SYMBOL.x ~ SYMBOL.y)
# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") %>% ggplot(aes(expression.x, expression.y)) + geom_point() + geom_smooth(method = "lm") + facet_grid(SYMBOL.y ~ SYMBOL.x)
# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") %>% ggplot(aes(expression.x, expression.y)) + geom_point() + geom_smooth(method = "lm") + facet_grid(SYMBOL.y ~ SYMBOL.x, scales = "free")
# full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") %>% filter(str_detect(sample, "SA")) %>% ggplot(aes(expression.x, expression.y)) + geom_point() + geom_smooth(method = "lm") + facet_grid(SYMBOL.y ~ SYMBOL.x, scales = "free")

motif4_tf_target_expression <- full_join(motif4_tf_names, all_motif4_gene_names, by = "sample") 
# motif4_tf_target_expression %>% group_by(SYMBOL.x, SYMBOL.y) %>% 
#      summarise(cor = tidy(cor.test(expression.x, expression.y))) %>% unnest(cols = cor) %>% 
#      mutate(SYMBOL.y = fct_reorder(SYMBOL.y, estimate, .desc = T)) %>% 
#      ggplot(aes(SYMBOL.y, estimate, fill = -log(p.value))) + geom_col() + 
#      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#      facet_grid(SYMBOL.x ~ ., scales = "free") + scale_fill_viridis() + 
#      xlab("Transcription factor") + ylab("Correlation with PRRC2A expression")

motif4_tf_tg_cor <- motif4_tf_target_expression %>% group_by(SYMBOL.x, SYMBOL.y) %>% 
  summarise(cor = tidy(cor.test(expression.x, expression.y))) %>% ungroup() %>% 
  split(.$SYMBOL.x) %>% map(function(x) {
    x %>%  unnest(cols = cor) %>% 
      mutate(SYMBOL.y = fct_reorder(SYMBOL.y, estimate, .desc = T)) %>% 
      ggplot(aes(SYMBOL.y, estimate, fill = -log(p.value))) + geom_col() + 
      theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
      facet_grid(SYMBOL.x ~ ., scales = "free") + scale_fill_viridis() + 
      xlab("Target gene") + ylab("Expression correlation")
  })

library(cowplot)

cowplot::plot_grid(motif4_tf_tg_cor$ELF1,
                   motif4_tf_tg_cor$ELF4,
                   motif4_tf_tg_cor$ELK1,
                   motif4_tf_tg_cor$ELK3, ncol = 1)



library(JASPAR2022)
library(TFBSTools)

opts <- list(species = "Homo sapiens",
             collection = "CORE",
             tax_group = "vertebrates",
             all_versions = F)
human_mats <- getMatrixSet(JASPAR2022, opts)
limd2 <- GRanges("chr17:61773942-61774230")
limd2_sentinel <- GRanges("chr17:61774174-61774175")
limd2_sentinel$name <- "cg05941027"
cpg_in_limd2 <- subsetByOverlaps(c(targeted_seq_GR, limd2_sentinel), limd2)

library(BSgenome.Hsapiens.UCSC.hg19)
limd2_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, cpg_in_limd2 + 6)

seqsInLimd2CpG <- searchSeq(toPWM(human_mats), limd2_seqs, min.score = "85%")
limd2_tf_match <- writeGFF3(seqsInLimd2CpG)

limd2_tfs <- limd2_tf_match$attributes %>% str_replace("TF=", "") %>% 
  str_replace(";.+", "") %>% unique()

tmp1 <- limd2_tfs %>% str_split("::", simplify = T)
tmp2 <- tmp1[tmp1[,2] != "" ,2 ]
prrc2a_tfs <- unique(c(tmp1[,1], tmp2))

library(org.Hs.eg.db)
limd2_tf_ensid <- select(org.Hs.eg.db, keys = limd2_tfs, 
                          columns = c("ENSEMBL", "GENENAME"),
                          keytype = "SYMBOL")

expr <- read_rds("~/data/adipocyte_methylation/transformed_counts_sa_va.rds")

limd2_ensid <- "ENSG00000136490"

limd2_expr <- expr %>% t() %>% as.data.frame() %>% rownames_to_column(var = "sample") %>%  
  filter(str_detect(sample, "VA")) %>%
  dplyr::select(any_of(c(limd2_ensid, limd2_tf_ensid$ENSEMBL)))

library(ggpubr)

limd2_expr %>% rowid_to_column() %>% pivot_longer(cols = -c(1, 2), names_to = "TF", values_to = "expr") %>% ggplot(aes(ENSG00000136490, expr)) + geom_point() + geom_smooth(method = "lm") + stat_cor() + facet_wrap(~ TF, scales = "free_y")

limd2_expr %>% rowid_to_column() %>% pivot_longer(cols = -c(1, 2), names_to = "ENSEMBL", values_to = "expr") %>% full_join(limd2_tf_ensid) %>% drop_na() %>% ggplot(aes(ENSG00000136490, expr)) + geom_point() + geom_smooth(method = "lm") + stat_cor() + facet_wrap(~ SYMBOL, scales = "free_y")

tfap2e <- "ENSG00000116819"
meis3 <- "ENSG00000105419"

# expr %>% t() %>% as_tibble() %>% dplyr::select(any_of(c(elk1, essra))) %>% ggplot(aes(ENSG00000126767, ENSG00000173153)) + geom_smooth(method = "lm") + geom_point() + stat_cor()

limd2_tf_expression <- limd2_expr %>% rowid_to_column() %>% 
  pivot_longer(cols = -c(1, 2), names_to = "ENSEMBL", values_to = "expr") %>% 
  full_join(limd2_tf_ensid)

limd2_cor_scatter <- limd2_tf_expression %>% filter(SYMBOL %in% c("TFAP2E", "MEIS3")) %>% 
     ggplot(aes(ENSG00000136490, expr)) + 
     geom_point() + geom_smooth(method = "lm") + stat_cor() + 
     facet_wrap( ~ SYMBOL, scales = "free_y") + 
     xlab("LIMD2 expression") + ylab("Associated TF expression")

library(broom)
library(viridis)

limd2_tf_expression %>% group_by(SYMBOL) %>% drop_na() %>%
     summarise(cor = tidy(cor.test(ENSG00000136490, expr))) %>% unnest(cols = cor) %>% 
     mutate(SYMBOL = fct_reorder(SYMBOL, estimate, .desc = T)) %>% 
     ggplot(aes(SYMBOL, estimate, fill = -log(p.value))) + geom_col() + 
     theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
     annotation_custom(ggplotGrob(limd2_cor_scatter), 
                       xmin = 18, xmax = 35, ymin = 0.15, ymax = 0.75) + 
  scale_fill_viridis() + 
  xlab("Transcription factor") + ylab("Correlation with LIMD2 expression")

prrc2a_sentinel_methylation <- replication_betaMethylation_subset["cg15809217", ] %>%
  data.frame %>% rownames_to_column(var = "Array_ID") %>% as_tibble()

elk1_expression <- expr[elk1, ] %>% data.frame %>% 
  rownames_to_column(var = "Sample.ID") %>% as_tibble()

prrc2a_sentinel_methylation <- rename(prrc2a_sentinel_methylation, Methylation = `.`)

elk1_expression <- rename(elk1_expression, ELK1 = `.`)

samples_arrayID <- replication_phe %>% dplyr::select(Sample.ID, Array_ID) %>% mutate(Sample.ID = str_replace(Sample.ID, " ", ".")) %>% full_join(prrc2a_sentinel_methylation) %>% full_join(elk1_expression) %>% drop_na() %>% filter(str_detect(Sample.ID, "SA"))

samples_arrayID %>% mutate(across(c(3, 4), ~ scale(.x))) %>% ggplot(aes(ELK1, Methylation)) + geom_point() + geom_smooth(method = "lm") + stat_cor()
