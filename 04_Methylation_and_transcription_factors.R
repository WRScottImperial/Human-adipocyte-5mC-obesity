##### Interactions between DNA methylation and transcription factors 
##### in human subcutaneous adipocytes

#############################################################################################

#############################################################################################

##### Overview:

### 1. Motif enrichment in sentinel region
### 2. Ratio of observed vs expected values of motif presence
### 3. Location of CG dinucleotide ariubd motifs in differentially methylated regions

#############################################################################################

library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(pbapply)
library(magrittr)

motif_instances_occurence_tidy <- 
  read_tsv("sentinel_with_motifs_occurence_tidy.tsv")

mosaNewBackground <- 
  read_csv("~/Downloads/mosa_1k_permutations_sliding_msd_dist_matches.csv") %>%
  column_to_rownames(var = "CG")

expected <- motif_instances_occurence_tidy %>%
  group_by(Motif) %>% nest() %>%
  mutate(expected = map(data, function(x){
    tmp_bg <- mosaNewBackground[x$PeakID, ]
    pbapply(tmp_bg, 2, function(xx) {
      table(Islands.UCSC[xx, "Relation_to_Island"] %>% 
              str_replace(pattern = "[SN]_", replacement = "") %>%
              factor(levels = c("Island", "OpenSea", "Shelf", "Shore")))
    }) # %>% { .[is.na(.)] <- 0 } # %>%
    # rowMeans() %>% round() %>% t() %>% as.data.frame()
  })) # %>% unnest(expected)

expected_1 <- expected %>%
  mutate(expected = map(expected, function(x){
    as.data.frame(t(round(rowMeans(x))))
  })) %>% unnest(expected)

observed <- motif_instances_occurence_tidy %>%
  group_by(Motif) %>% nest() %>%
  mutate(observed = map(data, function(x){
    table(Islands.UCSC[x$PeakID, "Relation_to_Island"] %>% 
            str_replace(pattern = "[SN]_", replacement = "")) %>%
      as.data.frame() %>% pivot_wider(names_from = "Var1", values_from = "Freq")
  })) %>% unnest(observed) %>% replace_na(list(Island = 0))

observed_mat <- observed[,3:6] %>% 
  as.matrix() %T>% { rownames(.) <- observed$Motif }

expected_mat <- expected_1[,3:6] %>% 
  as.matrix() %T>% { rownames(.) <- expected_1$Motif }

motif_num <- rowSums(observed_mat)

cgi_pvals <- t(sapply(1:7, function(i){
  sapply(1:4, function(j){
    obs <- observed_mat[i, j]
    expt <- expected_mat[i, j]
    obs_resid <- motif_num[i] - obs
    expt_resid <- motif_num[i] - expt
    fisher.test(matrix(c(obs, obs_resid, expt, expt_resid), ncol = 2))$p.value
  })
}))

rownames(cgi_pvals) <- rownames(observed_mat)
colnames(cgi_pvals) <- colnames(observed_mat)

chisq <- chisq.test(observed_mat)
obs_over_expected <- observed_mat / expected_mat

xl <- observed_mat %>% as.data.frame() %>% rownames_to_column(var = "motif") %>%
  pivot_longer(cols = -1, names_to = "cgi", values_to = "count")
yl <- obs_over_expected %>% as.data.frame() %>% 
  rownames_to_column(var = "motif") %>% 
  pivot_longer(cols = -1, names_to = "cgi", values_to = "ratio")
csgi_counts <-  plyr::join(xl, yl, by = c("motif", "cgi"))

cgi_plot <- csgi_counts %>% 
  ggplot(aes(cgi, motif, colour = ratio, size = count)) +
  geom_point() + scale_color_distiller(palette = "RdBu") + 
  scale_size(range = c(0, 20)) + theme_bw()


roadmap_assoc <- read_csv("~/Downloads/k450_Roadmap_Chromatin_States.csv")
roadmap_assoc <- column_to_rownames(roadmap_assoc, "cg")

roadmap_annotation <- c(
  "Active TSS",
  "Flanking Active TSS",
  "Transcr. at gene 5' and 3'",
  "Strong transcription",
  "Weak transcription",
  "Genic enhancers",
  "Enhancers",
  "ZNF genes & repeats",
  "Heterochromatin",
  "Bivalent/Poised",
  "Flanking Bivalent TSS/Enh",
  "Bivalent Enhancer",
  "Repressed PolyComb",
  "Weak Repressed PolyComb",
  "Quiescent/Low"
)

names(roadmap_annotation) <- str_c("E", 1:15)

expected_roadmap <- motif_instances_occurence_tidy %>%
  group_by(Motif) %>% nest() %>%
  mutate(expected = map(data, function(x){
    tmp_bg <- mosaNewBackground[x$PeakID, ]
    pbapply(tmp_bg, 2, function(xx) {
      table(roadmap_assoc[xx, "RM_E063_primary"] %>% 
              factor(levels = str_c("E", c(1:15))))
    }) # %>% { .[is.na(.)] <- 0 } # %>%
    # rowMeans() %>% round() %>% t() %>% as.data.frame()
  })) # %>% unnest(expected)

expected_roadmap_1 <- expected_roadmap %>%
  mutate(expected = map(expected, function(x){
    as.data.frame(t(round(rowMeans(x))))
  })) %>% unnest(expected)

observed_roadmap <- motif_instances_occurence_tidy %>%
  group_by(Motif) %>% nest() %>%
  mutate(observed = map(data, function(x){
    table(roadmap_assoc[x$PeakID, "RM_E063_primary"]) %>% 
      as.data.frame() %>% pivot_wider(names_from = "Var1", values_from = "Freq")
  })) %>% unnest(observed) %>% mutate_all(function(x) replace_na(x, 0))

observed_roadmap_filt <- observed_roadmap %>% select(c("Motif", "data",
                                                       str_c("E", c(1, 2, 4, 5, 7, 11, 12, 13, 14, 15))))

expected_roadmap_filt <- expected_roadmap_1 %>% select(c("Motif", "data",
                                                         str_c("E", c(1, 2, 4, 5, 7, 11, 12, 13, 14, 15))))

observed_roadmap_mat <- observed_roadmap_filt[,3:12] %>% 
  as.matrix() %T>% { rownames(.) <- observed_roadmap_filt$Motif }

expected_roadmap_mat <- expected_roadmap_filt[,3:12] %>% 
  as.matrix() %T>% { rownames(.) <- expected_roadmap_filt$Motif }

motif_num <- rowSums(observed_mat)

roadmap_pvals <- t(sapply(1:7, function(i){
  sapply(1:10, function(j){
    obs <- observed_roadmap_mat[i, j]
    expt <- expected_roadmap_mat[i, j]
    obs_resid <- motif_num[i] - obs
    expt_resid <- motif_num[i] - expt
    fisher.test(matrix(c(obs, obs_resid, expt, expt_resid), ncol = 2))$p.value
  })
}))

rownames(roadmap_pvals) <- rownames(observed_roadmap_mat)
colnames(roadmap_pvals) <- colnames(observed_roadmap_mat)

chisq_roadmap <- chisq.test(observed_roadmap_mat)
obs_over_expected_roadmap <- observed_roadmap_mat / (expected_roadmap_mat + 1)

xl <- observed_roadmap_mat %>% as.data.frame() %>% 
  rownames_to_column(var = "motif") %>%
  pivot_longer(cols = -1, names_to = "roadmap", values_to = "count")
yl <- obs_over_expected_roadmap %>% as.data.frame() %>% 
  rownames_to_column(var = "motif") %>% 
  pivot_longer(cols = -1, names_to = "roadmap", values_to = "ratio")
roadmap_counts <-  plyr::join(xl, yl, by = c("motif", "roadmap"))

roadmap_plot <- roadmap_counts %>% 
  mutate(roadmap = factor(roadmap, levels = 
                            str_c("E", c(1, 2, 4, 5, 7, 11, 12, 13, 14, 15)))) %>%
  ggplot(aes(roadmap, motif, colour = ratio, size = count)) +
  geom_point() + scale_color_distiller(palette = "RdBu") + 
  scale_size(range = c(0, 20)) + theme_bw() + 
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

cowplot::plot_grid(cgi_plot, roadmap_plot, rel_widths = c(1, 2))

motif1_bars <- roadmap_counts %>% filter(motif == "Motif_down_1") %>% 
  mutate(expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"),
         full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio),
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = expected, label = count)) +
  geom_col() + geom_text(aes(hjust = hjust_aes)) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) + theme_bw() + 
  scale_fill_manual(values = c("red3", "royalblue4")) +
  labs(x = "Roadmap state", y = "", fill = "Expectation direction") + coord_flip()

motif2_bars <- roadmap_counts %>% filter(motif == "Motif_down_4") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"), 
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = expected, label = count)) +
  geom_col() + geom_text(aes(hjust = hjust_aes)) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_fill_manual(values = c("red3", "royalblue4")) +
  labs(y = "Observed over expected ratio", x = "",
       fill = "Expectation direction") + coord_flip()

motif3_bars <- roadmap_counts %>% filter(motif == "Motif_down_5") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"), 
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = expected, label = count)) +
  geom_col() + geom_text(aes(hjust = hjust_aes)) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_fill_manual(values = c("red3", "royalblue4")) +
  labs(x = "", y = "", fill = "Expectation direction") + coord_flip()

prow <- plot_grid(
  motif1_bars + theme(legend.position="none"),
  motif2_bars + theme(legend.position="none"),
  motif3_bars + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

legend_b <- get_legend(
  motif1_bars + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

pdf("roadmap_enrichment_v1.pdf", width = 18, height = 6)
plot_grid(prow, legend_b, ncol = 1, rel_heights = c(1, .1))
dev.off()

motif1_lolli <- roadmap_counts %>% filter(motif == "Motif_down_1") %>% 
  mutate(expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"),
         full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio)) %>% 
  ggplot(aes(full_names, ratio - 1, color = expected, label = count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + geom_text(color = "white") + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) + theme_bw() + 
  scale_color_manual(values = c("red3", "royalblue4")) +
  labs(x = "Roadmap state", y = "", color = "Expectation direction") + coord_flip()

motif2_lolli <- roadmap_counts %>% filter(motif == "Motif_down_4") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected")) %>% 
  ggplot(aes(full_names, ratio - 1, color = expected, label = count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + geom_text(color = "white") + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_color_manual(values = c("red3", "royalblue4")) +
  labs(y = "Observed over expected ratio", x = "",
       color = "Expectation direction") + coord_flip()

motif3_lolli <- roadmap_counts %>% filter(motif == "Motif_down_5") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected")) %>% 
  ggplot(aes(full_names, ratio - 1, color = expected, label = count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + geom_text(color = "white") + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_color_manual(values = c("red3", "royalblue4")) +
  labs(x = "", y = "", color = "Expectation direction") + coord_flip()

prow_lolli <- plot_grid(
  motif1_lolli + theme(legend.position="none"),
  motif2_lolli + theme(legend.position="none"),
  motif3_lolli + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

pdf("roadmap_enrichment_v2.pdf", width = 18, height = 6)
plot_grid(prow_lolli, legend_b, ncol = 1, rel_heights = c(1, .1))
dev.off()

motif1_bars_v2 <- roadmap_counts %>% filter(motif == "Motif_down_1") %>% 
  mutate(expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"),
         full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio),
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = count)) +
  geom_col()  + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) + theme_bw() + 
  scale_fill_distiller(palette = "Blues", direction = -11) +
  labs(x = "Roadmap state", y = "", fill = "Count") + coord_flip()

motif2_bars_v2 <- roadmap_counts %>% filter(motif == "Motif_down_4") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"), 
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = count)) +
  geom_col() + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(y = "Observed over expected ratio", x = "",
       fill = "Count") + coord_flip()

motif3_bars_v2 <- roadmap_counts %>% filter(motif == "Motif_down_5") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"), 
         hjust_aes = ifelse(ratio > 1, 0, 1)) %>% 
  ggplot(aes(full_names, ratio - 1, fill = count)) +
  geom_col() + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_fill_distiller(palette = "Blues", direction = 1) +
  labs(x = "", y = "", fill = "Count") + coord_flip()

prow_v2 <- plot_grid(
  motif1_bars_v2 + theme(legend.position="none"),
  motif2_bars_v2 + theme(legend.position="none"),
  motif3_bars_v2 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

legend_b_v2 <- get_legend(
  motif1_bars_v2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

pdf("roadmap_enrichment_v3.pdf", width = 18, height = 6)
plot_grid(prow_v2, legend_b_v2, ncol = 1, rel_heights = c(1, .1))
dev.off()

motif1_lolli_v2 <- roadmap_counts %>% filter(motif == "Motif_down_1") %>% 
  mutate(expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected"),
         full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio)) %>% 
  ggplot(aes(full_names, ratio - 1, color =  count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) + theme_bw() + 
  scale_color_distiller(palette = "Blues", direction = 1) +
  labs(x = "Roadmap state", y = "", color = "Expectation direction") + coord_flip()

motif2_lolli_v2 <- roadmap_counts %>% filter(motif == "Motif_down_4") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected")) %>% 
  ggplot(aes(full_names, ratio - 1, color = count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_color_distiller(palette = "Blues", direction = 1) +
  labs(y = "Observed over expected ratio", x = "",
       color = "Expectation direction") + coord_flip()

motif3_lolli_v2 <- roadmap_counts %>% filter(motif == "Motif_down_5") %>% 
  mutate(full_names = fct_reorder(factor(roadmap_annotation[roadmap]), ratio), 
         expected = ifelse(ratio > 1,
                           "More than expected", 
                           "Less than expected")) %>% 
  ggplot(aes(full_names, ratio - 1, color = count)) +
  geom_segment(aes(y = 0,
                   yend = ratio - 1,
                   xend = full_names),
               color = "black") +
  geom_point(size = 6) + 
  scale_y_continuous(labels = function(y) y + 1) + 
  facet_grid(~ motif) +  theme_bw() + 
  scale_color_distiller(palette = "Blues", direction = 1) +
  labs(x = "", y = "", color = "Expectation direction") + coord_flip()

prow_lolli_v2 <- plot_grid(
  motif1_lolli_v2 + theme(legend.position="none"),
  motif2_lolli_v2 + theme(legend.position="none"),
  motif3_lolli_v2 + theme(legend.position="none"),
  align = 'vh',
  hjust = -1,
  nrow = 1
)

legend_b_v3 <- get_legend(
  motif1_lolli_v2 + 
    guides(color = guide_legend(nrow = 1)) +
    theme(legend.position = "bottom")
)

pdf("roadmap_enrichment_v4.pdf", width = 18, height = 6)
plot_grid(prow_lolli_v2, legend_b_v3, ncol = 1, rel_heights = c(1, .1))
dev.off()


motif_names <- c(str_c("Motif_down_", c(1:5)), str_c("Motif_up_", c(1, 7)))
motifsOccurence <- read_tsv("sentinel_with_motifs_v1.tsv")

colnames(motifsOccurence)[1] <- "PeakID"

motifOccurenceGR <- motifsOccurence %T>% {colnames(.)[22:28] <- motif_names} %>% 
  pivot_longer(cols = 22:28, names_to = "motif", values_to = "hit") %>%
  filter(!is.na(motif)) %>%
  mutate(hit = 
           str_replace_all(hit, 
                           pattern = "(\\d+\\([A-Za-z]+,[*+-],\\d+\\.\\d+\\)),", 
                           replacement = "\\1_")) %>% 
  separate(hit, into = c("occ1", "occ2", "occ3", "occ4"), "_") %>% 
  pivot_longer(cols = starts_with("occ"), 
               names_to = "ignore", values_to = "hit") %>% 
  filter(!is.na(hit)) %>% 
  tidyr::extract(col = "hit", 
                 into = c("offset", "sequence", "motif_strand"), 
                 "(\\d+)\\(([A-Za-z]+),([*+-]).+", convert = T) 


motif_files <- list.files("motif_centered_MOSA", full.names = T)

motif_names <- c(str_c("Motif_down_", c(1:5)), str_c("Motif_up_", c(1, 7)))

names(motif_files) <- motif_names

library(heatmaps)
library(BSgenome.Hsapiens.UCSC.hg19)

motifs_gr <- lapply(motif_files, function(x){
  motif <- read_tsv(x, col_names = F)
  makeGRangesFromDataFrame(motif, 
                           seqnames.field = "X2", 
                           start.field = "X3", 
                           end.field = "X4", keep.extra.columns = T)
})

pwms <- lapply(names(motif_matrices), function(x) {
  t(as.matrix(motif_matrices[[x]] * 1000 )) %T>% 
    { rownames(.) <- c("A", "C", "G", "T") } %>% 
    { PFMatrix(name = x, profileMatrix = .) } %>% toPWM()
})

motif_orient <- lapply(1:7, function(x) {
  tmpSites <- searchSeq(pwms[[x]], getSeq(BSgenome.Hsapiens.UCSC.hg19, 
                                          resize(motifs_gr[[x]], 25, 
                                                 fix = "center")))
  ifelse(sapply(tmpSites, function(y) start(y[1])) >= 9, "+", "-")
})

motifs_gr_orient <- lapply(1:7, function(x){
  tmp <- motifs_gr[[x]]
  strand(tmp) <- motif_orient[[x]]
  tmp
})

names(motifs_gr_orient) <- names(motifs_gr)

heatMotifs <- lapply(motifs_gr_orient, function(x){
  seq_motif_7 <- getSeq(BSgenome.Hsapiens.UCSC.hg19, 
                        x)
  PatternHeatmap(seq_motif_7, "CG", coords = c(-150, 150)) })
plotHeatmapList(heatMotifs)

genomic_cg_motif1 <- heatMotifs$Motif_down_1 %>% image() %>% colSums() %>% 
  { ksmooth(coords, ., bandwidth = 5) } %>% as.data.frame() %>% 
  ggplot(aes(x, y)) + 
  geom_rect(aes(xmin=0, xmax=9, ymin=min(y), ymax=max(y)), 
            color="transparent", fill="orange", alpha=0.1) +
  geom_line(size = 1.25, col = "royalblue4") + 
  facet_grid(~ "Motif 1")  +
  labs(x = "", y = "") + theme_bw() + 
  theme(axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm") )

genomic_cg_motif2 <- heatMotifs$Motif_down_2 %>% image() %>% colSums() %>% 
  { ksmooth(coords, ., bandwidth = 5) } %>% as.data.frame() %>% 
  ggplot(aes(x, y)) + 
  geom_rect(aes(xmin=0, xmax=11, ymin=min(y), ymax=max(y)), 
            color="transparent", fill="orange", alpha=0.1) +
  geom_line(size = 1.25, col = "royalblue4") + 
  facet_grid(~ "Motif 2") + 
  labs(x = "", y = "CG dinucleotide density") + theme_bw() + 
  theme(axis.text.x = element_blank(), plot.margin = unit(c(0, 0, 0, 0), "cm") )

genomic_cg_motif3 <- heatMotifs$Motif_down_4 %>% image() %>% colSums() %>% 
  { ksmooth(coords, ., bandwidth = 5) } %>% as.data.frame() %>% 
  ggplot(aes(x, y)) + 
  geom_rect(aes(xmin=0, xmax=9, ymin=min(y), ymax=max(y)), 
            color="transparent", fill="orange", alpha=0.1) +
  geom_line(size = 1.25, col = "royalblue4") + 
  facet_grid(~ "Motif 4") +
  labs(x = "Position relative to motif", 
       y = "") + theme_bw() +  theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
  scale_x_continuous(breaks = seq(-150, 150, 50))

plot_grid(genomic_cg_motif1, genomic_cg_motif2, genomic_cg_motif3, 
          ncol = 1, align = "vh")

all_pvals <- read_tsv("final_motifs_MOSA_stats_all/knownResults.txt")
down_pvals <- read_tsv("final_motifs_MOSA_stats_down/knownResults.txt")
up_pvals <- read_tsv("final_motifs_MOSA_stats_up/knownResults.txt")

discovery_pvals <- inner_join(all_pvals, down_pvals, 
                              by = "Motif Name",
                              suffix = c("_all", "_down")) %>% 
  inner_join(up_pvals, by = "Motif Name") %>% 
  select(`Motif Name`, contains("P-value"),
         "# of Target Sequences with Motif(of 691)")

log_discovery_pvals_mat <- discovery_pvals %>% select(starts_with("Log")) %>%
  as.matrix()

rownames(log_discovery_pvals_mat) <- discovery_pvals$`Motif Name`
colnames(log_discovery_pvals_mat) <- c("All", "Down", "Up")
log_discovery_pvals_mat <- log_discovery_pvals_mat * -1
log_discovery_pvals_mat <- log_discovery_pvals_mat[
  down_pvals$`Motif Name`[order(down_pvals$`P-value`)], 
]

sentinel_motif_number <- discovery_pvals$`# of Target Sequences with Motif(of 691)`
names(sentinel_motif_number) <- discovery_pvals$`Motif Name`
sentinel_motif_number <- sentinel_motif_number[
  down_pvals$`Motif Name`[order(down_pvals$`P-value`)]
]

motif_matrices_files <- c(str_c("homer_new/MOSA_down_arraybg_2/homerResults/motif", c(1:5), ".motif"),
                          str_c("homer_new/MOSA_up_arraybg_2/homerResults/motif", c(1, 7), ".motif"))
names(motif_matrices_files) <- motif_names

motif_matrices <- lapply(motif_matrices_files, read_table2, skip = 1, col_names = F)
foo <- lapply(names(motif_matrices), function(x) t(as.matrix(motif_matrices[[x]] * 1000 )) %T>% { rownames(.) <- c("A", "C", "G", "T") } %>% { PFMatrix(name = x, profileMatrix = .) } %>% toICM())
pdf("motif_1.pdf", width = 8, height = 4); sl1(); dev.off()
pdf("motif_2.pdf", width = 8, height = 4); seqLogo(foo[[2]]); dev.off()
pdf("motif_3.pdf", width = 8, height = 4); seqLogo(foo[[3]]); dev.off()
pdf("motif_4.pdf", width = 8, height = 4); seqLogo(foo[[4]]); dev.off()
pdf("motif_5.pdf", width = 8, height = 4); seqLogo(foo[[5]]); dev.off()
pdf("motif_6.pdf", width = 8, height = 4); seqLogo(foo[[6]]); dev.off()
pdf("motif_7.pdf", width = 8, height = 4); seqLogo(foo[[7]]); dev.off()

pdf("heatmap_enrichment_motifs.pdf", width = 6, height = 8)
ComplexHeatmap::Heatmap(log_discovery_pvals_mat[ ,c(2:3)],
                        cluster_columns = F, cluster_rows = F,
                        col = colorRamp2(breaks = c(0, 2, 30),
                                         c("royalblue4", "white", "red3") ),
                        show_row_names = F, name = "-log(p-value)",
                        right_annotation = rowAnnotation(
                          `Number of sentinels with motif` = 
                            anno_barplot(sentinel_motif_number,
                                         bar_width = .9), 
                          width = unit(2.5, "in")
                        ))
dev.off


all_motif_ranges <- lapply(names(motifs_gr), function(x) {
  tmp <- motifs_gr[[x]]
  tmp$motif <- x
  tmp
}) %>% purrr::reduce(c)
all_motif_ranges$genomic_cpg <- str_count(toupper(getSeq(BSgenome.Hsapiens.UCSC.hg19, all_motif_ranges)), "CG")
all_motif_ranges$array_cpg <- countOverlaps(all_motif_ranges, arrayRanges)
all_motif_ranges$cpg_coverage <- all_motif_ranges$array_cpg / all_motif_ranges$genomic_cpg
plot_grid(~ hist(motif_array_coverage), ~ hist(motif_genomic_coverage), ~ hist(motif_array_coverage / motif_genomic_coverage), ~ plot(motif_array_coverage, motif_genomic_coverage))

#############################################################################################
