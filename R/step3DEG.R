library(limma)
library(ggvenn)
library(ggVolcano)

load('step2.RData')

################################################################################
# meat
meta_combined <- as.data.frame(meta_combined)
rownames(meta_combined) <- meta_combined[,1]
meta_combined <- meta_combined[,-1]
meta_combined <- meta_combined[colnames(exp_all_combat), ]
meta_combined$disease <- factor(meta_combined$disease)


################################################################################
# DEG
design <- model.matrix(~ disease, data = meta_combined)
print(colnames(design))

fit <- lmFit(exp_all_combat, design)
fit <- eBayes(fit)

coef_to_use <- colnames(design)[2]

deg.results <- topTable(
  fit,
  number = Inf,
  coef = coef_to_use,
  sort.by = "P"
)

DEG <- deg.results[deg.results$adj.P.Val < 0.05 & abs(deg.results$logFC) > 0.2, ]

write.table(deg.results, "results.txt", sep="\t", quote=F)
write.table(DEG, "DEG.txt", sep="\t", quote=F)

################################################################################
# volcano
volcano_data2 = tibble::rownames_to_column(deg.results, "row")

data <- add_regulate(volcano_data2, log2FC_name = "logFC",
                     fdr_name = "adj.P.Val",log2FC = 0.2, fdr = 0.05)

ggvolcano(data, x = "log2FoldChange", y = "padj",
          label = "row", label_number = 10, output = FALSE, FDR_cut  = 0.05,  pointSize = 2,
          log2FC_cut  = 0.2, )


################################################################################
# venn
FRGs = read.csv('D:/Desktop/Bioinformatics_Practice_AD_and_ferroptosis/2step/ferroptosis_driver.csv')
venn_list <- list(
  DEGs = rownames(DEG),
  FRGs = FRGs$symbol
)

ggvenn(
  venn_list,
  fill_color = c("#FF9999", "#66B2FF"),
  fill_alpha = 0.7,
  stroke_color = "white",
  stroke_size = 2,
  text_size = 5,
  set_name_size = 6,
  set_name_color = c("#FF9999", "#66B2FF")
) + 
  theme_void()


intersect_genes = intersect(rownames(DEG), FRGs$symbol)
write.table(intersect_genes, "intersect_genes.txt", sep="\t", quote=F)
