library(data.table)
library(RColorBrewer)
library(tinyarray)
library(sva)

# 加载数据
load('step1.RData')

################################################################################
draw_batch_boxplot <- function(exp_matrix, batch_vector) {
  exp_matrix <- as.matrix(exp_matrix)
  
  if (!is.factor(batch_vector)) batch_vector <- factor(batch_vector)
  if (length(batch_vector) != ncol(exp_matrix)) {
    stop("batch_vector长度必须与exp_matrix列数相同")
  }
  
  n_batch <- length(levels(batch_vector))
  col_base <- brewer.pal(min(max(n_batch, 3), 9), "Set1")
  col_use <- if(n_batch > length(col_base)) colorRampPalette(col_base)(n_batch) else col_base[1:n_batch]
  col_batch <- setNames(col_use, levels(batch_vector))
  col_sample <- col_batch[batch_vector]
  
  boxplot(exp_matrix, col = col_sample, outline = FALSE, las = 2, cex.axis = 0.6)
  legend("topright", legend = names(col_batch), fill = col_batch, cex = 0.6, inset = 0.05)
}

################################################################################
merge_exp <- function(datasets_list, keep_all = FALSE) {
  if (length(datasets_list) < 2) stop("至少需要两个数据框")
  result <- datasets_list[[1]]
  for (i in 2:length(datasets_list)) {
    result <- merge(result, datasets_list[[i]], by = 0, all = keep_all)
    rownames(result) <- result[, 1]
    result <- result[, -1, drop = FALSE]
  }
  return(result)
}

add_prefix <- function(data_list) {
  lapply(names(data_list), function(n) {
    df <- data_list[[n]]
    colnames(df) <- paste0(n, "_", colnames(df))
    df
  }) |> setNames(names(data_list))
}

################################################################################
# 处理表达矩阵
datasets <- list(
  GSE5281 = exp_GSE5281, GSE28146 = exp_GSE28146, GSE173954 = exp_GSE173954,
  GSE29378 = exp_GSE29378, GSE36980 = exp_GSE36980, GSE48350 = exp_GSE48350,
  GSE84422 = exp_GSE84422
)
datasets <- add_prefix(datasets)
exp_all <- merge_exp(datasets)

################################################################################
# Metadata合并
merge_meta <- function(sample_list) {
  meta_list <- lapply(names(sample_list), function(gse) {
    df <- sample_list[[gse]]
    if(!"geo_accession" %in% colnames(df)) stop(paste(gse, "缺少geo_accession列"))
    df$geo_accession <- paste0(gse, "_", df$geo_accession)
    df$batch <- gse
    df
  })
  meta_combined <- rbindlist(meta_list, fill = TRUE)
  cat("合并后的疾病分布:\n")
  print(table(meta_combined$batch, meta_combined$disease))
  return(meta_combined)
}

sample_list <- list(
  GSE5281 = sample_GSE5281, GSE28146 = sample_GSE28146, GSE173954 = sample_GSE173954,
  GSE29378 = sample_GSE29378, GSE36980 = sample_GSE36980, GSE48350 = sample_GSE48350,
  GSE84422 = sample_GSE84422
)
meta_combined <- merge_meta(sample_list)
meta_combined$disease = ifelse(meta_combined$disease, 'AD', 'nonAD')


################################################################################
# 排序
common_samples <- intersect(colnames(exp_all), meta_combined$geo_accession)

# 创建带名称的分组向量
group_list_batch <- setNames(factor(meta_combined$batch), meta_combined$geo_accession)
group_list_disease <- setNames(factor(meta_combined$disease), meta_combined$geo_accession)

# 显式按共同样本排序
exp_all <- exp_all[, common_samples]
group_list_batch <- group_list_batch[common_samples]
group_list_disease <- group_list_disease[common_samples]


draw_batch_boxplot(exp_all, group_list_batch)
draw_pca(exp_all, group_list_batch)
draw_pca(exp_all, group_list_disease)


################################################################################
# 批次效应校正
cat("\n=== 批次与疾病分布 ===\n")
print(table(group_list_batch, group_list_disease))

# ComBat校正
exp_all_combat <- ComBat(
  dat = exp_all,
  batch = group_list_batch,
  mod = group_list_disease
)

draw_batch_boxplot(exp_all_combat, group_list_batch)
draw_pca(exp_all_combat, group_list_batch)
draw_pca(exp_all_combat, group_list_disease)


################################################################################
# save
save(exp_all_combat, meta_combined, file = 'step2.RData')