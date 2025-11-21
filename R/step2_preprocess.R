library(data.table)
library(ggtree)
library(ggplot2)
library(ape)
library(RColorBrewer)
library(tinyarray)
library(sva)

load('step1.RData')

################################################################################
# merge exp
merge_all_exp_loop <- function(datasets_list, keep_all = FALSE) {
  if (length(datasets_list) < 2) {
    stop("至少需要两个数据框")
  }
  
  result <- datasets_list[[1]]
  
  for (i in 2:length(datasets_list)) {
    result <- merge(result, datasets_list[[i]], by = 0, all = keep_all)
    rownames(result) <- result[, 1]
    result <- result[, -1, drop = FALSE]
  }
  
  return(result)
}

rename_and_prefix <- function(data_list) {
  lapply(names(data_list), function(name) {
    df <- data_list[[name]]
    colnames(df) <- paste0(name, "_", colnames(df))
    return(df)
  }) |> setNames(names(data_list))
}


datasets <- list(
  GSE5281 = exp_GSE5281,
  GSE28146 = exp_GSE28146,
  GSE173954 = exp_GSE173954,
  GSE193438 = exp_GSE193438,
  GSE236562 = exp_GSE236562,
  GSE29378 = exp_GSE29378,
  GSE36980 = exp_GSE36980,
  GSE48350 = exp_GSE48350,
  GSE67333 = exp_GSE67333,
  GSE84422 = exp_GSE84422
)

datasets = rename_and_prefix(datasets)
exp_all = merge_all_exp_loop(datasets)
boxplot(exp_all, las = 2, cex.axis = 0.6)



################################################################################
safe_combine_metadata <- function(sample_list) {
  meta_list <- list()
  
  for(gse in names(sample_list)) {
    df <- sample_list[[gse]]
    
    # 确保geo_accession列存在且正确命名
    if(!"geo_accession" %in% colnames(df)) {
      stop(paste("数据集", gse, "缺少geo_accession列"))
    }
    
    # 添加前缀到样本ID
    df$geo_accession <- paste0(gse, "_", df$geo_accession)
    df$batch <- gse
    
    meta_list[[gse]] <- df
  }
  
  # 合并并检查
  meta_combined <- rbindlist(meta_list, fill = TRUE)
  
  # 验证合并结果
  cat("合并后的疾病分布:\n")
  print(table(meta_combined$batch, meta_combined$disease))
  
  return(meta_combined)
}

sample_list <- list(
  GSE5281 = sample_GSE5281,
  GSE28146 = sample_GSE28146,
  GSE173954 = sample_GSE173954,
  GSE193438 = sample_GSE193438,
  GSE236562 = sample_GSE236562,
  GSE29378 = sample_GSE29378,
  GSE36980 = sample_GSE36980,
  GSE48350 = sample_GSE48350,
  GSE67333 = sample_GSE67333,
  GSE84422 = sample_GSE84422
)

# 创建metadata
meta_combined <- safe_combine_metadata(sample_list)

# 检查样本名称匹配情况
cat("表达矩阵样本数量:", ncol(exp_all), "\n")
cat("metadata样本数量:", nrow(meta_combined), "\n")

# 检查样本名称的前几个
cat("表达矩阵前5个样本:", head(colnames(exp_all), 5), "\n")
cat("metadata前5个样本:", head(meta_combined$geo_accession, 5), "\n")

# 查找交集
common_samples <- intersect(colnames(exp_all), meta_combined$geo_accession)
cat("共同样本数量:", length(common_samples), "\n")

# 创建分组向量
group_list_batch <- setNames(factor(meta_combined$batch), meta_combined$geo_accession)
group_list_disease <- setNames(factor(meta_combined$disease), meta_combined$geo_accession)

# 确保样本顺序一致
common_samples <- intersect(colnames(exp_all), names(group_list_batch))

cat("最终共同样本数量:", length(common_samples), "\n")

# 按相同顺序排列所有数据
exp_all_ordered <- exp_all[, common_samples]
batch_ordered <- group_list_batch[common_samples]
disease_ordered <- group_list_disease[common_samples]


################################################################################
# batch effect
exp_all_scaled = scale(exp_all_ordered)
boxplot(exp_all_scaled, las = 2, cex.axis = 0.6)

draw_pca(exp = exp_all_scaled, group_list = batch_ordered)

cat("批次与疾病分布:\n")
print(table(batch_ordered, disease_ordered))

mod_matrix <- model.matrix(~ disease_ordered)

# ComBat
exp_all_combat <- ComBat(
  dat = exp_all_scaled,
  batch = batch_ordered,
  mod = mod_matrix
)

draw_pca(exp = exp_all_combat, group_list = batch_ordered)

save(exp_all_combat, meta_combined, file = 'step2.RData')
