library(GEOquery)
library(AnnoProbe)
library(oligo)
library(stringr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(sva)
library(tinyarray)
library(data.table)
library(illuminaio)
library(AnnotationDbi)
library(org.Hs.eg.db)


filterEM2 <- function(probes_expr,probe2gene,method="mean"){
  colnames(probe2gene) <- c("probeid","symbol")
  probe2gene$probeid=as.character(probe2gene$probeid)
  probe2gene$symbol=trimws(probe2gene$symbol)
  # head(probe2gene)
  
  message(paste0('input expression matrix is ',nrow(probes_expr),' rows(genes or probes) and ',ncol(probes_expr),' columns(samples).\n'))
  message(paste0('input probe2gene is ',nrow(probe2gene),' rows(genes or probes)\n'))
  
  probe2gene=na.omit(probe2gene)
  # if one probe mapped to many genes, we will only keep one randomly.
  probe2gene=probe2gene[!duplicated(probe2gene$probeid),]
  # 这个地方是有问题的，随机挑选一个注释进行后续分析。
  probe2gene = probe2gene[probe2gene$probeid %in% rownames(probes_expr),]
  
  message(paste0('after remove NA or useless probes for probe2gene, ',nrow(probe2gene),' rows(genes or probes) left\n'))
  
  #probes_expr <- exprs(eSet);dim(probes_expr)
  probes_expr <- as.data.frame(probes_expr)
  genes_expr <- tibble::rownames_to_column(probes_expr,var="probeid") %>%
    merge(probe2gene,.,by="probeid")
  genes_expr <- genes_expr[-1]
  # probes_expr[1:4,1:4]
  
  # remove duplicates sympol:method=mead,also median, max ,min
  message("remove duplicates sympols, it will take a while")
  genes_expr <- aggregate(x=genes_expr[,2:ncol(genes_expr)],by=list(genes_expr$symbol),FUN=method,na.rm=T)
  genes_expr <- tibble::column_to_rownames(genes_expr,var = "Group.1")
  # genes_expr[1:4,1:4]
  message(paste0('output expression matrix is ',nrow(genes_expr),' rows(genes or probes) and ',ncol(genes_expr),' columns(samples).'))
  # probes_expr['AGAP6',]
  return(genes_expr)
}

read_CEL <- function(filedir){
  celFiles <- list.celfiles(filedir,listGzipped=T,
                            full.name=TRUE)
  exon_data = oligo::read.celfiles(celFiles)
  exon_data_rma = oligo::rma(exon_data)
  exp_probe = Biobase::exprs(exon_data_rma)
  colnames(exp_probe) = sub("\\.CEL\\.gz$", "", colnames(exp_probe))
  exp_probe[is.na(exp_probe)] = 0
  return(exp_probe)
}

get_exp <- function(filedir){
  geo = getGEO(filename = filedir, getGPL = FALSE)
  exp = exprs(geo)
  if (max(exp) > 50){
    exp[is.na(exp)] = 0
    exp = log2(exp)
  }
  return(exp)
}

get_probe2symbol <- function(filename, flag, number = 1) {
  geo <- getGEO(filename = filename, getGPL = FALSE)
  platform_ids <- names(geo@gpls)
  if (length(platform_ids) == 0) {
    stop("未找到任何GPL平台信息")
  }
  gpl_id <- platform_ids[number]
  message("检测到平台ID: ", gpl_id)
  GPL_table <- geo@gpls[[gpl_id]]@dataTable@table
  if (flag == 1){
    GPL_table = GPL_table[,c(1,11)]
    setDT(GPL_table)
    setnames(GPL_table, c("ProbeID", "Gene.Symbol"))
    probe2symbol <- GPL_table[
      !Gene.Symbol %in% c("", NA, "---", "--- /// ---")
    ][,
      .(Symbol = trimws(strsplit(Gene.Symbol, "///", fixed = TRUE)[[1]][1])), 
      by = ProbeID
    ]
  }
  if (flag == 2) {
    GPL_table = GPL_table[, c(1, 10)]
    setDT(GPL_table)
    setnames(GPL_table, c("ProbeID", "Gene.Symbol"))
    probe2symbol <- GPL_table[
      !is.na(Gene.Symbol) & Gene.Symbol != "" & !grepl("^---", Gene.Symbol)
    ][,
      .(Symbol = {
        # 取第一个基因块，然后提取其中的第二个字段
        first_block <- strsplit(Gene.Symbol, " /// ", fixed = TRUE)[[1]][1]
        trimws(strsplit(first_block, " // ", fixed = TRUE)[[1]][2])
      }), 
      by = ProbeID
    ]
  }
  if (flag == 3){
    GPL_table = GPL_table[, c(1, 7)]
    setDT(GPL_table)
    setnames(GPL_table, c("ProbeID", "Gene.Symbol"))
    probe2symbol = GPL_table
  }
  if (flag == 4){
    GPL_table = GPL_table[, c(1,8)]
    setDT(GPL_table)
    setnames(GPL_table, c("ProbeID", "Gene.Symbol"))
    probe2symbol <- GPL_table[
      !is.na(Gene.Symbol) & Gene.Symbol != "" & !grepl("^---", Gene.Symbol)
    ][,
      .(Symbol = {
        # 取第一个基因块，然后提取其中的第二个字段
        first_block <- strsplit(Gene.Symbol, " /// ", fixed = TRUE)[[1]][1]
        trimws(strsplit(first_block, " // ", fixed = TRUE)[[1]][2])
      }), 
      by = ProbeID
    ]
  }
  return(probe2symbol)
}



################################################################################
# GSE5281
## exp
exp_GSE5281 = get_exp('geo_data/GSE5281/GSE5281_series_matrix.txt.gz')

## GPL
probe2symbol = get_probe2symbol('./geo_data/GSE5281/GSE5281_family.soft.gz', 1)
exp_GSE5281 = filterEM2(exp_GSE5281, probe2symbol)
exp_GSE5281 = exp_GSE5281[,-1]

## group
geo = getGEO(filename = 'geo_data/GSE5281/GSE5281_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
p = p[,c(1,2)]
sample_GSE5281 = p[grepl("HIP", p$title), ]
sample_GSE5281$'disease' = grepl("affected", sample_GSE5281$title)
sample_GSE5281 = sample_GSE5281[,-1]

exp_GSE5281 = exp_GSE5281[, c(sample_GSE5281$geo_accession)]


################################################################################
# GSE28146
## exp
exp_GSE28146 = get_exp('geo_data/GSE28146/GSE28146_series_matrix.txt.gz')

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE28146/GSE28146_family.soft.gz', 1)
exp_GSE28146 = filterEM2(exp_GSE28146, probe2symbol)
exp_GSE28146 = exp_GSE28146[,-1]

## group
geo = getGEO(filename = 'geo_data/GSE28146/GSE28146_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE28146 = p[,c(1,2)]
sample_GSE28146$'disease' = !grepl("Control", sample_GSE28146$title)
sample_GSE28146 = sample_GSE28146[,-1]

exp_GSE28146 = exp_GSE28146[, c(sample_GSE28146$geo_accession)]


################################################################################
# GSE29378
## exp
exp_GSE29378 = get_exp('geo_data/GSE29378/GSE29378_series_matrix.txt.gz')

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE29378/GSE29378_family.soft.gz', 3)
exp_GSE29378 = filterEM2(exp_GSE29378, probe2symbol)
exp_GSE29378 = exp_GSE29378[,-1]

## group
geo = getGEO(filename = 'geo_data/GSE29378/GSE29378_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE29378 = p[,c(1,2)]
sample_GSE29378$'disease' = grepl("AD", sample_GSE29378$title)
sample_GSE29378 = sample_GSE29378[,-1]

exp_GSE29378 = exp_GSE29378[, c(sample_GSE29378$geo_accession)]


################################################################################
# GSE36980
## exp
exp_GSE36980 = get_exp('geo_data/GSE36980/GSE36980_series_matrix.txt.gz')

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE36980/GSE36980_family.soft.gz', 2)
exp_GSE36980 = filterEM2(exp_GSE36980, probe2symbol)
exp_GSE36980 = exp_GSE36980[,-1]

## group
geo = getGEO(filename = 'geo_data/GSE36980/GSE36980_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
p = p[,c(1,2)]
sample_GSE36980 = p[grepl('HI', p$title),]
sample_GSE36980$'disease' = !grepl("non-AD", sample_GSE36980$title)
sample_GSE36980 = sample_GSE36980[, -1]

exp_GSE36980 = exp_GSE36980[, c(sample_GSE36980$geo_accession)]


################################################################################
# GSE48350
## exp
exp_GSE48350 = read_CEL('geo_data/GSE48350/GSE48350_RAW/')

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE48350/GSE48350_family.soft.gz', 1)
exp_GSE48350 = filterEM2(exp_GSE48350, probe2symbol)
exp_GSE48350 = exp_GSE48350[,-1]

## group
geo = getGEO(filename = 'geo_data/GSE48350/GSE48350_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
p = p[,c(1,2)]
sample_GSE48350 = p[grepl('Hippocampus', p$title) | grepl('hippocampus', p$title),]
sample_GSE48350$'disease' = grepl("AD", sample_GSE48350$title)
sample_GSE48350 = sample_GSE48350[, -1]

exp_GSE48350 = exp_GSE48350[, c(sample_GSE48350$geo_accession)]


################################################################################
# GSE84422
## exp
exp_GSE84422 = get_exp('geo_data/GSE84422/GPL96/GSE84422-GPL96_series_matrix.txt.gz')

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE84422/GSE84422_family.soft.gz', 1)
exp_GSE84422 = filterEM2(exp_GSE84422, probe2symbol)
exp_GSE84422 = exp_GSE84422[, -1]

## group
geo = getGEO(filename = 'geo_data/GSE84422/GPL96/GSE84422-GPL96_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE84422 = p[p$`brain region:ch1` == 'Hippocampus' & p$`neuropathological category:ch1` != 'Probable AD' & p$`neuropathological category:ch1` != 'Possible AD', ]
sample_GSE84422 = sample_GSE84422[, c(2, 50)]
sample_GSE84422$'disease' = (sample_GSE84422$`neuropathological category:ch1`=='definite AD')
sample_GSE84422 = sample_GSE84422[,-2]

exp_GSE84422 = exp_GSE84422[, c(sample_GSE84422$geo_accession)]


################################################################################
# GSE173954
## exp
exp_GSE173954 = read_CEL('geo_data/GSE173954/GSE173954_RAW/')
colnames(exp_GSE173954) <- sub("^(GSM\\d+).*", "\\1", colnames(exp_GSE173954))

## GPL
probe2symbol = get_probe2symbol('geo_data/GSE173954/GSE173954_family.soft.gz', 4)
exp_GSE173954 = filterEM2(exp_GSE173954, probe2symbol)
exp_GSE173954 = exp_GSE173954[, -1]

## group
geo = getGEO(filename = 'geo_data/GSE173954/GSE173954_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE173954 = p[,c(1,2)]
sample_GSE173954$'disease' = !grepl('non-Alzheimer\'s Disease', sample_GSE173954$title)
sample_GSE173954 = sample_GSE173954[,-1]

exp_GSE173954 = exp_GSE173954[, c(sample_GSE173954$geo_accession)]


################################################################################
# GSE236562
## exp
exp_GSE236562 = read.delim('geo_data/GSE236562/GSE236562_raw_counts_GRCh38.p13_NCBI.tsv', row.names = 1)
exp_GSE236562 = log2(exp_GSE236562)
exp_GSE236562[exp_GSE236562 == '-Inf'] = 0

## GPL
probe2symbol = read.delim('geo_data/GSE236562/Human.GRCh38.p13.annot.tsv')
probe2symbol = probe2symbol[, c(1, 2)]
exp_GSE236562 = filterEM2(exp_GSE236562, probe2symbol)

## group 
geo = getGEO(filename = 'geo_data/GSE236562/GSE236562_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE236562 = p[!grepl('COVID', p$title) & grepl('HP', p$title), c(1,2)]
sample_GSE236562$'disease' = grepl('ADpos', sample_GSE236562$title)
sample_GSE236562 = sample_GSE236562[,-1]

exp_GSE236562 = exp_GSE236562[, c(sample_GSE236562$geo_accession)]


################################################################################
# GSE67333
## exp
exp_GSE67333 = read.delim('geo_data/GSE67333/GSE67333_raw_counts_GRCh38.p13_NCBI.tsv', row.names = 1)
exp_GSE67333 = log2(exp_GSE67333)
exp_GSE67333[exp_GSE67333 == -Inf] = 0


## GPL
probe2symbol = read.delim('geo_data/GSE67333/Human.GRCh38.p13.annot.tsv')
probe2symbol = probe2symbol[, c(1, 2)]
exp_GSE67333 = filterEM2(exp_GSE67333, probe2symbol)

## group
geo = getGEO(filename = 'geo_data/GSE67333/GSE67333_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE67333 = p[, c(1,2)]
sample_GSE67333$'disease' = grepl('AD', sample_GSE67333$title)
sample_GSE67333 = sample_GSE67333[,-1]

exp_GSE67333 = exp_GSE67333[, c(sample_GSE67333$geo_accession)]


################################################################################
# GSE193438
## exp
exp_GSE193438 = read.delim('geo_data/GSE193438/GSE193438_raw_counts_GRCh38.p13_NCBI.tsv', row.names = 1)
exp_GSE193438 = log2(exp_GSE193438)
exp_GSE193438[exp_GSE193438 == -Inf] = 0

## GPL
probe2symbol = read.delim('geo_data/GSE193438/Human.GRCh38.p13.annot.tsv')
probe2symbol = probe2symbol[, c(1, 2)]
exp_GSE193438 = filterEM2(exp_GSE193438, probe2symbol)

## group
geo = getGEO(filename = 'geo_data/GSE193438/GSE193438_series_matrix.txt.gz', getGPL = FALSE)
p = phenoData(geo)
p = p@data
sample_GSE193438 = p[p$source_name_ch1 == 'Hippocampus' & p$`disease:ch1` != 'LBD', c(2,38)]
sample_GSE193438$'disease' = grepl('AD', sample_GSE193438$`disease:ch1`)
sample_GSE193438 = sample_GSE193438[,-2]

exp_GSE193438 = exp_GSE193438[, c(sample_GSE193438$geo_accession)]


################################################################################
# save
remove(geo)
remove(p)
remove(probe2symbol)
remove(filterEM2)
remove(get_exp)
remove(get_probe2symbol)
remove(read_CEL)



save.image(file = "step1.RData")
