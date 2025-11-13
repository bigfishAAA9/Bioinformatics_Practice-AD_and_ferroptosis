library(GEOquery)

gse_number = c('GSE173955', 'GSE184942', 'GSE263319', 'GSE84422', 'GSE48350', 'GSE5281')
dest_dir <- "./geo_data"

geo_data = getGEO('GSE173954', destdir = dest_dir, getGPL = TRUE)
