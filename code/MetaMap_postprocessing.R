library(stringr)
setwd("CAD out") # the directory where MetaMap output files are saved
outputs = list()
for(f in dir()) {
  textlines = readLines(f)
  CUIs = str_extract(textlines, "C[0-9]{7}")
  CUIs = CUIs[!is.na(CUIs)]
  outputs[[f]] = unique(CUIs)
}
CUIs.all = unique(unlist(outputs))
CUI.occ = matrix(0, nrow = length(CUIs.all), ncol = length(outputs))
colnames(CUI.occ) = dir()
row.names(CUI.occ) = CUIs.all
for(f in dir()) {
  CUI.occ[outputs[[f]],f] = 1
}
writeLines(CUIs.all[rowSums(CUI.occ) > length(outputs)/2], "CUIs.txt") # majority voting and output
