library(tidyverse)
alltis = list.files('./data/processed/expression/change_w_age/')
genesx = unique(unname(unlist(lapply(alltis,function(x)rownames(readRDS(paste('./data/processed/expression/change_w_age/',x,sep='')))))))

rhomat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'rho']
})
rownames(rhomat) = genesx
saveRDS(rhomat,'./data/processed/expression/rhomat.rds')

padjmat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'padj']
})
rownames(padjmat) = genesx
saveRDS(padjmat,'./data/processed/expression/padjmat.rds')

pmat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'pval']
})
rownames(pmat) = genesx
saveRDS(pmat,'./data/processed/expression/pmat.rds')