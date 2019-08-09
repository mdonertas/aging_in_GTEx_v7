library(tidyverse)
library(RColorBrewer)
library(ggpubr)
alltis = list.files('./data/processed/expression/change_w_age/')
genesx = unique(unname(unlist(lapply(alltis,function(x)rownames(readRDS(paste('./data/processed/expression/change_w_age/',x,sep='')))))))

rhomat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'rho']
})
rownames(rhomat) = genesx
colnames(rhomat) = gsub('.rds','',colnames(rhomat))
saveRDS(rhomat,'./data/processed/expression/rhomat.rds')

padjmat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'padj']
})
rownames(padjmat) = genesx
colnames(padjmat) = gsub('.rds','',colnames(padjmat))
saveRDS(padjmat,'./data/processed/expression/padjmat.rds')

pmat = sapply(alltis,function(x){
  ch = readRDS(paste('./data/processed/expression/change_w_age/',x,sep=''))
  ch[genesx,'pval']
})
rownames(pmat) = genesx
colnames(pmat) = gsub('.rds','',colnames(pmat))
saveRDS(pmat,'./data/processed/expression/pmat.rds')

rhocor = cor(rhomat,use = 'pairwise',method = 'spearman')
annotfr = data.frame(tissue = sapply(strsplit(colnames(rhocor),'-'),function(x)x[1]))
rownames(annotfr) = colnames(rhocor)
pheatmap::pheatmap(rhocor, breaks = seq(-1,1,length.out = 101),color=colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(100),cellwidth = 10,cellheight = 10, annotation_row = annotfr, annotation_col = annotfr, cutree_rows = 6,cutree_cols = 6,filename = './results/correlations_all.pdf')

srhomat = rhomat
srhomat[padjmat>0.1] = NA
srhocor = cor(srhomat,use = 'pairwise',method = 'spearman')
annotfr = data.frame(tissue = sapply(strsplit(colnames(srhocor),'-'),function(x)x[1]))
rownames(annotfr) = colnames(srhocor)
mytissues = names(which(rowMeans(is.na(srhocor))<1))
srhocor = srhocor[mytissues,mytissues]
pheatmap::pheatmap(srhocor, breaks = seq(-1,1,length.out = 101),color=colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(100),cellwidth = 10,cellheight = 10, annotation_row = annotfr, annotation_col = annotfr, cutree_rows = 6,cutree_cols = 6,filename = './results/correlations_signif.pdf')

completerho = rhomat[complete.cases(rhomat),]

completerho_sc = apply(completerho,2,scale)

pcxx = prcomp(t(completerho_sc),scale. = T)
pcx = pcxx$x

as.data.frame(pcx) %>%
  mutate(min_tissue = rownames(pcx)) %>% 
  mutate(tissue = sapply(strsplit(min_tissue,'-'),function(x)x[1])) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue)) +
  geom_point(size = 3) +
  theme_pubr()


completerho = rhomat[complete.cases(rhomat),]
scompleterho = completerho[rownames(completerho)%in%names(which(rowSums(padjmat<=0.1,na.rm=T)>=1)),]
dim(scompleterho)

scompleterho_sc = apply(scompleterho,2,scale)

pcxx = prcomp(t(scompleterho_sc),scale. = T)
pcx = pcxx$x

as.data.frame(pcx) %>%
  mutate(min_tissue = rownames(pcx)) %>% 
  mutate(tissue = sapply(strsplit(min_tissue,'-'),function(x)x[1])) %>%
  ggplot(aes(x = PC1, y= PC2, color = tissue)) +
  geom_point(size = 3) +
  theme_pubr()
