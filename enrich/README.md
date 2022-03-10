# Select significant tissues by p-value and z-score

```R
library(dplyr)
library(tidyr)

f_pval = 'roadmap_bed-overlap-1_32-permn_1000-pval.tsv'
f_zscore = 'roadmap_bed-overlap-1_32-permn_1000-zscore.tsv'
f_overlap = 'roadmap_bed-overlap-1_32-permn_1000-overlap.tsv'
f_meta = '../../../Postgwas_v4/db_gwas/roadmap_meta.tsv'
out = 'roadmap_bed-overlap-1_32-summary-pval0.05'

pval = read.delim(f_pval)
zscore = read.delim(f_zscore)
num = read.delim(f_overlap)
meta = read.delim(f_meta)

pvalL = gather(pval, EID, pval, E003:E129)
zscoreL = gather(zscore, EID, zscore, E003:E129)
numL = gather(num, EID, num, E003:E129)

fdr = p.adjust(pvalL$pval)
pvalL$fdr = fdr # Nothing significant!

prom = c('1_TssA','2_TssFlnk','3_TssFlnkU','4_TssFlnkD')
enh = c('7_EnhG1','8_EnhG2','9_EnhA1','10_EnhA2','11_EnhWk')
prom_enh = c(prom,enh)

pvalL_sub = subset(pvalL, Status %in% prom_enh & pval < 0.05)
zscoreL_sub = subset(zscoreL, Status %in% prom_enh & zscore > 0)
numL_sub = subset(numL, Status %in% prom_enh)
perm_list = list(pvalL_sub,zscoreL_sub,numL_sub)

# Merge list with intersect
perm_merge = Reduce(function(x,y) merge(x,y, by=c('Status','EID')), perm_list)
perm_prom = subset(perm_merge, Status %in% prom)
perm_enh = subset(perm_merge, Status %in% enh)

# Summary & Merge prom and enh
perm_prom_summ = data.frame(
    EID = perm_prom$EID,
    Promoter = paste0(
        perm_prom$Status,
        " (p=",formatC(perm_prom$pval, ormat="e",digits=0)
        ,", n=",perm_prom$num,")"))
perm_enh_summ = data.frame(
    EID = perm_enh$EID,
    Enhancer = paste0(
        perm_enh$Status,
        " (p=",formatC(perm_enh$pval,format="e",digits=0)
        ,", n=",perm_enh$num,")"))
perm_prom_enh = merge(perm_prom_summ, perm_enh_summ, by='EID', all=T)

# Annotate by roadmap_meta
perm_merge_annot = merge(perm_merge, meta[,c('EID','ANATOMY','EDACC_NAME')], by='EID', all.x=T) %>% select('ANATOMY','EDACC_NAME','EID','Status','pval','fdr','zscore','num')
perm_summ = merge(perm_prom_enh, meta[,c('EID','ANATOMY','EDACC_NAME')], by='EID', all.x=T) %>% select('ANATOMY','EDACC_NAME','EID','Promoter','Enhancer')

# Save as a TSV file
write.table(perm_merge_annot, paste0(out,'_long.tsv'), sep='\t',row.names=F,quote=F)
write.table(perm_summ,p aste0(out,'.tsv'), sep='\t',row.names=F,quote=F)
```



# Wilcoxon test for tissue groups



## Calculate Wilcoxon Rank Sum Test

```R
library(dplyr)
library(tidyr)

test_df = read.delim('roadmap_highprob_snp_count.tsv')
groups = test_df$sort_simplified %>% unique
status = c("sum","promoter","X1_TssA","X2_TssFlnk","X3_TssFlnkU","X4_TssFlnkD","enhancer","active_enh","X7_EnhG1","X8_EnhG2","X9_EnhA1","X10_EnhA2","X11_EnhWk")

tissue_p_li = lapply(c(1:length(groups)), function(j) {
    test = test_df %>% mutate(group = ifelse(sort_simplified==groups[j],groups[j],"OTHER"))
    test_sub = test %>% filter(group == groups[j])
    eids = test_sub %>% select('EID') %>% unlist

    p_li = lapply(c(1:length(status)), function(i) {
        t = as.character(status[i])
        alter = test_sub %>% select(t) %>% unlist
        #subset1 = unlist(subset1[t])
        null = test %>% filter(group == "OTHER") %>% select(t) %>% unlist
        res <- wilcox.test(alter,null,alternative = "greater")
        c(t, res$p.value %>% as.numeric)
    })
  
    data.frame(
        Group = groups[j],
        EID_N = length(eids),
        EIDs = paste0(eids, collapse=','),
        Promoter_Enhancer = p_li[[1]][2],
        Promoter = p_li[[2]][2],
        TssA = p_li[[3]][2],
        TssFlnk = p_li[[4]][2],
        Tss_FlnkU = p_li[[5]][2],
        Tss_FlnkD = p_li[[6]][2],
        Enhancer = p_li[[7]][2],
        Active_enhancer = p_li[[8]][2],
        EnhA1 = p_li[[11]][2],
        EnhA2 = p_li[[12]][2],
        EnhG1 = p_li[[9]][2],
        EnhG2 = p_li[[10]][2],
        EnhWk = p_li[[13]][2]
    )
})

tissue_p = data.table::rbindlist(tissue_p_li) %>% arrange(Promoter_Enhancer)
write.csv(tissue_p,"roadmap_highprob_snp_pvals.csv")
```

## Draw plots by tissues

```R
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)

#target = 'HEART'
target = c('SKIN','LUNG','BONE','FETAL','FAT','SEXUAL')
f_snp_count = 'roadmap_highprob_snp_count.tsv'

dat = read.delim(f_snp_count)
#table(dat$sort_simplified)

# Data preparation
dat_long = dat %>% select(EID,sort_simplified,sort,full_name,sum,promoter,enhancer,active_enh,X7_EnhG1:X4_TssFlnkD) %>% gather(Status, SNP_n, sum:X4_TssFlnkD)
dat_long$Tissue = 'Others'
for(i in 1:length(target)) {
    dat_long$Tissue[dat_long$sort_simplified==target[i]] = target[i]
}
dat_long$Tissue = factor(dat_long$Tissue, levels=c('Others',target))
dat_long$Status = factor(dat_long$Status, levels=c("sum","promoter","X1_TssA","X2_TssFlnk","X3_TssFlnkU","X4_TssFlnkD","enhancer","active_enh","X9_EnhA1","X10_EnhA2","X7_EnhG1","X8_EnhG2","X11_EnhWk"))

# Draw scatter plot
my_col = c('grey50',brewer.pal(length(target),name="Set3"))
ggplot(dat_long, aes(Status, SNP_n))+
    geom_jitter(aes(colour=Tissue),size=0.5,alpha=0.5)+
    geom_boxplot(aes(fill=Tissue), outliers.shape=NA)+
    theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))+
    scale_color_manual(values=my_col)+
    scale_fill_manual(values=my_col)
targets = paste0(target, collapse=',')
ggsave(paste0('roadmap_jitter_',targets,'.svg'), width=7,height=5, unit='in')
```

