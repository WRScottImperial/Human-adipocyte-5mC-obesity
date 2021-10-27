#############################################################################################

##### Identification of sentinel DNA methylation sites associated with expression.
##### In combined subcutaneous and visceral adipocyte samples, using mixed effects controlling for sample relatedness.

#############################################################################################

##### load packages and set options

library("DESeq2")
library('biomaRt')
library('variancePartition')
library('edgeR')
library('BiocParallel')

options(stringsAsFactors=F)

# Load methylation betas, control probe PCs and annotation.

load("...beta_QN_rep.RData")
phe.cg=read.table("...phe.rep.txt", header=T, sep='\t')
load("...ctrlpca_rep.RData")
phe.cg=merge(phe.cg, ctrlprobes.scores, by.y='row.names', by.x='ID')
scr=1-colSums(is.na(beta))/nrow(beta)
phe.cg = phe.cg[phe.cg$Array_ID %in% names(scr[scr>0.98]),] 

# Load RNA seq counts, surrogate variables and annotation.

counts = round(read.csv(".../RawCounts.combinedFlowcells.csv",row.names=1))
phe.gn=read.table("...phe.rep.with.RNA.qc.txt", header=T, sep='\t')
phe.gn=phe.gn[which(phe.gn$Lib_Prep_Batch != "Failed"),]
rownames(phe.gn)=phe.gn$Sample.ID
SVAs = read.csv(".../vst_count_SVs.csv")
phe.gn=merge(phe.gn, SVAs, by.y='row.names', by.x='ID')

# Merge datasets and remove failed RNA seq samples/extreme outliers

phe.all = merge(phe.gn,phe.cg,by=c('Sample.ID','CaseControl','Group','Age','Sex','Sample.Type','Ethnicity'))  
rownames(phe.all)=phe.all$Sample.ID

outliers = c("P116B.SA","P176B.VA","P38B.VA","P10B.SA","P76B.SA") # 
phe.sel = phe.all[!phe.all$Sample.ID %in% outliers,]

beta = beta[,phe.all$Sample.ID]
counts = counts[,phe.all$Sample.ID]

# load sentinels 

cpgs=read.csv(".../Obese_lean_SA_Combined_sentinels.csv")
cpgs.assoc=cpgs$CG

# load Illumina manifest and transcript annotations

load(".../HumanMethylation450_15017482_v.1.1.RData")
anno=anno[,c("IlmnID","CHR","MAPINFO")]
colnames(anno)=c("ID","chr","pos")
load(".../Genemap_Hg19_restricted_to_RNAseq_counts.RData")

# Run models of 5mC ~ cis-gene expression using variancePartition: DREAM to allow for mixed-effects modeling.

filter=rowSums(counts >= 5) >= 0.2*ncol(counts)
filteredCounts = counts[filter,]

ds.all = DGEList(filteredCounts)
ds.all = calcNormFactors(ds.all)

# Loop for each sentinel 

for(cg in 1:length(cpgs.assoc)) {
  
  cpg = cpgs.assoc[cg]
  
  # 2. Run with Limma in SA and VA combined using Dream mixed effect model 
  
  phe.sel$sent = beta[cpgs.assoc[cg],]
  phe.all.sel = phe.sel[!is.na(phe.sel$sent),]
  ds.all.sel = ds.all[,phe.all.sel$Sample.ID]
  
  param = SnowParam(4, "SOCK", progressbar=TRUE)
  register(param)
  
  formula = ~ sent + Age + Sex + Ethnicity + RIN + SV1_SAVA + SV2_SAVA + PC1_cp + PC2_cp + (1|Participant.ID) #
  
  voom.dream.ds.all = voomWithDreamWeights(ds.all.sel, formula, phe.all.sel)
  fit.dream.ds.all = dream(voom.dream.ds.all, formula, phe.all.sel)
  
  resu = topTable(fit.dream.ds.all, coef="sent", n = Inf)
  
  # save results for each sentinel
  
  write.csv(resu, file=file.path(dir,paste(cpg,"_cisgene_association_dream.csv",sep=""),row.names=T))
  
}

# EOF





