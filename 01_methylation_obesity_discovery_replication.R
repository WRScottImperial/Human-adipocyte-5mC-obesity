#############################################################################################

##### Identification of DNA methylation sites associated with obesity
##### Separately in i. subcutaneous adipocytes and ii. visceral adipocytes

#############################################################################################

#############################################################################################

##### Overview:

### 1. Run tissue specific discovery association analyses (for each CG on obesity status) 
### 2. Run tissue specific replication association analyses (for each CG on obesity status) 
### 3. Use Inverse variance weighted meta-analysis to combine results
### 4. Identify significant CG sites passing significance threshold, and sentinel site for each locus

#############################################################################################

### 1. Discovery Analysis (450K array)

# load normalised betas, control probe principal components, and annotation 

load("...beta_QN_disc.RData")
phe=read.table("...phe.disc.txt", header=T, sep='\t')
load("...ctrlpca_disc.RData")
phe.all=merge(phe.all, ctrlprobes.scores, by.y='row.names', by.x='ID')

# Remove failed sample with to extremely low call-rate (~14%)

scr=1-colSums(is.na(beta))/nrow(beta)
phe.all = phe.all[phe.all$Array_ID %in% names(scr[scr>0.98]),]  #sample 137 (V;BS) excluded due

# Set models

MF=as.formula('beta.sel[i,] ~ phe$Group2 + as.factor(phe$Sex) + phe$Age + as.factor(phe$Eth3) + phe$PC1_cp + phe$PC2_cp + phe$PC3_cp + phe$PC4_cp')
models=c(MF)
modelNames=c("MF_final")

# Select tissue

tissue = "SA" # or "VA

# Run models

dir=".../K450_discovery"

for(m in 1:length(models)){
  
  # Subset betas by tissue
  
  lfla=models[[m]]; modName=modelNames[[m]];
  phe=phe.all[phe.all$Sample_Type ==tissue,]
  beta.sel=beta[,as.character(phe$Array_ID)]
  beta.sel=beta.sel[-which(rowMeans(is.na(beta.sel)) > 0.5), ]
  
  # Generate models and tabulate
  res=matrix(nrow=nrow(beta.sel), ncol=5)
  for(i in 1:nrow(beta.sel)){
    tryCatch({fit = summary(lm(lfla)) }, error = function(error) {return(NA)})
    if(!exists("fit")){
      res[i,]=NA
    }else{
      res[i,]=c(fit$coefficients[2,],fit$df[2])
    }
    if(i%%10000==0){ print(paste(modName,i,sep=" ")) }
  }
  rownames(res)=rownames(beta.sel); colnames(res)=c("beta","se","t","p","df")
  fname2 = paste("Obese_lean_",tissue,"_disc_",modName,".csv", sep="");
  write.csv(res,file=file.path(dir, fname2))
  
}

#############################################################################################

### 2. Replication Analysis (EPIC array)

# load normalised betas, control probe principal components, and annotation 

load("...beta_QN_rep.RData")
phe=read.table("...phe.rep.txt", header=T, sep='\t')
load("...ctrlpca_rep.RData")
phe.all=merge(phe.all, ctrlprobes.scores, by.y='row.names', by.x='ID')

# Remove failed samples with low call-rates (<98%)

scr=1-colSums(is.na(beta))/nrow(beta)
phe.all = phe.all[phe.all$Array_ID %in% names(scr[scr>0.98]),] 

# Set models

MF=as.formula('beta.sel[i,] ~ phe$Group2 + as.factor(phe$Sex) + phe$Age + as.factor(phe$Eth3) + phe$PC1_cp + phe$PC2_cp + phe$PC3_cp + phe$PC4_cp')
models=c(MF)
modelNames=c("MF_final")

# Select tissue

tissue = "SA" # or "VA

# Run models

dir=".../EPIC_Replication"

for(m in 1:length(models)){
  
  # Subset betas by tissue
  
  lfla=models[[m]]; modName=modelNames[[m]];
  phe=phe.all[phe.all$Sample_Type ==tissue,]
  beta.sel=beta[,as.character(phe$Array_ID)]
  beta.sel=beta.sel[-which(rowMeans(is.na(beta.sel)) > 0.5), ]
  
  # Generate models and tabulate
  res=matrix(nrow=nrow(beta.sel), ncol=5)
  for(i in 1:nrow(beta.sel)){
    tryCatch({fit = summary(lm(lfla)) }, error = function(error) {return(NA)})
    if(!exists("fit")){
      res[i,]=NA
    }else{
      res[i,]=c(fit$coefficients[2,],fit$df[2])
    }
    if(i%%10000==0){ print(paste(modName,i,sep=" ")) }
  }
  rownames(res)=rownames(beta.sel); colnames(res)=c("beta","se","t","p","df")
  fname2 = paste("Obese_lean_",tissue,"_rep_",modName,".csv", sep="");
  write.csv(res,file=file.path(dir, fname2))
  
}

#############################################################################################

### 3. Combined by meta analysis ( inverse variance weighted ).

require(qvalue)

# Load discovery and replication datasets

options(stringsAsFactors=F)

disc=read.csv(file=".../K450_discovery/Obese_lean_SA_disc_MF_final.csv")
rep=read.csv(file=".../EPIC_Replication/Obese_lean_SA_rep_MF_final.csv")

# Merge discovery and replicaiton shared 5mC sites 

all=merge(disc,rep,by="CG")

# Meta-analyse results

dir=".../Combined"

all$comb_beta = ((all$disc_beta*(1/all$disc_se^2))+(all$rep_beta*(1/all$rep_se^2)))/((1/all$disc_se^2)+(1/all$rep_se^2))
all$comb_se = sqrt( 1/ ((1/all$disc_se^2)+(1/all$rep_se^2)))
all$comb_t = all$comb_beta/all$comb_se
all$comb_p = 2*pt(-abs(all$comb_t),df=(all$disc_df + all$rep_df))

# save

write.csv(all,file=file.path(dir, ".../Obese_lean_SA_comb_MF_final.csv",row.names=F))

#############################################################################################

### 4. Identify significant CG sites, and sentinel site for each independent genomic locus (+/-5000-bp)

require(qvalue)
options(stringsAsFactors=F)

# load meta analysis results

all=read.csv(file.path(dir, ".../Obese_lean_SA_comb_MF_final.csv"))
all$direction = ifelse(sign(all$disc_beta)==sign(all$rep_beta),T,F)

# load problem probes lists

xr = read.csv(".../k450_Annotation/k450_problem_probes.csv")
snps = read.csv(".../k450_1000_genomes_probe_snps.csv")
cgs1 = xr$CG[xr$xReactiveReported==T]
cgs2 = snps$CG[!is.na(snps$probe_snps_1000g)]
problemCGs = unique(c(cgs1,cgs2))

# remove problem probe CGs

all = all[!all$CG %in% problemCGs,]

### tabulate discovery and replication thresholds as Y / N

# calculate discovery q values

disc_fdr = qvalue(all$disc_p)

# identify cg sites passing discovery and replicaiton criteria thresholds

signif = all$CG[disc_fdr$qvalues<0.01 & all$comb_p<1e-7 & all$direction == T ]
all$signif = ifelse(all$CG %in% signif,T,F)

# idenitfy sentinel cg sites (top ranking cg site wihtin +/-5000-bp)

anno=read.csv(".../k450_manifest.csv")
anno=anno[,c("CG","chr","pos","SYMBOL")]; 
allAnno=merge(all,anno,by="CG",all.x=T)

signif = allAnno[which(allAnno$signif==T),]
signif = signif[order(signif$comb_p),]

for(cg in 1:length(sent$CG)){
  cpg=sent[cg,]
  cpgChr = cpg$chr
  cpgPos = cpg$pos
  if(cg==1) {sentinel="Y"}
  if(cg>1)  {
    preceedingCpgs = sent[1:cg-1,] 
    keep=ifelse(length(which(preceedingCpgs$chr == cpgChr & abs(preceedingCpgs$pos-cpgPos)<5000))==0,T,NA)
    sentinel=c(sentinel,keep)
  }
}

signif = cbind(signif,sentinel)
all$sent = ifelse(all$CG %in% signif$CG[signif$sentinel==T],T,F)

# save
write.csv(all,file=file.path(dir, ".../Obese_lean_SA_comb_MF_final.csv",row.names=F))

#############################################################################################

