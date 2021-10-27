#############################################################################################

##### Identification of DNA methylation sites associated with obesity
##### Separately in i. subcutaneous adipocytes and ii. visceral adipocytes

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

write.csv(all,file=file=file.path(dir, ".../Obese_lean_SA_comb_MF_final.csv",row.names=F))

#############################################################################################

