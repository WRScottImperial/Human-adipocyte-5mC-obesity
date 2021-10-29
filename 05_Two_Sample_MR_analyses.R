#############################################################################################

##### Two sample Mendelian Randomisation analyses

#############################################################################################

##### Overview:

### 1. Prune Twins UK WAT mqtl SNPs (within +/- 500-kb) by LD R2 < 0.01 
### 2. Run Discovery and Sensitivity Two Sample Mendelian Randomisation analyses (in TwoSampleMR)
### 3. Identify significant exposure-outcome Two Sample MR results using most powerful tests 
### 4. Tabulate TwoSampleMR Discovery and Sensitivity results
### 5. Run sensitivity analyses using correlated mqtl cis-SNPs (LD R2 <0.8)
### 6. Combine all Two Sample MR discovery and sensitivity results

#############################################################################################

require(datatable)
require(TwoSampleMR)
require(MRInstruments)
library(qvalue)
options(stringsAsFactors=F)

#############################################################################################

### 1. Prune Twins UK WAT mqtl SNPs (within +/- 500-kb of each sentinel) by LD R2 < 0.01 

# load mqtls 

mqtls = fread(".../Twins_WAT_mQTLs_sentinels_snps_500kbp.txt")

# load intersected snps and subset mqtl SNPs (intersected SNPs are those present in each of the human gwas summary results ~N=2.5M) 

load(".../Intersected_SNPs_for_MR_analyses_all_GWAS.RData")
gwas_snps = as.character(unique(allGwasSnps))
MQTLS = mqtls[mqtls$rsid %in% gwas_snps,]

# Split MQTLS by discovery tissue  

mosaSent = read.csv(".../Obese_lean_SA_Combined_sentinels.csv")
movaSent = read.csv(".../Obese_lean_VA_Combined_sentinels.csv")
scMqtls = MQTLS[MQTLS$PHENOTYPE %in% as.character(mosaSent$CG),]
viMqtls = MQTLS[MQTLS$PHENOTYPE %in% as.character(movaSent$CG),]
mqtlsList = list(scMqtls,viMqtls)
names(mqtlsList) = c("sc","vi")

# clump mqtls by LD of 0.01 (uncorrelated discovery/sensitivity) and 0.8 (correlated sensitivity) using ieugwasr::ld_clump keeping strongest MQTL snp

lds = c(0.01,0.8)

for(l in 1:length(lds)){
  
  # select ld
  ld = lds[l]
  
  for(m in 1:length(mqtlsList)){
    
    # select tissue
    
    mQTLs = mqtlsList[[m]]
    
    # split mqtls by sentinel CG site
    
    mQTLs = split(mQTLs, mQTLs$PHENOTYPE)
    
    for(c in 1:length(mQTLs)){
      
      # get SNPs table for each sentinel 5mC site
      
      print(paste("LD=",ld," ",c,sep="")); 
      cisSnp=mQTLs[[c]];
      
      # return NA id no cis- SNPs 
      
      if(dim(cisSnp)[1]==0) { cisSnpClump = rep(NA,15) } else {
        
        # tabulate as TwoSample MR format
        
        cisSnp = format_data(cisSnp,type='exposure',phenotype_col="PHENOTYPE", snp_col="rsid", 
                             beta_col="BETA", se_col="SE", effect_allele_col="EA",
                             other_allele_col="NEA", pval_col="PVAL", samplesize_col="N",
                             chr_col="chr", pos_col="pos", eaf_col = "EAF")
        names(cisSnp)[8] = 'rsid'; names(cisSnp)[3] = 'pval'
        
        # clump SNPs based on selected LD
        
        clump = ieugwasr::ld_clump(cisSnp,clump_r2=ld)
        
      }
      
      # list and tabulate clumped mQTLs for each sentinel
      
      if(c==1){all_list=list(clump); all_table = clump}
      if(c>1){all_list=c(all_list,list(clump)); all_table = rbind(all_table,clump) }
      
    }
    
    # generate adjusted pvalues
    
    all_table$padj_bf = p.adjust(all_table$pval,method="bonferroni",n=length(all_table$pval))
    all_table$padj_fdr = qvalue::qvalue(all_table$pval)$qvalue
    
    if(m==1) {watMqtls = list(all_list,all_table)}
    if(m>1) {watMqtls = c(watMqtls,list(all_list,all_table)) }
    
  }

  names(watMqtls)=c("scMqtlsList","scMqtlsTable","viMqtlsList","viMqtlsTable")

  dir = ".../mQTLs/WAT/"

  save(watMqtls,file=file.path(dir,paste("/WAT_mQTLs_obese_lean_sent_",ld,".RData",sep="")))

}

#############################################################################################

### 2. Run Discovery and sensitivity two sample Mendelian Randomisation analyses on mqtls clumped at LD R2 < 0.01 (in TwoSampleMR)

# load clumped mqtls LD R2 0.01

load(".../mQTLs/WAT/WAT_mQTLs_obese_lean_sent_0.01.RData")

# load 5mC ~ Obesity and metabolic disease GWAS summary data (this is a list of BMI, WHR, T2D, T2DadjWHR, FG, FI, HbA1c summary tables) 

load("/rds/general/project/lms-scott-raw/live/Other_datasets/GWAS_summary_datasets/GWAS_results_tables_list_for_MR.RData")
gwas=allGwasRes

# select mqtls based on stringent bonferroni correction to reduce weak instrument bias 

mqtlsList = list(watMqtls$scMqtlsTable,watMqtls$viMqtlsTable)
names(mqtlsList) =c("sc_mqtls","v_mqtls")

for(m in 1:length(mqtlsList)){
  mqtl=mqtlsList[[m]]
  mqtl = mqtl[which(mqtl$padj_bf<0.05),]
  mqtlsList[[m]]=mqtl
}

# format exposure (mqtl) data for TwoSampleMMR

for(m in 1:length(mqtlsList)){
  mqtl=mqtlsList[[m]]
  MQTL <- format_data(
    mqtl, type="exposure", snp_col = "rsid", phenotype_col = "exposure", 
    beta_col = "beta.exposure",  se_col = "se.exposure", pval_col = "pval",
    effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure", samplesize_col = "samplesize.exposure", chr_col = "chr.exposure", pos_col = "pos.exposure"
  )
  MQTL$pair = paste(MQTL$SNP,MQTL$exposure,sep="_")
  mqtlsList[[m]] = MQTL
}

# format outcome (gwas) data for TwoSampleMR

for(g in 1:length(gwas)) {
  print(g)
  outcomeGwas <- format_data(
    gwas[[g]], snps=MQTL$SNP, type="outcome", snp_col = "rsid", beta_col = "beta",  se_col = "se", pval_col = "p",
    effect_allele_col = "EA", other_allele_col = "NonEA", eaf_col = "EAF", samplesize_col = "n",
    chr_col = "chr", pos_col = "pos", phenotype_col = "gwas"
  )
  if(g==1) { outcomeGwasList = list(outcomeGwas) };
  if(g>1)  { outcomeGwasList = c(outcomeGwasList, list(outcomeGwas)) }
}
names(outcomeGwasList)=names(gwas)

# select tissue

MQTL = mqtlsList$sc_mqtls # or MQTL = mqtlsList$vi_mqtls
tissue = "SA" # or tissue = "VA"

# harmonise snp directions and remove ambiguous palendromic SNPs using TwoSampleMR

for(ogl in 1:length(outcomeGwasList)){
  GWAS = outcomeGwasList[[ogl]]
  harmonisedGwas = harmonise_data(exposure_dat = MQTL, outcome_dat = GWAS)
  harmonisedGwas$snp = paste(harmonisedGwas$SNP,harmonisedGwas$effect_allele.exposure,harmonisedGwas$other_allele.exposure,sep="_")
  if(ogl==1) { harmonisedGwasList = list(harmonisedGwas) }
  if(ogl>1)  { harmonisedGwasList = c(harmonisedGwasList,list(harmonisedGwas)) }
}
names(harmonisedGwasList)=names(gwas)

# run TwoSampleMR 

for(hgl in 1:length(harmonisedGwasList)){
  
  dat = harmonisedGwasList[[hgl]]
  
  mrBasic = mr(dat)
  mrSingle = mr_singlesnp(dat)
  mrSingle$snp =dat$snp[match(mrSingle$SNP,dat$SNP)]
  mrSingle$EA =dat$effect_allele.exposure[match(mrSingle$SNP,dat$SNP)]
  mrSingle$NEA =dat$other_allele.exposure[match(mrSingle$SNP,dat$SNP)]
  mrHet = mr_heterogeneity(dat)
  mrHZP = mr_pleiotropy_test(dat) 
  steiger = directionality_test(dat)
  
  phe = list(mrBasic,mrSingle,mrHet,mrHZP,steiger)
  names(phe) = c("MR_results","MR_results_single_snps","MR_heterogenity","MR_hz_pleiotrophy","MR_steiger")
  
  if(hgl==1) { mrResList = list(phe) }
  if(hgl>1)  { mrResList = c(mrResList,list(phe)) }
  
}

names(mrResList)=names(gwas)

save(mrResList,file=".../mQTLs/WAT/TwoSampleMR/TwoSampleMR_results_SA_mqtls_gwas.RData")

#############################################################################################

### 3. Identify significant exposure-outcome Two Sample MR results using most powerful tests (Discovery MR)
# (i.e. single instrumental variable (SNP) with Wald ratio, and multiple instrumental variables (SNPs) with IVW)

mr_file = ".../mQTLs/WAT/TwoSampleMR/TwoSampleMR_results_SA_mqtls_gwas.RData"
load(mr_file)
  
for(m in 1:length(mrResList)){
  
  # subset IVW and Wald results
  res = mrResList[[m]]$MR_results
  res = res[res$method=="Wald ratio"|res$method=="Inverse variance weighted",]
  
  # adjusted p values
  res$padj_bf = p.adjust(res$pval,method="bonferroni",n=length(res$pval))
  res$padj_fdr = qvalue::qvalue(res$pval)$qvalue
  
  # evaluate concordance with methylation directions of effect (CG-obesity assoacition)
  if( grepl("SA",mr_file) ) { print("use mosa betas"); betas = match(as.character(res$exposure),mosaSent$CG) ; res = cbind(res,mosaSent[betas,c("CG","comb_beta","comb_se","comb_t","comb_p")]) }
  if( grepl("VA",mr_file) ) { print("use mova betas"); betas = match(as.character(res$exposure),movaSent$CG) ; res = cbind(res,movaSent[betas,c("CG","comb_beta","comb_se","comb_t","comb_p")]) }
  res$concord = ifelse(sign(res$b)==sign(res$comb_beta),TRUE,FALSE)
  
  # add steiger directionality results
  steiger = mrResList[[m]]$MR_steiger 
  res = cbind(res,steiger[match(res$id.exposure,steiger$id.exposure),c("snp_r2.exposure","snp_r2.outcome","correct_causal_direction","steiger_pval")])
  
  # Exponentiate T2D results for adds ration 
  if ( grepl("t2d",names(mrResList)[m]) )  {
    res$or <- exp(res$b)
    res$or_lci95 <- exp(res$b - res$se * 1.96)
    res$or_uci95 <- exp(res$b + res$se * 1.96)
  } 
  
  # order results
  res=res[order(res$pval),]
  
  mrResList[[m]]$MR_discovery = res
  
}

save(mrResList,file=mr_file)

# tabulate for all gwas phenotypes

for(mr in 1:length(mrResList)){
  
  res = mrResList[[mr]]$MR_discovery
  res = res[res$padj_fdr<0.01,]
  
  if(mr==1) { all = res }
  if(mr>1) { all = rbind(all,res) }
  
  write.csv(all,file=".../mQTLs/WAT/TwoSampleMR/TwoSampleMR_results_SA_mqtls_gwas_table.csv",row.names=F)
  
}

#############################################################################################

### 4. Combine TwoSampleMR Sensitivity results

load(mr_file)
  
for(m in 1:length(mrResList)){
  
  print(m)
  
  requireNamespace("plyr", quietly = TRUE)
  
  res=mrResList[[m]]$MR_results
  het=mrResList[[m]]$MR_heterogenity
  plt=mrResList[[m]]$MR_hz_pleiotrophy
  sin=mrResList[[m]]$MR_results_single_snps
  
  sin <- sin[grep("[:0-9:]", sin$SNP), ]
  sin$method <- "Wald ratio"
  names(sin)[names(sin) == "p"] <- "pval"
  names(res)[names(res) == "method"] <- "Method"
  names(sin)[names(sin) == "method"] <- "Method"
  
  if( dim(het)[1]>0 & dim(plt)[1]>0 ) {
    
    het <- het[, c("id.exposure", "id.outcome", "method", "Q", "Q_df", "Q_pval")]
    names(het)[names(het) == "method"] <- "Method"
    plt <- plt[, c("id.outcome", "id.exposure", "egger_intercept", "se", "pval")]
    plt$Method <- "MR Egger"
    names(plt)[names(plt) == "egger_intercept"] <- "intercept"
    names(plt)[names(plt) == "se"] <- "intercept_se"
    names(plt)[names(plt) == "pval"] <- "intercept_pval"
    res <- merge(res, het, by = c("id.outcome", "id.exposure", "Method"), all.x = T)
    res <- merge(res, plt, by = c("id.outcome", "id.exposure", "Method"), all.x = T) 
  } else {res = res}
  
  for(i in 1:nrow(res)){
    if(res$Method[i]=="Wald ratio") {res$rsid[i]=sin$SNP[sin$id.exposure==res$id.exposure[i]]; res$SNP[i]=sin$snp[sin$id.exposure==res$id.exposure[i]] } else {res$rsid[i] = NA; res$SNP[i] = NA}
  }
  for (i in unique(res$id.exposure)) {
    Methods <- unique(res$Method[res$id.exposure == i]); Methods <- Methods[Methods != "Wald ratio"]
    for (j in unique(Methods)) {
      res$rsid[res$id.exposure == i & res$Method == j] <- paste(sin$SNP[sin$id.exposure == i & sin$Method == "Wald ratio"], collapse = ";")
      res$SNP[res$id.exposure == i & res$Method == j] <- paste(sin$snp[sin$id.exposure == i & sin$Method == "Wald ratio"], collapse = ";")
    }
  }
  
  res=res[order(res$exposure),]
  mrResList[[m]]$MR_results_with_snps = res
  
}

save(mrResList,file=mr_file)

#############################################################################################

### 5. Run sensitivity analyses using correlated mqtl cis-SNPs (using LD R2 <0.8)

# load clumped mqtls LD R2 < 0.8

load(".../mQTLs/WAT/WAT_mQTLs_obese_lean_sent_0.8.RData")

# select correlated mqtls based on stringent bonferroni correction to reduce weak instrument bias 

mqtlsList = list(watMqtls$scMqtlsTable,watMqtls$viMqtlsTable)
names(mqtlsList) =c("sc_mqtls","v_mqtls")

for(m in 1:length(mqtlsList)){
  mqtl=mqtlsList[[m]]
  mqtl = mqtl[which(mqtl$padj_bf<0.05),]
  mqtlsList[[m]]=mqtl
}

# format exposure (mqtl) data 

for(m in 1:length(mqtlsList)){
  mqtl=mqtlsList[[m]]
  MQTL <- format_data(
    mqtl, type="exposure", snp_col = "rsid", phenotype_col = "exposure", 
    beta_col = "beta.exposure",  se_col = "se.exposure", pval_col = "pval",
    effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure", samplesize_col = "samplesize.exposure", chr_col = "chr.exposure", pos_col = "pos.exposure"
  )
  MQTL$pair = paste(MQTL$SNP,MQTL$exposure,sep="_")
  mqtlsList[[m]] = MQTL
}

# format outcome (gwas) data 

for(g in 1:length(gwas)) {
  print(g)
  outcomeGwas <- format_data(
    gwas[[g]], snps=MQTL$SNP, type="outcome", snp_col = "rsid", beta_col = "beta",  se_col = "se", pval_col = "p",
    effect_allele_col = "EA", other_allele_col = "NonEA", eaf_col = "EAF", samplesize_col = "n",
    chr_col = "chr", pos_col = "pos", phenotype_col = "gwas"
  )
  if(g==1) { outcomeGwasList = list(outcomeGwas) };
  if(g>1)  { outcomeGwasList = c(outcomeGwasList, list(outcomeGwas)) }
}
names(outcomeGwasList)=names(gwas)

# select tissue

MQTL = mqtlsList$sc_mqtls # or MQTL = mqtlsList$vi_mqtls
tissue = "SA" # or tissue = "VA"

# harmonise snp directions and remove ambiguous palindromic SNPs 

for(ogl in 1:length(outcomeGwasList)){
  GWAS = outcomeGwasList[[ogl]]
  harmonisedGwas = harmonise_data(exposure_dat = MQTL, outcome_dat = GWAS)
  harmonisedGwas$snp = paste(harmonisedGwas$SNP,harmonisedGwas$effect_allele.exposure,harmonisedGwas$other_allele.exposure,sep="_")
  if(ogl==1) { harmonisedGwasList = list(harmonisedGwas) }
  if(ogl>1)  { harmonisedGwasList = c(harmonisedGwasList,list(harmonisedGwas)) }
}
names(harmonisedGwasList)=names(gwas)

# Run correlated MR sensitivity analyses

for(g in 1:length(harmonisedGwasList)){
  
  dat=harmonisedGwasList[[g]]
  dat = split(dat,dat$exposure)
  
  for(t in 1:length(dat)){
    
    print(paste(names(harmonisedGwasList)[g],t,nrow(dat[[t]]),sep=" "))
    exposure = dat[[t]]$exposure[1]
    
    # if more than one mqtl use correlation matrix for correlated MR
    if(dim(dat[[t]])[1]>1) {
      cor = ieugwasr::ld_matrix(unique(dat[[t]]$SNP))
      MRInput <- MendelianRandomization::mr_input(bx = dat[[t]]$beta.exposure,bxse = dat[[t]]$se.exposure,
                                                  by = dat[[t]]$beta.outcome,byse = dat[[t]]$se.outcome,
                                                  snps = dat[[t]]$SNP,
                                                  corr = cor,
                                                  exposure=dat[[t]]$exposure, outcome=dat[[t]]$outcome,
                                                  effect_allele=dat[[t]]$effect_allele.exposure, 
                                                  other_allele=dat[[t]]$other_allele.outcom, 
                                                  eaf=dat[[t]]$eaf.outcome)
    } else {
      
      MRInput <- MendelianRandomization::mr_input(bx = dat[[t]]$beta.exposure,bxse = dat[[t]]$se.exposure,
                                                  by = dat[[t]]$beta.outcome,byse = dat[[t]]$se.outcome,
                                                  snps = dat[[t]]$SNP,
                                                  exposure=dat[[t]]$exposure, outcome=dat[[t]]$outcome,
                                                  effect_allele=dat[[t]]$effect_allele.exposure, 
                                                  other_allele=dat[[t]]$other_allele.outcom, 
                                                  eaf=dat[[t]]$eaf.outcome)
    }
    
    ivw <- MendelianRandomization::mr_ivw(MRInput)
    
    if(length(MRInput@snps)>2) {  egger <- MendelianRandomization::mr_egger(MRInput)
    } else ( egger = rep(NA,15) )
    
    SNPs = paste(rownames(ivw@Correlation),collapse=";")
    
    IVW = c(ivw@Exposure[1],ivw@Outcome[1],ivw@Estimate, ivw@StdError, ivw@CILower, ivw@CIUpper, ivw@Pvalue, ivw@Heter.Stat, ivw@SNPs  )
    
    if(class(egger)=="Egger") { EGGER = c(egger@Estimate, egger@StdError.Est, egger@CILower.Est, egger@CIUpper.Est, egger@Pvalue.Est,
                                          egger@Intercept, egger@StdError.Int, egger@CILower.Int, egger@CIUpper.Int, egger@Pvalue.Int, 
                                          egger@Pleio.pval, egger@Causal.pval, egger@Heter.Stat, egger@I.sq) 
    } else (EGGER = egger)
    
    all = c(IVW,EGGER,SNPs)          
    res = list(ivw,egger)
    names(res) = c("ivw","egger")
    
    if(t==1) {mr = list(res); exp = exposure; allmr = all}
    if(t>1){ mr = c(mr,list(res)); exp = c(exp,exposure); allmr = rbind(allmr,all) }
    
  }
  names(mr)=exp
  allmr=as.data.frame(allmr)
  colnames(allmr)=c("exposure","outcome","ivw_beta","ivw_se","ivw_CIlow","ivw_CIup","ivw_pval","ivw_het_stat","ivw_het_pval","nsnps",
                    "egger_beta","egger_se","egger_CIlow","egger_CIup","egger_pval",
                    "egger_intercept","egger_se_inter","egger_CIlow_inter","egger_CIup_inter","egger_pval_inter",
                    "egger_pval_pleio","egger_pval_causal","egger_het_stat","egger_het_pval","egger_het_I2","snps")
  allmr=allmr[order(allmr$ivw_pval),]
  
  # save results
  
  save(mr,allmr,file=paste(".../mQTLs/WAT/TwoSampleMR/Correlated/corr_MR_results_SA_mqtls_",names(harmonisedGwasList)[g],"_gwas.RData",sep=""))
  
}

# Tabulate correlated MR results 

mr_files=list.files(".../mQTLs/WAT/TwoSampleMR/Correlated/",full.names=T)
files = unlist(sapply(strsplit(basename(mr_files),split=".RData"),"[",1))

for(f in 1:length(mr_files)){
  
  load(mr_files[f])
  
  char = c("exposure","outcome","snps")
  numeric = allmr[,!colnames(allmr) %in% char]
  numeric = as.data.frame(apply(numeric, 2, as.numeric))
  
  allmr = cbind(allmr[,c("exposure","outcome")],numeric,allmr[,"snps"])
  allmr$ivw_padj_bf = p.adjust(as.numeric(allmr$ivw_pval),method="bonferroni",n=length(allmr$ivw_pval))
  allmr$ivw_padj_fdr = qvalue::qvalue(as.numeric(allmr$ivw_pval))$qvalue
  
  # evaluate concordance with methylation directions of effect (CG-obesity assoacition)
  if( grepl("SA",files[f]) ) { print("use mosa betas"); betas = mosaSent }
  if( grepl("VA",files[f]) ) { print("use mova betas"); betas = movaSent }
  
  allmr = cbind(allmr,betas[,c("CG","comb_beta","comb_se","comb_t","comb_p")])
  allmr$concord = ifelse(sign(allmr$ivw_beta)==sign(allmr$comb_beta),TRUE,FALSE)
  
  write.csv(allmr,file=".../mQTLs/WAT/TwoSampleMR/Correlated/corr_MR_results_SA_mqtls_all_gwas.csv",row.names=F)
  
}

#############################################################################################

### 6. Combine all Two Sample MR discovery and sensitivity results

# Load TwoSampleMR results (LD<0.01)

load(".../mQTLs/WAT/TwoSampleMR/TwoSampleMR_results_SA_mqtls_gwas.RData")
mrBase = mrResList

# Load correlated MR results (LD<0.8)

mrCor = read.csv(".../mQTLs/WAT/TwoSampleMR/Correlated/corr_MR_results_SA_mqtls_all_gwas.csv",full.names=T)

# Intersect results

for(mr in 1:length(mrBase)){

  base = mrBase[[mr]]$MR_discovery # discovery results
  snp = mrBase[[mr]]$MR_results_with_snps # sensitivity results
  
  base$MR_base_steiger_qval = qvalue::qvalue(base$steiger_pval)$qvalue
  
  cor = mrCor
  
  all = merge(base,snp[,c("id.exposure","rsid","SNP")],by="id.exposure")
  all = merge(all,cor,by=c("exposure","outcome","CG","comb_beta","comb_se","comb_t","comb_p"),all.x=T)
  all = all[order(all$pval),]
  
  colnames(all) = c("exposure","outcome","CG","cg_mo_beta","cg_mo_se","cg_mo_t",
                    "cg_mo_p","id.exposure","id.outcome","MR_base_method","MR_base_nsnp",
                    "MR_base_b","MR_base_se","MR_base_pval","MR_base_padj_bf","MR_base_padj_fdr",
                    "MR_base_concord","MR_base_snp_r2.exposure","MR_base_snp_r2.outcome",
                    "MR_base_correct_causal_direction","MR_base_steiger_pval","MR_base_steiger_qval",
                    "MR_base_rsid","MR_base_snps",
                    "MR_cor_ivw_beta","MR_cor_ivw_se","MR_cor_ivw_CIlow","MR_cor_ivw_CIup",
                    "MR_cor_ivw_pval","MR_cor_ivw_het_stat","MR_cor_ivw_het_pval",
                    "MR_cor_nsnps","MR_cor_egger_beta","MR_cor_egger_se","MR_cor_egger_CIlow",
                    "MR_cor_egger_CIup","MR_cor_egger_pval","MR_cor_egger_intercept",
                    "MR_cor_egger_se_inter","MR_cor_egger_CIlow_inter","MR_cor_egger_CIup_inter",
                    "MR_cor_egger_pval_inter","MR_cor_egger_pval_pleio","MR_cor_egger_pval_causal",
                    "MR_cor_egger_het_stat","MR_cor_egger_het_pval","MR_cor_egger_het_I2",
                    "MR_cor_snps","MR_cor_ivw_padj_bf","MR_cor_ivw_padj_fdr","MR_cor_concord")
  
  # select significant results
  
  dat = dat[dat$MR_base_padj_fdr<0.01 & dat$MR_base_steiger_qval<0.01 & dat$MR_base_concord==T,]
  
  # tabulate
  
  if(f==1){alldat = dat}
  if(f>1) {alldat = rbind(alldat,dat)}
  
  # save final results
  
  write.csv(all,file=".../mQTLs/WAT/TwoSampleMR/TwoSampleMR_all_discovery_sensitivity_results_SA_.csv",row.names = F)
  
}

#############################################################################################