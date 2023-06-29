

suppressMessages(library("plyr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("data.table", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("crayon", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggplot2", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("farver", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("labeling", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("optparse", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("dplyr", lib.loc = "/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("withr", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("backports", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("broom", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("rstudioapi", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("cli", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tzdb", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("tidyverse", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))
suppressMessages(library("svglite", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggeasy", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("sandwich", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("digest", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))
suppressMessages(library("ggforce", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))

suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)


collate_results_function = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
   
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Selected_vars ----
  
  
  SELECTED_VARS = unlist(strsplit(opt$SELECTED_VARS, split=","))
  
  cat("SELECTED_VARS_\n")
  cat(sprintf(as.character(SELECTED_VARS)))
  cat("\n")
  
  
  
 
  
  ##### LOOP TO READ ALL VARIABLES -----
  list_genes<-list()
 
 
  
  Condition_DEBUG <- 1
  
  for(i in 1:length(SELECTED_VARS))
  {
    
    SELECTED_VARS_sel<-SELECTED_VARS[i]
    
    cat("--------------->\t")
    cat(sprintf(as.character(SELECTED_VARS_sel)))
    cat("\n")
    
    
    ##### path6 ---
    
    path6<-paste(out,SELECTED_VARS_sel,'/', sep='')
    
    # cat("path6\n")
    # cat(sprintf(as.character(path6)))
    # cat("\n")
    
    
    setwd(path6)
    
    
    
    filename=paste(c('DTU_LogLM_HET_RESULTS_NOMINAL_',SELECTED_VARS_sel,'.rds'),collapse='')
    
   
    
    if (file.exists(filename)) {
      DTU_genes = readRDS(filename)
      DTU_genes$VAR<-SELECTED_VARS_sel
      
      
       # quit(status = 1)
      
      if(Condition_DEBUG == 1)
      {
        cat("DTU_genes_0\n")
        cat(str(DTU_genes))
        cat("\n")
      }
      
    
      
      list_genes[[i]]<-DTU_genes
      
      # ##########################################
      # quit(status = 1)
      
    }#file.exists(DTU_genes.rds)
      
    
    
    
  }# i in 1:length(SELECTED_VARS)
  
  # ##################################
  # quit(status = 1)
  
  
  Condition_DEBUG <- 1
  
 
  
  # list_DEF[[3]]<-Results_CIS_subset 
  
  
 
  
  Results_DEF<-data.frame()
  
 
 
  
  if(length(list_genes) >0)
  {

    Results_genes_nominal = as.data.frame(data.table::rbindlist(list_genes, fill=T), stringsAsFactors=F)



    if(Condition_DEBUG == 1)
    {
      cat("Results_genes_nominal_0\n")
      cat(str(Results_genes_nominal))
      cat("\n")
      cat(str(unique(Results_genes_nominal$VAR)))
      cat("\n")
    }
    
    
    setwd(out)
    
    write.table(Results_genes_nominal, file="DTU_results_ALLTOGETHER.tsv", sep="\t", quote = F, row.names = F)
    saveRDS(Results_genes_nominal, file="DTU_results_ALLTOGETHER.rds")


   
  }# length(list_genes) >0
 
  
  
 }

function_CIS= function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  ####Read Results_genes_nominal----
  
  setwd(out)
  
  Results_genes_nominal<-readRDS(file="DTU_results_ALLTOGETHER.rds")
  
  cat("Results_genes_nominal_0\n")
  cat(str(Results_genes_nominal))
  cat("\n")
  cat(str(unique(Results_genes_nominal$VAR)))
  cat("\n")
  
  #### VEP_CSQ ----
  
  VEP_CSQ<-as.data.frame(fread(file=opt$VEP_CSQ,sep=",") , stringsAsFactors=F)
  
  cat("VEP_CSQ_0\n")
  cat(str(VEP_CSQ))
  cat("\n")
  cat(str(unique(VEP_CSQ$VAR)))
  cat("\n")
  
  #### CIS consequences -----
  
  CIS_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3",
                 "INTRON","SPLICE","UPSTREAM")
  
  #### SELECT CIS AND VARIANTS ----
  
  VEP_CSQ_sel<-VEP_CSQ[which(VEP_CSQ$VAR%in%Results_genes_nominal$VAR),]
  
  cat("VEP_CSQ_sel_0\n")
  cat(str(VEP_CSQ_sel))
  cat("\n")
  cat(str(unique(VEP_CSQ_sel$VAR)))
  cat("\n")
  
  VEP_CSQ_sel_CIS<-VEP_CSQ_sel[which(VEP_CSQ_sel$VEP_DEF_LABELS%in%CIS_LABELS),]
  
  cat("VEP_CSQ_sel_CIS_0\n")
  cat(str(VEP_CSQ_sel_CIS))
  cat("\n")
  cat(str(unique(VEP_CSQ_sel_CIS$VAR)))
  cat("\n")
  
  
  Results_genes_nominal_CIS<-merge(Results_genes_nominal,
                                   VEP_CSQ_sel_CIS[,-which(colnames(VEP_CSQ_sel_CIS) == 'HGNC')],
                                   by=c("VAR","ensembl_gene_id"))
  
  cat("Results_genes_nominal_CIS_0\n")
  cat(str(Results_genes_nominal_CIS))
  cat("\n")
  cat(str(unique(Results_genes_nominal_CIS$VAR)))
  cat("\n")
  
  
  if(dim(Results_genes_nominal_CIS)[1] >0)
  {
    
    ##### Correction per cell type -----
    
    
    CT_levels<-unique(as.character(Results_genes_nominal_CIS$Cell_Type))
    
    cat("CT_levels\n")
    cat(str(CT_levels))
    cat("\n")
    
    list_CT<-list()
    
    for(iteration_CT in 1:length(CT_levels))
    {
      
      CT_levels_sel<-CT_levels[iteration_CT]
      
      cat("-------------------------------------------------------------------------------------------------------------------------------CT------>\t")
      cat(sprintf(as.character(CT_levels_sel)))
      cat("\n")
      
      Results_genes_nominal_CIS_CT_sel<-Results_genes_nominal_CIS[which(Results_genes_nominal_CIS$Cell_Type == CT_levels_sel),]
      
      cat("Results_genes_nominal_CIS_CT_sel\n")
      cat(str(Results_genes_nominal_CIS_CT_sel))
      cat("\n")
      
      
      VARS_array<-unique(Results_genes_nominal_CIS_CT_sel$VAR)
      
      list_VAR<-list()
      
      for(iteration_VARS_array in 1:length(VARS_array))
      {
        
        VAR_sel<-VARS_array[iteration_VARS_array]
        
        cat("-----------VAR------->\t")
        cat(sprintf(as.character(VAR_sel)))
        cat("\n")
        
        Results_genes_nominal_CIS_CT_sel_VAR_sel<-Results_genes_nominal_CIS_CT_sel[which(Results_genes_nominal_CIS_CT_sel$VAR%in%VAR_sel),]
        
        cat("Results_genes_nominal_CIS_CT_sel_VAR_sel\n")
        cat(str(Results_genes_nominal_CIS_CT_sel_VAR_sel))
        cat("\n")
        cat(str(unique(Results_genes_nominal_CIS_CT_sel_VAR_sel$VAR)))
        cat("\n")
        cat(str(unique(Results_genes_nominal_CIS_CT_sel_VAR_sel$HGNC)))
        cat("\n")
        
        #### SELECT TESTS TO CORRECT ----
        
        
        indx.select<-c(which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "VAR"),
                       which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "ensembl_gene_id"),which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "HGNC"),which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "transcript_id"),
                       which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "coefficient_Genotypes"),
                       which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "pvalue_Genotypes"),
                       which(colnames(Results_genes_nominal_CIS_CT_sel_VAR_sel)== "n_breakdown_string"))
        
        
        
        MT_set<-unique(Results_genes_nominal_CIS_CT_sel_VAR_sel[,indx.select])
        
        cat("MT_set_0\n")
        cat(str(MT_set))
        cat("\n")
        cat(str(unique(MT_set$ensembl_gene_id)))
        cat("\n")
        
        MT_set$ajusted.pvalue_Genotypes<-p.adjust(MT_set$pvalue_Genotypes, method = "BH")
        MT_set$ajusted.minuslogpvalue_Genotypes<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes),2)
        
        
        
        cat("MT_set_1\n")
        cat(str(MT_set))
        cat("\n")
        
        list_VAR[[iteration_VARS_array]]<-MT_set
        
      }#iteration_VARS_array in 1:length(VARS_array)
      
      if(length(list_VAR) > 0)
      {
        VAR_df = as.data.frame(data.table::rbindlist(list_VAR, fill=T), stringsAsFactors=F)
        
        VAR_df$Cell_Type<-CT_levels_sel
        
        cat("VAR_df\n")
        cat(str(VAR_df))
        cat("\n")
        
        list_CT[[iteration_CT]]<-VAR_df
        
        
      }#length(list_VAR) > 0
    }# iteration_CT in 1:length(CT_levels)
    
    if(length(list_CT) > 0)
    {
      CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
      
      cat("CT_df_0\n")
      cat(str(CT_df))
      cat("\n")
      
      
      check<-CT_df[which(CT_df$ajusted.minuslogpvalue_Genotypes >=1.3),]
      
      cat("check\n")
      cat(str(check))
      cat("\n")
      cat(sprintf(as.character(names(summary((unique(droplevels(interaction(check$VAR,check$HGNC)))))))))
      cat("\n")
      cat(sprintf(as.character(summary((unique(droplevels(interaction(check$VAR,check$HGNC))))))))
      cat("\n")
      
          
      setwd(out)
      
      if(dim(CT_df)[1] >0)
      {
        saveRDS(file="DTU_CIS_gene.rds", CT_df)  
      }
      
    }#length(list_CT) > 0
    
  }# dim(Results_genes_nominal_CIS)[1] >0
  
 
}

function_Block= function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
 
  
 
  #### Read ALL_dB file ----
  
  ALL_dB<-as.data.frame(fread(file=opt$ALL_dB,sep="\t") , stringsAsFactors=F)
  
  cat("ALL_dB_0\n")
  cat(str(ALL_dB))
  cat("\n")
  
  
  indx.dep<-c(which(colnames(ALL_dB) == "maf_origin"))
  
  ALL_dB_subset<-unique(ALL_dB[,-indx.dep])
  
  # cat("ALL_dB_subset\n")
  # cat(str(ALL_dB_subset))
  # cat("\n")
  
  
  ALL_dB_subset$Allelic_Series_ID<-paste(ALL_dB_subset$phenotype,ALL_dB_subset$block_no,sep='__')
  
  # cat("ALL_dB_subset_1\n")
  # cat(str(ALL_dB_subset))
  # cat("\n")
  
  
  indx.int<-c(which(colnames(ALL_dB_subset) == "VAR"),which(colnames(ALL_dB_subset) == "Allelic_Series_ID"))
  
  ALL_FINAL<-unique(ALL_dB_subset[,indx.int])
  
  cat("ALL_FINAL\n")
  cat(str(ALL_FINAL))
  cat("\n")
  
  cat("PCHiC_info_file:\n")
  cat(sprintf(as.character(opt$PCHiC_info)))
  cat("\n")
  
  PCHiC_info<-as.data.frame(fread(file=opt$PCHiC_info, sep =",", header =T), stringsAsFactors=F)
  
  cat("PCHiC_info:\n")
  cat(str(PCHiC_info))
  cat("\n")
  
  PCHiC_info_subset<-unique(PCHiC_info[,c(which(colnames(PCHiC_info) == "VAR"),
                                          which(colnames(PCHiC_info) == "ensembl_gene_id"))])
  
  cat("PCHiC_info_subset:\n")
  cat(str(PCHiC_info_subset))
  cat("\n")
  
  
  #### GENES_PER_BLOCKS ----
  
  GENES_PER_BLOCKS = as.data.frame(fread(file=opt$GENES_PER_BLOCKS, sep="\t", stringsAsFactors = F, header = T))
  
  cat("GENES_PER_BLOCKS\n")
  cat(str(GENES_PER_BLOCKS))
  cat("\n")
  
  GENES_PER_BLOCKS$phenotype<-gsub("^.+__","",GENES_PER_BLOCKS$BLOCK, perl=T)
  GENES_PER_BLOCKS$block_no<-gsub("__.+$","",GENES_PER_BLOCKS$BLOCK, perl=T)
  
  GENES_PER_BLOCKS$Allelic_Series_ID<-paste(GENES_PER_BLOCKS$phenotype,GENES_PER_BLOCKS$block_no,sep="__")
  
  GENES_PER_BLOCKS<-GENES_PER_BLOCKS[,-c(which(colnames(GENES_PER_BLOCKS) == "phenotype"),
                                         which(colnames(GENES_PER_BLOCKS) == "block_no"))]
  
  
  # cat("GENES_PER_BLOCKS\n")
  # cat(str(GENES_PER_BLOCKS))
  # cat("\n")
  # 
  
  AS_Gene_source<-merge(ALL_FINAL,
                        GENES_PER_BLOCKS,
                        by="Allelic_Series_ID")
  
  
  cat("AS_Gene_source\n")
  cat(str(AS_Gene_source))
  cat("\n")
  
  AS_Gene_source_subset<-unique(AS_Gene_source[,c(which(colnames(AS_Gene_source) == "VAR"),
                                                  which(colnames(AS_Gene_source) == "ensembl_gene_id"))])
  
  cat("AS_Gene_source_subset:\n")
  cat(str(AS_Gene_source_subset))
  cat("\n")
  
  
  Block_and_PCHiC_candidates<-rbind(PCHiC_info_subset,
                                    AS_Gene_source_subset)
  
  cat("Block_and_PCHiC_candidates:\n")
  cat(str(Block_and_PCHiC_candidates))
  cat("\n")
 
  
  
  ####Read Results_genes_nominal----
  
  setwd(out)
  
  Results_genes_nominal<-readRDS(file="DTU_results_ALLTOGETHER.rds")
  
  cat("Results_genes_nominal_0\n")
  cat(str(Results_genes_nominal))
  cat("\n")
  cat(str(unique(Results_genes_nominal$VAR)))
  cat("\n")
  
  
  Results_genes_nominal_Block_and_PCHiC_sel<-merge(Results_genes_nominal,
                                                   Block_and_PCHiC_candidates,
                                                   by=c("VAR","ensembl_gene_id"))
  
  cat("Results_genes_nominal_Block_and_PCHiC_sel\n")
  cat(str(Results_genes_nominal_Block_and_PCHiC_sel))
  cat("\n")
  cat(str(unique(Results_genes_nominal_Block_and_PCHiC_sel$VAR)))
  cat("\n")
  
  if(dim(Results_genes_nominal_Block_and_PCHiC_sel)[1] >0)
  {
    
    ##### Correction per cell type -----
    
    
    CT_levels<-unique(as.character(Results_genes_nominal_Block_and_PCHiC_sel$Cell_Type))
    
    cat("CT_levels\n")
    cat(str(CT_levels))
    cat("\n")
    
    list_CT<-list()
    
    for(iteration_CT in 1:length(CT_levels))
    {
      
      CT_levels_sel<-CT_levels[iteration_CT]
      
      cat("-------------------------------------------------------------------------------------------------------------------------------CT------>\t")
      cat(sprintf(as.character(CT_levels_sel)))
      cat("\n")
      
      Results_genes_nominal_Block_and_PCHiC_sel_CT_sel<-Results_genes_nominal_Block_and_PCHiC_sel[which(Results_genes_nominal_Block_and_PCHiC_sel$Cell_Type == CT_levels_sel),]
      
      cat("Results_genes_nominal_Block_and_PCHiC_sel_CT_sel\n")
      cat(str(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel))
      cat("\n")
      
      
      VARS_array<-unique(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel$VAR)
      
      list_VAR<-list()
      
      for(iteration_VARS_array in 1:length(VARS_array))
      {
        
        VAR_sel<-VARS_array[iteration_VARS_array]
        
        cat("-----------VAR------->\t")
        cat(sprintf(as.character(VAR_sel)))
        cat("\n")
        
        Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel<-Results_genes_nominal_Block_and_PCHiC_sel_CT_sel[which(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel$VAR%in%VAR_sel),]
        
        cat("Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel\n")
        cat(str(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel))
        cat("\n")
        cat(str(unique(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel$VAR)))
        cat("\n")
        cat(str(unique(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel$HGNC)))
        cat("\n")
        
        #### SELECT TESTS TO CORRECT ----
        
        
        indx.select<-c(which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "VAR"),
                       which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "ensembl_gene_id"),which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "HGNC"),which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "transcript_id"),
                       which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "coefficient_Genotypes"),
                       which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "pvalue_Genotypes"),
                       which(colnames(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel)== "n_breakdown_string"))
        
        
        
        MT_set<-unique(Results_genes_nominal_Block_and_PCHiC_sel_CT_sel_VAR_sel[,indx.select])
        
        cat("MT_set_0\n")
        cat(str(MT_set))
        cat("\n")
        cat(str(unique(MT_set$ensembl_gene_id)))
        cat("\n")
        
        MT_set$ajusted.pvalue_Genotypes<-p.adjust(MT_set$pvalue_Genotypes, method = "BH")
        MT_set$ajusted.minuslogpvalue_Genotypes<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes),2)
        
        
        
        cat("MT_set_1\n")
        cat(str(MT_set))
        cat("\n")
        
        list_VAR[[iteration_VARS_array]]<-MT_set
        
      }#iteration_VARS_array in 1:length(VARS_array)
      
      if(length(list_VAR) > 0)
      {
        VAR_df = as.data.frame(data.table::rbindlist(list_VAR, fill=T), stringsAsFactors=F)
        
        VAR_df$Cell_Type<-CT_levels_sel
        
        cat("VAR_df\n")
        cat(str(VAR_df))
        cat("\n")
        
        list_CT[[iteration_CT]]<-VAR_df
        
        
      }#length(list_VAR) > 0
    }# iteration_CT in 1:length(CT_levels)
    
    if(length(list_CT) > 0)
    {
      CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
      
      cat("CT_df_0\n")
      cat(str(CT_df))
      cat("\n")
      
      
      check<-CT_df[which(CT_df$ajusted.minuslogpvalue_Genotypes >=1.3),]
      
      cat("check\n")
      cat(str(check))
      cat("\n")
      cat(sprintf(as.character(names(summary((unique(droplevels(interaction(check$VAR,check$HGNC)))))))))
      cat("\n")
      cat(sprintf(as.character(summary((unique(droplevels(interaction(check$VAR,check$HGNC))))))))
      cat("\n")
      
      setwd(out)
      
      if(dim(CT_df)[1] >0)
      {
        saveRDS(file="DTU_Block_and_PCHic_gene.rds", CT_df)  
      }
      
    }#length(list_CT) > 0
    
  }# dim(Results_genes_nominal_Block_and_PCHiC_sel)[1] >0
  
  
}

put_together_CIS_and_Block= function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### Results_DEF ----
  
  Results_DEF<-data.frame()
  
  #### DTU_CIS----
  
  Condition_DEBUG <- 1
  
  setwd(out)
  
  DTU_CIS<-readRDS(file="DTU_CIS_gene.rds")
  
  cat("DTU_CIS_0\n")
  cat(str(DTU_CIS))
  cat("\n")
  cat(str(unique(DTU_CIS$VAR)))
  cat("\n")
  
  DTU_CIS.dt<-data.table(DTU_CIS, key=c("VAR","Cell_Type","ensembl_gene_id"))
  
  DTU_CIS_MAX<-as.data.frame(DTU_CIS.dt[,.SD[which.max(ajusted.minuslogpvalue_Genotypes)], by=key(DTU_CIS.dt)], stringsAsFactor=F)
  
  cat("DTU_CIS_MAX\n")
  cat(str(DTU_CIS_MAX))
  cat("\n")
  
  DTU_CIS_MAX.dt<-data.table(DTU_CIS_MAX, key=c("VAR","Cell_Type"))
  
  DTU_CIS_MAX_collapsed<-as.data.frame(DTU_CIS_MAX.dt[,.N,by=key(DTU_CIS_MAX.dt)], stringsAsFactors=F)
  
  colnames(DTU_CIS_MAX_collapsed)<-c("VAR","Cell_Type","n_tested_CIS")
  
  
  
  
  DTU_CIS_MAX_SIG<-DTU_CIS_MAX[which(DTU_CIS_MAX$ajusted.minuslogpvalue_Genotypes >= 1.3),]
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_CIS_MAX_SIG_1\n")
    cat(str(DTU_CIS_MAX_SIG))
    cat("\n")
  }
  
  
  DTU_CIS_MAX_SIG.dt<-data.table(DTU_CIS_MAX_SIG, key=c("VAR","Cell_Type"))
  
  DTU_CIS_MAX_SIG_collapsed<-as.data.frame(DTU_CIS_MAX_SIG.dt[,.N,by=key(DTU_CIS_MAX_SIG.dt)], stringsAsFactors=F)
  
  colnames(DTU_CIS_MAX_SIG_collapsed)<-c("VAR","Cell_Type","n_SIG_CIS")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_CIS_MAX_SIG_collapsed_0\n")
    cat(str(DTU_CIS_MAX_SIG_collapsed))
    cat("\n")
  }
  
  DTU_CIS_MAX<-merge(DTU_CIS_MAX_collapsed,
               DTU_CIS_MAX_SIG_collapsed,
               by=c("VAR","Cell_Type"),
               all=T)
  
  DTU_CIS_MAX$n_SIG_CIS[is.na(DTU_CIS_MAX$n_SIG_CIS)]<-0
  DTU_CIS_MAX$Perc_SIG_CIS<-round(100*(DTU_CIS_MAX$n_SIG_CIS/DTU_CIS_MAX$n_tested_CIS),1)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_CIS_MAX_2\n")
    cat(str(DTU_CIS_MAX))
    cat("\n")
  }
 
  DTU_CIS<-merge(DTU_CIS,
                     DTU_CIS_MAX,
                     by=c("VAR","Cell_Type"),
                     all=T)
  
  cat("DTU_CIS_1\n")
  cat(str(DTU_CIS))
  cat("\n")
  cat(str(unique(DTU_CIS$VAR)))
  cat("\n")
  
  
  
  
  
  indx.int<-c(which(colnames(DTU_CIS) == "VAR"),which(colnames(DTU_CIS) == "ensembl_gene_id"),which(colnames(DTU_CIS) == "HGNC"),which(colnames(DTU_CIS) == "transcript_id"),which(colnames(DTU_CIS) == "Cell_Type"),
              which(colnames(DTU_CIS) == "ajusted.minuslogpvalue_Genotypes"),which(colnames(DTU_CIS) == "coefficient_Genotypes"),which(colnames(DTU_CIS) == "n_breakdown_string"),
              which(colnames(DTU_CIS) == "n_SIG_CIS"),which(colnames(DTU_CIS) == "n_tested_CIS"),which(colnames(DTU_CIS) == "Perc_SIG_CIS"))
  
  DTU_CIS_subset<-unique(DTU_CIS[,indx.int])
  
  colnames(DTU_CIS_subset)[which(colnames(DTU_CIS_subset) == "ajusted.minuslogpvalue_Genotypes")]<-"CIS_gene_minuslogpvalue"
  colnames(DTU_CIS_subset)[which(colnames(DTU_CIS_subset) == "coefficient_Genotypes")]<-"CIS_Beta"
  colnames(DTU_CIS_subset)[which(colnames(DTU_CIS_subset) == "n_breakdown_string")]<-"CIS_n_breakdown_string"
  
  
  Results_DEF<-rbind(DTU_CIS_subset,Results_DEF)
  
  if(Condition_DEBUG == 1)
  {
    cat("Results_DEF_0\n")
    cat(str(Results_DEF))
    cat("\n")
  }
  
  # quit(status = 1)
  
  #### DTU_Block_and_PCHiC ----
  
  DTU_Block_and_PCHiC<-readRDS(file="DTU_Block_and_PCHic_gene.rds")
  
  cat("DTU_Block_and_PCHiC\n")
  cat(str(DTU_Block_and_PCHiC))
  cat("\n")
  cat(str(unique(DTU_Block_and_PCHiC$VAR)))
  cat("\n")
  
  DTU_Block_and_PCHiC.dt<-data.table(DTU_Block_and_PCHiC, key=c("VAR","Cell_Type","ensembl_gene_id"))
  
  DTU_Block_and_PCHiC_MAX<-as.data.frame(DTU_Block_and_PCHiC.dt[,.SD[which.max(ajusted.minuslogpvalue_Genotypes)], by=key(DTU_Block_and_PCHiC.dt)], stringsAsFactor=F)
  
  cat("DTU_Block_and_PCHiC_MAX\n")
  cat(str(DTU_Block_and_PCHiC_MAX))
  cat("\n")
  
  DTU_Block_and_PCHiC_MAX.dt<-data.table(DTU_Block_and_PCHiC_MAX, key=c("VAR","Cell_Type"))
  
  DTU_Block_and_PCHiC_MAX_collapsed<-as.data.frame(DTU_Block_and_PCHiC_MAX.dt[,.N,by=key(DTU_Block_and_PCHiC_MAX.dt)], stringsAsFactors=F)
  
  colnames(DTU_Block_and_PCHiC_MAX_collapsed)<-c("VAR","Cell_Type","n_tested_Block")
  
  
  
  
  DTU_Block_and_PCHiC_MAX_SIG<-DTU_Block_and_PCHiC_MAX[which(DTU_Block_and_PCHiC_MAX$ajusted.minuslogpvalue_Genotypes >= 1.3),]
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_Block_and_PCHiC_MAX_SIG_1\n")
    cat(str(DTU_Block_and_PCHiC_MAX_SIG))
    cat("\n")
  }
  
  
  DTU_Block_and_PCHiC_MAX_SIG.dt<-data.table(DTU_Block_and_PCHiC_MAX_SIG, key=c("VAR","Cell_Type"))
  
  DTU_Block_and_PCHiC_MAX_SIG_collapsed<-as.data.frame(DTU_Block_and_PCHiC_MAX_SIG.dt[,.N,by=key(DTU_Block_and_PCHiC_MAX_SIG.dt)], stringsAsFactors=F)
  
  colnames(DTU_Block_and_PCHiC_MAX_SIG_collapsed)<-c("VAR","Cell_Type","n_SIG_Block")
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_Block_and_PCHiC_MAX_SIG_collapsed_0\n")
    cat(str(DTU_Block_and_PCHiC_MAX_SIG_collapsed))
    cat("\n")
  }
  
  DTU_Block_and_PCHiC_MAX<-merge(DTU_Block_and_PCHiC_MAX_collapsed,
                     DTU_Block_and_PCHiC_MAX_SIG_collapsed,
                     by=c("VAR","Cell_Type"),
                     all=T)
  
  DTU_Block_and_PCHiC_MAX$n_SIG_Block[is.na(DTU_Block_and_PCHiC_MAX$n_SIG_Block)]<-0
  DTU_Block_and_PCHiC_MAX$Perc_SIG_Block<-round(100*(DTU_Block_and_PCHiC_MAX$n_SIG_Block/DTU_Block_and_PCHiC_MAX$n_tested_Block),1)
  
  
  if(Condition_DEBUG == 1)
  {
    cat("DTU_Block_and_PCHiC_MAX_2\n")
    cat(str(DTU_Block_and_PCHiC_MAX))
    cat("\n")
  }
  
  DTU_Block_and_PCHiC<-merge(DTU_Block_and_PCHiC,
                 DTU_Block_and_PCHiC_MAX,
                 by=c("VAR","Cell_Type"),
                 all=T)
  
  cat("DTU_Block_and_PCHiC_1\n")
  cat(str(DTU_Block_and_PCHiC))
  cat("\n")
  cat(str(unique(DTU_Block_and_PCHiC$VAR)))
  cat("\n")
  
  
  
  
  
  indx.int<-c(which(colnames(DTU_Block_and_PCHiC) == "VAR"),which(colnames(DTU_Block_and_PCHiC) == "ensembl_gene_id"),which(colnames(DTU_Block_and_PCHiC) == "HGNC"),which(colnames(DTU_Block_and_PCHiC) == "transcript_id"),which(colnames(DTU_Block_and_PCHiC) == "Cell_Type"),
              which(colnames(DTU_Block_and_PCHiC) == "ajusted.minuslogpvalue_Genotypes"),which(colnames(DTU_Block_and_PCHiC) == "coefficient_Genotypes"),which(colnames(DTU_Block_and_PCHiC) == "n_breakdown_string"),
              which(colnames(DTU_Block_and_PCHiC) == "n_SIG_Block"),which(colnames(DTU_Block_and_PCHiC) == "n_tested_Block"),which(colnames(DTU_Block_and_PCHiC) == "Perc_SIG_Block"))
  
  DTU_Block_and_PCHiC_subset<-unique(DTU_Block_and_PCHiC[,indx.int])
  
  colnames(DTU_Block_and_PCHiC_subset)[which(colnames(DTU_Block_and_PCHiC_subset) == "ajusted.minuslogpvalue_Genotypes")]<-"Block_PCHiC_gene_minuslogpvalue"
  colnames(DTU_Block_and_PCHiC_subset)[which(colnames(DTU_Block_and_PCHiC_subset) == "coefficient_Genotypes")]<-"Block_PCHiC_Beta"
  colnames(DTU_Block_and_PCHiC_subset)[which(colnames(DTU_Block_and_PCHiC_subset) == "n_breakdown_string")]<-"Block_PCHiC_n_breakdown_string"
  
  
  Results_DEF<-merge(Results_DEF,
                     DTU_Block_and_PCHiC_subset,
                     by=c("VAR","ensembl_gene_id","HGNC","transcript_id","Cell_Type"),
                     all=T)
  
  if(Condition_DEBUG == 1)
  {
    cat("Results_DEF_0\n")
    cat(str(Results_DEF))
    cat("\n")
  }
  
  setwd(out)
  
  write.table(file="DTU_LogRatio_FPKM_results.tsv", Results_DEF, sep="\t", row.names = F,quote = F)
  
  
  
}
printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----




main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--SELECTED_VARS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--VEP_CSQ"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC_info"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_dB"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENES_PER_BLOCKS"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
 
  collate_results_function(opt)
  function_CIS(opt)
  function_Block(opt)
  put_together_CIS_and_Block(opt)
  
}


###########################################################################

system.time( main() )
