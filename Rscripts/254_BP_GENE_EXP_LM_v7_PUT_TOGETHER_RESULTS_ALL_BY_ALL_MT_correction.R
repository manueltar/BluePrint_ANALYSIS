

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
  list_GW<-list()
  list_Block<-list()
  list_CIS<-list()
  
  
 
  
  Condition_DEBUG <- 0
  
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
    
    filename="DE_genome_wide.rds"
    
   
    
    if (file.exists(filename)) {
      DE_genome_wide = readRDS(filename)
      
       # quit(status = 1)
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_genome_wide_0\n")
        cat(str(DE_genome_wide))
        cat("\n")
      }
      
      DE_genome_wide_subset<-unique(DE_genome_wide[,-c(which(colnames(DE_genome_wide) == "ajusted.pvalue_Genotypes"),
                                                   which(colnames(DE_genome_wide) == "ajusted.minuslogpvalue_Genotypes"))])
      
      DE_genome_wide_subset$VAR<-SELECTED_VARS_sel
      if(Condition_DEBUG == 1)
      {
        cat("DE_genome_wide_subset_0\n")
        cat(str(DE_genome_wide_subset))
        cat("\n")
      }
      
      list_GW[[i]]<-DE_genome_wide_subset
      
      # ##########################################
      # quit(status = 1)
      
    }#file.exists(DE_genome_wide.rds)
    
    ### DTU block
    
    setwd(path6)
    
    filename="DE_Block.rds"
    
    
    
    if (file.exists(filename)) {
      DE_BLOCK = readRDS(filename)
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_BLOCK_0\n")
        cat(str(DE_BLOCK))
        cat("\n")
      }
      
      DE_BLOCK_subset<-unique(DE_BLOCK[,-c(which(colnames(DE_BLOCK) == "ajusted.pvalue_Genotypes"),
                                                       which(colnames(DE_BLOCK) == "ajusted.minuslogpvalue_Genotypes"))])
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_BLOCK_0\n")
        cat(str(DE_BLOCK))
        cat("\n")
      }
      
      DE_BLOCK_subset<-unique(DE_BLOCK[,-c(which(colnames(DE_BLOCK) == "ajusted.pvalue_Genotypes"),
                                                       which(colnames(DE_BLOCK) == "ajusted.minuslogpvalue_Genotypes"))])
      
      DE_BLOCK_subset$VAR<-SELECTED_VARS_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_BLOCK_subset_0\n")
        cat(str(DE_BLOCK_subset))
        cat("\n")
      }
      
      list_Block[[i]]<-DE_BLOCK_subset
      
    }
    
    # cat("DE_BLOCK_FINAL\n")
   
    
    ### CIS gene
    setwd(path6)
    
    filename="DE_CIS_gene.rds"
    
    
    
    if (file.exists(filename)) {
      DE_CIS_gene = readRDS(filename)
      
      # cat("DE_CIS_gene\n")
      # cat(str(DE_CIS_gene))
      # cat("\n")
      
      DE_CIS_gene_subset<-unique(DE_CIS_gene[,-c(which(colnames(DE_CIS_gene) == "ajusted.pvalue_Genotypes"),
                                           which(colnames(DE_CIS_gene) == "ajusted.minuslogpvalue_Genotypes"))])
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_0\n")
        cat(str(DE_CIS_gene))
        cat("\n")
      }
      
      DE_CIS_gene_subset<-unique(DE_CIS_gene[,-c(which(colnames(DE_CIS_gene) == "ajusted.pvalue_Genotypes"),
                                           which(colnames(DE_CIS_gene) == "ajusted.minuslogpvalue_Genotypes"))])
      
      DE_CIS_gene_subset$VAR<-SELECTED_VARS_sel
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_subset_0\n")
        cat(str(DE_CIS_gene_subset))
        cat("\n")
      }
      
      list_CIS[[i]]<-DE_CIS_gene_subset
      
      
   
    #  quit(status = 1)
      
      
    }
    
    
  }# i in 1:length(SELECTED_VARS)
  
  # ##################################
  # quit(status = 1)
  
  
  Condition_DEBUG <- 1
  
 
  
  # list_DEF[[3]]<-Results_CIS_subset 
  
  
 
  
  Results_DEF<-data.frame()
  
  if(length(list_CIS) >0)
  {
    
    Results_CIS_nominal = as.data.frame(data.table::rbindlist(list_CIS, fill=T), stringsAsFactors=F)
    
    
    if(Condition_DEBUG == 1)
    {
      cat("Results_CIS_nominal_0\n")
      cat(str(Results_CIS_nominal))
      cat("\n")
      cat(str(unique(Results_CIS_nominal$VAR)))
      cat("\n")
    }
    
    CT_levels<-unique(as.character(Results_CIS_nominal$Cell_Type))
    
    cat("CT_levels\n")
    cat(str(CT_levels))
    cat("\n")
    
    list_CT<-list()
    
    for(iteration_CT in 1:length(CT_levels))
    {
      
      CT_levels_sel<-CT_levels[iteration_CT]
      
      cat("-------------->\t")
      cat(sprintf(as.character(CT_levels_sel)))
      cat("\n")
      
      Results_CIS_nominal_CT_sel<-Results_CIS_nominal[which(Results_CIS_nominal$Cell_Type == CT_levels_sel),]
      
      cat("Results_CIS_nominal_CT_sel\n")
      cat(str(Results_CIS_nominal_CT_sel))
      cat("\n")
      
      
      #### SELECT TESTS TO CORRECT ----
      
      
      indx.select<-c(which(colnames(Results_CIS_nominal_CT_sel)== "VAR"),
                     which(colnames(Results_CIS_nominal_CT_sel)== "ensembl_gene_id"),which(colnames(Results_CIS_nominal_CT_sel)== "HGNC"),which(colnames(Results_CIS_nominal_CT_sel)== "transcript_id"),
                     which(colnames(Results_CIS_nominal_CT_sel)== "coefficient_Genotypes"),
                     which(colnames(Results_CIS_nominal_CT_sel)== "pvalue_Genotypes"),
                     which(colnames(Results_CIS_nominal_CT_sel)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_CIS_nominal_CT_sel[,indx.select])
      
      cat("MT_set_0\n")
      cat(str(MT_set))
      cat("\n")
      cat(str(unique(MT_set$ensembl_gene_id)))
      cat("\n")
      
      MT_set$ajusted.pvalue_Genotypes<-p.adjust(MT_set$pvalue_Genotypes, method = "BH")
      MT_set$ajusted.minuslogpvalue_Genotypes<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes),2)
      MT_set$Cell_Type<-CT_levels_sel
      
      
      cat("MT_set_1\n")
      cat(str(MT_set))
      cat("\n")
      
      list_CT[[iteration_CT]]<-MT_set
      
      
    }# iteration_CT in 1:length(CT_levels)
    
    if(length(list_CT) > 0)
    {
      CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
      
      cat("CT_df_0\n")
      cat(str(CT_df))
      cat("\n")
      
      CT_df.dt<-data.table(CT_df, key=c("VAR","Cell_Type"))
      
      CT_df_collapsed<-as.data.frame(CT_df.dt[,.N,by=key(CT_df.dt)], stringsAsFactors=F)
      
      colnames(CT_df_collapsed)<-c("VAR","Cell_Type","n_tested_CIS")
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_collapsed_0\n")
        cat(str(CT_df_collapsed))
        cat("\n")
      }
      
      
      
      DE_CIS_gene<-CT_df[which(CT_df$ajusted.minuslogpvalue_Genotypes >= 1.3),]
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_1\n")
        cat(str(DE_CIS_gene))
        cat("\n")
      }
      
      
      DE_CIS_gene.dt<-data.table(DE_CIS_gene, key=c("VAR","Cell_Type"))
      
      DE_CIS_gene_collapsed<-as.data.frame(DE_CIS_gene.dt[,.N,by=key(DE_CIS_gene.dt)], stringsAsFactors=F)
      
      colnames(DE_CIS_gene_collapsed)<-c("VAR","Cell_Type","n_SIG_CIS")
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_CIS_gene_collapsed_0\n")
        cat(str(DE_CIS_gene_collapsed))
        cat("\n")
      }
      
      CT_df<-merge(CT_df,
                    CT_df_collapsed,
                    by=c("VAR","Cell_Type"),
                    all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_2\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      CT_df<-merge(CT_df,
                    DE_CIS_gene_collapsed,
                    by=c("VAR","Cell_Type"),
                    all=T)
      
      CT_df$n_SIG_CIS[is.na(CT_df$n_SIG_CIS)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_3\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      CT_df$Perc_SIG_CIS<-round(100*(CT_df$n_SIG_CIS/CT_df$n_tested_CIS),1)
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_4\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      indx.int<-c(which(colnames(CT_df) == "VAR"),which(colnames(CT_df) == "ensembl_gene_id"),which(colnames(CT_df) == "HGNC"),which(colnames(CT_df) == "Cell_Type"),
                  which(colnames(CT_df) == "ajusted.minuslogpvalue_Genotypes"),which(colnames(CT_df) == "coefficient_Genotypes"),which(colnames(CT_df) == "n_breakdown_string"),
                  which(colnames(CT_df) == "n_SIG_CIS"),which(colnames(CT_df) == "n_tested_CIS"),which(colnames(CT_df) == "Perc_SIG_CIS"))
      
      CT_df_subset<-unique(CT_df[,indx.int])
      
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "ajusted.minuslogpvalue_Genotypes")]<-"CIS_gene_minuslogpvalue"
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "coefficient_Genotypes")]<-"CIS_Beta"
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "n_breakdown_string")]<-"CIS_n_breakdown_string"
      
      
      Results_DEF<-rbind(CT_df_subset,Results_DEF)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_0\n")
        cat(str(Results_DEF))
        cat("\n")
      }
      
     
      
    }#length(list_CT) > 0
    
    
    # ####################################################################################################
    # quit(status = 1)
    
  }# length(list_CIS) >0
 
  
  if(length(list_Block) >0)
  {
    
    Results_Block_nominal = as.data.frame(data.table::rbindlist(list_Block, fill=T), stringsAsFactors=F)
    
    
    
    if(Condition_DEBUG == 1)
    {
      cat("Results_Block_nominal_0\n")
      cat(str(Results_Block_nominal))
      cat("\n")
      cat(str(unique(Results_Block_nominal$VAR)))
      cat("\n")
    }
    
    
    CT_levels<-unique(as.character(Results_Block_nominal$Cell_Type))
    
    cat("CT_levels\n")
    cat(str(CT_levels))
    cat("\n")
    
    list_CT<-list()
    
    for(iteration_CT in 1:length(CT_levels))
    {
      
      CT_levels_sel<-CT_levels[iteration_CT]
      
      cat("-------------->\t")
      cat(sprintf(as.character(CT_levels_sel)))
      cat("\n")
      
      Results_Block_nominal_CT_sel<-Results_Block_nominal[which(Results_Block_nominal$Cell_Type == CT_levels_sel),]
      
      cat("Results_Block_nominal_CT_sel\n")
      cat(str(Results_Block_nominal_CT_sel))
      cat("\n")
      
      
      #### SELECT TESTS TO CORRECT ----
      
      
      indx.select<-c(which(colnames(Results_Block_nominal_CT_sel)== "VAR"),
                     which(colnames(Results_Block_nominal_CT_sel)== "ensembl_gene_id"),which(colnames(Results_Block_nominal_CT_sel)== "HGNC"),which(colnames(Results_Block_nominal_CT_sel)== "transcript_id"),
                     which(colnames(Results_Block_nominal_CT_sel)== "coefficient_Genotypes"),
                     which(colnames(Results_Block_nominal_CT_sel)== "pvalue_Genotypes"),
                     which(colnames(Results_Block_nominal_CT_sel)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_Block_nominal_CT_sel[,indx.select])
      
      cat("MT_set_0\n")
      cat(str(MT_set))
      cat("\n")
      cat(str(unique(MT_set$ensembl_gene_id)))
      cat("\n")
      
      MT_set$ajusted.pvalue_Genotypes<-p.adjust(MT_set$pvalue_Genotypes, method = "BH")
      MT_set$ajusted.minuslogpvalue_Genotypes<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes),2)
      MT_set$Cell_Type<-CT_levels_sel
      
      
      cat("MT_set_1\n")
      cat(str(MT_set))
      cat("\n")
      
      list_CT[[iteration_CT]]<-MT_set
      
      
    }# iteration_CT in 1:length(CT_levels)
    
    if(length(list_CT) > 0)
    {
      CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
      
      cat("CT_df_0\n")
      cat(str(CT_df))
      cat("\n")
      
      CT_df.dt<-data.table(CT_df, key=c("VAR","Cell_Type"))
      
      CT_df_collapsed<-as.data.frame(CT_df.dt[,.N,by=key(CT_df.dt)], stringsAsFactors=F)
      
      colnames(CT_df_collapsed)<-c("VAR","Cell_Type","n_tested_Block")
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_collapsed_0\n")
        cat(str(CT_df_collapsed))
        cat("\n")
      }
      
      
      
      DE_Block<-CT_df[which(CT_df$ajusted.minuslogpvalue_Genotypes >= 1.3),]
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_Block_1\n")
        cat(str(DE_Block))
        cat("\n")
      }
      
      
      DE_Block.dt<-data.table(DE_Block, key=c("VAR","Cell_Type"))
      
      DE_Block_collapsed<-as.data.frame(DE_Block.dt[,.N,by=key(DE_Block.dt)], stringsAsFactors=F)
      
      colnames(DE_Block_collapsed)<-c("VAR","Cell_Type","n_SIG_Block")
      
      
      if(Condition_DEBUG == 1)
      {
        cat("DE_Block_collapsed_0\n")
        cat(str(DE_Block_collapsed))
        cat("\n")
      }
      
      CT_df<-merge(CT_df,
                   CT_df_collapsed,
                   by=c("VAR","Cell_Type"),
                   all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_2\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      CT_df<-merge(CT_df,
                   DE_Block_collapsed,
                   by=c("VAR","Cell_Type"),
                   all=T)
      
      CT_df$n_SIG_Block[is.na(CT_df$n_SIG_Block)]<-0
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_3\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      CT_df$Perc_SIG_Block<-round(100*(CT_df$n_SIG_Block/CT_df$n_tested_Block),1)
      
      if(Condition_DEBUG == 1)
      {
        cat("CT_df_4\n")
        cat(str(CT_df))
        cat("\n")
      }
      
      indx.int<-c(which(colnames(CT_df) == "VAR"),which(colnames(CT_df) == "ensembl_gene_id"),which(colnames(CT_df) == "HGNC"),which(colnames(CT_df) == "Cell_Type"),
                  which(colnames(CT_df) == "ajusted.minuslogpvalue_Genotypes"),which(colnames(CT_df) == "coefficient_Genotypes"),
                  which(colnames(CT_df) == "n_breakdown_string"),
                  which(colnames(CT_df) == "n_SIG_Block"),which(colnames(CT_df) == "n_tested_Block"),which(colnames(CT_df) == "Perc_SIG_Block"))
      
      CT_df_subset<-unique(CT_df[,indx.int])
      
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "ajusted.minuslogpvalue_Genotypes")]<-"Block_PCHiC_minuslogpvalue"
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "coefficient_Genotypes")]<-"Block_PCHiC_Beta"
      colnames(CT_df_subset)[which(colnames(CT_df_subset) == "n_breakdown_string")]<-"Block_PCHiC_n_breakdown_string"
      
      
      Results_DEF<-merge(Results_DEF,
                         CT_df_subset,
                         by=c("VAR","ensembl_gene_id","HGNC","Cell_Type"),
                         all=T)
      
      if(Condition_DEBUG == 1)
      {
        cat("Results_DEF_1\n")
        cat(str(Results_DEF))
        cat("\n")
        
        # ########################################################
        # quit(status = 1)
      }
      
      
      
    }#length(list_CT) > 0
    
  }# length(list_Block) >0
  
 
  
  # if(length(list_GW) >0)
  # {
  #   
  #   Results_GW_nominal = as.data.frame(data.table::rbindlist(list_GW, fill=T), stringsAsFactors=F)
  #   
  #   
  #   
  #   if(Condition_DEBUG == 1)
  #   {
  #     cat("Results_GW_nominal_0\n")
  #     cat(str(Results_GW_nominal))
  #     cat("\n")
  #     cat(str(unique(Results_GW_nominal$VAR)))
  #     cat("\n")
  #   }
  #   
  #   
  #   CT_levels<-unique(as.character(Results_GW_nominal$Cell_Type))
  #   
  #   cat("CT_levels\n")
  #   cat(str(CT_levels))
  #   cat("\n")
  #   
  #   list_CT<-list()
  #   
  #   for(iteration_CT in 1:length(CT_levels))
  #   {
  #     
  #     CT_levels_sel<-CT_levels[iteration_CT]
  #     
  #     cat("-------------->\t")
  #     cat(sprintf(as.character(CT_levels_sel)))
  #     cat("\n")
  #     
  #     Results_GW_nominal_CT_sel<-Results_GW_nominal[which(Results_GW_nominal$Cell_Type == CT_levels_sel),]
  #     
  #     cat("Results_GW_nominal_CT_sel\n")
  #     cat(str(Results_GW_nominal_CT_sel))
  #     cat("\n")
  #     
  #     
  #     #### SELECT TESTS TO CORRECT ----
  #     
  #     
  #     indx.select<-c(which(colnames(Results_GW_nominal_CT_sel)== "VAR"),
  #                    which(colnames(Results_GW_nominal_CT_sel)== "ensembl_gene_id"),which(colnames(Results_GW_nominal_CT_sel)== "HGNC"),which(colnames(Results_GW_nominal_CT_sel)== "transcript_id"),
  #                    which(colnames(Results_GW_nominal_CT_sel)== "coefficient_Genotypes"),
  #                    which(colnames(Results_GW_nominal_CT_sel)== "pvalue_Genotypes"))
  #     
  #     
  #     
  #     MT_set<-unique(Results_GW_nominal_CT_sel[,indx.select])
  #     
  #     cat("MT_set_0\n")
  #     cat(str(MT_set))
  #     cat("\n")
  #     cat(str(unique(MT_set$ensembl_gene_id)))
  #     cat("\n")
  #     
  #     MT_set$ajusted.pvalue_Genotypes<-p.adjust(MT_set$pvalue_Genotypes, method = "BH")
  #     MT_set$ajusted.minuslogpvalue_Genotypes<-round(-1*log10(MT_set$ajusted.pvalue_Genotypes),2)
  #     MT_set$Cell_Type<-CT_levels_sel
  #     
  #     
  #     cat("MT_set_1\n")
  #     cat(str(MT_set))
  #     cat("\n")
  #     
  #     list_CT[[iteration_CT]]<-MT_set
  #     
  #     
  #   }# iteration_CT in 1:length(CT_levels)
  #   
  #   if(length(list_CT) > 0)
  #   {
  #     CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
  #     
  #     cat("CT_df_0\n")
  #     cat(str(CT_df))
  #     cat("\n")
  #     
  #     CT_df.dt<-data.table(CT_df, key=c("VAR","Cell_Type"))
  #     
  #     CT_df_collapsed<-as.data.frame(CT_df.dt[,.N,by=key(CT_df.dt)], stringsAsFactors=F)
  #     
  #     colnames(CT_df_collapsed)<-c("VAR","Cell_Type","n_tested_GW")
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("CT_df_collapsed_0\n")
  #       cat(str(CT_df_collapsed))
  #       cat("\n")
  #     }
  #     
  #     
  #     
  #     DE_GW<-CT_df[which(CT_df$ajusted.minuslogpvalue_Genotypes >= 1.3),]
  #     
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("DE_GW_1\n")
  #       cat(str(DE_GW))
  #       cat("\n")
  #     }
  #     
  #     
  #     DE_GW.dt<-data.table(DE_GW, key=c("VAR","Cell_Type"))
  #     
  #     DE_GW_collapsed<-as.data.frame(DE_GW.dt[,.N,by=key(DE_GW.dt)], stringsAsFactors=F)
  #     
  #     colnames(DE_GW_collapsed)<-c("VAR","Cell_Type","n_SIG_GW")
  #     
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("DE_GW_collapsed_0\n")
  #       cat(str(DE_GW_collapsed))
  #       cat("\n")
  #     }
  #     
  #     CT_df<-merge(CT_df,
  #                  CT_df_collapsed,
  #                  by=c("VAR","Cell_Type"),
  #                  all=T)
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("CT_df_2\n")
  #       cat(str(CT_df))
  #       cat("\n")
  #     }
  #     
  #     CT_df<-merge(CT_df,
  #                  DE_GW_collapsed,
  #                  by=c("VAR","Cell_Type"),
  #                  all=T)
  #     
  #     CT_df$n_SIG_GW[is.na(CT_df$n_SIG_GW)]<-0
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("CT_df_3\n")
  #       cat(str(CT_df))
  #       cat("\n")
  #     }
  #     
  #     CT_df$Perc_SIG_GW<-round(100*(CT_df$n_SIG_GW/CT_df$n_tested_GW),1)
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("CT_df_4\n")
  #       cat(str(CT_df))
  #       cat("\n")
  #     }
  #     
  #     indx.int<-c(which(colnames(CT_df) == "VAR"),which(colnames(CT_df) == "ensembl_gene_id"),which(colnames(CT_df) == "HGNC"),which(colnames(CT_df) == "Cell_Type"),
  #                 which(colnames(CT_df) == "ajusted.minuslogpvalue_Genotypes"),which(colnames(CT_df) == "coefficient_Genotypes"),
  #                 which(colnames(CT_df) == "n_SIG_GW"),which(colnames(CT_df) == "n_tested_GW"),which(colnames(CT_df) == "Perc_SIG_GW"))
  #     
  #     CT_df_subset<-unique(CT_df[,indx.int])
  #     
  #     colnames(CT_df_subset)[which(colnames(CT_df_subset) == "ajusted.minuslogpvalue_Genotypes")]<-"Genome_wide_minuslogpvalue"
  #     colnames(CT_df_subset)[which(colnames(CT_df_subset) == "coefficient_Genotypes")]<-"Genome_wide_Beta"
  #     
  #     Results_DEF<-merge(Results_DEF,
  #                        CT_df_subset,
  #                        by=c("VAR","ensembl_gene_id","HGNC","Cell_Type"),
  #                        all=T)
  #     
  #     if(Condition_DEBUG == 1)
  #     {
  #       cat("Results_DEF_GW\n")
  #       cat(str(Results_DEF))
  #       cat("\n")
  #       
  #       # ########################################################
  #       # quit(status = 1)
  #     }
  #   }#length(list_CT) > 0
  # }# length(list_GW) >0
 
  if(dim(Results_DEF)[1] >0)
  {
    
    cat("Results_DEF_BORN\n")
    cat(str(Results_DEF))
    cat("\n")
    cat(str(unique(Results_DEF$VAR)))
    cat("\n")
    
    
    # quit(status = 1)
    
    SIG_CIS<-Results_DEF[which(Results_DEF$CIS_gene_minuslogpvalue >= 1.3),]
    
    cat("SIG_CIS_0\n")
    cat(str(SIG_CIS))
    cat("\n")
    cat(str(unique(SIG_CIS$VAR)))
    cat("\n")
    
    SIG_CIS.dt<-data.table(SIG_CIS, key=c("VAR","Cell_Type"))
    
    SIG_CIS_collapsed<-data.frame(SIG_CIS.dt[,.(strings_ENSG_CIS=paste(ensembl_gene_id, collapse = ";"),
                                                strings_HGNC_CIS=paste(HGNC, collapse = ";"),
                                                strings_CIS_gene_minuslogpvalue=paste(CIS_gene_minuslogpvalue, collapse = ";"),
                                                n_SIG_CIS=n_SIG_CIS,
                                                n_tested_CIS=n_tested_CIS,
                                                Perc_SIG_CIS=Perc_SIG_CIS),
                                             by=key(SIG_CIS.dt)], stringsAsFactors=F)
    
    
    
    
    cat("SIG_CIS_collapsed_0\n")
    cat(str(SIG_CIS_collapsed))
    cat("\n")
    cat(str(unique(SIG_CIS_collapsed$VAR)))
    cat("\n")
    
    SIG_BLOCK<-Results_DEF[which(Results_DEF$Block_PCHiC_minuslogpvalue >= 1.3),]
    
    cat("SIG_BLOCK_0\n")
    cat(str(SIG_BLOCK))
    cat("\n")
    cat(str(unique(SIG_BLOCK$VAR)))
    cat("\n")
    
    SIG_BLOCK.dt<-data.table(SIG_BLOCK,  key=c("VAR","Cell_Type"))
    
    SIG_BLOCK_collapsed<-data.frame(SIG_BLOCK.dt[,.(strings_ENSG_Block=paste(ensembl_gene_id, collapse = ";"),
                                                    strings_HGNC_Block=paste(HGNC, collapse = ";"),
                                                    strings_Block_PCHiC_minuslogpvalue=paste(Block_PCHiC_minuslogpvalue, collapse = ";"),
                                                    n_SIG_Block=n_SIG_Block,
                                                    n_tested_Block=n_tested_Block,
                                                    Perc_SIG_Block=Perc_SIG_Block),
                                                 by=key(SIG_BLOCK.dt)], stringsAsFactors=F)
    
    
    
    cat("SIG_BLOCK_collapsed_0\n")
    cat(str(SIG_BLOCK_collapsed))
    cat("\n")
    cat(str(unique(SIG_BLOCK_collapsed$VAR)))
    cat("\n")
    
    # SIG_GW<-Results_DEF[which(Results_DEF$Genome_wide_minuslogpvalue >= 1.3),]
    # 
    # cat("SIG_GW_0\n")
    # cat(str(SIG_GW))
    # cat("\n")
    # cat(str(unique(SIG_GW$VAR)))
    # cat("\n")
    # 
    # SIG_GW.dt<-data.table(SIG_GW,  key=c("VAR","Cell_Type"))
    # 
    # SIG_GW_collapsed<-data.frame(SIG_GW.dt[,.(strings_ENSG_GW=paste(ensembl_gene_id, collapse = ";"),
    #                                           strings_HGNC_GW=paste(HGNC, collapse = ";"),
    #                                           strings_Genome_wide_minuslogpvalue=paste(Genome_wide_minuslogpvalue, collapse = ";"),
    #                                           n_SIG_GW=n_SIG_GW,
    #                                           n_tested_GW=n_tested_GW,
    #                                           Perc_SIG_GW=Perc_SIG_GW),
    #                                        by=key(SIG_GW.dt)], stringsAsFactors=F)
    # 
    # 
    # cat("SIG_GW_collapsed_0\n")
    # cat(str(SIG_GW_collapsed))
    # cat("\n")
    # cat(str(unique(SIG_GW_collapsed$VAR)))
    # cat("\n")
    
    DEF_SIG<-unique(merge(SIG_CIS_collapsed,
                          SIG_BLOCK_collapsed,
                          by=c("VAR","Cell_Type"),
                          all=T))
    
    #"strings_ENSG","strings_HGNC"
    
    # DEF_SIG<-unique(merge(DEF_SIG,
    #                       SIG_GW_collapsed,
    #                       by=c("VAR","Cell_Type"),
    #                       all=T))
    
    # ,"strings_ENSG","strings_HGNC"
    
    cat("DEF_SIG_0\n")
    cat(str(DEF_SIG))
    cat("\n")
    cat(str(unique(DEF_SIG$VAR)))
    cat("\n")
  
    setwd(out)
    
    cat("Results_DEF_PRE_Z_score_0\n")
    cat(str(Results_DEF))
    cat("\n")
    cat(str(unique(Results_DEF$VAR)))
    cat("\n")
    
    
    ## Explore per CT effect size values
    
    explore.dt<-data.table(Results_DEF, key="Cell_Type")
    
    exploration_df<-as.data.frame(explore.dt[,.(MIN_CIS=min(CIS_Beta[!is.na(CIS_Beta)]),
                                                Q1_CIS=summary(CIS_Beta[!is.na(CIS_Beta)])[2],
                                                MEDIAN_CIS=summary(CIS_Beta[!is.na(CIS_Beta)])[3],
                                                Q3_CIS=summary(CIS_Beta[!is.na(CIS_Beta)])[5],
                                                MAX_CIS=summary(CIS_Beta[!is.na(CIS_Beta)])[6],
                                                mean_CIS=mean(CIS_Beta[!is.na(CIS_Beta)]),
                                                sd_CIS=sd(CIS_Beta[!is.na(CIS_Beta)]),
                                                MIN_Block_PCHiC=min(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)]),
                                                Q1_Block_PCHiC=summary(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)])[2],
                                                MEDIAN_Block_PCHiC=summary(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)])[3],
                                                Q3_Block_PCHiC=summary(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)])[5],
                                                MAX_Block_PCHiC=summary(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)])[6],
                                                mean_Block_PCHiC=mean(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)]),
                                                sd_Block_PCHiC=sd(Block_PCHiC_Beta[!is.na(Block_PCHiC_Beta)])
                                                ),
                                             by=key(explore.dt)], stringsAsFactors=F)
    
    
    
    cat("exploration_df_PRE_Z_score\n")
    cat(str(exploration_df))
    cat("\n")
    
    exploration_subset_df<-exploration_df[,c(which(colnames(exploration_df) == "Cell_Type"),
                                             which(colnames(exploration_df) == "mean_CIS"),
                                             which(colnames(exploration_df) == "sd_CIS"),
                                             which(colnames(exploration_df) == "mean_Block_PCHiC"),
                                             which(colnames(exploration_df) == "sd_Block_PCHiC"))]
    
    cat("exploration_subset_df_PRE_Z_score\n")
    cat(str(exploration_subset_df))
    cat("\n")
    
    Results_DEF<-merge(Results_DEF,
                       exploration_subset_df,
                       by="Cell_Type",
                       all=T)
    
    
    
    cat("Results_DEF_PRE_Z_score_1\n")
    cat(str(Results_DEF))
    cat("\n")
    cat(str(unique(Results_DEF$VAR)))
    cat("\n")
    
    Results_DEF$CIS_Beta_Z_score<-((Results_DEF$CIS_Beta-Results_DEF$mean_CIS)/Results_DEF$sd_CIS)
    Results_DEF$Block_PCHiC_Beta_Z_score<-((Results_DEF$Block_PCHiC_Beta-Results_DEF$mean_Block_PCHiC)/Results_DEF$sd_Block_PCHiC)
   
    
    cat("Results_DEF_FINAL_POST_Z\n")
    cat(str(Results_DEF))
    cat("\n")
    cat(str(unique(Results_DEF$VAR)))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_DEF$CIS_Beta_Z_score)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_DEF$CIS_Beta_Z_score))))
    cat("\n")
    cat(sprintf(as.character(names(summary(Results_DEF$Block_PCHiC_Beta_Z_score)))))
    cat("\n")
    cat(sprintf(as.character(summary(Results_DEF$Block_PCHiC_Beta_Z_score))))
    cat("\n")
   
    # ###############################
    # quit(status = 1)
    # 
    
    write.table(file="ALL_BY_ALL_DE_LM_FPKM_results.tsv", Results_DEF, sep="\t", row.names = F,quote = F)
    write.table(file="ALL_BY_ALL_DE_LM_FPKM_results_SIG_collapsed.tsv", DEF_SIG, sep="\t", row.names = F,quote = F)
    
    
  }# dim(Results_DEF)[1] >0
  
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

  
}


###########################################################################

system.time( main() )
