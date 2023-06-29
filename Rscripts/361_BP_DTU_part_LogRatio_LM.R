

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


Filter_irrelevant_transcripts = function (option_list)
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
  cat(str(SELECTED_VARS))
  cat("\n")
  
  
  
  
  #### Read transposed expression ----
  
  if(length(SELECTED_VARS) >0)
  {
    
    if (file.exists(opt$BP_transcript_EXP)){
      
      BP_transcript_EXP<-readRDS(file=opt$BP_transcript_EXP)
      
      # as.data.frame(fread(file=opt$BP_transcript_EXP, sep=",", header=T) , stringsAsFactors=F)
      
      
      cat("BP_transcript_EXP\n")
      cat(str(BP_transcript_EXP))
      cat("\n")
      cat(sprintf(as.character(names(summary(2^BP_transcript_EXP$value)))))
      cat("\n")
      cat(sprintf(as.character(summary(2^BP_transcript_EXP$value))))
      cat("\n")
      
     
      
      ENST_array<-unique(BP_transcript_EXP$transcript_id)
      
      cat("ENST_array\n")
      cat(str(ENST_array))
      cat("\n")
      
      ENSG_array<-unique(BP_transcript_EXP$ensembl_gene_id)
      
      cat("ENSG_array_0\n")
      cat(str(ENSG_array))
      cat("\n")
      
      
      ### Read Transcripts_table----
      
      
      Transcripts_table<-as.data.frame(fread(file=opt$Transcripts_table, sep="\t", header=F) , stringsAsFactors=F)
      
      colnames(Transcripts_table)<-c("ensembl_gene_id","HGNC","transcript_id","transcript_name")
      
      cat("Transcripts_table\n")
      cat(str(Transcripts_table))
      cat("\n")
      
      
      #### MASTER LOOP ----
      
      
      list_ACCEPTED_ENST<-list()
      
      
      
      
     
      
      for(i in 1:length(ENSG_array))
      {
        Condition_DEBUG <- 0
        
        ENSG_array_sel<-ENSG_array[i]
        
        cat("----------------------------------------------------------------------------------------------->\t")
        cat(sprintf(as.character(i)))
        cat("\t")
        cat(sprintf(as.character(ENSG_array_sel)))
        cat("\n")
        
        TRANSCRIPTS_table_sel<-Transcripts_table[which(Transcripts_table$ensembl_gene_id == ENSG_array_sel),]
        
        
        if(Condition_DEBUG == 1)
        {
          cat("TRANSCRIPTS_table_sel\n")
          cat(str(TRANSCRIPTS_table_sel))
          cat("\n")
        }
        
        BP_transcript_EXP_sel<-BP_transcript_EXP[which(BP_transcript_EXP$transcript_id%in%TRANSCRIPTS_table_sel$transcript_id),]
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_transcript_EXP_sel_0\n")
          cat(str(BP_transcript_EXP_sel))
          cat("\n")
        }
        
        FLAG_genotypes_diversity<-length(unique(levels(droplevels(BP_transcript_EXP_sel$Genotype))))
        
        if(Condition_DEBUG == 1)
        {
          cat("FLAG_genotypes_diversity\n")
          cat(str(FLAG_genotypes_diversity))
          cat("\n")
        }
        
        if(FLAG_genotypes_diversity >1)
        {
          n_transcripts_per_gene<-length(unique(BP_transcript_EXP_sel$transcript_id))
          
          if(Condition_DEBUG == 1)
          {
            cat("n_transcripts_per_gene\n")
            cat(str(n_transcripts_per_gene))
            cat("\n")
          }
          
          if(n_transcripts_per_gene > 1)
          {
            
            BP_genes_DEF_sel.m<-unique(BP_transcript_EXP_sel[,c(which(colnames(BP_transcript_EXP_sel) == "ensembl_gene_id"),
                                                                which(colnames(BP_transcript_EXP_sel) == "transcript_id"),
                                                                which(colnames(BP_transcript_EXP_sel) == "BP_ID"),
                                                                  which(colnames(BP_transcript_EXP_sel) == "HGNC"),
                                                                  which(colnames(BP_transcript_EXP_sel) == "value"),
                                                                which(colnames(BP_transcript_EXP_sel) == "Cell_Type"),
                                                                which(colnames(BP_transcript_EXP_sel) == "Genotype"))])
            HGNC_sel<-unique(BP_genes_DEF_sel.m$HGNC)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_genes_DEF_sel.m_0\n")
              cat(str(BP_genes_DEF_sel.m))
              cat("\n")
            }
            
            
            
            BP_genes_DEF_sel.m$FPKM<-2^(BP_genes_DEF_sel.m$value)
            
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_genes_DEF_sel.m_1\n")
              cat(str(BP_genes_DEF_sel.m))
              cat("\n")
              
              cat("ensembl_gene_id\n")
              cat(str(unique(BP_genes_DEF_sel.m$ensembl_gene_id)))
              cat("\n")
              # ########################################
              # quit(status = 1)
            }
           
           
            
            cat("transcript_id\n")
            cat(str(unique(BP_genes_DEF_sel.m$transcript_id)))
            cat("\n") 
            
            BP_genes_DEF_sel.m_HET<-droplevels(BP_genes_DEF_sel.m[which(BP_genes_DEF_sel.m$Genotype != "HOM"),])
            
            
            if(Condition_DEBUG == 1)
            {
              
              cat("BP_genes_DEF_sel.m_HET_0\n")
              cat(str(BP_genes_DEF_sel.m_HET))
              cat("\n")
              cat(sprintf(as.character(names(summary(BP_genes_DEF_sel.m_HET$value)))))
              cat("\n")
              cat(sprintf(as.character(summary(BP_genes_DEF_sel.m_HET$value))))
              cat("\n")
              # quit(status = 1)
              
            }
           

            
            
            BP_genes_DEF_sel.m_HET.dt<-data.table(BP_genes_DEF_sel.m_HET,
                                                               key=c("Genotype","transcript_id","Cell_Type"))
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_genes_DEF_sel.m_HET.dt_2\n")
              cat(str(BP_genes_DEF_sel.m_HET.dt))
              cat("\n")
            }
            
           
            
            ### calculate median FPKM per genotype and cell type
            
            Summary_table<-as.data.frame(BP_genes_DEF_sel.m_HET.dt[, .(median=round(as.numeric(summary(FPKM)[3]),3)),
                                                                                by=key(BP_genes_DEF_sel.m_HET.dt)],stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_0\n")
              cat(str(Summary_table))
              cat("\n")
              
              #  quit(status = 1)
            }
            
            Summary_table.dt<-data.table(Summary_table,
                                         key=c("Genotype","Cell_Type"))
            
            #### calculate total GENE EXP per genotype as the sum of the medians
            
            Summary_table_TOTAL<-as.data.frame(Summary_table.dt[, .(TOTAL_GENE_EXP_median=sum(median)),
                                                                by=key(Summary_table.dt)],stringsAsFactors=F)
            
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_TOTAL_0\n")
              cat(str(Summary_table_TOTAL))
              cat("\n")
              
              #  quit(status = 1)
            }
            
            Summary_table<-merge(Summary_table,
                                 Summary_table_TOTAL,
                                 by=c("Genotype","Cell_Type"))
            
            #### Calculate Transcript Ratio
            
            Summary_table$Transcript_Ratio<-(Summary_table$median/Summary_table$TOTAL_GENE_EXP_median)
            
            
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_1\n")
              cat(str(Summary_table))
              cat("\n")
              
              # quit(status = 1)
            }
            
            
            Summary_table.dt<-data.table(Summary_table,
                                         key=c("transcript_id"))
            
            #### Obtain the max transcript ratio per transcript from all the genotypes
            
            
            Summary_table_FILTERED<-as.data.frame(Summary_table.dt[,.SD[which.max(Transcript_Ratio)],
                                                                   by=key(Summary_table.dt)],stringsAsFactors=F)
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_FILTERED_0\n")
              cat(str(Summary_table_FILTERED))
              cat("\n")
              
            }
            
            ### Filter if it is less than 10%
            
            Summary_table_FILTERED<-Summary_table_FILTERED[which(Summary_table_FILTERED$Transcript_Ratio >= 0.1),]
            
            
            if(Condition_DEBUG == 1)
            {
              
              cat("Summary_table_FILTERED_1\n")
              cat(str(Summary_table_FILTERED))
              cat("\n")
            }
            
            if(dim(Summary_table_FILTERED)[1] >1)
            {
              # INTERVAL_isoform_EXP_Filtered_sel_restricted.m_FILTERED<-INTERVAL_isoform_EXP_Filtered_sel_restricted.m[which(INTERVAL_isoform_EXP_Filtered_sel_restricted.m$transcript_id%in%Summary_table_FILTERED$transcript_id),]
              
              PASSED_ENST<-unique(as.character(Summary_table_FILTERED$transcript_id))
              
              cat("PASSED_ENST_0\n")
              cat(str(PASSED_ENST))
              cat("\n")
             
              
              list_ACCEPTED_ENST[[i]]<-PASSED_ENST
              
              
            }else{
              
              # After filtering the TR > 10% int at least 1 genotype the gene only has 1 isoform
            }#dim(Summary_table_FILTERED)[1] >1
          }#n_transcripts_per_gene > 1
        }#FLAG_genotypes_diversity >1
      }#i ENSG_array
      
      
      
      
      if(length(list_ACCEPTED_ENST) >0)
      {
        if(Condition_DEBUG == 1)
        {
          cat("list_ACCEPTED_ENST\n")
          cat(str(list_ACCEPTED_ENST))
          cat("\n")
        }
        
        FINAL_transcripts = unique(unlist(list_ACCEPTED_ENST))
        
        
        if(Condition_DEBUG == 1)
        {
          cat("FINAL_transcripts_0\n")
          cat(str(FINAL_transcripts))
          cat("\n")
        }
        
        FINAL_transcripts<-FINAL_transcripts[!is.null(FINAL_transcripts)]
        
        if(Condition_DEBUG == 1)
        {
          cat("FINAL_transcripts_1\n")
          cat(str(FINAL_transcripts))
          cat("\n")
        }
        
        Condition_DEBUG <- 1
        
        Transcripts_table_FINAL<-unique(Transcripts_table[which(Transcripts_table$transcript_id%in%FINAL_transcripts),])
        
        
        if(Condition_DEBUG == 1)
        {
          cat("Transcripts_table_FINAL:1\n")
          cat(str(Transcripts_table_FINAL))
          cat("\n")
        }
        
        
        path6<-paste(out,SELECTED_VARS,'/', sep='')
        
        # cat("path6\n")
        # cat(sprintf(as.character(path6)))
        # cat("\n")
        
        
        if (file.exists(path6)){
          
          
          
          
        } else {
          dir.create(file.path(path6))
          
        }
        
        setwd(path6)
        
        
        saveRDS(file=paste("DTU_PASS_Transcripts_",SELECTED_VARS,'.rds', sep=''), Transcripts_table_FINAL)
        
        
        
        # ###################################################################################################################
        # quit(status = 1)
        
        
        
      }# length(list_ACCEPTED_ENST) >0
        
        
      
      
    }#file.exists(opt$BP_transcript_EXP)
    
  }# length(SELECTED_VARS) >0
 
}


LogRatio_LM_model = function (option_list)
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
  cat(str(SELECTED_VARS))
  cat("\n")
  
  
  
  
  #### Read transposed expression ----
  
  if(length(SELECTED_VARS) >0)
  {
    
    if (file.exists(opt$BP_transcript_EXP)){
      
      BP_transcript_EXP<-readRDS(file=opt$BP_transcript_EXP)
      
      # as.data.frame(fread(file=opt$BP_transcript_EXP, sep=",", header=T) , stringsAsFactors=F)
      
      
      cat("BP_transcript_EXP\n")
      cat(str(BP_transcript_EXP))
      cat("\n")
      cat(sprintf(as.character(names(summary(2^BP_transcript_EXP$value)))))
      cat("\n")
      cat(sprintf(as.character(summary(2^BP_transcript_EXP$value))))
      cat("\n")
      
      
      
       path6<-paste(out,SELECTED_VARS,'/', sep='')
      
      # cat("path6\n")
      # cat(sprintf(as.character(path6)))
      # cat("\n")
      
      
      if (file.exists(path6)){
        
        
        
        
      } else {
        dir.create(file.path(path6))
        
      }
      
      setwd(path6)
      
      filtered_transcripts_file<-paste("DTU_PASS_Transcripts_",SELECTED_VARS,'.rds', sep='')
      
      if (file.exists(filtered_transcripts_file)){
        
        TRANSCRIPTS_table_FILTERED_FINAL<-readRDS(file=filtered_transcripts_file)
        
      
        cat("TRANSCRIPTS_table_FILTERED_FINAL\n")
        cat(str(TRANSCRIPTS_table_FILTERED_FINAL))
        cat("\n")
       
        
        ENST_array<-unique(TRANSCRIPTS_table_FILTERED_FINAL$transcript_id)
        
        cat("ENST_array\n")
        cat(str(ENST_array))
        cat("\n")
        
        ENSG_array<-unique(TRANSCRIPTS_table_FILTERED_FINAL$ensembl_gene_id)
        
        cat("ENSG_array_0\n")
        cat(str(ENSG_array))
        cat("\n")
        
        
        
        #### MASTER LOOP ----
        
        
      
        
        
        list_RESULT<-list()
        
        for(i in 1:length(ENSG_array))
        {
          Condition_DEBUG <- 0
          
          ENSG_array_sel<-ENSG_array[i]
          
          cat("---------------->\t")
          cat(sprintf(as.character(i)))
          cat("\t")
          cat(sprintf(as.character(ENSG_array_sel)))
          cat("\n")
          
          TRANSCRIPTS_table_sel<-TRANSCRIPTS_table_FILTERED_FINAL[which(TRANSCRIPTS_table_FILTERED_FINAL$ensembl_gene_id == ENSG_array_sel),]
          
          
          if(Condition_DEBUG == 1)
          {
            cat("TRANSCRIPTS_table_sel\n")
            cat(str(TRANSCRIPTS_table_sel))
            cat("\n")
          }
          
          BP_transcript_EXP_sel<-BP_transcript_EXP[which(BP_transcript_EXP$transcript_id%in%TRANSCRIPTS_table_sel$transcript_id),]
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_transcript_EXP_sel_0\n")
            cat(str(BP_transcript_EXP_sel))
            cat("\n")
          }
          
          FLAG_genotypes_diversity<-length(unique(levels(droplevels(BP_transcript_EXP_sel$Genotype))))
          
          cat("FLAG_genotypes_diversity\n")
          cat(str(FLAG_genotypes_diversity))
          cat("\n")
          
          if(FLAG_genotypes_diversity >1)
          {
           
            
           
              
            BP_transcript_EXP_sel.m<-unique(BP_transcript_EXP_sel[,c(which(colnames(BP_transcript_EXP_sel) == "ensembl_gene_id"),
                                                                which(colnames(BP_transcript_EXP_sel) == "transcript_id"),
                                                                which(colnames(BP_transcript_EXP_sel) == "BP_ID"),
                                                                which(colnames(BP_transcript_EXP_sel) == "HGNC"),
                                                                which(colnames(BP_transcript_EXP_sel) == "value"),
                                                                which(colnames(BP_transcript_EXP_sel) == "Cell_Type"),
                                                                which(colnames(BP_transcript_EXP_sel) == "Genotype"))])
            HGNC_sel<-unique(BP_transcript_EXP_sel.m$HGNC)
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_genes_DEF_sel.m_0\n")
              cat(str(BP_transcript_EXP_sel.m))
              cat("\n")
            }
            
            
            
            BP_transcript_EXP_sel.m$FPKM<-2^(BP_transcript_EXP_sel.m$value)
            
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_genes_DEF_sel.m_1\n")
              cat(str(BP_transcript_EXP_sel.m))
              cat("\n")
              
              cat("ensembl_gene_id\n")
              cat(str(unique(BP_transcript_EXP_sel.m$ensembl_gene_id)))
              cat("\n")
              # ########################################
              # quit(status = 1)
            }
            
            
            
            cat("transcript_id\n")
            cat(str(unique(BP_transcript_EXP_sel.m$transcript_id)))
            cat("\n") 
            
            BP_transcript_EXP_sel.m_HET<-droplevels(BP_transcript_EXP_sel.m[which(BP_transcript_EXP_sel.m$Genotype != "HOM"),])
            
            
            if(Condition_DEBUG == 1)
            {
              
              cat("BP_genes_DEF_sel.m_HET_0\n")
              cat(str(BP_transcript_EXP_sel.m_HET))
              cat("\n")
              cat(sprintf(as.character(names(summary(BP_transcript_EXP_sel.m_HET$value)))))
              cat("\n")
              cat(sprintf(as.character(summary(BP_transcript_EXP_sel.m_HET$value))))
              cat("\n")
              # quit(status = 1)
              
            }
            
            #### LM per CellType----
            
            Cell_Type_array<-unique(as.character(BP_transcript_EXP_sel.m_HET$Cell_Type))
            
            if(Condition_DEBUG == 1)
            {
              cat("Cell_Type_array\n")
              cat(str(Cell_Type_array))
              cat("\n")
              
              # quit(status = 1)
            }
            
            list_CT<-list()
            
            for(iteration_Cell_Type in 1:length(Cell_Type_array))
            {
              Cell_Type_array_sel<-Cell_Type_array[iteration_Cell_Type]
              
              if(Condition_DEBUG == 1)
              {
                cat("------------------>\t")
                cat(sprintf(as.character(Cell_Type_array_sel)))
                cat("\n")
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel<-droplevels(BP_transcript_EXP_sel.m_HET[which(BP_transcript_EXP_sel.m_HET$Cell_Type == Cell_Type_array_sel),])
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel))
                cat("\n")
                
                # quit(status = 1)
              }
              
              n_transcripts_per_gene<-length(unique(BP_transcript_EXP_sel.m_HET_CT_sel$transcript_id))
              
              if(Condition_DEBUG == 1)
              {
                cat("n_transcripts_per_gene\n")
                cat(str(n_transcripts_per_gene))
                cat("\n")
              }
              
              if(n_transcripts_per_gene > 1)
              {
              
              
              #### calculate log model ----
              
              #### Find lowest expression per transcript that is different from 0
              
              
              BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO<-BP_transcript_EXP_sel.m_HET_CT_sel[which(BP_transcript_EXP_sel.m_HET_CT_sel$FPKM >0),]
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO_0\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO))
                cat("\n")
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel_ZERO<-droplevels(BP_transcript_EXP_sel.m_HET_CT_sel[which(BP_transcript_EXP_sel.m_HET_CT_sel$FPKM <= 0),])
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel_ZERO_0\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO))
                cat("\n")
                cat(sprintf(as.character(names(summary(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$FPKM)))))
                cat("\n")
                cat(sprintf(as.character(summary(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$FPKM))))
                cat("\n")
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO.dt<-data.table(BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO,
                                                                     key=c("transcript_id"))
              
              
              
              Zero_imputation_values<-as.data.frame(BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO.dt[,.(min_FPKM=min(FPKM)),by=key(BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO.dt)], stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                cat("Zero_imputation_values_0\n")
                cat(str(Zero_imputation_values))
                cat("\n")
                
                
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel_ZERO<-merge(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO,
                                                          Zero_imputation_values,
                                                          by="transcript_id",
                                                          all.x=T)
              
              BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$FPKM<-0.65*BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$min_FPKM
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel_ZERO_1\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO))
                cat("\n")
                cat(sprintf(as.character(names(summary(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$FPKM)))))
                cat("\n")
                cat(sprintf(as.character(summary(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO$FPKM))))
                cat("\n")
                
                # ##################################################
                # quit(status = 1)
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel_ZERO<-BP_transcript_EXP_sel.m_HET_CT_sel_ZERO[,-which(colnames(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO)== "min_FPKM")]
              
              BP_transcript_EXP_sel.m_HET_CT_sel<-rbind(BP_transcript_EXP_sel.m_HET_CT_sel_ZERO,
                                                     BP_transcript_EXP_sel.m_HET_CT_sel_NO_ZERO)
              
              
              
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel_RBIND\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel))
                cat("\n")
                cat(sprintf(as.character(names(summary(BP_transcript_EXP_sel.m_HET_CT_sel$FPKM)))))
                cat("\n")
                cat(sprintf(as.character(summary(BP_transcript_EXP_sel.m_HET_CT_sel$FPKM))))
                cat("\n")
              }
              
              
              ### Calculate transcript_reference
              
              BP_transcript_EXP_sel.m_HET_CT_sel.dt<-data.table(BP_transcript_EXP_sel.m_HET_CT_sel,
                                                             key=c("transcript_id"))
              
              Reference_value<-as.data.frame(BP_transcript_EXP_sel.m_HET_CT_sel.dt[,.(mean_FPKM=mean(FPKM)),by=key(BP_transcript_EXP_sel.m_HET_CT_sel.dt)], stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                cat("Reference_value_0\n")
                cat(str(Reference_value))
                cat("\n")
                
                
              }
              
              max_value<-max(Reference_value$mean_FPKM)
              
              Reference_value_MAX<-Reference_value[which(Reference_value$mean_FPKM == max_value),]
              
              Reference_transcript<-Reference_value_MAX$transcript_id[1]
              
              if(Condition_DEBUG == 1)
              {
                cat("max_value\n")
                cat(sprintf(as.character(max_value)))
                cat("\n")
                
                cat("Reference_value_MAX_0\n")
                cat(str(Reference_value_MAX))
                cat("\n")
                
                cat("Reference_transcript\n")
                cat(sprintf(as.character(Reference_transcript)))
                cat("\n")
                
                
              }
              
              
              ### calculate summatory FPKM per sample
              
              BP_transcript_EXP_sel.m_HET_CT_sel.dt<-data.table(BP_transcript_EXP_sel.m_HET_CT_sel,
                                                             key=c("BP_ID"))
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel.dt_0\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel.dt))
                cat("\n")
              }
              
              
              Summary_table_GENE_EXP<-as.data.frame(BP_transcript_EXP_sel.m_HET_CT_sel.dt[, .(sum_GENE_EXP=sum(FPKM)),
                                                                                       by=key(BP_transcript_EXP_sel.m_HET_CT_sel.dt)],stringsAsFactors=F)
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Summary_table_GENE_EXP_0\n")
                cat(str(Summary_table_GENE_EXP))
                cat("\n")
                
                #  quit(status = 1)
              }
              
              
              BP_transcript_EXP_sel.m_HET_CT_sel<-merge(BP_transcript_EXP_sel.m_HET_CT_sel,
                                                     Summary_table_GENE_EXP,
                                                     by="BP_ID")
              
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_transcript_EXP_sel.m_HET_CT_sel_1\n")
                cat(str(BP_transcript_EXP_sel.m_HET_CT_sel))
                cat("\n")
              }
              
              BP_transcript_EXP_sel.m_HET_CT_sel.dt<-data.table(BP_transcript_EXP_sel.m_HET_CT_sel,
                                                             key=c("BP_ID","Genotype","transcript_id"))
              
              Ratio_df<-as.data.frame(BP_transcript_EXP_sel.m_HET_CT_sel.dt[, .(FPKM=FPKM,
                                                                             sum_GENE_EXP=sum_GENE_EXP,
                                                                             Ratio=FPKM/sum_GENE_EXP),by=key(BP_transcript_EXP_sel.m_HET_CT_sel.dt)],stringsAsFactors=F)
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_0\n")
                cat(str(Ratio_df))
                cat("\n")
                
                cat("distrib_ratios\n")
                cat(sprintf(as.character(names(summary(Ratio_df$Ratio)))))
                cat("\n")
                cat(sprintf(as.character(summary(Ratio_df$Ratio))))
                cat("\n")
                
              }
              
              #### scale ratios to reference transcripts
              
              if(ENSG_array_sel == "ENSG00000029534")
              {
                
                if(dim(Ratio_df[which(Ratio_df$transcript_id == "ENST00000347528"),])[1] >0)
                {
                  Reference_transcript<-"ENST00000347528"
                }
                
              }
              
              if(ENSG_array_sel == "ENSG00000132394")
              {
                if(dim(Ratio_df[which(Ratio_df$transcript_id == "ENST00000483569"),])[1] >0)
                {
                  Reference_transcript<-"ENST00000483569"
                }
                
              }
              
              if(ENSG_array_sel == "ENSG00000257923")
              {
                if(dim(Ratio_df[which(Ratio_df$transcript_id == "ENST00000437600"),])[1] >0)
                {
                  Reference_transcript<-"ENST00000437600"
                }
              }
              
              if(Condition_DEBUG == 1)
              {
                cat("Reference_transcript_AFTER\n")
                cat(sprintf(as.character(Reference_transcript)))
                cat("\n")
                
                
              }
              
              Ratio_df_reference<-droplevels(Ratio_df[which(Ratio_df$transcript_id == Reference_transcript),c(which(colnames(Ratio_df) == "BP_ID"),
                                                                                                              which(colnames(Ratio_df) == "transcript_id"),
                                                                                                              which(colnames(Ratio_df) == "Ratio"))])
              
              colnames(Ratio_df_reference)[which(colnames(Ratio_df_reference) == "Ratio")]<-"Reference_ratio_value"
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_reference_0\n")
                cat(str(Ratio_df_reference))
                cat("\n")
                
                cat("distrib_ratios\n")
                cat(sprintf(as.character(names(summary(Ratio_df_reference$Reference_ratio_value)))))
                cat("\n")
                cat(sprintf(as.character(summary(Ratio_df_reference$Reference_ratio_value))))
                cat("\n")
                
              }
              
              Ratio_df_reference<-unique(Ratio_df_reference[,-which(colnames(Ratio_df_reference) == "transcript_id")])
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_reference_1\n")
                cat(str(Ratio_df_reference))
                cat("\n")
              }
              
              Ratio_df<-merge(Ratio_df,
                              Ratio_df_reference,
                              by="BP_ID",
                              all.x=T)
              
              Ratio_df$scaled_ratio<-Ratio_df$Ratio/Ratio_df$Reference_ratio_value
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_1\n")
                cat(str(Ratio_df))
                cat("\n")
                
                cat("distrib_ratios\n")
                cat(sprintf(as.character(names(summary(Ratio_df$scaled_ratio)))))
                cat("\n")
                cat(sprintf(as.character(summary(Ratio_df$scaled_ratio))))
                cat("\n")
                
              }
              
              check<-Ratio_df[which(Ratio_df$scaled_ratio > 1),]
              
              if(dim(check)[1] >0)
              {
                check2<-Ratio_df[which(Ratio_df$BP_ID%in%check$BP_ID),]
                check2<-check2[order(check2$BP_ID),]
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("check2\n")
                  cat(str(check2))
                  cat("\n")
                }
                
              }
              
             
              
              
              #### calculate the log
              
              Ratio_df$logscaled_ratio<-log(Ratio_df$scaled_ratio)
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("Ratio_df_DEF\n")
                cat(str(Ratio_df))
                cat("\n")
                
                
                
              }
              
              ##### Fit the LM per transcript ----
              
              ENST_array<-unique(as.character(Ratio_df$transcript_id))
              
              
              if(Condition_DEBUG == 1)
              {
                
                cat("ENST_array_1\n")
                cat(str(ENST_array))
                cat("\n")
                
              }
              
              
              list_transcript<-list()
              
              Condition_DEBUG <- 0
              
              
              for(k in 1:length(ENST_array))
              {
                ENST_array_sel<-ENST_array[k]
                
                if(Condition_DEBUG == 1)
                {
                  cat("--->\t")
                  cat(sprintf(as.character(ENST_array_sel)))
                  cat("\n")
                  
                }
                
                Ratio_df_ENST_sel<-Ratio_df[which(Ratio_df$transcript_id%in%ENST_array_sel),] 
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("Ratio_df_ENST_sel_0\n")
                  cat(str(Ratio_df_ENST_sel))
                  cat("\n")
                  cat(str(unique(Ratio_df_ENST_sel$transcript_id)))
                  cat("\n")
                  
                  
                  # quit(status = 1)
                  
                }
                
                
                
                ACCEPTED_genotypes<-c("HOM_REF","HET")
                
                Ratio_df_ENST_sel_HET<-Ratio_df_ENST_sel[which(Ratio_df_ENST_sel$Genotype%in%ACCEPTED_genotypes),]
                
                Ratio_df_ENST_sel_HET<-droplevels(Ratio_df_ENST_sel_HET)
                
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("Ratio_df_ENST_sel_HET_0\n")
                  cat(str(Ratio_df_ENST_sel_HET))
                  cat("\n")
                  
                  
                  # quit(status = 1)
                  
                }
                
                
                n_summary<-as.numeric(summary(Ratio_df_ENST_sel_HET$Genotype))
                
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("n_summary_0\n")
                  cat(str(n_summary))
                  cat("\n")
                  
                  
                  # quit(status = 1)
                  
                }
                
                n_summary_names<-names(summary(Ratio_df_ENST_sel_HET$Genotype))
                
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("n_summary_names_0\n")
                  cat(str(n_summary_names))
                  cat("\n")
                  
                  
                  # quit(status = 1)
                  
                }
                
                RUN_OUT_OF_NAMES<-NULL
                
                for(iteration_n_summary_names in 1:length(n_summary_names))
                {
                  n_summary_names_sel<-n_summary_names[iteration_n_summary_names]
                  n_summary_sel<-n_summary[iteration_n_summary_names]
                  
                  n_summary_string_sel<-paste(n_summary_names_sel,n_summary_sel,sep="__")
                  
                  RUN_OUT_OF_NAMES[iteration_n_summary_names]<-n_summary_string_sel
                  
                }#iteration_n_summary_names
                
                
                
                RUN_OUT_OF_NAMES_DEF<-paste(RUN_OUT_OF_NAMES, collapse=";")
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("RUN_OUT_OF_NAMES_DEF_0\n")
                  cat(str(RUN_OUT_OF_NAMES_DEF))
                  cat("\n")
                  
                  
                  # quit(status = 1)
                  
                }
                
                
                ### FULL linear model of the logRatio (Genotype+covs) ----
                
                unselected_columns<-c("BP_ID","transcript_id","FPKM","sum_GENE_EXP","Ratio","Reference_ratio_value","scaled_ratio")
                
                Selected_columns<-c("logRatio","Genotype")
                
              
                
                Ratio_df_ENST_sel_HET_LM<-Ratio_df_ENST_sel_HET[,-which(colnames(Ratio_df_ENST_sel_HET)%in%unselected_columns)]
                
                row.names(Ratio_df_ENST_sel_HET_LM)<-Ratio_df_ENST_sel_HET$BP_ID
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("Ratio_df_ENST_sel_HET_LM_0\n")
                  cat(str(Ratio_df_ENST_sel_HET_LM))
                  cat("\n")
                  
                  
                  
                  
                }
                
                
                model<-lm(logscaled_ratio ~ ., data=Ratio_df_ENST_sel_HET_LM)
                
                
                
                A<-summary(model)
                
                if(Condition_DEBUG == 1)
                {
                  cat("A\n")
                  cat(str(A))
                  cat("\n")
                  
                  
                }
                
                # quit(status = 1)
                
                results<-A$coefficients
                
                if(Condition_DEBUG == 1)
                {
                  cat("results_0\n")
                  cat(str(results))
                  cat("\n")
                  cat(sprintf(as.character(colnames(results))))
                  cat("\n")
                  cat(sprintf(as.character(row.names(results))))
                  cat("\n")
                  cat(sprintf(as.character(results[,4])))
                  cat("\n")
                  
                  results_G<-results[which(row.names(results) == "Genotype.L"),]
                  
                  cat("results_G\n")
                  cat(str(results_G))
                  cat("\n")
                }
                
                pvalue_Genotypes_1<-as.numeric(results[which(row.names(results) == "Genotype.L"),4])
                coefficient_Genotypes_1<-as.numeric(results[which(row.names(results) == "Genotype.L"),1])
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("pvalue_Genotypes_1\n")
                  cat(str(pvalue_Genotypes_1))
                  cat("\n")
                  cat("coefficient_Genotypes_1\n")
                  cat(str(coefficient_Genotypes_1))
                  cat("\n")
                }
                
                A.df<-as.data.frame(cbind(ENST_array_sel,pvalue_Genotypes_1, coefficient_Genotypes_1, RUN_OUT_OF_NAMES_DEF))
                
                colnames(A.df)<-c("transcript_id","pvalue_Genotypes","coefficient_Genotypes","n_breakdown_string")
                
                A.df$pvalue_Genotypes<-as.numeric(A.df$pvalue_Genotypes)
                A.df$coefficient_Genotypes<-as.numeric(A.df$coefficient_Genotypes)
                
                list_transcript[[k]]<-A.df
                
                
                if(Condition_DEBUG == 1)
                {
                  
                  cat("A.df\n")
                  cat(str(A.df))
                  cat("\n")
                  
                  
                  # quit(status=1)
                }
                
                # ###############################
                # quit(status = 1)
                
               
              }# k ENST_array
              
              Condition_DEBUG <- 1
              
              
              if(length(list_transcript) >0)
              {
                
                Results_per_CT = as.data.frame(data.table::rbindlist(list_transcript, fill=T), stringsAsFactors=F)
                
                Results_per_CT$Cell_Type<-Cell_Type_array_sel
                
                list_CT[[iteration_Cell_Type]]<-Results_per_CT
                
                if(Condition_DEBUG == 1)
                {
                  cat("Results_per_CT_0\n")
                  cat(str(Results_per_CT))
                  cat("\n")
                  
                  # ###################################################################
                  # quit(status = 1)
                  
                  #quit(status = 1)
                }
              }# length(list_transcript) >0
              
              }#n_transcripts_per_gene > 1
            } #iteration_Cell_Type in 1:length(Cell_Type_array)
          }#FLAG_genotypes_diversity >1
          
          if(length(list_CT) >0)
          {
            
            Results_per_gene = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
            
            Results_per_gene$ensembl_gene_id<-ENSG_array_sel
            
            list_RESULT[[i]]<-Results_per_gene
            
            if(Condition_DEBUG == 1)
            {
              cat("Results_per_gene_0\n")
              cat(str(Results_per_gene))
              cat("\n")
              
              # ###################################################################
              # quit(status = 1)
              
              #quit(status = 1)
            }
          }# length(list_CT) >0
          
        }#i in 1:length(ENSG_array)
        
        Condition_DEBUG <- 1
        
        
        if(length(list_RESULT) >0)
        {
          
          Results_uncorrected = as.data.frame(data.table::rbindlist(list_RESULT, fill=T), stringsAsFactors=F)
          
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_uncorrected_0\n")
            cat(str(Results_uncorrected))
            cat("\n")
            cat(str(unique(Results_uncorrected$ensembl_gene_id)))
            cat("\n")
            cat(str(unique(Results_uncorrected$transcript_id)))
            cat("\n")
          }
          
          Results_uncorrected<-merge(Results_uncorrected,
                                     TRANSCRIPTS_table_FILTERED_FINAL,
                                     by=c("ensembl_gene_id","transcript_id"))
          
          if(Condition_DEBUG == 1)
          {
            cat("Results_uncorrected_1\n")
            cat(str(Results_uncorrected))
            cat("\n")
            cat(str(unique(Results_uncorrected$ensembl_gene_id)))
            cat("\n")
            cat(str(unique(Results_uncorrected$transcript_id)))
            cat("\n")
          }
          
          
          path6<-paste(out,SELECTED_VARS,'/', sep='')
          
          # cat("path6\n")
          # cat(sprintf(as.character(path6)))
          # cat("\n")
          
          
          if (file.exists(path6)){
            
            
            
            
          } else {
            dir.create(file.path(path6))
            
          }
          
          setwd(path6)
          
          saveRDS(Results_uncorrected,file=paste("DTU_LogLM_HET_RESULTS_NOMINAL_",SELECTED_VARS,'.rds', sep=''))
          
        }# length(list_RESULT) >0
        
      }# file.exists(filtered_transcripts_file)
    }#file.exists(opt$BP_transcript_EXP)
  }#length(SELECTED_VARS) >0
  
 
  
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
    make_option(c("--Transcripts_table"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BP_transcript_EXP"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
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
  
  Filter_irrelevant_transcripts(opt)
  LogRatio_LM_model(opt)

}


###########################################################################

system.time( main() )
