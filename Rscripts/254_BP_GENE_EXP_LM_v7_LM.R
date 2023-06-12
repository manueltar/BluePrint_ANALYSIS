

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


LM_model = function (option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
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
  
 
  
### Find if file exists -----
  
  path6<-paste(out,SELECTED_VARS,'/', sep='')
  
  # cat("path6\n")
  # cat(sprintf(as.character(path6)))
  # cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  setwd(path6)
  
  Condition_DEBUG <- 0

  if (file.exists("BP_for_LM.rds")){
    
    BP_genes_DEF<-readRDS(file="BP_for_LM.rds")
    
    
    if(Condition_DEBUG == 1)
    {
      cat("BP_genes_DEF_0\n")
      cat(str(BP_genes_DEF))
      cat("\n")
      
      # BP_genes_DEF_0
      # 'data.frame':   559 obs. of  8 variables:
      #   $ ensembl_gene_id: chr  "ENSG00000111252" "ENSG00000111252" "ENSG00000111252" "ENSG00000111252" ...
      # $ HGNC           : chr  "SH2B3" "SH2B3" "SH2B3" "SH2B3" ...
      # $ value          : num  13.1 12.6 9.7 13.1 12.7 ...
      # $ Genotype       : Ord.factor w/ 3 levels "HOM_REF"<"HET"<..: 1 1 1 1 1 1 1 1 1 1 ...
      # $ Cell_Type      : chr  "Monocyte" "Neutrophil" "Tcell" "Monocyte" ...
      # $ BP_ID          : chr  "S000GZ" "S000GZ" "S000GZ" "S000X1" ...
      # $ rsid           : chr  "rs540639423" "rs540639423" "rs540639423" "rs540639423" ...
      # $ VAR            : chr  "chr12_111844956_C_T" "chr12_111844956_C_T" "chr12_111844956_C_T" "chr12_111844956_C_T" ...
      # 
    }
    
    ENSG_array<-unique(BP_genes_DEF$ensembl_gene_id)
    
    # ENSG_array<-"ENSG00000126353"
    
    if(Condition_DEBUG == 1)
    {
      cat("ENSG_array_0\n")
      cat(str(ENSG_array))
      cat("\n")
    }
    
    FLAG_genotypes_diversity<-length(unique(levels(droplevels(BP_genes_DEF$Genotype))))

    cat("FLAG_genotypes_diversity\n")
    cat(str(FLAG_genotypes_diversity))
    cat("\n")

    if(FLAG_genotypes_diversity >1)
    {
      #### MASTER LOOP ----

      list_RESULT<-list()

      

      for(i in 1:length(ENSG_array))
      {
        ENSG_array_sel<-ENSG_array[i]

        cat("------------------------------>\t")
        cat(sprintf(as.character(i)))
        cat("\t")
        cat(sprintf(as.character(ENSG_array_sel)))
        cat("\n")
        
        # c(which(colnames(BP_genes_DEF) == "sample_id"),
        #   which(colnames(BP_genes_DEF) %in% ENSG_array_sel)

        BP_genes_DEF_sel.m<-unique(BP_genes_DEF[which(BP_genes_DEF$ensembl_gene_id == ENSG_array_sel),c(which(colnames(BP_genes_DEF) == "ensembl_gene_id"),
                                                                                               which(colnames(BP_genes_DEF) == "HGNC"),
                                                                                               which(colnames(BP_genes_DEF) == "value"),
                                                                                               which(colnames(BP_genes_DEF) == "Cell_Type"),
                                                                                               which(colnames(BP_genes_DEF) == "Genotype"))])
        HGNC_sel<-unique(BP_genes_DEF_sel.m$HGNC)
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_genes_DEF_sel.m_0\n")
          cat(str(BP_genes_DEF_sel.m))
          cat("\n")
        }



        BP_genes_DEF_sel.m$FPKM<-+2^(BP_genes_DEF_sel.m$value)


        if(Condition_DEBUG == 1)
        {
          cat("BP_genes_DEF_sel.m_1\n")
          cat(str(BP_genes_DEF_sel.m))
          cat("\n")

          
          cat("ensembl_gene_id\n")
          cat(str(BP_genes_DEF_sel.m$ensembl_gene_id))
          cat("\n")

          # ########################################
          # quit(status = 1)
        }

        
        BP_genes_DEF_sel.m_HET<-droplevels(BP_genes_DEF_sel.m[which(BP_genes_DEF_sel.m$Genotype != "HOM"),])


        if(Condition_DEBUG == 1)
        {

          cat("BP_genes_DEF_sel.m_HET_0\n")
          cat(str(BP_genes_DEF_sel.m_HET))
          cat("\n")


          # quit(status = 1)

        }
        
        #### LM per CellType----

        Cell_Type_array<-unique(as.character(BP_genes_DEF_sel.m_HET$Cell_Type))
        
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
          
          BP_genes_DEF_sel.m_HET_CT_sel<-droplevels(BP_genes_DEF_sel.m_HET[which(BP_genes_DEF_sel.m_HET$Cell_Type == Cell_Type_array_sel),])
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_genes_DEF_sel.m_HET_CT_sel\n")
            cat(str(BP_genes_DEF_sel.m_HET_CT_sel))
            cat("\n")
            
            # quit(status = 1)
          }
          
          n_summary<-as.numeric(summary(BP_genes_DEF_sel.m_HET_CT_sel$Genotype))
          
          
          if(Condition_DEBUG == 1)
          {
            
            cat("n_summary_0\n")
            cat(str(n_summary))
            cat("\n")
            
            
            # quit(status = 1)
            
          }
          
          n_summary_names<-names(summary(BP_genes_DEF_sel.m_HET_CT_sel$Genotype))
          
          
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
          
          SELECTED_TERMS<-c("Genotype")
          
          
          
          BP_genes_DEF_sel.m_HET_CT_sel<-BP_genes_DEF_sel.m_HET_CT_sel[,c(which(colnames(BP_genes_DEF_sel.m_HET_CT_sel) == "FPKM"),which(colnames(BP_genes_DEF_sel.m_HET_CT_sel)%in%SELECTED_TERMS))]
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_genes_DEF_sel.m_HET_CT_sel\n")
            cat(str(BP_genes_DEF_sel.m_HET_CT_sel))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          model<-lm(FPKM ~ ., data=BP_genes_DEF_sel.m_HET_CT_sel)
         
         
          A<-summary(model)
          
          if(Condition_DEBUG == 1)
          {
            cat("A\n")
            cat(str(A))
            cat("\n")
            
            
          }
          
          
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
          
          A.df<-as.data.frame(cbind(ENSG_array_sel,pvalue_Genotypes_1, coefficient_Genotypes_1,RUN_OUT_OF_NAMES_DEF))
          
          colnames(A.df)<-c("ensembl_gene_id","pvalue_Genotypes","coefficient_Genotypes","n_breakdown_string")
          
          A.df$pvalue_Genotypes<-as.numeric(A.df$pvalue_Genotypes)
          A.df$coefficient_Genotypes<-as.numeric(A.df$coefficient_Genotypes)
          A.df$Cell_Type<-Cell_Type_array_sel
            
          
          
          list_CT[[iteration_Cell_Type]]<-A.df
          
          if(Condition_DEBUG == 1)
          {
            
            cat("A.df\n")
            cat(str(A.df))
            cat("\n")
            
            
            # quit(status=1)
          }
        
         
          
        }# iteration_Cell_Type
        
        
        if(length(list_CT) > 0)
        {
          CT_df = as.data.frame(data.table::rbindlist(list_CT, fill=T), stringsAsFactors=F)
          CT_df$ensembl_gene_id<-ENSG_array_sel
          CT_df$HGNC<-HGNC_sel
          CT_df$pvalue_Genotypes<-as.numeric(CT_df$pvalue_Genotypes)
        
        
          
          
          if(Condition_DEBUG == 1)
          {
            cat("CT_df_0\n")
            cat(str(CT_df))
            cat("\n")
            
            # quit(status = 1)
          }
        }# length(list_CT) > 0

        list_RESULT[[i]]<-CT_df
      
      }#i ENSG_array


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
          cat(str(unique(Results_uncorrected$ensembl_gene_id)))
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

        saveRDS(file="LM_HET_RESULTS_NOMINAL.rds", Results_uncorrected)
      }# length(list_RESULT) >0

    }#FLAG_genotypes_diversity >1
    # 
    # # quit(status=1)
    
  } else {
    
    cat("NO_BP_file\n")
    cat("\n")
  }
  
  
  
  
  

 
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
  
  LM_model(opt)
  
  
  
  
}


###########################################################################

system.time( main() )
