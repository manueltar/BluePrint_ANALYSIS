

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
library("ggridges", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/")
suppressMessages(library("splitstackshape", lib.loc="/nfs/team151/software/manuel_R_libs_4_1//"))


suppressMessages(library("reshape2", lib.loc="/nfs/team151/software/manuel_R_libs_4_1/"))


opt = NULL

options(warn = 1)

fetch_data = function(option_list)
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
  
  ### Read GENES_table----
  
 
  GENES_table<-as.data.frame(fread(file=opt$GENES_table, sep="\t", header=F) , stringsAsFactors=F)
  
  colnames(GENES_table)<-c("chr","start","end","ensembl_gene_id","HGNC")
  cat("GENES_table\n")
  cat(str(GENES_table))
  cat("\n")
  
  
  GENES_table_subset<-GENES_table[,c(which(colnames(GENES_table) == "ensembl_gene_id"),
                                      which(colnames(GENES_table) == "HGNC"))]
  
 
  cat("GENES_table_subset\n")
  cat(str(GENES_table_subset))
  cat("\n")
  
  ### Read My_list_RNA_ids----
  
  
  Intersected_RNAids<-as.data.frame(readRDS(file=opt$Intersected_RNAids) , stringsAsFactors=F)
  
  cat("Intersected_RNAids\n")
  cat(str(Intersected_RNAids))
  cat("\n")
  cat(str(unique(Intersected_RNAids$VAR)))
  cat("\n")
  
  Intersected_RNAids_sel<-Intersected_RNAids[which(Intersected_RNAids$VAR %in% SELECTED_VARS),]
  
  
  cat("Intersected_RNAids_sel\n")
  cat(str(Intersected_RNAids_sel))
  cat("\n")
  cat(str(unique(Intersected_RNAids_sel$VAR)))
  cat("\n")
  
  if(dim(Intersected_RNAids_sel)[1] >0)
  {
    rm(Intersected_RNAids)
    
    HETS_sel<-unique(as.data.frame(cSplit(Intersected_RNAids_sel, splitCols = "string_Donors_RNA_HET",
                                          sep = "|", direction = "long", drop = TRUE),stringsAsFactors=F))
    
    
    cat("HETS_sel:0\n")
    cat(str(HETS_sel))
    cat("\n")
    
    
    
    HETS_sel<-HETS_sel$string_Donors_RNA_HET
    
    
    HETS_sel<-HETS_sel[!is.na(HETS_sel)]
    
    
    cat("HETS_sel1\n")
    cat(str(HETS_sel))
    cat("\n")
    
    
    HOMS_sel<-unique(as.data.frame(cSplit(Intersected_RNAids_sel, splitCols = "string_Donors_RNA_HOM",
                                          sep = "|", direction = "long", drop = TRUE),stringsAsFactors=F))
    
    
    cat("HOMS_sel:0\n")
    cat(str(HOMS_sel))
    cat("\n")
    
    
    HOMS_sel<-HOMS_sel$string_Donors_RNA_HOM
    
    
    HOMS_sel<-HOMS_sel[!is.na(HOMS_sel)]
    
    cat("HOMS_sel:1\n")
    cat(str(HOMS_sel))
    cat("\n")
    
    ############ FETCH FILE -----
    
    Condition_DEBUG <- 0
    
    rsid_sel<-unique(as.character(Intersected_RNAids_sel$rs))
    ref_sel<-unique(as.character(Intersected_RNAids_sel$ref))
    alt_sel<-unique(as.character(Intersected_RNAids_sel$alt))
    
    
    if(Condition_DEBUG == 1)
    {
      cat(sprintf(as.character(rsid_sel)))
      cat("\t")
      
      
      cat(sprintf(as.character(ref_sel)))
      cat("\t")
      
      
      cat(sprintf(as.character(alt_sel)))
      cat("\n")
    }
    
    kousik_path<-"/lustre/scratch123/hgi/mdt2/projects/hematopoiesis/Blueprint/Analysis/kk8/Manuel_variants/gene_expressio_results/"
    
    file_list <- list.files(path=kousik_path, include.dirs = FALSE)
    
    if(Condition_DEBUG == 1)
    {
      cat("file_list\n")
      cat(str(file_list))
      cat("\n")
    }
    
    
    indexes_gene_exp <- grep("gene_exp\\.txt$",file_list)
    
    if(Condition_DEBUG == 1)
    {
      cat("indexes_gene_exp\n")
      cat(str(indexes_gene_exp))
      cat("\n")
    }
    
    file_list_gene_exp<-file_list[indexes_gene_exp]
    
    if(Condition_DEBUG == 1)
    {
      cat("file_list_gene_exp\n")
      cat(str(file_list_gene_exp))
      cat("\n")
    }
    
    indexes_rsid <- grep(rsid_sel,file_list_gene_exp)
    
    if(Condition_DEBUG == 1)
    {
      cat("indexes_rsid\n")
      cat(str(indexes_rsid))
      cat("\n")
    }
    
    file_list_gene_exp_rsid<-file_list_gene_exp[indexes_rsid]
    
    if(Condition_DEBUG == 1)
    {
      cat("file_list_gene_exp_rsid\n")
      cat(str(file_list_gene_exp_rsid))
      cat("\n")
    }
    
    File_array<-unique(unlist(file_list_gene_exp_rsid))
    
    if(Condition_DEBUG == 1)
    {
      cat("File_array\n")
      cat(str(File_array))
      cat("\n")
    }
    
    ENSG_array<-gsub("gene_exp\\.txt$","",File_array)
    
    if(Condition_DEBUG == 1)
    {
      cat("ENSG_array_0\n")
      cat(str(ENSG_array))
      cat("\n")
    }
    
    ENSG_array<-gsub(paste(rsid_sel,"_",sep=''),"",ENSG_array)
    ENSG_array<-gsub("_$","",ENSG_array)
    
    
    if(Condition_DEBUG == 1)
    {
      cat("ENSG_array_1\n")
      cat(str(ENSG_array))
      cat("\n")
    }
    
    if(length(ENSG_array) >0)
    {
      list_genes<-list()
      
      for(i in 1:length(ENSG_array))
      {
        
        ENSG_array_sel<-ENSG_array[i]
        
        cat("----------------------->\t")
        cat(sprintf(as.character(ENSG_array_sel)))
        cat("\n")
        
        
        setwd(kousik_path)
        
        
        filename=paste(rsid_sel,ENSG_array_sel,"gene_exp.txt",sep='_')
        
        
        if (file.exists(filename)) {
          
          BP_EXP <- as.data.frame(fread(filename, sep="\t",header=T), stringsAsFactors = F)
          
          if(Condition_DEBUG == 1)
          {
            cat("BP_EXP_0\n")
            str(BP_EXP)
            cat("\n")
          }
          
          alt_sel_BP<-unique(BP_EXP$Gen)
          
          if(Condition_DEBUG == 1)
          {
            cat("--->\t")
            cat(sprintf(as.character(alt_sel_BP)))
            cat("\n")
          }
          
          if(alt_sel_BP == alt_sel)
          {
            HOM_REF_string<-paste(ref_sel,'/',ref_sel,sep='')
            
            if(Condition_DEBUG == 1)
            {
              cat("--->\t")
              cat(sprintf(as.character(HOM_REF_string)))
              cat("\n")
            }
            
            BP_EXP$Genotype[which(BP_EXP$Genotype == HOM_REF_string)]<-"HOM_REF"
            
            HET_string_1<-paste(alt_sel,'/',ref_sel,sep='')
            HET_string_2<-paste(ref_sel,'/',alt_sel,sep='')
            
            HET_string<-c(HET_string_1,HET_string_2)
            
            if(Condition_DEBUG == 1)
            {
              cat("--->\t")
              cat(sprintf(as.character(HET_string)))
              cat("\n")
            }
            
            BP_EXP$Genotype[which(BP_EXP$Genotype%in%HET_string)]<-"HET"
            
            HOM_string_1<-paste(alt_sel,'/',alt_sel,sep='')
            
            HOM_string<-c(HOM_string_1)
            
            if(Condition_DEBUG == 1)
            {
              cat("--->\t")
              cat(sprintf(as.character(HOM_string)))
              cat("\n")
            }
            
            BP_EXP$Genotype[which(BP_EXP$Genotype%in%HOM_string)]<-"HOM"
            
            BP_EXP$Genotype<-factor(BP_EXP$Genotype,
                                    levels=c("HOM_REF","HET","HOM"),
                                    ordered=T)
            if(Condition_DEBUG == 1)
            {
              cat("BP_EXP_1\n")
              str(BP_EXP)
              cat("\n")
              cat(sprintf(as.character(names(summary(BP_EXP$Genotype)))))
              cat("\n")
              cat(sprintf(as.character(summary(BP_EXP$Genotype))))
              cat("\n")
            }
            
            BP_EXP_subset<-BP_EXP[,c(which(colnames(BP_EXP) == "Value"),
                                     which(colnames(BP_EXP) == "Genotype"),
                                     which(colnames(BP_EXP) == "CellType"),
                                     which(colnames(BP_EXP) == "Doner"))]
            
            colnames(BP_EXP_subset)[which(colnames(BP_EXP_subset) == "Value")]<-"value"
            colnames(BP_EXP_subset)[which(colnames(BP_EXP_subset) == "CellType")]<-"Cell_Type"
            colnames(BP_EXP_subset)[which(colnames(BP_EXP_subset) == "Doner")]<-"BP_ID"
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_EXP_subset_0\n")
              str(BP_EXP_subset)
              cat("\n")
            }
            
            BP_EXP_subset<-BP_EXP_subset[,-which(colnames(BP_EXP_subset) == "Doner.1")]
            
            if(Condition_DEBUG == 1)
            {
              cat("BP_EXP_subset_0.5\n")
              str(BP_EXP_subset)
              cat("\n")
            }
            
            GENES_table_subset_sel<-GENES_table_subset[which(GENES_table_subset$ensembl_gene_id == ENSG_array_sel),]
            
            if(Condition_DEBUG == 1)
            {
              cat("GENES_table_subset_sel:\n")
              cat(str(GENES_table_subset_sel))
              cat("\n")
            }
            
            
            
            
            
            if(dim(GENES_table_subset_sel)[1] >0)
            {
              BP_EXP_subset$ensembl_gene_id<-GENES_table_subset_sel$ensembl_gene_id
              BP_EXP_subset$HGNC<-GENES_table_subset_sel$hgnc
              BP_EXP_subset$rsid<-rsid_sel
              BP_EXP_subset$VAR<-SELECTED_VARS
              
              if(Condition_DEBUG == 1)
              {
                cat("BP_EXP_subset_1\n")
                str(BP_EXP_subset)
                cat("\n")
              }
              
              list_genes[[i]]<-BP_EXP_subset
              
            }
            
            
            
            # ############################################################################################################################################################################
            # quit(status=1)
          }
        }#file.exists(filename)
        
        
      }#i in 1:length(ENSG_array)
      
      
      
      
      
      
      
      #### Recompose at the end ----
      
      if(length(list_genes) >0)
      {
        BP_genes_DEF = as.data.frame(data.table::rbindlist(list_genes, fill=T), stringsAsFactors=F)
        
        
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_genes_DEF\n")
          cat(str(BP_genes_DEF))
          cat("\n")
        }
        
        
        
        check<-BP_genes_DEF[is.na(BP_genes_DEF$value),]
        
        if(dim(check)[1]>0)
        {
          cat("check\n")
          cat(str(check))
          cat("\n")
          
          ensembl_gene_id_check<-unique(check$ensembl_gene_id)
          rsid_check<-unique(check$rsid)
          
          cat("WARNING_empty_values_for_genes!!!!\t")
          cat(sprintf(as.character(ensembl_gene_id_check)))
          cat("\t")
          cat(sprintf(as.character(rsid_check)))
          cat("\n")
          
          # ############################################################################################################################################################################
          # quit(status=1)
          
        }
        
        BP_genes_DEF<-BP_genes_DEF[!is.na(BP_genes_DEF$value),]  
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_genes_DEF_2\n")
          cat(str(BP_genes_DEF))
          cat("\n")
        }
        
        BP_genes_DEF<-merge(GENES_table_subset,
                            BP_genes_DEF,
                            by="ensembl_gene_id",
                            all.y=T)
        
        Condition_DEBUG <- 1
        
        
        if(Condition_DEBUG == 1)
        {
          cat("BP_genes_DEF_3\n")
          cat(str(BP_genes_DEF))
          cat("\n")
          cat(str(unique(BP_genes_DEF$HGNC)))
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
        
        saveRDS(file="BP_for_LM.rds", BP_genes_DEF)
        
        
      }# length(list_genes) >0
      
    }# length(ENSG_array) >0
    
  }# dim(Intersected_RNAids_sel)[1] >0

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
    make_option(c("--Intersected_RNAids"), type="character", default=NULL,
                metavar="type",
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--GENES_table"), type="character", default=NULL,
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
  
  fetch_data(opt)

  
}


###########################################################################

system.time( main() )
