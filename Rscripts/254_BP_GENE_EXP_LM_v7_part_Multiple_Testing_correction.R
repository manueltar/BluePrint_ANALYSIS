

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


MT_correction_in_cis_gene = function (option_list)
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
  
  #### VEP_CSQ ----
  
  VEP_CSQ<-as.data.frame(fread(file=opt$VEP_CSQ,sep=",") , stringsAsFactors=F)
  
  cat("VEP_CSQ_0\n")
  cat(str(VEP_CSQ))
  cat("\n")
  
  VEP_CSQ_sel<-VEP_CSQ[which(VEP_CSQ$VAR%in%SELECTED_VARS),]
  
  cat("VEP_CSQ_sel_0\n")
  cat(str(VEP_CSQ_sel))
  cat("\n")
  
  rm(VEP_CSQ)
  
  
  #### CIS consequences -----
  
  CIS_LABELS <-c("LOF","MISS","SYN","UTR5","UTR3",
                 "INTRON","SPLICE","UPSTREAM")
  
  
  VEP_CSQ_sel_CIS<-VEP_CSQ_sel[which(VEP_CSQ_sel$VEP_DEF_LABELS%in%CIS_LABELS),]
  
  cat("VEP_CSQ_sel_CIS_0\n")
  cat(str(VEP_CSQ_sel_CIS))
  cat("\n")
  
  #### SELECTED GENES ----
  
  ENSG_array<-vector()
  
  
  
  if(dim(VEP_CSQ_sel_CIS)[1] >0)
  {
    
    ENSG_sel<-unique(VEP_CSQ_sel_CIS$ensembl_gene_id)
    
    
    ENSG_array<-c(ENSG_array,ENSG_sel)
    
  }
  
  cat("ENSG_array\n")
  cat(str(ENSG_array))
  cat("\n")
  
  
  if(length(ENSG_array) >0)
  {
    
    
    
    
    ##### path6 ---
    
    path6<-paste(out,SELECTED_VARS,'/', sep='')
    
    # cat("path6\n")
    # cat(sprintf(as.character(path6)))
    # cat("\n")
    
    
    if (file.exists(path6)){
      
      
      
      
    } else {
      dir.create(file.path(path6))
      
    }
    
    setwd(path6)
    
    filename="LM_HET_RESULTS_NOMINAL.rds"
    
    
    if (file.exists(filename)) {
      Results_Nominal = readRDS(filename)
      
      cat("Results_Nominal\n")
      cat(str(Results_Nominal))
      cat("\n")
      
      Results_Nominal_sel<-Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),]
      
      cat("Results_Nominal_sel\n")
      cat(str(Results_Nominal_sel))
      cat("\n")
      
      # quit(status = 1)
      
      ##### Correction per cell type -----
      
      
      CT_levels<-unique(as.character(Results_Nominal_sel$Cell_Type))
      
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
        
        Results_Nominal_sel_CT_sel<-Results_Nominal_sel[which(Results_Nominal_sel$Cell_Type == CT_levels_sel),]
        
        cat("Results_Nominal_sel_CT_sel\n")
        cat(str(Results_Nominal_sel_CT_sel))
        cat("\n")
        
        
        #### SELECT TESTS TO CORRECT ----
        
        
        indx.select<-c(which(colnames(Results_Nominal_sel_CT_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_sel_CT_sel)== "HGNC"),which(colnames(Results_Nominal_sel_CT_sel)== "transcript_id"),
                       which(colnames(Results_Nominal_sel_CT_sel)== "coefficient_Genotypes"),
                       which(colnames(Results_Nominal_sel_CT_sel)== "pvalue_Genotypes"),
                       which(colnames(Results_Nominal_sel_CT_sel)== "n_breakdown_string"))
        
        
        
        MT_set<-unique(Results_Nominal_sel_CT_sel[,indx.select])
        
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
        
        path6<-paste(out,SELECTED_VARS,'/', sep='')
        
        # cat("path6\n")
        # cat(sprintf(as.character(path6)))
        # cat("\n")
        
        
        if (file.exists(path6)){
          
          
          
          
        } else {
          dir.create(file.path(path6))
          
        }
        
        setwd(path6)
        if(dim(CT_df)[1] >0)
        {
          saveRDS(file="DE_CIS_gene.rds", CT_df)  
        }
        
      }#length(list_CT) > 0
      
    }#file.exists(filename)
  }# length(ENSG_array) >0
}


MT_correction_in_block_Plus_PCHi_C = function (option_list)
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
  
  
  PCHiC_info<-as.data.frame(fread(file=opt$PCHiC_info, sep =",", header =T), stringsAsFactors=F)
  
  cat("PCHiC_info:\n")
  cat(str(PCHiC_info))
  cat("\n")
  
  #### GENES_PER_BLOCKS ----
  
  GENES_PER_BLOCKS = as.data.frame(fread(file=opt$GENES_PER_BLOCKS, sep="\t", stringsAsFactors = F, header = T))
  
  # cat("GENES_PER_BLOCKS\n")
  # cat(str(GENES_PER_BLOCKS))
  # cat("\n")
  
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
  
  
  #### SELECTED GENES ----
  
  ENSG_array<-vector()
  
  
  AS_Gene_source_sel<-AS_Gene_source[which(AS_Gene_source$VAR%in%SELECTED_VARS),]
  
  cat("AS_Gene_source_sel\n")
  cat(str(AS_Gene_source_sel))
  cat("\n")
  
  if(dim(AS_Gene_source_sel)[1] >0)
  {
    
    ENSG_sel<-unique(AS_Gene_source_sel$ensembl_gene_id)
    
    
    ENSG_array<-c(ENSG_array,ENSG_sel)
    
  }
  
  
  PCHiC_info_sel<-PCHiC_info[which(PCHiC_info$VAR%in%SELECTED_VARS),]
  
  if(dim(PCHiC_info_sel)[1] >0)
  {
    cat("PCHiC_info_sel\n")
    cat(str(PCHiC_info_sel))
    cat("\n")
    
    ENSG_sel<-unique(PCHiC_info_sel$ensembl_gene_id)
    
    
    ENSG_array<-c(ENSG_array,ENSG_sel)
  }
  
  if(SELECTED_VARS == "chr9_135920196_C_T")
  {
    ENSG_array<-unique(c(ENSG_array,"ENSG00000165702"))
  }
  
  
  cat("ENSG_array\n")
  cat(str(ENSG_array))
  cat("\n")
  
  
  
 
  ##### path6 ---
  
  path6<-paste(out,SELECTED_VARS,'/', sep='')
  
  # cat("path6\n")
  # cat(sprintf(as.character(path6)))
  # cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  setwd(path6)
  
  filename="LM_HET_RESULTS_NOMINAL.rds"
  
  
  if (file.exists(filename)) {
    Results_Nominal = readRDS(filename)
    
    cat("Results_Nominal\n")
    cat(str(Results_Nominal))
    cat("\n")
    
    Results_Nominal_sel<-Results_Nominal[which(Results_Nominal$ensembl_gene_id%in%ENSG_array),]
    
    cat("Results_Nominal_sel\n")
    cat(str(Results_Nominal_sel))
    cat("\n")
    
    
    ##### Correction per cell type -----
    
    
    CT_levels<-unique(as.character(Results_Nominal_sel$Cell_Type))
    
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
      
      Results_Nominal_sel_CT_sel<-Results_Nominal_sel[which(Results_Nominal_sel$Cell_Type == CT_levels_sel),]
      
      cat("Results_Nominal_sel_CT_sel\n")
      cat(str(Results_Nominal_sel_CT_sel))
      cat("\n")
      
      
      #### SELECT TESTS TO CORRECT ----
      
      
      indx.select<-c(which(colnames(Results_Nominal_sel_CT_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_sel_CT_sel)== "HGNC"),which(colnames(Results_Nominal_sel_CT_sel)== "transcript_id"),
                     which(colnames(Results_Nominal_sel_CT_sel)== "coefficient_Genotypes"),
                     which(colnames(Results_Nominal_sel_CT_sel)== "pvalue_Genotypes"),
                     which(colnames(Results_Nominal_sel_CT_sel)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_Nominal_sel_CT_sel[,indx.select])
      
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
      
      cat("CT_df_0_Block\n")
      cat(str(CT_df))
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
      if(dim(CT_df)[1] >0)
      {
        saveRDS(file="DE_Block.rds", CT_df)  
      }
      
    }#length(list_CT) > 0
    
  }#file.exists(filename)
  
  
  
}



MT_correction_genome_wide = function (option_list)
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

    ##### path6 ---
  
  path6<-paste(out,SELECTED_VARS,'/', sep='')
  
  # cat("path6\n")
  # cat(sprintf(as.character(path6)))
  # cat("\n")
  
  
  if (file.exists(path6)){
    
    
    
    
  } else {
    dir.create(file.path(path6))
    
  }
  
  setwd(path6)
   
  
  
  #### Read accepted transcripts for the variant ----
  
  filename="LM_HET_RESULTS_NOMINAL.rds"
  
  
  if (file.exists(filename)) {
    Results_Nominal = readRDS(filename)
    
    cat("Results_Nominal\n")
    cat(str(Results_Nominal))
    cat("\n")
    
    ##### Correction per cell type -----
    
    
    CT_levels<-unique(as.character(Results_Nominal$Cell_Type))
    
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
      
      Results_Nominal_CT_sel<-Results_Nominal[which(Results_Nominal$Cell_Type == CT_levels_sel),]
      
      cat("Results_Nominal_CT_sel\n")
      cat(str(Results_Nominal_CT_sel))
      cat("\n")
      
      
      #### SELECT TESTS TO CORRECT ----
      
      
      indx.select<-c(which(colnames(Results_Nominal_CT_sel)== "ensembl_gene_id"),which(colnames(Results_Nominal_CT_sel)== "HGNC"),which(colnames(Results_Nominal_CT_sel)== "transcript_id"),
                     which(colnames(Results_Nominal_CT_sel)== "coefficient_Genotypes"),
                     which(colnames(Results_Nominal_CT_sel)== "pvalue_Genotypes"),
                     which(colnames(Results_Nominal_CT_sel)== "n_breakdown_string"))
      
      
      
      MT_set<-unique(Results_Nominal_CT_sel[,indx.select])
      
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
      
      cat("CT_df_0_Genome_wide\n")
      cat(str(CT_df))
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
      if(dim(CT_df)[1] >0)
      {
        saveRDS(file="DE_genome_wide.rds", CT_df)  
      }
      
    }#length(list_CT) > 0
    
    
    
    
  }# file.exists(filename
  
 
  
  
 
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
  
  MT_correction_in_cis_gene(opt)
  MT_correction_in_block_Plus_PCHi_C(opt)
  MT_correction_genome_wide(opt)
  
}



###########################################################################

system.time( main() )
