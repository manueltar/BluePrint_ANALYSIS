#0;136;0c#!/bin/bash>

 
#### Rscript
 
Rscript=/software/R-4.1.0/bin/Rscript

SELECTED_VARS=$(echo 'chr10_9032555_C_G,chr1_101744646_C_G,chr1_158613314_G_A,chr1_158638260_G_A,chr4_1008212_C_T,chr14_74218569_C_T,chr12_110331967_G_A,chr2_169313518_C_T,chr12_112092462_G_A,chr12_22610292_G_C,chr12_53309769_A_G,chr9_135857646_C_T,chr7_50444152_G_T,chr14_23587046_A_C,chr14_51132622_T_C,chr6_7143859_C_A,chr14_74630170_A_G,chr1_202129205_G_A,chr8_130641322_C_T,chr15_65174494_A_G,chr15_65174438_C_A,chr16_24761046_G_A,chr16_333719_G_C,chr16_67250992_C_T,chr16_67690688_C_T,chr16_85595360_G_C,chr16_86016328_C_T,chr16_89094897_C_T,chr17_16949211_C_A,chr17_27197056_G_T,chr17_27778073_C_T,chr17_38764524_T_A,chr17_47856909_A_G,chr17_56339594_A_C,chr13_28604007_T_C,chr17_58602131_G_A,chr17_7106378_G_A,chr18_42041131_T_G,chr18_60880701_T_C,chr18_60920854_C_T,chr18_67856078_G_A,chr19_10676941_G_A,chr8_41589736_T_G,chr6_34947254_A_G,chr1_91606142_A_G,chr9_114663385_T_C,chr19_11210157_C_T,chr19_3203962_C_T,chr1_93482787_G_A,chr19_35776481_C_T,chr20_25409287_A_G,chr20_37544151_C_T,chr20_55990370_A_T,chr2_144084356_C_T,chr2_144162105_A_G,chr19_15653669_T_C,chr22_18252442_G_A,chr17_56603493_C_T,chr22_28761148_C_T,chr22_39362450_C_T,chr2_24091099_C_T,chr1_92981236_T_G,chr2_31476771_G_C,chr2_46293826_C_T,chr2_74920648_G_A,chr3_128317978_C_T,chr3_128322617_G_A,chr3_17098399_A_G,chr3_184091102_T_G,chr3_46354444_C_T,chr3_71355240_G_C,chr6_41924998_C_T,chr5_1041433_C_T,chr5_1093511_G_A,chr6_41952511_T_G,chr5_75563535_G_A,chr5_35476470_G_T,chr22_50949811_T_C,chr7_100309180_A_G,chr1_92925654_G_C,chr6_82476412_C_T,chr6_82501719_C_T,chr7_100083971_G_A,chr1_29217311_G_A,chr7_100314474_C_T,chr7_101499930_G_A,chr7_139875230_C_T,chr2_219020958_C_T,chr16_155132_C_T,chr8_130429059_G_A,chr7_99760955_G_C,chr1_198680015_G_A,chr8_90995426_C_T,chr12_111844956_C_T,chr15_64349614_G_A,chr9_135874752_G_A,chr9_135920196_C_T,chr9_20576825_C_A,chr9_4533390_C_G,chr9_5079248_T_C')
 



MASTER_ROUTE=$1
mem=$2
pc=$3
queue=$4

output="/nfs/users/nfs_m/mt19/Scripts/Wraper_scripts/302_CREATION_GENE_EXP_BP.sh"

touch $output
echo -n "" > $output

echo "#!/bin/bash"  >> $output

# # # # CAREFUL!!!!

#rm -rf $MASTER_ROUTE
#mkdir -p $MASTER_ROUTE
#SELECTED_VARS=$(echo 'chr17_38764524_T_A')
#SELECTED_VARS=$(echo 'chr3_17098399_A_G')

a=($(echo "$SELECTED_VARS" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        VAR_sel=${i}
        echo "$VAR_sel"

 	VAR_ROUTE=$(echo "$MASTER_ROUTE""/""$VAR_sel""/")

#        rm -rf $VAR_ROUTE
	#        mkdir -p $VAR_ROUTE

      Rscript_fetch_transcripts=/nfs/users/nfs_m/mt19/Scripts/R/360_BP_Transcript_EXP_fetch_data.R

      type=$(echo "fetch_transcripts""_""$VAR_sel")
      outfile_fetch_transcripts=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_fetch_transcripts
      echo -n "" > $outfile_fetch_transcripts
       name_fetch_transcripts=$(echo "$type""_job")

      Intersected_RNAids=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/INTERVAL/RNA_Seq/INTERVAL_CARRIERS_RNA_Seq_SEQUENCED_DEF.rds")
      Transcripts_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_Transcripts_table.txt")


      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_fetch_transcripts -M $DE_mem -J $name_fetch_transcripts -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_fetch_transcripts \\" >> $output
      echo "--Intersected_RNAids $Intersected_RNAids \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--Transcripts_table $Transcripts_table \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

      ### DTU_logRatio_model

      Rscript_DTU_logRatio_model=/nfs/users/nfs_m/mt19/Scripts/R/361_BP_DTU_part_LogRatio_LM.R

      type=$(echo "DTU_logRatio_model""_""$VAR_sel")
      outfile_DTU_logRatio_model=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_DTU_logRatio_model
      echo -n "" > $outfile_DTU_logRatio_model
       name_DTU_logRatio_model=$(echo "$type""_job")

      Transcripts_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_Transcripts_table.txt")
      BP_transcript_EXP=$(echo "$VAR_ROUTE""BP_Transcript_EXP_for_LM.rds")

      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_DTU_logRatio_model -M $DE_mem -w\"done($name_fetch_transcripts)\" -J $name_DTU_logRatio_model -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_DTU_logRatio_model -M $DE_mem -J $name_DTU_logRatio_model -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_DTU_logRatio_model \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--BP_transcript_EXP $BP_transcript_EXP \\" >> $output
      echo "--Transcripts_table $Transcripts_table \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output


       Rscript_fetch_data=/nfs/users/nfs_m/mt19/Scripts/R/254_BP_GENE_EXP_LM_v7_fetch_data.R

      type=$(echo "fetch_data""_""$VAR_sel")
      outfile_fetch_data=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_fetch_data
      echo -n "" > $outfile_fetch_data
       name_fetch_data=$(echo "$type""_job")

       Intersected_RNAids=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/INTERVAL/RNA_Seq/INTERVAL_CARRIERS_RNA_Seq_SEQUENCED_DEF.rds")
       GENES_table=$(echo "/nfs/users/nfs_m/mt19/RareVar_Dragana/Homo_sapiens.GRCh37.87_GENES_table.txt")


      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_fetch_data -M $DE_mem -J $name_fetch_data -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_fetch_data \\" >> $output
      echo "--Intersected_RNAids $Intersected_RNAids \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_table $GENES_table \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

       
      ### LM #########################################################################################################################################################################
      
      Rscript_LM=/nfs/users/nfs_m/mt19/Scripts/R/254_BP_GENE_EXP_LM_v7_LM.R


      type=$(echo "BP_DE_LM""_""$VAR_sel")
      outfile_LM=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_LM
      echo -n "" > $outfile_LM
      name_LM=$(echo "$type""_job")
      
           
      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


#      echo "bsub -G team151 -o $outfile_LM -M $DE_mem -w\"done($name_fetch_data)\" -J $name_LM -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "bsub -G team151 -o $outfile_LM -M $DE_mem -J $name_LM -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_LM \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output



      
       ### MT_correction #################################################################################################################################################################
      
      Rscript_MT_correction=/nfs/users/nfs_m/mt19/Scripts/R/254_BP_GENE_EXP_LM_v7_part_Multiple_Testing_correction.R

      type=$(echo "MT_correction_GENE_EXP_BP""_""$VAR_sel")
      outfile_MT_correction=$(echo "$VAR_ROUTE""outfile""_""$type"".out")
      touch $outfile_MT_correction
      echo -n "" > $outfile_MT_correction
        name_MT_correction=$(echo "$type""_job")


      ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
      GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
      PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
      VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")

      DE_mem=$(expr $mem \* 1)
      DE_pc=$(expr $pc \* 1)


      echo "$mem""->""$DE_mem"
      echo "$pc""->""$DE_pc"


      echo "bsub -G team151 -o $outfile_MT_correction -M $DE_mem -w\"done($name_LM)\" -J $name_MT_correction -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#      echo "bsub -G team151 -o $outfile_MT_correction -M $DE_mem -J $name_MT_correction -R\"select[model==Intel_Platinum]\"  -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
      echo "\"$Rscript $Rscript_MT_correction \\" >> $output
      echo "--SELECTED_VARS $VAR_sel \\" >> $output
      echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
      echo "--PCHiC_info $PCHiC_info \\" >> $output
      echo "--VEP_CSQ $VEP_CSQ \\" >> $output
      echo "--ALL_dB $ALL_dB \\" >> $output
      echo "--type $type --out $MASTER_ROUTE\"" >> $output

       
      
      if [ $VAR_sel == "chr10_9032555_C_G" ]; then 
      	  MT_correction_string=$(echo "done($name_MT_correction)")
      else
          MT_correction_string=$(echo "&& done($name_MT_correction)")

      fi

   
      echo "->>>$MT_correction_string"
      arr[${#arr[@]}]="$MT_correction_string"
      

done



done_string=$(echo "\""""${arr[@]}"""\"")
echo "$done_string"

 
#### Rscript_PUT_TOGETHER_RESULTS_ALL_BY_ALL ----
 
Rscript_PUT_TOGETHER_RESULTS_ALL_BY_ALL=/nfs/users/nfs_m/mt19/Scripts/R/254_BP_GENE_EXP_LM_v7_PUT_TOGETHER_RESULTS_ALL_BY_ALL_MT_correction.R

type=$(echo "PUT_TOGETHER_RESULTS_GENE_EXP_BP_ALL_BY_ALL")
outfile_PUT_TOGETHER_RESULTS_ALL_BY_ALL=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_PUT_TOGETHER_RESULTS_ALL_BY_ALL
echo -n "" > $outfile_PUT_TOGETHER_RESULTS_ALL_BY_ALL
name_PUT_TOGETHER_RESULTS_ALL_BY_ALL=$(echo "$type""_job")

DE_mem=$(expr $mem \* 2)
DE_pc=$(expr $pc \* 2)


echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"


#echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_ALL_BY_ALL -M $DE_mem -w$done_string -J $name_PUT_TOGETHER_RESULTS_ALL_BY_ALL -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "bsub -G team151 -o $outfile_PUT_TOGETHER_RESULTS_ALL_BY_ALL -M $DE_mem -J $name_PUT_TOGETHER_RESULTS_ALL_BY_ALL -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PUT_TOGETHER_RESULTS_ALL_BY_ALL \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output

#### Rscript_PUT_TOGETHER_DTU ----
 
Rscript_PUT_TOGETHER_DTU=/nfs/users/nfs_m/mt19/Scripts/R/362_BP_DTU_PUT_TOGETHER_RESULTS_MT_correction.R

type=$(echo "PUT_TOGETHER_DTU")
outfile_PUT_TOGETHER_DTU=$(echo "$MASTER_ROUTE""outfile""_""$type"".out")
touch $outfile_PUT_TOGETHER_DTU
echo -n "" > $outfile_PUT_TOGETHER_DTU
name_PUT_TOGETHER_DTU=$(echo "$type""_job")

DE_mem=$(expr $mem \* 2)
DE_pc=$(expr $pc \* 2)



ALL_dB=$(echo "/lustre/scratch123/hgi/mdt1/teams/soranzo/projects/Explore_Patricks_tables/ALL_db.tsv")
GENES_PER_BLOCKS=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/mt19/Allelic_Series_candidates_Allelic_Series_Generation.tsv")
PCHiC_info=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Variant_csv_tables/PCHIC_ChicagoScore_graphs.csv")
VEP_CSQ=$(echo "/lustre/scratch126/humgen/teams/soranzo/users/ALL_dB/Allelic_Series_csv_tables_AS_DEFINITIVE/VEP_CSQ.csv")


echo "$mem""->""$DE_mem"
echo "$pc""->""$DE_pc"

echo "-------------------->$PCHiC_info"

echo "bsub -G team151 -o $outfile_PUT_TOGETHER_DTU -M $DE_mem -w$done_string -J $name_PUT_TOGETHER_DTU -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
#echo "bsub -G team151 -o $outfile_PUT_TOGETHER_DTU -M $DE_mem -J $name_PUT_TOGETHER_DTU -R\"select[model==Intel_Platinum]\" -R\"select[mem>=$DE_mem] rusage[mem=$DE_mem] span[hosts=1]\" -n$DE_pc -q $queue -- \\" >> $output
echo "\"$Rscript $Rscript_PUT_TOGETHER_DTU \\" >> $output
echo "--SELECTED_VARS $SELECTED_VARS \\" >> $output
echo "--ALL_dB $ALL_dB \\" >> $output
echo "--GENES_PER_BLOCKS $GENES_PER_BLOCKS \\" >> $output
echo "--PCHiC_info $PCHiC_info \\" >> $output
echo "--VEP_CSQ $VEP_CSQ \\" >> $output
echo "--type $type --out $MASTER_ROUTE\"" >> $output


bash $output

  
