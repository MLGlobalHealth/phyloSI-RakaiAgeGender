#!/bin/sh
#PBS -lselect=1:ncpus=1:mem=10GB
#PBS -lwalltime=05:00:00
#PBS -j oe 

# The key driver of this analysis is the STEP parameter
# which should be passed through the qsub command.
# qsub -v STEP="net" runall_TSI_seroconv2.sh
# If unset default to "sim"
if [ -z "$STEP" ]
then
    echo "Intended use:\n"
    echo 'qsub -v STEP="xxx" runall_TSI_pairs.sh'
    exit 1
fi

${RES:=1} 
${REDO:=0}
echo "running '${STEP:=sim}' analysis"

# This includes all code necessary to run PHSC pipeline to produce TSI estimates
DEEPDATA="/rds/general/project/ratmann_pangea_deepsequencedata/live"
DEEPANALYSES="/rds/general/project/ratmann_deepseq_analyses/live"
HOME="/rds/general/user/ab1820/home"
XIAOYUE="/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI/"

# User specific paths
software_path="$HOME/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software"
phyloscanner_path="$HOME/git/phyloscanner"
hivtsipath="$HOME/git/HIV-phyloTSI"

# analysis specific paths & args
out_dir_base="$DEEPANALYSES/PANGEA2_RCCS_MRC_UVRI_TSI"
out_dir_rel="$out_dir_base/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T"

controller="$software_path/$PBS_JOBNAME" #current script location
inputsamples="$out_dir_base/220331_RCCSUVRI_phscinput_samples_with_bf.rds"
CLUSIZE='50'
DATE='2022-08-22'

echo "Check that DATE, CLUSIZE, out_dir_rel and inputsamples are correctly specified"

cwd=$(pwd)
echo $cwd
module load anaconda3/personal
source activate phylo_alignments

case $STEP in

    # In this analysis we avoid the first step of computing similarities, as there exist already
    net)
    echo "---- initialise analysis ----"
    Rscript $software_path/TSI_initialise.R \
        --controller $controller \
        --include_least_recent_only TRUE \
        --out_dir_base $out_dir_base
    ;;

    ali)
    echo "---- compute alignments ----"
    if [ "$REDO" = "0"]; then
        Rscript $software_path/make_deep_sequence_alignments.R \
            --out_dir_base $out_dir_base \
            --pkg_dir $software_path \
            --prog_dir $phyloscanner_path \
            --windows_start 550 \
            --windows_end 9500 \
            --sliding_width 25 \
            --n_control 0 \
            --cluster_size $CLUSIZE \
            --reference ConsensusGenomes.fasta \
            --mafft " --globalpair --maxiterate 1000 " \
            --rm_vloops FALSE \
            --controller $controller \
            --walltime_idx $RES \
            --tsi_analysis FALSE
    else
        Rscript $software_path/make_deep_sequence_alignments.R \
            --out_dir_base $out_dir_base \
            --pkg_dir $software_path \
            --prog_dir $phyloscanner_path \
            --windows_start 550 \
            --windows_end 9500 \
            --sliding_width 25 \
            --n_control 0 \
            --cluster_size $CLUSIZE \
            --reference ConsensusGenomes.fasta \
            --mafft " --globalpair --maxiterate 1000 " \
            --rm_vloops FALSE \
            --controller $controller \
            --walltime_idx $RES \
            --date $DATE \
            --tsi_analysis FALSE
    fi
    ;;
        
    btr)
    echo "----- build trees ----"
    Rscript $software_path/make_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --controller $controller \
        --walltime_idx $RES
    ;;

    ctr)
    echo "----- check trees ----"
    Rscript $software_path/check_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --iqtree_method "GTR+F+R6" \
        --env_name "phylostan" \
        --date $DATE \
        --walltime_idx $RES
    ;;


    # atm modified from here...
    atr)
    conda activate phylostan
    Rscript $software_path/analyse_trees.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --prog_dir $phyloscanner_path \
        --normRefFileName "$DEEPDATA/normalisation_ByPosition.csv" \
        --outgroupName "REF_B.FR.83.HXB2_LAI_IIIB_BRU_K03455" \
        --ratioBlacklistThreshold 0.01 \
        --distanceThreshold "0.02 0.05"   \
        --maxReadsPerHost 100 \
        --minReadsPerHost 30  \
        --multinomial \
        --noProgressBars  \
        --postHocCountBlacklisting  \
        --relaxedAncestry \
        --zeroLengthAdjustment \
        --date $DATE \
        --controller $controller \
        --env_name "phylostan" \
        --verbose TRUE
    ;;
    
    # ... to here

    tsi)
    echo "----- Run HIV-TSI -----"
    Rscript $software_path/TSI_run_predictions.R \
        --out_dir_base $out_dir_base \
        --pkg_dir $software_path \
        --relationship_dir $out_dir_rel \
        --TSI_dir $hivtsipath \
        --date $DATE \
        --controller $controller \
        --env_name 'hivphylotsi'
    ;;

    dti)
    echo "----- get dates of infection -----"
    Rscript $software_path/TSI_estimate_dates.R \
        --out_dir_base $out_dir_base \
        --relationship_dir $out_dir_rel \
        --pkg_dir $software_path \
        --date $DATE \
        --input_samples $inputsamples \
        --controller $controller 
    ;;

    *)
    echo "no R script run. STEP does not match any task.\n" 
    ;;
esac
