#!/bin/bash

# use as:
# phyloSI-RakaiAgeGende-run_phyloflows-hpc.bash JOBNAME="jobname"
# Flags:
#    STAN_MODEL : stan model to be run
#    JOBNAME : name of the job 
#    INDIR : directory where github directory is located
#    OUTDIR : directory where output will be saved
#    ENVNAME : name of the conda environment to be used

for ARGUMENT in "$@"
do
   KEY=$(echo $ARGUMENT | cut -f1 -d=)

   KEY_LENGTH=${#KEY}
   VALUE="${ARGUMENT:$KEY_LENGTH+1}"

   echo "export $KEY=$VALUE"
   export "$KEY"="$VALUE"
done

# default options
STAN_MODEL="${STAN_MODEL:-gp_230602}"
JOBNAME="${JOBNAME:-central}"
INDIR="${INDIR:-/home/andrea/git/phyloflows}"
OUTDIR="${OUTDIR:-/home/andrea/git/phyloflows/output}"
ENVNAME="${ENVNAME:-phylowSI-RakaiAgeGender}"

echo "Selected options:"
echo "STAN_MODEL = $STAN_MODEL"
echo "JOBNAME = $JOBNAME"
echo "INDIR = $INDIR"
echo "OUTDIR = $OUTDIR"
echo "ENVNAME = $ENVNAME"

# if ENVNAME is not-empty, need to add a line to both scripts below
if [ ! -z "$ENVNAME" ]
then
    ENVNAME="source activate $ENVNAME"
fi
# 

mkdir $OUTDIR

cat > $OUTDIR/bash_$STAN_MODEL-$JOBNAME.pbs <<EOF
#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
$ENVNAME
  
JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP  
PWD=\$(pwd)

INDIR=$INDIR
OUTDIR=$OUTDIR
STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME
  
# main directory
CWD=\$PWD/\$STAN_MODEL-\$JOBNAME

mkdir \$CWD
mkdir \$CWD/figures
  
Rscript \$INDIR/src/transmission_flows/run_stan.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME
  
cp -R --no-preserve=mode,ownership \$PWD/* \$OUTDIR
  
cd \$OUTDIR
qsub bash_$STAN_MODEL-$JOBNAME-postprocessing.pbs

EOF

cat > $OUTDIR/bash_$STAN_MODEL-$JOBNAME-postprocessing.pbs <<EOF
#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=480gb
#PBS -j oe
module load anaconda3/personal
$ENVNAME

INDIR=$INDIR
OUTDIR=$OUTDIR
STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME
  
# main directory
CWD=\$OUTDIR/\$STAN_MODEL-\$JOBNAME

# directories for figure and table
mkdir \$CWD/figures
mkdir \$CWD/tables

Rscript \$INDIR/src/transmission_flows/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
Rscript \$INDIR/src/transmission_flows/postprocessing_figures.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
Rscript \$INDIR/src/transmission_flows/postprocessing_figure_time_trends_sources.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
Rscript \$INDIR/src/transmission_flows/postprocessing-figure_contribution_sexual_contact.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
Rscript \$INDIR/src/transmission_flows/postprocessing_figure_counterfactual.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 

EOF
  
cd $OUTDIR
qsub bash_$STAN_MODEL-$JOBNAME.pbs
