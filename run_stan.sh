#!/bin/sh

STAN_MODEL="gp_220108"
JOBNAME="cutoff_2014"
INDIR="/rds/general/user/mm3218/home/git/phyloflows"
OUTDIR="/rds/general/user/mm3218/home/projects/2021/phyloflows"

mkdir $OUTDIR

cat > $OUTDIR/bash_$STAN_MODEL-$JOBNAME.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=48:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
  
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
  
Rscript \$INDIR/scripts/process_data.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME
Rscript \$INDIR/scripts/run_stan.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME
  
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

INDIR=$INDIR
OUTDIR=$OUTDIR
STAN_MODEL=$STAN_MODEL
JOBNAME=$JOBNAME
  
# main directory
CWD=\$OUTDIR/\$STAN_MODEL-\$JOBNAME

# directories for figure and table
mkdir \$CWD/figures
mkdir \$CWD/tables

Rscript \$INDIR/scripts/postprocessing_assess_mixing.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 
Rscript \$INDIR/scripts/postprocessing_figures.R -indir \$INDIR -outdir \$CWD -stan_model \$STAN_MODEL -jobname \$JOBNAME 

EOF
  
cd $OUTDIR
qsub bash_$STAN_MODEL-$JOBNAME.pbs