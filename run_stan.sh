#!/bin/sh

JOBID=$$
STAN_MODEL="gp_211117"
TAG="MRC_FALSE_OnlyHTX_TRUE_threshold_0.5_jobname_firstruncutoff"
DATADIR="/rds/general/user/mm3218/home/projects/2021/phyloflows"/$TAG
INDIR="/rds/general/user/mm3218/home/git/phyloflows"

mkdir $DATADIR

cat > $DATADIR/bash_$STAN_MODEL-$JOBID.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
  
JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP  
PWD=\$(pwd)
DATADIR=$DATADIR
INDIR=$INDIR
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
TAG=$TAG
  
# main directory
mkdir \$PWD/\$STAN_MODEL-\$JOBID
  
Rscript \$INDIR/scripts/run_stan.R -indir \$INDIR -datadir \$DATADIR -outdir \$PWD -stan_model \$STAN_MODEL -JOBID \$JOBID -lab \$TAG
  
cp -R --no-preserve=mode,ownership \$PWD/* \$DATADIR
  
cd \$DATADIR
qsub bash_$STAN_MODEL-$JOBID-postprocessing.pbs

EOF

cat > $DATADIR/bash_$STAN_MODEL-$JOBID-postprocessing.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:ompthreads=1:mem=480gb
#PBS -j oe
module load anaconda3/personal

JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP  
PWD=\$(pwd)
DATADIR=$DATADIR
INDIR=$INDIR
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
TAG=$TAG

# directories for figure and table
mkdir \$PWD/figures
mkdir \$PWD/tables

Rscript \$INDIR/scripts/postprocessing_results.R -indir \$INDIR -datadir \$DATADIR -outdir \$PWD -stan_model \$STAN_MODEL -JOBID \$JOBID -lab \$TAG

cp -R --no-preserve=mode,ownership \$PWD/* \$DATADIR/\$STAN_MODEL-\$JOBID

EOF
  
cd $DATADIR
qsub bash_$STAN_MODEL-$JOBID.pbs