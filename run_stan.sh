#!/bin/sh

JOBID=$$
STAN_MODEL="gp_211207"
TAG="MRC_FALSE_OnlyHTX_TRUE_threshold_0.5"
CWD="/rds/general/user/mm3218/home/projects/2021/phyloflows"/$TAG
INDIR="/rds/general/user/mm3218/home/git/phyloflows"

mkdir $CWD

cat > $CWD/bash_$STAN_MODEL-$JOBID.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
  
JOB_TEMP=\${EPHEMERAL}/\${PBS_JOBID}
mkdir -p \$JOB_TEMP
cd \$JOB_TEMP  
PWD=\$(pwd)
CWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
TAG=$TAG
  
# main directory
mkdir \$PWD/\$STAN_MODEL-\$JOBID
  
# directories for fits, figures and tables
mkdir \$PWD/\$STAN_MODEL-\$JOBID/fits
mkdir \$PWD/\$STAN_MODEL-\$JOBID/data
mkdir \$PWD/\$STAN_MODEL-\$JOBID/figure
mkdir \$PWD/\$STAN_MODEL-\$JOBID/table
  
Rscript \$INDIR/scripts/run_stan.R -indir \$INDIR -outdir \$PWD -stan_model \$STAN_MODEL -JOBID \$JOBID -lab \$TAG
  
cp -R --no-preserve=mode,ownership \$PWD/* \$CWD
  
cd \$CWD
qsub bash_\$STAN_MODEL-\$JOBID-postprocessing.pbs

EOF

cat > $CWD/bash_$STAN_MODEL-$JOBID-postprocessing.pbs <<EOF
  
#!/bin/sh
#PBS -l walltime=40:59:00
#PBS -l select=1:ncpus=10:ompthreads=1:mem=240gb
#PBS -j oe
module load anaconda3/personal
  
CWD=$CWD
PWD=$CWD
INDIR=$INDIR/inst
STAN_MODEL=$STAN_MODEL
JOBID=$JOBID
TAG=$TAG
  
# main directory
mkdir \$PWD/\$STAN_MODEL-\$JOBID
  
# directories for fits, figures and tables
mkdir \$PWD/\$STAN_MODEL-\$JOBID/figure
mkdir \$PWD/\$STAN_MODEL-\$JOBID/table
  
# to write ....
EOF

  
cd $CWD
qsub bash_$STAN_MODEL-$JOBID.pbs