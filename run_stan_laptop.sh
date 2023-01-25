#!/bin/bash
STAN_MODEL="gp_220108"
JOBNAME="cutoff_2014"
INDIR="$PWD"

# ========== Configure output directory ============
OUTDIR="/Users/shozendan/Imperial/phyloflows-output"
# ==================================================

# Create output directories if they don't exist
mkdir -p "$OUTDIR"
mkdir -p "$OUTDIR/${STAN_MODEL}-${JOBNAME}"
mkdir -p "$OUTDIR/${STAN_MODEL}-${JOBNAME}/figures"

# Create main script
cat > "$OUTDIR/bash_${STAN_MODEL}-${JOBNAME}.sh" <<EOF
#!/bin/bash
Rscript "$INDIR/scripts/run_stan.R" --indir "$INDIR" --outdir "$OUTDIR/${STAN_MODEL}-${JOBNAME}" --stan_model "$STAN_MODEL" --jobname "$JOBNAME"
EOF

# Create post-processing script
cat > "$OUTDIR/bash_${STAN_MODEL}-${JOBNAME}-postprocessing.sh" <<EOF
#!/bin/bash
mkdir -p "$OUTDIR/${STAN_MODEL}-${JOBNAME}/figures"
mkdir -p "$OUTDIR/${STAN_MODEL}-${JOBNAME}/tables"
Rscript "$INDIR/scripts/postprocessing_assess_mixing.R" --indir "$INDIR" --outdir "$OUTDIR/${STAN_MODEL}-${JOBNAME}" --stan_model "$STAN_MODEL" --jobname "$JOBNAME"
Rscript "$INDIR/scripts/postprocessing_figures.R" --indir "$INDIR" --outdir "$OUTDIR/${STAN_MODEL}-${JOBNAME}" --stan_model "$STAN_MODEL" --jobname "$JOBNAME"
EOF

# Execute the generated scripts in order
bash "$OUTDIR/bash_${STAN_MODEL}-${JOBNAME}.sh" && bash "$OUTDIR/bash_${STAN_MODEL}-${JOBNAME}-postprocessing.sh"

