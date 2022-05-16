# Define functions----
#
#	function to make PBS header
make.PBS.header <- function(hpc.walltime=47, hpc.select=1, hpc.nproc=1, hpc.mem= "6gb", hpc.load= "module load anaconda3/personal\nsource activate high_res_contacts", hpc.q="pqcovid19c", hpc.array=1 )
{
  pbshead <- "#!/bin/sh"
  tmp <- paste("#PBS -l walltime=", hpc.walltime, ":59:00", sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  tmp <- paste("#PBS -l select=", hpc.select, ":ncpus=", hpc.nproc,":ompthreads=", hpc.nproc,":mem=", hpc.mem, sep = "")
  pbshead <- paste(pbshead, tmp, sep = "\n")
  pbshead <- paste(pbshead, "#PBS -j oe", sep = "\n")
  if (hpc.array > 1)
  {
    pbshead	<- paste(pbshead, "\n#PBS -J 1-", hpc.array, sep = '')
  }
  if (!is.na(hpc.q))
  {
    pbshead <- paste(pbshead, paste("#PBS -q", hpc.q), sep = "\n")
  }
  pbshead	<- paste(pbshead, hpc.load, sep = "\n")
  pbshead
}

# yu's hpc input args age x gender ----
if (1)
{
  args <- list()
  # determine what to run
  args$run_analysis <- list(
    run_Poisson_BSGP_inland = 1,
    run_Poisson_BSGP_fishing = 1,
    run_Poisson_HSGP_inland = 1,
    run_Poisson_HSGP_fishing = 1,
    run_negbin_BSGP_inland = 1,
    run_negbin_BSGP_fishing = 1,
    run_negbin_HSGP_inland = 1,
    run_negbin_HSGP_fishing = 1
  )

  args$seed <- 18L
  args$on_hpc <- TRUE

  # define which datase
  args$data_set <- 'RCCS_reported_partnership_220505.csv'
  args$n_knots <- 30 # for B spline
  args$spline_degree <- 3
  args$hsgp_boundary_inflation <- 1.5 # for HSGP model
  args$hsgp_m <- 20 # for HSGP model

  # files
  args$pkg_dir <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership"
  args$out_dir_base <- "/rds/general/user/yc2819/home/github/phyloflows/sexual_partnership/results"
}

# Run Poisson BSGP Inland----
if (args$run_analysis$run_Poisson_BSGP_inland)
{
  out.dir <- paste0("BSGP_inland_knots-",args$n_knots, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S"))

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'inland',"\n")
  cmd <- paste0(cmd,"spline_degree=",args$spline_degree,"\n")
  cmd <- paste0(cmd,"n_knots=",args$n_knots,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_BSGP_Poisson.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --spline_degree $spline_degree',
                ' --n_knots $n_knots'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_Poisson.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}

# Run Poisson BSGP Fishing----
if (args$run_analysis$run_Poisson_BSGP_fishing)
{
  out.dir <- paste0("BSGP_fishing_knots-",args$n_knots, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S"))

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'fishing',"\n")
  cmd <- paste0(cmd,"spline_degree=",args$spline_degree,"\n")
  cmd <- paste0(cmd,"n_knots=",args$n_knots,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_BSGP_Poisson.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --spline_degree $spline_degree',
                ' --n_knots $n_knots'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_Poisson.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}


# Run Poisson HSGP Inland----
if (args$run_analysis$run_Poisson_HSGP_inland)
{
  out.dir <- paste0("HSGP_inland_c-",args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S")
  )

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'inland',"\n")
  cmd <- paste0(cmd,"hsgp_boundary_inflation=",args$hsgp_boundary_inflation,"\n")
  cmd <- paste0(cmd,"hsgp_m=",args$hsgp_m,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_HSGP_Poisson.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --hsgp_boundary_inflation $hsgp_boundary_inflation',
                ' --hsgp_m $hsgp_m'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_Poisson.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}
# Run Poisson HSGP Fishing----
if (args$run_analysis$run_Poisson_HSGP_fishing)
{
  out.dir <- paste0("HSGP_fishing_c-",args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S")
  )

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'fishing',"\n")
  cmd <- paste0(cmd,"hsgp_boundary_inflation=",args$hsgp_boundary_inflation,"\n")
  cmd <- paste0(cmd,"hsgp_m=",args$hsgp_m,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_HSGP_Poisson.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --hsgp_boundary_inflation $hsgp_boundary_inflation',
                ' --hsgp_m $hsgp_m'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_Poisson.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}


# Run NegBin BSGP Inland----
if (args$run_analysis$run_negbin_BSGP_inland)
{
  out.dir <- paste0("NegBin_BSGP_inland_knots-",args$n_knots, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S"))

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'inland',"\n")
  cmd <- paste0(cmd,"spline_degree=",args$spline_degree,"\n")
  cmd <- paste0(cmd,"n_knots=",args$n_knots,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_BSGP_NegBin.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --spline_degree $spline_degree',
                ' --n_knots $n_knots'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_NegBin.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}

# Run NegBin BSGP Fishing----
if (args$run_analysis$run_negbin_BSGP_fishing)
{
  out.dir <- paste0("NegBin_BSGP_fishing_knots-",args$n_knots, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S"))

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'fishing',"\n")
  cmd <- paste0(cmd,"spline_degree=",args$spline_degree,"\n")
  cmd <- paste0(cmd,"n_knots=",args$n_knots,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_BSGP_NegBin.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --spline_degree $spline_degree',
                ' --n_knots $n_knots'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_NegBin.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}


# Run NegBin HSGP Inland----
if (args$run_analysis$run_negbin_HSGP_inland)
{
  out.dir <- paste0("NegBin_HSGP_inland_c-",args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S")
  )

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'inland',"\n")
  cmd <- paste0(cmd,"hsgp_boundary_inflation=",args$hsgp_boundary_inflation,"\n")
  cmd <- paste0(cmd,"hsgp_m=",args$hsgp_m,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_HSGP_NegBin.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --hsgp_boundary_inflation $hsgp_boundary_inflation',
                ' --hsgp_m $hsgp_m'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_NegBin.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}
# Run NegBin HSGP Fishing----
if (args$run_analysis$run_negbin_HSGP_fishing)
{
  out.dir <- paste0("NegBin_HSGP_fishing_c-",args$hsgp_boundary_inflation*100,"_m-",args$hsgp_m, "_",
                    format(Sys.time(),"%y-%m-%d-%H-%M-%S")
  )

  out.dir <- file.path(args$out_dir_base, out.dir)
  if (!dir.exists(out.dir))
  {
    dir.create(out.dir)
  }
  cat("\noutput directory is ",out.dir)

  cmd <- ''
  cmd <- paste0(cmd,"CWD=$(pwd)\n")
  cmd <- paste0(cmd,"echo $CWD\n")
  tmpdir.prefix <- paste0('csim_',format(Sys.time(),"%y-%m-%d-%H-%M-%S"))
  tmpdir <- paste0("$CWD/",tmpdir.prefix)
  cmd <- paste0(cmd,"mkdir -p ",tmpdir,'\n')
  cmd <- paste0(cmd,"pkg_dir=",args$pkg_dir,"\n")
  cmd <- paste0(cmd,"seed=",args$seed,"\n")
  cmd <- paste0(cmd,"data_set=",args$data_set,"\n")
  cmd <- paste0(cmd,"community=",'fishing',"\n")
  cmd <- paste0(cmd,"hsgp_boundary_inflation=",args$hsgp_boundary_inflation,"\n")
  cmd <- paste0(cmd,"hsgp_m=",args$hsgp_m,"\n")
  cmd <- paste0(cmd,"out_dir_base=",tmpdir,"\n")
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_HSGP_NegBin.R'),
                ' --pkg_dir $pkg_dir',
                ' --out_dir_base $out_dir_base',
                ' --seed $seed',
                ' --data_set $data_set',
                ' --community $community',
                ' --hsgp_boundary_inflation $hsgp_boundary_inflation',
                ' --hsgp_m $hsgp_m'
  )
  cmd <- paste0(cmd, tmp, '\n')
  tmp <- paste0('Rscript ', file.path('$pkg_dir','scripts','make_partnerships_rate_postprocessing_NegBin.R'),
                ' --seed $seed',
                ' --pkg_dir $pkg_dir',
                ' --out_dir $out_dir_base'
  )
  cmd <- paste0(cmd, tmp, '\n')

  if (!args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R "', tmpdir,'"/* ',out.dir,'\n')
  }
  if (args$on_hpc)
  {
    cmd <- paste0(cmd, 'cp -R --no-preserve=mode,ownership "', tmpdir,'"/* ',out.dir,'\n')
  }
  cmd <- paste0(cmd, 'chmod -R g+rw ', out.dir,'\n')
  cmd <- paste0(cmd,"cd $CWD\n")

  if (!args$on_hpc)
  {
    cmd <- paste( cmds, collapse = '\n\n')
  }
  if (args$on_hpc)
  {
    pbshead <- make.PBS.header(	hpc.walltime = 10,
                                hpc.select = 1,
                                hpc.nproc = 4,
                                hpc.mem = "50gb",
                                hpc.q = NaN,
                                hpc.load = "module load anaconda3/personal\nsource activate high_res_contacts\nexport TBB_CXX_TYPE=gcc\nexport CXXFLAGS+=-fPIE"
    )
    cmd <- paste(pbshead,cmd ,sep = '\n')

  }

  jobfile <- gsub(':','',paste("csim",paste(strsplit(date(),split = ' ')[[1]],collapse = '_',sep = ''),'sh', sep = '.'))
  jobfile <- file.path(args$out_dir_base, jobfile)
  cat("\nWrite job script to file ", jobfile)
  cat(cmd, file = jobfile)

  if (args$on_hpc)
  {
    cmd <- paste("qsub", jobfile)
    cat(cmd)
    cat(system(cmd, intern = TRUE))
  }
}

