# The phyloscanner run - find transmission networks 
cat("\n===== find_chains_from_phylogenetics.R =====\n")

# Preamble This script aims to find transmission networks and the most likely transmission chains using the phyloscanner outputs.

# Load the required packages
library(data.table) |> suppressPackageStartupMessages()
library(tidyverse) |> suppressPackageStartupMessages()
library(dplyr) |> suppressPackageStartupMessages()
library(glue) |> suppressPackageStartupMessages()
library(igraph) |> suppressPackageStartupMessages()
library(RBGL) |> suppressPackageStartupMessages()
library(phyloscannerR) |> suppressPackageStartupMessages()
library(here) |> suppressPackageStartupMessages()

#
# Define input arguments that can be changed by users
#
option_list <- list(
optparse::make_option( 
    "--confidential",
    type = "logical",
    default = FALSE,
    help = "Flag on whether to use the confidential data (if access is granted) [Defaults to TRUE]", 
    dest = 'confidential'
),
  optparse::make_option(
    "--phylo-pairs-dir",
    type = "character",
    # default = "/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/deep_sequence_phylogenies_primary/data_for_likely_transmission_pairs/phyloscanner-results",
    default = NA_character_,
    help = "Absolute file path to base directory where the phyloscanner outputs are stored [default]",
    dest = 'phylo.pairs.dir'
  ),
  optparse::make_option(
    "--outdir",
    type = "character",
    default = NULL,
    help = "Absolute file path to output directory to save the scripts outputs in", 
    dest = 'outdir'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# Paths
#

gitdir <- here()

source(file.path(gitdir, 'config.R'))

# if output directory is null, set it to proc
if( is.null(args$outdir) )
    args$outdir <- dir.zenodo.phyloproc

# Chose confidential or randomized files:
if(args$confidential)
{
    path.meta.data  <-  path.meta.confidential
    path.sequence.dates <- path.collection.dates.confidential 
}else{
    path.meta.data  <- path.meta.randomized
    path.sequence.dates <-  path.collection.dates.randomized
}

if(is.na(args$phylo.pairs.dir))
{
    args$phylo.pairs.dir <- dir.zenodo.pairs.phsc
}

dir.exists(args$phylo.pairs.dir) |> stopifnot()
file.exists(path.meta.data) |> stopifnot()

#
# Helpers
#

catn <- function(x) cat('\n----', x ,'----\n')

read.infiles <- function(DT)
{
    with(DT,{
        dca <- data.table()
        dwina <- data.table()
        for(pty in unique(PTY_RUN) )
        {
            cat(pty, '\n')
            load(F[PTY_RUN == pty])
            dc[, PTY_RUN := pty]
            dwin[, PTY_RUN := pty]
            rbind(dca, dc)
            rbind(dwina, dwin)
        }
        list(dca, dwina)
    }) -> out
    out
}
get.sample.collection.dates <- function(select_aid=NULL, get_first_visit=FALSE)
{
    ddates <- dseqdates

    if(!is.null(select_aid))
        ddates <- ddates[aid %in% select_aid]

    if(get_first_visit)
        ddates <- ddates[, .(date_collection=min(visit_dt)),by='aid']
    ddates
}

get.infection.range.from.testing <- function()
{
    # get maximum and minimum dates

    # for missing 1st pos, use date collection instead (also store contradicting first pos)
    tmp <- meta[, get.sample.collection.dates(aid, get_first_visit = TRUE)]
    meta <- merge(meta, tmp, by='aid', all.x=TRUE)

    # check whether anybody tested positive before 15yo
    meta[, date15th := date_birth + as.integer(365.25*15)]
    stopifnot(meta[date15th > date_first_positive, .N == 0])

    contradict_firstpos_datecoll <<- meta[date_collection < date_first_positive]
    drange <- meta[, .(AID=aid, 
                       MIN=pmax(date_last_negative,  date15th, na.rm=TRUE),
                       MAX=pmin(date_first_positive - 30, date_collection - 30, na.rm=TRUE))]
    drange[, lowb15lastneg := MAX - as.integer(365.25*15) ]
    drange[  lowb15lastneg > MIN, MIN:=lowb15lastneg ]
    drange[, lowb15lastneg := NULL]
    drange
}

.format.controls <- function(DT)
{
        DT[, `:=` (CNTRL1=FALSE, CNTRL2=FALSE)]
        DT[ host.1 %like% 'CNTRL-', `:=` (CNTRL1=TRUE, host.1=gsub('CNTRL-', '', host.1))]
        DT[ host.2 %like% 'CNTRL-', `:=` (CNTRL2=TRUE, host.2=gsub('CNTRL-', '', host.2))]
}

.reorder.host.labels <- function(DT)
{
        tmp <- subset(DT, host.1 > host.2)

        cols1 <- c('host.1','host.2','paths12','paths21','nodes1','nodes2','CNTRL1','CNTRL2')
        cols2 <- c('host.2','host.1','paths21','paths12','nodes2','nodes1','CNTRL2','CNTRL1')
        cols1 <- cols1[cols1 %in% names(DT)]
        cols2 <- cols2[cols2 %in% names(DT)]
        setnames(tmp, cols1, cols2)

        .gs <- function(x)
                gsub('xx','21',gsub('21','12',gsub('12','xx',x)))

        cols <- c('close.and.contiguous.and.directed.cat',
                  'close.and.adjacent.and.directed.cat',
                  'close.and.contiguous.and.ancestry.cat',
                  'close.and.adjacent.and.ancestry.cat', 
                  'type')
        cols <- cols[ cols %in% names(DT)]

        tmp[, (cols) := lapply(.SD, .gs) , .SDcols=cols]

        DT <- rbind(DT[!(host.1>host.2)], tmp)
        return(DT)
}

fill.in.missing.reverse.direction <- function(DCA)
{

    # the all.x and all.y options allow us to know which entries weren't previously in DCA
    cols <- c('PTY_RUN', 'host.1', 'host.2', 'categorisation', 'categorical.distance', 'CNTRL1', 'CNTRL2')
    merge(
        DCA[type == '12', ..cols],
        DCA[type == '21', ..cols],
        by=cols, all.x=TRUE, all.y=TRUE
    ) -> tmp
    tmp.12 <- merge(tmp, DCA[type == '12'], by=cols, all.x=TRUE)
    tmp.12[is.na(type), type := '12']
    tmp.21 <- merge(tmp, DCA[type == '21'], by=cols, all.x=TRUE)
    tmp.21[is.na(type), type := '21']
    tmp <- rbind(tmp.12, tmp.21)

    # check
    idx <- tmp[!is.na(n),uniqueN(n), by=cols][ V1 > 1, ..cols]
    stopifnot(nrow(idx)==0)

    # and update: for each "by" group, n & n_eff should be constant.
    # some scores will still be NA
    tmp[, n:=na.omit(n)[1], by=cols] 
    tmp[, n.eff:=na.omit(n.eff)[1], by=cols] 
    tmp[, score:= k.eff/n.eff]
    
    out <- rbind(DCA[! type %in% c('12', '21')], tmp)
    setkey(out, PTY_RUN, host.1, host.2, categorisation, categorical.distance, type)
    out
}

update.category.counts <- function(DEXCLUDE, DCA)
{
    # DEXCLUDE <- copy(dexclude); DCA <- copy(DCA1)
    idx.12 <- DEXCLUDE[EXCLUDE == '12', .(host.1, host.2)]
    idx.21 <- DEXCLUDE[EXCLUDE == '21', .(host.1, host.2)]
    setkey(DCA, host.1, host.2); setkey(dwina, host.1, host.2);
    setkey(idx.12, host.1, host.2);setkey(idx.21, host.1, host.2);
    tmp.12 <- merge(DCA, idx.12, by=c('host.1', 'host.2'))
    tmp.21 <- merge(DCA, idx.21, by=c('host.1', 'host.2'))

    update.category.counts.by.unsupported.direction <- function(ctg, DT, dir)
    {
        .update <- function(TMP)
        {
            with(TMP,{
                # Find values associated with 'impossible' label
                # cat(unique(host.1), unique(host.2),'\n')
                new.dir <- setdiff(c('12', '21'), dir)
                which.dir <- type == dir
                which.new.dir <- type == new.dir

                TMP$k[which.new.dir] <- k1 <- sum(TMP$k)
                TMP$k.eff[which.new.dir] <- k1.eff <- sum(TMP$k.eff)
                TMP$k[!which.new.dir] <- 0
                TMP$k.eff[!which.new.dir] <- 0
                TMP$score <- TMP$k.eff / TMP$n.eff

                if(! any(which.new.dir))
                {
                    tmp_row <- TMP[1, ]
                    tmp_row$type <- new.dir
                    tmp_row$k <- k1
                    tmp_row$k.eff <- k1.eff
                    tmp_row$score <- tmp_row$k.eff / tmp_row$n.eff
                    TMP <- rbind(tmp_row, TMP)
                    return(TMP)
                }

                return(TMP)
            })
        }
        cols <- setdiff(names(DT),'categorisation')
        out <- DT[categorisation == ctg, .update(.SD) , by=c('host.1', 'host.2', 'PTY_RUN', 'categorisation'), .SDcols=cols]
        return(out)
    }
    
    
    dca.update.12 <- lapply(categories.to.update, update.category.counts.by.unsupported.direction, tmp.12, '12')
    dca.update.21 <- lapply(categories.to.update, update.category.counts.by.unsupported.direction, tmp.21, '21')
    dca.update <- rbind(rbindlist(dca.update.12), rbindlist(dca.update.21))
    # for some reason some cols are duplicated
    dca.update <- subset(dca.update, select=!duplicated(names(dca.update)))
    setcolorder(dca.update,names(dca))
    dca.update[! k %in% c('12', '21'), sum(k)]

    # lapply(dca.update.12, function(DT) DT[type=='21', sum(k)])
    # lapply(dca.update.21, function(DT) DT[type=='21', sum(k)])

    setkey(dca.update, host.1, host.2, PTY_RUN, categorisation)
    setkey(DCA, host.1, host.2, PTY_RUN, categorisation)
    DCA <- rbind(DCA[!dca.update], dca.update)
}


########
# Main #
########

dseqdates <- readRDS(path.sequence.dates)

catn('Load phyloscanner outputs')
stopifnot(dir.exists(args$phylo.pairs.dir))
infiles	<- data.table(F=list.files(args$phylo.pairs.dir, pattern='*workspace.rda$', full.names=TRUE))
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
setkey(infiles, PTY_RUN)

# 'dwin' : pairwise relationships between the host in each tree
# 'dc' :  summarises pairwise relationships between hosts across ALL trees. 
# 'a' suffix stands for aggregated

# current solution only loads F once, instead of twice
out <- read.infiles(DT=infiles) 
dca   <- out[['dca']]
dwina <- out[['dwina']]

# Change the format
.format.controls(dca)
.format.controls(dwina)

# Sort, so host.1 < host.2 
dca   <- .reorder.host.labels(dca)
dwina <- .reorder.host.labels(dwina)
setkey(dwina, PTY_RUN, host.1, host.2)
setkey(dca, PTY_RUN, host.1, host.2)

# subset to phylogenies with strongest data on linkage 
idx <- dca[ categorisation %in% 'close.and.adjacent.cat',
           .(n.eff=unique(n.eff)),
           by=c('host.1', 'host.2', 'PTY_RUN')]

idx <- idx[ ,
           .(PTY_RUN=PTY_RUN[which.max(n.eff)]),
           by=c('host.1', 'host.2')]
setcolorder(idx, c('PTY_RUN','host.1', 'host.2'))
dca <- merge(idx, dca)
dwina <- merge(idx, dwina)

# 
if( file.exists(path.meta.data) )
{   # Now, whenever a pairwise transmission is not supported by the serohistory, we need to "switch" the counts supporting that direction.
    cat('Using serohistory + demographic data to reweight evidence of direction\n')
    
    # load meta data on infection times
    meta_env <- new.env()
    load(path.meta.data, envir=meta_env)
    meta <- subset(meta_env$meta_data,
                   select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
    meta <- unique(meta[!is.na(aid)])
    # idx <- meta[, uniqueN(comm) == 2, by=aid][V1==TRUE, aid]
    stopifnot(meta[, uniqueN(aid) == .N])

    # get pairs for which one direction of transmission is unsupported according to demographic and serohistory data.
    drange <- get.infection.range.from.testing()
    dexclude <- dca[, .(host.1, host.2)] |> unique()
    dexclude <- merge(dexclude, drange[, .(host.2=AID, MIN.H2=MIN, MAX.H2=MAX)] , by='host.2', all.x=TRUE)
    dexclude <- merge(dexclude, drange[, .(host.1=AID, MIN.H1=MIN, MAX.H1=MAX)] , by='host.1', all.x=TRUE)
    dexclude[ MAX.H1 <= MIN.H2, EXCLUDE := '21']
    dexclude[ MAX.H2 <= MIN.H1, EXCLUDE := '12']
    dexclude[, table(EXCLUDE, useNA = 'always' )/.N * 100 ]

    # check that for each '12' entry, there is a '21' entry in the count data.table
    dca <- fill.in.missing.reverse.direction(dca)

    # Now, modify dca and dwina accordingly, by removing the counts supporting that direction
    categories.to.update <- c(
        "close.and.adjacent.and.ancestry.cat",
        "close.and.adjacent.and.directed.cat",
        "close.and.contiguous.and.ancestry.cat",
        "close.and.contiguous.and.directed.cat")
    dca_sero_only <- dca[categorisation %in% categories.to.update, ]

    # label.to.update <- function(ctg)
    #   fcase(ctg %like% 'directed', NA_character_, ctg %like% 'ancestry.cat', 'complex.or.no.ancestry')
    dca_sero_only <- update.category.counts(dexclude, dca_sero_only)
    dca_sero_only[, categorisation:=paste0(categorisation, '.sero')]
    dca <- rbind(dca, dca_sero_only)

    setkey(dca, host.1, host.2)
}

# find pairs according to classification rule and thresholds.
# classification rule o: Oliver Ratmann's
dir_group <- dca[categorisation %like% 'close.and.adjacent.and.directed.cat', unique(categorisation)]
idx <- dir_group %like% 'sero'
dir_group <- ifelse(any(idx), dir_group[idx], dir_group[1])

control <- list(linked.group='close.and.adjacent.cat',
                linked.no='not.close.or.nonadjacent',
                linked.yes='close.and.adjacent', 
                dir.group = dir_group,
                conf.cut=0.6, 
                neff.cut=3,
                weight.complex.or.no.ancestry=0.5)
# Find pairs
tmp <- find.pairs.in.networks(dwina, dca, control=control, verbose=TRUE)
dpl <- copy(tmp$network.pairs)
dc <- copy(tmp$relationship.counts)
dw <- copy(tmp$windows)

# Find chains
dc[, unique(CATEGORISATION)]

tmp <- find.networks(dc, control=control, verbose=TRUE)
dnet <- copy(tmp$transmission.networks)
dchain <- copy(tmp$most.likely.transmission.chains)

cat("\n---- Save in output directory (default:dir.zenodo.phyloproc) ----\n")
suff <- ''
if(args$confidential == FALSE)
    suff <- '_randomized'
filename <- paste0('Rakai_phscnetworks_ruleo_sero',suff,'.rda')
filename <- file.path(args$outdir, filename)
if(! file.exists(filename) | config$overwrite.existing.files)
{
    catn(paste0("saving ", filename))
    save(dpl, dc, dw, dnet, dchain, file=filename)
}else{
    catn(paste0("File ", filename, " already exists"))
}
cat('\n\n=====  end of script =====\n\n')
