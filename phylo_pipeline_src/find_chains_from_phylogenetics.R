# The phyloscanner run - find transmission networks 

# Preamble This script aims to find transmission networks and the most likely transmission chains using the phyloscanner outputs.

# Load the required packages
library(data.table)
library(tidyverse)
library(dplyr)
library(glue)
library(igraph)
library(RBGL)
library(phyloscannerR)
library(here)

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
    "--classification_rule",
    type = "character",
    default = 'o',
    help = "Rules for classifying linked and directed pairs. It takes values o or m. 
    o: a pair is linked if the linkage score > 0.6, and directed if the ancestral / (ancestral + descedant) > 0.6.
    m: a pair is linked if the linkage score > 0.5, and directed if the ancestral / all > 0.33. 
    b: both [default]",
    dest = 'classif_rule'
  ),
  optparse::make_option(
    "--out_dir_base",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to base directory where all the tree outputs are stored [default]",
    dest = 'out.dir'
  ),
  optparse::make_option(
    "--phylo_dir",
    type = "character",
    default = '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI/211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd',
    help = "Absolute file path to base directory where the phyloscanner outputs are stored [default]",
    dest = 'phylo.dir'
  ),
  optparse::make_option(
    "--date",
    type = 'character',
    default = '2022-02-04',
    metavar = '"YYYY-MM-DD"',
    help = 'As of date to extract data from.  Defaults to today.',
    dest = 'date'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))

#
# Paths
#

indir <- here()
indir.data <- file.path(indir, 'data')

usr <- Sys.info()[['user']]
if(usr=='andrea')
{
    indir.deepsequencedata <- '/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live'
    indir.deepanalyses.xiaoyue <- '/home/andrea/HPC/project/ratmann_xiaoyue_jrssc2022_analyses/live'
}


#
# Helpers
#

catn <- function(x) cat('\n----', x ,'----\n')


get.sample.collection.dates <- function(select_aid=NULL, get_first_visit=FALSE)
{
    # # get collection dates 
    # path.sdates.rccs <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS','200316_pangea_db_sharing_extract_rakai.csv')
    # path.sdates.mrc <- file.path(indir.deepsequencedata, 'PANGEA2_MRC','200319_pangea_db_sharing_extract_mrc.csv')

    # files <- c(path.sdates.rccs, path.sdates.mrc)
    # cols <- c('pt_id', 'pangea_id', 'visit_dt')
    # ddates <- rbindlist(lapply(files, fread, select=cols))
    # ddates <- unique(ddates)
    # ddates <- merge(ddates, aik, by.x='pt_id', by.y='PT_ID')
    # ddates[, pt_id := NULL]
    # stopifnot(ddates[, uniqueN(pangea_id)==.N])
    # setnames(ddates, 'AID', 'aid')
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


#
# test
#

if(args$confidential)
{
    # paths that should not be available to everyone, but were used for the analysis
    path.meta.data <- .fp('D', 'RCCS_R15_R18/Rakai_Pangea2_RCCS_Metadata_20221128.RData')
    path.sequence.dates <- file.path(indir.confidential, "sequences_collection_dates.rds")
}else{
    # randomized version of the paths used for the analysis
    path.meta.data <- file.path(indir, 'data',"Rakai_Pangea2_RCCS_Metadata_randomized.RData" )
    path.sequence.dates <- file.path(indir.data, "sequences_collection_dates_randomized.rds")
}

dseqdates <- readRDS(path.sequence.dates)

stopifnot("args$phylo.dir does not exist: make sure you specify the correct path"=file.exists(args$phylo.dir))
stopifnot("path.meta.data does not exist: make sure you specify the correct path"=file.exists(path.meta.data))

catn('Load phyloscanner outputs')
stopifnot(dir.exists(args$phylo.dir))
infiles	<- data.table(F=list.files(args$phylo.dir, pattern='*workspace.rda$', full.names=TRUE))
infiles[, PTY_RUN:= as.integer(gsub('^ptyr([0-9]+)_.*','\\1',basename(F)))]
setkey(infiles, PTY_RUN)

# 'dwin' : pairwise relationships between the host in each tree
# 'dc' :  summarises pairwise relationships between hosts across ALL trees. 
# 'a' suffix stands for aggregated

# there must be a slightly faster way to avoid double loading, but nvm
dca <-   infiles[, { cat(PTY_RUN,'\n'); load(F); dc }, by='PTY_RUN']
dwina <- infiles[, { cat(PTY_RUN,'\n'); load(F); dwin }, by='PTY_RUN']

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
meta_data

# 
if( file.exists(path.meta.data) )
{   # Now, whenever a pairwise transmission is not supported by the serohistory, we need to "switch" the counts supporting that direction.
    cat('Using serohistory + demographic data to reweight evidence of direction\n')
    
    # Load anon. keys
    # file.anonymisation.keys <- file.path(indir.deepanalyses.xiaoyue, 'PANGEA2_RCCS1519_UVRI/important_anonymisation_keys_210119.csv')
    # aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

    # load meta data on infection times
    meta_env <- new.env()
    load(path.meta.data, envir=meta_env)
    meta <- subset(meta_env$meta_data,
                   select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
    meta <- unique(meta[!is.na(aid)])
    # idx <- meta[, uniqueN(comm) == 2, by=aid][V1==TRUE, aid]
    stopifnot(meta[, uniqueN(aid) == .N])

    # IFNEEDED save ids herefor which we do not have a first positive date, so we can ask Joseph.

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
    categories.to.update <- c("close.and.adjacent.and.ancestry.cat", "close.and.adjacent.and.directed.cat",
                              "close.and.contiguous.and.ancestry.cat", "close.and.contiguous.and.directed.cat")
    dca_sero_only <- dca[categorisation %in% categories.to.update, ]

    # label.to.update <- function(ctg)
    #   fcase(ctg %like% 'directed', NA_character_, ctg %like% 'ancestry.cat', 'complex.or.no.ancestry')
    dca_sero_only <- update.category.counts(dexclude, dca_sero_only)
    dca_sero_only[, categorisation:=paste0(categorisation, '.sero')]
    dca <- rbind(dca, dca_sero_only)

    # deleted an if(0)

    setkey(dca, host.1, host.2)
}

if(0)
{
    # checking that redistribution of count has beem done in the correct direction
    drange
    tmp <- dca_sero_only[type %in% c('12', '21') & categorisation %like% 'sero', .(host.1, host.2, type, score)]
    tmp <- merge(tmp, drange[, .(host.1=AID, MIN.1=MIN, MAX.1=MAX)], by='host.1')
    tmp <- merge(tmp, drange[, .(host.2=AID, MIN.2=MIN, MAX.2=MAX)], by='host.2')
    tmp[MAX.1 < MIN.2]
    stopifnot(tmp[MAX.1 < MIN.2, mean(score[type=='12'], na.rm=TRUE) == 1 & mean(score[type=='21'], na.rm=TRUE) == 0])
    stopifnot(tmp[MAX.2 < MIN.1, mean(score[type=='12'], na.rm=TRUE) == 0 & mean(score[type=='21'], na.rm=TRUE) == 1])
}

if(0)
{
    # get pairs whose direction was changed according to serohistory
    class <- 'close.and.adjacent.and.directed.cat.sero'
    dca_class <- dca[categorisation == class | categorisation == gsub( '.sero', '', class)]
    # in the majority of cases, scores sum to 0
    summed_scores <- dca_class[, sum( score, na.rm=T), by=c('host.1', 'host.2', 'categorisation')]
    summed_scores[, table(V1)/.N*100 ]
    # in 8 cases, the sum of the scores went from 1 to 0... (why?)
    summed_scores |>
        dcast(host.1 + host.2 ~ categorisation) |> 
        with( table(close.and.adjacent.and.directed.cat, close.and.adjacent.and.directed.cat.sero)) 
    
    # subset to cases where summed scores were always one ? 
    idx <- summed_scores[ , all(V1==1) , by=c('host.1', 'host.2') ][V1 == TRUE, .(host.1, host.2)]

    idx <- dca_class[idx][, uniqueN(score), by=c('host.1', 'host.2', 'type')][V1>1, .(host.1, host.2)]
    stopifnot( uniqueN(idx) == nrow(idx)/2 )
    idx <- unique(idx)

    # when changes, mass is moving towards 0
    dca_class[idx][categorisation %like% 'sero', table(score)]

}

# find pairs according to classification rule and thresholds.
# classification rule o: Oliver Ratmann's

if(args$classif_rule=='o'|args$classif_rule=='b')
{
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

    # filename <- file.path(args$phylo.dir, paste0('Rakai_phscnetworks_allpairs_ruleo_sero.rda'))
    # if(! file.exists(filename))
    # {
    #     save(dpl, dc, dw, file=filename)
    # }

    # Find chains
    tmp <- find.networks(dc, control=control, verbose=TRUE)
    dnet <- copy(tmp$transmission.networks)
    dchain <- copy(tmp$most.likely.transmission.chains)

    suff <- ''
    if(args$confidential == FALSE)
        suff <- 'randomized'
    filename <- paste0('Rakai_phscnetworks_ruleo_sero_randomized.rda')
    filename <- file.path(indir.data, filename)
    if(! file.exists(filename) )
    {
        catn("saving Rakai_phscnetworks_ruleo_sero.rda")
        save(dpl, dc, dw, dnet, dchain, file=filename)
    }else{
        catn("File Rakai_phscnetworks_ruleo_sero.rda already exists")
    }
}
