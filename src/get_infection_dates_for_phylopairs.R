# AIMS:
# compute generation intervals for pairs in a network.

################
# DEPENDENCIES #
################ 
library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)

################
#    PATHS     #
################

usr <- Sys.info()[['user']]
if(usr == 'andrea')
{
    base.hpc <- '/home/andrea/HPC'
    # out.dir <- '/home/andrea/HPC/Documents/Box/2022/genintervals/'
    out.dir <- '/home/andrea/HPC/ab1820/home/projects/2022/genintervals'
    indir.deepsequence_xiaoyue   <- file.path(base.hpc, 'project/ratmann_xiaoyue_jrssc2022_analyses/live/PANGEA2_RCCS1519_UVRI')
    indir.deepsequencedata <- file.path(base.hpc, 'project/ratmann_pangea_deepsequencedata/live')
    indir.deepsequence_analyses   <- file.path(base.hpc, 'project/ratmann_deepseq_analyses/live')
    indir <- '/home/andrea/git/phyloflows'

}else{
    indir.deepsequence_xiaoyue  #  <- TODO
    indir.deepsequencedata #  <- TODO
    indir # <- TODO
}

.fp <- function(C, x)
{
    if(C=='X')
        indir <- indir.deepsequence_xiaoyue
    if(C=='A')
        indir <- indir.deepsequence_analyses
    if(C=='D')
        indir <- indir.deepsequencedata
    file.path(indir, x)
}

file.path.meta <- .fp('D', 'RCCS_R15_R18/Rakai_Pangea2_RCCS_Metadata_20221128.RData')
# file.path.chains.data <- '~/Downloads/Rakai_phscnetworks_ruleo_sero.rda'
file.path.chains.data.old <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')
file.path.chains.data.old <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks_ruleo_sero.rda')
file.path.round.timeline <- .fp('D', 'RCCS_data_estimate_incidence_inland_R6_R18/220903/RCCS_round_timeline_220905.RData')
file.anonymisation.keys <- .fp('X','important_anonymisation_keys_210119.csv')
file.path.tsiestimates <- .fp('A', 'PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T/aggregated_TSI_with_estimated_dates.csv')

file.exists(file.path.meta,
            file.path.chains.data.old,
            file.path.chains.data,
            file.path.round.timeline,
            file.anonymisation.keys,
            file.path.tsiestimates) |> all() |> stopifnot()

################
#   OPTIONS    #
################
# eg: Rscript get_generation_intervals.R --rerun TRUE --RH-infectiousness 1 --RH-duration 1
# eg: Rscript get_generation_intervals.R --rerun TRUE --RH-infectiousness 5 --RH-duration 1/6
option_list <- list(
    optparse::make_option(
        "--rerun",
        type = "logical",
        default = FALSE,
        help = "Flag on whether to re-run rejection sampling [Defaults to FALSE]", 
        dest = 'rerun'
    ),
    optparse::make_option(
        "--only-hetero",
        type = "logical",
        default = TRUE,
        help = "Excludes heterosexual couples from transmission network", 
        dest = 'include.only.heterosexual.pairs'
    ),
    optparse::make_option(
        "--get-round-probabilities",
        type = "logical",
        default = TRUE,
        help = "Computes probability of infection occurring in different rounds", 
        dest = 'get.round.probabilities'
    ),
    optparse::make_option(
        "--only-rccs",
        type = "logical",
        default = TRUE,
        help = "Excludes non-participants of the RCCS", 
        dest = 'include.only.rccs'
    ), optparse::make_option( "--RH-infectiousness",
        type = "numeric",
        default = 1,
        help = "Relative Hazard of transmission during the acute versus chronic infection phase [Defaults to 5]", 
        dest = 'RH.infectiousness'
    ),
    optparse::make_option(
        "--RH-duration",
        type = "numeric",
        default = 1,
        help = "Duration of acute phase, in years [Defaults to 2 months]", 
        dest = 'RH.duration'
    ),
    optparse::make_option( "--sensitivity-no-refinement",
        type = "logical",
        default = FALSE,
        help = "Relative Hazard of transmission during the acute versus chronic infection phase [Defaults to 5]", 
        dest = 'sensitivity.no.refinement'
    )
)

args <-  optparse::parse_args(optparse::OptionParser(option_list = option_list))
print(args)

################
#    HELPERS   #
################

source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'src/utils', 'gi_analysis_functions.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
# source(file.path(indir, 'functions', 'statistics_functions.R'))
# source(file.path(indir, 'functions', 'stan_utils.R'))
# source(file.path(indir, 'functions', 'check_potential_TNet.R'))
find_palette_round()

plot.rectangles <- function(DT, values=c('', 'TSI', 'INTERSECT'), idx=data.table())
{
    .p <- function(...) 
    {
        out <- paste(..., sep='.')
        out <- gsub('^\\.', '', out)
    }
    .c <- function(x)
        fcase(x=='', 'blue', x=='TSI', 'red', x=='INTERSECT','orange')
    .s <- 'SOURCE'; .r <- 'RECIPIENT'

    plot_dt <- double.merge(chain, DT)
    setkey(plot_dt, SOURCE,RECIPIENT)
    setkey(idx, SOURCE,RECIPIENT)
    values[.p(values, 'MIN', .s) %in% names(DT)]

    if(nrow(idx))
        plot_dt <- plot_dt[idx]


    add.rect <- function(pre)
    {
        geom_rect( aes_string(xmin=.p(pre,'MIN',.s),
                              xmax=.p(pre,'MAX',.s),
                              ymin=.p(pre,'MIN',.r),
                              ymax=.p(pre,'MAX',.r)), fill=NA, color=.c(pre) )
    }
    
    g <- ggplot(plot_dt)

    for(v in values)
        g <- g + add.rect(v)

    g <- g + 
        geom_abline(aes(intercept=0, slope=1), linetype='dashed', color='black')  +
        theme_bw() + 
        labs(x='DOI source', y='DOI recipient')
    g
}

################
#     MAIN     #
################

# get.extra.pairs.from.serohistory <- 1
threshold.likely.connected.pairs <- 0.5
get.sero.extra.pairs <- TRUE
build.network.from.pairs <- TRUE

# Load anon. keys
aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

# load meta data
meta_env <- new.env()
load(file.path.meta, envir=meta_env)
meta <- subset(meta_env$meta_data,
               select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
meta <- unique(meta[!is.na(aid)])
stopifnot(meta[, uniqueN(aid) == .N])

# file.path.chains.data <- file.path.chains.data.old
# load chains 
chains_env <- new.env()
load(file.path.chains.data, envir=chains_env)
dchain <- as.data.table(chains_env$dchain)
chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)
if(get.sero.extra.pairs)
    chain <- get.extra.pairs.from.serohistory(dchain, meta)
chain <- subset(chain, select=c('SOURCE', 'RECIPIENT', 'IDCLU'))
setkey(chain, IDCLU, SOURCE)

if(build.network.from.pairs)
{   # study dc instead..
    dpl <- setDT(chains_env$dpl)
    dc <- setDT(chains_env$dc)
    stopifnot(dpl[, .N, by=c('H1', 'H2')][, all(N==1)])
    stopifnot(dpl[,all(H1 < H2)])

    idx <- dpl[SCORE > threshold.likely.connected.pairs, .(H1, H2, SCORE)]
    tmp <- dc[ CATEGORISATION %like% 'close.and.adjacent.and.directed.cat.sero' , .(H1, H2, TYPE, SCORE)]
    stopifnot(tmp[, .N > 0])
    tmp <- tmp[, {z <- which(SCORE>.5); list(TYPE=TYPE[z], SCORE_DIR=SCORE[z])}, by=c('H1', 'H2')]
    dlinkdir <- merge(idx, tmp, by=c('H1', 'H2'), all.x=TRUE)
    stopifnot(tmp[,.N,by=c('H1', 'H2')][, all(N==1)])
    stopifnot(dlinkdir[,.N,by=c('H1', 'H2')][, all(N==1)])

    # get range of dates
    drange <- get.infection.range.from.testing()

    # if add unsupperted but with strong SCORE_DIR
    if(get.sero.extra.pairs)
    {
        tmp <- dpl[SCORE <= threshold.likely.connected.pairs, .(H1, H2, SCORE)]
        tmp <- merge(tmp, drange[,.(H1=AID, MIN.1=MIN, MAX.1=MAX)], by='H1')
        tmp <- merge(tmp, drange[,.(H2=AID, MIN.2=MIN, MAX.2=MAX)], by='H2')
        tmp <- rbind(
            tmp[MAX.1 < MIN.2, `:=` (SCORE_DIR=1, TYPE='12')],
            tmp[MAX.2 < MIN.1, `:=` (SCORE_DIR=1, TYPE='21')]
        )[!is.na(TYPE), .(H1, H2, SCORE, SCORE_DIR, TYPE), ] 
        dlinkdir <- rbind(dlinkdir, tmp)
    }

    # assign SCORE_DIR and check...
    dlinkdir <- merge(dlinkdir, drange[,.(H1=AID, MIN.1=MIN, MAX.1=MAX)], by='H1')
    dlinkdir <- merge(dlinkdir, drange[,.(H2=AID, MIN.2=MIN, MAX.2=MAX)], by='H2')
    dlinkdir[MAX.1 < MIN.2 & is.na(TYPE), `:=` (TYPE='12', SCORE_DIR=1) ]
    dlinkdir[MAX.2 < MIN.1 & is.na(TYPE), `:=` (TYPE='21', SCORE_DIR=1) ]
    dlinkdir[MAX.1 < MIN.2, stopifnot(all(TYPE=='12'))]
    dlinkdir[MAX.2 < MIN.1, stopifnot(all(TYPE=='21'))]
    dlinkdir[, `:=` (MAX.1=NULL, MAX.2=NULL, MIN.1=NULL, MIN.2=NULL)]

    # Assign direction of transmission
    dnewpairs <- rbind(
        dlinkdir[TYPE == '12', .(SOURCE=H1, RECIPIENT=H2, SCORE, SCORE_DIR)],
        dlinkdir[TYPE == '21', .(SOURCE=H2, RECIPIENT=H1, SCORE, SCORE_DIR)]
    ) |> unique() 
    stopifnot(dnewpairs[,.N,by=c('SOURCE', 'RECIPIENT')][, all(N==1)])

    filename <- file.path(out.dir, paste0('221124_study_networks_all.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, ' recipients')]
    p1 <- plot.chains(dnewpairs, size.threshold = 3, ttl='Before Subsetting', sbttl=lab)
    dev.off()

    # Remove non-RCCS participants
    tmp <- nrow(dnewpairs)
    rccs_ids <- meta[, unique(aid)]
    dnewpairs <- dnewpairs[ SOURCE %in% rccs_ids & RECIPIENT %in% rccs_ids ]
    cat('Excluding', tmp - nrow(dnewpairs), 'of', tmp, 'pairs outside of RCCS\n')
    cat(nrow(dnewpairs), 'pairs remaining\n')
    # plot
    filename <- file.path(out.dir, paste0('221124_study_networks_norccs.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, 'recipients')]
    p2 <- plot.chains(dnewpairs, size.threshold = 3, ttl='After removing non-RCCS', sbttl=lab)
    dev.off()

    # Remove homosexual pairs
    dsex <- meta[, .(aid, sex)] |> unique()
    idx <- merge(dnewpairs, dsex[, .(SOURCE=aid, SEX.SOURCE=sex)], all.x=T, by='SOURCE')
    idx <- merge(idx, dsex[, .(RECIPIENT=aid, SEX.RECIPIENT=sex)], all.x=T, by='RECIPIENT')
    # idx[, table(SEX.SOURCE,SEX.RECIPIENT, useNA = 'ifany'),]
    dhomosexualpairs <- idx[ (SEX.SOURCE==SEX.RECIPIENT) , ]
    dhomosexualpairs[, IDCLU:=NULL]
    idx <- idx[!(SEX.SOURCE==SEX.RECIPIENT), .(SOURCE, RECIPIENT)]

    tmp <- nrow(dnewpairs) - nrow(idx)
    cat('Excluding', tmp, 'of', nrow(dnewpairs), 'pairs of homosexual or unknown sex\n')
    setkey(dnewpairs, SOURCE,RECIPIENT)
    setkey(idx,SOURCE,RECIPIENT)
    dnewpairs <- dnewpairs[idx]
    cat(nrow(dnewpairs), 'pairs remaining\n')
    filename <- file.path(out.dir, paste0('221124_study_networks_norccs_nohomosex.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, 'recipients')]
    p3 <- plot.chains(dnewpairs, size.threshold = 3, ttl='After removing non-RCCS + homosexuals', sbttl=lab)
    dev.off()

    # Solve multiple sources by assigning recipient with strongest linkage support
    dnewpairs <- dnewpairs[, {
        z <- which.max(SCORE);
        list(SOURCE=SOURCE[z], SCORE=SCORE[z], SCORE_DIR=SCORE_DIR[z])
    }, by='RECIPIENT']

    # get an idea of how many can be inland
    meta_env <- new.env()
    load(file.path.meta, envir=meta_env)
    cols <- c('aid', 'comm', 'round', 'sample_date')
    dcomms <- subset(meta_env$meta_data, select=cols)
    names(dcomms) <- toupper( names(dcomms) )
    idx <- dcomms[!is.na(AID), .(COMM=paste0(sort(unique(COMM)), collapse='-'), uniqueN(COMM)), by='AID']
    dnewpairs <- double.merge(dnewpairs, idx[, .(AID, COMM)])
    dnewpairs[COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland', {
        cat(.N);
        print(table(COMM.SOURCE, COMM.RECIPIENT));
        NULL
    }]
    dnewpairs[COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland', uniqueN(RECIPIENT)]
    dnewpairs[COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland', table(COMM.RECIPIENT)] |> knitr::kable()
    dnewpairs[COMM.RECIPIENT %like% 'inland', table(COMM.SOURCE)] |> knitr::kable()
    dnewpairs[COMM.RECIPIENT %like% 'inland', uniqueN(RECIPIENT)]
    
    filename <- file.path(out.dir, paste0('221124_study_networks_norccs_nohomosex_sourceattr.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, ' recipients')]
    p3a <- plot.chains(dnewpairs, size.threshold = 3, ttl='After removing non-RCCS + homosexuals + src attribution', sbttl=lab)
    dev.off()

    setcolorder(dnewpairs, c('SOURCE','RECIPIENT'))
    idclus <- dnewpairs |> graph_from_data_frame() |> components()
    idclus <- data.table(RECIPIENT=names(idclus$membership), IDCLU=unname(idclus$membership))
    dnewpairs <- merge(dnewpairs, idclus,by='RECIPIENT')
    setcolorder(dnewpairs, c('SOURCE','RECIPIENT'))

    if(0)
    {
    meta_env <- new.env()
    load(file.path.meta, envir=meta_env)
    cols <- c('aid', 'comm', 'round', 'sample_date')
    dcomms <- subset(meta_env$meta_data, select=cols)
    names(dcomms) <- toupper( names(dcomms) )
    idx <- dcomms[!is.na(AID), .(COMM=paste0(unique(COMM), collapse='-'), uniqueN(COMM)), by='AID']
    dnewpairs <- double.merge(dnewpairs, idx[, .(AID, COMM)])
    
    dnewpairs <- dnewpairs[COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland']
    filename <- file.path(out.dir, paste0('221124_study_networks_norccs_nohomosex_onlyinland.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, 'recipients')]
    p4 <- plot.chains(dnewpairs, size.threshold = 2, ttl='After removing non-RCCS + homosexuals + non-inland', sbttl=lab)
    dev.off()

    tmp1 <- graph_from_data_frame(dnewpairs[, .(SOURCE,RECIPIENT)], directed=TRUE, vertices=NULL) |> components()
    comps <- tmp1$membership
    comps <- data.table(SOURCE=names(comps), COMP=unname(comps))
    tmp1 <- merge(comps, dnewpairs)
    setcolorder(tmp1, c('SOURCE','RECIPIENT'))
    tmp1[, SIZE:=.N, by=COMP]

    # If a recipient has multiple sources, pick the recipient with highest transmission score
    dnewpairs <- tmp1[,  list(SOURCE=SOURCE[which.max(SCORE)]) ,by=RECIPIENT]
    filename <- file.path(out.dir, paste0('221124_study_networks_norccs_nohomosex_onlyinland_maxscore.png'))
    png(filename, width = 600, height = 600)
    lab <- dnewpairs[, uniqueN(SOURCE), by=RECIPIENT][, paste0('Multiple source for ',sum(V1!=1), '/', .N, 'recipients')]
    p5 <- plot.chains(dnewpairs, size.threshold = 2, ttl='After removing non-RCCS + homosexuals + non-inland', sbttl=lab)
    dev.off()
    }
}

if(0)
{
    # compare new and old chains
    tmp_env <- new.env()
    load(file.path.chains.data.old, envir = tmp_env)
    
    dchain_old <- as.data.table(tmp_env$dchain)
    chain_old <- keep.likely.transmission.pairs(dchain_old, threshold.likely.connected.pairs)
    # if(! args$sensitivity.no.refinement)
    #    chain_old <- get.extra.pairs.from.serohistory(dchain_old, meta)
    chain_old <- subset(chain_old, select=c('SOURCE', 'RECIPIENT', 'IDCLU'))
    setkey(chain_old, IDCLU, SOURCE)

    tmp <- merge(chain[, .(SOURCE, RECIPIENT, NEW=TRUE)], chain_old[,  .(SOURCE, RECIPIENT, OLD=TRUE)], all.x=TRUE, all.y=TRUE)
    tmp[, table(!is.na(NEW), !is.na(OLD))] |> knitr::kable()
    setkey(tmp, SOURCE, RECIPIENT)


    # THERE ARE NO PAIRS THAT ARE SWAPPED...
    merge(tmp[NEW==TRUE, .(SOURCE, RECIPIENT)], tmp[OLD==TRUE, .(SOURCE=RECIPIENT, RECIPIENT=SOURCE)])
    merge(tmp[OLD==TRUE, .(SOURCE, RECIPIENT)], tmp[NEW==TRUE, .(SOURCE=RECIPIENT, RECIPIENT=SOURCE)])

}

# exclude non-RCCS
if( args$include.only.rccs & ! build.network.from.pairs )
{
    tmp <- nrow(chain)
    rccs_ids <- meta[, unique(aid)]
    chain <- chain[ SOURCE %in% rccs_ids & RECIPIENT %in% rccs_ids ]
    cat('Excluding', tmp - nrow(chain), 'of', tmp, 'pairs outside of RCCS\n')
    cat(nrow(chain), 'pairs remaining\n')

}

# exclude homosexuals
if(args$include.only.heterosexual.pairs & ! build.network.from.pairs)
{
    dsex <- meta[, .(aid, sex)] |> unique()
    idx <- merge(chain, dsex[, .(SOURCE=aid, SEX.SOURCE=sex)], all.x=T, by='SOURCE')
    idx <- merge(idx, dsex[, .(RECIPIENT=aid, SEX.RECIPIENT=sex)], all.x=T, by='RECIPIENT')
    # idx[, table(SEX.SOURCE,SEX.RECIPIENT, useNA = 'ifany'),]
    dhomosexualpairs <- idx[ (SEX.SOURCE==SEX.RECIPIENT) , ]
    dhomosexualpairs[, IDCLU:=NULL]
    idx <- idx[!(SEX.SOURCE==SEX.RECIPIENT), .(SOURCE, RECIPIENT)]

    tmp <- nrow(chain) - nrow(idx)
    cat('Excluding', tmp, 'of', nrow(chain), 'pairs of homosexual or unknown sex\n')
    setkey(chain, SOURCE,RECIPIENT)
    setkey(idx,SOURCE,RECIPIENT)
    chain <- chain[idx]
    cat(nrow(chain), 'pairs remaining\n')
}

if(build.network.from.pairs)
{
    chain <- copy(dnewpairs)
    drange <- get.infection.range.from.testing()
    check <- merge(dnewpairs, drange[, .(SOURCE=AID, MIN.SOURCE=MIN, MAX.SOURCE=MAX)], by='SOURCE')
    check <- merge(check, drange[, .(RECIPIENT=AID, MIN.RECIPIENT=MIN, MAX.RECIPIENT=MAX)], by='RECIPIENT')
    check[ MAX.RECIPIENT < MIN.SOURCE, stopifnot(.N==0)]
}

if(1)
{
    meta_env <- new.env()
    load(file.path.meta, envir=meta_env)
    cols <- c('aid', 'comm', 'round', 'sample_date')
    dcomms <- subset(meta_env$meta_data, select=cols)
    names(dcomms) <- toupper( names(dcomms) )

    tmp <- dcomms[!is.na(AID), .(COMMS=paste0(unique(sort(COMM)), collapse='-')) , by='AID']
    tmp <- double.merge(chain, tmp)
    tmp <- tmp[ COMMS.RECIPIENT %like% 'inland' & ! COMMS.SOURCE == '']
    tmp[ , table(COMMS.SOURCE)]
    tmp[ , uniqueN(RECIPIENT)]
}

dancestors <- get.ancestors.from.chain(chain)

# get plausible infection ranges.
drange <- get.infection.range.from.testing()
chain <- check.inconsistent.testing(drange, switch_if_no_other_src = TRUE)
drange <- shrink.intervals(drange)

# update using TSI estimates
drange_tsi <- get.infection.range.from.tsi(file.path.tsiestimates)
drange_tsi <- check.inconsistent.testing(drange_tsi, switch_if_no_other_src = FALSE)
drange_tsi <- shrink.intervals(drange_tsi)
setnames(drange_tsi, c('MIN', 'MAX'), c('TSI.MIN', 'TSI.MAX'))

# check no contradictions after shrinkning
idx <- double.merge(chain, drange_tsi, by_col = "AID")[TSI.MAX.RECIPIENT - TSI.MIN.SOURCE < 0, unique(IDCLU)]
stopifnot(nrow(idx) == 0)

# intersect the serohistory and tsi ranges at the individal-level. If intersection is empty, priorities tsi.
tmp <- merge(drange, drange_tsi)
tmp[,{
    rng_intersect <- (MAX>TSI.MIN & TSI.MAX>MIN)

    fifelse(rng_intersect,
            yes=pmax(MIN, TSI.MIN),
            no=MIN
    ) -> r_min
    fifelse(rng_intersect,
            yes=pmin(MAX, TSI.MAX),
            no=MAX
    ) -> r_max
    list(AID, INTERSECT.MIN=r_min, INTERSECT.MAX=r_max)
    }] -> tmp1

if(0)
{   # VISUALISE
    for(i in 1:100)
    {
        idx <- chain[i, .(SOURCE,RECIPIENT)]
        p <- plot.rectangles(tmp, idx=idx, values = c("", "TSI"))
        print(p)
        Sys.sleep(1.5)
        if(i == 10)
            plot.rect <- copy(p)
    }
}

cols <- c('MIN', 'MAX')
setnames(tmp1, paste0('INTERSECT.', cols), cols)
drange <- rbind( drange[! (AID %in% tmp1$AID)],tmp1)


# Specify relative likelihood of infectiousness
# _____________________________________________

cat('specifying relative infectiousness...\n')

# Bellan 2015: 5.3 EHM_\text{acute} d_\text{acute}=1.7months
with(args,
     data.table( START = c(0, RH.duration), VALUE = c(RH.infectiousness, 1))
) -> dinfectiousness
dinfectiousness[, END := c(START[-1], Inf)]
setcolorder(dinfectiousness, c('START', 'END', 'VALUE'))


# Load tranmsission pairs and run MC
# __________________________________

dpairs <- double.merge(chain, drange, by_col = "AID")
dpairs[MAX.RECIPIENT - MIN.SOURCE < 0, stopifnot(.N == 0)]
# dpairs[, sort(MAX.SOURCE-MIN.SOURCE)]
# dpairs[, sort(MAX.RECIPIENT-MIN.RECIPIENT)]

# get transmission cluster ids. 
dclus <- get.transmission.cluster.ids(dpairs, check=FALSE)

# initialise data.table summarising geometric properties
dclus[, { g <- GROUP;
         out <- dclus[GROUP==g, .(SOURCE=SOURCE, RECIPIENT=RECIPIENT)]
         out <- unique(out)
         N_out <- out[, uniqueN(c(SOURCE,RECIPIENT))]
         out <- list(out)
         out <- list(IDS=out, N_IDS=N_out)
         }, by=GROUP] -> dcohords

# specify pdf through infectiousness

# RUN MCMC
# ________

range_gi <- c(.5, 1:10)
tmp <- dinfectiousness[!is.infinite(END), paste(round(END, 2), sep='')]
tmp <- gsub('\\.', '', tmp)

filename_net <- 'networks_GICentroids'
suffix <- ''
if( ! is.null(dinfectiousness))
    suffix <- paste0('_d', tmp, '_w', paste0( dinfectiousness$VALUE, collapse=''))
if( build.network.from.pairs  )
    suffix <- paste0(suffix, '_netfrompairs')
if( get.sero.extra.pairs )
    suffix <- paste0(suffix, '_seropairs')
filename_net <- file.path(out.dir, paste0(filename_net, suffix, '.rds'))

if(args$get.round.probabilities)
{
    
    tmp_env <- new.env()
    load(file.path.round.timeline, envir = tmp_env)
    df_round_inland <- tmp_env$df_round_inland
    df_round_fishing <- tmp_env$df_round_fishing

    start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
    stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
    start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
    stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

    stopifnot(start_first_period_inland < stop_first_period_inland)
    stopifnot(stop_first_period_inland < start_second_period_inland)
    stopifnot(start_second_period_inland < stop_second_period_inland)

    cols <- c("ROUND", "MIN_SAMPLE_DATE", "MAX_SAMPLE_DATE", "INDEX_TIME", "INDEX_ROUND", "round")
    df_round_i <- make.df.round(df_round_inland) |> subset(select=cols)
    df_round_f <- make.df.round(df_round_fishing) |> subset(select=cols)

    cols2 <- cols[! cols %like% 'SAMPLE_DATE|INDEX']
    df_round_gi <- merge(df_round_i, df_round_f, by=cols2, all.x=TRUE)
    df_round_gi[, `:=` (
                        MIN_SAMPLE_DATE=pmin(MIN_SAMPLE_DATE.x, MIN_SAMPLE_DATE.y, na.rm=TRUE),
                        MAX_SAMPLE_DATE=pmax(MAX_SAMPLE_DATE.x, MAX_SAMPLE_DATE.y, na.rm=TRUE),
                        INDEX_TIME=INDEX_TIME.x,
                        INDEX_ROUND=INDEX_ROUND.x
                        )]
    df_round_gi <- subset(df_round_gi, select= ! names(df_round_gi) %like% '.x$|.y$' )
    df_round_gi <- subset(df_round_gi, select= cols)

    df_round_gi[, MAX_SAMPLE_DATE := pmax(MAX_SAMPLE_DATE, shift(MIN_SAMPLE_DATE, -1), na.rm=TRUE)]
    df_round_gi


    rm(tmp_env, df_round_inland, df_round_fishing, df_round_i, df_round_f)
}

if(file.exists(filename_net) & ! args$rerun )
{
    dcohords <- readRDS(filename_net)
}else{

    # get A LOT of Uniform Samples for MCMC
    provide_samples <- lapply( dcohords[, seq(max(N_IDS))], function(x){runif(1e6, 0, 1)})
    n_iter <- nrow(dcohords)

    dcohords[ , {
        cat(GROUP, '/', n_iter, '\n')
        DT <- IDS[[1]]
        out <- get.volume.genints.multid(DT,
                                         N_IDS=N_IDS,
                                         range_gi=range_gi,
                                         verbose=T,
                                         importance_weights=dinfectiousness,
                                         get_volume=T,
                                         provide_samples=provide_samples,
                                         df_round=df_round_gi)
        list(out, NAME=c('CENTROID', 'VOLUME', 'GI', 'RPROBS'))
    } , by=GROUP] -> tmp

    Reduce('merge',
           list(tmp[NAME == 'CENTROID', .(GROUP, CENTROIDS=out)],
                tmp[NAME == 'VOLUME', .(GROUP, VOLUME=unlist(out))],
                tmp[NAME == 'GI', .(GROUP, GENINTS=out)]),
                tmp[NAME == 'RPROBS', .(GROUP, RPROBS=out)]
           ) -> tmp
    dcohords <- merge( dcohords[, .(GROUP, IDS, N_IDS)], tmp, by='GROUP')

    rm(provide_samples, n_iter)
    saveRDS(dcohords, filename_net)
}

if(1)   
{   # compute 95% ranges
    dpred_ranges <- dcohords[, CENTROIDS[[1]] ,by=GROUP]
    stopifnot( dpred_ranges[, .N == uniqueN(ID)] )
    stopifnot( dpred_ranges[, all(IL <= M & M <= IU) ] )
    stopifnot( dpred_ranges[, all(CL <= IL & IU <= CU) ] )

    tmp <- merge(dpred_ranges, drange, by.x='ID', by.y='AID')
    stopifnot( tmp[, all(IL >= MIN)] )
    stopifnot( tmp[, all(IU <= MAX)] )
}

if(nrow(df_round_gi))
{   # compute probability of assigning source to different period
    dprobs_roundallocation <- dcohords[ , RPROBS[[1]], by='GROUP']
    dprobs_roundallocation <- merge(
                                    df_round_gi[, .(ROUND, INDEX_TIME)],
                                    dprobs_roundallocation,
                                    by='ROUND', all.y=TRUE)
    dprobs_roundallocation <- dprobs_roundallocation[, 
            .(P=sum(P)),
            by=c('SOURCE','RECIPIENT','INDEX_TIME')]
    
    setkey(dprobs_roundallocation, RECIPIENT, SOURCE, INDEX_TIME)

    dprobs_roundallocation[is.na(INDEX_TIME), INDEX_TIME := 0 ]

    stopifnot(dprobs_roundallocation[, abs(sum(P) - 1) <= .00001, by=c('SOURCE', 'RECIPIENT')][,all(V1)])

    dprobs_roundallocation <- dcast(dprobs_roundallocation, SOURCE+RECIPIENT ~ INDEX_TIME)
    setnames(dprobs_roundallocation, c('0', '1', '2'), c('BeforeR10', 'R10_R15', 'R16_18'))
    cols <- c('BeforeR10', 'R10_R15', 'R16_18')
    dprobs_roundallocation[, (cols) := lapply(.SD, function(x){x[is.na(x)] <- 0; x} ) , .SDcols=cols]
    
    filename <- file.path(out.dir, paste0('probs_roundallocations', suffix, '.rds'))
    saveRDS(dprobs_roundallocation, filename)
}

# prepare input for Bayesian model
centroids <- dcohords[, CENTROIDS[[1]], by='GROUP']
dresults <- prepare.pairs.input.for.bayesian.model(centroids)
# fsetdiff(dresults[, .(SOURCE, RECIPIENT)], chain[, .(SOURCE, RECIPIENT)])
# dresults[, table(DIRECTION)]
dhomosexualpairs[, DIRECTION := 'phyloscanner']
dhomosexualpairs <- dhomosexualpairs[! RECIPIENT %in% dresults$RECIPIENT]
dhomosexualpairs[, `:=` (SCORE=NULL, SCORE_DIR=NULL)]
rbind(
    dresults,
    dhomosexualpairs,
    fill=TRUE
) -> dresults
dresults <- get.community.type.at.infection.date(dresults)

if(get.sero.extra.pairs)
{
    if(nrow(additional_pairs_from_serohistory))
    {
        setkey(dresults, SOURCE, RECIPIENT)
        idx <- additional_pairs_from_serohistory[, .(SOURCE, RECIPIENT)]
        dresults[idx, DIRECTION := 'serohistory']

        # 
        idx2 <- idx[! RECIPIENT %in% dresults$RECIPIENT]
        dadditional <- meta[, .(SOURCE=aid, SEX.SOURCE=sex, DIRECTION='serohistory')]
        dadditional <- merge(idx2, dadditional)
        tmp1 <- meta[, .(RECIPIENT=aid, SEX.RECIPIENT=sex)]
        dadditional <- merge(dadditional, tmp1, by='RECIPIENT')
        setcolorder(dadditional, c('SOURCE','RECIPIENT', 'SEX.SOURCE', 'SEX.RECIPIENT', 'DIRECTION'))
        rbind(
            dresults,
            dadditional,
            fill=TRUE
        ) -> dresults
    }
}

# get sample collection dates
tmp <- get.sample.collection.dates(get_first_visit=TRUE)
names(tmp) <- toupper(gsub('_','.',names(tmp)))
dresults <- double.merge(dresults, tmp)
stopifnot(dresults[, all(CU < DATE.COLLECTION.RECIPIENT, na.rm=TRUE)])

dresults[!is.na(M) &
         SEX.SOURCE != SEX.RECIPIENT &
         COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland', {
             cat('There are', .N, 'pairs with median DOI estimate(', sum(is.na(ROUND.M)), ') outside of range.\n'); .SD
         } ][is.na(ROUND.M), .(SOURCE, RECIPIENT, M)] -> idx
idx <- double.merge(idx, meta[, .(AID=aid, DFP=date_first_positive)])

setkey(dprobs_roundallocation, SOURCE,RECIPIENT)
setkey(idx, SOURCE,RECIPIENT)
dprobs_roundallocation[idx] |> knitr::kable()




# Save
filename <-file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_toshare', suffix, '.rds')) 
saveRDS(dresults, filename)    

if(0)
{   # study all FF pairs to be sent to Griffin
    idx <- unique(dchain[, .(H1, H2)])
    idx <- merge(idx, meta[, .(H1=aid, SEX.H1=sex)], by='H1')
    idx <- merge(idx, meta[, .(H2=aid, SEX.H2=sex)], by='H2')
    idx <- tmp[SEX.H1 == SEX.H2 & SEX.H1 == 'F', .(H1, H2)]

    cols <- c('SOURCE', 'RECIPIENT', 'SEX.SOURCE', 'SEX.RECIPIENT', 'ROUND.M', 'DIRECTION', 'COMM.SOURCE', 'COMM.RECIPIENT')
    tmp12 <- merge(dresults, idx, by.x=c('SOURCE', 'RECIPIENT'), by.y=c('H1', 'H2'))
    tmp21 <- merge(dresults, idx, by.x=c('SOURCE', 'RECIPIENT'), by.y=c('H2', 'H1'))
    
    tmpUN <- merge( idx, tmp12, by.x=c('H1', 'H2'), by.y=c('SOURCE', 'RECIPIENT'), all.x=TRUE)[is.na(SEX.SOURCE), .(H1, H2) ]
    tmpUN <- merge(tmpUN, tmp21, by.x=c('H1', 'H2'), by.y=c('RECIPIENT', 'SOURCE'), all.x=TRUE)[is.na(SEX.RECIPIENT), .(H1, H2) ]
    
    tmpUN[, `:=`(SOURCE=H1, RECIPIENT=H2, SEX.SOURCE='F', SEX.RECIPIENT='F', ROUND.M=NA_character_, DIRECTION='phyloscanner_unclear')]
    tmpUN[, `:=`(H1=NULL, H2=NULL)]
    stopifnot(nrow(tmpUN) + nrow(tmp12) + nrow(tmp21) == nrow(idx))

    all_ff_pairs <- rbind(tmp12, 
                          tmp21, 
                          tmpUN,
                          fill=TRUE)

    load(file.path.meta, envir=meta_env)
    dcomms <- subset(meta_env$meta_data, select=c('aid', 'comm', 'round'))
    dcomms <- unique(dcomms[!is.na(aid),])

    dcomms <- dcomms[ is.na(comm), comm :='neuro']
    dcomms <- dcomms[, list(comm=fifelse(uniqueN(comm)==1, yes=comm[1], no=NA_character_ )), by='aid']

    all_ff_pairs[is.na(COMM.SOURCE), COMM.SOURCE:=dcomms[aid == SOURCE, comm], by='SOURCE']
    all_ff_pairs[is.na(COMM.RECIPIENT), COMM.RECIPIENT:=dcomms[aid == RECIPIENT, comm], by='RECIPIENT']

    filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('221117_all_ff_pairs.csv'))
    fwrite(all_ff_pairs, filename)
}

# make table
cols <- c('SOURCE', 'RECIPIENT', 'SEX.SOURCE', 'SEX.RECIPIENT', 'M', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT', 'DIRECTION')
dtable <- subset(dresults,select=cols)
dtable[, table(SEX.SOURCE,SEX.RECIPIENT)]
setnames(dtable, 'M', 'INFECTION_DATE')
filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_table_toshare', suffix, '.csv'))
# [1] "/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/pairsdata_table_toshare_d1_w11_netfrompairs_seropairs.csv"
fwrite(dtable, filename)

##############################
#   EXTENDED DATA FIGURES   # 
##############################

library(ggpubr)

# For plots, only select pairs which appear in the analysis
dresults <- dresults[SEX.SOURCE!=SEX.RECIPIENT]
dresults <- dresults[ COMM.SOURCE != 'neuro' & COMM.RECIPIENT != 'neuro'] 
dresults <- dresults[COMM.RECIPIENT == 'inland' & COMM.SOURCE == 'inland', ]
dresults <- dresults[! is.na(ROUND.M) ]
cat(dresults[, .N], 'source-recipient pairs selected\n')
# dresults[, table(COMM.SOURCE)]

naturemed_reqs()

if(1)
{
    # plot final predictions v1
    {
        tmp_aids <- dresults[!is.na(M), .(SOURCE,RECIPIENT, M)][, RECIPIENT[order(M)] ]

        drange_topo <- get.infection.range.from.testing()
        drange_test <- copy(drange_topo)
        shrink.intervals(drange_topo)
        
        drange_final <- copy(drange)
        lvls <- c('Demographic and testing', '+ directionality', '+ phyloTSI')
        drange_test[, TYPE:=lvls[1]]
        drange_topo[, TYPE:=lvls[2]]
        drange_final[, TYPE:=lvls[3]]
        drange_tmp <- rbind(drange_topo, drange_test, drange_final)
        drange_tmp[, TYPE := ordered(TYPE, levels=lvls)]

        if(0)
        {
            p_predictions <- ggplot(dresults[!is.na(M)], aes(y=ordered(RECIPIENT, levels= tmp_aids) )) + 
                geom_point(aes(x=M), color='red')  +
                geom_vline(data=df_round_gi, aes(xintercept=MIN_SAMPLE_DATE), linetype='dotted') +
                geom_vline(data=df_round_gi, aes(xintercept=MAX_SAMPLE_DATE), linetype='dotted') +
                geom_vline(data=df_round_gi[round==15], aes(xintercept=MAX_SAMPLE_DATE), color='red') +
                geom_rect(data = df_round_gi, aes(y=NULL, ymin = first(tmp_aids), ymax = last(tmp_aids), 
                                                xmin = MIN_SAMPLE_DATE, xmax = MAX_SAMPLE_DATE,
                                                fill = as.ordered(round)), alpha = 0.5) + 
                geom_linerange(data=drange_tmp[AID %in% tmp_aids], aes(xmin=MIN, xmax=MAX, y=AID, color=TYPE)) +
                geom_point(aes(x=M), color='red')  +
                scale_fill_manual(values=palette_round) +
                scale_color_manual(values=c('grey80', 'grey30', 'black') )  +
                guides(color=guide_legend(nrow=1, byrow=TRUE, override.aes=list(size=3))) + 
                guides(fill=guide_legend(nrow=1, byrow=TRUE, override.aes=list(size=1))) + 
                theme_bw() + 
                theme(panel.border = element_blank(), 
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      axis.line = element_line(colour = "black")) +
                theme(legend.position='bottom', legend.box="vertical", 
                      axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
                labs(x='Date of infection', y='Recipient', fill='Round', pch='median doi', color='Infection range')

        }

        tmp <- unique(dresults[, .(RECIPIENT, CL, M, CU)])

        cols <- c('RENAME_ID','pred_doi_mid',  'pred_doi_min', 'pred_doi_max')
        drange_tsi2 <- fread(file.path.tsiestimates, select = cols) 
        cols_new <- c('RECIPIENT', 'MID','MIN', 'MAX')
        setnames(drange_tsi2, cols, cols_new)
        drange_tsi2[, `:=` (RECIPIENT=gsub('-fq[0-9]', '', RECIPIENT)) ] 
        cols_new <- setdiff(cols_new, 'RECIPIENT')
        drange_tsi2[, (cols_new) := lapply(.SD, as.Date) , .SDcols=cols_new]
        drange_tsi2 <- drange_tsi2[RECIPIENT %in% tmp$RECIPIENT]

        tmp <- merge(tmp, drange_tsi2, by='RECIPIENT', all.x=TRUE)
        tmp[, INTERSECT := pmax(CL, MIN) <= pmin(CU, MAX)]
        tmp <- tmp[!is.na(MIN)]


        p_predictions2 <- ggplot(tmp) + 
            geom_rect(data=df_round_gi, aes(xmin=MIN_SAMPLE_DATE, xmax=MAX_SAMPLE_DATE, 
                                         ymin=MIN_SAMPLE_DATE, ymax=MAX_SAMPLE_DATE, fill=as.ordered(round))) + 
            geom_rect(data=df_round_gi, aes(xmin=as.Date(-Inf), xmax=MAX_SAMPLE_DATE , 
                                        ymin=MIN_SAMPLE_DATE, ymax=MAX_SAMPLE_DATE, 
                                         fill=as.ordered(round)), alpha=.1) + 
            geom_rect(data=df_round_gi, aes(ymin=as.Date(-Inf), ymax=MAX_SAMPLE_DATE , 
                                        xmin=MIN_SAMPLE_DATE, xmax=MAX_SAMPLE_DATE, 
                                        fill=as.ordered(round)), alpha=.1) + 
            geom_abline(aes(slope=1, intercept=0), linetype='dashed', color='red' ) + 
            geom_point(aes(x=MID, y=M), size=.5) + 
            geom_errorbar(aes(x=MID ,ymin=CL, ymax=CU), alpha=.1) + 
            geom_errorbarh(aes(y=M ,xmin=MIN, xmax=MAX), alpha=.1) + 
            scale_fill_manual(values=palette_round) +
            guides(fill=guide_legend(ncol=1, override.aes=list(size=1))) + 
            theme_bw() + 
            labs(x="infection time estimates using phyloTSI on deep-sequence data",
                 y="refined infection time estimates accounting further for sero-history and transmission direction",
                 fill='Round', alpha='') + 
            theme(legend.position='right', legend.direction='vertical', 
                  panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black"))
        force(p_predictions2)

        # get MAE for people with last negative test
        # Median absolute error for 
        dsero <- meta[!is.na(date_first_positive) & !is.na(date_last_negative), 
                      .(AID=aid, MIN=as.Date(date_last_negative), MAX=as.Date(date_first_positive))]
        dsero <- dsero[, list( midpoint=mean(c(MIN, MAX)) ), by='AID']

        tmp1 <- merge(tmp, dsero, by.x='RECIPIENT', by.y='AID')
        tmp1[, final := as.numeric(abs(M - midpoint)/365.25)]
        tmp1[, phyloTSI := as.numeric(abs(MID - midpoint)/365.25)]
        tmp2 <- melt( tmp1[, .(RECIPIENT, final, phyloTSI)],
                     id.vars='RECIPIENT', variable.name='METHOD', value.name='AE') 
        means <- tmp2[, .(MEAN=mean(AE)), by='METHOD']

        .rm <- function(x)
            fifelse(x %like% 'final|Final',
                    yes='refined infection time estimates accounting\nfurther for serohistory and transmission direction',
                    no='phyloTSI on\ndeep sequence data')

        change_histogram_colors <- function(p)
        {
            p + scale_color_manual(values=c('red', 'blue')) + 
                scale_fill_manual(values=c('red', 'blue')) + 
                scale_alpha_manual(values=c(.4, .5))
        }

        plot_mae <- ggplot(tmp2, aes(color=.rm(METHOD), fill=.rm(METHOD))) + 
            geom_histogram(aes(AE, alpha=.rm(METHOD)), center=.25, binwidth=.5,
                           position='identity' ) + 
            geom_point(data=means, aes(y=-150/60, x=MEAN, color=.rm(METHOD)), pch=2, size=2)+ 
            # geom_vline(data=means, aes(xintercept=MEAN, color=METHOD)) +
            theme_bw() +
            theme(legend.position='bottom') + 
            scale_x_continuous(expand = c(0.01, .1), breaks = 1:20) +
            scale_y_continuous(limits=c(-150/30, 150), expand = c(0, 0)) +
            labs(fill='', color='', alpha='',
                 x='absolute difference between infection time estimates and\nthe midpoint of the seroconversion interval',
                 y='source-recipient pairs for whom recipient has a last negative test') 

        plot_mae <- change_histogram_colors(plot_mae)
        force(plot_mae)

        # get ages and rounds if TSI criterion is used
        tsi_predictions <- copy(tmp)
        tsi_predictions <- tsi_predictions[ , .(ID=RECIPIENT, CL=MIN, M=MID, CU=MAX)]

        dresults[ ! RECIPIENT %in% tsi_predictions$ID]
        dresults_tsi <- prepare.pairs.input.for.bayesian.model(tsi_predictions)

        plot.age.comparison <- function(part)
        {
            stopifnot(part %in% c('SOURCE', 'RECIPIENT'))
            part_age_col <- fifelse(part=='SOURCE', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT')

            part_age_col %in% names(dresults)
            part_age_col %in% names(dresults_tsi)

            cols <- c(part, part_age_col)
            .f <- function(DT) subset(DT, select=cols)
            dage_comparison <- list(.f(dresults_tsi), .f(dresults))
            

            stopifnot(part %in% names(dage_comparison[[1]]))
            stopifnot(part %in% names(dage_comparison[[2]]))

            dage_comparison <- merge(dage_comparison[[1]], dage_comparison[[2]], all.y=TRUE, by=part)

            if(part=='RECIPIENT')
            {
                dage_comparison <- melt(dage_comparison, id.vars='RECIPIENT',
                                        value.name='AGE_INFECTION.RECIPIENT',
                                        variable.name='METHOD')
            }else{
                dage_comparison <- melt(dage_comparison, id.vars='SOURCE',
                                        value.name='AGE_TRANSMISSION.SOURCE',
                                        variable.name='METHOD')
            }
            dage_comparison[METHOD %like% '.x$', METHOD := 'phyloTSI' ]
            dage_comparison[METHOD %like% '.y$', METHOD := 'final' ]

            if(part=='RECIPIENT')
            {
                medians2 <- dage_comparison[, .(MEDIAN=median(AGE_INFECTION.RECIPIENT, na.rm=TRUE)), by='METHOD']
            }else{
                medians2 <- dage_comparison[, .(MEDIAN=median(AGE_TRANSMISSION.SOURCE, na.rm=TRUE)), by='METHOD']
            }

            setnames(dage_comparison, part_age_col,'X')

            plot_age_comparison <- ggplot(dage_comparison, aes(color=.rm(METHOD), fill=.rm(METHOD))) + 
                geom_histogram(aes(x=X),
                               center=.5, binwidth=1,
                               position='identity', alpha=.5) +
                geom_point(data=medians2, aes(y=-25/60, x=MEDIAN, color=.rm(METHOD)), size=2, pch=2)+ 
                theme_bw() +
                theme(legend.position='bottom') + 
                scale_x_continuous(expand = c(0.01, .1), breaks=seq(0, 100, 5)) +
                scale_y_continuous(limits=c(-25/30, 25), expand = c(0, 0)) +
                labs(fill='', color='', alpha='',
                     y = "source-recipient pairs",
                     x = paste("estimated age of the phylogenetically likely\n",tolower(part),"at time of infection"))

            plot_age_comparison <- change_histogram_colors(plot_age_comparison)
            plot_age_comparison
        }
        plot_age_comparison_source <- plot.age.comparison(part='SOURCE')
        plot_age_comparison_recipient <- plot.age.comparison(part='RECIPIENT')
        plot_age_comparison_source
    }


    # plot seroconverters
    {
        cols <- c('RENAME_ID', 'pred_doi_min', 'pred_doi_max', 'visit_dt')
        drange_tsi2 <- fread(file.path.tsiestimates, select = cols) 
        cols_new <- c('AID', 'MIN', 'MAX', 'visit_dt')
        setnames(drange_tsi2, cols, cols_new)
        cols_new <- setdiff(cols_new, 'AID')
        drange_tsi2[, (cols_new) := lapply(.SD, as.Date) , .SDcols=cols_new]
        drange_tsi2[, `:=` (AID=gsub('-fq[0-9]', '', AID)) ] 
        dsero <- meta[!is.na(date_first_positive) & !is.na(date_last_negative), 
                      .(AID=aid, MIN=date_last_negative, MAX=date_first_positive, TYPE='Testing')]

        # note: 157 seroconverters do not appear in the tsi predictions although they appear in db sharing....
        # maybe those are the 'poorly' sequenced? But then why do they have an AID????
        idx <- dsero[! AID %in% drange_tsi2$AID , AID]
        
        dvisits <- drange_tsi2[, .(AID, visit_dt)]
        drange_tsi2 <- drange_tsi2[, -'visit_dt']
        drange_tsi2[, TYPE := 'HIV-PhyloTSI']

        # intersect
        idx <- intersect(drange_tsi2$AID, dsero$AID)
        tmp <- rbind( drange_tsi2[AID %in% idx, ], dsero[AID %in% idx, ])

        setorder(dvisits, visit_dt)
        dvisits[, AID := ordered(AID, levels=dvisits$AID)]
        tmp[, AID := ordered(AID, levels=dvisits$AID )]
        setkey(tmp, AID)
        setorder(tmp, AID)
        
        # median lenghts of intervals...
        tmp[, as.numeric(median(MAX - MIN)/365.25), by='TYPE']
        # proportion of intersecting...
        tmp[, max(MIN) <= min(MAX), by='AID'][, paste0(round(mean(V1)*100, 2), '%') ]
        

        p_tsisero <- ggplot(tmp, aes(y=AID)) +
            geom_point(data=dvisits[AID %in% tmp$AID], aes(x=visit_dt)) + 
            geom_linerange(aes(xmin=MIN, xmax=MAX, color=TYPE), alpha=.5) + 
            scale_color_manual(values=c('red', 'blue'))  +
            theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black")) +
            theme(legend.position='bottom', axis.ticks.y=element_blank(), axis.text.y=element_blank()) +
            labs(x='Date of infection', y='Seroconverters', fill='Round', pch='median doi', color='Infection range')
        p_tsisero
    }

    # Plot schema
    for(i in 101:200)
    {
        i <- 129
        cat(i, '\n')
        idx <- dresults[i, .(SOURCE,RECIPIENT)]
        # IDX <- copy(idx)

        tmp <- merge(idx, drange_test[, .(AID, MIN.SOURCE=MIN, MAX.SOURCE=MAX)], by.x='SOURCE', by.y='AID')
        tmp <- merge(tmp, drange_test[, .(AID, MIN.RECIPIENT=MIN, MAX.RECIPIENT=MAX)], by.x='RECIPIENT', by.y='AID')
        cols <- names(tmp)[names(tmp) %like% '^MIN|^MAX'] 
        setnames(tmp, cols, paste0('TESTING.', cols))

        tmp1 <- merge(idx,  drange_tsi2[, .(AID, MIN.SOURCE=MIN, MAX.SOURCE=MAX)], by.x='SOURCE', by.y='AID')
        tmp1 <- merge(tmp1, drange_tsi2[, .(AID, MIN.RECIPIENT=MIN, MAX.RECIPIENT=MAX)], by.x='RECIPIENT', by.y='AID')
        cols <- names(tmp1)[names(tmp1) %like% '^MIN|^MAX'] 
        setnames(tmp1, cols, paste0('TSI.', cols))
        tmp1

        tmp2 <- merge(idx,  drange_tsi[, .(AID, MIN.SOURCE=TSI.MIN, MAX.SOURCE=TSI.MAX)], by.x='SOURCE', by.y='AID')
        tmp2 <- merge(tmp2, drange_tsi[, .(AID, MIN.RECIPIENT=TSI.MIN, MAX.RECIPIENT=TSI.MAX)], by.x='RECIPIENT', by.y='AID')
        cols <- names(tmp2)[names(tmp2) %like% '^MIN|^MAX'] 
        setnames(tmp2, cols, paste0('TSI.', cols))

        trap <- tmp[,list(
                          X=c(TESTING.MIN.RECIPIENT, TESTING.MIN.SOURCE, TESTING.MIN.SOURCE, TESTING.MAX.RECIPIENT),
                          Y=c(TESTING.MIN.RECIPIENT, TESTING.MIN.RECIPIENT, TESTING.MAX.RECIPIENT, TESTING.MAX.RECIPIENT)
        )]
        trap1 <- tmp1[,list(
                          X=c(TSI.MIN.SOURCE, TSI.MIN.SOURCE, TSI.MAX.RECIPIENT),
                          Y=c(TSI.MIN.SOURCE, TSI.MAX.RECIPIENT, TSI.MAX.RECIPIENT)
        )]
        trap2 <- tmp2[,list(
                          X=c(TSI.MIN.SOURCE, TSI.MIN.SOURCE, TSI.MAX.RECIPIENT),
                          Y=c(TSI.MIN.RECIPIENT, TSI.MAX.RECIPIENT, TSI.MAX.RECIPIENT)
        )]
        triangle <- data.table(X=c(tmp1$TSI.MIN.SOURCE, tmp1$TSI.MIN.SOURCE, tmp$TESTING.MAX.RECIPIENT),
                               Y=c(tmp1$TSI.MIN.SOURCE, tmp$TESTING.MAX.RECIPIENT, tmp$TESTING.MAX.RECIPIENT))

        .f <- function(x) sort(unique(x))
        center.of.mass <- triangle[, lapply(.SD,.f ),]
        center.of.mass <- center.of.mass[,.( 
                                         X = weighted.mean(X, c(2,1)),
                                         Y = weighted.mean(X, c(1,2))
                                         )]

        .d <- as.Date(-Inf)


        p_schema <- ggplot(data=tmp) +
            geom_rect(data=tmp,
                      aes(xmin=TESTING.MIN.SOURCE, xmax=TESTING.MAX.SOURCE,
                          ymin=TESTING.MIN.RECIPIENT, ymax=TESTING.MAX.RECIPIENT),
                      fill=NA, color='blue', linetype='dotted'
                      ) +
            geom_polygon(data=trap, aes(x=X, y=Y),
                         alpha=.2, fill='blue') + 
            geom_linerange(data=tmp, aes(y=as.Date(-Inf),
                                         xmin=TESTING.MIN.SOURCE,
                                         xmax=TESTING.MAX.RECIPIENT),
                           color='blue', size=3, alpha=.4) +
            geom_linerange(data=tmp, aes(x=as.Date(-Inf),
                                         ymin=TESTING.MIN.RECIPIENT,
                                         ymax=TESTING.MAX.RECIPIENT),
                           color='blue', size=3, alpha=.4) +
            geom_rect(data=tmp1,
                      aes(xmin=TSI.MIN.SOURCE, xmax=TSI.MAX.SOURCE,
                          ymin=TSI.MIN.RECIPIENT, ymax=TSI.MAX.RECIPIENT),
                      fill=NA, color='red', linetype='dotted'
                      ) +
            geom_linerange(data=tmp1, aes(y=as.Date(-Inf),
                                         xmin=TSI.MIN.SOURCE,
                                         xmax=TSI.MAX.RECIPIENT),
                           color='red', size=3, alpha=.4) +
            geom_linerange(data=tmp1, aes(x=as.Date(-Inf) + 1000,
                                         ymin=TSI.MIN.RECIPIENT,
                                         ymax=TSI.MAX.RECIPIENT),
                           color='red', size=3, alpha=.4) +
            geom_polygon(data=trap1, aes(x=X, y=Y),
                         alpha=.2, fill='red') + 
            geom_polygon(data=triangle, aes(x=X, y=Y),
                         alpha=.2, color='purple', size=2) + 
            geom_point(data=center.of.mass, aes(x=X, y=Y), 
                       color='purple', size=5, pch=4) + 
            geom_abline(slope=1, linetype='dashed') +
            theme_bw() +
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  axis.line = element_line(colour = "black")) +
            labs(x='Source date of infection', y='Recipient date of infection') +
            reqs + coord_flip()

        print(p_schema)
        #filename <- file.path(out.dir, 'supp_triangle_schema.pdf')
        #ggsave(filename, p_schema, width=10, height=10, units = 'cm')
        Sys.sleep(1)
    }

    tmp <- ggarrange_nature(plot_mae + theme(legend.spacing.x=unit(.3, 'cm')),
                            plot_age_comparison_source, plot_age_comparison_recipient,
                            common.legend = TRUE, legend = 'bottom',
                            labels = c('b','c','d'), ncol=3)
    tmp <- ggarrange_nature(p_predictions2, tmp, ncol=1, heights = c(6,4), labels=c('a', '')) 
    tmp
    filename <- file.path(out.dir, 'edf_tsis_v_v2.pdf')
    ggsave_nature(filename, tmp, add_reqs=FALSE)

    
    if(0)
    {
        # vertical layout
        p_tmp <- ggarrange(p_schema , p_tsisero, common.legend = TRUE, legend='bottom')
        p1 <- ggarrange(p_predictions, p_tmp, ncol=1, heights=c(.6, .4))
        # horizontal layour
        p_tmp_h <- ggarrange(p_schema+coord_flip() + reqs , p_tsisero+reqs, ncol=1, 
                             common.legend = TRUE, legend='bottom',
                             heights=c(.4, .6),
                             labels=c('b','c'), font.label=list(size=8))
        p_hor <- ggarrange(p_predictions+reqs, p_tmp_h,
                           ncol=2, widths=c(.6, .4),
                           labels=c('a', ''), font.label=list(size=8))

        filename <- file.path(out.dir, 'edf_tsis_v.pdf')
        ggsave(filename, p_hor, width=17, height=18, units = 'cm')


        p_tmp_h <- ggarrange(plot_mae + reqs, plot_age_comparison+reqs,
                             common.legend = TRUE, legend='right',
                             labels=c('b','c'), font.label=list(size=8))

        .v_leg <- theme(legend.position='right', legend.direction='vertical', legend.box='vertical')
        p_hor2 <- ggarrange(p_predictions2 + reqs +
                            .v_leg + 
                            guides(fill=guide_legend(ncol=1, byrow=TRUE, override.aes=list(size=1))),
                        p_tmp_h, 
                            ncol=1, heights=c(6.5,3.5), 
                            labels=c('a', ''), font.label=list(size=8))

        filename <- file.path(out.dir, 'edf_tsis_v_v2.pdf')
        ggsave(filename, p_hor2, width=17, height=18, units = 'cm')
    }

}


if(0)
{
    # Compare against previous ones
    setkey(centroids, ID)
    setkey(dpairs, RECIPIENT)
    dpairs2 <- merge(dpairs,  centroids[, .(DOI2.RECIPIENT=M, RECIPIENT=ID)], by='RECIPIENT')
    dpairs2 <- merge(dpairs2, centroids[, .(DOI2.SOURCE=M, SOURCE=ID)], by='SOURCE')

    drounds <- load.rounds()[COMM=='aggregated']
    drounds[, max_sample_date := shift(min_sample_date, -1)]
    drounds[round=='R018', max_sample_date := as.Date('2018-06-01') ]

    dclass <- data.table(
                         uniform=dpairs2[, classify.transmission.round.given.centroid(DOI.RECIPIENT)],
                         weighted=dpairs2[, classify.transmission.round.given.centroid(DOI2.RECIPIENT)] 
    )
# knitr::kable(table(dclass))

    plot.scatter.comparison.vs.uniform(dpairs2)

    p1 <- plot.phylopair.dates.scores(dpairs2, doi.center.var = "DOI" )
    p2 <- plot.phylopair.dates.scores(dpairs2, doi.center.var = "DOI2" )
    p <- ggpubr::ggarrange(p1, p2)
    filename <- file.path(out.dir,'compare_uniform_weighted_dois',suffix,'.png')
    ggsave(filename, p, width=12, height=8)
}

if(1)   # study generation intervals
{       
    dgen <- dcohords[, GENINTS[[1]], by='GROUP']
    dgen <- dgen[, mean(V1), by=GI]
    dgen[GI %in% c(0.5, 1), {
            cat('We estimate', round(100*V1, 2), '% generation intervals to be less than', GI, 'years\n')
            V1
    }, by='GI']
}


if(0)   # Check whether MCMC performs well in 2d case: YES
{
        dcohords[ N_IDS == 2, {
                cat(GROUP, '\n')
                DT <- IDS[[1]]
                out <- get.volume.genints.multid(DT, N_IDS=N_IDS, range_gi=range_gi, verbose=T)
                list(out, NAME=c('CENTROID', 'VOLUME', 'GI'))
        } , by=GROUP] -> tmp
        Reduce('merge',
                list(tmp[NAME == 'CENTROID', .(GROUP, CENTROIDS=out)],
                     tmp[NAME == 'VOLUME', .(GROUP, VOLUME=unlist(out))],
                     tmp[NAME == 'GI', .(GROUP, GENINTS=out)])
        ) -> tmp

        cols <- c('GROUP', 'CENTROIDS', 'VOLUME', 'GENINTS')
        tmp0 <- subset(dcohords[ N_IDS == 2], select=cols)

        dcomp <- rbind(tmp[, TYPE:='MCMC'], tmp0[, TYPE:='EXACT'])
        # compuate volume errors
        dcomp[, (VOLUME[TYPE=='MCMC'] - VOLUME[TYPE=='EXACT'])/VOLUME[TYPE=='EXACT'] , by=GROUP]
        # compute GENINTS errors
        dcomp[GROUP == 1, (GENINTS[[1]]$PR - GENINTS[[2]]$PR),  by='GROUP']
        dcomp[, {
                z0 <- GENINTS[[1]]$PR
                z1 <- GENINTS[[2]]$PR
                max((z0 - z1)/z1)
        },  by='GROUP']
        # compute CENTROIDS errors
        dcomp[, {
                z0 <- unlist(CENTROIDS[[1]])
                z1 <- unlist(CENTROIDS[[1]])
                norm(z0-z1, type='2')
        }, by='GROUP']
}

if(0)   
{   # interesting case with source RK-B110064
    range_gi
    tmp0 <- dcohords[ N_IDS == 6 & GROUP==129, GENINTS[[1]] ]
    tmp0[, list(X=GI-.5, Y=diff(c(0, PR))), by=c('SOURCE','RECIPIENT')] |> 
            ggplot(aes(X, Y, color=RECIPIENT)) +
                    geom_line() + 
                    theme_bw()

    setkey(dpairs, SOURCE, RECIPIENT)
    idx <- dcohords[GROUP == 129, IDS[[1]]]
    dpairs[idx, ..cols0]

    cols <- names(dpairs)[names(dpairs) %like% 'REC|date_first_positive.SOURCE']
    cols0 <-names(dpairs)[names(dpairs) %like% '^SOURCE|^RECIPIENT|MAX|MIN']
    dpairs[RECIPIENT %like% 'G109883' , ]
    dpairs[SOURCE %like% 'B110064' , ..cols0]
    dpairs[RECIPIENT %like% 'G109883', ..cols]
    dpairs[RECIPIENT %like% 'G109883', ..cols0]
}

if(0)   
{   # Look at what is the minimum upper bound for generation time
    dgen[PR == 0, range(GI)]
    dgen[PR == 0 & GI >= 8 ] 

    idx <- dgen[PR == 0 & GI >= 8, unique(.(SOURCE=SOURCE, RECIPIENT=RECIPIENT))]

    setkey(dpairs, SOURCE, RECIPIENT)
    cols0 <-names(dpairs)[names(dpairs) %like% '^SOURCE|^RECIPIENT|MAX|MIN']
    dpairs[idx, ..cols0]

    dpairs[ , range(MIN.RECIPIENT - MAX.SOURCE)/365.25]
}


# sensitivity analysis
generation.interval.sensitiviy.analysis <- function()
{

    files <- list.files(file.path(indir.deepsequencedata, 'RCCS_R15_R18'), 
                        pattern='pairsdata.*rds$', full.names = TRUE)
    lresults <- lapply(files, readRDS)
    lapply(lresults, function(DT)
           {    # Subset Results to source-recipient pairs in the analysis
               DT <- DT[SEX.SOURCE!=SEX.RECIPIENT]
               DT <- DT[ COMM.SOURCE != 'neuro' & COMM.RECIPIENT != 'neuro'] 
               DT <- DT[COMM.RECIPIENT == 'inland', ]
               DT <- DT[! is.na(ROUND.M) ]
               cat(DT[, .N], 'source-recipient pairs selected\n')
               DT
           }) -> lresults

    # get cols of interest
    cols <- c('SOURCE', 'RECIPIENT', 'CL', 'IL', 'M', 'IU', 'CU', 'ROUND.M')
    lresults <- lapply(lresults, subset, select=cols)
    names(lresults) <- fifelse(basename(files) %like% 'w11', 'uniform', 'non-uniform')

    dsens <- merge(lresults[[1]], lresults[[2]], by=c('SOURCE', "RECIPIENT"))
    names(dsens) <- gsub('\\.x$', '.uniform', names(dsens))
    names(dsens) <- gsub('\\.y$', '.bellan', names(dsens))
    names(dsens)
#  [1] "SOURCE"     "RECIPIENT"  "CL.uniform" "IL.uniform" "M.uniform" 
#  [6] "IU.uniform" "CU.uniform" "CL.bellan"  "IL.bellan"  "M.bellan"  
# [11] "IU.bellan"  "CU.bellan" 
    .gs <- function(x) gsub('R0', 'R ', x)
    table(dsens[, .(ROUND.M.uniform, ROUND.M.bellan)]) |> 
        xtable( label='t:sens_analysis_bellan_genints') 


    ggplot(dsens) +
        geom_errorbar(aes(xmin=CL.uniform, xmax=IU.uniform, y=M.bellan)) + 
        geom_errorbar(aes(ymin=CL.bellan, ymax=IU.bellan, x=M.uniform)) + 
        geom_abline
        theme_bw() +
        labs()
}
