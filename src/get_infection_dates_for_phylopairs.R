# AIMS:
# compute generation intervals for pairs in a network.

################
# DEPENDENCIES #
################ 
library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)
library(igraph)

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
file.path.chains.data <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks_ruleo_sero.rda')
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
# eg: Rscript get_infection_dates_for_phylopairs.R --rerun TRUE --RH-infectiousness 1 --RH-duration 1
# eg: Rscript get_infection_dates_for_phylopairs.R --rerun TRUE --RH-infectiousness 5 --RH-duration  0.1666667
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

################
#     MAIN     #
################

# settings 
threshold.likely.connected.pairs <- 0.5
threshold.direction <- 0.5
get.sero.extra.pairs <- FALSE
build.network.from.pairs <- TRUE
postpone.samesex.removal <- TRUE

# initialise overleaf substitute-expressions.
overleaf_expr <- list()

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

if(! build.network.from.pairs)
{
    chain <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)
    cat('There are', chain[, .N], 'pairs with linkage score greater than the threshold',threshold.likely.connected.pairs,'.\n')
    if(get.sero.extra.pairs)
        chain <- get.extra.pairs.from.serohistory(dchain, meta)
    chain <- subset(chain, select=c('SOURCE', 'RECIPIENT', 'IDCLU'))
    setkey(chain, IDCLU, SOURCE)

}else{

    # prepare helpers
    dcomms <- get.communities.where.participated()

    # study dc instead..
    cat('\n Building Network from pairs...\n\n')

    # get couples and chains from run
    dpl <- setDT(chains_env$dpl)
    dc <- setDT(chains_env$dc)
    serohistory_impact <- summarise_serohistory_impact_on_pairs(dc)
    stopifnot(dpl[, .N, by=c('H1', 'H2')][, all(N==1)])
    stopifnot(dpl[,all(H1 < H2)])

    # Merge linkage and direction score and subset to linkage score > threshold
    idx <- dpl[SCORE > threshold.likely.connected.pairs, .(H1, H2, SCORE)]
    dir_scores <- dc[ CATEGORISATION %like% 'close.and.adjacent.and.directed.cat.sero' , .(H1, H2, TYPE, SCORE)]
    dir_scores <- dir_scores[, {z <- which(SCORE>threshold.direction); list(TYPE=TYPE[z], SCORE_DIR=SCORE[z])}, by=c('H1', 'H2')]
    dlinkdir <- merge(idx, dir_scores, by=c('H1', 'H2'), all.x=TRUE)
    stopifnot(dir_scores[,.N,by=c('H1', 'H2')][, all(N==1)])
    stopifnot(dlinkdir[,.N,by=c('H1', 'H2')][, all(N==1)])
    

    # get range of dates
    drange <- get.infection.range.from.testing()

    # if add unsupperted but with strong SCORE_DIR
    if(get.sero.extra.pairs)
    {
        sero_extra_pairs <- dpl[SCORE <= threshold.likely.connected.pairs, .(H1, H2, SCORE)]
        sero_extra_pairs <- merge(sero_extra_pairs, drange[,.(H1=AID, MIN.1=MIN, MAX.1=MAX)], by='H1')
        sero_extra_pairs <- merge(sero_extra_pairs, drange[,.(H2=AID, MIN.2=MIN, MAX.2=MAX)], by='H2')
        sero_extra_pairs <- rbind(
            sero_extra_pairs[MAX.1 < MIN.2, `:=` (SCORE_DIR=1, TYPE='12')],
            sero_extra_pairs[MAX.2 < MIN.1, `:=` (SCORE_DIR=1, TYPE='21')]
        )[!is.na(TYPE), .(H1, H2, SCORE, SCORE_DIR, TYPE), ] 
        sero_extra_pairs <- unique(sero_extra_pairs)
        dlinkdir <- rbind(dlinkdir, sero_extra_pairs)
    }

    # check direction is consistent with the serohistory (btw we subset to RCCS here)
    dlinkdir <- merge(dlinkdir, drange[,.(H1=AID, MIN.1=MIN, MAX.1=MAX)], by='H1')
    dlinkdir <- merge(dlinkdir, drange[,.(H2=AID, MIN.2=MIN, MAX.2=MAX)], by='H2')
    dlinkdir[MAX.1 < MIN.2 & is.na(TYPE), `:=` (TYPE='12', SCORE_DIR=1) ]
    dlinkdir[MAX.2 < MIN.1 & is.na(TYPE), `:=` (TYPE='21', SCORE_DIR=1) ]
    dlinkdir[MAX.1 < MIN.2, stopifnot(all(TYPE=='12'))]
    dlinkdir[MAX.2 < MIN.1, stopifnot(all(TYPE=='21'))]
    dlinkdir[, `:=` (MAX.1=NULL, MAX.2=NULL, MIN.1=NULL, MIN.2=NULL)]

    # Assign source and recipient labels
    dnewpairs <- rbind(
        dlinkdir[TYPE == '12', .(SOURCE=H1, RECIPIENT=H2, SCORE, SCORE_DIR)],
        dlinkdir[TYPE == '21', .(SOURCE=H2, RECIPIENT=H1, SCORE, SCORE_DIR)]
    ) |> unique() 
    stopifnot(dnewpairs[,.N,by=c('SOURCE', 'RECIPIENT')][, all(N==1)])

    dnewpairs[, cat('There are ', .N, 'pairs with score>threshold and',
        sum(SCORE < threshold.likely.connected.pairs),
        'pairs with only one possible direction of tranmsission\n')]

    # table.pairs.in.inland(dnewpairs)



    if( exists("overleaf_expr") )
    {
        # Pairs involving inland participants
        # table.pair.in.inland(dnewpairs)
        double.merge(dnewpairs, dcomms[, .(AID, COMM, SEX)], by_col = 'AID') |>
            subset(COMM.RECIPIENT %like% 'inland' | COMM.SOURCE %like% 'inland') -> tmp_overleaf
        tmp_overleaf[, table(COMM.SOURCE, COMM.RECIPIENT)]
        overleaf_expr[['N_pairs_1_in_study']] <- tmp_overleaf[, .N]

        # now get how many have direction swapped
        rbind(
            serohistory_impact[, .(SOURCE=H1, RECIPIENT=H2, CHANGED_DIR)],
            serohistory_impact[, .(RECIPIENT=H1, SOURCE=H2, CHANGED_DIR)]  
        ) |>
            merge(x = tmp_overleaf, all.x=TRUE, by=c('SOURCE', 'RECIPIENT')) |>
            with(sum(CHANGED_DIR)) -> tmp_overleaf2
        overleaf_expr[['N_swapped']] <- tmp_overleaf2
        
        # number of unique recipients
        overleaf_expr[['N_unique_rec']] <- tmp_overleaf[, uniqueN(RECIPIENT) ]
    }

    # Remove same-sex pairs
    if(! postpone.samesex.removal)
    {
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
    }

    # Solve multiple sources by assigning recipient with strongest linkage support (same for homosexual pairs)
    dnewpairs <- dnewpairs[, {
        z <- which.max(SCORE);
        list(SOURCE=SOURCE[z], SCORE=SCORE[z], SCORE_DIR=SCORE_DIR[z])
    }, by='RECIPIENT']


    if( exists('overleaf_expr') )
    {
        # Pairs involving at both inland participants
        tmp <- table.pairs.in.inland(dnewpairs) |> as.vector()
        names(tmp) <- c('FF', 'MF', 'FM',  'MM')
        overleaf_expr[['N_unique_rec_participants']] <- sum(tmp)
        overleaf_expr[['N_unique_rec_participants_hetero']] <- tmp['MF']+ tmp['FM']
        overleaf_expr[['N_unique_rec_participants_dirs']] <- tmp
    }

    # remove homosexual pairs
    if(! postpone.samesex.removal)
    {
        # solve multiple sources for same-sex too
        dhomosexualpairs <- dnewpairs[ ! RECIPIENT %in% dnewpairs$RECIPIENT, {
            z <- which.max(SCORE);
            list(SOURCE=SOURCE[z], SCORE=SCORE[z], SCORE_DIR=SCORE_DIR[z])
        }, by='RECIPIENT']

    }else{

        cat("Postponing removal of same sex pairs after solving for double sources.\n")
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
    }


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

        tmp1 <- graph_from_data_frame(dnewpairs[, .(SOURCE,RECIPIENT)],
                                      directed=TRUE,
                                      vertices=NULL) |> components()
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

cat('the DOI algorithm is run for a total number of pairs equal to:')
double.merge(chain, meta[, .(AID=aid, SEX=sex)])[, table(SEX.RECIPIENT, SEX.SOURCE)] |>
    knitr::kable() |> print()

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

if(1) # Count number of removed individuals/pairs
{
    tmp <- double.merge(chain, dcomms[, .(AID, COMM, SEX)]) |>
        subset( COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland', select = c('SOURCE','RECIPIENT')) 

    noprobs <- drange_tsi$AID
    overleaf_expr[['N_noTSI']] <- tmp[ ! RECIPIENT %in% noprobs, .N ]
}

# check no contradictions after shrinkning
idx <- double.merge(chain, drange_tsi, by_col = "AID")[TSI.MAX.RECIPIENT - TSI.MIN.SOURCE < 0, unique(IDCLU)]
stopifnot(nrow(idx) == 0)

# 
if(! args$sensitivity.no.refinement)
{
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
        list(AID, RANGE_INTERSECT=rng_intersect, INTERSECT.MIN=r_min, INTERSECT.MAX=r_max)
        }] -> tmp1

    cat('For', tmp1[RANGE_INTERSECT==FALSE, .N], 'individuals the 2 ranges did not intersect\n')
    tmp1[, RANGE_INTERSECT := NULL]


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

}else{

    cat('SENSITIVITY ANALYSIS: Using phyloTSI\n')

    # are there pairs with 

    # keep coherent TSI predictions, and use serohistory if incoherencies
    aid_to_predict <- chain[, unique(c(SOURCE, RECIPIENT))]
    tmp_tsi <- drange_tsi[ AID %in% aid_to_predict]
    tmp_ser <- drange[AID %in% aid_to_predict]
    cat('In', sum(!aid_to_predict %in% tmp$AID), 'cases, incoherent predictions, using serohistory...\n')
    rbind(
          tmp_tsi[, .(AID, MIN=as.Date(TSI.MIN), MAX=as.Date(TSI.MAX))],
          tmp_ser[ ! (AID %in% tmp_tsi$AID) ]
    )  -> drange

    idx <- double.merge(chain[, .(SOURCE,RECIPIENT)], 
                        drange,
                        by_col = "AID")[MAX.RECIPIENT - MIN.SOURCE < 0,]
    cat('In', nrow(idx), 'pairs, the ranges were inconsistent, use serohistory for',
        idx[, uniqueN(c(SOURCE, RECIPIENT))],' participants instead...\n')
    idx <- idx[, c(SOURCE,RECIPIENT)]

    drange <- rbind( drange[! AID %in% idx], tmp_ser[AID %in% idx])
}


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

# RUN MCMC
# ________

range_gi <- c(.5, 1:10)
tmp <- dinfectiousness[!is.infinite(END), paste(round(END, 2), sep='')]
tmp <- gsub('\\.', '', tmp)

filename_drange <- 'networks_individualDOIrange'
filename_net <- 'networks_GICentroids'
filename_overleaf <- 'pairsinfo_overleaf'
suffix <- ''
if( ! is.null(dinfectiousness))
    suffix <- paste0('_d', tmp, '_w', paste0( dinfectiousness$VALUE, collapse=''))
if( build.network.from.pairs  )
    suffix <- paste0(suffix, '_netfrompairs')
if( threshold.likely.connected.pairs != .5 )
    suffix <- paste0(suffix, '_thr', gsub( '0\\.', '', threshold.likely.connected.pairs))
if( threshold.direction != .5 )
    suffix <- paste0(suffix, '_dir', gsub( '0\\.', '', threshold.direction))
if( get.sero.extra.pairs )
    suffix <- paste0(suffix, '_seropairs')
if( args$sensitivity.no.refinement)
    suffix <- paste0(suffix, '_sensnoref')
if( postpone.samesex.removal )
    suffix <- paste0(suffix, '_postponessrem')
cat('\n Chosen suffix: ', suffix, '\n')

filename_drange <- file.path(out.dir, paste0(filename_drange, suffix,  '.rds'))
filename_net <- file.path(out.dir, paste0(filename_net, suffix, '.rds'))
filename_overleaf <- file.path(out.dir, paste0(filename_overleaf, suffix, '.rds'))

if( ! is.null(dinfectiousness) & build.network.from.pairs & get.sero.extra.pairs & ! args$sensitivity.no.refinement)
{
    out <- list(drange, chain)
    cat('Saving', filename_drange, '...\n')
    saveRDS(out, filename_drange)
}

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

if( file.exists(filename_net) & ! args$rerun )
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

# Final results
filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_toshare', suffix, '.rds')) 

if(  file.exists(filename) & ! args$rerun == TRUE )
{
    dresults <- readRDS(filename)
}else{

    # prepare input for Bayesian model
    centroids <- dcohords[, CENTROIDS[[1]], by='GROUP']
    dresults <- prepare.pairs.input.for.bayesian.model(centroids)
    dresults[, table(DIRECTION, useNA='ifany')]

    dhomosexualpairs[, DIRECTION := 'phyloscanner']
    dhomosexualpairs[, `:=` (SCORE=NULL, SCORE_DIR=NULL, SEX.SOURCE=NULL, SEX.RECIPIENT=NULL)]
    dhomosexualpairs <-  double.merge(dhomosexualpairs,meta[, .(AID=aid, SEX=sex)])
    rbind(
        dresults,
        dhomosexualpairs,
        fill=TRUE
    ) -> dresults
    dresults <- get.community.type.at.infection.date(dresults)

    # denote which directions were fixed through serohistory
    dresults <- rbind(
        serohistory_impact[, .(SOURCE=H1, RECIPIENT=H2, CHANGED_DIR)],
        serohistory_impact[, .(RECIPIENT=H1, SOURCE=H2, CHANGED_DIR)]  
    ) |> merge(x = dresults, all.x=TRUE)
    dresults[ CHANGED_DIR == TRUE, DIRECTION := 'sero-adjusted']
    dresults[, CHANGED_DIR := NULL]

    if(get.sero.extra.pairs & ! build.network.from.pairs) 
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
    if(! args$sensitivity.no.refinement)
        stopifnot(dresults[, all(CU < DATE.COLLECTION.RECIPIENT, na.rm=TRUE)])

    if(exists('overleaf_expr'))
    {
        overleaf_expr
        idx_inland <- dresults[ COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland']
        idx_inland_hetero <- idx_inland[ SEX.SOURCE != SEX.RECIPIENT ]
        overleaf_expr[['N_final_inland']] <- idx_inland_hetero[, .N]

        overleaf_expr[['N_final_inland_outstudy']] <- idx_inland_hetero[ is.na(ROUND.M), .N ]
        overleaf_expr[['N_final']] <- idx_inland_hetero[ !is.na(ROUND.M), .N ]
    }

    dresults[!is.na(M) &
             SEX.SOURCE != SEX.RECIPIENT &
             COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland', {
                 cat('There are', .N, 'pairs with median DOI estimate(', sum(is.na(ROUND.M)), 'outside of range).\n'); .SD
             } ][is.na(ROUND.M), .(SOURCE, RECIPIENT, M)] -> idx
    idx <- double.merge(idx, meta[, .(AID=aid, DFP=date_first_positive)])
    setkey(dprobs_roundallocation, SOURCE,RECIPIENT)
    setkey(idx, SOURCE,RECIPIENT)

    cols <- c('SOURCE', 'RECIPIENT', 'M')
    saveRDS(dresults, filename)    
    dresults_new <- readRDS(filename)

    saveRDS(overleaf_expr, filename_overleaf)
}

# make table
cols <- c('SOURCE', 'RECIPIENT',
          'SEX.SOURCE', 'SEX.RECIPIENT',
          'M', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT',
          'DIRECTION')

double.merge( dresults[, -c('COMM.SOURCE', 'COMM.RECIPIENT')], dcomms[, .(AID, COMM)]) |>
    subset( COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland',
        select=c('SOURCE', 'RECIPIENT')
    ) -> idx

dtable <- merge(idx, dresults, all.x=TRUE, by=c('SOURCE', 'RECIPIENT'))
dtable <- subset(dtable,select=cols)
setnames(dtable, 'M', 'INFECTION_DATE')
filename <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', paste0('pairsdata_table_toshare', suffix, '.csv'))
fwrite(dtable, filename)
