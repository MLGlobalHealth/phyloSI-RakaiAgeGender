################
# DEPENDENCIES #
################ 
library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)
library(ggpubr)

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
# file.path.chains.data <- .fp('X','211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks_ruleo_sero.rda')
file.anonymisation.keys <- .fp('X','important_anonymisation_keys_210119.csv')
file.path.tsiestimates <- .fp('A', 'PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_001_rla_T_zla_T/aggregated_TSI_with_estimated_dates.csv')
file.path.round.timeline <- .fp('D', 'RCCS_data_estimate_incidence_inland_R6_R18/220903/RCCS_round_timeline_220905.RData')
path.network <-file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'pairsdata_toshare_d1_w11_netfrompairs_postponessrem.rds')
path.range <-file.path(out.dir, 'networks_individualDOIrange_d1_w11_netfrompairs_seropairs.rds')

file.exists(file.path.meta,
            # file.path.chains.data,
            file.anonymisation.keys,
            file.path.tsiestimates,
            file.path.round.timeline 
            ) |> all() |> stopifnot()

threshold.likely.connected.pairs <- 0.5
################
#    HELPERS   #
################

source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'confidential_data_src/utils', 'functions_tsi_attribution.R'))
source(file.path(indir, 'functions', 'plotting_functions.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
# source(file.path(indir, 'functions', 'statistics_functions.R'))
# source(file.path(indir, 'functions', 'stan_utils.R'))
# source(file.path(indir, 'functions', 'check_potential_TNet.R'))
find_palette_round()
naturemed_reqs()

##############################
#   EXTENDED DATA FIGURES   # 
##############################

#
# Load data
#

dresults <- readRDS(path.network)    
tmp <- readRDS(path.range)
drange <- tmp[[1]]
chain <- tmp[[2]]
aik <- fread(file.anonymisation.keys, 
             header = TRUE, select=c('PT_ID', 'AID'))

# get meta 
meta_env <- new.env()
load(file.path.meta, envir=meta_env)
meta <- subset(meta_env$meta_data,
               select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
meta <- unique(meta[!is.na(aid)])
stopifnot(meta[, uniqueN(aid) == .N])

dcomms <- get.communities.where.participated()

# For plots, only select pairs which appear in the analysis

if ("SEX.RECIPIENT.x" %in% names(dresults)) {

    dresults[, SEX.RECIPIENT := fcoalesce(SEX.RECIPIENT, SEX.RECIPIENT.x, SEX.RECIPIENT.y )]
    dresults[, SEX.SOURCE := fcoalesce(SEX.SOURCE, SEX.SOURCE.x, SEX.SOURCE.y )]

    dresults[, `:=` (SEX.RECIPIENT.x = NULL, SEX.RECIPIENT.y = NULL, SEX.SOURCE.x = NULL, SEX.SOURCE.y = NULL  )]
}

dresults[, {cat(.N); table(COMM.SOURCE, COMM.RECIPIENT)} ]

only.inland.participants <- 1
if( only.inland.participants )
{
    idx_ever_inland <- double.merge( dresults[, .(SOURCE,RECIPIENT)], dcomms[, .(AID, COMM)] ) |>
        subset(COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland') |>
        subset(select=c('SOURCE', 'RECIPIENT'))

    dresults <- dresults[idx_ever_inland]
}else{
    cat("Not subsetting exclusively to participants.\n")
}

# subset to heterosexual pairs only 
dresults <- dresults[SEX.SOURCE != SEX.RECIPIENT ]
dresults[, .(.N, sum(!is.na(M)), sum(!is.na(ROUND.M))) ]

cat(dresults[, .N], 'source-recipient pairs selected\n')

if(1)
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

if(1)
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

    dresults_tsi <- prepare.pairs.input.for.bayesian.model(tsi_predictions)

    plot.age.comparison <- function(part)
    {
        stopifnot(part %in% c('SOURCE', 'RECIPIENT'))

        part_age_col <- fifelse(part=='SOURCE', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT')
        stopifnot(part_age_col %in% names(dresults) & part_age_col %in% names(dresults_tsi))

        cols <- c(part, part_age_col)
        .f <- function(DT) subset(DT, select=cols)
        dage_comparison <- list(.f(dresults_tsi), .f(dresults))

        stopifnot(part %in% names(dage_comparison[[1]]) & part %in% names(dage_comparison[[2]]))


        if(part=='RECIPIENT')
        {
            dage_comparison <- merge(dage_comparison[[1]], dage_comparison[[2]], all.y=TRUE, by=part)
            dage_comparison <- melt(dage_comparison, id.vars='RECIPIENT',
                                    value.name='AGE_INFECTION.RECIPIENT',
                                    variable.name='METHOD')
            dage_comparison[METHOD %like% '.x$', METHOD := 'phyloTSI' ]
            dage_comparison[METHOD %like% '.y$', METHOD := 'final' ]
        }else{

            dage_comparison[[1]][, METHOD := 'phyloTSI']
            dage_comparison[[2]][, METHOD := 'final']
            dage_comparison <- Reduce(rbind,dage_comparison)
        }

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
    force(plot_age_comparison_source)
    force(plot_age_comparison_recipient)
}


# plot seroconverters
if(0)
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
}

# Plot schema
if(0)
{   # BETTER TO GO ON WITH FAKE DATA
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
filename <- file.path(out.dir, 'edf_tsis_v_v2_230116.pdf')
cat('Saving', filename, '...\n')
ggsave_nature(filename, tmp, add_reqs=FALSE)

if(0)   # study generation intervals
{       
    dgen <- dcohords[, GENINTS[[1]], by='GROUP']
    dgen <- dgen[, mean(V1), by=GI]
    dgen[GI %in% c(0.5, 1), {
            cat('We estimate', round(100*V1, 2), '% generation intervals to be less than', GI, 'years\n')
            V1
    }, by='GI']
}



# sensitivity analysis
generation.interval.sensitiviy.analysis <- function()
{

    files <- list.files(file.path(indir.deepsequencedata, 'RCCS_R15_R18'), 
                        pattern='pairsdata.*rds$', full.names = TRUE)
    files <- grep('netfrompairs_seropairs', files, value=TRUE)

    lresults <- lapply(files, readRDS)
    lapply(lresults, function(DT)
           {    # Subset Results to source-recipient pairs in the analysis
               DT <- DT[SEX.SOURCE!=SEX.RECIPIENT]
               DT <- DT[ COMM.SOURCE == 'inland' & COMM.RECIPIENT == 'inland'] 
               # DT <- DT[COMM.RECIPIENT == 'inland', ]
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

    dresults[!is.na(M) & COMM.SOURCE=='inland'& COMM.SOURCE=='inland', .N]

    ggplot(dsens) +
        geom_errorbar(aes(xmin=CL.uniform, xmax=IU.uniform, y=M.bellan)) + 
        geom_errorbar(aes(ymin=CL.bellan, ymax=IU.bellan, x=M.uniform)) + 
        #geom_abline() +
        theme_bw() +
        labs()
}

cat('End of script\n')
