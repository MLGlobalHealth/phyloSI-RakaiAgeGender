################
# DEPENDENCIES #
################ 

library(data.table)
library(ggplot2)
library(lubridate)
library(xtable)
library(ggpubr)
library(here)

################
#    PATHS     #
################

gitdir <- here::here()
source(file.path(gitdir, 'config.R'))

if(usr == 'andrea')
{
    out.dir <- '/home/andrea/HPC/ab1820/home/projects/2022/genintervals'
}

# file.anonymisation.keys <- file.path(indir.deepanalyses.xiaoyue,'important_anonymisation_keys_210119.csv')

# file.pairs |> readRDS() |> nrow()
# file.pairs.nonrefined |> readRDS() |> nrow()

file.exists(path.meta.confidential,
            file.anonymisation.keys,
            path.tsiestimates,
            file.path.round.timeline 
            ) |> all() |> stopifnot()

threshold.likely.connected.pairs <- 0.5

################
#    HELPERS   #
################

source(file.path(gitdir.R, 'utils.R'))
source(file.path(gitdir.R.conf, 'functions_tsi_attribution.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R'))
source(file.path(gitdir.R.flow, 'summary_functions.R'))
find_palette_round()
naturemed_reqs()

##############################
#   EXTENDED DATA FIGURES   # 
##############################

#
# Load data
#

only.inland.participants <- TRUE

df_round_gi <- get.round.dates(file.path.round.timeline)

# get pairs and remove those for which we do not have predictions 
dresults <- readRDS(file.pairs) |> 
    subset(!is.na(M))
# aik <- fread(file.anonymisation.keys, header = TRUE, select=c('PT_ID', 'AID'))

# get meta 
meta_env <- new.env()

load(path.meta.confidential, envir=meta_env)

meta <- subset(meta_env$meta_data,
               select=c('aid', 'sex', 'date_birth', 'date_first_positive', 'date_last_negative'))
meta <- unique(meta[!is.na(aid)])
stopifnot(meta[, uniqueN(aid) == .N])

dcomms <- get.communities.where.participated(path.meta.confidential)

if( only.inland.participants )
{
    cat("only select inland pairs\n") 
    idx_ever_inland <- double.merge( dresults[, .(SOURCE,RECIPIENT)], dcomms[, .(AID, COMM)] ) |>
        subset(COMM.SOURCE %like% 'inland' & COMM.RECIPIENT %like% 'inland') |>
        subset(select=c('SOURCE', 'RECIPIENT'))

    dresults <- dresults[idx_ever_inland]
}else{
    cat("Not subsetting exclusively to participants.\n")
}
dresults[, {cat(.N); table(COMM.SOURCE, COMM.RECIPIENT)} ]
cat(dresults[, .N], 'source-recipient pairs selected\n')

if(1)
{
    # removed code 4 beautiful p_predictions
    rec_preds <- unique(dresults[, .(RECIPIENT, CL, M, CU)])

    # get phyloTSI estimates (+range) for recipients 
    drange_tsi2 <- get.infection.range.from.tsi(path.tsiestimates, path.sequence.dates=path.collection.dates.confidential, exclude_mid=FALSE, chain_subset=FALSE) |> 
        setnames('AID', 'RECIPIENT') |> 
        set(i=NULL, 'RENAME_ID', NULL)
    # transform iDate to Date...
    cols <- c('MIN', 'MID', 'MAX')
    drange_tsi2[, (cols) := lapply(.SD, as.Date) , .SDcols=cols]
    drange_tsi2 <- drange_tsi2[RECIPIENT %in% rec_preds$RECIPIENT]

    rec_preds <- merge(rec_preds, drange_tsi2, by='RECIPIENT', all.x=TRUE)
    rec_preds[, INTERSECT := pmax(CL, MIN) <= pmin(CU, MAX)]
    rec_preds <- rec_preds[!is.na(MIN)]

    # Figure 5A 
    p_predictions2 <- ggplot(rec_preds) + 
        geom_rect(data=df_round_gi, aes(xmin=MIN_SAMPLE_DATE, xmax=MAX_SAMPLE_DATE, 
                                     ymin=MIN_SAMPLE_DATE, ymax=MAX_SAMPLE_DATE, fill=as.ordered(round))) + 
        geom_rect(data=df_round_gi, aes(xmin=as.Date(-Inf), xmax=MAX_SAMPLE_DATE, ymin=MIN_SAMPLE_DATE, ymax=MAX_SAMPLE_DATE, 
                                     fill=as.ordered(round)), alpha=.1) + 
        geom_rect(data=df_round_gi, aes(ymin=as.Date(-Inf), ymax=MAX_SAMPLE_DATE, 
                                    xmin=MIN_SAMPLE_DATE, xmax=MAX_SAMPLE_DATE, 
                                    fill=as.ordered(round)), alpha=.1) + 
        geom_abline(aes(slope=1, intercept=0), linetype='dashed', color='red' ) + 
        geom_point(aes(x=MID, y=M), size=.4) + 
        geom_errorbar(aes(x=MID ,ymin=CL, ymax=CU), alpha=.1) + 
        geom_errorbarh(aes(y=M ,xmin=MIN, xmax=MAX), alpha=.1) + 
        scale_fill_manual(values=palette_round) +
        guides(fill=guide_legend(ncol=1, override.aes=list(size=1))) + 
        theme_bw() + 
        theme(legend.position='right', legend.direction='vertical', 
              panel.border = element_blank(), 
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black")) +
        labs(x="infection time estimates using phyloTSI on deep-sequence data",
             y="refined infection time estimates accounting further for sero-history and transmission direction",
             fill='Round', alpha='') 
    force(p_predictions2)

    # get MAE for people with last negative test
    # Median absolute error for 
    dsero <- meta |> 
        subset(!is.na(date_first_positive) & !is.na(date_last_negative), 
            select=c('aid', 'date_last_negative', 'date_first_positive')) |>
            setnames('aid', 'AID')
    stopifnot(dsero[, uniqueN(AID) == .N ])
    dsero <- dsero[, .(midpoint=mean(c(date_last_negative, date_first_positive))), by='AID']

    dsero_rec <- merge(rec_preds, dsero, by.x='RECIPIENT', by.y='AID') |> 
        subset(!is.na(M))
    dsero_rec[, final := as.numeric(abs(M - midpoint)/365.25)]
    dsero_rec[, phyloTSI := as.numeric(abs(MID - midpoint)/365.25)]

    tmp2 <- melt( dsero_rec[, .(RECIPIENT, final, phyloTSI)],
                 id.vars='RECIPIENT', variable.name='METHOD', value.name='AE') 
    means <- tmp2[, .(MEAN=mean(AE)), by='METHOD']

    .rm <- function(x){
        fifelse(x %like% 'final|Final',
                yes='refined infection time estimates accounting\nfurther for serohistory and transmission direction',
                no='phyloTSI on\ndeep sequence data')
    }

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
    tsi_predictions <- copy(rec_preds)
    tsi_predictions <- tsi_predictions[ , .(ID=RECIPIENT, CL=MIN, M=MID, CU=MAX)]
    dresults_tsi <- prepare.pairs.input.for.bayesian.model(tsi_predictions, CHAIN=dresults)

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

p_edf_doi <- ggarrange_nature(plot_mae + theme(legend.spacing.x=unit(.3, 'cm')),
                        plot_age_comparison_source, plot_age_comparison_recipient,
                        common.legend = TRUE, legend = 'bottom',
                        labels = c('b','c','d'), ncol=3)
p_edf_doi <- ggarrange_nature(p_predictions2, p_edf_doi, ncol=1, heights = c(6,4), labels=c('a', '')) 
filename <- file.path(out.dir, 'edf_tsis_v_v2_230116.pdf')
cat('Saving', filename, '...\n')
ggsave_nature(filename, p_edf_doi, add_reqs=FALSE)

cat('End of script\n')

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
