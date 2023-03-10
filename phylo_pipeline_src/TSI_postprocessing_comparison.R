# We need to compare the predictions obtained with Tanya's run and with ours
# As such, I want to be able to produce two plots systematically:
# - the first one comparing the predictors and predictions between methods
# - the second one testing the predictions on the known seroconversion intervals

# I need to load the aggregated predictions for our analysis
# And Tanyas


# Packages
#_________

require(ggplot2)
require(data.table)
require(gridExtra)


# Args
#_____
cat("arguments...\n")

option_list <- list(
  optparse::make_option(
    "--relationship_dir",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to directory containing analyse_trees.R results", 
    dest = 'rel.dir'
  ),
  optparse::make_option(
    "--TSI_dir",
    type = "character",
    default = '~/git/HIV-phyloTSI-main/',
    help = "Absolute file path to HIV-phylo-TSI-main repository", 
    dest = 'TSI.dir'
  ),
  optparse::make_option(
    "--input_samples",
    type = "character",
    default = NA_character_,
    help = "Absolute file path to input samples rds containing PANGEA_IDs and RENAME_IDs", 
    dest = 'phsc.samples'
  )
)

args <- optparse::parse_args(optparse::OptionParser(option_list = option_list))
args[['help']] <- NULL

################
# Testing
################

user <- Sys.info()[['user']]
if(user == 'andrea')
{
    out.dir <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/'
    args <- list(
                 out.dir=out.dir,
                 rel.dir= file.path(out.dir, "2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_00"), 
                 phsc.samples=file.path(out.dir, "220331_RCCSUVRI_phscinput_samples_with_bf.rds"),
                 date='2022-08-22'
    )
}

user <- Sys.info()[['user']]
if(0){
  args <- list(
    out.dir = "/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/",
    pkg.dir = "/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software",
    rel.dir="/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/2022_08_22_phsc_phscTSI_sd_42_sdt_002_005_dsl_100_mr_30_mlt_T_npb_T_og_REF_BFR83HXB2_LAI_IIIB_BRU_K03455_phcb_T_rtt_00",
    phsc.samples="/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS_MRC_UVRI_TSI/220331_RCCSUVRI_phscinput_samples_with_bf.rds",
    TSI.dir="/rds/general/user/ab1820/home/git/HIV-phyloTSI",
    controller='/rds/general/user/ab1820/home/git/Phyloscanner.R.utilities/misc_data_analysis_RCCS1519/software/runall_TSI_pairs2.sh'
  )
        if(user == 'andrea')
        {
                .f <- function(x) 
                {
                        x <- gsub('/rds/general/project','/home/andrea/Documents/Box',x)
                        x <- gsub('/rds/general/user/ab1820/home','/home/andrea',x)
                        x <- gsub('HIV-phyloTSI','HIV-phyloTSI-main',x)
                }
                args <- lapply(args, .f)
        }
}

print(args)
tmp <- file.exists(unlist(args))
stopifnot( all(tmp) )


# Paths
#________ 

usr <- Sys.info()[['user']]
if (usr == 'andrea')
{
        # If local
        indir.deepsequence.analyses <- '~/Documents/Box/ratmann_deepseq_analyses/live'
        indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata'
        tanya.rakai.dir <- '~/git/HIV-phyloTSI-main/RakExample_Tanya'
        indir.deepsequence.analyses.xiaoyue <-  indir.deepsequence.analyses

        # args$rel.dir <-  file.path(indir.deepsequence.analyses, 'seroconverters2/phscrel')
        # tmp <- dirname(args$rel.dir)
        # args$phsc.samples <- list.files(tmp, pattern='phscinput_samples', full.names=T)
}else{
        # if HPC
        indir.deepsequence.analyses <- '/rds/general/project/ratmann_deepseq_analyses/live'
        indir.deepsequence.analyses.xiaoyue <- "/rds/general/project/ratmann_xiaoyue_jrssc2022_analyses/live"
        indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live'
        tanya.rakai.dir <- file.path(args$TSI.dir, 'RakExample_Tanya')
}

# Analysis results:
dtsi.all.path <- list.files(args$rel.dir, pattern='tsi', full.names=TRUE)
dtsi.aggregated.path <- file.path(args$rel.dir, "aggregated_TSI_with_estimated_dates.csv")

# Meta data
indir.deepsequence.analyses.old <- file.path(indir.deepsequence.analyses.xiaoyue, 'PANGEA2_RCCS1519_UVRI')
file.db.sharing <- file.path(indir.deepsequencedata,"/PANGEA2_RCCS/200316_pangea_db_sharing_extract_rakai.csv")
file.anonymisation.keys <- file.path(indir.deepsequence.analyses.old,'important_anonymisation_keys_210119.csv')
file.seroconverters <- file.path(indir.deepsequencedata, 'PANGEA2_RCCS', '220329_TSI_seroconverters.csv')

tmp <-file.exists(indir.deepsequence.analyses.old, file.db.sharing,
                  file.anonymisation.keys, file.seroconverters)

print(tmp)
stopifnot(all(tmp))




# Microproccesing + helper functions
#___________________________________

predictors.hiv.phylo.tsi <-  c("gag_lrtt", "pol_lrtt", "gp120_lrtt", "gag_maf3c", "gp41_maf3c", "gp41_maf12c", "gag_tips", "gp41_tips", "gp120_tips")

aik <- fread(file.anonymisation.keys, header=TRUE)
aik[, V1:=NULL]

.aid2studyid <- function(x, rm_prefix=FALSE, with_fq=FALSE)
{
    x <- data.table(RENAME_ID = x)
    
    cols <- c('RENAME_ID', 'UNIT_ID')
    tmp <- unique(dsamples[, ..cols])
    if(! with_fq)
            tmp[, RENAME_ID:=gsub('-fq[0-9]$', '', RENAME_ID)]

    x <- merge(x, tmp, by='RENAME_ID', all.x=TRUE)$UNIT_ID
    if(rm_prefix)
            x <- gsub('^.*?-','',x)
    x
}

.sid2pid.and.date <- function(x)
{
    y1 <- gsub('^(.*?)-(.*?)$', '\\1', x)
    y2 <- gsub('^(.*?)-(.*?)$', '\\2', x)

    y2 <- as.numeric(y2)
    year <- y2 %/% 1 
    date <- as.Date(paste0('01-01-', year), '%d-%m-%Y')
    tmp1 <- round(y2 %% 1 *365)
    date <- date + tmp1

    list(PID=y1,date=date)
}

preprocess.tanya <- function()
{
    # Tanya's analysis
    maf_tanya <- fread(file.path(tanya.rakai.dir, 'rak_maf.csv'), header=TRUE)
    ps_tanya <- fread(file.path(tanya.rakai.dir, 'rak_ps.csv'), header=TRUE)
    out_tanya <- fread(file.path(tanya.rakai.dir, 'rak_out.csv'), header=TRUE)

    # extract PID and sample date from host_id
    newcols <- c('PT_ID', 'sample_date')
    maf_tanya[, (newcols):=.sid2pid.and.date(sid)]
    ps_tanya[, (newcols):=.sid2pid.and.date(sid)]
    out_tanya[, (newcols):=.sid2pid.and.date(host.id)]
    out_tanya
    
    # extract predictors and predictions
    tmp <- grep('sqrt|cc', colnames(out_tanya), value=TRUE)
    cols <- c("PT_ID", "sample_date", predictors.hiv.phylo.tsi, tmp)
    predictors.tanya <- out_tanya[, .SD, .SDcols=cols]
    predictors.tanya
}

preprocess.ours <- function()
{
    # Load all outputs from HIVphyloTSI
    # dtsi <- fread(dtsi.aggregated.path)
    # dtsi[, V1 := NULL]
    tmp <- lapply(dtsi.all.path, fread)
    dtsi.all <- rbindlist(tmp, use.names=TRUE)
    dtsi.all[, AID:=gsub('-fq[0-9]','',host.id)]
    dtsi.all <- merge(aik, dtsi.all, by='AID')

    if(0)
    {
            # It seems like lrtt vary for different runs of the same ID
            # On the other hand the mafs dont. This makes sense
            dtsi.all.counts <- dtsi.all[, .(N=.N, 
                                            Nmaf=uniqueN(genome_maf3c),
                                            Nlrtt=uniqueN(genome_lrtt)), by='host.id']
    }

    # extract predictors and predictions
    tmp <- grep('sqrt|cc', colnames(dtsi.all), value=TRUE)
    cols <- c("AID", "PT_ID", 'host.id', predictors.hiv.phylo.tsi, tmp)
    predictors.ours <- dtsi.all[, .SD, .SDcols=cols]
    predictors.ours <- predictors.ours[grepl('^RK',PT_ID)]
    predictors.ours[, PT_ID := gsub('RK-', '',PT_ID)] 
    predictors.ours
}

compare.predictors <- function(predictors.ours, predictors.tanya)
{
    # merge data and sort by differences in preds
    cols <- colnames(predictors.ours)
    cols <- cols[cols %in% c(colnames(predictors.tanya), 'visit_dt', 'sample_date')]
    
    # If sample dates are included in our predictor datasets, rename
    if( any(c('visit_dt','sample_date') %in% cols) )
    {
      setnames(predictors.ours, 'visit_dt', 'sample_date', skip_absent=TRUE)
      cols <- gsub( 'visit_dt', 'sample_date', cols) 
      predictors.merged <- merge(predictors.ours[,..cols],
                                 predictors.tanya[,..cols],
                                 by=c('PT_ID', 'sample_date'))
      # predictors.ours[, .N, by=c('PT_ID', 'sample_date')]
      
    }else{
      predictors.merged <- merge(predictors.ours[,..cols], predictors.tanya[,..cols], by='PT_ID')
    }
    
    predictors.merged[, diff := (RF_pred_sqrt.x - RF_pred_sqrt.y)]
    
    # plot comparisons
    .plot.pred.comparison <- function(tx)
    {
            ty <- gsub('x$', 'y', tx)
            t <- c(tx,ty, 'diff')

            t %in% names(predictors.merged)
            
            dplot <- predictors.merged[, .SD, .SDcols=t]
            colnames(dplot) <- c('x','y','c')
            gg <- ggplot(dplot, aes(x=x, y=y, color=c)) +
                    geom_abline(slope=1, color='black') + 
                    geom_smooth(method='lm', formula= y~x, color="grey50", linetype='dotted', se=FALSE)+ 
                    geom_point() +
                    guides(colour='none') +
                    # scale_colour_gradient(low="blue", high="red", midpoint=0)
                    scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 0) +
                    labs(x="our", y="tanya's",                             title=gsub('.x$','',tx)) +
                    theme_bw() +
                    expand_limits(x = 0, y = 0) +
                    theme(plot.title = element_text(hjust = 0.5))
            gg
    }
    cols <- colnames(predictors.merged)
    tmp <- grep('\\.x$',cols, value=T)
    tmp
    gg_list <- lapply(tmp, .plot.pred.comparison)


    if(1)
    {
            do.call("grid.arrange", c(gg_list, ncol=3)) -> gg1
            filename=file.path(args$rel.dir , "inputs_outputs_comparison.png")
            ggsave(filename , gg1, w=10, h=14)
    }

    
    .plot.compare.95ci <- function(predictors.merged)
    {
            cols <- grep('RF',colnames(predictors.merged), value=T) 
            tmp <- predictors.merged[, .SD, .SDcols=cols]
            gg <- ggplot(tmp, aes(x=RF_pred_sqrt.x, y=RF_pred_sqrt.y)) + 
                    geom_abline(slope=1, linetype='dotted', color='red') + 
                    geom_point() + 
                    geom_crossbar(aes(ymin=RF_cc025.y, ymax=RF_cc975.y)) + 
                    geom_linerange(aes(xmin=RF_cc025.x, xmax=RF_cc975.x))
            gg
    }
    gg2 <- .plot.compare.95ci(predictors.merged)
    
    return(list(gg1=gg1, gg2=gg2, DT=predictors.merged))
}

plot.timeline.all.seroconverters <- function(DTSI, DSERO, key_value='visit_dt', show_preduncertainty=FALSE, shift_visit_dt=FALSE)
{
    # check uncertainty
    cols <- names(DTSI)[names(DTSI) %like% 'study_id|doi|visit_dt']
    tmp <- subset(DTSI, select=cols)

    tmp <- merge(DSERO, tmp, by=c('study_id', 'visit_dt'))
    cat(nrow(tmp), 'seroconverters selected in `plot.timeline.all.seroconverters`\n')
    tmp[, date_fpos := visit_dt - 365.25 * knownTSI_min ]
    tmp[, date_lneg := visit_dt - 365.25 * knownTSI_max ]
    tmp[, `:=` (knownTSI_min=NULL, knownTSI_max=NULL )]

    # key_value = c('pred_doi_mid', 'midpoint')
    if( all(key_value %in% names(tmp)) )
    {
        if(length(key_value) == 1)
            tmp[, KEY := .SD, .SDcols=key_value ]
        if(length(key_value) == 2)
        {
            tmp1 <- tmp[, .SD, .SDcols=key_value]
            tmp[, KEY:=(tmp1[[1]] - tmp1[[2]])]
        }
    }else{
        stop('provide a key_value')
    }
    setkey(tmp, KEY)
    tmp[, study_id := ordered(study_id, levels=tmp$study_id)]

    # get intersections
    dintersect <- tmp[ ! (pred_doi_max < date_lneg | pred_doi_min > date_fpos) ]
    # dintersect[, .N] / nrow(tmp)
    coh <- tmp[,  mean(date_lneg <= pred_doi_mid & pred_doi_mid <= date_fpos) ] 

    is.idate <- function(x) inherits(x, 'IDate')
    cols <- tmp[, sapply(.SD, is.idate) , ]
    cols <- names(which(cols))
    tmp[, (cols) := lapply(.SD, as.Date) , .SDcols=cols]

    cols <- names(tmp)[names(tmp) %like% 'doi|midpoint|date|dt']
    lims <- c()
    if(shift_visit_dt)
    {
        tmp[, (cols) := lapply(.SD, function(x) x-visit_dt) , .SDcols=cols]
        lims <- tmp[KEY > - 365.25, study_id[1]]
        lims <- c(lims, tmp[KEY > 365.25, study_id[1]])
    }
    

    gg <- ggplot(tmp, aes(y=study_id, yend=study_id))

    if(show_preduncertainty)
        gg <- gg + geom_segment( aes(x=pred_doi_min, xend=pred_doi_max), alpha=.5, color='blue' ) 
    for(l in lims)
        gg <- gg +
            geom_hline(aes(yintercept=l), color='red', linetype='dotted') 

    gg <- gg +
        geom_segment( aes(x=date_lneg, xend=date_fpos), alpha=.5, color='green' ) + 
        geom_point( aes(x=visit_dt)) + 
        geom_point( aes(x=pred_doi_mid), color='blue', pch=23) + 
        geom_point( aes(x=midpoint), color='green', pch=23) + 
        labs(x='dates', y='study_participants') + 
        theme_bw() + 
        theme(legend.position='bottom',
              axis.text.y=element_blank(),
              axis.ticks.x=element_blank())

    .year.diff <- function(x) as.numeric(x/365.25)
    .mae <- function(x,y)
        mean(abs(.year.diff(x-y)))
    .mse <- function(x,y)
        mean(.year.diff(x-y)^2)

    mae <- tmp[, .mae(midpoint, pred_doi_mid)]
    mse <- tmp[, .mse(midpoint, pred_doi_mid)]
    list(plot=gg, mae=mae, mse=mse, coh=coh)
}

plot.sero.predictions <- function(data, title='')
{       
    # data <- copy(dtsi)
    # uniqueN(data) == nrow(data)
    data[, in_sero_interval := (RF_pred_linear < knownTSI_max) & (RF_pred_linear > knownTSI_min)]
    data <- data[!is.na(knownTSI) & !is.na(RF_pred_linear)]

    # get lab: correlation + predictions falling in seroconversion interval
    lab <- data[knownTSI < 5, cor(knownTSI, RF_pred_linear)]
    lab <- round(lab, 2)
    lab <- paste0('correlation = ', lab, '\n')
    tmp <- data[, round(mean(in_sero_interval), 2)]
    lab <- paste0(lab, 'coherent predictions = ', tmp)
    y <- data[, 0.9* max(knownTSI-RF_pred_linear)]

    cat(data[,.N], ' datapoints plotted \n' )
    gg <- ggplot(data, aes(x=knownTSI, y=knownTSI - RF_pred_linear)) + 
            geom_point(aes(color=in_sero_interval)) +
            geom_hline(aes(yintercept=0), color='red', linetype='dotted') + 
            theme_bw() +
            theme(legend.position='bottom') +
            geom_label(aes(x=3, y=y, label=lab)) + 
            labs(x="Known TSI (sample collection date - midpoint)", y="Known TSI - estimated TSI", color='Prediction in seroconversion interval') +
            ggtitle(title)

    return(gg)
    data[knownTSI - RF_pred_linear <= 1, mean(in_sero_interval)]
    data[knownTSI - RF_pred_linear <= 1 & ! in_sero_interval, ]

}

show.corrs <- function(data, title='')
{      
    data <- data[!is.na(knownTSI) & !is.na(RF_pred_linear)]
    lab <- round( cor(data$knownTSI, data$RF_pred_linear), 2 )
    lab <- paste0('corr =', lab)
    tmp <- data[, 0.9* max(RF_pred_linear)]

    cat(data[,.N], ' datapoints plotted \n' )
    ggplot(data, aes(x=knownTSI, y= RF_pred_linear)) + 
            geom_point() + 
            geom_hline(aes(yintercept=0), color='red', linetype='dotted') + 
            geom_label(aes(x=2.5, y=tmp, label=lab)) + 
            theme_bw() + 
            theme(legend.position='bottom') +
            labs(x="Known TSI" , y="Estimated TSI") +
            ggtitle(title)
}

.extract.repeated.only <- function(DT, cols)
{
    cols <- c('study_id')
    tmp <- DT[, .N, by=cols]
    tmp <- tmp[ N > 1, ..cols]
    merge(tmp, DT, by=cols)
}

########
# Main #
########

# get samples
dsamples <- readRDS(args$phsc.samples)
dsamples[, pangea_id:= gsub("^.*?_", "",PANGEA_ID)]

# get TSI estimates
dtsi <- fread(dtsi.aggregated.path, drop='V1')
dtsi[, study_id:= .aid2studyid(RENAME_ID, with_fq=TRUE),]

if('visit_dt' %in% names(dtsi))
        ddates <- dtsi[, .(RENAME_ID, visit_dt)]

predictors.tanya <- preprocess.tanya() 
predictors.ours <- preprocess.ours()
predictors.ours <- merge(predictors.ours, ddates, by.x='host.id', by.y='RENAME_ID')

# Save 1st plot of interest
x <- compare.predictors(predictors.ours, predictors.tanya)

# tmp <- x[['DT']]
# cols <- grep('PT_ID|maf',names(tmp), value=T)
# tmp[, ..cols]
# tmp

# 2nd plot: Extract all seroconverters
        
dsero <- fread(file.seroconverters)
dsero <- dsero[pangea_id %in% dsamples$pangea_id]
cols <- c('date_last_negative', 'date_first_positive')
new_cols <- c('knownTSI_max', 'knownTSI_min')
dsero[, (new_cols) := lapply(.SD , function(x) (visit_dt - x)/365), .SDcols=cols ]
dsero <- unique(dsero[, .(study_id, visit_dt, midpoint, knownTSI, knownTSI_min, knownTSI_max)])

dsero[ , grepl('^RK-',  study_id) ] |> all()
dtsi[! grepl('^RK-',  study_id) ] 

cat('Among all seroconveters from Rakai:\n\t')
p <- plot.timeline.all.seroconverters(DSERO=dsero, DTSI=dtsi)
cat('\t Coherency of:', scales::label_percent()(p$coh), '\n')
cat('\t MAE of:', round(p$mae, 2), 'years \n')
cat('\t MSE of:', round(p$mse, 2), 'years \n')
filename <- file.path(args$rel.dir, 'seroconv_timeline_all_rakai.png')
ggsave(filename,p$plot, width=10, height=14)

cat('Among seroconveters with visit_dt==sample collection date:\n\t')
p <- plot.timeline.all.seroconverters(DSERO=dsero[knownTSI_min == 0], DTSI=dtsi)
cat('\t Coherency of:', scales::label_percent()(p$coh), '\n')
cat('\t MAE of:', round(p$mae, 2), 'years \n')
cat('\t MSE of:', round(p$mse, 2), 'years \n')

cat('Among seroconveters with midpoint <= 10yr prior to visit_dt:\n\t')
p <- plot.timeline.all.seroconverters(DSERO=dsero[knownTSI <= 10], DTSI=dtsi)
cat('\t Coherency of:', scales::label_percent()(p$coh), '\n')
cat('\t MAE of:', round(p$mae, 2), 'years \n')
cat('\t MSE of:', round(p$mse, 2), 'years \n')



# TODO: some of the inds which Tanya treated as seroconverters 
# are not seroconverters in our eyes:
# if(0)
# {
#         ?write.csv
#         dtsi[, study_id[! study_id %in% dsero$study_id] ] |>
#                 write.table('220428_tanyaseroc_not_ourserconv.csv', 
#                             row.names=F, sep=',', col.names=F)
# }

dtsi <- merge(dtsi, dsero, by=c('study_id', 'visit_dt'))
dtsi <- unique(dtsi)
dtsi[, AID := gsub('-fq[0-9]', '', RENAME_ID)]
dtsi[, PT_ID := gsub('RK-|MRC-', '', study_id)]
dtsi <- dtsi[PT_ID %in% predictors.tanya$PT_ID,]


ttl=gsub('^.*?live/','',args$rel.dir)
ttl=gsub('/',' ',ttl)
ttl=gsub('_phsc_.*?$',' ',ttl)

g <- plot.sero.predictions(dtsi, title=ttl)
g

filename <- file.path(args$rel.dir, 'sero_predictions_vs_known.png')
ggsave(filename,g, width=10, height=8)

if(0)
{
        # Check blood samples collected on same date
        tmp <- dtsi[, list(rep=uniqueN(visit_dt) != .N, mult=.N >1) , by='AID']
        tmp[mult==TRUE,
            cat(sum(rep),'(', round(mean(rep*100),2),
                '%) participants with multiple associated blood samples had reported same collection date \n' ) ]

}



if(0)
{

  names(dpreds)
  .plot.comparison <- function(col)
  {

    # Extract participants with multiple associated predictions (ie. multiple blood samples)
    idx <- dpreds[, .N, by=AID ][N==2, AID ]
    rgx <- paste0('^', col, '$|AID|visit_dt')
    
    # Extract columns of interest
    cols <- grep(rgx, names(dpreds), value=T)
    tmp <- dpreds[AID %in% idx, ..cols]
    cols <- setdiff(cols, 'AID')

    # extract associated predictions + visit dates
    .f <- function(x, i) x[i]
    .f1 <- function(i){
      tmp1 <- tmp[, lapply(.SD, .f,i), .SDcols=cols, by='AID']
      idx <- !names(tmp1) %like% 'AID'
      names(tmp1)[idx] <- paste0(names(tmp1)[idx], '_' ,i)
      tmp1
    }  
    tmp <- lapply(c(1, 2), .f1)
    tmp <- Reduce(merge, tmp)
    tmp[, delta_t := as.numeric(visit_dt_2 - visit_dt_1)/365]
    tmp
    
    ttl <- fcase(
          col %like% 'sqrt', 'Square rooted TSI median',
          col %like% 'linear', 'Linear TSI median',
          col %like% 'pred_doi', 'Median date of infection '
    ) 

    g <- ggplot(tmp, aes_string(x=paste0(col, '_1'), y=paste0(col, '_2'))) +
      geom_point() + 
      geom_abline(intercept = 0, slope = 1, color='red', linetype='dotted') + 
      theme_bw() +
      labs(title=ttl)
    
    return(list(tmp=tmp, g=g))
  }
  names(dpreds)
  dpreds[,.N,by='AID'][N==2, .N]



  library(ggpubr)

  ggarrange(
          .plot.comparison('RF_pred_sqrt')$g,
          .plot.comparison('RF_pred_linear')$g,
          .plot.comparison('pred_doi_mid')$g,
          nrow=3
  ) -> p

  filename=file.path(args$rel.dir , "twosamples_predictions_comparison.png")
  ggsave(p,file=filename, w=8, h=14)


  
  .plot.comparison('RF_pred_sqrt')$tmp[, mean(RF_pred_sqrt_1 > RF_pred_sqrt_2)]
  .plot.comparison('RF_pred_linear')$tmp[, mean(RF_pred_linear_1 > RF_pred_linear_2)]
  .plot.comparison('pred_doi_mid')$tmp[, mean(pred_doi_mid_1 > pred_doi_mid_2)]
  
  .d <- function(a,b) cor(a,b)
  .d1 <- function(a, b) mean(abs(a-b))      # mean absolute error    
  .d2 <- function(a, b) mean((a-b)^2)
  tmp <- .plot.comparison('RF_pred_linear')$tmp
  
  cat('Correlations\n')
  tmp[, .d(RF_pred_linear_1, RF_pred_linear_2)]
  tmp[, .d(RF_pred_linear_1, RF_pred_linear_2 - delta_t)]

  tmp[, .d1(RF_pred_linear_1, RF_pred_linear_2)]
  tmp[, .d1(RF_pred_linear_1, RF_pred_linear_2 - delta_t)]
  tmp[, .d2(RF_pred_linear_1, RF_pred_linear_2)]
  tmp[, .d2(RF_pred_linear_1, RF_pred_linear_2 - delta_t)]
  ggplot(tmp, aes(x=RF_pred_linear_1, y=RF_pred_linear_2)) + 
          geom_abline(slope=1, linetype='dashed', color='red') +
          geom_label()
          geom_point() 

  ggplot(tmp, aes(x=RF_pred_linear_1, y=RF_pred_linear_2 - delta_t)) + 
          geom_abline(slope=1, linetype='dashed', color='red') +
          geom_point() 

  
  cat('Mean absolute "Error" accounting for time distance between visit dates:\n', 
      'NO : ', tmp[, .d1(RF_pred_linear_1, RF_pred_linear_2)],'\n',
      'YES: ', tmp[, .d1(RF_pred_linear_1, RF_pred_linear_2 - delta_t)], '\n')
  
  cat('Mean squared "Error" accounting for time distance between visit dates:\n', 
      'NO : ', tmp[, .d2(RF_pred_linear_1, RF_pred_linear_2)],'\n',
      'YES: ', tmp[, .d2(RF_pred_linear_1, RF_pred_linear_2 - delta_t)], '\n')
  
}

