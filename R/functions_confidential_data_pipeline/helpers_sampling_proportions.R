cat("\n_ Sourcing 'helpers_sampling_proportions _\n")

get_hivids <- function( agerange=c(15,49))
{

    meta_data <- fread(file.path.metadata) #additional meta_data

    # find age
    meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
    meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
    meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]

    # restrict age
    meta_data <- meta_data[AGEYRS %between% agerange ]

    # find community
    meta_data[, COMM := 'inland']
    meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']

    # find hiv status
    meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]

    # subset round 14 to 30 continuously surveyed communities
    meta_data <- meta_data[!(round=='14' & !COMM %in% c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754))]
    
    # find art use
    meta_data[, ART := artslfuse == 'yes']

    # keep variable of interest
    meta_data[, round := paste0('R0', round)]
    colnames(meta_data) <- toupper(colnames(meta_data))
    meta_data <- meta_data[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV, ART)]

    # set 15.1 to be 15S
    meta_data[ROUND == 'R015.1', ROUND := 'R015S']


    #
    # Quest
    #

    # keep variable of interest, merge communities
    rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, arvmed, cuarvmed)]
    rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
    colnames(rinc) <- toupper(colnames(rinc))
    rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

    # get ART status
    rinc[!ROUND %in% c('R016', 'R017', 'R018'), ART := ARVMED ==1]
    rinc[!ROUND %in% c('R016', 'R017', 'R018') & is.na(ARVMED), ART := F]
    rinc[ROUND == 'R016', ART := ARVMED ==1 | CUARVMED ==1]
    rinc[ROUND == 'R016' & (is.na(ARVMED) | is.na(CUARVMED)), ART := F]
    rinc[ROUND %in% c('R017', 'R018'), ART := CUARVMED ==1]
    rinc[ROUND %in% c('R017', 'R018') & is.na(CUARVMED), ART := F]

    # art was not reported in round 10
    rinc[ROUND == 'R010', ART := NA]

    # add meta data from Kate
    tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
    tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
    # rinc <- merge(rinc, meta_data, by=c('STUDY_ID', 'ROUND', 'SEX', 'AGEYRS','ART'), all.x=TRUE)

    rinc <- rbind(
        rinc[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)], 
        tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)]
    )

    #
    # HIV
    #

    hiv <- fread(file.path.hiv)
    hiv[, round := gsub(' ', '', round)] # remove space in string

    # get hiv status
    rhiv <- hiv[, .(study_id, round, hiv)]
    rhiv[, round := gsub(" ", '', round, fixed = T)]
    colnames(rhiv) <- toupper(colnames(rhiv))

    # add meta data from Joseph 
    hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

    # add meta data from Kate
    tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], hivs[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
    tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
    hivs <- rbind(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)])

    # SET ROUND 15S IN INLAND AS 15
    hivs[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
    hivs[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
    hivs <- hivs[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]
    hivs <- hivs[AGEYRS %between% agerange ]
    return(hivs)

}

extract_posterior_suppression_samples <- function(DT=NULL, ps = c(M=0.5, CL=0.025, CU=0.975))
{
    by_cols <- c('SEX', 'COMM', 'ROUND')

    if(is.null(DT)){
        n_infected <- add_susceptible_infected(
            eligible_count_smooth = fread(file.eligible.count),
            proportion_prevalence = fread(file.prevalence.prop),
            participation = fread(file.participation),
            nonparticipants.male.relative.infection  =  1, 
            nonparticipants.female.relative.infection  =  1
        ) 
    }else{
        n_infected <- copy(DT)
    }
    n_infected <- subset(n_infected, select=c(by_cols, 'AGEYRS','INFECTED'))

    # extract suppression samples and prettify
    samples_suppression <- readRDS(file.treatment.cascade) |> 
        subset(COMM == 'inland', select=c(by_cols, 'AGEYRS','iterations', 'PROP_SUPPRESSED_POSTERIOR_SAMPLE')) |>
        within(ROUND <- gsub('R0', '', ROUND))

    # get number unsuppressed among pop
    samples_suppression <- merge( samples_suppression, n_infected, by=c(by_cols, 'AGEYRS'))
    samples_suppression[, N_UNSUPPRESSED := INFECTED * (1-PROP_SUPPRESSED_POSTERIOR_SAMPLE) ]
    samples_suppression[, AGEGP := ageyrs2agegp(AGEYRS)]
    samples_suppression <- samples_suppression[, .(
        N_UNSUPPRESSED = sum(N_UNSUPPRESSED)),
        by=c(by_cols, 'AGEGP','iterations')]

    # get quantiles
    samples_suppression <- samples_suppression[, {
        z <- as.list(quantile(N_UNSUPPRESSED, probs = ps))
        names(z) <- paste0('INFECTED_NON_SUPPRESSED_',names(ps))
        z
    }, by=c(by_cols, 'AGEGP') ]
    samples_suppression[, ROUND := paste0('R0', ROUND)]
    samples_suppression
}

load.viralsuppression.test.results <- function(QUEST=quest)
{
    # VL_DETECTABLE = 400 # Not needed
    VIREMIC_VIRAL_LOAD = 1000 # WHO standards

    # Load data: exclude round 20 as incomplete
    dall <- fread(path.tests)
    dall <- dall[ROUND %in% c(15:18)]
    # dall <- dall[ROUND == round]

    # rename variables according to Oli's old script + remove 1 unknown sex
    setnames(dall, c('HIV_VL', 'COMM'), c('VL_COPIES', 'FC') )
    dall[, HIV_AND_VL := ifelse( HIV_STATUS == 1 & !is.na(VL_COPIES), 1, 0)]
    dall <- dall[! SEX=='']

    # remove HIV+ individuals with missing VLs  
    DT <- subset(dall, HIV_STATUS==0 | HIV_AND_VL==1)

    # set ARVMED to 0 for HIV-
    set(DT, DT[, which(ARVMED==1 & HIV_STATUS==0)], 'ARVMED', 0) 

    # define VLC as VL_COPIES for HIV+ and as 0 for HIV-
    DT[,  VLC := fifelse(HIV_STATUS == 0, yes=0, no=VL_COPIES)]

    # define suppressed VL as VLS and unsuppressed as VLNS (according to WHO criteria)	
    DT[, `:=` (
        VLS = as.integer(VLC<VIREMIC_VIRAL_LOAD),
        VLNS = as.integer(VLC>=VIREMIC_VIRAL_LOAD)
    )]

    # reset undetectable to VLC 0
    setkey(DT, ROUND, FC, SEX, AGEYRS)

    # keep within census eligible age
    DT <- subset(DT, AGEYRS > 14 & AGEYRS < 50)

    # keep infected
    DT <- DT[HIV_STATUS ==1]

    # use age from quest
    DT[, round := paste0('R0', ROUND)]
    set(DT, NULL, 'AGEYRS', NULL)
    # DT <- merge(DT, hiv, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
    DT <- merge(DT, QUEST, by.x = c('round', 'STUDY_ID'), by.y = c('round', 'study_id'))
    setnames(DT, 'ageyrs', 'AGEYRS')
}

prettify_labs <- function(DT){
    nms <- names(DT)

    if (! 'ROUND_LAB' %in% nms & 'ROUND' %in% nms )
        DT[, ROUND_LAB := gsub("R0", "Round ", ROUND)]

    if (! 'SEX_LAB' %in% nms)
        DT[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")] 

    return(DT)
}

.make.plot.with.binconf <- function(DT, x, n, xvar=ROUND_LAB, .ylab, ylims = c(0,1) )
{
    .xvar <- enexpr(xvar)
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot <- prettify_labs(dplot)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]

    if(deparse(.xvar) == 'ROUND_LAB'){
        dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]
    }
    dplot[, AGEGP_LAB := paste(AGEGP, 'years')]

    ggplot(dplot, aes(x=eval(.xvar), color=SEX_LAB, pch=AGEGP_LAB, linetype=AGEGP_LAB, y=P )) + 
        geom_hline(yintercept = 1, linetype='dashed', color='grey50') +
        geom_point(position=position_dodge(width=.8) ) + 
        geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width =.8)) +
        scale_y_continuous(limits = ylims, expand=c(0,0),labels=scales::percent) +
        scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        labs( 
            x = NULL,
            y = .ylab,
            linetype = "Age group", 
            pch = "Age group", 
            color = NULL,
        ) +
        theme_bw() + 
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        rotate_x_axis(30) +
        reqs
}

.binconf.ratio.plot <- function(DT, x, n, .ylab, xaes=expr(ROUND_LAB), .yrange=NA)
{
    x <- enexpr(x)
    n <- enexpr(n)

    dplot <- copy(DT)
    dplot[, (c('P', 'CL', 'CU')) := binconf(x=eval(x), n=eval(n), return.df=TRUE) ]

    if('ROUND_LAB' %in% names(DT))
        dplot[, ROUND_LAB := gsub('^R0','Round', ROUND_LAB)]

    by_col <- fifelse('PERIOD' %in% names(DT), yes='PERIOD', no='ROUND' )
    dplot[, (c('LOG_RATIO_P', 'LOG_RATIO_CL', 'LOG_RATIO_CU')):= {
        P.both <- P[which(SEX == 'Both')]
        list( 
            LOG_RATIO_P = log(P/P.both),
            LOG_RATIO_CL = log(CL/P.both),
            LOG_RATIO_CU = log(CU/P.both)
        )
    }, by=by_col]

    dplot |> subset(SEX_LAB != 'Both' & AGEGP != 'All') |> 
        ggplot(aes(x=!!xaes, color=SEX_LAB, pch=AGEGP, linetype=AGEGP, y=LOG_RATIO_P )) + 
        geom_hline(yintercept = 0, linetype='dashed', color='grey50') +
        geom_point(position=position_dodge(width=.8) ) + 
        geom_linerange(aes(ymin=LOG_RATIO_CL, ymax=LOG_RATIO_CU), position=position_dodge(width =.8)) +
        scale_y_continuous(limits = .yrange, expand=c(0,0)) +
        scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
        # coord_cartesian(ylim = ylims) +
        labs( 
            x = NULL,
            y = .ylab,
            linetype = "Age group", 
            pch = "Age group", 
            color = NULL,
        ) +
        theme_bw() + 
        theme(legend.position = "bottom", strip.background = element_blank()) + 
        rotate_x_axis(30) +
        reqs
}

# plot.hist.numerators.denominators <- function(DT, DRANGE, filltext=NA_character_)
# {
#     # can show agegroups as different bar fills! How to do this? 
#     dplot <- copy(DT)
#     dplot[, variable_lab := fifelse(variable == 'N_EVERSEQ', 'Ever-sequenced', 'Never deep-sequenced') ]
#     lims <- dplot[ , levels(interaction(variable_lab, SEX_LAB))[c(1,3)]] 
# 
#     ggplot(dplot, aes(x=ROUND_LAB, pch=AGEGP, y=value, fill=interaction(variable_lab, SEX_LAB) )) +
#         geom_col(data=dplot[variable == 'INFECTED_NON_SUPPRESSED_M'], position=position_dodge(width=.9), width=.9, color='black') + 
#         geom_linerange(data=DRANGE, position = position_dodge(width =.9),
#             aes(y=NULL, ymin=INFECTED_NON_SUPPRESSED_CL, ymax=INFECTED_NON_SUPPRESSED_CU, fill=NULL)) + 
#         # geom_col_pattern(
#         #     data=dplot[variable == 'N_EVERSEQ'],
#         #     aes(pattern=AGEGP),
#         #     position=position_dodge(width=.9), color='black', width=.9,
#         #     pattern_fill = 'white', pattern_color='white', pattern_density = .2, pattern_spacing=.008) +
#         facet_grid(SEX_LAB~.) + 
#         theme_bw() + 
#         scale_y_continuous(expand = expansion(c(0,0.1)) ) +
#         # scale_pattern_manual(values = c('15-24' = 'circle' , '25-34'= 'none', '35-49' = 'stripe' )) + 
#         # guides(pattern = guide_legend(override.aes = list(fill = "#BDC5D0"), pattern_density = .1)) +
#         guides(fill = guide_legend(override.aes = list(pattern = "none"))) +
#         scale_fill_manual(
#             values=fillpalette,
#             labels=filllabels,
#             limits = lims,
#             na.value = 'white'
#         ) +  
#         theme(legend.position = "bottom", strip.background = element_blank()) + 
#         labs(
#             x=NULL,
#             y="Estimated number of individuals with\nunsuppressed HIV in the population",
#             fill=NULL,
#             pattern='Age',
#             color=NULL
#         ) +
#         rotate_x_axis(30)
# }


plot.hist.numerators.denominators <- function(DT, DRANGE, filltext=NA_character_)
{
    # can show agegroups as different bar fills! How to do this? 
    # DT <- copy(dplot); DRANGE <- copy(tab_seq_unsup)
    dplot <- copy(DT)
    dplot[, variable_lab := fifelse(variable == 'N_EVERSEQ', 'Yes', 'No') ]
    dplot[, AGEGP_LAB := paste(AGEGP, "years")]
    
    ggplot(dplot, aes(x=ROUND_LAB, pch=AGEGP, y=value, fill=AGEGP_LAB )) +
        geom_col_pattern(
            data=dplot[variable_lab == 'No'],
            color='black',
            aes(group=AGEGP_LAB, pattern=variable_lab),
            position=position_dodge(width=.8), width=.8) +
        geom_col_pattern(
            data=dplot[variable_lab == 'Yes'],
            color='black',
            pattern_fill = "black",
            pattern_angle = 45,
            pattern_density = 0.1,
            pattern_spacing = 0.025,
            pattern_key_scale_factor = 0.6,
            aes(pattern=variable_lab, group=AGEGP_LAB),
            position=position_dodge(width=.8),  width=.8) +
        geom_linerange(data=DRANGE, position = position_dodge(width =.8),
            aes(y=NULL, ymin=INFECTED_NON_SUPPRESSED_CL, ymax=INFECTED_NON_SUPPRESSED_CU, fill=NULL)) + 
        facet_grid(SEX_LAB~.) + 
        theme_bw() + 
        scale_y_continuous(expand = expansion(c(0,0.1)) ) +
        scale_pattern_manual(values = c(Yes = "stripe", No = "none"), breaks=c('Yes', 'No')) +
        scale_fill_viridis_d( ) +
        guides(
            fill = guide_legend(override.aes = list(pattern="none")),
            pattern = guide_legend(override.aes = list(fill = "#BDC5D0", pattern_density = .001, pattern_spacing = .01))
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom", 
            strip.background = element_blank(),
            panel.grid.minor= element_blank(),
            panel.grid.major= element_blank()
        ) +
        labs(
            x=NULL,
            y="Estimated number of individuals with\nunsuppressed HIV in the population",
            fill="Age group",
            pattern= "",
            NULL
        ) +
        rotate_x_axis(30)
}

plot.hist.numerators.denominators.2 <- function(DT, range=FALSE, filltext=NA_character_)
{
    # DT <- copy(detectionprob)
    by_cols <- c('ROUND', 'SEX', 'COMM', 'AGEGP')
    dplot <- melt( DT, id.vars = by_cols, measure.vars = c('count', 'INCIDENT_CASES')) |> suppressWarnings()
    dplot <- prettify_labs(dplot)
 
    if(range){
        drange <- subset(DT, select=c('INCIDENT_CASES_LB', 'INCIDENT_CASES_UB', 'ROUND', 'AGEGP', 'SEX'))
        drange <- prettify_labs(drange)
        drange[, AGEGP_LAB := paste(AGEGP, "years")]
    }
    
    # relabel
    dplot[, variable_lab := fifelse(variable == 'count', 'Yes', 'No') ]
    dplot[, variable_lab := factor(variable_lab, levels = c('Yes', 'No'))]
    dplot[, AGEGP_LAB := paste(AGEGP, "years")]

    ggplot(dplot, aes(x=ROUND_LAB,  y=value, fill=AGEGP_LAB )) +
        geom_col_pattern(
            data=dplot[variable_lab == 'No'],
            color='black',
            aes(group=AGEGP_LAB, pattern=variable_lab),
            position=position_dodge(width=.8), width=.8) +
        geom_col_pattern(
            data=dplot[variable_lab == 'Yes'],
            color='black',
            pattern_fill = "black",
            pattern_angle = 45,
            pattern_density = 0.1,
            pattern_spacing = 0.025,
            pattern_key_scale_factor = 0.6,
            aes(pattern=variable_lab, group=AGEGP_LAB),
            position=position_dodge(width=.8),  width=.8) + {
            if(range) 
            geom_linerange(data=drange,position = position_dodge(width =.8), color='black', size=.5,
                aes(ymin=INCIDENT_CASES_LB, ymax=INCIDENT_CASES_UB, y=NULL))
        }  +
        facet_grid(SEX_LAB~.) +
        scale_pattern_manual(values = c(Yes = "stripe", No = "none"), breaks=c('Yes', 'No')) +
        scale_y_continuous(expand = expansion(c(0,0.1)) ) +
        scale_fill_viridis_d() +
        guides(
            fill = guide_legend(override.aes = list(pattern="none")),
            pattern = guide_legend(override.aes = list(fill = "#BDC5D0", pattern_density = .001, pattern_spacing = .01))
        ) +
        theme_bw() +
        theme(
            legend.position = "bottom", 
            strip.background = element_blank(),
            panel.grid.minor= element_blank(),
            panel.grid.major= element_blank()
        ) +
        labs(
            x=NULL,
            y="Estimated number of\ninfection events",
            fill="Age group",
            pattern= "",
            NULL
        ) +
        rotate_x_axis(30)
}
