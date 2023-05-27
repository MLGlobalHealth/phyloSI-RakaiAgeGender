cat("Start of get_main_figure4.R")

{
    library(rstan)
    library(data.table)
    library(ggplot2)
    library(ggpubr)
    library(gridExtra)
    library(matrixStats)
    library(dplyr)
    library(lubridate)
    library(ggnewscale)
    library(patchwork)
} |> suppressPackageStartupMessages()

gitdir <- here::here()
source(file.path(gitdir, 'config.R'))
# load functions
source(file.path(gitdir.R.flow, 'postprocessing_summary_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_plot_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_utils_functions.R'))
source(file.path(gitdir.R.flow, 'postprocessing_statistics_functions.R'))


jobname <- 'central3'
stan_model <- 'gp_221201d'
suffix <- paste(stan_model, jobname, sep='-')

# if on HPC, we know the username
if(dir.exists('/rds/')){
    if(usr == 'ab1820')
        outdir <- paste0("/rds/general/user/",usr,"/home/projects/2022/phyloflows/", suffix)
    if(usr == 'mm3218')
        outdir <- paste0("/rds/general/user/",usr,"/home/projects/2021/phyloflows/", suffix)
} else {
    if(usr == 'andrea'){
        outdir <- file.path(outdir.phyloflows, suffix) 
    }
}


outfile <- file.path(outdir, paste0(stan_model,'-', jobname))

# paths
path.to.stan.output = paste0(outfile, "-stanout_", jobname, ".rds")
.outfile.figures <- file.path(outdir, 'figures',suffix)
.outdir.table <-  file.path(outdir, 'tables', suffix)
if(0)
    path.to.suboutput <- '~/Downloads/subsample.rds'

# load data
path.to.stan.data <- paste0(outfile, "-stanin_",jobname,".RData")
load(path.to.stan.data)

# samples: select the thinned sample if specified
if(exists('path.to.suboutput') )
{
    stopifnot("Thinned chains not found"=file.exists(path.to.suboutput))
    .outfile.figures <- '~/Downloads' 
    .outdir.table <- '~/Downloads' 
    samples <- readRDS( file=path.to.suboutput)
}else{
    fit <- readRDS(path.to.stan.output)
    samples <- rstan::extract(fit)
}

outfile.figures <- .outfile.figures
outdir.table <- .outdir.table

## Function to extract legend
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 


# need to be able to get to: 
# a -output-counterfactual_budget_age_Unsuppressed_inland.pdf
# b -output-counterfactual_budget_incidence_panel_plot1_Unsuppressed_inland.pdf
# c -output-counterfactual_budget_incidence_panel_plot2_Unsuppressed_inland.pdf
# d -output-counterfactual_budget_incidence_panel_plot3_Unsuppressed_inland 

## UNSUPPRESSED ## 

log_offset_round <- find_log_offset_by_round(
    stan_data,
    eligible_count_round,
    df_estimated_contact_rates,
    use_number_susceptible_offset,
    use_contact_rates_prior)


# log offset formula (per year per susceptible)
log_offset_formula_persusceptible <- 'log_INFECTED_NON_SUPPRESSED'
if(use_contact_rates_prior)
  log_offset_formula_persusceptible = paste0(log_offset_formula_persusceptible, ' + log_CONTACT_RATES')


cat("\nPlot relative incidence infection if different groups of male are targeted\n")


# find incidence under the factual scenario by sex and age
# THIS SEEMS COMPLETELY WRONG? 
incidence_factual <- find_summary_output_by_round(
    samples,
    'log_beta',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_INFECTION.RECIPIENT'),
    transform = 'exp',
    log_offset_formula = log_offset_formula_persusceptible,
    log_offset_round = log_offset_round)

expected_contribution_age_source <- find_summary_output_by_round(samples,
    'log_lambda_latent',
    c('INDEX_DIRECTION', 'INDEX_ROUND', 'AGE_TRANSMISSION.SOURCE'),
    transform = 'exp',
    standardised.vars = c('INDEX_DIRECTION', 'INDEX_ROUND'))

# identify the main spreaders, sources that contributes the most to incidence
spreaders <- find_spreaders(expected_contribution_age_source, outdir.table)

counterfactuals_p_a <- make_counterfactual_target(samples,
    spreaders,
    log_offset_round,
    stan_data,
    eligible_count_smooth,
    eligible_count_round,
    treatment_cascade,
    proportion_prevalence,
    participation,
    only_participant = T,
    art_up_to_female = F,
    outdir.table)

counterfactuals_a_a <- make_counterfactual_target(samples,
    spreaders,
    log_offset_round,
    stan_data,
    eligible_count_smooth,
    eligible_count_round,
    treatment_cascade,
    proportion_prevalence,
    participation,
    only_participant = F,
    art_up_to_female = F,
    outdir.table)

if(exists("treatment_cascade_samples"))
{
    counterfactuals_a_f <- make_counterfactual(
        samples,
        log_offset_round,
        stan_data, 
        eligible_count_smooth,
        eligible_count_round, 
        treatment_cascade_samples,
        proportion_prevalence,
        participation,
        only_participant = F,
        art_up_to_female = 1,
        s959595 = NULL,
        s909090 = NULL,
        outdir.table)

    #  generate counterfactual treating all men half way to as much as female are diagnosed/treated/suppressed
    counterfactuals_a_f05 <- make_counterfactual(
        samples,
        log_offset_round,
        stan_data, 
        eligible_count_smooth,
        eligible_count_round, 
        treatment_cascade_samples,
        proportion_prevalence,
        participation,
        only_participant = F,
        art_up_to_female = 0.5,
        s959595 = NULL,
        s909090 = NULL,
        outdir.table)

    # generate counterfactual treating all men 95 95 95
    counterfactuals_a_959595 <- make_counterfactual(
        samples,
        log_offset_round,
        stan_data, 
        eligible_count_smooth,
        eligible_count_round, 
        treatment_cascade_samples,
        proportion_prevalence,
        participation,
        only_participant = F,
        art_up_to_female = NULL,
        s959595 = 1,
        s909090 = NULL,
        outdir.table)

    # generate counterfactual treating all men 90 90 90
    counterfactuals_a_909090 <- make_counterfactual(samples,
        log_offset_round,
        stan_data, 
        eligible_count_smooth,
        eligible_count_round, 
        treatment_cascade_samples,
        proportion_prevalence,
        participation,
        only_participant = F,
        art_up_to_female = NULL,
        s959595 = NULL,
        s909090 = 1,
        outdir.table)

    # plot
    naturemed_reqs()
    fig4_list <- plot_counterfactual(
        counterfactuals_a_f,
        counterfactuals_a_f05,
        counterfactuals_a_959595,
        counterfactuals_a_909090,
        incidence_factual,
        "Unsuppressed",
        outfile.figures, 
        nm_reqs = TRUE
    )
}

# remember to add labels
p_a <-fig4_list[['a']]
p_b <-fig4_list[['b']] 
p_c <-fig4_list[['c']]
p_d <-fig4_list[['d']]

# Prepare legend grabbing it from 
cols <- c('#F1A661', '#C55300', '#749F82')
p_a_labs <- c(
    'Remaining virally unsuppressed in R18',
    'Closing the suppression gap\nin men relative to women',
    'Closing half the suppression gap\nin men relative to women',
    '95-95-95 in men',
    'Already virally suppressed in R18'
)

p_a_legend <- p_a  +
    labs(fill='Counterfactual interventions', color='Counterfactual interventions') +
    theme(legend.position='bottom') +
    scale_fill_manual(values = c('grey50', cols[c(2,1,3)], 'grey80'), labels=p_a_labs)

p_d_legend <- p_d +
    theme(legend.position='bottom') + 
    guides(fill='none', color='none')

legend_a <- get_legend(p_a_legend)
legend_d <- get_legend(p_d_legend)

legends <- ggarrange_nature(
    as_ggplot(legend_a),
    as_ggplot(legend_d),
    ncol=2)

# consider removing legend for a 

p_bcd <- ggarrange_nature(
    p_b , p_c , p_d,
    ncol=1,
    heights=c(3.3, 3.3, 4.6),
    labels = c('b', 'c', 'd'))

p_a_with_legend <- ggarrange_nature(
    p_a , legends,
    ncol = 1,
    heights = c(4, 2),
    labels = c('a', '')
)

p_abcd <- ggarrange_nature(
    p_a_with_legend, p_bcd,
    ncol=2,
    widths=c(3.2, 3.4)
)

ggsave_nature(p_abcd,
    filename=paste0(outfile.figures, 'MainFigure4.pdf'),
    w=18, h=15)
