if(TRUE)
{
    library(data.table)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(knitr)
    require(lubridate)
    library(rstan)
    library(gridExtra)
    library(lognorm)
    library(ggExtra)
    library(Hmisc)
    require(ggplotify)
}


gitdir <- here::here()

source(file.path(gitdir, 'config.R' ))
source(file.path(gitdir.R, 'utils.R'))
source(file.path(gitdir.R.flow, 'summary_functions.R'))
source(file.path(gitdir.R.flow, 'plotting_functions.R'))
source(file.path(gitdir.R.flow, 'statistics_functions.R'))
source(file.path(gitdir.R.flow, 'stan_utils.R'))

naturemed_reqs()
find_palette_round()


# indicators -- sensitivity analyses

#
# Define start time, end time and cutoff
#

load(file.path.round.timeline)
df_round_inland[, `:=` (min_sample_date = as.Date(min_sample_date), max_sample_date = as.Date(max_sample_date))]

start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

stopifnot(start_first_period_inland < stop_first_period_inland)
stopifnot(stop_first_period_inland < start_second_period_inland)
stopifnot(start_second_period_inland < stop_second_period_inland)

# 2nd line needed for ROUND_SPANYRS
df_round <- make.df.round(df_round_inland)
df_period <- make.df.period(start_first_period_inland, stop_first_period_inland, 
                            start_second_period_inland, stop_second_period_inland, 
                            df_round)


# load pairs
pairs.all <- read_pairs(file.pairs)


#
# Find phylo pairs 

#

pairs <- select.pairs.for.analysis(pairs.all, 
    only.one.community='inland', 
    use_30com_pairs=FALSE, 
    only.transmission.after.start.observational.period=TRUE, 
    only.transmission.before.stop.observational.period=TRUE, 
    remove.pairs.from.rounds=NULL
)

# keep only pairs with source-recipient with a time of infection
pairs <- pairs[!is.na(AGE_TRANSMISSION.SOURCE) & !is.na(AGE_INFECTION.RECIPIENT)]
pairs[, DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT := DATE_INFECTION.RECIPIENT < start_second_period_inland]
tab <- pairs[, list(count = .N), by = c('DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT', 'COMM.RECIPIENT', 'SEX.RECIPIENT')]
print_table(tab[order(DATE_INFECTION_BEFORE_CUTOFF.RECIPIENT, COMM.RECIPIENT, SEX.RECIPIENT)])

# outfile.figures <- '~/Downloads/EDF5/inkscape'


p <- plot_pairs(pairs.all[BOTH_PARTICIPATED == TRUE & SEX.RECIPIENT != SEX.SOURCE], outfile.figures, nm_reqs = TRUE)
p1 <- plot_pairs_all(pairs.all[BOTH_PARTICIPATED==TRUE], outfile.figures, nm_reqs=TRUE)
p2 <- plot_transmission_events_over_time(pairs.all[BOTH_PARTICIPATED==TRUE], outdir = NULL, nm_reqs=TRUE)

# extract legend
legend <- cowplot::get_legend(p2)
p2 <- p2 + theme(legend.position = 'none')

require(patchwork)
top <- (p1 + p2 + plot_layout(widths = c(1.2,3))) / 
    legend + 
    plot_layout( heights = c(100,1)) + 
    plot_annotation(tag_levels='a') & theme(plot.tag = element_text(size=8, face='bold', family='sans')) +
    reqs

p <- as.ggplot(p) + reqs
all <- ggarrange(top + reqs, p + reqs , heights=c(1.8,3), ncol=1)
ggsave_nature('~/Downloads/EDF5/extended_data_figure_pairs.pdf',all, w=18, h=19)
