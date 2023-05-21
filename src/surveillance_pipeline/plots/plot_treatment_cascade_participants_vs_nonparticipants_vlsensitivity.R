library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)
library(stringr)
library(ggsci)
library(gridExtra)
library(ggpubr)

# directory of the repository
gitdir <- here()
source(file.path(gitdir, "config.R"))

# directory to save the figure
outdir <- file.path('/home/andrea/HPC/project/ratmann_pangea_deepsequencedata/live','PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')
if(usr == 'melodiemonod'){
  outdir <- file.path('~/Box\ Sync/2021/ratmann_deepseq_analyses/live/', 'PANGEA2_RCCS', 'treatment_cascade_by_gender_loc_age')
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# check files exist
file.exists(c(
  file.treatment.cascade.prop.participants,
  file.treatment.cascade.prop.nonparticipants,
  file.treatment.cascade.prop.participants.vl200 ,
  file.treatment.cascade.prop.nonparticipants.vl200))  |> all() |> stopifnot()


##############

# PRE-PROCESS NON PARTICIPANTS

##############

# load files
cascade.np <- as.data.table(read.csv(file.treatment.cascade.prop.nonparticipants))
cascade.np200 <- as.data.table(read.csv(file.treatment.cascade.prop.nonparticipants.vl200))

# combine
cascade.np[, type := 'Viral load threshold 1,000 mL']
cascade.np200[, type := 'Viral load threshold 200 mL']
cascade.all.var <- rbind(cascade.np, cascade.np200)

# select variable of interst for cascade
cascade <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,
                               PROP_DIAGNOSED_M, PROP_DIAGNOSED_CL, PROP_DIAGNOSED_CU,
                               PROP_ART_COVERAGE_M, PROP_ART_COVERAGE_CL, PROP_ART_COVERAGE_CU,
                               SUPPRESSION_RATE_M, SUPPRESSION_RATE_CL, SUPPRESSION_RATE_CU)]
setnames(cascade, c('PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'),
         c('PROP_ART_COVERAGE_GIVEN_DIAGNOSED_M', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CL', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CU'))
setnames(cascade, c('SUPPRESSION_RATE_M', 'SUPPRESSION_RATE_CL', 'SUPPRESSION_RATE_CU'),
         c('PROP_SUPPRESSION_GIVEN_ART_M', 'PROP_SUPPRESSION_GIVEN_ART_CL', 'PROP_SUPPRESSION_GIVEN_ART_CU'))

# melt
cascade <- melt.data.table(cascade, id.vars = c('ROUND', 'SEX', 'AGEYRS', 'COMM', 'type'))

# find quantiles index
cascade[, q := gsub('_', '', str_sub(variable, start= -2))]
cascade[, variable := gsub('(.+)_.*', '\\1',variable)]

# cast
cascade <- dcast.data.table(cascade, ROUND + SEX + AGEYRS + COMM + type + variable ~ q, value.var = 'value')


###########################

# PLOT PROPORTION OF UNSUPPRESSED

###########################

# select variable of interst for cascade
uns <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,PROP_UNSUPPRESSED_EMPIRICAL,
                           PROP_UNSUPPRESSED_M, PROP_UNSUPPRESSED_CL, PROP_UNSUPPRESSED_CU)]
uns[,SEX_LABEL := 'Men']
uns[SEX == 'F',SEX_LABEL := 'Women']
uns[,COMM_LABEL := 'Inland communities']
uns[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
uns[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

# for participants
tab <- uns[COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = type), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL, col = type, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = type, linetype= 'Fit')) + 
  facet_grid(ROUND_LABEL~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with unsuppressed viral load', 
       col = '', fill = '', linetype = '', shape = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))

# save
file.name <- file.path(outdir, 'non_participants_smooth_unsuppressed_vl200vs1000_221208.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(file = file.name, w = 7, h = 12)  
} else {
  cat("File:", file.name, "already exists...\n")
}


##############

# PRE-PROCESS PARTICIPANTS

##############

# load files
cascade <- as.data.table(read.csv(file.treatment.cascade.prop.participants))
cascade.200 <- as.data.table(read.csv(file.treatment.cascade.prop.participants.vl200))

# combine
cascade[, type := 'Viral load threshold 1,000 mL']
cascade.200[, type := 'Viral load threshold 200 mL']
cascade.all.var <- rbind(cascade, cascade.200)

# select variable of interst for cascade
cascade <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,
                               PROP_DIAGNOSED_M, PROP_DIAGNOSED_CL, PROP_DIAGNOSED_CU,
                               PROP_ART_COVERAGE_M, PROP_ART_COVERAGE_CL, PROP_ART_COVERAGE_CU,
                               SUPPRESSION_RATE_M, SUPPRESSION_RATE_CL, SUPPRESSION_RATE_CU)]
setnames(cascade, c('PROP_ART_COVERAGE_M', 'PROP_ART_COVERAGE_CL', 'PROP_ART_COVERAGE_CU'),
         c('PROP_ART_COVERAGE_GIVEN_DIAGNOSED_M', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CL', 'PROP_ART_COVERAGE_GIVEN_DIAGNOSED_CU'))
setnames(cascade, c('SUPPRESSION_RATE_M', 'SUPPRESSION_RATE_CL', 'SUPPRESSION_RATE_CU'),
         c('PROP_SUPPRESSION_GIVEN_ART_M', 'PROP_SUPPRESSION_GIVEN_ART_CL', 'PROP_SUPPRESSION_GIVEN_ART_CU'))

# melt
cascade <- melt.data.table(cascade, id.vars = c('ROUND', 'SEX', 'AGEYRS', 'COMM', 'type'))

# find quantiles index
cascade[, q := gsub('_', '', str_sub(variable, start= -2))]
cascade[, variable := gsub('(.+)_.*', '\\1',variable)]

# cast
cascade <- dcast.data.table(cascade, ROUND + SEX + AGEYRS + COMM + type + variable ~ q, value.var = 'value')


###########################

# PLOT PROPORTION OF UNSUPPRESSED

###########################

# select variable of interst for cascade
uns <- cascade.all.var[, .(ROUND, SEX, AGEYRS, COMM, type,PROP_UNSUPPRESSED_EMPIRICAL,
                           PROP_UNSUPPRESSED_M, PROP_UNSUPPRESSED_CL, PROP_UNSUPPRESSED_CU)]
uns[,SEX_LABEL := 'Men']
uns[SEX == 'F',SEX_LABEL := 'Women']
uns[,COMM_LABEL := 'Inland communities']
uns[COMM == 'fishing',COMM_LABEL := 'Fishing communities']
uns[,ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]

# for participants
tab <- uns[COMM == 'inland' & ROUND != 'R015S']
ggplot(tab, aes(x = AGEYRS)) + 
  geom_ribbon(aes(ymin = PROP_UNSUPPRESSED_CL, ymax = PROP_UNSUPPRESSED_CU, fill = type), alpha = 0.5) + 
  geom_point(aes(y = PROP_UNSUPPRESSED_EMPIRICAL, col = type, shape = 'Data'), alpha = 0.5) + 
  geom_line(aes(y = PROP_UNSUPPRESSED_M, col = type, linetype= 'Fit')) + 
  facet_grid(ROUND_LABEL~SEX_LABEL) + 
  theme_bw() + 
  labs(x = 'Age', y = 'Proportion of HIV-positive participants with unsuppressed viral load', 
       col = '', fill = '', linetype = '', shape = '') +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white")) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))  + 
  guides(shape = guide_legend(order = 1), linetype = guide_legend(order = 2), 
         color = guide_legend(order = 3),fill = guide_legend(order = 3))

file.name <- file.path(outdir, 'participants_smooth_unsuppressed_vl200vs1000_221208.pdf')
if (! file.exists(file.name) || config$overwrite.existing.files) {
  cat("Saving file:", file.name, "\n")
  ggsave(file = file.name, w = 7, h = 12)  
} else {
  cat("File:", file.name, "already exists...\n")
}



