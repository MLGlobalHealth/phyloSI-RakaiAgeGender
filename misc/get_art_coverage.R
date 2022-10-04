## find art coverage by combining proportion of suppressed viral load and self-reported art coverage. 

library(data.table)
library(ggplot2)
require(lubridate)
library(dplyr)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'

outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'artcoverage_by_gender_loc_age')

# summarised posterior samples
file.unsuppressedviralload.s <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_vl_1000_220803.csv'))
file.selfreportedart.s <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_selfreportedart_estimates_220906.csv'))

# posterior samples
file.unsuppressedviralload <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0('RCCS_nonsuppressed_proportion_posterior_samples_vl_1000_220818.csv'))
file.selfreportedart <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_selfreportedart_posterior_samples_220906.csv'))

####################################

# COMBINE SUMMARISED POSTERIOR SAMPLES (i.e., median, cl, cu)

####################################

uns <- as.data.table(read.csv(file.unsuppressedviralload.s))
sre <- as.data.table(read.csv(file.selfreportedart.s))

# use viral laods measurements for round 15 
sre <- sre[ROUND != 'R015']

# combine
art <- rbind(uns, sre)
art <- art[order(COMM, ROUND, SEX, AGEYRS)]

# restrict age 
art <- art[AGEYRS > 14 & AGEYRS < 50]

# check
stopifnot(nrow(art[COMM == 'inland']) == art[, length(unique(AGEYRS))] * art[, length(unique(SEX))] * art[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(art[COMM == 'fishing']) == art[, length(unique(AGEYRS))] * art[, length(unique(SEX))] * art[COMM == 'fishing', length(unique(ROUND))])

# plot
tmp <- copy(art)
tmp[, ROUND_LABEL := paste0('Round ', gsub('R0(.+)', '\\1', ROUND))]
tmp[, SEX_LABEL := 'Female']
tmp[SEX== 'M', SEX_LABEL := 'Male']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']
tmp <- tmp[!(ROUND == 'R015S' & COMM == 'inland')]
ggplot(tmp, aes(x = AGEYRS)) + 
  geom_line(aes(y = 1-PROP_UNSUPPRESSED_M, col = SEX_LABEL)) + 
  geom_ribbon(aes(ymin = 1-PROP_UNSUPPRESSED_CL, ymax = 1-PROP_UNSUPPRESSED_CU, fill = SEX_LABEL), alpha = 0.7) +
  facet_grid( ROUND_LABEL~COMM_LABEL) +
  labs(x = 'Age', y = 'ART coverage among HIV-positive participants', col = '', fill = '') +
  theme_bw() +
  scale_fill_manual(values = c('Male'='lightblue3','Female'='lightpink1')) + 
  scale_color_manual(values = c('Male'='royalblue3','Female'='deeppink')) + 
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1))) + 
  scale_y_continuous(labels = scales::percent, limits= c(0,1))
ggsave(file=file.path(outdir, paste0('smooth_artcoverage.png')), w=8, h=10)

# save
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_estimates_220906.csv'))
write.csv(art, file = file.name, row.names = F)


####################################

# COMBINE SUMMARISE POSTERIOR SAMPLE

####################################

uns <- as.data.table(read.csv(file.unsuppressedviralload))
sre <- as.data.table(read.csv(file.selfreportedart))

# use viral laods measurements for round 15 
sre <- sre[ROUND != 'R015']

# combine
art <- rbind(uns, sre)
art <- art[order(COMM, ROUND, SEX, AGEYRS)]

# restrict age 
art <- art[AGEYRS > 14 & AGEYRS < 50]

# check
stopifnot(nrow(unique(art[COMM == 'inland', .(AGEYRS, SEX, ROUND)])) == art[, length(unique(AGEYRS))] * art[, length(unique(SEX))] * art[COMM == 'inland', length(unique(ROUND))])
stopifnot(nrow(unique(art[COMM == 'fishing', .(AGEYRS, SEX, ROUND)])) == art[, length(unique(AGEYRS))] * art[, length(unique(SEX))] * art[COMM == 'fishing', length(unique(ROUND))])

# save
file.name <- file.path(indir.deepsequencedata, 'RCCS_data_estimate_incidence_inland_R6_R18/220903/', paste0('RCCS_artcoverage_posterior_samples_220906.csv'))
write.csv(art, file = file.name, row.names = F)

