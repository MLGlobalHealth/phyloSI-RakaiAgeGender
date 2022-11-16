library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggpubr)

# directory repository
indir.repository <- '~/git/phyloflows'

# directory to save the figure
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

# files
eligible <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')
participation <- file.path(indir.repository, 'data', 'RCCS_participation_221116.csv')

# load
ncen <- as.data.table(read.csv(eligible))
rinc <- as.data.table(read.csv(participation))

#  plot census eligible 
tmp <- copy(ncen)
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09") & COMM == 'inland' & SEX== 'M']
pE.M <- ggplot(tmp1, aes(x = AGEYRS)) +
  geom_bar(aes(y = ELIGIBLE_NOT_SMOOTH, fill = SEX_LABEL), stat = 'identity', alpha = 0.7) +
  geom_line(data = select(tmp1[ROUND == '10'], -'ROUND_LABEL'), aes(y = ELIGIBLE, linetype = 'Round 10'), color = 'paleturquoise4') +
  geom_line(aes(y = ELIGIBLE, col  = SEX_LABEL)) +
  labs(y = 'Census eligible individuals', x = 'Age', col = '', fill = '', linetype = '') +
  facet_grid(ROUND_LABEL~SEX_LABEL) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 12), 
        legend.text = element_text(size=12),
        strip.text.y = element_blank()) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 800))  

tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09") & COMM == 'inland' & SEX== 'F']
pE.F <- ggplot(tmp1, aes(x = AGEYRS)) +
  geom_bar(aes(y = ELIGIBLE_NOT_SMOOTH, fill = SEX_LABEL), stat = 'identity') +
  geom_line(data = select(tmp1[ROUND == '10'], -'ROUND_LABEL'), aes(y = ELIGIBLE, linetype = 'Round 10'), color = 'pink4') +
  geom_line(aes(y = ELIGIBLE, col  = SEX_LABEL)) +
  labs(y = 'Census eligible individuals', x = 'Age', col = '', fill = '', linetype = '') +
  facet_grid(ROUND_LABEL~SEX_LABEL) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 12), 
        legend.text = element_text(size=12),
        strip.text.y = element_blank()) + 
  scale_x_continuous(expand = c(0,0)) + 
  scale_y_continuous(expand = expansion(mult = c(0,0.1)), limits = c(0, 800))  


# plot participation
tmp <- select(ncen, c('AGEYRS', 'ROUND', 'SEX', 'COMM', 'ELIGIBLE_NOT_SMOOTH', 'ELIGIBLE'))
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rpr <- merge(rinc,tmp , by =  c('AGEYRS', 'ROUND', 'SEX', 'COMM'))
rpr[, PARTICIPATION_NOT_SMOOTH := PARTICIPANT / ELIGIBLE_NOT_SMOOTH]
rpr[PARTICIPATION_NOT_SMOOTH > 1, PARTICIPATION_NOT_SMOOTH := 1]

tmp <- copy(rpr)
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

# plot participation rate
tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09") & COMM=='inland' & SEX == "M"]
pP.M <- ggplot(tmp1, aes(x = AGEYRS)) +
  geom_bar(aes(y = PARTICIPATION_NOT_SMOOTH, fill = SEX_LABEL), stat = 'identity') +
  geom_line(data = select(tmp1[ROUND == '10'], -'ROUND_LABEL'), aes(y = PARTICIPATION, linetype = 'Round 10'), color = 'paleturquoise4') +
  geom_line(aes(y = PARTICIPATION, col = SEX_LABEL)) +
  labs(y = 'RCCS participation rate', x = 'Age', col = '', fill = '', linetype ='') +
  facet_grid(ROUND_LABEL~SEX_LABEL) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = 12), 
        legend.text = element_text(size=12),
        panel.spacing.y = unit(0.7, "lines")) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))

tmp1 <- tmp[!ROUND %in% c("06", "07", "08", "09") & COMM=='inland' & SEX == "F"]
pP.F <- ggplot(tmp1, aes(x = AGEYRS)) +
  geom_bar(aes(y = PARTICIPATION_NOT_SMOOTH, fill = SEX_LABEL), stat = 'identity') +
  geom_line(data = select(tmp1[ROUND == '10'], -'ROUND_LABEL'), aes(y = PARTICIPATION, linetype = 'Round 10'), color = 'pink4') +
  geom_line(aes(y = PARTICIPATION, col = SEX_LABEL)) +
  labs(y = 'RCCS participation rate', x = 'Age', col = '', fill = '', linetype ='') +
  facet_grid(ROUND_LABEL~SEX_LABEL) +
  scale_fill_manual(values = c('Men'='lightblue3','Women'='lightpink1')) + 
  scale_color_manual(values = c('Men'='royalblue3','Women'='deeppink')) +
  theme_bw() +
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_blank(), 
        legend.text = element_text(size=12),
        panel.spacing.y = unit(0.7, "lines")) + 
  scale_y_continuous(labels = scales::percent, limits = c(0,1), expand = c(0,0)) + 
  scale_x_continuous(expand = c(0,0))

# combine figure
p <- ggarrange(pE.F, pE.M, pP.F, pP.M, legend = 'bottom', common.legend = T, nrow = 1, labels = c('a', 'b', 'c', 'd'), 
               widths = c(rep(0.22, 3),0.25))
ggsave(p, file = file.path(outdir, 'CensusEligibleCount_Participation_221116.pdf'), w = 9, h = 12)

