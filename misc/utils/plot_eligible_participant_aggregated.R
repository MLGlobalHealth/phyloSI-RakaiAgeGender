library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(ggnewscale)
#remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)

# repository directory
indir.repository <- '~/git/phyloflows'

# directory to save the figure
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/'
outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')

# files
participants <- file.path(indir.repository, 'data', 'RCCS_participation_221116.csv')
eligible <- file.path(indir.repository, 'data', 'RCCS_census_eligible_individuals_221116.csv')
file.path.round.timeline <- file.path(indir.repository, 'data', 'RCCS_round_timeline_220905.RData')

# load
rinc <- as.data.table(read.csv(participants))
ncen <- as.data.table(read.csv(eligible))
load(file.path.round.timeline)

# merge
tmp <- select(ncen, c('AGEYRS', 'ROUND', 'SEX', 'COMM', 'ELIGIBLE_NOT_SMOOTH', 'ELIGIBLE'))
rinc[, ROUND := gsub('R0(.+)', '\\1', ROUND)]
rpr <- merge(rinc,tmp , by =  c('AGEYRS', 'ROUND', 'SEX', 'COMM'))

# plot count participation rate over age pretty
tmp <- copy(rpr)
tmp[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp <- tmp[!(ROUND == '15S' & COMM == 'inland')]
tmp[, SEX_LABEL := 'Women']
tmp[SEX== 'M', SEX_LABEL := 'Men']
tmp[, COMM_LABEL := 'Fishing\n communities']
tmp[COMM == 'inland', COMM_LABEL := 'Inland\n communities']

tmp1 <- tmp[, list(PARTICIPANT = sum(PARTICIPANT), ELIGIBLE = sum(ELIGIBLE)), by = c('ROUND', 'SEX_LABEL', 'COMM_LABEL', 'COMM')]
tmp1[, NON_PARTICIPANT := ELIGIBLE - PARTICIPANT]
tmp1 <- melt.data.table(tmp1, id.vars = c('SEX_LABEL', 'COMM_LABEL', 'ROUND', 'COMM'))
tmp1 <- tmp1[variable != 'ELIGIBLE']
tmp1 <- tmp1[!ROUND %in% c("06", "07", "08", "09")]
tmp1[, ROUND_LABEL := paste0('Round ', ROUND)]
tmp1[, VARIABLE_LABEL := ifelse(variable == 'PARTICIPANT', 'Participant', 'Non-participant')]
tmp1[, COMM_LABEL := gsub('\n', '', COMM_LABEL)]

df_round <- data.table(ROUND = tmp1[, unique(ROUND)])
df_round[, ROUND_INDEX := 1:length(ROUND)]
df_round_inland[, ROUND := gsub('R0(.+)', '\\1', round)]
df_round_inland <- df_round_inland[ROUND %in% df_round$ROUND]
df_round_fishing[, ROUND := gsub('R0(.+)', '\\1', round)]
df_round_location <- rbind(df_round_inland[ROUND != '15'], df_round_fishing[ROUND %in% c('15', '15S')])
df_round_location[, midpoint_sample_date := min_sample_date + (max_sample_date - min_sample_date)/2]
set(df_round_location, NULL, 'COMM', NULL)
df_round <- merge(df_round, df_round_location, by = 'ROUND')

tmp1 <- merge(tmp1, df_round, by = 'ROUND')
tmp1 <- tmp1[order(COMM_LABEL, SEX_LABEL, variable)]
df_round_inland[, ROUND_LABEL := paste0('Round ', ROUND)]

# with round in the background
palette_round_inland <<- grDevices::colorRampPalette(c("#264653", "#2A9D8F", "#E9C46A", "#F4A261", "#E76F51"))(9)

# for readibility
df_round_inland[round == 'R012', max_sample_date := as.Date('2008-05-15')]
df_round_inland[round == 'R015', max_sample_date := as.Date('2013-06-01')]
df_round_inland[round == 'R016', max_sample_date := as.Date('2015-01-20')]

p <- ggplot(tmp1[COMM == 'inland']) +
  geom_rect(data = df_round_inland, aes(ymin = 200, ymax = 25000, xmin = min_sample_date,
                                        xmax = max_sample_date, col = ROUND_LABEL), alpha = 0.3, fill = 'white', size = 0.8) +
  scale_color_manual(values= palette_round_inland)+
  guides(color = guide_legend(byrow = T, ncol = 2, order = 4, override.aes = list(pattern = "none", col = 'white')))+
  new_scale_color() +
  geom_bar_pattern(aes(x = midpoint_sample_date, y = value, fill = SEX_LABEL, pattern = VARIABLE_LABEL), 
                   colour = 'grey60',
                   stroke = 1,
                   stat = 'identity', width=300, 
                   pattern_fill = "white",
                   pattern_colour = 'grey50',
                   pattern_angle = 45,
                   pattern_density = 0.01,
                   pattern_spacing = 0.035,
                   pattern_key_scale_factor = 0.6, 
                   size = 0.3) +
  labs(y = 'Census eligible population') +
  scale_fill_manual(values = c('Women' = 'lightpink', 'Men' = 'lightblue'))+
  scale_color_manual(values = c('Women' = 'lightpink', 'Men' = 'lightblue'))+
  scale_pattern_manual(values = c('Participant' = "none", 'Non-participant' = "stripe")) +
  theme_bw() + 
  theme(legend.position = 'bottom', 
        strip.background = element_rect(colour="white", fill="white"),
        strip.text = element_text(size = rel(1)), 
        panel.grid.major = element_blank(), 
        panel.grid.minor.x = element_blank(), 
        axis.text.x = element_text(angle = 30, hjust =1), 
        axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        legend.box = 'vertical') + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.009))) + 
  scale_x_date(expand = expansion(mult = c(0.005, 0.005))) + 
  # scale_x_continuous(breaks = tmp1[, unique(ROUND_INDEX)], label =  tmp1[, unique(ROUND_LABEL)] )  + 
  guides(alpha = guide_legend(byrow = T, nrow = 2, order = 1, override.aes = list(pattern = "none", col = 'white')),
         fill = guide_legend(byrow = T, nrow = 2, order = 2, override.aes = list(pattern = "none", col = 'white')),
         color = guide_legend(byrow = T, nrow = 2, order = 2, override.aes = list(pattern = "none", col = 'white')),
         pattern = guide_legend(byrow = T, nrow = 2, override.aes = list(fill = "white", col = 'black'), order = 3))
ggsave(p, file = file.path(outdir, 'Participants_aggregated_age_221116.pdf'), w = 3.2, h = 6.2)

