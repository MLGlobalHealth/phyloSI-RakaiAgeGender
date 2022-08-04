library(data.table)
library(ggplot2)


# laptop
indir <- '~/git/phyloflows'
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'
outdir <- file.path(indir.deepsequence.analyses, 'RCCS_R15_R20', 'vl_suppofinfected_by_gender_loc_age')

# file paths
file.nonsuppressed.prop.self.reported <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "RCCS_nonsuppressed_proportion_arvmed_220801.csv")
file.nonsuppressed.prop.vl.1000 <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0("RCCS_nonsuppressed_proportion_vl_", 1000, "_220803.csv"))
file.nonsuppressed.prop.vl.500 <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0("RCCS_nonsuppressed_proportion_vl_", 500, "_220803.csv"))
file.nonsuppressed.prop.vl.loess.1000 <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', paste0("RCCS_nonsuppressed_proportion_loess_vl_", 1000, "_220803.csv"))

#load self reported unsuppressed proportion
proportion_unsuppressed.sr <- as.data.table(read.csv(file.nonsuppressed.prop.self.reported))
proportion_unsuppressed.sr[, type := 'self-reported']
proportion_unsuppressed.sr<- proportion_unsuppressed.sr[ROUND != 'R015S']
proportion_unsuppressed.sr[, method := 'GP']

# load viral load unsuppressed proportion with viremic viral load 500 GP
proportion_unsuppressed.vl.500 <- as.data.table(read.csv(file.nonsuppressed.prop.vl.500))
proportion_unsuppressed.vl.500[, type := paste0('Viremic viral load\n(VL > ', '500', ')')]
proportion_unsuppressed.vl.500[, method := 'GP']

# load viral load unsuppressed proportion with viremic viral load 1000 GP
proportion_unsuppressed.vl.1000 <- as.data.table(read.csv(file.nonsuppressed.prop.vl.1000))
proportion_unsuppressed.vl.1000[, type := paste0('Viremic viral load\n(VL > ', '1000', ')')]
proportion_unsuppressed.vl.1000[, method := 'GP']

# load viral load unsuppressed proportion with viremic viral load 1000 loess
proportion_unsuppressed.loess.vl.1000 <- as.data.table(read.csv(file.nonsuppressed.prop.vl.loess.1000))
proportion_unsuppressed.loess.vl.1000 <- melt.data.table(proportion_unsuppressed.loess.vl.1000, id.vars = c('ROUND', 'COMM', 'SEX', 'AGEYRS', 'N', 'HIV_N', 'VLNS_N', 'ARV_N', 'AGE', 
                                                                   'HIV_AND_VLNS', 'ROW_ID', 'EMPIRICAL_VLNS_IN_HIV', 'PROP_NON_SUPPRESSED_EMPIRICAL'))
setnames(proportion_unsuppressed.loess.vl.1000, c('variable', 'value'), c('method', 'M'))
proportion_unsuppressed.loess.vl.1000[, type := paste0('Viremic viral load\n(VL > ', '1000', ')')]

# merge
proportion_unsuppressed <- do.call('rbind', list(proportion_unsuppressed.sr, proportion_unsuppressed.vl.500, proportion_unsuppressed.vl.1000, proportion_unsuppressed.loess.vl.1000, fill=TRUE))

#
# compare viral viremic viral load 500 and 1000 with GP
p <- ggplot(proportion_unsuppressed[!grepl('self', type) & method == 'GP']) + 		
  geom_ribbon(aes(x=AGEYRS, ymin=CL, ymax=CU, fill = type), alpha=0.2) +
  geom_line(aes(x=AGEYRS, y=M, colour=type)) +
  geom_point(aes(x=AGEYRS, y=PROP_NON_SUPPRESSED_EMPIRICAL, colour=type), alpha = 0.5) +
  scale_x_continuous( expand=c(0,0) ) + 
  scale_y_continuous(label=scales:::percent) +
  facet_grid(ROUND~SEX+COMM) +
  theme_bw() +
  ggsci::scale_color_jama() + 
  ggsci::scale_fill_jama() + 
  labs(x='\nage at visit (years)', 
       y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n')
ggsave(p, file=file.path(outdir, paste0('suppofinfected_by_gender_loc_age_vl500_vs_vl1000.png')), w=10, h=9)

#
# compare GP and loess smoth with viral viremic viral load 1000 with GP
VIREMIC_VIRAL_LOAD <- 1000
tmp <- proportion_unsuppressed[!grepl(paste0('self|500'), type) & !grepl('50', method)]
p <- ggplot(tmp) + 		
  geom_line(aes(x=AGEYRS, y=M, colour=method)) +
  geom_point(data = tmp[method == 'GP'], aes(x=AGEYRS, y=PROP_NON_SUPPRESSED_EMPIRICAL), alpha = 0.5, color = 'grey') +
  scale_x_continuous( expand=c(0,0) ) + 
  scale_y_continuous(label=scales:::percent) +
  facet_grid(ROUND~SEX+COMM) +
  theme_bw() +
  ggsci::scale_color_jama() + 
  ggsci::scale_fill_jama() + 
  labs(x='\nage at visit (years)', 
       y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n')
ggsave(p, file=file.path(outdir, paste0('suppofinfected_by_gender_loc_age_loess_vs_gp_vl_', VIREMIC_VIRAL_LOAD, '.png')), w=10, h=9)


# plot
proportion_unsuppressed <- proportion_unsuppressed[grepl(paste0('self|', VIREMIC_VIRAL_LOAD), type)]
p <- ggplot(proportion_unsuppressed) + 		
  geom_ribbon(aes(x=AGEYRS, ymin=CL, ymax=CU, fill = type), alpha=0.2) +
  geom_line(aes(x=AGEYRS, y=M, colour=type)) +
  geom_point(aes(x=AGEYRS, y=PROP_NON_SUPPRESSED_EMPIRICAL, colour=type), alpha = 0.5) +
  scale_x_continuous( expand=c(0,0) ) + 
  scale_y_continuous(label=scales:::percent) +
  facet_grid(ROUND~SEX+COMM) +
  theme_bw() +
  labs(x='\nage at visit (years)', 
       y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n')
ggsave(p, file=file.path(outdir, paste0('suppofinfected_by_gender_loc_age_vlvsselfreported_', VIREMIC_VIRAL_LOAD, '.png')), w=10, h=9)

# plot vl col gender
p <- ggplot(proportion_unsuppressed.vl) + 		
  geom_ribbon(aes(x=AGEYRS, ymin=CL, ymax=CU, group = SEX), alpha=0.2) +
  geom_line(aes(x=AGEYRS, y=M, colour=SEX)) +
  geom_point(aes(x=AGEYRS, y=PROP_NON_SUPPRESSED_EMPIRICAL, colour=SEX), alpha = 0.5) +
  scale_x_continuous( expand=c(0,0) ) + 
  scale_y_continuous(label=scales:::percent) +
  facet_grid(ROUND~COMM) +
  scale_colour_manual(values=c('M'='royalblue3','F'='deeppink2')) +
  theme_bw() +
  labs(x='\nage at visit (years)', 
       y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n')
ggsave(p, file=file.path(outdir, paste0('suppofinfected_by_gender_loc_age_vl_', VIREMIC_VIRAL_LOAD, '.png')), w=7, h=8)

# plot vl col round
p <- ggplot(proportion_unsuppressed.vl) + 		
  geom_ribbon(aes(x=AGEYRS, ymin=CL, ymax=CU, fill = ROUND), alpha=0.2) +
  geom_line(aes(x=AGEYRS, y=M, colour=ROUND)) +
  geom_point(aes(x=AGEYRS, y=PROP_NON_SUPPRESSED_EMPIRICAL, colour=ROUND), alpha = 0.5) +
  scale_x_continuous( expand=c(0,0) ) + 
  scale_y_continuous(label=scales:::percent) +
  facet_grid(SEX~COMM) +
  theme_bw() +
  labs(x='\nage at visit (years)', 
       y='HIV+ individuals with unsuppressed viral load\n(95% credibility interval)\n')
ggsave(p, file=file.path(outdir, paste0('suppofinfected_by_gender_loc_age_vl_', VIREMIC_VIRAL_LOAD, '_2.png')), w=7, h=8)


