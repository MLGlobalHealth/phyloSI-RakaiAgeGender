library(data.table)
library(ggplot2)


# laptop
indir <- '~/git/phyloflows'
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence.analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live'
outdir <- file.path(indir.deepsequence.analyses, 'RCCS_R15_R20', 'vl_suppofinfected_by_gender_loc_age')

# file paths
file.nonsuppressed.prop.self.reported <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', "RCCS_nonsuppressed_proportion_arvmed_220801.csv")
file.nonsuppressed.prop.vl <- file.path(indir.deepsequencedata, 'RCCS_R15_R20', "RCCS_nonsuppressed_proportion_vl_220803.csv")

#load self reported unsuppressed proportion
proportion_unsuppressed.sr <- as.data.table(read.csv(file.nonsuppressed.prop.self.reported))
proportion_unsuppressed.sr[, type := 'self-reported']
proportion_unsuppressed.sr<- proportion_unsuppressed.sr[ROUND != 'R015S']

# load viral load unsuppressed proportion
proportion_unsuppressed.vl <- as.data.table(read.csv(file.nonsuppressed.prop.vl))
proportion_unsuppressed.vl[, type := 'Viremic viral load\n(VL > 500)']

# merge
proportion_unsuppressed <- rbind(proportion_unsuppressed.sr, proportion_unsuppressed.vl)

# plot
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
ggsave(p, file=file.path(outdir, 'suppofinfected_by_gender_loc_age_vlvsselfreported.png'), w=10, h=9)

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
ggsave(p, file=file.path(outdir, 'suppofinfected_by_gender_loc_age_vl.png'), w=7, h=8)

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
ggsave(p, file=file.path(outdir, 'suppofinfected_by_gender_loc_age_vl2.png'), w=7, h=8)


