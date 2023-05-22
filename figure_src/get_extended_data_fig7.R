library(data.table)
library(ggplot2)

jobname <- 'central3'
stan_model <- 'gp_221201d'

# indir repository
indir <- here::here() # please modify accordingly
outdir <- here::here()
gitdir <- indir
source( file.path( indir, 'config.R' ) )  # loads the paths

outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
if (!dir.exists((outfile))) dir.create((outfile))
if (!dir.exists(file.path(outdir, 'figures'))) dir.create(file.path(outdir, 'figures'))
if (!dir.exists((outfile.figures))) dir.create((outfile.figures))

# paths to data
# file.unsuppressed.share <- file.path(indir, 'fit', paste0('RCCS_unsuppressed_share_sex_221208.csv')) # please modify accordingly
# file.path.round.timeline <- file.path(indir, 'data', 'RCCS_round_timeline_220905.RData') # please modify accordingly
# file.expected_contribution_age_source <- file.path(indir, 'data', paste0(stan_model, '-', jobname, '-output-log_lambda_latentby_direction_round_age_transmission.sourcestandardisedby_direction_round.rds')) # please modify accordingly
# file.contribution.sexual.contacts <- file.path(indir, 'data', 'inland-R015_age-dist_ma_cntct_area_1549.rds')
# load functions
source(file.path(indir, 'R', 'functions_transmission_flow', 'postprocessing_summary_functions.R'))
source(file.path(indir, 'R', 'functions_transmission_flow', 'postprocessing_plot_functions.R'))
source(file.path(indir, 'R', 'functions_transmission_flow', 'postprocessing_utils_functions.R'))
source(file.path(indir, 'R', 'functions_transmission_flow', 'postprocessing_statistics_functions.R'))
source(file.path(indir, 'R', 'functions_transmission_flow', 'summary_functions.R'))

# get fig A ----
# map direction
df_direction <- get.df.direction()
df_community <- get.df.community()
load(file.path.round.timeline)
df_round <- make.df.round(df_round_inland)

# load data
unsuppressed_share <- fread(file.unsuppressed.share.fig) # share of unsuppressed count by sex
expected_contribution_age_source <- as.data.table(readRDS(file.expected_contribution_age_source))
# load the sexual contact data
contribution_sexual_contact <- as.data.table(readRDS(file.age_dist_ma_cntct_area))

# clean
df_contribution_sexual_contact <- clean_contribution_sexual_contact_yu(contribution_sexual_contact)
unsuppressed_share_age <- prepare_unsuppressed_share(unsuppressed_share, c('SEX', 'AGEYRS'), 'SEX')

# plot A----
# contribution to incidence
cti <- copy(expected_contribution_age_source)

# contribution to sexual contact
cts <- copy(df_contribution_sexual_contact)

# contribution to individual with unsuppressed virus
ctu <- copy(unsuppressed_share_age)

# select one round
Round <- 'R018'
cti <- cti[ROUND == Round]
ctu <- ctu[ROUND == Round]

label.cont.inc <- 'Contribution to HIV incidence'
label.cont.uns <- 'Contribution to individuals with unsuppressed virus'
label.cont.cont <- 'Contribution to sexual contacts'
# plot
communities <- cti[, unique(COMM)]
Comm <- communities

cti.c <- cti[COMM == Comm]
ctu.c <- ctu[COMM == Comm]
setnames(cti.c, 'AGE_TRANSMISSION.SOURCE', 'AGEYRS')
setnames(ctu.c, 'LABEL_GENDER_SOURCE', 'SEX')
ctu.c[, label := 'Age contribution with unsuppressed HIV']
setnames(cti.c, 'LABEL_GENDER_SOURCE', 'SEX')
cti.c[, gender := ifelse(LABEL_GENDER_RECIPIENT == 'Women','male', 'female')]
ctu.c[, gender2 := ifelse(LABEL_GENDER_RECIPIENT == 'Women', 'men', 'women')]
cts[, gender := ifelse(LABEL_GENDER_RECIPIENT == 'Women', 'male', 'female')]
ctu.c[, gender := ifelse(LABEL_GENDER_RECIPIENT == 'Women','male', 'female')]


cti.c[, label := paste0('Contribution to ', gender ,  '\ntransmitting partners   ' )]
ctu.c[, label := paste0('Contribution to ', gender2, '\nwith unsuppressed HIV   ')]
cts[, label := paste0('Contribution to ', gender, '\nsexual contacts   ')]

pltA <- rbind(cti.c[, list(AGEYRS,M,CL,CU,LABEL_SOURCE, label, gender)],
              ctu.c[, list(AGEYRS,M,CL,CU,LABEL_SOURCE, label, gender)],
              cts[, list(AGEYRS,M,CL,CU,LABEL_SOURCE, label, gender)])
# plt
setkey(pltA, LABEL_SOURCE)
pltA[, facet_title := ifelse(gender == 'female', 'Women', 'Men')]
pA <- ggplot(pltA, aes(x = AGEYRS, y = M)) +
  geom_line(aes(y = M, col = factor(label, levels = unique(label)), linetype = factor(label, levels = unique(label)))) +
  geom_ribbon(aes(ymin = CL, ymax = CU, fill = factor(label, levels = unique(label))), alpha = .3) +
  facet_grid(.~factor(facet_title, levels = c('Women', 'Men'))) +
  scale_color_manual(values = c("#ae017e",'grey50', '#cc4c02',"lightblue4",'grey50' ,'#41ab5d')) +
  scale_fill_manual(values = c("#fa9fb5",'grey70','#ec7014',"lightblue3",'grey50', '#78c679' )) +
  scale_x_continuous(breaks = seq(15, 50, 5)) +
  theme_bw() +
  labs(x = 'Female age                                                                                                           Male age',
       y = 'Percent') +
  theme(
    panel.grid.major = element_line(linewidth = 0.3),

    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.2),
    strip.background = element_blank(),
    axis.text = element_text(size = 5, family = 'sans'),
    legend.text =  element_text(size = 5, family = 'sans'),
    strip.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank(),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.title = element_text(size = 7, family = 'sans')
  ) +
  theme(
    legend.position = c(0.69, 1.05),
    legend.justification = c("right", "top"),
    legend.box.just = "right",
    legend.background = element_blank()) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0, 0.08)) +
  scale_linetype_manual(values = c(1,1, 2, 1, 1, 2)) +
  guides(color = guide_legend(nrow = 3, byrow = F))

# conditional contact intensities----
pd.a <- as.data.table(readRDS(file.age_dist_cntct_area))

pd.a <- pd.a[cont.age %in% c(15, 20, 25, 30, 35, 40)]
pd.a$part.sex <- ifelse(pd.a$part.sex == 'F', 'Women', 'Men')
pd.a$cont.sex <- ifelse(pd.a$part.sex == 'Women',  'Men', 'Women')
pd.a <- pd.a[, label := ifelse(cont.age %in% c(15, 20, 25), '15, 20, 25 years', '30, 35, 40 years')]

pd.a[, plt.bar := (part.age - cont.age) %in% c( -10, -5, 0,  5,  15)]
pd.a[, plt.age := cont.age]
pd.a[, cont.age := paste0('Sexual contact intensities to ', cont.sex, ' aged ', cont.age)]
pd.a[, type := 'Contribution to sexual contacts']

# load data----
dta <- readRDS(file.expected_contribution_sliced_age_source)
dta <- dta[ROUND == 'R018'] # subset to gender of the source (i.e.., Female sources or Male sources)
setnames(dta, c('LABEL_GENDER_SOURCE', 'LABEL_GENDER_RECIPIENT', 'AGE_TRANSMISSION.SOURCE', 'AGE_INFECTION.RECIPIENT'), c('part.sex', 'cont.sex', 'part.age', 'cont.age'))
dta <- dta[cont.age %in% c(15, 20, 25, 30, 35, 40)]
dta <- dta[, label := ifelse(cont.age %in% c(15, 20, 25), '15, 20, 25 years', '30, 35, 40 years')]
dta[, plt.age := cont.age]

dta[, plt.bar := (part.age - cont.age) %in% c( -9, -4, 1,  6,  16)]
dta[, type := 'Contribution to transmitting partners ']
dta <- dta[, cont.age := ifelse(cont.age %in% c(15, 20, 25), paste0('Sources trainsmitting to ', cont.sex, ' aged ', cont.age), paste0('Sources trainsmitting to ', cont.sex, ' aged ', cont.age))]

plt <- rbind(pd.a[, list(part.age,cont.age,CL,CU,M,part.sex,cont.sex,label, plt.bar, plt.age,type)],
             dta[, list(part.age,cont.age,CL,CU,M,part.sex,cont.sex,label, plt.bar, plt.age,type)]
             , use.names = TRUE, fill = TRUE)
plt <- plt[order(label)]
plt$cont.age <- as.factor(plt$cont.age)

# plot B new----
plt[,label := paste0(cont.sex, ' aged ', label)]
# try to get all labels
# if part.sex is Men, then change women's age 15 -> 30 etc
plt[cont.sex == 'Men', plt.age := ifelse(plt.age %in% c(15, 20, 25), plt.age + 15, plt.age)]
# if part.sex is Women, then change men's age 30 -> 15 etc
plt.xaxis <- ifelse(plt$part.sex == 'Men', 'Male', 'Female')
p1b <- ggplot(plt[grepl('15', label)], aes(x = part.age, y = M, col = factor(plt.age, levels = unique(plt.age)), linetype = type)) +
  geom_line(data =  plt[grepl('15', label)], aes(x = part.age, y = M)) +
  #
  geom_pointrange(data = plt[grepl('15', label) &  (plt.bar == TRUE )],
                  aes(ymin = CL, ymax = CU), linetype = 1, size = .5, fatten = 2) +

  facet_grid(.~paste0(label)) +
  scale_x_continuous(breaks = seq(15, 70, 5)) +
  scale_color_manual(values = c(
    '#f46d43',
    '#3288bd',
    '#88419d',
    '#f46d43',
    '#3288bd',
    '#88419d'
  )) +

  scale_linetype_manual(values = c(2,1)) +
  labs(x = 'Female age                                                                                                           Male age',
       y = paste0('Percent'),
       col = 'Age of contact',
       linetype = 'type of contact'
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size = 0.3),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.3),
    axis.ticks = element_line(size = 0.2),
    strip.background = element_blank(),
    axis.text = element_text(size = 5, family = 'sans'),
    legend.text =  element_text(size = 5, family = 'sans'),
    strip.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank(),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.title = element_text(size = 7, family = 'sans')
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))

p2b <- ggplot(plt[grepl('35', label)], aes(x = part.age, y = M, col = factor(plt.age, levels = unique(plt.age)), linetype = type)) +
  geom_line(data =  plt[grepl('35', label)], aes(x = part.age, y = M)) +
  geom_pointrange(data = plt[grepl('35', label) &  (plt.bar == TRUE )],
                  aes(ymin = CL, ymax = CU), linetype = 1, size = .5, fatten = 2) +
  facet_grid(.~paste0(label)) +
  scale_x_continuous(breaks = seq(15, 70, 5)) +
  scale_color_manual(values = c(
    '#f46d43',
    '#3288bd',
    '#88419d',
    '#f46d43',
    '#3288bd',
    '#88419d'
  )) +

  scale_linetype_manual(values = c(2,1)) +
  labs(x = 'Female age                                                                                                           Male age',
       y = paste0('Percent'),
       col = 'Age of contact',
       linetype = 'type of contact'
  ) +
  guides(color = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(
    legend.position = 'bottom',
    panel.grid.minor = element_blank(),
    panel.background = element_rect('white'),
    panel.grid.major = element_line(size = 0.3),
    panel.border = element_rect(fill = NA, colour = "black", size = 0.3),
    axis.ticks = element_line(size = 0.2),
    strip.background = element_blank(),
    axis.text = element_text(size = 5, family = 'sans'),
    legend.text =  element_text(size = 5, family = 'sans'),
    strip.text = element_text(size = 7, family = 'sans'),
    legend.title = element_blank(),
    axis.line = element_line(colour = 'black', size = 0.2),
    axis.title = element_text(size = 7, family = 'sans')
  ) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, .05)), limits = c(0,NA))

pB <- ggpubr::ggarrange(p1b, p2b, ncol = 1, common.legend = T, legend = 'bottom')

# combine two plots----
p <- ggpubr::ggarrange(pA, pB, ncol = 1, heights  = c(1.35, 2.8), labels = c('a', 'b'), font.label = list(color = "black", size = 8))
ggsave(file = file.path(outfile.figures, 'extended-data-fig_age-dist_EDF7.pdf'), p, width = 18, height = 20, units = 'cm', dpi = 310, limitsize = FALSE)
