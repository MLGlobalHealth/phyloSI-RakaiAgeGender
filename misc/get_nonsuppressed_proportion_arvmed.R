library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
indir.repository <- '~/git/phyloflows'

outdir <- file.path(indir.deepsequence_analyses, 'preliminary', 'IncidentCases')

file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')
file.path.quest <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')

hiv <- as.data.table(read.csv(file.path.hiv))
quest <- as.data.table(read.csv(file.path.quest))
community.keys <- as.data.table(read.csv(file.community.keys))

file.stan.model <- file.path(indir.repository, 'misc', 'stan_models', 'one_dimensional_gp_binomial.stan')
model = rstan::stan_model(file.stan.model)
ps <- c(0.5, 0.025, 0.975)
p_labs <- c('M','CL','CU')

# keep variable of interest
rin <- quest[, .(ageyrs, round, study_id, sex, comm_num, intdate, arvmed, cuarvmed)]

# Set to date format
rin[, intdate := as.Date(intdate, format = '%d-%b-%y')]

# find  community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# remove R016 for fishing 
rinc[, sum(!is.na(arvmed)), by = c('round', 'comm')]
rinc <- rinc[!(comm == 'fishing' & round %in% c('R016')) ]

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# census eligible age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# keep HIV positive
hivs <- hivs[HIV == 'P']

# compare arvmed and cuarvmed
hivs[, sum(!is.na(ARVMED)), by = c('ROUND', 'COMM')]
hivs[, sum(!is.na(CUARVMED)), by = c('ROUND', 'COMM')]
hivs[, mean(na.omit(ARVMED == CUARVMED)), by = c('ROUND', 'COMM')]
hivs[, mean(na.omit(ARVMED == CUARVMED)), by = c('ROUND')]

# get ART status
hivs[, ART := ARVMED ==1]
hivs[is.na(ARVMED), ART := F]
hivs[, table(ROUND)]

# find non-suppresed for participant
rart <- hivs[, list(COUNT = sum(ART == F), TOTAL_COUNT = length(ART)), by = c('ROUND', 'SEX', 'COMM', 'AGEYRS')]

# fit loess smooth through 
rart <- rart[order(ROUND, COMM, SEX, AGEYRS)]
AGEYRSPREDICT = 15:49
rartfit <- rart[, {
  fit = rstan::sampling(model,data=list(N = length(AGEYRS),
                                        x = sort(AGEYRS),
                                        x_hat = AGEYRSPREDICT,
                                        N_hat = length(AGEYRSPREDICT),
                                        COUNT = COUNT,
                                        TOTAL_COUNT = TOTAL_COUNT,
                                        x_in_x_hat = which(AGEYRS %in% AGEYRSPREDICT)),
                        iter=1000,warmup=500,chains=1,
                        control = list(max_treedepth = 15, adapt_delta = 0.99))
  fit_samples = rstan::extract(fit)
  tmp1 = as.data.table( reshape2::melt(fit_samples[['mu']]))
  setnames(tmp1, 2, 'AGEYRS_INDEX')
  tmp1 = tmp1[, list(q= quantile(value, prob=ps, na.rm = T), q_label=p_labs), by=c('AGEYRS_INDEX')]
  tmp1 = dcast(tmp1, AGEYRS_INDEX ~ q_label, value.var = "q")
  tmp1[, AGEYRS := AGEYRSPREDICT[AGEYRS_INDEX]]
  
  PROP_NON_SUPPRESSED_EMPIRICAL <- vector(mode = 'numeric', length = length(AGEYRSPREDICT))
  PROP_NON_SUPPRESSED_EMPIRICAL[AGEYRSPREDICT %in% AGEYRS ] = COUNT / TOTAL_COUNT
  PROP_NON_SUPPRESSED_EMPIRICAL[!AGEYRSPREDICT %in% AGEYRS ] = NA
  
  list(M = tmp1[, M], CL = tmp1[, CL], CU = tmp1[, CU], AGEYRS = tmp1[, AGEYRS], PROP_NON_SUPPRESSED_EMPIRICAL = PROP_NON_SUPPRESSED_EMPIRICAL)
}, by = c('ROUND', 'SEX', 'COMM')]


if(0){#plot
  ggplot(rartfit, aes(x = AGEYRS)) + 
    geom_point(aes(y = PROP_NON_SUPPRESSED_EMPIRICAL)) +
    geom_line(aes(y = M, col = SEX)) + 
    geom_ribbon(aes(ymin = CL, ymax = CU, fill = SEX), alpha = 0.5) + 
    facet_grid(COMM~ROUND)
}

# no data for R016 in fishing for now use the same value than for round R016
# rartfit2 <- copy(rartfit[ROUND == 'R017' & COMM == 'fishing'])
# rartfit2[, ROUND := 'R016']
# rartfit <- rbind(rartfit, rartfit2)
# 
# # no data for R015 for now use the same value than for round R016
# rartfit2 <- copy(rartfit[ROUND == 'R016'])
# rartfit2[, ROUND := 'R015']
# rartfit <- rbind(rartfit, rartfit2)

# check
n = rartfit[, length(unique(ROUND))] * rartfit[, length(unique(AGEYRS))] * rartfit[, length(unique(COMM))] * rartfit[, length(unique(SEX))]
stopifnot(nrow(rartfit) == n)

# save
write.csv(rartfit, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', "RCCS_nonsuppressed_proportion_arvmed_220801.csv"), row.names = F)


