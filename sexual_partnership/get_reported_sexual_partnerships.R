library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)

.year.diff <- function(x, y)
{
  if (!is.Date(x)) {x <- as.Date(x, format = '%Y-%m-%d')}
  if (!is.Date(y)) {y <- as.Date(y, format = '%Y-%m-%d')}
  lubridate::time_length(difftime(x, y),"years")
}

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'

file.meta.R1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.census.eligible.individuals.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220411.csv')

community.keys <- as.data.table(read.csv(file.community.keys))
meta <- as.data.table(read.csv(file.meta.R1518))
census_eligible_count <- as.data.table(read.csv(file.census.eligible.individuals.count))

code_na = 97:99

#####################################
# find self-reported age preference #
#####################################

meta[, pt_id := paste0('RK-', study_id)]
meta[, visit_date := as.Date(intdate, '%d-%b-%y')]
meta[, birth_date := as.Date(paste0(gsub('.{2}$', '', birthdat), birthyr), '%d-%b-%Y')]
meta[is.na(birth_date) | nchar(birthyr) == 2 | birthyr < 1900, birth_date := visit_date - ageyrs*365]

# clean
tmp <- meta[, .(pt_id, birth_date, visit_date, sex, comm_num, round,
                days1, weeks1, months1, years1,
                days2, weeks2, months2, years2,
                days3, weeks3, months3, years3,
                days4, weeks4, months4, years4,
                sexyear, eversex, sexp1yr,
                rldyslt1, rlwkslt1, rlmoslt1, rlyrslt1,
                rldyslt2, rlwkslt2, rlmoslt2, rlyrslt2,
                rldyslt3, rlwkslt3, rlmoslt3, rlyrslt3,
                rldyslt4, rlwkslt4, rlmoslt4, rlyrslt4,
                rltnage1, rltnage2, rltnage3, rltnage4,
                rltnyrs1, rltnyrs2, rltnyrs3, rltnyrs4)]
tmp <- unique(tmp)

tmp <- as.data.table( reshape2::melt(tmp, id.vars = c('pt_id', 'birth_date', 'visit_date', 'sex', 'comm_num', 'round', 'sexyear', 'sexp1yr', 'eversex')) )
tmp[, relation := as.numeric(sub('.*(?=.$)', '', variable, perl=T))]
tmp[, variable := gsub('.{2}$', '', variable)]
age_preference <- as.data.table(reshape2::dcast(tmp, pt_id + birth_date + visit_date + sex + comm_num + round + sexyear + sexp1yr + eversex + relation ~ variable, value.var = 'value'))

# find start and end date of the relationship
age_preference[, days_since_start := 0]
age_preference[!day %in% code_na, days_since_start := days_since_start + day]
age_preference[!week %in% code_na, days_since_start := days_since_start + week*7]
age_preference[!month %in% code_na, days_since_start := days_since_start + month*31]
age_preference[!year %in% code_na, days_since_start := days_since_start + year*365]
age_preference[day %in% code_na & week %in% code_na & month %in% code_na & year %in% code_na, days_since_start := NA]
age_preference[, date_start := visit_date - days_since_start]
age_preference[, age_start := .year.diff(date_start, birth_date) ]

age_preference[, days_since_stop := 0]
age_preference[!rldysl %in% code_na, days_since_stop := days_since_stop + rldysl]
age_preference[!rlwksl %in% code_na, days_since_stop := days_since_stop + rlwksl*7]
age_preference[!rlmosl %in% code_na, days_since_stop := days_since_stop + rlmosl*31]
age_preference[!rlyrsl %in% code_na, days_since_stop := days_since_stop + rlyrsl*365]
age_preference[rldysl %in% code_na & rlwksl %in% code_na & rlmosl %in% code_na & rlyrsl %in% code_na, days_since_stop := NA]
age_preference[, date_stop := visit_date - days_since_stop]
age_preference[, age_stop := .year.diff(date_stop, birth_date) ]
age_preference <- age_preference[(age_stop > age_start) | is.na(age_stop) | is.na(age_start)]

# create time since end relationship
age_preference[, time_since_date_stop := .year.diff(visit_date, date_stop)]

# find age index
age_preference[, age_index := .year.diff(visit_date, birth_date)]

# find age partner
age_preference[, age_partner_classification := NA_character_]
age_preference[rltnag == 1, age_partner_classification := 'older']
age_preference[rltnag == 2, age_partner_classification := 'younger']
age_preference[rltnag == 3, age_partner_classification := 'same']

age_preference[!rltnyr %in% code_na & age_partner_classification == 'younger', age_partner := age_index - rltnyr]
age_preference[!rltnyr %in% code_na & age_partner_classification == 'older', age_partner := age_index + rltnyr]
age_preference[age_partner_classification == 'same', age_partner := age_index ]

# check that individual do not appear more than 4 times in each round (for the 4 relationshups)
tmp <- age_preference[, list(count = .N), by = c('pt_id', 'round')]
stopifnot(all(tmp[, count <= 4]))

# total number of relationships reported in the last 12 months
ra <- merge(age_preference, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
ra[, part.comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
setnames(ra, 'sex', 'part.sex')
ra[, part.age := floor(age_index)]
ra <- ra[sexyear == 1, list(count = length(unique(pt_id))), by = c('sexp1yr', 'round', 'part.comm', 'part.age', 'part.sex')]
tmp <- ra[round == 'R015']
tmp[sexp1yr == 93, sum(count)]
tmp[sexp1yr == 0, sum(count)]

# remove missing variables
ragem <- age_preference[sexyear %in% 1:2] # remove participant who did not say wether or not they had sexual intercourse in the last year
ragem <- ragem[!(sexyear == 1 & is.na(date_stop))] # remove participant who had sex in the last year and who did not report the stop data of any of their relations
ragem <- ragem[!(sexyear == 1 & !rltnag %in% 1:3)]  # remove participant who ever had sex in the last year and who did not report the age clasification of any of their relations
ragem <- ragem[!(sexyear == 1 & rltnag %in% 1:2 & rltnyr %in% code_na)]  # remove participant who ever had sex in the last year and who did not report the age diff of any of their relations (if the classification was not same)
# ragem <- ragem[!(is.na(rltnag) & is.na(rltnyr))]  # remove participant who were not asked age of any of their relations
ragem <- ragem[!(round %in% c('R016', 'R017', 'R018'))] # same as above

# change coding of sexyear if one relationship ended within a year
ragem[sexyear == 2 & time_since_date_stop <= 1, sexyear := 1]

# keep age index within census eligible age
ragem <- ragem[!is.na(age_index)]
ragem <- ragem[age_index >= 15 & age_index < 50]
ragem <- ragem[!(sexyear == 1 &  is.na(age_partner))]
ragem <- ragem[!(sexyear == 1 & (age_partner < 15 | age_partner >= 50))]

# plot
if(0){
  tmp <- unique(age_preference[, .(pt_id, round, sexyear)])
  tmp <- tmp[, list(YES = sum(sexyear == 1), 
                               NO = sum(sexyear == 2), 
                               NR = sum(sexyear %in% 8:9), 
                               N = length(unique(pt_id))), by = c('round')]
  tmp <- as.data.table(reshape2::melt(tmp, id.vars = c('N', 'round')))
  tmp[, round_label := paste0('Round ', round, '\nN = ', N)]
  ggplot(tmp, aes(x = variable, y = value)) + 
    geom_bar(stat = 'identity') + 
    facet_wrap(~round_label, nrow = 1) +
    labs(x = "Did the participant had sexual intercouse in the last 12 months", y = 'participants count') + 
    ggtitle('All participants')
  ggsave('~/Downloads/sexyear.png', w = 9, h = 5)
  
  #####
  
  tmp <- ra[round == 'R015']
  tmp[, part.age.group := '35-49']
  tmp[part.age < 35, part.age.group := '25-34']
  tmp[part.age < 25, part.age.group := '15-24']
  
  tmp1 <- tmp[, list(count = sum(count)), by = c('sexp1yr', 'part.age.group', 'part.sex')]
  tmp1[sexp1yr == 93, sum(count)]
  tmp1[sexp1yr == 0, sum(count)]
  tmp1 <- tmp1[sexp1yr != 93]
  ggplot(tmp1, aes(x = sexp1yr, y = count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid(part.sex~part.age.group, scale = 'free_y') + 
    labs(x='Number of different sexual partners in the last 12 months', y = 'count participants') 
  ggsave('~/Downloads/number_sexual_partners_12months_long.png', w = 9, h = 5)
  
  ggplot(tmp1, aes(x = sexp1yr, y = count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid(part.sex~part.age.group, scale = 'free_y') + 
    labs(x='Number of different sexual partners in the last 12 months', y = 'count participants')+
    coord_cartesian(xlim = c(0, 15))
  ggsave('~/Downloads/number_sexual_partners_12months_short.png', w = 9, h = 5)
  
  tmp <- age_preference[round == 'R015']
  tmp[, part.age := floor(age_index)]
  tmp[, missing_values := !(sexyear %in% 1:2)] # remove participant who did not say wether or not they had sexual intercourse in the last year
  tmp[missing_values == F, missing_values := sexyear == 1 & is.na(date_stop)] # remove participant who had sex in the last year and who did not report the stop data of any of their relations
  tmp[missing_values == F, missing_values := sexyear == 1 & !rltnag %in% 1:3]  # remove participant who ever had sex in the last year and who did not report the age clasification of any of their relations
  tmp[missing_values == F, missing_values := sexyear == 1 & rltnag %in% 1:2 & rltnyr %in% code_na] 
  tmp[sexyear == 2, sexp1yr := 0]
  tmp1 <- tmp[, missing_values := !(sum(!missing_values) >= unique(sexp1yr)), by = c('pt_id', 'part.age', 'sex')]
  tmp1 <- tmp1[, list(count = length(unique(pt_id))), by = c('missing_values', 'part.age', 'sex')]
  tmp1[, missing_values_label := 'Missing values in age partner questions']
  tmp1[missing_values == F,  missing_values_label := 'No missing values in age partner questions']
  tmp1[, total_count := sum(count), by = c('part.age', 'sex')]
  tmp1[, proportion_count := count / total_count]
  ggplot(tmp1, aes(x = part.age, y = count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid(missing_values_label~sex, scale = 'free') + 
    labs(x='Age participant', y = 'count participants')
  ggsave('~/Downloads/missing_values_age_sex.png', w = 9, h = 5)
  ggplot(tmp1[missing_values == T], aes(x = part.age, y = proportion_count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid( .~sex, scale = 'free') + 
    labs(x='Age participant', y = 'Proprotion of participants with missing values')
  ggsave('~/Downloads/prop_missing_values_age_sex.png', w = 9, h = 5)
  
  #####
  
  tmp <- unique(ragem[, .(pt_id, round, sexyear)])
  tmp <- tmp[, list(YES = sum(sexyear == 1), 
                    NO = sum(sexyear == 2), 
                    N = length(unique(pt_id))), by = c('round')]
  tmp <- as.data.table(reshape2::melt(tmp, id.vars = c('N', 'round')))
  tmp[, round_label := paste0('Round ', round, '\nN = ', N)]
  ggplot(tmp, aes(x = variable, y = value)) + 
    geom_bar(stat = 'identity') + 
    facet_wrap(~round_label, nrow = 1) +
    labs(x = "Did the participant had sexual intercouse in the last 12 months", y = 'participants count') + 
    ggtitle('Participants without missing data')
  ggsave('~/Downloads/sexyear_womissing.png', w = 5, h = 5)
  
  ggplot(ragem, aes(x = time_since_date_stop)) + 
    geom_histogram(bins = 50) + 
    facet_grid(round~., scale = 'free') + 
    scale_y_log10() + 
    geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkred') + 
    theme_bw() + 
    labs(x = 'Time since end of the relationghip (years)', y = 'relationships count')
  ggsave('~/Downloads/count_time_since_date_stop.png', w = 6, h = 5)
}

# keep relationship that happened within a year for individuals who had sex within a year
rage <- ragem[!(sexyear == 1 & time_since_date_stop > 1)]

# keep one relationship (for one count) for individuals who did not have sex within a year
tmp <- rage[sexyear == 2, list(relation = min(relation)), by = c('pt_id', 'round')]
tmp <- merge(tmp, rage[sexyear == 2], by = c('pt_id', 'round', 'relation'))
tmp[, age_partner := NA]
rage <- rbind(rage[sexyear == 1], tmp)


###############
# FORMAT DATA #
###############

# create community variable
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rage <- merge(rage, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
rage[, `Community index` := comm]

# aggregate by 1 year age band
rage[, part.age := floor(age_index)]
rage[, cont.age := floor(age_partner)]

# rename comm and round
setnames(rage, c('comm', 'round'), c('part.comm', 'part.round'))

# find number of relationships reported in each age-age cat
tmp <- as.data.table(expand.grid(part.age = rage[, sort(unique(part.age))],
                                 cont.age = rage[, sort(unique(cont.age))],
                                 sex = rage[, sort(unique(sex))],
                                 part.comm = rage[, sort(unique(part.comm))],
                                 part.round = rage[, sort(unique(part.round))]))
ragea <- rage[, list(y = .N), by = c('part.age', 'cont.age', 'sex', 'part.comm', 'part.round')]

stopifnot(ragea[!is.na(cont.age), sum(y)] == nrow(rage[sexyear == 1]))
stopifnot(ragea[is.na(cont.age), sum(y)] == nrow(rage[sexyear == 2])) # they will be remove with the merge, ok

ragea <- merge(ragea, tmp, by = c('part.age', 'cont.age', 'sex', 'part.comm', 'part.round'), all.y = T)
ragea[is.na(y), y := 0]

# find number of participants in each age category 
tmp1 <- rage[, list(N = length(unique(pt_id))), by = c('part.age', 'sex', 'part.comm', 'part.round')]
stopifnot(tmp1[part.round == 'R015', sum(N)] == rage[part.round == 'R015', length(unique(pt_id))])

ragea <- merge(ragea, tmp1, by = c('part.age', 'sex', 'part.comm', 'part.round'), all.x = T)
ragea[is.na(N), N := 0]

# find census eligible count across round and communities by age and sex
tmp <- census_eligible_count[, list(T = sum(ELIGIBLE)), by = c('AGEYRS', 'SEX', 'ROUND', 'COMM')]
tmp[, ROUND := paste0('R0', ROUND)]
ragea <- merge(ragea, tmp, by.x = c('cont.age', 'sex', 'part.round', 'part.comm'), by.y = c('AGEYRS', 'SEX', 'ROUND', 'COMM'))

# find sex of participant and partner
ragea[, part.sex := ifelse(sex == 'M', 'M', 'F')]
ragea[, cont.sex := ifelse(part.sex == 'M', 'F', 'M')]
set(ragea, NULL, 'sex', NULL)

# input homosexual partnership with  0
tmp <- ragea[part.sex == 'F']
tmp[, `:=`(part.sex = 'M', y = 0, N = NULL)]
tmp <- merge(tmp, unique(ragea[part.sex == 'M', .(part.age, part.round, part.comm, N)]), by = c('part.age', 'part.round', 'part.comm'))

tmp1 <- ragea[part.sex == 'M']
tmp1[, `:=`(part.sex = 'F', y = 0, N = NULL)]
tmp1 <- merge(tmp1, unique(ragea[part.sex == 'F', .(part.age, part.round, part.comm, N)]), by = c('part.age', 'part.round', 'part.comm'))

ragea <- rbind(ragea, rbind(tmp, tmp1))

# find offset
ragea[, U := T * N]

# order
ragea <- ragea[order(part.round, part.comm, part.sex, cont.sex, part.age, cont.age)]

# check that the number of participants is correct
tmp <- unique(ragea[, .(part.comm, part.sex, part.round, part.age, N)])
tmp[, list(total_N = sum(N)), by = c('part.comm', 'part.sex', 'part.round')]
rage[part.comm == 'fishing' & part.round == 'R015' & sex == 'F', length(unique(pt_id))]

# plot
if(0){
  tmp <- ragea[part.round == 'R015']
  tmp[, part.sex := factor(part.sex, levels = c('M', 'F'))]
  tmp[, cont.sex := factor(cont.sex, levels = c('F', 'M'))]

  ggplot(tmp[y > 0 & part.comm == 'inland'], aes(y = cont.age, x = part.age)) +
    geom_point(aes(alpha = y)) +
    theme_bw() +
    labs(x = 'Participant age', y='Partner age') +
    facet_grid(cont.sex~part.sex, label = 'label_both')+
    ggtitle('Participant from inland communities')
  ggsave('~/Downloads/reported_relationships_inland.png', w = 6, h = 5.5)
  
  ggplot(tmp[y > 0 & part.comm == 'fishing'], aes(y = cont.age, x = part.age)) +
    geom_point(aes(alpha = y)) +
    theme_bw() +
    labs(x = 'Participant age', y='Partner age') +
    facet_grid(cont.sex~part.sex, label = 'label_both')+
    ggtitle('Participant from fishing communities')
  ggsave('~/Downloads/reported_relationships_fishing.png', w = 6, h = 5.5)

  tmp <- unique(ragea[, .(part.comm, part.sex, part.round, part.age, N)])
  tmp[, list(total_N = sum(N)), by = c('part.round', 'part.comm', 'part.sex')][order(part.round, part.comm, part.sex)]
  ggplot(tmp, aes(y = N, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm) +
    labs(x = 'Participant age', y = 'Number of participant')
  ggsave('~/Downloads/number_participant.png', w = 6, h = 5.5)
  
  tmp <- ragea[, list(total_y = sum(y)), by = c('part.comm', 'part.sex', 'part.round', 'part.age')]
  tmp[, list(total_y = sum(total_y)), by = c('part.round', 'part.comm', 'part.sex')]
  ggplot(tmp, aes(y = total_y, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm)  +
    labs(x = 'Participant age', y = 'Number of relationships reported')
  ggsave('~/Downloads/number_relationships.png', w = 6, h = 5.5)
  
  tmp <- ragea[, list(total_intensity = sum(y / N)), by = c('part.comm', 'part.sex', 'part.round', 'part.age')]
  tmp[, list(total_intensity = sum(total_intensity)), by = c('part.comm', 'part.sex', 'part.round')]
  ggplot(tmp, aes(y = total_intensity, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm)  +
    labs(x = 'Participant age', y = 'Number of relationships reported / Number of participants')
  ggsave('~/Downloads/number_relationships_by_N.png', w = 6, h = 5.5)
  
  tmp <- ragea[, list(total_intensity = sum(y / U)), by = c('part.comm', 'part.sex', 'part.round', 'part.age')]
  tmp[, list(total_intensity = sum(total_intensity)), by = c('part.comm', 'part.sex', 'part.round')]
  ggplot(tmp, aes(y = total_intensity, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm)  +
    labs(x = 'Participant age', y = 'Number of relationships reported / U')
  ggsave('~/Downloads/number_relationships_by_U.png', w = 6, h = 5.5)
  
  tmp <- census_eligible_count[, list(T = sum(ELIGIBLE)), by = c( 'SEX', 'ROUND', 'COMM')]
  tmp <- tmp[ROUND %in% c('15', '15S')]
  ggplot(tmp, aes(x = SEX)) + 
    geom_bar(aes(y = T), stat = 'identity') + 
    facet_grid(ROUND~COMM)  + 
    theme_bw() +
    labs(x = 'Sex', y = 'Census eligible count')
  ggsave('~/Downloads/census_eligible_count.png', w = 6, h = 5.5)
  
  tmp <- census_eligible_count[, list(T = sum(ELIGIBLE)), by = c( 'SEX', 'ROUND', 'COMM', 'AGEYRS')]
  tmp <- tmp[ROUND %in% c('15', '15S')]
  ggplot(tmp, aes(x = AGEYRS)) + 
    geom_line(aes(y = T, col = SEX), stat = 'identity') + 
    facet_grid(ROUND~COMM)  + 
    theme_bw() +
    labs(x = 'Age', y = 'Census eligible count')
  ggsave('~/Downloads/census_eligible_count_by_age.png', w = 6, h = 5.5)
}

# save
tmp <- ragea[part.round == 'R015']
write.csv(tmp, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220505.csv'), row.names = F)

tmp <- ra[round == 'R015']
write.csv(tmp, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_sexp1yr_220514.csv'), row.names = F)

