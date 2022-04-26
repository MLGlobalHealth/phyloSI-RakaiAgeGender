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
                rldyslt1, rlwkslt1, rlmoslt1, rlyrslt1,
                rldyslt2, rlwkslt2, rlmoslt2, rlyrslt2,
                rldyslt3, rlwkslt3, rlmoslt3, rlyrslt3,
                rldyslt4, rlwkslt4, rlmoslt4, rlyrslt4,
                rltnage1, rltnage2, rltnage3, rltnage4,
                rltnyrs1, rltnyrs2, rltnyrs3, rltnyrs4)]
tmp <- unique(tmp)

tmp <- as.data.table( reshape2::melt(tmp, id.vars = c('pt_id', 'birth_date', 'visit_date', 'sex', 'comm_num', 'round')) )
tmp[, relation := as.numeric(sub('.*(?=.$)', '', variable, perl=T))]
tmp[, variable := gsub('.{2}$', '', variable)]
age_preference <- as.data.table(reshape2::dcast(tmp, pt_id + birth_date + visit_date + sex + comm_num + round + relation ~ variable, value.var = 'value'))

# find age index
age_preference[, days_since_start := 0]
age_preference[!day %in% code_na, days_since_start := days_since_start + day]
age_preference[!week %in% code_na, days_since_start := days_since_start + week*7]
age_preference[!month %in% code_na, days_since_start := days_since_start + month*31]
age_preference[!year %in% code_na, days_since_start := days_since_start + year*365]
age_preference[day %in% code_na & week %in% code_na & month %in% code_na & year %in% code_na, days_since_start := 10*365]
age_preference[, date_start := visit_date - days_since_start]
age_preference[, age_start := .year.diff(date_start, birth_date) ]

age_preference[, days_since_stop := 0]
age_preference[!rldysl %in% code_na, days_since_stop := days_since_stop + rldysl]
age_preference[!rlwksl %in% code_na, days_since_stop := days_since_stop + rlwksl*7]
age_preference[!rlmosl %in% code_na, days_since_stop := days_since_stop + rlmosl*31]
age_preference[!rlyrsl %in% code_na, days_since_stop := days_since_stop + rlyrsl*365]
age_preference[, date_stop := visit_date - days_since_stop]
age_preference[, age_stop := .year.diff(date_stop, birth_date) ]
age_preference <- age_preference[age_stop > age_start]

age_preference[, age_index := mean(c(age_start, age_stop)), by = c('pt_id', 'relation')]
age_preference[, relation_date := mean(c(date_start, date_stop)), by = c('pt_id', 'relation')]

# find age partner
age_preference[, age_partner_classification := NA_character_]
age_preference[rltnag == 1, age_partner_classification := 'older']
age_preference[rltnag == 2, age_partner_classification := 'younger']
age_preference[rltnag == 3, age_partner_classification := 'same']

age_preference[!rltnyr %in% code_na & age_partner_classification == 'younger', age_partner := age_index - rltnyr]
age_preference[!rltnyr %in% code_na & age_partner_classification == 'older', age_partner := age_index + rltnyr]
age_preference[age_partner_classification == 'same', age_partner := age_index ]

# clean
rage <- age_preference[!is.na(age_index) & !is.na(age_partner)]
rage <- rage[age_index >= 15 & age_index < 50]
rage <- rage[age_partner >= 15 & age_partner < 50]
rage <- rage[relation_date > as.Date('2010-01-01') & relation_date < as.Date('2020-01-01')]

# create variable
rage[, is_before_cutoff_date := relation_date < as.Date('2014-01-01')]
rage <- merge(rage, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
rage[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rage[, `Community index` := comm]

# aggregate by 1 year age band
rage[, part.age := floor(age_index)]
rage[, cont.age := floor(age_partner)]

# rename comm and round
setnames(rage, c('comm', 'round'), c('part.comm', 'part.round'))

# find number of participants in each age category and number of relationships reported in each age-age cat
tmp <- as.data.table(expand.grid(part.age = rage[, sort(unique(part.age))],
                                 cont.age = ragea[, sort(unique(cont.age))],
                                 sex = rage[, sort(unique(sex))],
                                 part.comm = rage[, sort(unique(part.comm))],
                                 part.round = rage[, sort(unique(part.round))]))
ragea <- rage[, list(y = .N), by = c('part.age', 'cont.age', 'sex', 'part.comm', 'part.round')]
ragea <- merge(ragea, tmp, by = c('part.age', 'cont.age', 'sex', 'part.comm', 'part.round'), all.y = T)
ragea[is.na(y), y := 0]
tmp1 <- rage[, list(N = length(unique(pt_id))), by = c('part.age', 'sex', 'part.comm', 'part.round')]
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

# find offset
ragea[, U := T * N]

ragea <- rbind(ragea, rbind(tmp, tmp1))
ragea <- ragea[order(part.round, part.comm, part.sex, cont.sex, part.age, cont.age)]

# plot
if(0){
  tmp <- ragea[part.round == 'R015']
  tmp[, part.sex := factor(part.sex, levels = c('M', 'F'))]
  tmp[, cont.sex := factor(cont.sex, levels = c('F', 'M'))]

  ggplot(tmp[y > 0 & part.comm == 'inland'], aes(y = cont.age, x = part.age)) +
    geom_point(aes(alpha = y)) +
    theme_bw() +
    labs(x = 'Participant s age', y='Partner s age') +
    facet_grid(cont.sex~part.sex, label = 'label_both')

  ggplot(tmp[y > 0 & part.comm == 'fishing'], aes(y = cont.age, x = part.age)) +
    geom_point(aes(alpha = y)) +
    theme_bw() +
    labs(x = 'Participant s age', y='Partner s age') +
    facet_grid(cont.sex~part.sex, label = 'label_both')

  tmp <- unique(ragea[, .(part.comm, part.sex, part.round, part.age, N, U)])
  ggplot(tmp, aes(y = N, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm)
  ggplot(tmp, aes(y = U, x = part.age, col = part.sex)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.round~part.comm)
}

# save
tmp <- ragea[part.round == 'R015']
write.csv(tmp, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220426.csv'), row.names = F)
