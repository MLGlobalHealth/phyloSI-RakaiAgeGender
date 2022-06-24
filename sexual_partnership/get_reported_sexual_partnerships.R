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
file.path.hiv <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'HIV_R15_R18_VOIs_220129.csv')

file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.census.eligible.individuals.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220620.csv')

community.keys <- as.data.table(read.csv(file.community.keys))
meta <- as.data.table(read.csv(file.meta.R1518))
census_eligible_count <- as.data.table(read.csv(file.census.eligible.individuals.count))
hiv <- as.data.table(read.csv(file.path.hiv))

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
                sexyear, sexp1yr, sexp1out,
                cuarvmed, 
                rldyslt1, rlwkslt1, rlmoslt1, rlyrslt1,
                rldyslt2, rlwkslt2, rlmoslt2, rlyrslt2,
                rldyslt3, rlwkslt3, rlmoslt3, rlyrslt3,
                rldyslt4, rlwkslt4, rlmoslt4, rlyrslt4,
                rltnage1, rltnage2, rltnage3, rltnage4,
                rltnyrs1, rltnyrs2, rltnyrs3, rltnyrs4, 
                rltnhh1, rltnhh2, rltnhh3, rltnhh4, 
                rltncm1, rltncm2, rltncm3, rltncm4)]
tmp <- unique(tmp)

tmp <- as.data.table( reshape2::melt(tmp, id.vars = c('pt_id', 'birth_date', 'visit_date', 'sex', 'comm_num', 'round', 'sexyear', 'sexp1yr', 'sexp1out', 'cuarvmed')) )
tmp[, relation := as.numeric(sub('.*(?=.$)', '', variable, perl=T))]
tmp[, variable := gsub('.{2}$', '', variable)]
age_preference <- as.data.table(reshape2::dcast(tmp, pt_id + birth_date + visit_date + sex + comm_num + round + sexyear + sexp1out + sexp1yr + cuarvmed + relation ~ variable, value.var = 'value'))

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
age_preference[, yrs_since_date_stop := days_since_stop / 365]

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


####################
# apply restrictions #
####################

# keep non-missing age index
rag <- age_preference[!is.na(age_index)]

# keep age index within census eligible age
rag <- rag[age_index >= 15 & age_index < 50]

# remove participant who did not say whether or not they had sexual intercourse in the last year
rag <- rag[sexyear %in% 1:2]
rag[, table(sexyear)]

# change coding of sexyear if not coherent with sexp1yr
rag[, table(sexp1yr)]
rag[, sexyear := ifelse(any(sexp1yr == 0), 2, sexyear), by = c('pt_id', 'round')]
rag[, sexp1yr := ifelse(all(sexyear == 2), 0, sexp1yr), by = c('pt_id', 'round')]
rag[, table(sexp1yr)]

# if no sex with anyone then no sex with anyone outside of the community 
rag[sexp1yr == 0, sexp1out := 0]
rag[sexyear == 1 & sexp1yr > 0  & round == 'R015' & sexp1out == 98, sexp1out := 0]

# keep partnership only within the last year
rag[sexyear == 2 & !is.na(age_partner), age_partner := NA]
rag[(sexyear == 1 & yrs_since_date_stop > 1) | (sexyear == 1 & is.na(yrs_since_date_stop)), age_partner := NA]

# keep age partner within census eligible age
rag[(sexyear == 1 & (age_partner < 15)), age_partner := NA]
rag[(sexyear == 1 & (age_partner > 69)), age_partner := NA]

# aggregate by 1 year age band
rag[, part.age := floor(age_index)]
rag[, cont.age := floor(age_partner)]


#####################################
# find meta patner-specific variable #
#####################################

# dcast cont.age
ragem <- dcast.data.table(rag, comm_num + round + sex + pt_id + sexp1yr + sexp1out + cuarvmed + part.age ~ relation, value.var = 'cont.age')
setnames(ragem, c('1', '2', '3', '4'), paste0('partner_age_', 1:4))

# dcast rltnh
rag[, cont.same.household := NA_character_]
rag[rltnh == 1, cont.same.household := 'YES']
rag[rltnh == 2, cont.same.household := 'NO']
tmp <- dcast.data.table(rag, comm_num + round + sex + pt_id + sexp1yr + sexp1out + cuarvmed + part.age ~ relation, value.var = 'cont.same.household')
setnames(tmp, c('1', '2', '3', '4'), paste0('partner_living_same_household_', 1:4))
ragem <- merge(ragem, tmp, by = c('comm_num', 'round', 'sex', 'pt_id', 'sexp1yr', 'sexp1out', 'cuarvmed', 'part.age'))

# dcast rltnc
rag[, cont.same.comm := NA_character_]
rag[rltnc == 1, cont.same.comm := 'YES']
rag[rltnc == 2, cont.same.comm := 'NO']
rag[rltnc == 7, cont.same.comm := 'DK']
tmp <- dcast.data.table(rag, comm_num + round + sex + pt_id + sexp1yr + sexp1out + cuarvmed + part.age ~ relation, value.var = 'cont.same.comm')
setnames(tmp, c('1', '2', '3', '4'), paste0('partner_living_same_community_', 1:4))
ragem <- merge(ragem, tmp, by = c('comm_num', 'round', 'sex', 'pt_id', 'sexp1yr', 'sexp1out', 'cuarvmed', 'part.age'))


#####################################
# find HIV patient-specific variable #
#####################################

# report ART
# ragem[, part.report.art := NA_character_ ]
# ragem[cuarvmed %in% c(0,2,8), part.report.art := 'NO']
# ragem[cuarvmed == 1, part.report.art := 'YES']
# set(ragem, NULL, 'cuarvmed', NULL)

# get hiv status
hiv[, pt_id := paste0('RK-', study_id)]
rhiv <- hiv[, .(pt_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
setnames(rhiv, 'hiv', 'part.hiv')
ragem <- merge(ragem, rhiv, by = c('pt_id', 'round'), all.x = T)
ragem[is.na(part.hiv) & round == 'R015',length(unique(pt_id))] #


#########################################
# Format meta patient-specific variable #
#########################################

# create community variable
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
rage <- merge(ragem, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')
set(rage, NULL, 'COMM_NUM_A', NULL)

# rename comm and round
setnames(rage, c('comm', 'comm_num', 'round', 'sex'), c('part.comm', 'part.comm_num', 'part.round', 'part.sex'))

# find census eligible count across round and communities by age and sex
tmp <- census_eligible_count[, list(part.T = sum(ELIGIBLE)), by = c('AGEYRS', 'SEX', 'ROUND', 'COMM_NUM')]
tmp[, ROUND := paste0('R0', ROUND)]
rage <- merge(rage, tmp, by.x = c('part.age', 'part.sex', 'part.round', 'part.comm_num'), by.y = c('AGEYRS', 'SEX', 'ROUND', 'COMM_NUM'))

# Z as a factor
rage[, Z := as.character(sexp1yr)]
rage[Z == '93', Z := '>3']
set(rage, NULL, 'sexp1yr', NULL)

# Z as a factor
rage[, sexp1out := as.character(sexp1out)]
rage[sexp1out == '93', sexp1out := '>3']

# save
tmp <- rage[part.round == 'R015']
write.csv(tmp, file = file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220620.csv'), row.names = F)


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

  ggplot(age_preference, aes(x = yrs_since_date_stop)) + 
    geom_histogram(bins = 50) + 
    facet_grid(round~., scale = 'free') + 
    scale_y_log10() + 
    geom_vline(xintercept = 1, linetype = 'dashed', col = 'darkred') + 
    theme_bw() + 
    labs(x = 'Time since end of the relationship (years)', y = 'relationships count') + 
    ggtitle('Among all participants')
  ggsave('~/Downloads/count_yrs_since_date_stop.png', w = 6, h = 5)
  
  #####
  
  tmp <- ragem[round == 'R015']
  tmp[, part.age.group := '35-49']
  tmp[part.age < 35, part.age.group := '25-34']
  tmp[part.age < 25, part.age.group := '15-24']
  
  tmp1 <- tmp[, list(count = length(unique(pt_id))), by = c('sexp1yr', 'part.age.group', 'sex')]
  tmp1[sexp1yr == 93, sum(count)]
  tmp1[sexp1yr == 0, sum(count)]
  tmp1 <- tmp1[sexp1yr != 93]
  ggplot(tmp1, aes(x = sexp1yr, y = count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid(sex~part.age.group, scale = 'free_y') + 
    labs(x='Number of different sexual partners in the last 12 months', y = 'count participants') 
  ggsave('~/Downloads/number_sexual_partners_12months_long.png', w = 9, h = 5)
  
  ggplot(tmp1, aes(x = sexp1yr, y = count))+ 
    geom_bar(stat = 'identity') + 
    facet_grid(sex~part.age.group, scale = 'free_y') + 
    labs(x='Number of different sexual partners in the last 12 months', y = 'count participants')+
    coord_cartesian(xlim = c(0, 15))
  ggsave('~/Downloads/number_sexual_partners_12months_short.png', w = 9, h = 5)
  
  
  ####
  tmp <- rage[part.round == 'R015' & Z != '>3']
  tmp[, y := as.numeric(!is.na(partner_age_1)) + as.numeric(!is.na(partner_age_2)) + as.numeric(!is.na(partner_age_3)) + as.numeric(!is.na(partner_age_4)) ]
  tmp <- tmp[, list(age_specific_relationship_reported = sum(y), relationship_reported = sum(as.numeric(Z))), by = c('part.comm', 'part.sex', 'part.age')]
  tmp[, list(age_specific_relationship_reported = sum(age_specific_relationship_reported)), by = c( 'part.comm', 'part.sex')]
  tmp <- melt.data.table(tmp, id.vars = c('part.comm', 'part.sex', 'part.age'))
  ggplot(tmp, aes(y = value, x = part.age, col = variable)) +
    geom_line() +
    theme_bw() +
    facet_grid(part.sex~part.comm)  +
    labs(x = 'Participant age', y = 'Number of relationships reported')
  ggsave('~/Downloads/number_relationships.png', w = 6, h = 5.5)
  
  ####
  
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

