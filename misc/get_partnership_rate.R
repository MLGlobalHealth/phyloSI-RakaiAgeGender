library(data.table)	
library(dplyr)	
library(ggplot2)
library(scales)
library(lubridate)
library(INLA)
library(reshape2)
library(RColorBrewer)

.year.diff <- function(x, y)
{
  if (!is.Date(x)) {x <- as.Date(x, format = '%Y-%m-%d')}
  if (!is.Date(y)) {y <- as.Date(y, format = '%Y-%m-%d')}
  lubridate::time_length(difftime(x, y),"years")
}

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
indir.repo = "~/git/phyloflows/"
outdir <- file.path(indir.deepsequence_analyses, 'preliminary', 'PartnerPreference')

file.meta.R1518 <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'quest_R15_R18_VoIs_220129.csv')
file.community.keys <- file.path(indir.deepsequence_analyses,'community_names.csv')
file.census.eligible.individuals.count <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_census_eligible_individuals_220411.csv')

community.keys <- as.data.table(read.csv(file.community.keys))
meta <- as.data.table(read.csv(file.meta.R1518))
census_eligible_count <- as.data.table(read.csv(file.census.eligible.individuals.count))

code_na = 97:99

source(file = file.path(indir.repo, 'misc', "functions", "get_partnership_rate-functions.R"))


#####################################
# find self-reported age preference #
#####################################

meta[, pt_id := paste0('RK-', study_id)]
meta[, visit_date := as.Date(intdate, '%d-%b-%y')]
meta[, birth_date := as.Date(paste0(gsub('.{2}$', '', birthdat), birthyr), '%d-%b-%Y')]
meta[is.na(birth_date) | nchar(birthyr) == 2 | birthyr < 1900, birth_date := visit_date - ageyrs*365]

# clean
tmp <- meta[, .(pt_id, birth_date, visit_date, sex, comm_num,
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

tmp <- as.data.table( reshape2::melt(tmp, id.vars = c('pt_id', 'birth_date', 'visit_date', 'sex', 'comm_num')) )
tmp[, relation := as.numeric(sub('.*(?=.$)', '', variable, perl=T))]
tmp[, variable := gsub('.{2}$', '', variable)]
age_preference <- as.data.table(reshape2::dcast(tmp, pt_id + birth_date + visit_date + sex + comm_num + relation ~ variable, value.var = 'value'))

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

# plot
if(0){
  rageac <- rage[, list(count = .N), by = c('age_index', 'sex', 'age_partner')]
  ggplot(rageac[sex == 'M'], aes(x = age_index, y = age_partner)) + 
    geom_point(aes(size = count)) + 
    theme_bw() + 
    labs(x = 'Participant s age (Male)', y='Partner s age (Female)') + 
    ggtitle('Partnership reported by males')
  ggsave(paste0(outdir, '-CountRelationships_Male_220412.png'), w = 7, h = 7)
  
  ggplot(rageac[sex == 'F'], aes(y = age_index, x = age_partner)) + 
    geom_point(aes(size = count)) + 
    theme_bw() + 
    labs(x = 'Partner s age (Male)', y='Participant s age (Female)') + 
    ggtitle('Partnership reported by females')
  ggsave(paste0(outdir, '-CountRelationships_Female_220412.png'), w = 7, h = 7)
}

# aggregate by 1 year age band
rage[, part.age := floor(age_index)]
rage[, cont.age := floor(age_partner)]

# find number of participants in each age category andnumber of relationships reported in each age-age cat
tmp1 <- rage[, list(N = length(unique(pt_id))), by = c('part.age', 'sex')]
ragea <- rage[, list(y = .N), by = c('part.age', 'cont.age', 'sex')]
tmp <- as.data.table(expand.grid(part.age = ragea[, sort(unique(part.age))], 
                                 cont.age = ragea[, sort(unique(cont.age))], 
                                 sex = ragea[, sort(unique(sex))]))
ragea <- merge(ragea, tmp, by = c('part.age', 'cont.age', 'sex'), all.y = T)
ragea <- merge(ragea, tmp1, by = c('part.age', 'sex'))
ragea[is.na(y), y := 0]

# find census eligible count across round and communities by age and sex
tmp <- census_eligible_count[, list(T = sum(ELIGIBLE)), by = c('AGEYRS', 'SEX')]
ragea <- merge(ragea, tmp, by.x = c('cont.age', 'sex'), by.y = c('AGEYRS', 'SEX'))
ragea[, U := T * N]

# find sex of participant and partner
ragea[, part.sex := ifelse(sex == 'M', 'M', 'F')]
ragea[, cont.sex := ifelse(part.sex == 'M', 'F', 'M')]
set(ragea, NULL, 'sex', NULL)

# input data without observation 
tmp <- ragea[part.sex == 'F']
tmp[, `:=`(part.sex = 'M', y = 0, N = NULL)]
tmp <- merge(tmp, unique(ragea[part.sex == 'M', .(part.age, N)]), by = 'part.age')

tmp1 <- ragea[part.sex == 'M']
tmp1[, `:=`(part.sex = 'F', y = 0, N = NULL)]
tmp1 <- merge(tmp1, unique(ragea[part.sex == 'F', .(part.age, N)]), by = 'part.age')

ragea <- rbind(ragea, rbind(tmp, tmp1))
ragea <- ragea[order(cont.sex, part.sex, cont.age, part.age)]

# get smooth estimate using both female and male reporting
AGEYRS <- 15:49
contact.matrix <- obtain.contact.matrix(ragea,  AGEYRS)
sragea <- contact.matrix[[2]]

# plot
sragea[, `Participant sex` := ifelse(part.sex == 'M', 'Male', 'Female')]
sragea[, `Partner sex` := ifelse(cont.sex == 'M', 'Male', 'Female')]
sragea[, `Participant sex` := factor(`Participant sex`, levels = c('Male', 'Female'))]
sragea[, `Partner sex` := factor(`Partner sex`, levels = c('Female', 'Male'))]

plot_crude_estimate(sragea)
ggsave(file = paste0(outdir, 'CrudeEstimates_220422.png'), w = 7, h = 6)
plot_smooth_estimate(sragea)
ggsave(file = paste0(outdir, 'SmoothEstimates_220422.png'), w = 7, h = 6)

# age profile of the partner
age_profile <- sragea[part.sex!=cont.sex]
age_profile[, total_c := sum(c), by = c('part.age', 'part.sex')]
age_profile[, prop := c / total_c]
age_profile <- age_profile[, list(M = as.numeric(matrixStats::weightedMedian(cont.age, prop)), 
                                  CL = as.numeric(modi::weighted.quantile(cont.age, prop, 0.1)), 
                                  CU = as.numeric(modi::weighted.quantile(cont.age, prop, 0.9))), by = c('part.age', 'part.sex')]
age_profile[, `Participant sex` := ifelse(part.sex == 'M', 'Male', 'Female')]
age_profile[, `Participant sex` := factor(`Participant sex`, levels = c('Male', 'Female'))]

plot_age_profile(age_profile)
ggsave(file = paste0(outdir, 'AgeProfile_220422.png'), w = 7, h = 4)

# save
write.csv(sragea, file.path(dirname(dirname(outdir)), 'RCCS_partnership_rate_220422.csv'), row.names = F)
# ragea <- as.data.table(read.csv(file.path(dirname(dirname(outdir)), 'RCCS_partnership_rate_220412.csv')))

# plot ratio
rageas <- copy(ragea)
rageas[, m_crude := y / N]
rageas[, female.age := ifelse(sex == 'F', part.age, cont.age)]
rageas[, male.age := ifelse(sex == 'M', part.age, cont.age)]
rageas <- dcast.data.table(rageas, female.age + male.age ~ sex, value.var = 'm_crude')
rageas[, ratio := ((F + 1) / (M + 1))]

ggplot(rageas, aes(x = male.age, y = female.age)) + 
  geom_raster(aes(fill = ratio)) + 
  labs(x = 'Male age', y = 'Female age', 
       fill = 'Females report more (>1) or less (<1) than males') + 
  scale_fill_gradient2(low = 'blue', mid = 'beige', high = 'red', midpoint = 1) +
  scale_x_continuous(expand = c(0,0))+
  scale_y_continuous(expand = c(0,0)) + 
  theme(legend.position = 'bottom')
ggsave(paste0(outdir, 'Differencemcrude_Female_220412.png'), w= 6, h = 7)
