library(data.table)
library(ggplot2)

source('~/git/phyloflows/misc/functions/get_census_eligible_count-functions.R')

outdir <- "~/Box\ Sync/2021/phyloflows/incidence"

infile					<- "~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/rakai_R1516_elibility.rda"
eligible <- RakaiFull.phylogeography.220310.data.eligibility.participation.sequenced(infile)

tmp		<- data.table(	COMM_NUM=	c("1","2","3","4","5","6","7","8","9","14","15","16","18","19","22","23","24","25","29","32","33","34","35","36","38","40","44","45","46","51","52","53","54","55", "56","57","58","59","60","61","62","65","67","74","77","81","84","89","94","95","103","106","107","108","109","120","177", "183", "256", "370","391","401","451", "468","602", "754", "755", "760", "770","771","772","773","774","776"),
                    COMM_TYPE=	c("T","A","A","T","A","A","A","A","A", "A", "A", "T", "A", "A", "T", "A", "T", "A", "A", "A", "T", "A", "A", "A", "F", "A", "A", "A", "A", "T", "A", "A", "A", "A",  "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A",  "A","A",   "A",  "T",  "A",  "A",  "A",  "A",   "A",   "A",   "A",  "A", "A",  "A",    "A",  "A",  "A",    "A",   "A", "F",  "F",  "A",  "A",  "F",   "T"))
set(tmp, NULL, 'COMM_TYPE', tmp[, as.character(factor(COMM_TYPE, levels=c('A','T','F'), labels=c('agrarian','trading','fisherfolk')))])
# set(tmp, NULL, 'COMM_NUM', tmp[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM_RAW))))])
eligible		<- merge(eligible, tmp, by='COMM_NUM')

file.community.keys <- file.path('~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/','community_names.csv')
community.keys <- as.data.table(read.csv(file.community.keys))
community.keys[, COMM := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']
eligible <- merge(eligible, community.keys, by = 'COMM_NUM_RAW')

eligible <- eligible[AGEYRS >= 15 & AGEYRS <= 49]
setnames(eligible, 'VISIT', 'ROUND')
eligible <- eligible[, list(ELIGIBLE = .N, PARTICIPANTS = sum(PARTICIPATED)), by = c('ROUND', 'COMM', 'AGEYRS', 'SEX')]
write.csv(eligible, "~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/RCCS_census_eligible_220311.csv", row.names = F)

# plot
if(0){
  tmp <- melt.data.table(eligible, id.vars = c('ROUND', 'COMM', 'AGEYRS', 'SEX'))
  ggplot(tmp[ROUND == 16], aes(x = AGEYRS, y = value)) + 
    geom_bar(aes(fill = variable), stat = 'identity', position = "dodge") + 
    labs(y = 'Count', x = 'Age (year)') + 
    facet_grid(SEX~COMM, label = 'label_both') +
    theme_bw() + 
    theme(legend.position = 'bottom') 
  ggsave(file.path(outdir, 'census_eligible_count_round16.png'), w = 9, h = 6)
  
  tmp <- eligible[, list(ELIGIBLE = sum(ELIGIBLE), PARTICIPANTS = sum(PARTICIPANTS)), by = c('ROUND', 'COMM')]
  tmp <- melt.data.table(tmp, id.vars = c('ROUND', 'COMM'))
  ggplot(tmp, aes(x = as.factor(ROUND), y = value)) + 
    geom_bar(aes(fill = variable), stat = 'identity', position = "dodge") + 
    labs(y = 'Count', x = 'ROUND') + 
    facet_grid(.~COMM, label = 'label_both') + 
    theme_bw()
  ggsave(file.path(outdir, 'census_eligible_count.png'), w = 7, h = 5)
  
}

# load Adam's incidence estimates
infile					<- "~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/Rakai_incpredictions_220310.csv"

incidence <- as.data.table(read.csv(infile))
colnames(incidence) <- toupper(colnames(incidence))
setnames(incidence, 'AGE', 'AGEYRS')
incidence[, COMM := 'inland']
incidence[, SEX := substring(SEX, 1, 1)]

setnames(eligible, 'ROUND', 'ROUND_CENSUS')

di <- merge(incidence, eligible[ROUND_CENSUS == 16], by = c('COMM', 'AGEYRS', 'SEX'))
di[, INFECTIONS := INCIDENCE * ELIGIBLE]
di[, INFECTIONS_LB := LB * ELIGIBLE]
di[, INFECTIONS_UB := UB * ELIGIBLE]
write.csv(di, "~/Box Sync/2019/ratmann_pangea_deepsequencedata/live//RCCS_R15_R18/RCCS_incident_cases_220311.csv", row.names = F)


ggplot(di, aes(x = AGEYRS, y = INCIDENCE)) + 
  geom_bar(stat = 'identity', position = "dodge") + 
  labs(y = 'Incidence rate per 1 PY', x = 'Age (year)') + 
  facet_grid(SEX~., label = 'label_both') +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  ggtitle("Fishing community in round 16")
ggsave(file.path(outdir, 'incidence_rate.png'), w = 7, h = 5)


ggplot(di, aes(x = AGEYRS, y = INFECTIONS)) + 
  geom_bar(aes(fill = MODEL), stat = 'identity', position = "dodge") + 
  labs(y = 'Expected number of infection', x = 'Age (year)') + 
  facet_grid(SEX~., label = 'label_both') +
  theme_bw() + 
  theme(legend.position = 'bottom') +
  ggtitle("Fishing community in round 16")
ggsave(file.path(outdir, 'infections.png'), w = 7, h = 5)


