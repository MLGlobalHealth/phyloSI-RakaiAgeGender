library(data.table)
library(ggplot2)

# laptop
if(dir.exists('~/Box\ Sync/2021/ratmann_deepseq_analyses/'))
{
  indir <- '~/git/phyloflows'
  indir.deepsequence_analyses <- '~/Box\ Sync/2021/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI/'
  indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'
  outdir <- '~/Box\ Sync/2021/phyloflows/'

  jobname <- 'test_new'
  stan_model <- 'gp_220905a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/home/andrea'))
{
  indir <-'~/git/phyloflows'
  indir.deepsequence_analyses   <- '~/Documents/Box/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '~/Documents/Box/ratmann_pangea_deepsequencedata/'
  outdir <- '~/Documents/Box/2021/phyloflows'

  jobname <- 'test'
  stan_model <- 'gp_220901a'
  outdir <- file.path(outdir, paste0(stan_model, '-', jobname))
  dir.create(outdir)
}

if(dir.exists('/rds/general'))
{
  indir.deepsequence_analyses   <- '/rds/general/project/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- '/rds/general/project/ratmann_pangea_deepsequencedata/live/'
}

if(dir.exists('/Users/dalma/'))
{
  indir <-'C:/Users/dalma/Desktop/IMPERIAL PROJECT/Code/IMP'
  indir.deepsequence_analyses   <-  'R:/projects/ratmann_deepseq_analyses/live/PANGEA2_RCCS1519_UVRI'
  indir.deepsequencedata <- 'R:/projects/ratmann_pangea_deepsequencedata/live'
  outdir <- 'C:/Users/dalma/Desktop/IMPERIAL PROJECT/Results'
  
 
  outdir <- file.path(outdir, 'Exploratory_data_analysis')
  dir.create(outdir)
}

args_line <-  as.list(commandArgs(trailingOnly=TRUE))
print(args_line)
if(length(args_line) > 0)
{
  stopifnot(args_line[[1]]=='-indir')
  stopifnot(args_line[[3]]=='-outdir')
  stopifnot(args_line[[5]]=='-stan_model')
  stopifnot(args_line[[7]]=='-jobname')
  indir <- args_line[[2]]
  outdir <- args_line[[4]]
  stan_model <- args_line[[6]]
  jobname <- args_line[[8]]
}

# TODO: set an output directory you want
# outfile <- file.path(outdir, paste0(stan_model,'-', jobname))
# outfile.figures <- file.path(outdir, 'figures', paste0(stan_model,'-', jobname))
# outdir.table <- file.path(outdir, 'tables', paste0(stan_model,'-', jobname))
# if(!dir.exists(dirname(outfile.figures))) dir.create(dirname(outfile.figures))
# if(!dir.exists(dirname(outdir.table))) dir.create(dirname(outdir.table))

# indicators
include.only.heterosexual.pairs <- T
threshold.likely.connected.pairs <- 0.5
use.tsi.estimates <- F
remove.inconsistent.infection.dates <- F
remove.young.individuals <- T
remove.missing.community.recipient <- T

remove.neuro.individuals <- T
only.transmission.after.start.observational.period <- T
only.transmission.before.stop.observational.period <- T
only.transmission.same.community <- F

# file paths
file.path.chains.data <- file.path(indir.deepsequence_analyses,'211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd/Rakai_phscnetworks.rda')


# from misc/
file.path.meta <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'Rakai_Pangea2_RCCS_Metadata_household_20220913.RData')

# load functions
source(file.path(indir, 'functions', 'utils.R'))
source(file.path(indir, 'functions', 'summary_functions.R'))
source(file.path(indir, 'misc', 'functions', 'household_helpers.R'))

# load anonymous aid TODO: add file.anonymisation.keys path
#aik <- fread(file.anonymisation.keys)
#aik$X <- NULL TODO: do I need it? I think no!



# load chains
load(file.path.chains.data)
dchain <- as.data.table(dchain)
dchain

load(file.path.meta)
meta_data

dchain <- keep.likely.transmission.pairs(dchain,threshold.likely.connected.pairs)
#NB:1020 rows now
###### EXPLORATORY ANALYSIS ############
dchain[, range(SCORE_LINKED)]
hist(dchain$SCORE_LINKED)
plot(dchain$SCORE_DIR_12,dchain$SCORE_DIR_21) #NB: some of them don't sum to one! Actually, they do after the new threshold =0.5
nrow(dchain[SCORE_DIR_12+ SCORE_DIR_21 >1.01 | SCORE_DIR_12+ SCORE_DIR_21<0.99]) #96 cases with problems, with old threshold =0.4
#pairs.all <- pairs.get.meta.data(dchain, meta_data, aik)

#Extract only useful info
tmp<-dchain[,c('CLU_SIZE','SOURCE','RECIPIENT')]
unique(tmp) #check is ok
tmpsource<-merge(tmp,meta_data,by.x='SOURCE',by.y='aid')
meta_data[,uniqueN(comm),by='study_id'][,all(V1==1)]
index<-meta_data[,uniqueN(comm),by='study_id'][(V1!=1),study_id]
meta_data[study_id %in% index,]
#some individuals change community, it is not constant over time
meta_data[,uniqueN(comm),by='hh_num'][,mean(V1==1)]
#now, this is shit. Does the numeration of households changes over different rounds?
#Probably, some people have two households. They have family in inland and they work at fishing. Is there
#data about the location of the household? solution we found:
tmp<-meta_data[,list(nf=sum(comm=='fishing'),ni=sum(comm=='inland')),by='hh_num']
#tmp$hh_comm<-fifelse(tmp$ni>tmp$nf,'inland','fishing')
tmp[,hh_comm:=fifelse(ni>nf,'inland','fishing')]
meta_data<-merge(meta_data,tmp[,c('hh_num','hh_comm')],by='hh_num')
tmp<-meta_data[,uniqueN(hh_num),by=c('study_id','sex')]
tmp[,mean(V1==1)]
tmp[,table(sex)] #There are more females than males, why? Selection bias?
tmp[V1>1,table(sex)] #still more women moving hhh rather than men. They change it after getting married?
#93% of people stay in the same household
View(tmp)
meta_data[,unique(sex),by='study_id'][,table(V1)]

dchain[, mean(unique(c(SOURCE, RECIPIENT)) %in% meta_data$aid) ] #NB: not all the sources and recipients
#in dchain have a corresponding person in meta_data

tmp<-dchain[,c('CLU_SIZE','SOURCE','RECIPIENT')]
tmp2<-meta_data[,c('study_id','aid','date_birth','date_first_positive','date_last_negative','sex')]
tmp2<-unique(tmp2)
tmpsource<-merge(tmp,tmp2,by.x='SOURCE',by.y='aid')
tmprecipient<-merge(tmp,tmp2,by.x='RECIPIENT',by.y='aid')
cols <- setdiff(names(tmpsource), c('RECIPIENT', 'CLU_SIZE', 'SOURCE'))
newcols <- paste0(cols, '_source')
setnames(tmpsource,cols,newcols)
newcols <- paste0(cols, '_recipient')
setnames(tmprecipient,cols,newcols)

#now, we define a dataset of transmissions with only the features which remain constant over time
d_trans_constant<-merge(tmpsource,tmprecipient,by=c('SOURCE','RECIPIENT','CLU_SIZE'))
#NB: we lose 283 rows when merging the transmission dataset with the meta_data, because of missing AIDs
#in either the source or the recipient
d_trans_constant<-d_trans_constant[sex_recipient!=sex_source,]
#NB: only 505, is that correct?
#function which returns whether source and recipient have lived in the same hh(even in different rounds)
d_trans_constant<-d_trans_constant[date_last_negative_source<date_first_positive_recipient| is.na(date_last_negative_source)==TRUE | is.na(date_first_positive_recipient)==TRUE,]
#we have just excluded inconsistent infection dates
count_same_hh<-function(dchain,meta_data){
  source<-dchain$SOURCE
  recipient<-dchain$RECIPIENT
  meta_useful<-meta_data[,c('hh_num','aid')]
  same_hh<-c()
  same_hh_comm<-rep(NA,length(source))
  for (i in 1:length(source)){
    hh_source<-meta_useful[aid==source[i],hh_num]
    hh_recipient<-meta_useful[aid==recipient[i],hh_num]
    hh_int <- intersect(hh_source,hh_recipient)
    same_hh[i]<-length(hh_int)
    if(same_hh[i]>0)
      same_hh_comm[i]<-unique(meta_data[hh_num==hh_int[1],hh_comm])
  }
  return (list(same_hh=same_hh,same_hh_comm=same_hh_comm))
}

d_trans_constant$same_hh<-count_same_hh(d_trans_constant,meta_data)$same_hh
source1<-d_trans_constant[same_hh==2,SOURCE][1]
recipient1<-d_trans_constant[same_hh==2,RECIPIENT][1]
meta_data[aid==source1|aid==recipient1,] #they are both houses in a fishing community, nice, they probably moved together
source2<-d_trans_constant[same_hh==2,SOURCE][2]
recipient2<-d_trans_constant[same_hh==2,RECIPIENT][2]
meta_data[aid==source2|aid==recipient2,] #they are both houses in a inland community, nice, they probably moved together
#NB: I have already modified the function above to account for it, and to add same_hh_comm!

d_trans_constant$same_hh_comm<-count_same_hh(d_trans_constant,meta_data)$same_hh_comm
table(d_trans_constant$same_hh)
#nice!
table(d_trans_constant$same_hh_comm) #TODO: compare with how many people are in inland and fishing community
plot(as.factor(d_trans_constant[same_hh>0,]$same_hh_comm))
population_census<-fread("R:/projects/ratmann_pangea_deepsequencedata/live/RCCS_R15_R18/RCCS_census_eligible_individuals_220830.csv")
#data about all the population, according to this census, divided by age and round, it indicates the 
#number of eligible people. Probably pretty useless indeed.
#Now I find the meta_data only for those people represented in pairs. Then, for every person I find their
#community using the mode and assess whether they are moving fishing<->inland or no
vectoraid<-unique(c(d_trans_constant$SOURCE,d_trans_constant$RECIPIENT))
tmp<-meta_data[aid %in% vectoraid,]
find_comm_and_moving<-function(meta_data){
  aid<-unique(meta_data$aid)
  comm<-c()
  moving_fi<-rep(FALSE,length(aid))
  for (i in c(1:length(aid))){
    comm[i]<-mode(meta_data[aid==aid[i],comm])
    comm[i]<-names(sort(summary(as.factor(meta_data[aid==aid[i],]$comm)),decreasing = TRUE))[1]
    if(length(names(sort(summary(as.factor(meta_data[aid==aid[i],]$comm)),decreasing = TRUE)))>1)
      moving_fi[i]<-TRUE
  }
  return(list(aid=aid,comm=comm,moving_fi=moving_fi))
}
res<-find_comm_and_moving(tmp)
#looks like it works :)
tmp2<-data.frame(aid=res$aid,comm=res$comm,moving_fi=res$moving_fi)
tmpsource<-merge(d_trans_constant,tmp2,by.x='SOURCE',by.y='aid')
tmprecipient<-merge(d_trans_constant,tmp2,by.x='RECIPIENT',by.y='aid')
cols <- c('comm','moving_fi')
newcols <- paste0(cols, '_source')
setnames(tmpsource,cols,newcols)
newcols <- paste0(cols, '_recipient')
setnames(tmprecipient,cols,newcols)
tmp<-merge(tmpsource,tmprecipient,by=colnames(d_trans_constant))
d_trans_constant<-copy(tmp)
#looks like there are many inconsistencies between hh_comm and the comm of different people

#From now on, we imagine our estimates of same_hh are correct, and build our analysis on that
#In case, we will just change them later and use the following code

############Let's start with the QUESTIONS################
#1) How many m->f in HH? and out of HH?
num_mf_hh<-nrow(d_trans_constant[same_hh>=1 & sex_source=='M' & sex_recipient=='F',])
num_fm_hh<-nrow(d_trans_constant[same_hh>=1 & sex_source=='F' & sex_recipient=='M',])
num_mf_outofhh<-nrow(d_trans_constant[same_hh==0 & sex_source=='M' & sex_recipient=='F',])
num_fm_outofhh<-nrow(d_trans_constant[same_hh==0 & sex_source=='F' & sex_recipient=='M',])
num_hh<-num_mf_hh+num_fm_hh
num_outofhh<-num_mf_outofhh+num_fm_outofhh
num_mf_hh/num_hh #59% of infections in HH are M->F
num_mf_outofhh/num_outofhh #55.8% of infections out of HH are M->F
diff_age<-(d_trans_constant$date_birth_recipient-d_trans_constant$date_birth_source)/365.25 #pretty naive approach,
#do we prefer to find the difference in years with lubridate and the remainders of days and weeks?
hist(as.numeric(diff_age))
diff_age_hh<-(d_trans_constant[same_hh>=1,]$date_birth_recipient-d_trans_constant[same_hh>=1,]$date_birth_source)/365.25
hist(as.numeric(diff_age_hh))
diff_age_outofhh<-(d_trans_constant[same_hh==0,]$date_birth_recipient-d_trans_constant[same_hh==0,]$date_birth_source)/365.25
hist(as.numeric(diff_age_outofhh)) 


###########Here, we should start with the PLOTS ################
#I think it is better to modify the variable same_hh: instead of 0,1,2, only 0,1 (better FALSE,TRUE) (we do not really care
#if the common HHs are 2, there is not much difference)
d_trans_constant[same_hh==2,]$same_hh<-rep(1,length(d_trans_constant[same_hh==2,]$same_hh))

plot(as.factor(d_trans_constant$same_hh):as.factor(d_trans_constant$comm_source))
plot(as.factor(d_trans_constant$same_hh):as.factor(d_trans_constant$sex_source))
ggplot(data = d_trans_constant, aes(x = as.logical(same_hh), fill = sex_source)) +
  geom_bar()+
  labs(title='Number of infections in and outside of the household',
       subtitle='Grouped by sex',
       x='Same household',
       y='Number of infections',
       fill='Sex of the source')

ggplot(data = d_trans_constant, aes(x = as.logical(same_hh), fill = comm_source)) +
  geom_bar()+
  labs(title='Number of infections in and outside of the household',
       subtitle = 'Grouped by community of the source',
       x='Same household',
       y='Number of infections',
       fill='Community of the source')

ggplot(data = d_trans_constant, aes(x = as.logical(same_hh), y = (date_birth_recipient-date_birth_source)/365.25)) +
  geom_boxplot()+
  labs(title='Age difference between source and recipient in and outside of the household',
       x='Same household',
       y='Age difference between source and recipient')
#Does it work? add better colors!
ggplot(data = d_trans_constant, aes(x = comm_source, y = (date_birth_recipient-date_birth_source)/365.25)) +
  geom_boxplot()+
  labs(title='Age difference between source and recipient',
       subtitle='By community of the source',
       x='Community of the source',
       y='Age difference between source and recipient')
ggplot(data = d_trans_constant, aes(x = sex_source, y = (date_birth_recipient-date_birth_source)/365.25)) +
  geom_boxplot()+
  labs(title='Age difference between source and recipient',
       subtitle='By sex of the source',
       x='Sex of the source',
       y='Age difference between source and recipient')
#as expected, older men tend to infect younger women

if(0)
{
        # NEED TO GET PAIRS pairs.all

        if(only.transmission.after.start.observational.period){
          cat('\nFor inland excluding recipients infected before ', as.character(start_observational_period_inland), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_inland & COMM.RECIPIENT == 'inland')]
          
          cat('\nFor fishing excluding recipients infected before ', as.character(start_observational_period_fishing), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT < start_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
          
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
        if(only.transmission.before.stop.observational.period){
          cat('\nFor inland excluding recipients infected after ', as.character(stop_observational_period_inland), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_inland & COMM.RECIPIENT == 'inland')]
          
          cat('\nFor fishing excluding recipients infected after ', as.character(stop_observational_period_fishing), '\n')
          cat('Removing ', nrow(pairs.all[DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing']), ' pairs\n')
          pairs.all <- pairs.all[!(DATE_INFECTION.RECIPIENT > stop_observational_period_fishing & COMM.RECIPIENT == 'fishing')]
          
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
        if(remove.missing.community.recipient){
          cat('\nExcluding recipients without community \n')
          cat('Removing ', nrow(pairs.all[is.na(COMM.RECIPIENT)]), ' pairs\n')
          pairs.all <- pairs.all[!is.na(COMM.RECIPIENT)]
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        } 
        if(only.transmission.same.community ){
          cat('\nExcluding transmission events between communities (I->F or F->I) \n')
          cat('Removing ', nrow(pairs.all[COMM.SOURCE != COMM.RECIPIENT]), ' pairs\n')
          pairs.all <- pairs.all[COMM.SOURCE == COMM.RECIPIENT]
          cat('resulting in a total of ', nrow(pairs.all),' pairs\n\n')
        }
}


