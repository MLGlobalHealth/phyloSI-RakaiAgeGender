library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(rstan)
library(haven)

# directory of the repository
gitdir <- here::here()

# load file paths
source(file.path(gitdir, 'config.R'))

# outdir for figures
if(usr == 'alexb' | usr == 'melodiemonod' | usr == 'andrea') {
  outdir <- file.path(indir.deepsequence_analyses, 'PANGEA2_RCCS', 'participants_count_by_gender_loc_age')
}else{
  outdir #<- #TODO
}
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

file.exists(c(
  file.community.keys ,
  file.path.hiv,
  file.path.quest,
  file.characteristics_sequenced_ind_R14_18,
  file.path.metadata))  |> all() |> stopifnot()

# load functions 
source( file.path(gitdir.R,'naturemed_reqs.R') )
source( file.path(gitdir.R,'utils.R') )
naturemed_reqs()

# load files
community.keys <- fread(file.community.keys)

# rounds of interest
df_round <- rbind(data.table(COMM = 'inland', ROUND = paste0('R0', 10:18)),
                  data.table(COMM = 'fishing', ROUND = paste0('R0', c(10:15,'15S', 16:18))))

# community
community.keys[, comm := ifelse(strsplit(as.character(COMM_NUM_A), '')[[1]][1] == 'f', 'fishing', 'inland'), by = 'COMM_NUM_A']

###############################################

# GET META DATA

###############################################

#
# Meta data
#

meta_data <- fread(file.path.metadata) #additional meta_data

# find age
meta_data[, date_birth := as.Date(paste0(birthyr, '-', birthmo, '-', '01'), format = '%Y-%m-%d')]
meta_data[, AGEYRS := round(lubridate::time_length(difftime(sample_date, date_birth),"years"))]
meta_data[is.na(AGEYRS), AGEYRS := round(lubridate::time_length(difftime(firstposvd, date_birth),"years"))]
meta_data[is.na(AGEYRS)]

# restrict age
meta_data <- meta_data[AGEYRS > 14 & AGEYRS < 50]

# find community
meta_data[, COMM := 'inland']
meta_data[LakeVictoria_FishingCommunity == 'yes', COMM := 'fishing']

# find hiv status
meta_data[, HIV := ifelse(is.na(firstposvd), 'N', 'P')]

# subset round 14 to 30 continuously surveyed communities
idr14not30comm <- meta_data[(round=='14' & !COMM %in% c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754)), unique(study_id)]
meta_data <- meta_data[!(round=='14' & !COMM %in% c(1, 2, 4, 5, 6, 7, 8, 16, 19, 22, 24, 29, 33, 34, 40, 56, 57, 58, 62, 74, 77, 89, 94, 106, 107, 108, 120, 391, 602, 754))]

# find art use
meta_data[, ART := artslfuse == 'yes']

# keep variable of interest
meta_data[, round := paste0('R0', round)]
colnames(meta_data) <- toupper(colnames(meta_data))
meta_data <- meta_data[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, HIV, ART)]

# set 15.1 to be 15S
meta_data[ROUND == 'R015.1', ROUND := 'R015S']


#
# Quest
#

rin <- fread(file.path.quest, select=c('ageyrs', 'round', 'study_id', 'sex', 'comm_num', 'arvmed', 'cuarvmed'))

# find  community
rinc <- merge(rin, community.keys, by.x = 'comm_num', by.y = 'COMM_NUM_RAW')

# to upper
colnames(rinc) <- toupper(colnames(rinc))

# restrict age
rinc <- rinc[AGEYRS > 14 & AGEYRS < 50]

# get ART status
rinc[!ROUND %in% c('R016', 'R017', 'R018'), ART := ARVMED ==1]
rinc[!ROUND %in% c('R016', 'R017', 'R018') & is.na(ARVMED), ART := F]
rinc[ROUND == 'R016', ART := ARVMED ==1 | CUARVMED ==1]
rinc[ROUND == 'R016' & (is.na(ARVMED) | is.na(CUARVMED)), ART := F]
rinc[ROUND %in% c('R017', 'R018'), ART := CUARVMED ==1]
rinc[ROUND %in% c('R017', 'R018') & is.na(CUARVMED), ART := F]

# art was not reported in round 10
rinc[ROUND == 'R010', ART := NA]

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], rinc[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
rinc <- rbind(rinc[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART)])

# SET ROUND 15S IN INLAND AS 15
rincp <- copy(rinc)
rincp[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
rincp[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
rincp <- rincp[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

#
# HIV
#

hiv <- fread(file.path.hiv)

hiv[, round := gsub(' ', '', round)] # remove space in string

# get hiv status
rhiv <- hiv[, .(study_id, round, hiv)]
rhiv[, round := gsub(" ", '', round, fixed = T)]
colnames(rhiv) <- toupper(colnames(rhiv))

# add meta data from Joseph
hivs <- merge(rhiv, rinc, by = c('STUDY_ID', 'ROUND'))

# add meta data from Kate
tmp <- anti_join(meta_data[, .(STUDY_ID, ROUND)], hivs[, .(STUDY_ID, ROUND)], by = c('STUDY_ID', 'ROUND'))
tmp <- merge(tmp, meta_data, by = c('STUDY_ID', 'ROUND'))
hivs <- rbind(hivs[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)], tmp[, .(STUDY_ID, SEX, ROUND, COMM, AGEYRS, ART, HIV)])

# SET ROUND 15S IN INLAND AS 15
hivs[, PARTICIPATED_TO_ROUND_RO15 := any(ROUND == 'R015'), by= 'STUDY_ID']
hivs[ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == F, ROUND := 'R015']
hivs <- hivs[!(ROUND == 'R015S' & COMM == 'inland' & PARTICIPATED_TO_ROUND_RO15 == T)]

# restric age
hivs <- hivs[AGEYRS > 14 & AGEYRS < 50]


#################################

# GET PARTICIPANTS #

#################################


# find participant
part <- rincp[, list(PARTICIPANT = .N,  TYPE = 'Total'), by = c('COMM', 'ROUND')]
part <- rbind(part, rincp[SEX == 'F', list(PARTICIPANT = .N,  TYPE = 'Female'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M', list(PARTICIPANT = .N,  TYPE = 'Male'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Female, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Female, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'F' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Female, 35-49'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS < 25, list(PARTICIPANT = .N,  TYPE = 'Male, 15-24'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS > 24 & AGEYRS < 35, list(PARTICIPANT = .N,  TYPE = 'Male, 25-34'), by = c('COMM', 'ROUND')])
part <- rbind(part, rincp[SEX == 'M' & AGEYRS > 34, list(PARTICIPANT = .N,  TYPE = 'Male, 35-49'), by = c('COMM', 'ROUND')])

# keep round of interest
part <- merge(part, df_round, by = c('COMM', 'ROUND'))


########################################

# GET HIV-POSITIVE AMONG PARTICIPANTS #

########################################

# keep rounds of interest and only positives
hivsp <- merge(hivs, df_round, by = c('COMM', 'ROUND'))
hivsp <- subset(hivsp,HIV=='P')

# group into analysis age groups
hivsp[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]

# save individual-level data
hivi <- copy(hivsp)

# keep first round only so age is unique per participant
hivsp[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
setkey(hivsp,STUDY_ID,r)
hivsp <- hivsp[,.SD[1],by = STUDY_ID]

# count unique participants with positive test during rounds
hivp <- hivsp[,list(HIV = length(unique(STUDY_ID))), by = c('COMM','SEX','AGEGP')]

tot1 <- hivsp[,list(SEX = 'Total', AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM')]
tot2 <- hivsp[,list(AGEGP = 'Total', HIV = length(unique(STUDY_ID))), by = c('COMM','SEX')]
hivp <- rbind(hivp,tot1,tot2)

hivp[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
hivp[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
hivp[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]
hivp <- hivp[order(COMM,SEX,AGEGP),]


######################################################

# GET HIV-POSITIVE AND USING ART AMONG PARTICIPANTS #

######################################################

# keep HIV positive and find art status
sart <- hivs[HIV == 'P']

# age groups
sart[, AGEGP:= cut(AGEYRS,breaks=c(15,25,35,50),include.lowest=T,right=F,
                   labels=c('15-24','25-34','35-49'))]

# keep round of interest
sart <- merge(sart, df_round, by = c('COMM', 'ROUND'))

# save individual-level record
arti <- copy(sart)

# keep first visit only
sart[, r:= as.numeric(gsub('R','',gsub('S','.1',ROUND)))]
setkey(sart,STUDY_ID,r)
sart <- sart[,.SD[1],by = STUDY_ID]

art <- sart[,list(NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM','SEX','AGEGP')]

tot1 <- sart[,list(SEX = 'Total', AGEGP = 'Total', NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM')]
tot2 <- sart[,list(AGEGP = 'Total', NO_ART = length(unique(STUDY_ID[ART==F]))), by = c('COMM','SEX')]
art <- rbind(art,tot1,tot2)

art[, COMM:= factor(COMM,levels=c('Total','inland','fishing'),labels=c('Total','Inland','Fishing'))]
art[, SEX:= factor(SEX,levels=c('Total','F','M'),labels=c('Total','Female','Male'))]
art[, AGEGP:= factor(AGEGP,levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))]



##############################

# NOTE: this is what I did, rest was from Alex.

################################

# I think i can focus in inland? 
# probably try to do as much with rincp? maybe not: probs all participants?
# arti is a subset of rincp

# keep age within 15-49
sequ <- as.data.table(readRDS(file.characteristics_sequenced_ind_R14_18))
sequ <- sequ[AGEYRS > 14 & AGEYRS < 50]
sequ[,  EVER.SEQ:=T]
sequ.ever <- subset(sequ,select=c('STUDY_ID', 'EVER.SEQ')) |> 
    unique()

# get unsuppressed at time x / not on ART at time x
## TODO I GET AN ERROR HERE
stopifnot( arti[, all(HIV == "P")] )
tmp <- arti[ COMM=='inland' & AGEYRS %between% c(15, 50)]
tmp <- merge(x=tmp, y=sequ.ever, by='STUDY_ID', all.x=TRUE)

arti[ART == FALSE, unique(ROUND)]
arti[, table(ART),by=ROUND]

tmp[is.na(EVER.SEQ), EVER.SEQ := FALSE ]
tmp <- tmp[
    ART == FALSE, 
    list( 
        N_EVERSEQ = sum(EVER.SEQ == TRUE), 
        N_TOT = .N ),
    by=c('ROUND', 'SEX', 'AGEGP')]

tmp[, (c('P_EVERSEQ', 'CL', 'CU')) := Hmisc::binconf(x=N_EVERSEQ, n=N_TOT, return.df=TRUE) ]
tmp[, ROUND_LAB := gsub("R0", "Round ", ROUND)]
tmp[, SEX_LAB := fifelse(SEX=="F", "Women", "Men")]

# Make Figure

naturemed_reqs()
p_everseq_givenart <- ggplot(tmp, aes(x=ROUND_LAB, color=SEX_LAB, pch=AGEGP, linetype=AGEGP)) + 
    geom_point(aes(y=P_EVERSEQ), position=position_dodge(width=.6) ) +
    geom_linerange(aes(ymin=CL, ymax=CU), position=position_dodge(width =.6)) +
    # facet_grid(ROUND ~ .) +
    scale_y_continuous(limits = c(0,1), expand=c(0,0),labels=scales::percent) +
    scale_color_manual(values=c(Women="#F4B5BD", Men="#85D4E3" )) + 
    labs( 
        x="Interview round",
        y="Proportion of ART negative participants eventually deep-sequenced", 
        linetype = "Age group", 
        pch = "Age group", 
        color = "Gender"
    ) +
    theme_bw() + 
    theme(legend.position = "bottom", strip.background = element_blank()) + 
    reqs

file.name <- file.path(outdir, "prop_arteverdeepseq_byroundagesex.pdf") 
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  ggsave_nature(filename=file.name, p=p_everseq_givenart, w=15, h=12)
}else{
  cat("File:", file.name, "already exists...\n")
}



# write table

# get totals
tab <- cube(tmp,
    j=list( N_TOT=sum(N_TOT), N_EVERSEQ=sum(N_EVERSEQ) ),
    by=c('ROUND_LAB','SEX_LAB', 'AGEGP')) |>
    subset(!is.na(ROUND_LAB))
tab[is.na(SEX_LAB), SEX_LAB:='Both' ]
tab[is.na(AGEGP), AGEGP:='All' ]
tab[, P_LAB := scales::label_percent(accuracy = .01)(N_EVERSEQ/N_TOT)]

setkey(tab, ROUND_LAB, SEX_LAB, AGEGP)
setcolorder(tab, c("ROUND_LAB", "SEX_LAB", "AGEGP",  "N_TOT", "N_EVERSEQ","P_LAB") )

fixed_nms <- c('ROUND_LAB', "SEX_LAB")
nms_toshift <- setdiff(names(tab), fixed_nms)
shifted_idx <- c(4, 1, 2, 3)
n_groups <- nrow(tab)/4
new_idx <- rep(1:n_groups-1, times=rep(4, n_groups))*4
new_idx <- new_idx + shifted_idx
tab <- cbind(
    tab |> subset(select=fixed_nms),
    tab[new_idx] |> subset(select= nms_toshift)
)


tab2 <- copy(tab)
cols <- c("ROUND_LAB", "SEX_LAB")
tab2[, (cols) := lapply(.SD, function(vec) {
    idx <- vec == shift(vec)
    vec[idx == TRUE] <- ""
    vec
}), .SDcols = cols]
rownames(tab2) <- NULL


file.name <- file.path(outdir, 'table_proportionART_eventuallysequenced.tex')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  write.to.tex(tab2, file=file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}


# if we want to split them
file.name <- file.path(outdir, 'table_proportionART_eventuallysequenced_part1.tex')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  write.to.tex(tab2[1:(.N/2)],file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}

file.name <- file.path(outdir, 'table_proportionART_eventuallysequenced_part2.tex')
if(! file.exists(file.name) | config$overwrite.existing.files )
{
  cat("Saving file:", file.name, '\n')
  write.to.tex(tab2[((.N/2)+1):.N],file.name)
}else{
  cat("File:", file.name, "already exists...\n")
}

