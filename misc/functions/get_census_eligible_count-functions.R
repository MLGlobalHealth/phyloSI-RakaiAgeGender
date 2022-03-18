hivc.db.Date2numeric<- function( x )
{
  if(!class(x)%in%c('Date','character'))	return( x )
  x	<- as.POSIXlt(x)
  tmp	<- x$year + 1900
  x	<- tmp + round( x$yday / ifelse((tmp%%4==0 & tmp%%100!=0) | tmp%%400==0,366,365), d=3 )
  x
}

RakaiFull.phylogeography.220310.data.eligibility.participation.sequenced<- function(infile)
{
  require(data.table)
  
  # infile					<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/rakai_elibility.rda"
  # outfile.base			<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/"

  
  load(infile)
  
  #	subset to data of interest
  de	<- as.data.table(eldat)
  de	<- subset(de, status%in%c('_Participated','Away','Blood refusal','Missing data','Other','Refused','urine sample'))
  de	<- subset(de, visit<17)
  de	<- subset(de, select=c(PERM_ID, CURR_ID, visit, RESIDENT, MOBILITY, REGION, COMM_NUM, HH_NUM, SEX, STUDY_ID, status, date, AGEYRS, inmigrant))
  setnames(de, colnames(de), toupper(colnames(de)))
  #	define PARTICIPATED as "participated, missing data ok"
  #	TODO take out missing data
  de[, PARTICIPATED:= as.integer(STATUS%in%c('_Participated'))]
  de[, DATE:= hivc.db.Date2numeric(DATE)]
  
  #	merge communities that are very close / identical
  setnames(de, 'COMM_NUM', 'COMM_NUM_RAW')
  set(de, NULL, 'COMM_NUM', de[, gsub('^107$|^16$','16m',gsub('^776$|^51$','51m',gsub('^4$|^24$','24m',gsub('^1$|^22$','22m',COMM_NUM_RAW))))])
  
  #	fixup missing permanent IDs with Study IDs
  stopifnot( nrow(subset(de, is.na(PERM_ID) & is.na(STUDY_ID)))==0 )
  tmp		<- de[, which(is.na(PERM_ID))]
  set(de, tmp, 'PERM_ID', de[tmp, STUDY_ID])
  #	find study participants with multiple permanent IDs and fixup
  tmp		<- subset(de, !is.na(STUDY_ID))[, list(N_PERM=length(unique(PERM_ID))), by='STUDY_ID']
  tmp		<- subset(tmp, N_PERM>1, STUDY_ID)
  tmp[, DUMMY:=1]
  de		<- merge(de, tmp, by='STUDY_ID', all.x=1)
  tmp		<- de[, which(DUMMY==1)]
  set(de, tmp, 'PERM_ID', de[tmp, STUDY_ID])
  de[, DUMMY:=NULL]
  
  #	we care about wether individuals participated in any of the visits 15, 15.1, 16
  #	the best we can do is look at permanent ids, even though they are not unique
  tmp		<- de[, list(PARTICIPATED_ANY_VISIT=as.integer(any(PARTICIPATED==1))), by='PERM_ID']
  de		<- merge(de, tmp, by='PERM_ID')
  #	to those of ever participated, find study id
  tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
  stopifnot( nrow(subset(tmp, PARTICIPATED==1 & is.na(STUDY_ID)))==0 )
  tmp		<- tmp[, list(DUMMY= STUDY_ID[PARTICIPATED==1][1]), by='PERM_ID']
  de		<- merge(de, tmp, by='PERM_ID', all.x=TRUE)
  tmp		<- de[, which(is.na(STUDY_ID) & !is.na(DUMMY))]
  set(de, tmp, 'STUDY_ID', de[tmp,DUMMY])
  de[, DUMMY:=NULL]
  setnames(de, 'STUDY_ID','RID')
  #	there are a few individuals with Study ID who did not "participate" when missing not included
  #subset(de, PARTICIPATED_ANY_VISIT==0 & !is.na(STUDY_ID))
  
  return(de)

  # #	prepare HIV status
  # tmp		<- RakaiCirc.epi.get.info.170208()
  # rd		<- tmp$rd
  # rneuro	<- tmp$rn
  # ra		<- tmp$ra
  # #	prepare self report ART
  # #	find those who report to be on ART at their first visit in 15-17
  # rart	<- merge(ra, ra[VISIT>=14 & VISIT<17, list(VISIT=min(VISIT)), by='RID'], by=c('RID','VISIT'))
  # set(rart, rart[, which(is.na(ARVMED))], 'ARVMED', 2L)
  # setnames(rart, 'ARVMED', 'SELFREPORTART_AT_FIRST_VISIT')
  # rart	<- subset(rart, select=c(RID, SELFREPORTART_AT_FIRST_VISIT))
  # #set(rart, NULL,  'SELFREPORTART', rart[, factor(SELFREPORTART, levels=c(0,1,2), labels=c('no','yes','unknown'))])
  # #	select meta-data closest to first pos date
  # tmp		<- rd[, list(VISIT=VISIT[which.min(abs(DATE-FIRSTPOSDATE))]), by='RID']
  # tmp2	<- rd[, list(PANGEA=as.integer(any(PANGEA==1))), by='RID']
  # rd		<- merge(rd, tmp, by=c('RID','VISIT'))
  # rd[, HIV:= 1L]
  # rd[, PANGEA:=NULL]
  # rd		<- merge(rd, tmp2, by='RID')
  # #
  # ra		<- unique(subset(ra, !is.na(FIRSTPOSDATE), select=c(RID, SEX, VISIT, VISDATE, FIRSTPOSDATE, ARVMED, COMM_NUM, COMM_NUM_A)))
  # tmp		<- ra[, list(VISIT=VISIT[which.min(abs(VISDATE-FIRSTPOSDATE))]), by='RID']
  # ra		<- merge(ra, tmp, by=c('RID','VISIT'))
  # ra[, HIV_1517:= 1L]
  # #rd		<- subset(rd, BIRTHYR>2010-50 & (is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>2010)))
  # 
  # #	add HIV status
  # tmp		<- rd[FIRSTPOSDATE<2015+2/12, list(HIV=max(HIV), PANGEA=max(PANGEA)), by='RID']
  # de		<- merge(de, tmp, by='RID', all.x=1)
  # tmp		<- unique(subset(ra, FIRSTPOSDATE<2015+2/12, select=c(RID, HIV_1517, FIRSTPOSDATE)))
  # de		<- merge(de, tmp, by='RID', all.x=1)
  # 
  # #	add ART status
  # de		<- merge(de, rart, by='RID', all.x=1)
  # 
  # #	get individuals with at least 750nt overlap with another individual at 20X
  # infile	<- "~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/full_run/close_pairs_170704_cl35_withmedian.rda"
  # load(infile)
  # rn	<- subset(rn, NEFF>3)
  # rn	<- data.table(RID=rn[, unique(c(ID1, ID2))],MIN_PNG_OUTPUT=1)
  # #	add to rn individuals with any sequence data
  # load('~/Dropbox (SPH Imperial College)/2015_PANGEA_DualPairsFromFastQIVA/RakaiAll_input_170301/Rakai_phyloscanner_170301_b75_part2.rda')
  # tmp		<- unique(subset(dc, !is.na(SID), select=RID))
  # tmp[, BAM_OUTPUT:=1L]
  # rn		<- merge(tmp, rn, all=1, by='RID')
  # set(rn, rn[, which(is.na(MIN_PNG_OUTPUT))],'MIN_PNG_OUTPUT',0L)
  # #	add
  # de		<- merge(de, rn, by='RID', all.x=1)
  # rd		<- merge(rd, rn, by='RID', all.x=1)
  # rneuro	<- merge(rneuro, rn, by='RID', all.x=1)
  
  #	tmp above has 4074 individuals
  #	of those, 3758 are in 'de', ie participated in rounds 15-17
  
  #	for every individual that ever participated, keep first record of when participated
  #	for every individual that never participated, keep last record
  tmp		<- subset(de, PARTICIPATED_ANY_VISIT==0)[, list(VISIT=max(VISIT)), by='PERM_ID']
  tmp		<- merge(de, tmp, by=c('PERM_ID','VISIT'))
  tmp2	<- subset(de, PARTICIPATED_ANY_VISIT==1)[, list(VISIT=min(VISIT[PARTICIPATED==1])), by='PERM_ID']
  de		<- merge(de, tmp2, by=c('PERM_ID','VISIT'))
  de		<- rbind(de, tmp)
  
  #	fill in missing HIV etc
  # tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(HIV))]
  # set(de, tmp, 'HIV', 0L)
  # tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(HIV_1517))]
  # set(de, tmp, 'HIV_1517', 0L)
  # tmp		<- de[, which(PARTICIPATED_ANY_VISIT==1 & is.na(BAM_OUTPUT))]
  # set(de, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
  # tmp		<- rd[, which(is.na(BAM_OUTPUT))]
  # set(rd, tmp, c('BAM_OUTPUT','MIN_PNG_OUTPUT'), 0L)
  # 
  # de[, table(HIV, HIV_1517)]
  # 672 with HIV_1517==1 and HIV==0
  #	0 with HIV_1517==0 and HIV==1
  # use HIV_1517 below
  
  # rneuro[, sum(MIN_PNG_OUTPUT)]			# 224
  # rd[, sum(MIN_PNG_OUTPUT)]				# 2746
  # de[, sum(MIN_PNG_OUTPUT, na.rm=TRUE)]	# 2689
  
  #	calculate number for whom we have deep seq output
  # tmp		<- subset(rd, MIN_PNG_OUTPUT==1)[, list(DEEP_SEQ=length(RID)), by=c('COMM_NUM')]
  # tmp		<- subset(tmp, DEEP_SEQ>0)
  # tmp[, DEEP_SEQ:=NULL]
  # tmp[, COMM_ANY_MIN_PNG_OUTPUT:= 1]
  # de		<- merge(de, tmp, by='COMM_NUM', all.x=1)
  # rd		<- merge(rd, tmp, by='COMM_NUM', all.x=1)
  # set(de, de[, which(is.na(COMM_ANY_MIN_PNG_OUTPUT))], 'COMM_ANY_MIN_PNG_OUTPUT', 0L)
  # set(rd, rd[, which(is.na(COMM_ANY_MIN_PNG_OUTPUT))], 'COMM_ANY_MIN_PNG_OUTPUT', 0L)
  # tmp		<- unique(subset(rd, select=c(COMM_NUM,COMM_NUM_A)))
  # de		<- merge(de, tmp, by='COMM_NUM')
  
  #	subset to only those communities with seq data
  #	4 fishing, 34 inland
  #def		<- copy(de)
  #de		<- subset(de, !is.na(DUMMY))
  #de[, DUMMY:=NULL]
  #rd		<- subset(rd, !is.na(DUMMY))
  #rd[, DUMMY:=NULL]
  #
  
  #	add age category at midpoint of observation period, 2013.25
  set(de, de[, which(is.na(DATE))], 'DATE', hivc.db.Date2numeric(as.Date("2011-09-02")))
  de[, AGE_AT_MID:= 2013.25 - (hivc.db.Date2numeric(DATE)-AGEYRS)]
  de[, AGE_AT_MID_C:= cut(AGE_AT_MID, breaks=c(10,25,35,65), labels=c('15-24','25-34','35+'), right=FALSE)]
  stopifnot( nrow(subset(de, is.na(AGE_AT_MID_C)))==0 )

  #
  # #	prepare inmigrant -- identify inmigrants from fishing communities and from external
  # #
  # infile.migrant		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/RakaiPangeaMetaData_v2.rda"
  # infile.migrant2		<- "~/Dropbox (SPH Imperial College)/Rakai Pangea Meta Data/Data for Fish Analysis Working Group/migrants_withMissingGPS.csv"
  # load(infile.migrant)
  # inmigrant	<- as.data.table(inmigrant)
  # inmigrant2	<- as.data.table(read.csv(infile.migrant2))
  # #	plot fisherfolk	to figure out how much of a radius we need
  # if(0)
  # {
  #   zf		<- data.table(longitude=c(31.763,31.7968,31.754,31.838), latitude=c(-0.915, -0.6518, -0.703, -0.497), ID= c('Kasensero','Bukyanju','NearBwende','Fish4'))
  #   make_circles <- function(centers, radius, nPoints = 100){
  #     # centers: the data frame of centers with ID
  #     # radius: radius measured in kilometer
  #     #
  #     meanLat <- mean(centers$latitude)
  #     # length per longitude changes with lattitude, so need correction
  #     radiusLon <- radius /111 / cos(meanLat/57.3)
  #     radiusLat <- radius / 111
  #     circleDF <- data.frame(ID = rep(centers$ID, each = nPoints))
  #     angle <- seq(0,2*pi,length.out = nPoints)
  # 
  #     circleDF$lon <- unlist(lapply(centers$longitude, function(x) x + radius * cos(angle)))
  #     circleDF$lat <- unlist(lapply(centers$latitude, function(x) x + radius * sin(angle)))
  #     return(circleDF)
  #   }
  #   zc <- make_circles(zf, 0.01)
  #   ggmap(zm) +
  #     geom_point(data=zf, aes(x=longitude, y=latitude, pch=ID), stroke=1.5, alpha=0.8) +
  #     geom_polygon(data=zc, aes(lon, lat, group = ID), color = "red", alpha = 0)
  #   #	radius of length 0.01 should catch
  #   tmp		<- inmigrant[, list( 	DIST_KASENSERO= sqrt( (inmig_lon- 31.763)^2 + (inmig_lat - (-0.915))^2),
  #                              DIST_BUKYANJU= sqrt( (inmig_lon- 31.7968)^2 + (inmig_lat - (-0.6518))^2),
  #                              DIST_NEARBWENDE= sqrt( (inmig_lon- 31.754)^2 + (inmig_lat - (-0.703))^2),
  #                              DIST_FISH4= sqrt( (inmig_lon- 31.838)^2 + (inmig_lat - (-0.497))^2)
  #   ), by=c('RCCS_studyid','visit')]
  #   tmp		<- melt(tmp, id.vars=c('RCCS_studyid','visit'))
  #   ggplot(subset(tmp, value<0.3), aes(x=value)) +
  #     geom_histogram(binwidth=0.01) +
  #     facet_grid(variable~.)
  #   zfd		<- merge(inmigrant, subset(tmp, value<0.01, c(RCCS_studyid, visit)), by=c('RCCS_studyid','visit'))
  # }
  # #	so fishing sites are MALEMBO DIMU KASENSERO NAMIREMBE but there are spelling mistakes
  # #	clean up inmigrant
  # #
  # #	inmigrant[, unique(sort(inmig_place))]
  # set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub('DDIMO|DDIMU|DIMO|DIMU','DIMU',inmig_place)])
  # set(inmigrant, inmigrant[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
  # set(inmigrant, NULL, 'inmig_place', inmigrant[, gsub("KASEMSERO","KASENSERO",inmig_place)])
  # set(inmigrant, inmigrant[, which(grepl('KASENSERO',inmig_place))], 'inmig_place', 'KASENSERO')
  # set(inmigrant2, inmigrant2[, which(grepl('MALEMBO',inmig_place))], 'inmig_place', 'MALEMBO')
  # #	define from_fishing and from_outside and from_inland
  # inmigrant[, INMIG_LOC:= 'inland' ]
  # set(inmigrant, inmigrant[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE|KYABASIMBA',inmig_place))], 'INMIG_LOC','fisherfolk')
  # set(inmigrant, inmigrant[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
  # set(inmigrant, inmigrant[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
  # set(inmigrant, inmigrant[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')
  # inmigrant2[, INMIG_LOC:= 'inland' ]
  # set(inmigrant2, inmigrant2[, which(grepl('MALEMBO|DIMU|KASENSERO|NAMIREMBE|KYABASIMBA',inmig_place))], 'INMIG_LOC','fisherfolk')
  # set(inmigrant2, inmigrant2[, which(inmig_admin0!='Uganda')], 'INMIG_LOC','external')
  # set(inmigrant2, inmigrant2[, which(inmig_admin1!='Rakai')], 'INMIG_LOC','external')
  # set(inmigrant2, inmigrant2[, which(is.na(inmig_admin1))], 'INMIG_LOC','unknown')
  # inmigrant2[, table(INMIG_LOC)]
  # #	  external fisherfolk     inland    unknown
  # #       1          1          6          4
  # #	  now resolved 8 more locations, great
  # setnames(inmigrant2, 'INMIG_LOC', 'INMIG_LOC2')
  # inmigrant2	<- subset(inmigrant2, select=c('RCCS_studyid','visit','INMIG_LOC2'))
  # inmigrant	<- merge(inmigrant, inmigrant2, by=c('RCCS_studyid','visit'), all.x=TRUE)
  # tmp			<- inmigrant[, which(!is.na(INMIG_LOC2) & INMIG_LOC2!='unknown')]
  # set(inmigrant, tmp, 'INMIG_LOC', inmigrant[tmp, INMIG_LOC2])
  # inmigrant[, table(INMIG_LOC)]
  # #	  external fisherfolk     inland    unknown
  # #	   763         41       1577        373
  # inmigrant[, INMIGRATE_DATE:= hivc.db.Date2numeric(date)]
  # set(inmigrant, NULL, c('date','inmigrant','INMIG_LOC2'), NULL)
  # setnames(inmigrant, colnames(inmigrant), gsub('RCCS_STUDYID','RID',toupper(colnames(inmigrant))))
  # #
  # #	inmigrants done
  # #
  # 
  # 
  # #
  # #	add inmigrants to eligibility / participation
  # #
  # di			<- subset(inmigrant, select=c(RID, INMIGRATE_DATE, INMIG_LOC))
  # di			<- merge(di, de, by='RID', all.x=TRUE)
  # #	ignore inmigrants not seen in 15-16
  # di			<- subset(di, !is.na(STATUS))
  # #	only inmigration date closest to first visit and before the first visit
  # tmp			<- di[, list(FIRSTVISITDATE=min(DATE)), by='RID']
  # di			<- merge(di, tmp, by='RID')
  # tmp			<- di[INMIGRATE_DATE<=FIRSTVISITDATE, list(INMIGRATE_DATE= INMIGRATE_DATE[which.min(FIRSTVISITDATE-INMIGRATE_DATE)]), by='RID']
  # di			<- merge(di, tmp, by=c('RID','INMIGRATE_DATE'))
  # di			<- unique(di, by=c('RID','VISIT'))
  # #	select those inmigrated in the last 2 years
  # tmp			<- subset(di, (FIRSTVISITDATE-INMIGRATE_DATE)<=2, c(RID, INMIG_LOC))
  # setnames(tmp, 'INMIG_LOC', 'INMIG_2YRS_LOC')
  # tmp[, INMIG_2YRS:=1]
  # de			<- merge(de, tmp, by='RID', all.x=TRUE)
  # set(de, de[, which(is.na(INMIG_2YRS))], 'INMIG_2YRS', 0)
  # set(de, de[, which(is.na(INMIG_2YRS_LOC))], 'INMIG_2YRS_LOC', 'resident')
  # 
  # #	some fixup
  # set(de, de[, which(PARTICIPATED_ANY_VISIT==0 & MIN_PNG_OUTPUT==1)], 'PARTICIPATED_ANY_VISIT', 1L)
  # set(de, de[, which(SELFREPORTART_AT_FIRST_VISIT==1 & MIN_PNG_OUTPUT==1)], 'SELFREPORTART_AT_FIRST_VISIT', 0L)


  # de[, table(HIV_1517, SELFREPORTART_AT_FIRST_VISIT, useNA='if')]
  #	        SELFREPORTART
  #	HIV_1517     0         1     2  <NA>
  #		    0    20715     8     7    10
  #			1     3884  1256     2     0
  #			<NA>     0     0     0 11763

  # #	now calculate participation by community and gender
  # des		<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','PARTICIPATED_ANY_VISIT')]
  # set(des, NULL, 'PARTICIPATED_ANY_VISIT', des[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
  # des		<- dcast.data.table(des, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX~PARTICIPATED_ANY_VISIT, value.var='N')

  #	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
  #	among those that participated in rounds 15 - 16
  # tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
  # tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)),
  #                     HIV_1516_NO=length(which(HIV_1517==0)),
  #                     SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
  #                     DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]
  # des		<- merge(des, tmp, by=c('COMM_NUM','SEX'), all=TRUE)

  #	add community type
  # des[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
  
  #	now do the same by community and gender and age category
  # desa	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','AGE_AT_MID_C','PARTICIPATED_ANY_VISIT')]
  # set(desa, NULL, 'PARTICIPATED_ANY_VISIT', desa[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
  # desa	<- dcast.data.table(desa, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX+AGE_AT_MID_C~PARTICIPATED_ANY_VISIT, value.var='N')
  # tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
  # tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)),
  #                     HIV_1516_NO=length(which(HIV_1517==0)),
  #                     SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
  #                     DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX','AGE_AT_MID_C')]
  # desa	<- merge(desa, tmp, by=c('COMM_NUM','SEX','AGE_AT_MID_C'), all=TRUE)
  # desa[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
  # #	check for missing and replace by 0
  # for(x in c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','SLART_AT_FIRST_VISIT','DEEP_SEQ_1516'))
  #   set(desa, which(is.na(desa[[x]])), x, 0)
  
  #	now do the same by community and gender and migration category
  # desm	<- de[, list(N=length(CURR_ID)), by=c('COMM_NUM','COMM_NUM_A','COMM_ANY_MIN_PNG_OUTPUT','SEX','INMIG_2YRS_LOC','PARTICIPATED_ANY_VISIT')]
  # set(desm, NULL, 'PARTICIPATED_ANY_VISIT', desm[, factor(PARTICIPATED_ANY_VISIT, levels=c('1','0'), labels=c('PART_EVER','PART_NEVER'))])
  # desm	<- dcast.data.table(desm, COMM_NUM+COMM_NUM_A+COMM_ANY_MIN_PNG_OUTPUT+SEX+INMIG_2YRS_LOC~PARTICIPATED_ANY_VISIT, value.var='N')
  # tmp		<- subset(de, PARTICIPATED_ANY_VISIT==1)
  # tmp		<- tmp[, list(	HIV_1516_YES=length(which(HIV_1517==1)),
  #                     HIV_1516_NO=length(which(HIV_1517==0)),
  #                     SLART_AT_FIRST_VISIT=length(which(SELFREPORTART_AT_FIRST_VISIT==1)),
  #                     DEEP_SEQ_1516=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX','INMIG_2YRS_LOC')]
  # desm	<- merge(desm, tmp, by=c('COMM_NUM','SEX','INMIG_2YRS_LOC'), all=TRUE)
  # desm[, COMM_TYPE:= as.character(factor(substr(COMM_NUM_A,1,1)=='f',levels=c(TRUE,FALSE),labels=c('fisherfolk','inland')))]
  # #	check for missing and replace by 0
  # for(x in c('PART_EVER','PART_NEVER','HIV_1516_YES','HIV_1516_NO','SLART_AT_FIRST_VISIT','DEEP_SEQ_1516'))
  #   set(desm, which(is.na(desm[[x]])), x, 0)
  
  
  
  #	now calculate HIV pos by end of round 16, and sequenced by end of round 16, by community and gender
  #	among those that ever participated
  #	excluding those that died before 2010, or that reached age 60 before 2010
  # tmp		<- subset(rd, is.na(EST_DATEDIED) | (!is.na(EST_DATEDIED) & EST_DATEDIED>=2010))
  # tmp		<- subset(tmp, (2010-BIRTHYR)<60)
  # rds		<- tmp[, list(HIV_EVER_YES=length(which(HIV==1)), DEEP_SEQ_EVER=length(which(MIN_PNG_OUTPUT==1)) ), by=c('COMM_NUM','SEX')]
  # rds[, sum(DEEP_SEQ_EVER)]
  # # 	2700 (among those in rd who are not too old and did not die before 2010)
  # 
  # save(des, desa, desm, rds, de, rd, di, inmigrant, file='~/Dropbox (SPH Imperial College)/Rakai Fish Analysis/180322_sampling_by_gender_age.rda')
  return(de)
}
