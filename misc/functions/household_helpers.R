.p <- function(x)
        round(100*x, 2)

get.marriage.from.quest <- function(DT)
{
        cols <- grep('mar', names(DT), value=T)
        cols <- c( 'study_id', 'round', cols)
        # tmp <- DT[, ..cols]

        .translate.numbers2bool <- function(x)
        {
                y <- fcase(
                           x == 1, TRUE,
                           x == 2, FALSE,
                           x %in% c(0, 8), NA
                )
                y
        }

        .translate.polymar <- function(x)
        {
                # assuming "a few (01-02)" means 2 wives
                # assuming "a lot (03+)" means 3 wives
                x <- as.integer(x)
                x[ x== 92] <- 2L
                x[ x== 93] <- 3L
                x[ x== 98] <- NA_integer_

                cat('There were', sum(x > 10, na.rm=TRUE),
                    '(', mean(x > 10, na.rm=TRUE),') times in which more than 10 wives were reported \n')
                x
        }

        cols0 <- grep('marr', names(DT), value=T)
        DT[, (cols0) := lapply(.SD, .translate.numbers2bool) , .SDcols=cols0]
        DT[, polymar := .translate.polymar(polymar)]

        DT
}

preprocess.flow <- function(DT)
{

        # Fix missing study_id issues and impute via curr_id
        # __________________________________________________
        
        flow[,
                cat("There are ",  .p(mean(study_id == '')), '% NA study_ids\n')
             , ]

        cat("Use curr_id to attribute some of these...\n ")
        stopifnot(flow[study_id != '', uniqueN(curr_id), by='study_id'][, all(V1)])

        tmp <- unique(flow[, .(study_id, curr_id)])
        idx <- tmp[, .N , by=curr_id][N == 2, curr_id ]
        dict <- tmp[ curr_id %in% idx & study_id != '', ]

        impute_currid2studyid <- function(DT){
                DT <- merge(DT, dict, by='curr_id', all.x=TRUE)
                DT[ !is.na(study_id.y), study_id.x := study_id.y ]
                DT[, study_id.y := NULL ]
                setnames(DT, 'study_id.x', 'study_id')
                DT <- unique(DT)
                DT
        }

        flow <- impute_currid2studyid(flow)

        # For attributed study_id's appearing multiple times at same round,
        # throw away those that were not attributed
        # (using id as marker for previous study_id version)
        idx <- flow[study_id != '', uniqueN(.SD), by=c('study_id', 'round')][V1 > 1, study_id]
        rbind(
                flow[!(study_id %in% idx), ],
                flow[study_id %in% idx & id != '', ]
        ) -> flow

        # check consistency after study_id attribution
        stopifnot(flow[study_id != '', uniqueN(.SD), by=c('study_id', 'round')][, all(V1==1)])
        stopifnot(
                flow[ study_id != '' , uniqueN(sex), by='study_id'][, all(V1)],
                flow[ study_id != '' , uniqueN(birthdat), by='study_id'][, all(V1)],
                # flow[, list( uniqueN(round), .N ), by='study_id']
                flow[study_id != '', uniqueN(round == .N) , by='study_id'][, all(V1)]
        )

        # Extract columns of interest
        cols <- c(grep('^study_id|round|sex|hh|member_num', names(flow), value=T))
        flow <- flow[study_id != '', ..cols]
        flow[ , study_id := paste0('RK-', study_id) ]
        flow
}


