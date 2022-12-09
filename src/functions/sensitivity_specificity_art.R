require(xtable)
require(stringr)

get.dtable <- function(DT)
{
    tmp <- DT[COMM=='inland' & ! is.na(ART) & ! is.na(VLNS)] 
    tmp[, AGE_GROUP := '35-49']
    tmp[AGEYRS < 35, AGE_GROUP := '25-34']
    tmp[AGEYRS < 25,AGE_GROUP := '15-24' ]

    dtable <- tmp[, .N , by=c('ART', 'VLNS', 'ROUND', 'SEX', 'AGE_GROUP')]
    dtable <- dtable |> dcast(ROUND + SEX + AGE_GROUP ~ ART + VLNS)
    dtable[is.na(dtable)] <- 0


    # rename
    names(dtable) <- gsub( 'TRUE', 'ARV_REPORT', names(dtable))
    names(dtable) <- gsub( 'FALSE', 'ARV_NAIVE', names(dtable))
    names(dtable) <- gsub( '0', 'SUPP', names(dtable))
    names(dtable) <- gsub( '1', 'NOTSUPP', names(dtable))

    dtable
}

find_sensitivity_specificity_art <- function(rprev)
{
    if(0) # all of this is in the table
    {
        # find percentage of participant who did not report art but had not viremic viral load
        tmp <- rprev[COMM == 'inland' & ART == F & !is.na(VLNS), 
                     list(X = length(STUDY_ID[VLNS == 0]), 
                          N = length(STUDY_ID)), by = 'ROUND']
        tmp[, PROP := round(X / N * 100, 2)]
        tmp

        # find percentage of participant who report art and had not viremic viral load
        tmp <- rprev[COMM == 'inland' & ART == T & !is.na(VLNS),
                     list(X = length(STUDY_ID[VLNS == 0]), 
                          N = length(STUDY_ID)), by = 'ROUND']

        tmp[, PROP := round(X / N * 100, 2)]
        tmp

        if(0)
        { # plot

            df <- copy(rprev)
            df[, AGE_GROUP := '35-49']
            df[AGEYRS < 35, AGE_GROUP := '25-34']
            df[AGEYRS < 25,AGE_GROUP := '15-24' ]

            # find percentage of participant who did not report art but had not viremic viral load
            tmp <- df[COMM == 'inland' & ART == F & !is.na(VLNS),
                      list(X = length(STUDY_ID[VLNS == 0]), 
                           N = length(STUDY_ID)), by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp[, PROP := X / N ]
            tmp[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower, by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper, by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp1 <- df[COMM == 'inland' & ART == F & !is.na(VLNS),
                       list(X = length(STUDY_ID[VLNS == 0]), 
                            N = length(STUDY_ID)), by = 'ROUND']

            tmp1[, PROP := X / N ]
            tmp1[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower, by = c('ROUND')]
            tmp1[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper, by = c('ROUND')]

            ggplot(tmp, aes(x = SEX, fill = AGE_GROUP)) + 
                geom_hline(data = tmp1, aes(yintercept= PROP), col = 'grey30') +
                geom_hline(data = tmp1, aes(yintercept= CL), linetype = 'dashed', col = 'grey30') +
                geom_hline(data = tmp1, aes(yintercept= CU), linetype = 'dashed', col = 'grey30') +
                geom_bar(aes(y =PROP),stat = 'identity', position = position_dodge()) + 
                geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(width = 0.9), width = 0.2) + 
                theme_bw() + 
                facet_grid(ROUND~.) + 
                labs(x = 'sex', y = 'p(suppressed == T | art == no)') + 
                scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01)))

            # find percentage of participant who report art and had not viremic viral load
            tmp1 <- df[COMM == 'inland' & ART == T & !is.na(VLNS),
                       list(X = length(STUDY_ID[VLNS == 0]), 
                            N = length(STUDY_ID)),
                        by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp1[, PROP := X / N ]
            tmp1[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower,
                 by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp1[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper,
                 by = c('ROUND', 'AGE_GROUP', 'SEX')]
            tmp <- df[COMM == 'inland' & ART == T & !is.na(VLNS),
                      list(X = length(STUDY_ID[VLNS == 0]), 
                           N = length(STUDY_ID)), by = 'ROUND']
            tmp[, PROP := X / N ]
            tmp[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower,
                by = c('ROUND')]
            tmp[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper,
                by = c('ROUND')]

            ggplot(tmp1, aes(x = SEX, fill = AGE_GROUP)) + 
                geom_hline(data = tmp, aes(yintercept= PROP), col = 'grey30') +
                geom_hline(data = tmp, aes(yintercept= CL), linetype = 'dashed', col = 'grey30') +
                geom_hline(data = tmp, aes(yintercept= CU), linetype = 'dashed', col = 'grey30') +
                geom_bar(aes(y =PROP),stat = 'identity', position = position_dodge()) + 
                geom_errorbar(aes(ymin = CL, ymax = CU), position = position_dodge(width = 0.9), width = 0.2) + 
                theme_bw() + 
                facet_grid(ROUND~.) + 
                labs(x = 'age group', y = 'p(suppressed == T | art == yes)') + 
                scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01)))
        }
    }

    tmp <- get.dtable(rprev)
    .agr.coull <- function(yes,no)
    {
        agr.coull <- binom::binom.confint(yes, yes+no, methods = 'agresti-coull')
        agr.coull$lower[which(agr.coull$lower < 0)] <- 0
        agr.coull$upper[which(agr.coull$upper > 1)] <- 1
        agr.coull
    }

    cbind(tmp, tmp[, {
        z.spec <- .agr.coull(ARV_REPORT_SUPP, ARV_REPORT_NOTSUPP)
        z.sens <- .agr.coull(ARV_NAIVE_NOTSUPP, ARV_NAIVE_SUPP)
        list(
             SPEC_M=z.spec$mean,
             SPEC_U=z.spec$upper,
             SPEC_L=z.spec$lower,
             SENS_M=z.sens$mean,
             SENS_U=z.sens$upper,
             SENS_L=z.sens$lower
        )
    }]) -> tmp
    tmp
}

make_table_sensitivity_specificity_art <- function(rprev, break.lines=FALSE)
{
    # Helpers
    .agr.coull.format <- function(yes,no, newline=FALSE)
    {
        agr.coull <- binom::binom.confint(yes, yes+no, methods = 'agresti-coull')

        .r <- function(x) 
        {
            x[which(x < 0)] <- 0
            x[which(x > 1)] <- 1
            format(round(100*x, 1), nsmall=1)
        }

        if(newline==TRUE)
        {
            out <- with(agr.coull, 
                        sprintf("%s%%\n[%s-%s]", .r(mean), .r(lower), .r(upper)))
        }else{
            out <- with(agr.coull, 
                        sprintf("%s%% [%s-%s]", .r(mean), .r(lower), .r(upper)))
        }
        out
    }

    # get confusion matrix and totals
    dtable <- get.dtable(rprev)
    cols <- names(dtable)[names(dtable) %like% 'ARV']
    tmp <- dtable[, c(AGE_GROUP='ALL', lapply(.SD, sum)) 
                  , .SDcols=cols,
                  by=c('ROUND', 'SEX')]
    tmp1 <- tmp[, c(SEX='ALL', AGE_GROUP='ALL', lapply(.SD, sum)) 
                , .SDcols=cols, by='ROUND']
    dtable <- rbind(dtable, tmp, tmp1)
    # reorder rows so that Totals are at the top
    dtable[, AGE_GROUP := paste0('Z', AGE_GROUP)]
    dtable[, AGE_GROUP := gsub('ZA', 'A', AGE_GROUP)]
    setkey(dtable, ROUND, SEX, AGE_GROUP)
    dtable[, AGE_GROUP := gsub('Z', '', AGE_GROUP)]

    # get agresti coull intervals for sensitivity and specificity
    dtable[, `:=` (
        SPECIFICITY = .agr.coull.format(ARV_REPORT_SUPP,
                                        ARV_REPORT_NOTSUPP,
                                        newline=TRUE),
        SENSITIVITY = .agr.coull.format(ARV_NAIVE_NOTSUPP,
                                        ARV_NAIVE_SUPP,
                                        newline=TRUE)
    ), ]
    dtable[SPECIFICITY %like% 'NaN', SPECIFICITY := '']
    dtable[SENSITIVITY %like% 'NaN', SENSITIVITY := '']

    if(break.lines)
    {   # split CIs in second row
        cols <- setdiff(names(dtable), c('SPECIFICITY', 'SENSITIVITY'))
        .f <- function(x) stringr::str_split(x, '\n') |> unlist()
        dtable <- dtable[, .(
            SENSITIVITY = .f(SENSITIVITY),
            SPECIFICITY = .f(SPECIFICITY)
        ) , by=cols] 
    }

    dtable[ AGE_GROUP=='ALL' & SEX=='ALL', LABEL := 'Total']
    dtable[ AGE_GROUP=='ALL' & SEX=='F', LABEL := 'Female']
    dtable[ AGE_GROUP=='ALL' & SEX=='M', LABEL := 'Male']
    dtable[ AGE_GROUP!='ALL', LABEL := AGE_GROUP]

    if(break.lines)
    {
        second.entry.to.na <- function(x){
            x[ (1:length(x)) %%2 == 0]  <- NA; x
        }
        cols <- c("SEX", "AGE_GROUP", 
                  "ARV_NAIVE_SUPP","ARV_NAIVE_NOTSUPP",
                  "ARV_REPORT_SUPP","ARV_REPORT_NOTSUPP")
        dtable[, (cols) := lapply(.SD, second.entry.to.na) , .SDcols=cols]
    }

    # add empty rows only containing label 'Age' after the 'Female'/'Male' labels
    empty_row  <- dtable[1, lapply(.SD, function(x) NA) , .SDcols=names(dtable)]
    empty_row[, LABEL:='Age']
    idx <- dtable[, LABEL %like% 'ale']
    idx <- c(0, cumsum(idx[-length(idx)])) 
    if(break.lines)
        idx <- idx %/% 2
    dtable[, DUMMY := idx]
    dtable <- dtable[, rbind(.SD, empty_row), by=DUMMY]
    dtable <- dtable[1:(.N-1)]
    dtable[, DUMMY := NULL]

    # fix rounds after transform
    rnds <- dtable[!is.na(ROUND), unique(ROUND)]
    dtable[is.na(ROUND), ROUND:=rep(rnds, each=.N/length(rnds))]

    # fix labels now:
    idx <- which(dtable$LABEL == shift(dtable$LABEL, 1))
    dtable[idx, LABEL := NA]

    # integer as integers
    cols <- names(dtable)[names(dtable) %like% 'ARV']
    dtable[, (cols):=lapply(.SD, as.integer) , .SDcols=cols]
    
    make.table.for.round <- function(i)
    {
        .gs <- function(x)
        {
            # makes correct align for Age levels
            x <- gsub('(Age)', '\\\\\\multicolumn\\{1\\}\\{r|\\}\\{Age\\}', x) 
            x <- gsub('(15-24)', '\\\\\\multicolumn\\{1\\}\\{r|\\}\\{15-24\\}', x)  
            x <- gsub('(25-34)', '\\\\\\multicolumn\\{1\\}\\{r|\\}\\{25-34\\}', x)  
            x <- gsub('(35-49)', '\\\\\\multicolumn\\{1\\}\\{r|\\}\\{35-49\\}', x)  

            # puts CI on two different lines 95%\n[80-100]
            if(break.lines)
            {
                x <- gsub("([0-9]*.[0-9]\\\\\\%)",
                          "\\\\shortstack{\\1\\\\\\\\", x)
                x <- gsub("(\\[.*\\])",
                          "\\1}", x)
            }
            x
        }

        table.cols <- c("LABEL", "ARV_NAIVE_SUPP", "ARV_NAIVE_NOTSUPP", "ARV_REPORT_SUPP", "ARV_REPORT_NOTSUPP", "SPECIFICITY", "SENSITIVITY")
        SD <- dtable[ROUND == i, .SD, .SDcols=table.cols]


        filename <- paste0('~/Downloads/ART_supp_sensitivity_table',i,'.tex')
        SD |> 
            xtable() |>
            print(include.rownames=FALSE,
                  include.colnames=FALSE,
                  only.contents=TRUE,
                  booktabs=FALSE,
                  comment=FALSE,
                  type='latex',
                  file=filename,
                  hline.after = NULL,
                  compress=FALSE
                  ) -> x

        writeLines(text=.gs(readLines(filename)), filename)

    }

    tex.table <- lapply(rnds, make.table.for.round)
    # filename <- file.path('~/Downloads/ART_supp_sensitivity_table.rds')
    # saveRDS(tex.table, filename)
}
