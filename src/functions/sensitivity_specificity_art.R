require(xtable)
require(stringr)

get.dtable <- function(DT, quest)
{
    tmp <- DT[COMM=='inland' & ! is.na(ART) & ! is.na(VLNS)] 
    
    # use age from quest
    set(tmp, NULL, 'AGEYRS', NULL)
    tmp <- merge(tmp, quest, by.x = c('ROUND', 'STUDY_ID'), by.y = c('round', 'study_id'))
    setnames(tmp, 'ageyrs', 'AGEYRS')
    
    # find age group
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


find_sensitivity_specificity_art <- function(rprev, quest)
{ 

    tmp <- get.dtable(rprev, quest)
    
    # save XX needed for paper
    .r <- function(x) formatC(x*100, digits=2, format='f')
    tmp[, 
        sum(ARV_REPORT_SUPP)/ (sum(ARV_REPORT_NOTSUPP) + sum(ARV_REPORT_SUPP))
        ] |> .r() -> aggregated_sensitivity
    filename <- file.path('~/Downloads', 'arv_sensitivity_suppression_aggregate.rds')
    cat('Saving aggregated sensitivity for Sexpr in',filename,'\n')
    saveRDS(aggregated_sensitivity, filename)

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
    
    #
    # extrapolate linearly
    
    # age to predict
    AGEYRSPREDICT = 15:49
    
    # take midpoont of age group
    tmp[, AGEYRS := mean(c(as.numeric(gsub('(.+)-.*', '\\1', AGE_GROUP)), as.numeric(gsub('.*-(.+)', '\\1', AGE_GROUP)))), by = 'AGE_GROUP']
    
    # plug value of round 16 for missing value in round 15
    tmp[ROUND == 'R015' & SEX == "M" & AGE_GROUP == '15-24', SPEC_M := tmp[ROUND == 'R016' & SEX == "M" & AGE_GROUP == '15-24', SPEC_M]]
    
    # take linear interpolation
    tmp1 <- tmp[, {
        
        linearapprox.SENS_M <- approx(AGEYRS, SENS_M, xout=AGEYRSPREDICT)$y
        linearapprox.SENS_M[1:(which.min(is.na(linearapprox.SENS_M)) - 1)] = SENS_M[1]
        linearapprox.SENS_M[which.max(is.na(linearapprox.SENS_M)):length(linearapprox.SENS_M)] = SENS_M[length(SENS_M)]
        
        linearapprox.SPEC_M  <- approx(AGEYRS, SPEC_M , xout=AGEYRSPREDICT)$y
        linearapprox.SPEC_M [1:(which.min(is.na(linearapprox.SPEC_M )) - 1)] = SPEC_M [1]
        linearapprox.SPEC_M [which.max(is.na(linearapprox.SPEC_M )):length(linearapprox.SPEC_M )] = SPEC_M [length(SPEC_M )]
        
        list(AGEYRS = AGEYRSPREDICT, SENS_M = linearapprox.SENS_M, SPEC_M = linearapprox.SPEC_M)
    }, by = c('ROUND', 'SEX')]
    
    if(0){ # plot

        ggplot(tmp) + 
            geom_point(aes(x = AGEYRS, y =1-SENS_M), alpha = 0.4) + 
            geom_line(data = tmp1, aes(x = AGEYRS, y = 1-SENS_M, col= 'linear')) + 
            theme_bw() + 
            facet_grid(ROUND~SEX) + 
            labs(x = 'age group', y = 'p(suppressed == T | art == no)') + 
            scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01)))
        ggsave('~/Downloads/sens_linear_new.png', w = 6, h = 8)
        
        ggplot(tmp) + 
            geom_point(aes(x = AGEYRS, y =SPEC_M), alpha = 0.4) + 
            geom_line(data = tmp1, aes(x = AGEYRS, y = SPEC_M, col= 'linear')) + 
            theme_bw() + 
            facet_grid(ROUND~SEX) + 
            labs(x = 'age group', y = 'p(suppressed == T | art == yes)') + 
            scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0, 0.01)))
        ggsave('~/Downloads/spe_linear_new.png', w = 6, h = 8)
        
    }
    
    return(tmp1)
}


make_table_sensitivity_specificity_art <- function(rprev, quest, break.lines=FALSE)
{

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
    dtable <- get.dtable(rprev, quest)
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
