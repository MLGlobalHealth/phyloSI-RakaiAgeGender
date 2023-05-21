.int2date <- function(x) {
    as.Date(x, origin = "1970-01-01")
}

`%larger%` <- function(a, b) {
    if (length(a) != 2 | length(b) != 2) {
        stop("Error in %larger%: a and b should both be 2-dimensional")
    }
    diff(a) > diff(b)
}

plot.rectangles <- function(DT, values = c("", "TSI", "INTERSECT"), idx = data.table()) {
    .p <- function(...) {
        out <- paste(..., sep = ".")
        out <- gsub("^\\.", "", out)
    }
    .c <- function(x) {
        fcase(x == "", "blue", x == "TSI", "red", x == "INTERSECT", "orange")
    }
    .s <- "SOURCE"
    .r <- "RECIPIENT"

    plot_dt <- double.merge(chain, DT)
    setkey(plot_dt, SOURCE, RECIPIENT)
    setkey(idx, SOURCE, RECIPIENT)
    values[.p(values, "MIN", .s) %in% names(DT)]

    if (nrow(idx)) {
        plot_dt <- plot_dt[idx]
    }


    add.rect <- function(pre) {
        geom_rect(aes_string(
            xmin = .p(pre, "MIN", .s),
            xmax = .p(pre, "MAX", .s),
            ymin = .p(pre, "MIN", .r),
            ymax = .p(pre, "MAX", .r)
        ), fill = NA, color = .c(pre))
    }

    g <- ggplot(plot_dt)

    for (v in values) {
        g <- g + add.rect(v)
    }

    g <- g +
        geom_abline(aes(intercept = 0, slope = 1), linetype = "dashed", color = "black") +
        theme_bw() +
        labs(x = "DOI source", y = "DOI recipient")
    g
}

check.range.consistency <- function(drange)
{
    inconsistent <- TRUE
    while(inconsistent)
    {
        inconsistent <- FALSE
        check <- merge(chain, drange[, .(SOURCE=AID, MIN.SOURCE=MIN, MAX.SOURCE=MAX)], by='SOURCE')
        check <- merge(check, drange[, .(RECIPIENT=AID, MIN.RECIPIENT=MIN, MAX.RECIPIENT=MAX)], by='RECIPIENT')
        if( ! args$confidential ) 
        {
            cat('\n(randomized data: artificially fix min-max dates of infection)\n')
            idx <- check[ MAX.RECIPIENT < MIN.SOURCE,  ]
            cat(idx[, .N], '\n')
            if(idx[, .N > 0])
                inconsistent <- TRUE
            idx[ MAX.RECIPIENT < MIN.SOURCE, MIN.SOURCE := MAX.RECIPIENT - 365 ]
            idx <- idx[, list(MIN=min(MIN.SOURCE), MAX=unique(MAX.SOURCE)), by=SOURCE]
            setnames(idx, 'SOURCE', 'AID')
            drange <- rbind(drange[!AID %in% idx$AID], idx)
            stopifnot(drange[, ! any(is.na(MIN) | is.na(MAX)) ])
        }
    }
    check[ MAX.RECIPIENT < MIN.SOURCE, stopifnot(.N==0)]
    return(drange)
}

get.communities.where.participated <- function() {
    # For each AID reports the community type(s) where each participant had participated
    meta_env <- new.env()
    load(path.meta, envir = meta_env)
    cols <- c("aid", "comm", "round", "sample_date", "sex")
    dcomms <- subset(meta_env$meta_data, select = cols)
    names(dcomms) <- toupper(names(dcomms))
    dcomms <- dcomms[!is.na(AID),
        .(COMM = paste0(sort(unique(COMM)), collapse = "-"), uniqueN(COMM), SEX = unique(SEX)),
        by = "AID"
    ]
    dcomms[, uniqueN(AID) == .N]
    dcomms
}

table.pairs.in.inland <- function(DT) {
    double.merge(DT, dcomms[, .(AID, COMM, SEX)], by_col = "AID")[
        COMM.SOURCE %like% "inland" & COMM.RECIPIENT %like% "inland"
    ][, {
        cat(.N, "\n")
        table(SEX.SOURCE, SEX.RECIPIENT)
    }]
}

build.phylo.network.from.pairs <- function(path.chains = path.chains.data) {
    # prepare helpers
    dcomms <<- get.communities.where.participated()

    # get couples and chains from run
    cat("\n Building Network from pairs...\n\n")

    chains_env <- new.env()
    load(path.chains, envir = chains_env)
    dpl <- setDT(chains_env$dpl)
    dc <- setDT(chains_env$dc)
    serohistory_impact <- summarise_serohistory_impact_on_pairs(dc)
    stopifnot(dpl[, .N, by = c("H1", "H2")][, all(N == 1)])
    stopifnot(dpl[, all(H1 < H2)])

    # Merge linkage and direction score and subset to linkage score > threshold
    idx <- dpl[SCORE > threshold.likely.connected.pairs, .(H1, H2, SCORE)]
    dir_scores <- dc[CATEGORISATION %like% "close.and.adjacent.and.directed.cat.sero", .(H1, H2, TYPE, SCORE)]
    dir_scores <- dir_scores[,
        {
            z <- which(SCORE > threshold.direction)
            list(TYPE = TYPE[z], SCORE_DIR = SCORE[z])
        },
        by = c("H1", "H2")
    ]
    dlinkdir <- merge(idx, dir_scores, by = c("H1", "H2"), all.x = TRUE)
    stopifnot(dir_scores[, .N, by = c("H1", "H2")][, all(N == 1)])
    stopifnot(dlinkdir[, .N, by = c("H1", "H2")][, all(N == 1)])

    # get range of dates
    drange <- get.infection.range.from.testing()

    # if add unsupported but with strong SCORE_DIR
    if (get.sero.extra.pairs) {
        sero_extra_pairs <- dpl[SCORE <= threshold.likely.connected.pairs, .(H1, H2, SCORE)]
        sero_extra_pairs <- merge(sero_extra_pairs, drange[, .(H1 = AID, MIN.1 = MIN, MAX.1 = MAX)], by = "H1")
        sero_extra_pairs <- merge(sero_extra_pairs, drange[, .(H2 = AID, MIN.2 = MIN, MAX.2 = MAX)], by = "H2")
        sero_extra_pairs <- rbind(
            sero_extra_pairs[MAX.1 < MIN.2, `:=`(SCORE_DIR = 1, TYPE = "12")],
            sero_extra_pairs[MAX.2 < MIN.1, `:=`(SCORE_DIR = 1, TYPE = "21")]
        )[!is.na(TYPE), .(H1, H2, SCORE, SCORE_DIR, TYPE), ]
        sero_extra_pairs <- unique(sero_extra_pairs)
        dlinkdir <- rbind(dlinkdir, sero_extra_pairs)
    }

    # check direction is consistent with the serohistory (btw we subset to RCCS here)
    # WE CAN ONLY REALLY DO THIS IF THE DATA ARE NOT RANDOMIZED...
    if (basename(path.chains) %like% "randomized") {
        dlinkdir <- merge(dlinkdir, drange[, .(H1 = AID, MIN.1 = MIN, MAX.1 = MAX)], by = "H1")
        dlinkdir <- merge(dlinkdir, drange[, .(H2 = AID, MIN.2 = MIN, MAX.2 = MAX)], by = "H2")
        dlinkdir[MAX.1 < MIN.2 & is.na(TYPE), `:=`(TYPE = "12", SCORE_DIR = 1)]
        dlinkdir[MAX.2 < MIN.1 & is.na(TYPE), `:=`(TYPE = "21", SCORE_DIR = 1)]
        dlinkdir[MAX.1 < MIN.2, stopifnot(all(TYPE == "12"))]
        dlinkdir[MAX.2 < MIN.1, stopifnot(all(TYPE == "21"))]
        dlinkdir[, `:=`(MAX.1 = NULL, MAX.2 = NULL, MIN.1 = NULL, MIN.2 = NULL)]
    }

    # Assign source and recipient labels
    dnewpairs <- rbind(
        dlinkdir[TYPE == "12", .(SOURCE = H1, RECIPIENT = H2, SCORE, SCORE_DIR)],
        dlinkdir[TYPE == "21", .(SOURCE = H2, RECIPIENT = H1, SCORE, SCORE_DIR)]
    ) |> unique()
    stopifnot(dnewpairs[, .N, by = c("SOURCE", "RECIPIENT")][, all(N == 1)])

    dnewpairs[, cat(
        "There are ", .N, "pairs with score>threshold and",
        sum(SCORE < threshold.likely.connected.pairs),
        "pairs with only one possible direction of tranmsission\n"
    )]

    if (exists("overleaf_expr")) {
        # Pairs involving inland participants
        # table.pair.in.inland(dnewpairs)
        double.merge(dnewpairs, dcomms[, .(AID, COMM, SEX)], by_col = "AID") |>
            subset(COMM.RECIPIENT %like% "inland" | COMM.SOURCE %like% "inland") -> tmp_overleaf
        tmp_overleaf[, table(COMM.SOURCE, COMM.RECIPIENT)]
        overleaf_expr[["N_pairs_1_in_study"]] <- tmp_overleaf[, .N]

        # now get how many have direction swapped
        rbind(
            serohistory_impact[, .(SOURCE = H1, RECIPIENT = H2, CHANGED_DIR)],
            serohistory_impact[, .(RECIPIENT = H1, SOURCE = H2, CHANGED_DIR)]
        ) |>
            merge(x = tmp_overleaf, all.x = TRUE, by = c("SOURCE", "RECIPIENT")) |>
            with(sum(CHANGED_DIR)) -> tmp_overleaf2
        overleaf_expr[["N_swapped"]] <- tmp_overleaf2

        # number of unique recipients
        overleaf_expr[["N_unique_rec"]] <- tmp_overleaf[, uniqueN(RECIPIENT)]
    }

    # Remove same-sex pairs
    if (!postpone.samesex.removal) {
        dsex <- meta[, .(aid, sex)] |> unique()
        idx <- merge(dnewpairs, dsex[, .(SOURCE = aid, SEX.SOURCE = sex)], all.x = T, by = "SOURCE")
        idx <- merge(idx, dsex[, .(RECIPIENT = aid, SEX.RECIPIENT = sex)], all.x = T, by = "RECIPIENT")
        # idx[, table(SEX.SOURCE,SEX.RECIPIENT, useNA = 'ifany'),]
        dhomosexualpairs <- idx[(SEX.SOURCE == SEX.RECIPIENT), ]
        dhomosexualpairs[, IDCLU := NULL]
        idx <- idx[!(SEX.SOURCE == SEX.RECIPIENT), .(SOURCE, RECIPIENT)]

        tmp <- nrow(dnewpairs) - nrow(idx)
        cat("Excluding", tmp, "of", nrow(dnewpairs), "pairs of homosexual or unknown sex\n")
        setkey(dnewpairs, SOURCE, RECIPIENT)
        setkey(idx, SOURCE, RECIPIENT)
        dnewpairs <- dnewpairs[idx]
        cat(nrow(dnewpairs), "pairs remaining\n")
    }

    # Solve multiple sources by assigning recipient with strongest linkage support (same for homosexual pairs)
    dnewpairs <- dnewpairs[,
        {
            z <- which.max(SCORE)
            list(SOURCE = SOURCE[z], SCORE = SCORE[z], SCORE_DIR = SCORE_DIR[z])
        },
        by = "RECIPIENT"
    ]


    if (exists("overleaf_expr")) {
        # Pairs involving at both inland participants
        tmp <- table.pairs.in.inland(dnewpairs) |> as.vector()
        names(tmp) <- c("FF", "MF", "FM", "MM")
        overleaf_expr[["N_unique_rec_participants"]] <- sum(tmp)
        overleaf_expr[["N_unique_rec_participants_hetero"]] <- tmp["MF"] + tmp["FM"]
        overleaf_expr[["N_unique_rec_participants_dirs"]] <- tmp
    }

    # remove homosexual pairs
    if (!postpone.samesex.removal) {
        # solve multiple sources for same-sex too
        dhomosexualpairs <- dnewpairs[!RECIPIENT %in% dnewpairs$RECIPIENT,
            {
                z <- which.max(SCORE)
                list(SOURCE = SOURCE[z], SCORE = SCORE[z], SCORE_DIR = SCORE_DIR[z])
            },
            by = "RECIPIENT"
        ]
    } else {
        cat("Postponing removal of same sex pairs after solving for double sources.\n")
        dsex <- meta[, .(aid, sex)] |> unique()
        idx <- merge(dnewpairs, dsex[, .(SOURCE = aid, SEX.SOURCE = sex)], all.x = T, by = "SOURCE")
        idx <- merge(idx, dsex[, .(RECIPIENT = aid, SEX.RECIPIENT = sex)], all.x = T, by = "RECIPIENT")
        # idx[, table(SEX.SOURCE,SEX.RECIPIENT, useNA = 'ifany'),]
        dhomosexualpairs <- idx[(SEX.SOURCE == SEX.RECIPIENT), ]
        dhomosexualpairs[, IDCLU := NULL]
        idx <- idx[!(SEX.SOURCE == SEX.RECIPIENT), .(SOURCE, RECIPIENT)]

        tmp <- nrow(dnewpairs) - nrow(idx)
        cat("Excluding", tmp, "of", nrow(dnewpairs), "pairs of homosexual or unknown sex\n")
        setkey(dnewpairs, SOURCE, RECIPIENT)
        setkey(idx, SOURCE, RECIPIENT)
        dnewpairs <- dnewpairs[idx]
        cat(nrow(dnewpairs), "pairs remaining\n")
    }

    # return to the environment the data.tables of interest.
    if (get.sero.extra.pairs) 
        sero_extra_pairs <<- sero_extra_pairs
    
    dhomosexualpairs <<- dhomosexualpairs

    if (exists("overleaf_expr"))
        overleaf_expr <<- overleaf_expr

    if (exists("serohistory_impact"))
        serohistory_impact <<- serohistory_impact 

    setcolorder(dnewpairs, c("SOURCE", "RECIPIENT"))
    idclus <- dnewpairs |>
        graph_from_data_frame() |>
        components()
    idclus <- data.table(RECIPIENT = names(idclus$membership), IDCLU = unname(idclus$membership))
    dnewpairs <- merge(dnewpairs, idclus, by = "RECIPIENT")
    setcolorder(dnewpairs, c("SOURCE", "RECIPIENT"))

    return(dnewpairs)
}

summarise_serohistory_impact_on_pairs <- function(DC, categ = "close.and.adjacent.and.directed") {
    # DC <- copy(dc); categ = 'close.and.adjacent.and.directed'

    # check that both cat and cat.sero categories are available.
    categ <- gsub(".sero", "", categ)
    stopifnot(DC[, unique(CATEGORISATION) %like% categ |> sum() == 2])

    # subset of interest
    dc_cat <- DC[CATEGORISATION %like% categ]
    setkey(dc_cat, H1, H2)
    cat("in how many cases were the scores changed?\n")
    dc_changed <- dc_cat[, uniqueN(SCORE) > 1, by = c("H1", "H2", "TYPE")]
    dc_changed <- dc_changed[, .(CHANGED_SCORE = any(V1)), by = c("H1", "H2")]
    dc_changed[, table(CHANGED_SCORE)] |>
        knitr::kable() |>
        print()

    cat("in how many cases were the directions changed?\n")
    idx <- dc_changed[CHANGED_SCORE == TRUE, .(H1, H2)]
    dc_changed_dir <- dc_cat[idx]
    # dc_changed_dir[is.na(SCORE)]

    # for each pair and class type, evaluate whether score in dir 1->2 is larger
    dc_changed_dir <- dc_changed_dir[,
        {
            dir12 <- TYPE %like% "12"
            dir21 <- TYPE %like% "21"
            list(DIR12 = SCORE[dir12] >= SCORE[dir21])
        },
        by = c("H1", "H2", "CATEGORISATION")
    ]

    # for each pair, look whether categorisations disagree.
    dc_changed_dir <- dc_changed_dir[, .(CHANGED_DIR = uniqueN(DIR12) == 2), by = c("H1", "H2")]
    dc_changed_dir[, table(CHANGED_DIR)]

    # now let's output a DT with the pairs, indicating whether
    # - serohistory changed scores
    # - serohistory changed directions!
    .f <- function(DT1, DT2) merge(DT1, DT2, all.x = TRUE, all.y = TRUE)
    out <- list(unique(dc_cat[, .(H1, H2)]), dc_changed, dc_changed_dir) |> Reduce(f = .f)
    out[is.na(CHANGED_DIR), CHANGED_DIR := FALSE]
    out[, table(CHANGED_SCORE, CHANGED_DIR, useNA = "ifany")]

    return(out)
}

get.extra.pairs.from.serohistory <- function(DCHAIN, META) {
    tmp <- DCHAIN[SCORE_LINKED <= threshold.likely.connected.pairs, .(H1, H2)]
    tmp <- merge(tmp,
        META[, .(
            H1 = aid,
            SEX = sex,
            DFP = date_first_positive,
            DLN = date_last_negative
        )],
        by = "H1"
    )
    tmp <- merge(tmp,
        META[, .(
            H2 = aid,
            SEX = sex,
            DFP = date_first_positive,
            DLN = date_last_negative
        )],
        by = "H2"
    )
    names(tmp) <- gsub(".x", ".H1", names(tmp))
    names(tmp) <- gsub(".y", ".H2", names(tmp))
    tmp <- tmp[DLN.H2 >= DFP.H1 | DLN.H1 >= DFP.H2]
    tmp[DLN.H2 >= DFP.H1, `:=`(EST_DIR = "12", SOURCE = H1, RECIPIENT = H2)]
    tmp[DLN.H1 >= DFP.H2, `:=`(EST_DIR = "21", SOURCE = H2, RECIPIENT = H1)]
    cat("There are", tmp[, .N], "pairs with only one possible direction implied by serohistory\n")
    tmp1 <- tmp[SEX.H1 != SEX.H2]
    cat(tmp1[, .N], "of those are among heterosexuals\n")

    idx <- tmp[, .(H1, H2, EST_DIR, SOURCE, RECIPIENT)]
    additional_pairs_from_serohistory <- merge(DCHAIN, idx, by = c("H1", "H2"))
    additional_pairs_from_serohistory[, `:=`(H1 = NULL, H2 = NULL)]
    additional_pairs_from_serohistory <<- additional_pairs_from_serohistory

    chain <- rbind(chain, additional_pairs_from_serohistory)
    unique(chain)
}

update.ranges.with.tsi.estimates <- function() {
    cols <- c("RENAME_ID", "pred_doi_min", "pred_doi_mid", "pred_doi_max", "visit_dt")
    dtsi <- fread(path.tsiestimates, select = cols)
    dtsi[, `:=`(AID = gsub("-fq[0-9]", "", RENAME_ID), RENAME_ID = NULL)]

    # extend short TSI ranges to the left so that they are 1 yr wide
    dtsi[, LENGTH_TSI := pred_doi_max - pred_doi_min, ]
    dtsi[LENGTH_TSI < 365, pred_doi_min := pred_doi_min - (365 - LENGTH_TSI)]
    dtsi[, LENGTH_TSI := pred_doi_max - pred_doi_min, ]

    # 80% of midpoint predictions are inside the plausible range
    initial_avg_length <- drange[, mean(MAX - MIN)]

    tmp <- merge(drange, dtsi, all.x = TRUE)
    tmp[, LENGTH_DEMO := MAX - MIN]
    tmp[, mean(pred_doi_mid >= MIN & pred_doi_mid <= MAX)]
    # 15 (1.5%) have non-intersecting range.
    # Make new intervals 1 yr wide towards direction pointed by TSI
    tmp[pred_doi_max < MIN, `:=`(NEW.MIN = MIN, NEW.MAX = MIN + pmin(365, LENGTH_DEMO))]
    tmp[pred_doi_min > MAX, `:=`(NEW.MIN = MAX - pmin(365, LENGTH_DEMO), NEW.MAX = MAX)]

    # for intersecting ranges, take intersection, and extend intersection to 1 yr if "allowed" by demographic
    tmp[!(pred_doi_min > MAX | MIN > pred_doi_max), `:=`(
        INTERSECT.MIN = pmax(MIN, pred_doi_min),
        INTERSECT.MAX = pmin(MAX, pred_doi_max)
    )]
    tmp[, LENGTH := as.integer(INTERSECT.MAX - INTERSECT.MIN)]
    tmp[is.na(NEW.MIN) & LENGTH < 365 & pred_doi_max > MAX, `:=`(
        NEW.MAX = MAX,
        NEW.MIN = MAX - pmin(365, LENGTH_DEMO)
    )]
    tmp[is.na(NEW.MIN) & LENGTH < 365 & pred_doi_min < MIN, `:=`(
        NEW.MAX = MIN + pmax(365, LENGTH_DEMO),
        NEW.MIN = MIN
    )]
    tmp[, NEW.MIN := fcoalesce(NEW.MIN, INTERSECT.MIN, MIN)]
    tmp[, NEW.MAX := fcoalesce(NEW.MAX, INTERSECT.MAX, MIN)]
    # tmp[, .(LENGTH_DEMO, NEW.MAX - NEW.MIN) ][V2 < 365]

    # only keep those consistent with direction of transmission
    tmp <- tmp[, .(AID, MIN = NEW.MIN, MAX = NEW.MAX)]
    tmp1 <- double.merge(chain[, .(SOURCE, RECIPIENT)], tmp, by_col = "AID")

    idx <- tmp1[MAX.RECIPIENT < MIN.SOURCE, unique(c(SOURCE, RECIPIENT))]
    rbind(
        tmp[!AID %in% idx],
        drange[AID %in% idx]
    ) -> tmp
    tmp

    final_avg_length <- tmp[, mean(MAX - MIN)]

    initial_avg_length / 365.25
    final_avg_length / 365.25
    tmp
}

double.merge <- function(DT1, DT2, by_col = "AID") {
    .ps <- function(x) paste0(x, ".SOURCE")
    .pr <- function(x) paste0(x, ".RECIPIENT")

    tmp <- copy(DT2)
    setnames(tmp, names(tmp), .ps(names(tmp)))
    out <- merge(DT1, tmp, by.x = "SOURCE", by.y = .ps(by_col))

    tmp <- copy(DT2)
    setnames(tmp, names(tmp), .pr(names(tmp)))
    out <- merge(out, tmp, by.x = "RECIPIENT", by.y = .pr(by_col))

    out
}

get.infection.range.from.testing <- function() 
{
    # get maximum and minimum dates

    # for missing 1st pos, use date collection instead (also store contradicting first pos)
    tmp <- meta[, get.sample.collection.dates(select_aid=aid, get_first_visit = TRUE)]
    meta <- merge(tmp, meta)

    # get 15th birthdate and check whether anyone was 100% infected previously
    meta[, date15th := date_birth + as.integer(365.25 * 15)]
    stopifnot(meta[date15th > date_first_positive, .N == 0])
    # don't min date more than 15 years prior first positive
    # meta[, mean(lowb15lastneg > date_last_negative), ]

    contradict_firstpos_datecoll <<- meta[date_collection < date_first_positive]
    drange <- meta[, .(
        AID = aid,
        MIN = pmax(date_last_negative, date15th, na.rm = TRUE),
        MAX = pmin(date_first_positive - 30, date_collection - 30, na.rm = TRUE)
    )]
    drange[, lowb15lastneg := MAX - as.integer(365.25 * 15)]
    drange[lowb15lastneg > MIN, MIN := lowb15lastneg]
    drange[, lowb15lastneg := NULL]
    drange
}

get.infection.range.from.tsi <- function(path, path.sequence.dates, chain_subset = TRUE) 
{
    seqdates <- readRDS(path.sequence.dates)
    names(seqdates) <- toupper(names(seqdates)) 

    if(path.tsiestimates %like% 'csv$'){
        dtsi <- fread(path.tsiestimates) 
        names(dtsi) <- toupper(names(dtsi))
    }else if(path.tsiestimates %like% 'rds$'){
        dtsi <- readRDS(path.tsiestimates)
    }
    dtsi <- merge(dtsi, seqdates, by=c('AID','RENAME_ID'))

    drange_tsi <- dtsi[, .(
        AID=AID,
        RENAME_ID=RENAME_ID,
        MIN=VISIT_DT - as.integer(365.25*RF_PRED_MAX_LINEAR),
        MAX=VISIT_DT - as.integer(365.25*RF_PRED_MIN_LINEAR)
    )]

    if (chain_subset) {
        idx <- chain[, unique(c(SOURCE, RECIPIENT))]
        drange_tsi <- drange_tsi[AID %in% idx, ]
    }
    drange_tsi
}

get.infection.range.from.tsi.deprecated <- function(path, chain_subset = TRUE) {
    cols <- c("RENAME_ID", "pred_doi_min", "pred_doi_max")
    drange_tsi <- fread(path, select = cols)
    setnames(drange_tsi, cols, c("AID", "MIN", "MAX"))
    drange_tsi[, `:=`(AID = gsub("-fq[0-9]", "", AID))]
    if (chain_subset) {
        idx <- chain[, unique(c(SOURCE, RECIPIENT))]
        drange_tsi <- drange_tsi[AID %in% idx, ]
    }
    drange_tsi
}

get.ancestors.from.chain <- function(CHAIN) {
    # start off with only recipients as
    out <- as.list(CHAIN$SOURCE)
    names(out) <- CHAIN$RECIPIENT

    while (TRUE) {
        .u <- unlist
        .add.one.level.ancestor <- function(x) {
            unique(c(.u(out[x]), .u(out[.u(out[x])])))
        }

        tmp <- lapply(names(out), .add.one.level.ancestor)
        names(tmp) <- names(out)
        tmp
        check <- length(unlist(tmp)) == length(unlist(out))
        out <- copy(tmp)
        if (check) break
    }
    sapply(out, length) |> table()
    out
}

check.inconsistent.testing <- function(DRANGE, switch_if_no_other_src = FALSE, remove_others_in_network = FALSE) {
    # DRANGE <- copy(drange_tsi); dancestors <- get.ancestors.from.chain(chain)
    chain_tmp <- double.merge(chain, DRANGE)

    # get the maximum earliest infection date across ancestors
    tmp <- lapply(dancestors, function(vec) DRANGE[AID %in% vec, max(MIN, na.rm = TRUE)])
    tmp <- data.table(RECIPIENT = names(tmp), MIN.ANCESTOR = as.Date(unlist(tmp), origin = "1970-01-01"))

    chain_tmp <- merge(chain_tmp, tmp, by = "RECIPIENT", all.x = TRUE)
    problematic_pair <- chain_tmp[MAX.RECIPIENT <= MIN.SOURCE]
    cat("There are", nrow(problematic_pair), "recipients whose maximum infection date preceeds the minium infection date of the source\n")
    problematic_net <- setdiff(chain_tmp[MAX.RECIPIENT <= MIN.ANCESTOR], problematic_pair)
    cat("Additionaly, there are", nrow(problematic_net), "recipients whose maximum infection date preceeds the the minimum infection date of one of its 'indirect' ancestors.\n")

    if (switch_if_no_other_src) {
        if (nrow(problematic_net)) {
            cat("\tremoving recipients inconsistent with 'deeper' ancestors...\n")
            setkey(problematic_net, SOURCE, RECIPIENT)
            setkey(chain, SOURCE, RECIPIENT)
            idx <- problematic_net[, .(SOURCE, RECIPIENT)]
            chain <- chain[!idx]
        }

        if (nrow(problematic_pair) != 0) {
            cat("\tswithing direction of tranmsission...\n")
            idx <- problematic_pair[, .(SOURCE, RECIPIENT)]
            setkey(chain, SOURCE, RECIPIENT)
            setkey(idx, SOURCE, RECIPIENT)
            # tmp <- chain[idx][, .(SOURCE=RECIPIENT,RECIPIENT=SOURCE, IDCLU)]
            idx

            # have to check that sources in pairs who'll swith to recipients don't already appear as recipients
            setkey(chain, SOURCE, RECIPIENT)
            cond <- !idx$SOURCE %in% chain$RECIPIENT
            idx_rm <- idx[!cond]
            idx <- idx[cond]

            chain <- chain[!idx_rm]

            chain_sub <- copy(chain[idx])
            setnames(chain_sub, c("SOURCE", "RECIPIENT"), c("RECIPIENT", "SOURCE"))
            setcolorder(chain_sub, c("SOURCE", "RECIPIENT"))

            chain <- rbind(chain[!idx], chain_sub)
        }
        return(chain)
    }

    if (!switch_if_no_other_src) {
        # reason beeing that there is no
        tmp <- rbind(problematic_pair, problematic_net)
        if (nrow(tmp) != 0) {
            cat(paste0(
                "Removing inconsistent ranges from DRANGE ",
                ifelse(remove_others_in_network, "and individuals in these networks", ""), " ...\n"
            ))
            nrow(tmp)
            idx <- tmp[, unique(c(SOURCE, RECIPIENT))]
            while (TRUE & remove_others_in_network) {
                l <- length(idx)
                idx <- chain[SOURCE %in% idx | RECIPIENT %in% idx, unique(c(SOURCE, RECIPIENT))]
                if (length(idx) == l) {
                    break
                }
            }
            cat("\t", DRANGE[AID %in% idx, .N], "individuals in total\n")
            DRANGE <- DRANGE[!AID %in% idx]
        }
        return(DRANGE)
    }
}


.shrink.based.on.source <- function(CHAIN, DRANGE) {
    rec <- CHAIN$RECIPIENT
    setkey(CHAIN, RECIPIENT)
    src <- CHAIN[rec, SOURCE]

    setkey(DRANGE, AID)
    min.source <- DRANGE[src, MIN]
    length(min.source)
    length(rec)
    nrow(DRANGE)
    DRANGE[rec, MIN.SOURCE := min.source]

    cond <- DRANGE[!is.na(MIN.SOURCE), sum(MIN < MIN.SOURCE)]
    if (cond) {
        DRANGE[, MIN := pmax(MIN, MIN.SOURCE, na.rm = TRUE), by = "AID"]
    }

    DRANGE[, MIN.SOURCE := NULL]

    return(list(C = cond, D = DRANGE))
}

.shrink.based.on.recipient <- function(CHAIN, DRANGE) {
    src <- unique(CHAIN$SOURCE)
    setkey(CHAIN, SOURCE)

    tmp <- merge(CHAIN, DRANGE[, .(RECIPIENT = AID, MAX)], by = "RECIPIENT")
    tmp <- tmp[, .(MAX.RECIPIENT = min(MAX)), by = "SOURCE"]

    # cols <- names(DRANGE)[names(DRANGE) %like% 'RECIPIENT']

    DRANGE <- merge(DRANGE, tmp, all.x = TRUE, by.x = "AID", by.y = "SOURCE")

    cond <- DRANGE[!is.na(MAX.RECIPIENT), sum(MAX > MAX.RECIPIENT)]
    if (cond) {
        DRANGE[, MAX := pmin(MAX, MAX.RECIPIENT, na.rm = TRUE), by = "AID"]
    }

    DRANGE[, MAX.RECIPIENT := NULL]

    return(list(C = cond, D = DRANGE))
}

shrink.intervals <- function(DT) {
    cat("Initial mean width:", DT[, round(as.numeric(mean(MAX - MIN) / 365.25), 2)], "\n")

    # This function return the avg int. length, but it should change DT (drange/drange_tsi in the background)
    cond <- TRUE
    iter <- 1
    while (cond) {
        cat(iter, "\n")
        # these functions modify drange, but should NOT modify dchain
        out0 <- .shrink.based.on.source(chain, DT)
        DT <- out0[["D"]]
        cond0 <- out0[["C"]]
        out1 <- .shrink.based.on.recipient(chain, DT)
        DT <- out1[["D"]]
        cond1 <- out1[["C"]]
        cat(cond0, cond1, "\n")
        cond <- cond0 | cond1
        iter <- iter + 1
    }
    cat("Final mean width:", DT[, round(as.numeric(mean(MAX - MIN) / 365.25), 2)], "\n")
    DT
}

intersect.serohistory.and.predictions.withrefinement <- function(DRANGE, DRANGE_TSI)
{
        # intersect the serohistory and tsi ranges at the individal-level. If intersection is empty, priorities tsi.
    tmp <- merge(DRANGE, DRANGE_TSI)
    tmp[,{
        rng_intersect <- (MAX>TSI.MIN & TSI.MAX>MIN)

        fifelse(rng_intersect,
                yes=pmax(MIN, TSI.MIN),
                no=MIN
        ) -> r_min
        fifelse(rng_intersect,
                yes=pmin(MAX, TSI.MAX),
                no=MAX
        ) -> r_max
        list(AID, RANGE_INTERSECT=rng_intersect, INTERSECT.MIN=r_min, INTERSECT.MAX=r_max)
        }] -> tmp1

    cat('For', tmp1[RANGE_INTERSECT==FALSE, .N], 'individuals the 2 ranges did not intersect\n')
    tmp1[, RANGE_INTERSECT := NULL]


    if(0)
    {   # VISUALISE
        for(i in 1:100)
        {
            idx <- chain[i, .(SOURCE,RECIPIENT)]
            p <- plot.rectangles(tmp, idx=idx, values = c("", "TSI"))
            print(p)
            Sys.sleep(1.5)
            if(i == 10)
                plot.rect <- copy(p)
        }
    }

    cols <- c('MIN', 'MAX')
    setnames(tmp1, paste0('INTERSECT.', cols), cols)
    DRANGE <- rbind( DRANGE[! (AID %in% tmp1$AID)],tmp1)

    return(DRANGE)
}

intersect.serohistory.and.predictions.withoutrefinement <- function(DRANGE, DRANGE_TSI)
{
    # keep coherent TSI predictions, and use serohistory if incoherencies
    aid_to_predict <- chain[, unique(c(SOURCE, RECIPIENT))]
    tmp_tsi <- DRANGE_TSI[ AID %in% aid_to_predict]
    tmp_ser <- DRANGE[AID %in% aid_to_predict]
    cat('In', sum(!aid_to_predict %in% tmp$AID), 'cases, incoherent predictions, using serohistory...\n')
    rbind(
          tmp_tsi[, .(AID, MIN=as.Date(TSI.MIN), MAX=as.Date(TSI.MAX))],
          tmp_ser[ ! (AID %in% tmp_tsi$AID) ]
    )  -> DRANGE

    idx <- double.merge(chain[, .(SOURCE,RECIPIENT)], 
                        DRANGE,
                        by_col = "AID")[MAX.RECIPIENT - MIN.SOURCE < 0,]
    cat('In', nrow(idx), 'pairs, the ranges were inconsistent, use serohistory for',
        idx[, uniqueN(c(SOURCE, RECIPIENT))],' participants instead...\n')
    idx <- idx[, c(SOURCE,RECIPIENT)]

    DRANGE <- rbind( DRANGE[! AID %in% idx], tmp_ser[AID %in% idx])
}

specify.relative.infectiousness <- function(ARGS)
{
    with(ARGS,
         data.table( START = c(0, RH.duration), VALUE = c(RH.infectiousness, 1))
    ) -> dinfectiousness
    dinfectiousness[, END := c(START[-1], Inf)]
    setcolorder(dinfectiousness, c('START', 'END', 'VALUE'))
    return(dinfectiousness)
}

set.mcmc.outputs.suffix <- function()
{
    suffix <- ''
    if( ! is.null(dinfectiousness))
        suffix <- paste0('_d', tmp, '_w', paste0( dinfectiousness$VALUE, collapse=''))
    # if( build.network.from.pairs  )
        suffix <- paste0(suffix, '_netfrompairs')
    if( threshold.likely.connected.pairs != .5 )
        suffix <- paste0(suffix, '_thr', gsub( '0\\.', '', threshold.likely.connected.pairs))
    if( threshold.direction != .5 )
        suffix <- paste0(suffix, '_dir', gsub( '0\\.', '', threshold.direction))
    if( get.sero.extra.pairs )
        suffix <- paste0(suffix, '_seropairs')
    if( args$sensitivity.no.refinement)
        suffix <- paste0(suffix, '_sensnoref')
    if( postpone.samesex.removal )
        suffix <- paste0(suffix, '_postponessrem')
    if( ! args$confidential)
        suffix <- paste0(suffix, '_randomized')

    cat('\n Chosen MCMC outputs suffix: ', suffix, '\n')
    return(suffix)
}

get.round.dates <- function(path=path.round.timeline)
{
    round_env <- new.env()
    load(path, envir = round_env)
    df_round_inland <- round_env$df_round_inland
    df_round_fishing <- round_env$df_round_fishing

    start_first_period_inland <- df_round_inland[round == 'R010', min_sample_date] # "2003-09-26"
    stop_first_period_inland <- df_round_inland[round == 'R015', max_sample_date] # "2013-07-05"
    start_second_period_inland <-df_round_inland[round == 'R016', min_sample_date] #  "2013-07-08"
    stop_second_period_inland <- df_round_inland[round == 'R018', max_sample_date] #  "2018-05-22"

    stopifnot(start_first_period_inland < stop_first_period_inland)
    stopifnot(stop_first_period_inland < start_second_period_inland)
    stopifnot(start_second_period_inland < stop_second_period_inland)

    cols <- c("ROUND", "MIN_SAMPLE_DATE", "MAX_SAMPLE_DATE", "INDEX_TIME", "INDEX_ROUND", "round")
    df_round_i <- make.df.round(df_round_inland) |> subset(select=cols)
    df_round_f <- make.df.round(df_round_fishing) |> subset(select=cols)

    cols2 <- cols[! cols %like% 'SAMPLE_DATE|INDEX']
    df_round_gi <- merge(df_round_i, df_round_f, by=cols2, all.x=TRUE)
    df_round_gi[, `:=` (
                        MIN_SAMPLE_DATE=pmin(MIN_SAMPLE_DATE.x, MIN_SAMPLE_DATE.y, na.rm=TRUE),
                        MAX_SAMPLE_DATE=pmax(MAX_SAMPLE_DATE.x, MAX_SAMPLE_DATE.y, na.rm=TRUE),
                        INDEX_TIME=INDEX_TIME.x,
                        INDEX_ROUND=INDEX_ROUND.x
                        )]
    df_round_gi <- subset(df_round_gi, select= ! names(df_round_gi) %like% '.x$|.y$' )
    df_round_gi <- subset(df_round_gi, select= cols)

    df_round_gi[, MAX_SAMPLE_DATE := pmax(MAX_SAMPLE_DATE, shift(MIN_SAMPLE_DATE, -1), na.rm=TRUE)]

    rm(round_env, df_round_inland, df_round_fishing, df_round_i, df_round_f)

    df_round_gi
}

classify.transmission.round.given.centroid <- function(idx) {
    # Requires definition of drounds in code
    .f <- function(r) {
        tmp <- idx %between% drounds[round == r, c(min_sample_date, max_sample_date)]
        out <- rep(NA_character_, length(idx))
        out[tmp] <- r
        out
    }
    drounds[, shift(min_sample_date, -1)]
    rs <- lapply(drounds$round, .f)
    rs <- Reduce("fcoalesce", rs)
}

group_nest_dt <- function(dt, ..., .key = "data") {
    stopifnot(is.data.table(dt))
    by <- substitute(list(...))
    dt <- dt[, list(list(.SD)), by = eval(by)]
    setnames(dt, old = "V1", new = .key)
    dt
}

unlist.centroids <- function(DT) {
    centroids <- DT[, unlist(CENTROIDS)]
    centroids <- data.table(study_id = names(centroids), doi_centroid = centroids)
    centroids[, doi_centroid := .int2date(doi_centroid)]
    centroids
}

# Think there might be a bug due to original choice not really being "unif"
plot.scatter.comparison.vs.uniform <- function(DT) {
    ggplot(data = DT) +
        geom_rect(data = drounds, aes(
            xmin = min_sample_date, xmax = max_sample_date, fill = round,
            ymin = min_sample_date, ymax = max_sample_date
        ), alpha = .5) +
        geom_point(aes(x = DOI.RECIPIENT, y = DOI2.RECIPIENT)) +
        geom_abline(slope = 1, linetype = "dashed", color = "red") +
        # geom_abline(slope=1, intercept=-365.25, linetype='dashed', color='grey50') +
        geom_hline(data = drounds, aes(yintercept = max_sample_date), linetype = "dotted") +
        geom_vline(data = drounds, aes(xintercept = max_sample_date), linetype = "dotted") +
        labs(
            title = "Comparison of attributed recipient doi",
            x = "DOI under uniform pdf", y = "DOI under custom pdf"
        ) +
        theme_bw() +
        theme() -> p
    filename <- file.path(out.dir, "scatter_uniform_weighted_centroids_recipient", suffix, ".png")
    ggsave(filename, p, width = 8, height = 8)
    p
}

compare.to.uniform <- function() {
    filename_net0 <- file.path(out.dir, "networks_GICentroids_d1_w11.rds")

    .f <- function(path) {
        tmp <- readRDS(path)
        tmp <- Reduce("c", tmp[, CENTROIDS])
        tmp <- unlist(tmp)
        tmp <- data.table(ID = names(tmp), DATE = .int2date(tmp))
        tmp
    }
    out <- lapply(c(filename_net, filename_net0), .f)
    out <- merge(out[[1]], out[[2]], by = "ID")
    setnames(out, c("DATE.x", "DATE.y"), c("DATE.custom", "DATE.uniform"))

    ggplot(out, aes(x = DATE.custom, DATE.uniform)) +
        geom_point() +
        theme_bw()
}
update.meta.pairs.after.doi.attribution <- function(path, outdir, overwrite = F) {
    if (!file.exists(path) | overwrite) {
        cat("\tRunning infection dates estimation...\n\n\n")
        output <- study_infection_pairs_dates(CHAIN = chain)
        saveRDS(output, path)
    } else {
        cat("\tLoading previously estimated dates of infection...\n")
        output <- readRDS(path)
    }

    # extract info
    dinfrange <- output$pairs
    dage <- output$age

    # plot
    g <- plot.phylopair.dates.scores(dinfrange,
        add.dots = TRUE,
        daterange.vars = c("MIN", "MAX"),
        doi.center.var = "DOI",
        title = "Midpoint/Network Attribution"
    )
    filename <- paste0(outdir, "phylopairs_dates_scores_after_networkattribution.png")
    ggsave(filename, g, w = 15, h = 15, units = "cm")

    # fix meta data
    aik2 <- copy(aik)
    names(aik2) <- tolower(names(aik2))

    dage[, AID := .pt2aid(study_id, aik2)]

    meta_data <- merge(meta_data[, -c("date_infection", "age_infection")],
        dage[, .(date_infection = DOI, age_infection = age_at_infection, aid = AID)],
        by = "aid", all.x = T
    )

    cols <- grep("^RECIPIENT$|^SOURCE$|SCORE|CLU", names(dinfrange), value = TRUE)
    chain <- dinfrange[, ..cols, ]
    cols <- c("SOURCE", "RECIPIENT")
    chain[, (cols) := lapply(.SD, function(x) .pt2aid(x, aik2)), .SDcols = cols, by = cols]

    pairs.all <- pairs.get.meta.data(chain, meta_data, aik)

    list(pairs.all = pairs.all, chain = chain, meta_data = meta_data)
}

study_infection_pairs_dates <- function(CHAIN) {
    cat("== Estimating network-coherent infection dates == \n")
    .is.mrc <- function(x) {
        grepl("^MRC", x)
    }

    .aid.2.studyid <- function(DT, by_var) {
        tmp <- merge(DT, aik, all.x = TRUE, by.x = by_var, by.y = "AID")
        stopifnot(tmp[, all(!is.na(PT_ID))])
        tmp[, (by_var) := NULL]
        setnames(tmp, "PT_ID", by_var)
        tmp
    }

    .double.merge <- function(DT1, DT2, by_var = "study_id") {
        cols <- setdiff(names(DT2), c("AID", "study_id"))
        cols.SOURCE <- paste0(cols, ".SOURCE")
        cols.RECIPIENT <- paste0(cols, ".RECIPIENT")
        tmp <- merge(DT1, DT2, by.x = "SOURCE", by.y = by_var, all.x = TRUE)
        setnames(tmp, cols, cols.SOURCE)
        tmp <- merge(tmp, DT2, by.x = "RECIPIENT", by.y = by_var, all.x = TRUE)
        setnames(tmp, cols, cols.RECIPIENT)
        tmp
    }

    # Get dates and ages data
    .select <- function(cols) {
        unique(meta_data[, ..cols])
    }

    ddates <- .select(c("study_id", "date_last_negative", "date_infection", "date_first_positive"))
    dage <- .select(c("study_id", "age_first_positive", "date_birth"))

    .check <- function(DT) {
        DT[, .N, by = "study_id"][, all(N == 1)]
    }
    stopifnot(.check(dage) & .check(ddates))

    dage <- merge(dage, aik, all.x = TRUE, by.x = "study_id", by.y = "PT_ID")
    ddates <- merge(ddates, aik, all.x = TRUE, by.x = "study_id", by.y = "PT_ID")

    chain2 <- copy(CHAIN)

    # remove uninteresting cols
    chain2[, CNTRL_ANY := CNTRL1 | CNTRL2]
    cols <- grep("CATEGORISATION|PTY_RUN|CNTRL[0-9]", names(chain2), value = TRUE)
    cols1 <- grep("^LINK|EST_DIR", names(chain2), value = TRUE)
    chain2[, (cols) := NULL]
    chain2[SCORE_DIR_12 == SCORE_DIR_21, cat("- There are ", .N, "pairs where direction of transmission is exactly equally likely\n")]
    # if(0) chain2[SCORE_DIR_12 > SCORE_DIR_21 , unique(.SD), .SDcols=cols1]
    chain2[, (cols1) := NULL]
    chain2[, all(round(SCORE_DIR_12 + SCORE_DIR_21, 2) == 1)]
    chain2[, SCORE_DIR_SR := pmax(SCORE_DIR_12, SCORE_DIR_21)]
    cols <- grep("SCORE_NW|SCORE_DIR_[0-9][0-9]", names(chain2), value = TRUE)
    chain2[, (cols) := NULL]

    # round entries for readability
    cols <- names(which(unlist(lapply(chain2, is.numeric))))
    chain2[, (cols) := lapply(.SD, function(x) round(x, 2)), .SDcols = cols]

    # Translate AIDs to study_ids
    chain2 <- .aid.2.studyid(chain2, by_var = "SOURCE")
    chain2 <- .aid.2.studyid(chain2, by_var = "RECIPIENT")

    # merge with dates
    chain2 <- .double.merge(DT1 = chain2, DT2 = ddates[, -"AID"])
    stopifnot(
        chain2[is.na(date_infection.RECIPIENT), all(.is.mrc(RECIPIENT))],
        chain2[is.na(date_infection.SOURCE), all(.is.mrc(SOURCE))]
    )
    cat("- CHECK: All individuals in infection pairs without ddates are in the MRC\n")
    cat("\t(removing them...)\n")
    chain2 <- chain2[!(.is.mrc(SOURCE) | .is.mrc(RECIPIENT))]

    if (0) {
        plot.phylopair.dates.scores(chain2, title = "Testing history")
        plot.phylopair.dates.scores(chain2, only.coherent = TRUE)
        plot.phylopair.dates.scores(chain2, only.contradict = TRUE)
        plot.phylopair.dates.scores(chain2, only.crossing = TRUE)
        plot.phylopair.dates.scores(chain2, only.contradict = TRUE, only.rect.pairs = FALSE)
    }

    tmp <- chain2[, date_last_negative.SOURCE >= date_first_positive.RECIPIENT]
    cat("- There are ", sum(tmp == TRUE, na.rm = TRUE), " couples with inconsistent first positive and last negative dates...\n")

    # Use age to give a min bound
    # setting date of last negative to be the minimum of 15th birthday and age first pos
    chain2 <- .double.merge(DT1 = chain2, DT2 = dage[, -c("AID", "age_first_positive")])
    cols <- grep("date_birth", names(chain2), value = TRUE)
    cols1 <- gsub("birth", "15yr", cols)
    chain2[, (cols) := lapply(.SD, function(x) x + as.integer(365.25 * 15)), .SDcols = cols]
    setnames(chain2, cols, cols1)
    cols <- grep("date_last", names(chain2), value = TRUE)
    chain2[, date_last_negative.SOURCE := pmax(date_last_negative.SOURCE, date_15yr.SOURCE, na.rm = T), by = "SOURCE"]
    chain2[, date_last_negative.RECIPIENT := pmax(date_last_negative.RECIPIENT, date_15yr.RECIPIENT, na.rm = T), by = "RECIPIENT"]

    # check
    # cols <- grep('date_first|date_last', names(chain2), value=TRUE)
    # chain2[, lapply(.SD,is.na) , .SDcols=cols][, lapply(.SD, sum) , .SDcols=cols]
    cols_src <- grep("SOURCE", names(chain2), value = TRUE)
    cols_rcp <- grep("RECIPIENT", names(chain2), value = TRUE)
    stopifnot(chain2[, all(date_first_positive.RECIPIENT >= date_last_negative.RECIPIENT, na.rm = T)])
    stopifnot(chain2[, all(date_first_positive.SOURCE >= date_last_negative.SOURCE, na.rm = T)])

    if (0) # date of infection does may not fall in range if we use "1 yr prior" DoI.
        {
            chain2[date_infection.RECIPIENT < date_last_negative.RECIPIENT, ..cols_rcp]
            chain2[date_infection.SOURCE < date_last_negative.SOURCE, ..cols_src]
            plot.phylopair.dates.scores(chain2[date_infection.SOURCE < date_last_negative.SOURCE])
            plot.phylopair.dates.scores(chain2, only.coherent = TRUE)
            plot.phylopair.dates.scores(chain2, only.contradict = TRUE)
            plot.phylopair.dates.scores(chain2, only.crossing = TRUE)
            plot.phylopair.dates.scores(chain2, only.contradict = TRUE, only.rect.pairs = FALSE)
        }

    cat("- Change direction of tranmsission for couple with contradicting times of infection...\n")

    cols <- grep("SOURCE|RECIPIENT", names(chain2), value = TRUE)
    .swap.src.rec <- function(x) {
        y <- copy(x)
        src <- grep("SOURCE", x)
        rcp <- grep("RECIPIENT", x)
        y[src] <- gsub("SOURCE", "RECIPIENT", x[src])
        y[rcp] <- gsub("RECIPIENT", "SOURCE", x[rcp])
        y
    }
    cols1 <- .swap.src.rec(cols)
    # chain2[, (uniqueN(RECIPIENT) == .N) ]
    idx <- chain2[date_last_negative.SOURCE >= date_first_positive.RECIPIENT, SOURCE]
    idx0 <- idx[!idx %in% chain2$RECIPIENT]
    idx1 <- idx[idx %in% chain2$RECIPIENT]
    cat("\tonly if the SOURCE does not appear as a RECIPIENT elsewhere(", length(idx1), "\n")
    chain2[date_last_negative.SOURCE >= date_first_positive.RECIPIENT & SOURCE %in% idx0,
        (cols1) := (.SD),
        .SDcols = cols
    ]
    cat("\t(exclude those others)")
    chain2 <- chain2[date_last_negative.SOURCE < date_first_positive.RECIPIENT]

    # individuals involved in multiple events:
    # SOURCEs should be infected

    cat("- Attribute date of infections based on transmission network structure...\n")

    shrink.ranges.doi.coherently.to.network <- function(CHAIN) {
        cat("Shrinking date of infection dates based on network relationships...\n")
        get_range_doi_by_studyid <- function(DT) {
            cols <- c("ID", "MIN", "MAX")
            cols0 <- c("SOURCE", "date_last_negative.SOURCE", "date_first_positive.SOURCE")
            cols1 <- c("RECIPIENT", "date_last_negative.RECIPIENT", "date_first_positive.RECIPIENT")
            tmp0 <- DT[, unique(.SD), .SDcols = cols0]
            tmp1 <- DT[, unique(.SD), .SDcols = cols1]
            setnames(tmp0, cols0, cols)
            setnames(tmp1, cols1, cols)
            tmp0[, uniqueN(ID) == .N]
            tmp1[, uniqueN(ID) == .N]
            tmp <- unique(rbind(tmp0, tmp1))
            stopifnot(tmp[, uniqueN(ID) == .N])
            return(tmp)
        }

        get_doi <- function(s, col) {
            stopifnot(length(col) == 1)
            dpairs_doi[ID %in% s, ..col][[1]]
        }

        shrink_max_doi_source <- function(PAIR, DATE) {
            # for each source, the MAX doi must preceed the MAX doi of the recipients
            tmp <- PAIR[, list(MAX2 = min(get_doi(RECIPIENT, "MAX"))), by = SOURCE]
            DATE <- merge(DATE, tmp, all.x = T, by.y = "SOURCE", by.x = "ID")

            if (DATE[MAX2 < MAX, .N == 0]) {
                DATE[, MAX2 := NULL]
                return(list(DATE = DATE, done = TRUE))
            }

            DATE[MAX2 < MAX, MAX := MAX2]
            DATE[, MAX2 := NULL]
            list(DATE = DATE, done = FALSE)
        }

        shrink_min_doi_recipient <- function(PAIR, DATE) {
            # for each recipient, the MIN doi must succeed the MIN doi of the source
            tmp <- PAIR[, list(MIN2 = get_doi(SOURCE, "MIN")), by = RECIPIENT]
            DATE <- merge(DATE, tmp, all.x = T, by.y = "RECIPIENT", by.x = "ID")

            if (DATE[MIN2 > MIN, .N == 0]) {
                DATE[, MIN2 := NULL]
                return(list(DATE = DATE, done = TRUE))
            }

            DATE[MIN2 > MIN, MIN := MIN2]
            DATE[, MIN2 := NULL]
            list(DATE = DATE, done = FALSE)
        }

        dpairs <- chain2[, .(SOURCE, RECIPIENT)]
        dpairs_doi <- get_range_doi_by_studyid(CHAIN)

        if (dpairs_doi[, any(is.na(MAX))]) {
            warning("Setting unknown date first positives as max date. Need changing\n")
            tmp <- dpairs_doi[, max(MAX, na.rm = T)]
            dpairs_doi[is.na(MAX), MAX := tmp]
        }

        done <- FALSE
        count <- 0
        while (done == FALSE) {
            count <- count + 1
            cat("\tIteration number:", count, "...\n")
            tmp <- shrink_max_doi_source(dpairs, dpairs_doi)
            dpairs_doi <- tmp$DATE
            done <- tmp$done
            tmp <- shrink_min_doi_recipient(dpairs, dpairs_doi)
            dpairs_doi <- tmp$DATE
            done <- done & tmp$done
        }

        return(list(P = dpairs, D = dpairs_doi))
    }

    tmp <- shrink.ranges.doi.coherently.to.network(chain2)
    dpairs <- copy(tmp$P)
    dpairs_doi <- copy(tmp$D)
    rm(tmp)
    setnames(dpairs_doi, "ID", "study_id")

    if (0) {
        tmp <- .double.merge(DT1 = chain2, dpairs_doi)
        idx <- tmp[, .N, by = "SOURCE"][N == 2, SOURCE]
        plot.phylopair.dates.scores(tmp[SOURCE %in% idx[3]], daterange.vars = c("MIN", "MAX"), add.dots = FALSE)
        plot.phylopair.dates.scores(tmp[SOURCE %in% idx[3]], add.dots = FALSE)
        plot.phylopair.dates.scores(tmp, daterange.vars = c("MIN", "MAX"), add.dots = FALSE)
        plot.phylopair.dates.scores(tmp, daterange.vars = c("MIN", "MAX"), add.dots = FALSE)
        plot.phylopair.dates.scores(tmp, daterange.vars = c("MIN", "MAX"), only.contradict = T, add.dots = FALSE)
        plot.phylopair.dates.scores(tmp, daterange.vars = c("MIN", "MAX"), only.coherent = T, add.dots = FALSE)
    }

    # now need to think about how to attribute the date of infection...
    # for simple case: midpoint attribution
    # for sources with multiple recipients: generalisation averaging pfds
    # for chains of A -> B -> C ... attribute the last DoI, and then "backpropagate" with midpoint assignment

    assign_dates_of_infection <- function(dpairs, dpairs_doi) {
        # attribute "depth" of transmission chain
        depth <- 1
        dpairs[, DEPTH := depth]
        tmp <- dpairs[DEPTH == depth, unique(SOURCE)]
        while (length(tmp)) {
            depth <- depth + 1
            dpairs[RECIPIENT %in% tmp, DEPTH := depth]
            tmp <- dpairs[DEPTH == depth, unique(SOURCE)]
        }
        dpairs[, table(DEPTH)]

        traceback_depth <- function(src, rcp, arrows = F, df = T) {
            tmp <- dpairs[SOURCE == src & RECIPIENT == rcp]
            if (tmp$DEPTH == 1) {
                if (df) {
                    return(data.frame(SOURCE = scr, RECIPIENT = rcp, DEPTH = tmp$DEPTH))
                }
                out <- c(src, rcp)
                if (arrows) {
                    out <- paste0(out, collapse = " -> ")
                }
                return(out)
            }

            tmp <- dpairs[SOURCE == rcp & DEPTH == tmp$DEPTH - 1, .(SOURCE, RECIPIENT)]

            if (df) {
                return(rbind(
                    data.frame(SOURCE = src, RECIPIENT = rcp, DEPTH = tmp$DEPTH),
                    traceback_depth(tmp$SOURCE, tmp$RECIPIENT, df = T)
                ))
            }

            out <- c(rcp, traceback_depth(tmp$SOURCE, tmp$RECIPIENT))
            if (arrows) {
                out <- paste0(out, collapse = " -> ")
            }
            return(out)
        }

        # get doi source | recipients bounds, and then assign midpoint for recps.
        assign.date.infection.sources <- function(src, d, precision = 50) {
            # src <- dpairs[,.N, by='SOURCE'][N==2, SOURCE[1]]
            # cat(src, '\n')
            tmp <- dpairs[SOURCE == src & DEPTH == d, -"DEPTH"]
            tmp <- .double.merge(tmp, dpairs_doi)

            # if any source has a known date of infection,
            # assign midpoint of plausible interval
            if (tmp[!is.na(DOI.RECIPIENT), .N]) {
                x_max <- tmp[!is.na(DOI.RECIPIENT), min(DOI.RECIPIENT)]
                return(tmp[, mean.Date(c(MIN.SOURCE[1], x_max))])
            }

            # there may be some analytic formula, but lets just use numeric integration
            x_range <- tmp[, seq(from = MIN.SOURCE[1], to = MAX.SOURCE[1], length.out = precision)]
            rcp_cross_section_length <- function(x) {
                tmp1 <- tmp[, .(L = as.numeric(MAX.RECIPIENT - max(MIN.RECIPIENT, x)) / 365.25), by = "RECIPIENT"]
                tmp1[, prod(L)]
            }
            x_pdf <- sapply(x_range, rcp_cross_section_length)
            x_cdf <- cumsum(x_pdf) / sum(x_pdf)
            x_range[which.min(abs(x_cdf - .5))]
        }

        assign.date.infection.recipient <- function(rcp, d) {
            # rcp <- idx[1]
            # cat(rcp, '\n')
            src <- dpairs[RECIPIENT == rcp, SOURCE]
            dsrc <- dpairs_doi[study_id == src]
            drcp <- dpairs_doi[study_id == rcp]
            x_min <- dsrc$DOI
            x_min <- max(x_min, drcp$MIN)
            x_mid <- mean.Date(c(x_min, drcp$MAX))
            x_mid
        }


        dpairs_doi[, DOI := NA_Date_]
        depth <- dpairs[, max(DEPTH)] + 1

        while (depth > 0 & dpairs_doi[, any(is.na(DOI))]) {
            depth <- depth - 1
            cat("Assigning dates of infection for pairs at depth:", depth, "...\n")

            dpairs_doi[study_id %in% dpairs[DEPTH == depth, SOURCE] & is.na(DOI),
                DOI := assign.date.infection.sources(study_id, depth),
                by = "study_id"
            ]
            cat(dpairs_doi[!is.na(DOI), .N], "\n")

            dpairs_doi[study_id %in% dpairs[DEPTH == depth, RECIPIENT] & is.na(DOI),
                DOI := assign.date.infection.recipient(study_id, depth),
                by = "study_id"
            ]
            cat(dpairs_doi[!is.na(DOI), .N], "\n")

            .check <- function(src, rcp) {
                tmp_s0 <- dpairs_doi[study_id == src, DOI >= MIN & DOI <= MAX]
                tmp_r0 <- dpairs_doi[study_id == rcp, DOI >= MIN & DOI <= MAX]
                tmp_s1 <- dpairs_doi[study_id == src, DOI]
                tmp_r1 <- dpairs_doi[study_id == rcp, DOI]
                tmp_s0 & tmp_r0 & (tmp_s1 <= tmp_r1)
            }
            tmp <- dpairs[, .check(SOURCE, RECIPIENT), by = c("SOURCE", "RECIPIENT")]
            stopifnot(tmp[, all(V1 == T, na.rm = T)])
        }
        dpairs_doi
    }

    # prepare output

    dpairs_doi <- assign_dates_of_infection(dpairs, dpairs_doi)

    dage <- merge(dage[, .(study_id, date_birth)],
        dpairs_doi[, .(study_id, DOI)],
        by = "study_id"
    )
    dage[, age_at_infection := .year.diff(DOI, date_birth)]

    list(pairs = .double.merge(DT1 = chain2, dpairs_doi), age = dage)
}



get.transmission.cluster.ids <- function(DT, check = F) {
    tmp <- subset(DT, select = c("SOURCE", "RECIPIENT"))
    tmp[, GROUP := 1:.N, ]
    setkey(tmp, GROUP, SOURCE)

    ids.all <- tmp[, unique(c(SOURCE, RECIPIENT))]
    for (idx in ids.all)
    {
        grps <- tmp[SOURCE == idx | RECIPIENT == idx, GROUP]
        if (length(grps) > 1) {
            tmp[GROUP %in% grps, GROUP := min(grps)]
        }
    }

    # rename
    setkey(tmp, GROUP)
    tmp1 <- data.table(GROUP = unique(tmp$GROUP))
    tmp1[, NEWGROUP := 1:.N]
    tmp <- merge(tmp, tmp1, by = "GROUP")
    tmp <- tmp[, .(SOURCE, RECIPIENT, GROUP = NEWGROUP)]

    if (check) {
        # does everyone appear in 1 group only?
        tmp1 <- rbind(tmp[, .(id = SOURCE, GROUP)], tmp[, .(id = RECIPIENT, GROUP)]) |> unique()
        .f <- function(x) {
            tmp1[id == x, uniqueN(GROUP)]
        }
        stopifnot(all(sapply(ids.all, .f) == 1))
    }

    tmp
}

plot.polygon <- function(df, centroid = FALSE, todate = FALSE) {
    tmp <- as.data.table(df)
    names(tmp) <- c("X", "Y")

    if (centroid) {
        tmp1 <- get.polygon.centroid(df)
    }

    if (todate) {
        tmp <- tmp[, lapply(.SD, .int2date)]
        if (centroid) {
            tmp1 <- .int2date(tmp1)
        }
    }

    g <- ggplot(tmp, aes(x = X, y = Y)) +
        geom_polygon(alpha = .5, fill = "blue") +
        geom_abline(slope = 1, color = "red", linetype = "dashed") +
        theme_bw() +
        theme(legend.position = "bottom")

    if (centroid) {
        g <- g + geom_point(aes(x = tmp1[1], y = tmp1[2]))
    }
    g
}

get.polygon.cohordinates <- function(DT, usenames = FALSE) {
    # DT <- copy(TMP)
    x.min <- DT[, MIN.SOURCE]
    x.max <- DT[, MAX.SOURCE]
    y.min <- DT[, MIN.RECIPIENT]
    y.max <- DT[, MAX.RECIPIENT]

    vertices1 <- rbind(c(x.min, y.min), c(x.max, y.min), c(x.max, y.max), c(x.min, y.max))

    # count how many vertices are above the y=x line
    tmp <- vertices1[, 2] > vertices1[, 1]
    N_above <- sum(tmp)
    N_above_x <- uniqueN(vertices1[, 1][tmp])


    # based on this, define the other
    if (N_above == 0) {
        return(c(NA, NA))
    }

    if (N_above == 4) {
        return(vertices1)
    }

    if (N_above == 1) {
        vertices2 <- rbind(c(x.min, y.max), c(x.min, x.min), c(y.max, y.max))
    }

    if (N_above == 3) {
        vertices2 <- rbind(c(x.min, y.max), c(x.min, y.min), c(y.min, y.min), c(x.max, x.max), c(x.max, y.max))
    }

    if (N_above == 2) {
        if (N_above_x == 2) {
            vertices2 <- rbind(c(x.min, y.max), c(x.max, y.max), c(x.max, x.max), c(x.min, x.min))
        }
        if (N_above_x == 1) {
            vertices2 <- rbind(c(x.min, y.min), c(x.min, y.max), c(y.max, y.max), c(y.min, y.min))
        }
    }

    if (usenames) {
        colnames(vertices2) <- c(DT[, SOURCE], DT[, RECIPIENT])
    }

    if (check.if.anticlockwise(as.data.table(vertices2)) == "clockwise") {
        vertices2 <- vertices2[nrow(vertices2):1, ]
    }

    vertices2
}


get.polygon.area <- function(COHORDS, year.normalise = TRUE) {
    if (is.na(COHORDS[[1]][1])) {
        return(0)
    }

    if (!all(first(COHORDS) != last(COHORDS))) {
        COHORDS <- rbind(COHORDS, first(COHORDS))
    }

    x_cohords <- COHORDS[, 1]
    y_cohords <- COHORDS[, 2]
    formula1 <- x_cohords * shift(y_cohords, n = -1)
    formula2 <- y_cohords * shift(x_cohords, n = -1)
    area <- abs(sum(formula1 - formula2, na.rm = TRUE)) / 2
    if (year.normalise) {
        area <- area / (365.25^2)
    }

    return(area)
}

get.polygon.centroid <- function(COHORDS, usenames = FALSE, as.dt = FALSE) {
    # COHORDS <- copy(DT)
    N <- nrow(COHORDS)
    # from: https://www.johndcook.com/blog/2020/04/05/center-of-mass-and-vectorization/
    if (!all(first(COHORDS) != last(COHORDS))) {
        COHORDS <- rbind(last(COHORDS), COHORDS)
    }

    x_cohords <- COHORDS[, 1] |>
        unlist() |>
        unname()
    y_cohords <- COHORDS[, 2] |>
        unlist() |>
        unname()

    x <- x_cohords[1:N]
    y <- y_cohords[1:N]
    x_dash <- x_cohords[2:(N + 1)]
    y_dash <- y_cohords[2:(N + 1)]

    const <- x * y_dash - x_dash * y
    A <- sum(const) / 2
    c_x <- 1 / (6 * A) * (x + x_dash) %*% (const)
    c_y <- 1 / (6 * A) * (y + y_dash) %*% (const)

    out <- c(c_x, c_y)


    if (as.dt) {
        out <- data.table(V1 = c_x, V2 = c_y)
    }

    if (usenames == TRUE) {
        names(out) <- names(COHORDS)
    }

    return(out)
}

# intersect.polygon.45degline(COHORDS, intercept=365.25, lower=TRUE, keepnames=FALSE) |> plot.polygon()
# intersect.polygon.45degline(COHORDS, intercept=365.25, keepnames=FALSE) |> plot.polygon()

intersect.polygon.45degline <- function(COHORDS, intercept = 0, lower = FALSE, keepnames = TRUE) {
    # intercept=365

    # Note this exploits the fact that the polygon's only has vertical, horizontal and 45deg lines
    # and it expects the line to be of the type y = x + intercept
    # Also, we assume the COHORDS are already sorted anti-clockwise
    tmp <- copy(COHORDS)
    # plot.polygon(COHORDS)
    names(tmp) <- c("X", "Y")

    # plot.polygon(tmp)

    # SETUP
    x_cohords <- tmp[, X]
    y_cohords <- tmp[, Y]
    x_ranges <- tmp[, .(MIN_X = min(X), MAX_X = max(X)), by = Y]
    y_ranges <- tmp[, .(MIN_Y = min(Y), MAX_Y = max(Y)), by = X]
    setkey(x_ranges, Y)
    setkey(y_ranges, X)

    # GET INTERSECTIONS
    # to get duplicated cohordinates
    .f <- function(x) unique(x[duplicated(x)])

    # get the intersections with the horizontal segments
    y_cohord_horsegs <- .f(y_cohords)
    h_intersect <- data.table(X = y_cohord_horsegs - intercept, Y = y_cohord_horsegs)
    h_intersect <- merge(h_intersect, x_ranges, by = "Y")
    h_intersect <- h_intersect[X <= MAX_X & X >= MIN_X, ]

    # get the intersections with the vertical segments
    x_cohord_versegs <- .f(x_cohords)
    v_intersect <- data.table(X = x_cohord_versegs, Y = x_cohord_versegs + intercept)
    v_intersect <- merge(v_intersect, y_ranges, by = "X")
    v_intersect <- v_intersect[Y <= MAX_Y & Y >= MIN_Y, ]

    intersections <- rbind(h_intersect[, .(X, Y)], v_intersect[, .(X, Y)])

    # NOW SELECT VERTICES: incorrect
    if (nrow(intersections) == 0) {
        checkwhetherlower <- unique(tmp[, Y - X <= intercept])
        stopifnot(length(checkwhetherlower) == 1)

        if (keepnames) {
            names(tmp) <- names(COHORDS)
        }

        if ((lower & checkwhetherlower) | (!lower & !checkwhetherlower)) {
            return(tmp)
        } else {
            # just want to return an empty dt
            return(tmp[1 > 2])
        }
    }


    if (lower) {
        setorder(intersections, -"X")
        idx <- tmp[, Y - X <= intercept]
    } else {
        setorder(intersections, "X")
        idx <- tmp[, Y - X >= intercept]
    }

    tmp[, `:=`(GROUP = rleid(idx), KEEP = idx)]
    subpos <- tmp[KEEP == FALSE, unique(GROUP)]
    stopifnot(length(subpos) == 1 | (subpos %in% range(tmp$GROUP) & length(subpos == 2)))
    intersections[, `:=`(GROUP = subpos[1], KEEP = TRUE)]
    tmp <- rbind(intersections, tmp)[KEEP == TRUE]
    setorder(tmp, GROUP)
    tmp[, `:=`(GROUP = NULL, KEEP = NULL)]

    intersections
    if (keepnames) {
        names(tmp) <- names(COHORDS)
    }
    tmp
}

check.if.anticlockwise <- function(COHORDS) {
    # https://stackoverflow.com/questions/1165647/how-to-determine-if-a-list-of-polygon-points-are-in-clockwise-order#:~:text=Here's%20a%20simple%20one%20that,the%20curve%20is%20counter%2Dclockwise.
    tmp <- copy(COHORDS)
    names(tmp) <- c("X", "Y")

    tmp <- rbind(last(tmp), tmp)
    tmp[, X_shift := shift(X, 1)]
    tmp[, Y_shift := shift(Y, 1)]
    tmp <- tmp[, sum((X_shift - X) * (Y_shift + Y), na.rm = TRUE)]
    fifelse(sign(tmp) == 1, yes = "anticlockwise", no = "clockwise")
}

load.rounds <- function() {
    rounds_env <- new.env()
    load(path.round.timeline, envir = rounds_env)
    rbind(
        rounds_env$df_round_fishing,
        rounds_env$df_round_inland
    ) -> dround
    tmp <- dround[, list(
        COMM = "aggregated",
        min_sample_date = min(min_sample_date),
        max_sample_date = min(max_sample_date)
    ), by = "round"]
    dround <- rbind(dround, tmp)
    setkey(dround, COMM, min_sample_date)
    dround
}

get.volume.genints.multid <- function(DT, N_IDS = NULL,
                                      n_samples = 10e5,
                                      provide_samples = NULL,
                                      range_gi = NULL,
                                      verbose = FALSE,
                                      importance_weights = NULL,
                                      get_volume = FALSE,
                                      df_round = data.frame()) {
    # DT <- dcohords[ GROUP == 1, IDS[[1]]]; df_round = df_round_gi; N_IDS=3; importance_weights=dinfectiousness


    if (nrow(df_round)) {
        round_dates <- df_round[, .(ROUND, MIN_SAMPLE_DATE, MAX_SAMPLE_DATE)]
        cols <- names(round_dates)[names(round_dates) %like% "DATE"]
        round_dates[, (cols) := lapply(.SD, as.integer), .SDcols = cols]
    }

    if (!is.null(provide_samples)) {
        n_samples <- nrow(provide_samples)
    }

    # Checking
    .s <- function(str) stop(paste0("Error in `get.volume.genints.multid`: ", str))
    if (is.null(N_IDS)) {
        .s("Provide N_IDS")
    }
    if (ncol(DT) != 2) {
        .s("DT should only have 2 columns")
    }
    if (!is.null(importance_weights) &
        !all(c("START", "END", "VALUE") %in% names(importance_weights))) {
        .s("importance_weights should have columns: START, END, VALUE")
    }

    # rank by how many times they appear in source?
    tmp <- copy(DT)
    codes <- LETTERS[(27 - N_IDS):26]
    ids <- unique(unlist(tmp))
    tbl <- table(factor(tmp[[1]], levels = ids))
    ids <- names(sort(tbl, decreasing = TRUE))

    codes <- data.table(ID = ids, CODE = codes)
    setkey(codes, ID)
    setkey(tmp, SOURCE)
    tmp <- tmp[codes[, .(ID = ID, C1 = CODE)], on = .(SOURCE = ID)]
    setkey(tmp, RECIPIENT)
    tmp <- tmp[codes[, .(ID = ID, C2 = CODE)], on = .(RECIPIENT = ID)]
    dineq <- tmp[!is.na(SOURCE) & !is.na(RECIPIENT), .(INEQ = paste(C1, C2, sep = "<"), C1, C2)]

    # To do MCMC integration, we need to recover the ranges...
    setkey(dpairs, SOURCE, RECIPIENT)
    tmp <- dpairs[DT]
    tmp <- rbind(
        tmp[, .(ID = SOURCE, MIN = MIN.SOURCE, MAX = MAX.SOURCE)],
        tmp[, .(ID = RECIPIENT, MIN = MIN.RECIPIENT, MAX = MAX.RECIPIENT)]
    ) |> unique()
    dranges <- codes[tmp, on = "ID"]
    cols <- c("MIN", "MAX")
    new_cols <- paste0("NEW_", cols)
    dranges[, (new_cols) := lapply(.SD, as.numeric), .SDcols = cols]

    # MCMC sampling
    if (!is.null(provide_samples)) {
        .f <- function(d) {
            spls <- provide_samples[[d]]
            dranges[d, NEW_MIN + (NEW_MAX - NEW_MIN) * spls]
        }
        dsamples <- lapply(seq(nrow(dranges)), .f)
        names(dsamples) <- dranges$CODE
        setDT(dsamples)
        dsamples[, I := 1:.N, ]
        setcolorder(dsamples, "I")
    } else {
        dsamples <- dranges[, .(I = 1:n_samples, V = runif(n_samples, min = NEW_MIN, max = NEW_MAX)), by = CODE]
        dsamples <- dcast(dsamples, I ~ CODE, value.var = "V")
    }

    ### now we can ask questions of interest.
    # use importance weights here...

    # Formulate region of interest with inequalities
    impsamples <- dineq[, subset(dsamples, select = c("I", C1, C2)), by = INEQ]
    oldcols <- setdiff(names(impsamples), c("I", "INEQ"))
    newcols <- c("S", "R")
    setnames(impsamples, oldcols, newcols)
    impsamples[, IN := S < R]
    impsamples[, IN := all(IN), by = I]
    impsamples[IN == TRUE, GI := (R - S) / 365.25, by = "I"]

    # attribute importance weights based on GI
    for (rw in seq_len(nrow(importance_weights)))
    {
        rng <- importance_weights[rw][, c(START, END, VALUE)]
        impsamples[GI > rng[1] & GI <= rng[2], W1d := rng[3]]
    }

    # attribute round of infection of recipient
    if (nrow(df_round)) {
        setkey(impsamples, R)
        for (rnd in round_dates$ROUND)
        {
            # rnd <- round_dates$ROUND[1]
            rng <- round_dates[ROUND == rnd, c(MIN_SAMPLE_DATE, MAX_SAMPLE_DATE)]
            impsamples[R > rng[1] & R < rng[2], ROUND := rnd]
        }
    }

    if (get_volume) {
        ### Volume
        if (verbose) {
            cat("Computing volumes...\n")
        }

        volume <- dranges[, prod(NEW_MAX - NEW_MIN) / 365.25^N_IDS]
        prop <- impsamples[, all(IN), by = I][, mean(V1)]
        volume <- volume * prop

        if (verbose) {
            cat("Acceptance rate:", round(100 * prop, 2), "%\n")
        }
    }
    impsamples <- impsamples[IN == TRUE]

    # get quantiles
    if (verbose) {
        cat("Computing centroids...\n")
    }

    if (is.null(importance_weights)) {
        centroid <- dsamples[impsamples[, I], lapply(.SD, median), .SDcols = setdiff(names(dsamples), "I")]
    } else {
        tmp0 <- impsamples[, .(WEIGHTS = prod(W1d)), by = I][WEIGHTS != 0]
        dsamples <- merge(dsamples, tmp0, on = "I")

        dsamples

        .get.ranges <- function(nm, p = c(.10, .25, .5, .75, .9)) {
            cols <- c(nm, "WEIGHTS")
            tmp <- subset(dsamples, select = cols)
            tmp <- setorderv(tmp, nm)
            tmp[, CS := cumsum(WEIGHTS)]
            tmp[, Q := CS / last(CS)]

            .get.quantile <- function(x) {
                tmp[, .int2date(unlist(.SD[which.min(abs(Q - x))])), .SDcols = nm]
            }

            out <- lapply(p, .get.quantile)
            names(out) <- c("CL", "IL", "M", "IU", "CU")
            out
        }

        cols <- setdiff(names(dsamples), c("I", "WEIGHTS"))
        centroid <- lapply(cols, .get.ranges)
        centroid <- rbindlist(centroid)
        centroid[, CODE := cols]
        centroid <- merge(codes, centroid, by = "CODE")
        centroid[, CODE := NULL]

        impsamples <- merge(impsamples, dsamples[, .(I, WEIGHTS)], by = "I")
        impsamples[, TOT_WEIGHTS := sum(WEIGHTS), by = INEQ]
    }

    nms <- codes[, ID[chmatch(names(centroid), CODE)]]
    names(centroid) <- nms

    # generation intervals
    if (verbose) {
        cat("Computing GI probs...\n")
    }

    if (!is.null(importance_weights) & !is.null(range_gi)) {
        .get.gi.probs <- function(gi) {
            idx <- unique(impsamples[, INEQ])
            out_tmp <- data.table(INEQ = idx, V2 = 0)
            out <- impsamples[R < S + gi * 365.25, sum(WEIGHTS) / unique(TOT_WEIGHTS), by = INEQ]
            out <- merge(out_tmp, out, all.x = TRUE)
            out <- out[, .(V1 = fcoalesce(V1, V2)), by = "INEQ"]
            setnames(out, "V1", paste0("GI", gi))
            out
        }
        tmp1 <- lapply(range_gi, .get.gi.probs)
        tmp1 <- Reduce("merge", tmp1)

        tmp1 <- melt(tmp1, "INEQ", variable.name = "GI", value.name = "")
        tmp1[, GI := as.numeric(gsub("^GI", "", GI))]
        setkey(codes, CODE)
        dineq <- dineq[, list(SOURCE = codes[C1, ID], RECIPIENT = codes[C2, ID], INEQ)]
        genints <- merge(dineq, tmp1)[, -"INEQ"]
    } else if (!is.null(range_gi)) {
        .f <- function(gi) {
            out <- impsamples[, .(R < S + gi * 365.25)]
            names(out) <- paste0("GI", gi)
            out
        }
        tmp1 <- lapply(range_gi, .f)
        tmp1 <- Reduce("cbind", tmp1)
        impsamples <- cbind(impsamples, tmp1)

        cols <- names(impsamples)[names(impsamples) %like% "^G"]
        impsamples <- impsamples[, lapply(.SD, mean), .SDcols = cols, by = INEQ]
        impsamples <- melt(impsamples, "INEQ", variable.name = "GI", value.name = "PR")
        impsamples[, GI := as.numeric(gsub("GI", "", GI))]
        setkey(impsamples, INEQ)
        impsamples[, `:=`(
            SOURCE = gsub("^(.*?)<.*?$", "\\1", INEQ),
            RECIPIENT = gsub(".*?<(.*?)$", "\\1", INEQ)
        )]
        impsamples <- merge(impsamples, codes[, .(CODE, ID_S = ID)], by.x = "SOURCE", by.y = "CODE")
        impsamples <- merge(impsamples, codes[, .(CODE, ID_R = ID)], by.x = "RECIPIENT", by.y = "CODE")
        cols <- c("SOURCE", "RECIPIENT")
        impsamples[, (cols) := NULL]
        old_cols <- paste0("ID_", c("S", "R"))
        setnames(impsamples, old_cols, cols)
        impsamples[, INEQ := NULL]
        genints <- impsamples
    }

    if (verbose) {
        cat("Computing probability of recipient in different round...\n")
    }

    if (nrow(df_round)) {
        round_probs <- impsamples[, .(W = sum(WEIGHTS)), by = c("INEQ", "ROUND")]
        round_probs[, TOT_W := sum(W), by = "INEQ"]
        round_probs[, P := W / TOT_W]
        round_probs <- merge(round_probs, dineq, by = "INEQ")
        cols <- c("SOURCE", "RECIPIENT", "ROUND", "P")
        round_probs <- subset(round_probs, select = cols)
        setkeyv(round_probs, cols)
    }

    if (verbose) {
        cat("Finalizing...\n")
    }

    out <- list(CENTROID = centroid)

    if (get_volume) {
        out <- append(out, list(VOLUME = volume))
    }

    if (!is.null(range_gi)) {
        out <- append(out, list(GI = genints))
    }

    if (nrow(df_round)) {
        out <- append(out, list(RPROBS = round_probs))
    }

    out
}

if (0) # 2D, uniform pdf
    {
        # get the cohordinates of each 2D plausible region
        setkey(dpairs, SOURCE, RECIPIENT)
        dcohords[N_IDS == 2,
            {
                pair.data <- dpairs[IDS]
                # cat(pair.data[, SOURCE], '-', pair.data[, RECIPIENT], '\n')
                out <- data.table(get.polygon.cohordinates(pair.data, usenames = TRUE))
                list(COHORDS = list(out))
            },
            by = GROUP
        ] -> tmp
        dcohords <- merge(dcohords, tmp, all.x = T)

        # get the volume of every 2D plausible region
        dcohords[N_IDS == 2, VOLUME := get.polygon.area(COHORDS[[1]]), by = GROUP]

        # get the centroid of every 2D plausible region
        tmp <- dcohords[N_IDS == 2,
            {
                get.polygon.centroid(COHORDS[[1]], as.dt = TRUE, usenames = TRUE)
            },
            by = GROUP
        ]
        tmp <- group_nest_dt(tmp, GROUP, .key = "CENTROIDS")
        dcohords <- merge(dcohords, tmp, all.x = TRUE)

        # compute probability of generation interval less than certain threshold.
        cat("compute generation intervals\n")
        range_gi <- c(.5, 1:10)
        dcohords[N_IDS == 2,
            {
                cat(GROUP, "\t")
                .f <- function(i) {
                    z0 <- intersect.polygon.45degline(COHORDS[[1]], intercept = i * 365.25, lower = TRUE)
                    z1 <- get.polygon.area(z0)
                }
                z1 <- sapply(range_gi, .f)
                out <- list(GI = range_gi, PR = z1 / VOLUME)
            },
            by = GROUP
        ] -> tmp
        cat("\n")
        tmp <- group_nest_dt(tmp, GROUP, .key = "GENINTS")
        dcohords <- merge(dcohords, tmp, all.x = TRUE)
    }

prepare.pairs.input.for.bayesian.model <- function(DT) {
    # get predictions for time of infection
    merge(
        chain[, .(SOURCE, RECIPIENT)],
        DT,
        by.x = "RECIPIENT", by.y = "ID"
    ) -> dresults
    dresults

    # get sex and age at infection
    tmp <- meta[, .(AID = aid, SEX = sex, DB = date_birth)]

    merge(
        dresults,
        tmp,
        by.y = "AID", by.x = "SOURCE"
    ) -> dresults
    dresults[, AGE_INFECTION := as.numeric(round((M - DB) / 365.25, 1))]
    dresults[, DB := NULL]
    setnames(dresults, c("AGE_INFECTION", "SEX"), c("AGE_TRANSMISSION.SOURCE", "SEX.SOURCE"))

    merge(
        dresults,
        tmp,
        by.y = "AID", by.x = "RECIPIENT"
    ) -> dresults
    dresults[, AGE_INFECTION := as.numeric(round((M - DB) / 365.25, 1))]
    dresults[, DB := NULL]
    setnames(dresults, c("AGE_INFECTION", "SEX"), c("AGE_INFECTION.RECIPIENT", "SEX.RECIPIENT"))

    # get round
    for (r in df_round_gi[, ROUND]) {
        dresults[M %between% df_round_gi[ROUND == r, c(MIN_SAMPLE_DATE, MAX_SAMPLE_DATE)], ROUND.M := r]
    }

    dresults[, GROUP := NULL]
    cols <- c("SOURCE", "RECIPIENT", "SEX.SOURCE", "SEX.RECIPIENT", "CL", "IL", "M", "IU", "CU", "AGE_TRANSMISSION.SOURCE", "AGE_INFECTION.RECIPIENT", "ROUND.M")
    cols <- intersect(cols, names(dresults))
    setcolorder(dresults, cols)

    cols <- c("SOURCE", "RECIPIENT")
    chain_tmp <- chain[SCORE > threshold.likely.connected.pairs]
    # chain_tmp <- keep.likely.transmission.pairs(as.data.table(dchain), threshold.likely.connected.pairs)
    idx <- merge(dresults[, ..cols], chain_tmp[, ..cols], by = cols)
    idx[, DIRECTION := "phyloscanner"]
    dresults <- merge(dresults, idx, all.x = TRUE, by = cols)
    dresults[, DIRECTION := fcoalesce(DIRECTION, "serohistory")]

    dresults
}

get.community.type.at.infection.date <- function(DT, comm_number = TRUE) {
    # DT <- copy(dresults); comm_number=TRUE
    meta_env <- new.env()
    load(path.meta, envir = meta_env)

    cols <- c("aid", "comm", "round", "sample_date")
    dcomms <- subset(meta_env$meta_data, select = cols)
    names(dcomms) <- toupper(names(dcomms))
    dcomms[ROUND == "neuro", COMM := "neuro"]
    dcomms <- dcomms[!is.na(AID)]

    if (comm_number) 
    {
        .fp <- function(x) file.path(indir.deepsequencedata, "RCCS_R15_R18", x)
        path.quest <- .fp("quest_R15_R18_VoIs_220129.csv")
        path.metadata <- .fp("Rakai_Pangea2_RCCS_Metadata__12Nov2019.csv")
        paths <- c(path.quest, path.metadata)
        cols <- c("study_id", "round", "comm_num", "community_number")
        dcomms2 <- lapply(paths, fread, select = cols)
        dcomms2 <- lapply(dcomms2, function(DT) setnames(DT, names(DT), toupper(cols[1:3])))
        dcomms2 <- rbindlist(dcomms2)
        dcomms2 <- dcomms2[!is.na(ROUND)]
        dcomms2[, ROUND := fifelse(ROUND %like% "R0", yes = ROUND, no = paste0("R0", ROUND))]
        dcomms2[ROUND %like% "15.1", ROUND := "R015S"]
        dcomms2 <- unique(dcomms2)
        dcomms2[, uniqueN(COMM_NUM), by = c("STUDY_ID", "ROUND")][, stopifnot(all(V1 == TRUE))]
        dcomms2[, STUDY_ID := paste0("RK-", STUDY_ID)]
        aids <- aik[, AID[match(dcomms2$STUDY_ID, PT_ID)]]
        dcomms2[, AID := aids]
        dcomms2 <- dcomms2[!is.na(AID), -"STUDY_ID"]
        dcomms <- merge(dcomms, dcomms2, by = c("AID", "ROUND"), all.x = TRUE)
        dcomms[is.na(COMM_NUM), stopifnot(all(ROUND == "neuro"))]
    }

    .f <- function(x) paste0(unique(sort(as.character(x))), collapse = "_")

    cols <- names(dcomms)[names(dcomms) %like% 'COMM']
    dcomms_collapse <- dcomms[,lapply(.SD, .f),.SDcols=cols, by="AID"]
    # dcomms_collapse <- dcomms[, .(COMM = .f(COMM), COMM_NUM = .f(COMM_NUM)), by = "AID"]

    # For participants with date of infection, set
    # community as the comm with visit date closest to estimated infection time
    .f <- function(x) which.min(abs(fcoalesce(x, Inf)))
    DT[!is.na(M),
        {
            idx <- c(SOURCE, RECIPIENT)
            tmp <- dcomms[AID %in% idx, ]
            out0 <- tmp[, .(C = COMM[.f(as.numeric(SAMPLE_DATE - M))]), by = "AID"]
            CS <- out0[AID == SOURCE, C]
            CR <- out0[AID == RECIPIENT, C]
            out_list <- list(SOURCE = SOURCE, COMM.SOURCE = CS, COMM.RECIPIENT = CR)
            if (comm_number) {
                out1 <- tmp[, .(CN = COMM_NUM[.f(as.numeric(SAMPLE_DATE - M))]), by = "AID"]
                CNS <- out1[AID == SOURCE, CN]
                CNR <- out1[AID == RECIPIENT, CN]
                out_list <- c(out_list, list(COMM_NUM.SOURCE = CNS, COMM_NUM.RECIPIENT = CNR))
            }
            out_list
        },
        by = "RECIPIENT"
    ] -> dt_with_inf_date

    # find inf date also for remaining participants without estimated inf. date
    dt_without_inf_date <- double.merge(DT[is.na(M), .(SOURCE, RECIPIENT)], dcomms_collapse)

    # check
    dt_with_inf_date[, .N, by = c("RECIPIENT", "SOURCE")][, stopifnot(all(N == 1))]
    dt_without_inf_date[, .N, by = c("RECIPIENT", "SOURCE")][, stopifnot(all(N == 1))]
    tmp <- rbind(dt_without_inf_date, dt_with_inf_date)

    DT <- merge(DT, tmp, by = c("SOURCE", "RECIPIENT"), all.x = TRUE)
    DT[SEX.SOURCE == SEX.RECIPIENT]
    return(DT)
}

count.number.intersecting.ranges.tsi.sero.wrt.chain <- function(CHAIN, DCOMPARE) {
    # CHAIN <- copy(chain);     DCOMPARE <- copy(compare_ranges)
    # individuals in the source-recipient pairs
    idx <- CHAIN[, unique(c(SOURCE, RECIPIENT))]

    # individuals with null TSI plausible region &
    # individuals where TSI and serohistory ranges do not intersect
    idx.tsi.null <- DCOMPARE[, idx[!idx %in% AID]]
    idx.no.marg.intersect <- DCOMPARE[RNG.INTERSECT == FALSE, AID]
    idx2 <- c(idx.tsi.null, idx.no.marg.intersect)

    .g <- function(x) CHAIN[SOURCE %in% x | RECIPIENT %in% x, ]
    .g(idx.no.marg.intersect)
    .f <- function(x) CHAIN[SOURCE %in% x | RECIPIENT %in% x, .N]
    length(idx.tsi.null)
    length(idx.no.marg.intersect)
    cat(
        "There are", length(idx2), "individuals with null marginal intersections (", .f(idx2), "):\n\t-",
        length(idx.tsi.null), "(", .f(idx.tsi.null), ")had null TSI plausible region\n\t-",
        length(idx.no.marg.intersect), "(", .f(idx.no.marg.intersect), ")additional inds had null marginal interesections.\n"
    )
    cat("\tAccounting for", CHAIN[SOURCE %in% idx2 | RECIPIENT %in% idx2, .N], "S-R pairs\n")

    return(NULL)
}

get.round.dates <- function(file = path.round.timeline) {
    tmp_env <- new.env()
    load(file, envir = tmp_env)
    df_round_inland <- tmp_env$df_round_inland
    df_round_fishing <- tmp_env$df_round_fishing

    start_first_period_inland <- df_round_inland[round == "R010", min_sample_date] # "2003-09-26"
    stop_first_period_inland <- df_round_inland[round == "R015", max_sample_date] # "2013-07-05"
    start_second_period_inland <- df_round_inland[round == "R016", min_sample_date] #  "2013-07-08"
    stop_second_period_inland <- df_round_inland[round == "R018", max_sample_date] #  "2018-05-22"

    stopifnot(start_first_period_inland < stop_first_period_inland)
    stopifnot(stop_first_period_inland < start_second_period_inland)
    stopifnot(start_second_period_inland < stop_second_period_inland)

    cols <- c("ROUND", "MIN_SAMPLE_DATE", "MAX_SAMPLE_DATE", "INDEX_TIME", "INDEX_ROUND", "round")
    df_round_i <- make.df.round(df_round_inland) |> subset(select = cols)
    df_round_f <- make.df.round(df_round_fishing) |> subset(select = cols)

    cols2 <- cols[!cols %like% "SAMPLE_DATE|INDEX"]
    df_round_gi <- merge(df_round_i, df_round_f, by = cols2, all.x = TRUE)
    df_round_gi[, `:=`(
        MIN_SAMPLE_DATE = pmin(MIN_SAMPLE_DATE.x, MIN_SAMPLE_DATE.y, na.rm = TRUE),
        MAX_SAMPLE_DATE = pmax(MAX_SAMPLE_DATE.x, MAX_SAMPLE_DATE.y, na.rm = TRUE),
        INDEX_TIME = INDEX_TIME.x,
        INDEX_ROUND = INDEX_ROUND.x
    )]
    df_round_gi <- subset(df_round_gi, select = !names(df_round_gi) %like% ".x$|.y$")
    df_round_gi <- subset(df_round_gi, select = cols)

    df_round_gi[, MAX_SAMPLE_DATE := pmax(MAX_SAMPLE_DATE, shift(MIN_SAMPLE_DATE, -1), na.rm = TRUE)]
    rm(tmp_env, df_round_inland, df_round_fishing, df_round_i, df_round_f)

    df_round_gi
}

plot.chains <- function(DPAIRS, size.threshold = 1, ttl = NA, sbttl = NA) {
    library(igraph)
    tmp <- copy(DPAIRS)

    if (size.threshold > 1) {
        tmp1 <- graph_from_data_frame(tmp, directed = TRUE, vertices = NULL) |>
            components()
        comps <- tmp1$membership
        comps <- data.table(SOURCE = names(comps), COMP = unname(comps))
        tmp <- merge(comps, tmp)
        tmp[, SIZE := .N, by = "COMP"]
        tmp <- tmp[SIZE >= size.threshold, .(SOURCE, RECIPIENT)]
    }

    .gs <- function(x) gsub("AID", "", x)
    graph_from_data_frame(tmp[, .(SOURCE = .gs(SOURCE), RECIPIENT = .gs(RECIPIENT))], directed = TRUE, vertices = NULL) |>
        plot(
            vertex.label = NA,
            vertex.label.cex = .7,
            vertex.label.distance = .5,
            vertex.size = 1,
            edge.arrow.size = .4,
            sub = sbttl,
            main = ttl
        ) -> p
    return(p)
}

load.meta.data <- function(path) 
{
    meta_env <- new.env()
    load(path, envir = meta_env)
    meta <- subset(meta_env$meta_data,
        select = c("aid", "sex", "date_birth", "date_first_positive", "date_last_negative")
    )
    meta <- unique(meta[!is.na(aid)])
    stopifnot(meta[, uniqueN(aid) == .N])
    return(meta)
}


get.household.data.marco <- function() {
    path.flow <- file.path(indir.deepsequencedata, "RCCS_R15_R18", "FlowR15_R18_VoIs_220129.csv")
    cols <- c("study_id", "round", "region", "comm_num", "hh_num", "member_num")
    flow <- fread(path.flow, select = cols)
    names(flow) <- toupper(names(flow))
    flow <- unique(flow[STUDY_ID != "", STUDY_ID := paste0("RK-", STUDY_ID)])
    flow <- merge(flow, aik, by.x = "STUDY_ID", by.y = "PT_ID")
    flow[, COMM_ID := paste(COMM_NUM, HH_NUM, sep = "_")]
    flow <- subset(flow, select = c("AID", "ROUND", "COMM_ID", "COMM_NUM", "HH_NUM"))

    # pairs F -> M  with big age differences
    dmother <- dresults[SEX.SOURCE == "F" & SEX.RECIPIENT == "M" & AGE_TRANSMISSION.SOURCE > AGE_INFECTION.RECIPIENT + 10, ]
    tmp <- double.merge(dmother[, .(RECIPIENT, SOURCE, AGE_INFECTION.RECIPIENT, AGE_TRANSMISSION.SOURCE, ROUND.M)], flow)
    cols <- names(tmp)[names(tmp) %like% "ROUND"]
    .f <- function(x) as.integer(gsub("R0|S", "", x))
    tmp[, (cols) := lapply(.SD, .f), .SDcols = cols]
    tmp[is.na(ROUND.SOURCE)]

    tmp[, length(intersect(COMM_ID.SOURCE, COMM_ID.RECIPIENT)), by = c("SOURCE", "RECIPIENT", "ROUND.M")]
    tmp[, length(intersect(COMM_NUM.RECIPIENT, COMM_NUM.SOURCE)), by = c("SOURCE", "RECIPIENT")]
    tmp
}

find_ff_pairs_for_Griffin <- function() { # study all FF pairs to be sent to Griffin
    idx <- unique(dchain[, .(H1, H2)])
    idx <- merge(idx, meta[, .(H1 = aid, SEX.H1 = sex)], by = "H1")
    idx <- merge(idx, meta[, .(H2 = aid, SEX.H2 = sex)], by = "H2")
    idx <- tmp[SEX.H1 == SEX.H2 & SEX.H1 == "F", .(H1, H2)]

    cols <- c("SOURCE", "RECIPIENT", "SEX.SOURCE", "SEX.RECIPIENT", "ROUND.M", "DIRECTION", "COMM.SOURCE", "COMM.RECIPIENT")
    tmp12 <- merge(dresults, idx, by.x = c("SOURCE", "RECIPIENT"), by.y = c("H1", "H2"))
    tmp21 <- merge(dresults, idx, by.x = c("SOURCE", "RECIPIENT"), by.y = c("H2", "H1"))

    tmpUN <- merge(idx, tmp12, by.x = c("H1", "H2"), by.y = c("SOURCE", "RECIPIENT"), all.x = TRUE)[is.na(SEX.SOURCE), .(H1, H2)]
    tmpUN <- merge(tmpUN, tmp21, by.x = c("H1", "H2"), by.y = c("RECIPIENT", "SOURCE"), all.x = TRUE)[is.na(SEX.RECIPIENT), .(H1, H2)]

    tmpUN[, `:=`(SOURCE = H1, RECIPIENT = H2, SEX.SOURCE = "F", SEX.RECIPIENT = "F", ROUND.M = NA_character_, DIRECTION = "phyloscanner_unclear")]
    tmpUN[, `:=`(H1 = NULL, H2 = NULL)]
    stopifnot(nrow(tmpUN) + nrow(tmp12) + nrow(tmp21) == nrow(idx))

    all_ff_pairs <- rbind(tmp12,
        tmp21,
        tmpUN,
        fill = TRUE
    )

    load(path.meta, envir = meta_env)
    dcomms <- subset(meta_env$meta_data, select = c("aid", "comm", "round"))
    dcomms <- unique(dcomms[!is.na(aid), ])

    dcomms <- dcomms[is.na(comm), comm := "neuro"]
    dcomms <- dcomms[, list(comm = fifelse(uniqueN(comm) == 1, yes = comm[1], no = NA_character_)), by = "aid"]

    all_ff_pairs[is.na(COMM.SOURCE), COMM.SOURCE := dcomms[aid == SOURCE, comm], by = "SOURCE"]
    all_ff_pairs[is.na(COMM.RECIPIENT), COMM.RECIPIENT := dcomms[aid == RECIPIENT, comm], by = "RECIPIENT"]

    filename <- file.path(indir.deepsequencedata, "RCCS_R15_R18", paste0("221117_all_ff_pairs.csv"))
    fwrite(all_ff_pairs, filename)
}
