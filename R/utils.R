print_table <- function(table) print(knitr::kable(table))

`%which.like%` <- function(x, rgx){
    x[ data.table::`%like%`(x,rgx)]
}

.lu <- function(x){length(unique(x))}

.read <- function(x){
  if(grepl('.csv$', x)){return(as.data.table(read.csv(x)))}
  if(grepl('.rds$|.RDS$',x)){return(as.data.table(readRDS(x)))}
}

.aid2pt <- function(x, aik){
  if ( .lu(x) < length(x)){stop('Error: avoid repeated entries')}
  colnames(aik) <- tolower(colnames(aik))
  x <- data.table(aid=x)
  x <- merge(x, aik, all.x=T, by='aid')
  x$pt_id
}

.pt2aid <- function(x, aik){
  if (.lu(x) < length(x)){stop('Error: avoid repeated entries')}
  x <- data.table(pt_id=x)
  x <- merge(x, aik, all.x=T, by='pt_id')
  x$aid
}

.remove.spaces <- function(x) {
  x <- gsub(' ', '', x)
  x[which(x=='')] <- NA
  return(x)
}

.dates <- function(x)
{
  tmp <- as.integer(gsub('^[0-9]{2}-[A-Za-z]{3}-([0-9]{2})$', '\\1', x))
  pre <- fifelse(tmp >= 30, '19', '20')
  x <- gsub('-([0-9]{2})$', paste0('-',pre, '\\1'),x)
}

.vars.with.multiple.values <- function(dt, name)
{
  tmp <- dt[, .N, by=name][N>1, ]
  return(tmp)
}

print.which.NA <- function(dt,regex='', name='dt')
{
  cols <- colnames(dt); cols <- grep(regex, cols, value=T)
  .f <- function(x){sum(is.na(x))}
  tmp <- dt[, lapply(.SD, .f), .SDcols=cols]
  cols <- cols[which(tmp[1,] != 0)]
  cat('\nData.Table: ', name , ' | Numbers of rows: ', dt[, .N])
  tmp <- tmp[,
             {
               n <- names(.SD);
               cat('\n-------------------------------------- \nColumns with NA entries : # NA entries \n-------------------------------------- \n')
               lapply(seq_along(.SD),
                      FUN=function(i){ cat(n[[i]], ': ', .SD[[i]], '\n'); 0})
             }
             , .SDcols=cols]
}

.year.diff <- function(x, y)
{
  if (!is.Date(x)) {x <- as.Date(x, format = '%Y-%m-%d')}
  if (!is.Date(y)) {y <- as.Date(y, format = '%Y-%m-%d')}
  lubridate::time_length(difftime(x, y),"years")
}

.month.diff <- function(x, y)
{
  if (!is.Date(x)) {x <- as.Date(x, format = '%Y-%m-%d')}
  if (!is.Date(y)) {y <- as.Date(y, format = '%Y-%m-%d')}
  lubridate::time_length(difftime(x, y),"months")
}

ageyrs2agegp = function(x){
    out <- cut(x, breaks=c(15,25,35,50), include.lowest=T,right=F, labels=c('15-24','25-34','35-49')) 
    out <- factor(out, levels=c('Total','15-24','25-34','35-49'),labels=c('Total','15-24','25-34','35-49'))
    return(out)
}


print.which.NA <- function(dt,regex='SOURCE|RECIPIENT')
{
  cols <- colnames(dt); cols <- grep(regex, cols, value=T)
  .f <- function(x){sum(is.na(x))}
  tmp <- dt[, lapply(.SD, .f), .SDcols=cols]
  cols <- cols[which(tmp[1,] != 0)]
  tmp <- tmp[,
             {
               n <- names(.SD);
               cat('\n-------------------------------------- \nColumns with NA entries : # NA entries \n-------------------------------------- \n')
               lapply(seq_along(.SD),
                      FUN=function(i){ cat(n[[i]], ': ', .SD[[i]], '\n'); 0})
             }
             , .SDcols=cols]
}

.assign.code.meta <- function(x)
        fcase( x %like% 'Rakai_Pangea2_RCCS_Metadata_20221128.RData', 1)
.assign.code.chain <- function(x)
        fcase( x %like% '211220_phsc_phscrelationships_02_05_30_min_read_100_max_read_posthoccount_im_mrca_fixpd', 1)


write.to.tex <- function(DT, file){

    # write to file the "body" of the table
    tab_latex <- xtable::xtable(DT, type='latex') 
    print(tab_latex, 
        file=file,
        hline.after=c(),
        include.rownames=FALSE,
        include.colnames=FALSE,
        only.contents=TRUE,
        comment=FALSE)

    # however, we want to be able to quickly copy-paste the tabular env
    # so print a skeleton
    
    skeleton <- xtable::xtable(DT[0], type='latex') %>% print(hline.after=c(), comment=FALSE) 

    cmd_input <- sprintf("\n \\\\input{%s} \n \\\\end{tabular}", file.path("TODO", basename(file)))
    cmd <- sub('\\n \\\\end\\{tabular\\}', cmd_input, skeleton) 
    cat('\n', cmd, '\n')
    cmd
}
