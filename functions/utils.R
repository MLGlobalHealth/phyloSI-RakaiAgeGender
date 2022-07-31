print_table <- function(table) print(knitr::kable(table))

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

.f <- function(x){round(x, 1)}

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