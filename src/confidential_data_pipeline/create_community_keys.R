# create new community keys file by adding the communities discountinued from round 15

library(data.table)
library(here)

gitdir <- here()
source(file.path(gitdir, "config.R"))

c(  file.community.keys.incomplete,
    file.community.keys.aggregated.incomplete) |> file.exists() |> all() |> stopifnot()

#
# Create community keys
#

# load files communities used until 2023-07-05
community.keys <- fread(file.community.keys.incomplete)

# communities discountinued in round 15, inland, that joseph share in Feb 2023
tmp <- data.table(COMM_NUM_RAW = c(3, 18, 25, 54, 59, 81, 95, 103, 109, 177, 451, 760), 
                  COMM_NUM_A = rep('iii', 12))
community.keys <- rbind(community.keys, tmp)


# `save
file.name <- file.community.keys
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  write.csv(community.keys, file = file.name, row.names = F)
}else{
  cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")

#
# Create community keys aggregated
#


# load files communities used until 2023-07-05
community.keys.aggregated <- fread(file.community.keys.aggregated.incomplete)

# communities discountinued in round 15, inland, that joseph share in Feb 2023
max_index <- community.keys.aggregated[grepl('i-', comm_id), max(as.numeric(gsub('i-(.+)', '\\1', comm_id)))]
tmp <- data.table(comm_num = c(3, 18, 25, 54, 59, 81, 95, 103, 109, 177, 451, 760), 
                  comm_id = paste0('i-', (max_index+1):(max_index+12)))
community.keys.aggregated <- rbind(community.keys.aggregated, tmp)

# `save
file.name <- file.community.keys.aggregated
if(! file.exists(file.name))
{
  cat("\n Saving output file", file.name, "\n")
  write.csv(community.keys.aggregated, file = file.name, row.names = F)
}else{
  cat("\n Output file", file.name, "already exists\n")
}
cat("\n Done \n")

