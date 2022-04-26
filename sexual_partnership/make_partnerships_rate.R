library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)

indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'

infile.reported.partnerships <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220426.csv')

# cont.age: contacted age
# part.age: participant age
# part.sex: participant sex
# cont.sex: contacted sex
# y: number of partnership reported by participant
# N: Number of participants
# T: population count of contacted
# U = N x T

# round: participant round
# comm: participant community = inland, fishing

# load data
reported.partnerships <- as.data.table(read.csv(infile.reported.partnerships))

