library(data.table)
library(dplyr)
library(ggplot2)
library(scales)
library(lubridate)
library(reshape2)
library(INLA)

indir <- '~/git/phyloflows/sexual_partnership'
indir.deepsequencedata <- '~/Box\ Sync/2019/ratmann_pangea_deepsequencedata/live/'

infile.reported.partnerships <- file.path(indir.deepsequencedata, 'RCCS_R15_R18', 'RCCS_reported_partnership_220505.csv')

source(file.path(indir, 'functions', 'INLA-functions.R'))

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

# arrange
age <- reported.partnerships[, sort(unique(part.age))]
reported.partnerships <- reported.partnerships[order(part.round, part.comm, cont.sex, part.sex, cont.age, part.age)]

reported.partnerships[part.sex == 'F' & cont.sex == 'F']
# R015 in inland
tmp <- reported.partnerships[part.round == 'R015' & part.comm == 'inland']
contact.matrix <- obtain.contact.matrix(tmp, age)[[2]]

# plot
contact.matrix[, `Participant sex` := ifelse(part.sex == 'M', 'Male', 'Female')]
contact.matrix[, `Partner sex` := ifelse(cont.sex == 'M', 'Male', 'Female')]
contact.matrix[, `Participant sex` := factor(`Participant sex`, levels = c('Male', 'Female'))]
contact.matrix[, `Partner sex` := factor(`Partner sex`, levels = c('Female', 'Male'))]

plot_smooth_estimate(contact.matrix)
plot_smooth_estimate_sd(contact.matrix)



