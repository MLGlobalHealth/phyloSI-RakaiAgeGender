find_sensitivity_specificity_art <- function(rprev){
  
  # find percentage of participant who did not report art but had not viremic viral load
  tmp <- rprev[COMM == 'inland' & ART == F & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
                                                                N = length(STUDY_ID)), by = 'ROUND']
  tmp[, PROP := round(X / N * 100, 2)]
  tmp
  
  # find percentage of participant who report art and had not viremic viral load
  tmp <- rprev[COMM == 'inland' & ART == T & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
                                                                N = length(STUDY_ID)), by = 'ROUND']
  tmp[, PROP := round(X / N * 100, 2)]
  tmp
  
  if(0){ # plot
    df <- copy(rprev)
    df[, AGE_GROUP := '35-49']
    df[AGEYRS < 35, AGE_GROUP := '25-34']
    df[AGEYRS < 25,AGE_GROUP := '15-24' ]
    
    # find percentage of participant who did not report art but had not viremic viral load
    tmp <- df[COMM == 'inland' & ART == F & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
                                                               N = length(STUDY_ID)), by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp[, PROP := X / N ]
    tmp[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower, by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper, by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp1 <- df[COMM == 'inland' & ART == F & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
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
    tmp1 <- df[COMM == 'inland' & ART == T & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
                                                                N = length(STUDY_ID)), by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp1[, PROP := X / N ]
    tmp1[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower, by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp1[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper, by = c('ROUND', 'AGE_GROUP', 'SEX')]
    tmp <- df[COMM == 'inland' & ART == T & !is.na(VLNS), list(X = length(STUDY_ID[VLNS == 0]), 
                                                               N = length(STUDY_ID)), by = 'ROUND']
    tmp[, PROP := X / N ]
    tmp[, CL := binom::binom.confint(X, N, methods = 'agresti-coull')$lower, by = c('ROUND')]
    tmp[, CU := binom::binom.confint(X, N, methods = 'agresti-coull')$upper, by = c('ROUND')]
    
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

make_table_sensitivity_specificity_art <- function(rprev){
  
}
