# Function for plotting contact intensities patterns and marginal age contact intensities posterior distribution----

# plot four panels contact intensities patterns----
plot_contact_intensity_matrix = function(d22)
{
  set(d22, NULL,"gender_label", d22[, paste0(part.sex,' (Participant)')])
  set(d22, NULL,"alter_gender_label", d22[, paste0(cont.sex,' (Contacted)')])

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    facet_grid(alter_gender_label ~ gender_label) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

plot_contact_intensity_matrix_emprical = function(d22)
{
  set(d22, NULL,"gender_label", d22[, paste0(part.sex,' (Participant)')])
  set(d22, NULL,"alter_gender_label", d22[, paste0(cont.sex,' (Contacted)')])

  ggplot(d22) +
    geom_tile(aes(x = part.age, y = cont.age, fill = M)) +
    coord_equal() +
    facet_grid(alter_gender_label ~ label + gender_label) +
    scale_fill_continuous(type = "viridis", na.value = 'transparent') +
    # scale_fill_gradient2(low = "darkblue", high = "darkgreen", guide = "colorbar", na.value = 'white'
    # ) +
    theme_bw() +
    theme(legend.position = 'bottom')
}

# plot posterior age marginal contact intensities distirbution----
plot_age_marginal_contact_intensity_MF = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex) &
                         grepl('F', cont.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                           colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_marginal_contact_intensity_FM = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex) &
                         grepl('M', cont.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted", alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_intensity_MF = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('M', part.sex) &
                         grepl('F', cont.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
        pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                           colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}

plot_age_contact_intensity_FM = function(pred_cntct)
{
  pred_cntct <- subset(pred_cntct,
                       grepl('F', part.sex) &
                         grepl('M', cont.sex)
  )
  if (!is.integer(pred_cntct$part.age)) {
    pred_cntct$age1 = as.numeric(pred_cntct$part.age - 0.5)
    pred_cntct$age2 = as.numeric(pred_cntct$part.age + 0.5)
    pred_cntct$age1 = ifelse(pred_cntct$age1 < 10, paste0('0',pred_cntct$age1), pred_cntct$age1)
    pred_cntct$age2 = ifelse(pred_cntct$age2 < 10, paste0('0',pred_cntct$age2), pred_cntct$age2)
    pred_cntct[, part.age := paste0("participants' age:", age1, '-', age2)]
    ncol_nb = 6
  }else{
    pred_cntct$part.age = ifelse(pred_cntct$part.age < 10, paste0('0',pred_cntct$part.age), pred_cntct$part.age)
    pred_cntct$part.age = as.character(pred_cntct$part.age)
    pred_cntct[, part.age := paste0("participants' age:", part.age)]
    ncol_nb = 8
  }
  p_m <- ggplot(pred_cntct, aes(x = cont.age,col = part.age)) +
    pammtools:::geom_stepribbon(aes(ymin = CL, ymax = CU, fill = part.age), linetype = "dotted",alpha = 0.2,
                                colour = "transparent", show.legend = FALSE) +
    geom_step(aes(y = M)) +
    geom_step(aes(y = cntct_intensity_empirical), col = "black", linetype = "dashed") +
    facet_wrap(~ part.age, ncol = ncol_nb, scales = "free_y") +
    theme_bw() +
    guides(fill = "none", scales = "free_y") + theme(legend.position = "none") +
    xlab("Age of contacted individuals") +
    ylab("Estimated contact intensities with 95% ci")
  return(p_m)
}
