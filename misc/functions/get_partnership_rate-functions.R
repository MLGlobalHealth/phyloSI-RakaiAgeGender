obtain.contact.matrix = function(dataset, age)
{
  polymod.tab = dataset
  
  # Get ages and number of age classes
  n.age <- length(age)
  
  # Add node and record IDs to polymod.tab.list
  polymod.tab <- within(polymod.tab, {
    record.id <- with(construct.recordID(n = n.age), c( rec ))
    node.id   <- with(construct.nodeID  (n = n.age), c(node))
  })
  
  # Construct structure matrix R
  ord <- 2
  R <- construct.Rmat(n = n.age, order = ord, sex = TRUE)$R
  
  # Run model
  polymod.mod <- inla(
    y ~ 1 + f(node.id, model = "generic0", Cmatrix = R,
              # log-Gamma(1, 0.0001) prior on log-precision
              hyper = list(prec = list(prior = "loggamma", param = c(1, 0.0001))),
              rankdef = 3*ord^2, constr = TRUE, diagonal = 0.001),
    E = U,
    family = "nbinomial",
    data = polymod.tab,
    # Normal(0, 0.001) prior on intercept
    control.fixed = list(mean.intercept = 0, prec.intercept = 0.001),
    # Normal(0, 0.001) prior on log-dispersion parameter
    control.family = list(
      hyper = list(theta = list(prior = "gaussian", param = c(0, 0.001)))),
    # Integration strategy eb for faster computation
    control.inla = list(int.strategy = "eb"),
    # Compute linear predictor
    control.predictor = list(compute = TRUE, link = 1),
    control.compute = list(
      waic = TRUE,    # Compute WAIC
      cpo = TRUE,     # Compute cross-validated predictive measures
      config = TRUE)) # Enable posterior sampling
  
  # Show summary
  summary(polymod.mod)
  
  #
  # Post processing
  #
  
  # Add expected c and m to polymod.tab
  polymod.tab <- within(polymod.tab, {
    # Compute expected linear predictor
    linpred_mean <- polymod.mod$summary.linear.predictor[, "mean"]
    linpred_sd <- polymod.mod$summary.linear.predictor[, "sd"]
    linpred_l95 <- polymod.mod$summary.linear.predictor[, "0.025quant"]
    linpred_u95 <- polymod.mod$summary.linear.predictor[, "0.975quant"]
    # c = exp(linear predictor) = contact rate
    # Divide c by 1e6 to go back to original scale (deflates contact rate c)
    c <- exp(linpred_mean)/1e6
    c.sd <- exp(linpred_sd)/(1e6^2)
    c.l95 <- exp(linpred_l95)/(1e6)
    c.u95 <- exp(linpred_u95)/(1e6)
    # m = T*c = contact intensity
    m <- T*c
    # Remove linpred
    rm(linpred_mean, linpred_sd, linpred_l95, linpred_u95)
  })
  return(list(polymod.mod, polymod.tab))
}

# Function to construct node IDs in matrix form
construct.nodeID <- function(n) {
    # Not sex specific
    nn <- n*(n + 1)/2
    node <- matrix(0, nrow = n, ncol = n)
    LT <- lower.tri(node, diag = TRUE)
    node[ LT] <- 1:nn
    node[!LT] <- t(node)[!LT]
    return(list(node = node))
}

# Function to construct record IDs in matrix form
construct.recordID <- function(n) {
    # Not sex specific
    rec <- matrix(1:n^2, nrow = n, ncol = n)
    return(list(rec = rec))
  
}

# Function to construct R-matrix
construct.Rmat <- function(n, order = 2, sex = TRUE) {
  if (sex) {
    # Sex specific
    D.MM <- D.FF <- construct.Dmat(n, order, tri = TRUE )$D # MM and FF
    D.FM         <- construct.Dmat(n, order, tri = FALSE)$D # FM
    D <- bdiag(D.MM, D.FM, D.FF) # Block diagonal difference matrix D
    R <- crossprod(D)            # Structure matrix R
    R <- as(R, "dsCMatrix")      # R is symmetric sparse (column compressed) matrix
    return(list(D = D, R = R))
  } else {
    # Not sex specific
    D <- construct.Dmat(n, order, tri = TRUE)$D
    R <- crossprod(D)
    R <- as(R, "dsCMatrix")
    return(list(D = D, R = R))
  }
}

# Function to construct D-matrix
construct.Dmat <- function(n, order, tri = FALSE) {
  require(Matrix)
  In <- Diagonal(n)            # n x n matrix
  D0 <- diff(In, diff = order) # Differences for single vector
  D1 <- kronecker(In, D0)      # Differences in horizontal direction
  D2 <- kronecker(D0, In)      # Differences in   vertical direction
  if (tri) {
    # For lower triangular matrix (MM and FF)
    ri1.mat <- matrix(1:((n - order)*n), nrow = n - order, ncol = n)
    ri2.mat <- matrix(1:(n*(n - order)), nrow = n, ncol = n - order)
    ci.mat  <- matrix(1:n^2, nrow = n, ncol = n)
    ri1 <- ri1.mat[row(ri1.mat) >=  col(ri1.mat)         ] # Rows to keep in D1
    ri2 <- ri2.mat[row(ri2.mat) >= (col(ri2.mat) + order)] # Rows to keep in D2
    ci  <-  ci.mat[row(ci.mat)  >=  col(ci.mat)          ] # Columns to keep in D1 & D2
    D <- rbind(D1[ri1, ci], D2[ri2, ci])
    return(list(D0 = D0, D1 = D1, D2 = D2, ri1 = ri1, ri2 = ri2, ci = ci, D = D))
  } else {
    # For full matrix (FM)
    D <- rbind(D1, D2)
    return(list(D0 = D0, D1 = D1, D2 = D2, D = D))
  }
}

plot_crude_estimate <- function(data, age){
  
  # Get ages and number of age classes
  n.age <- length(age)
  sex <- unique(data$sex)
  
  # crude estimate - age
  z <- with(data, log(1 + y/(U)))
  z.range <- range(z, na.rm = TRUE)
  
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  image(age, age, matrix(z[0*n.age^2 + 1:n.age^2], n.age, n.age), zlim = z.range, col = cols, useRaster = TRUE, add = TRUE)
  box(lwd = 0.5)
  
  if(sex == 'M'){
    mtext("Participant's age (Male)", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
    mtext("Partner's age (Female)",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
    mtext(paste0("Reported by male participants"),     side = 3, adj = 0.5, line = -1.5, outer = TRUE)
  }
  
  if(sex == 'F'){
    mtext("Partner's age (Male)", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
    mtext("Participant's age (Female)",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
    mtext(paste0("Reported by female participants"),     side = 3, adj = 0.5, line = -1.5, outer = TRUE)
  }
  
}

plot_smooth_estimate <- function(fit, age){
  
  # Get ages and number of age classes
  n.age <- length(age)
  sex <- unique(fit$sex)
  
  # Smooth estimate
  z <- log(fit$c)
  z.range <- range(z)
  
  
  plot.new()
  plot.window(
    xlim = c(min(age) - 0.5, max(age) + 0.5),
    ylim = c(min(age) - 0.5, max(age) + 0.5))
  image  (age, age, matrix(        z[0*n.age^2 + 1:n.age^2],  n.age, n.age), zlim = z.range, col = cols, add = TRUE)
  contour(age, age, matrix(1e6*exp(z[0*n.age^2 + 1:n.age^2]), n.age, n.age), levels = euro.levs, labcex = 0.5, lwd = 0.2, add = TRUE)
  axis(side = 1, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  axis(side = 2, at = ticks, labels = ifelse(ticks %% 20 != 0, NA, ticks), lwd = 0.5)
  box(lwd = 0.5)
  
  
  if(sex == 'M'){
    mtext("Participant's age (Male)", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
    mtext("Partner's age (Female)",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
    mtext(paste0("Reported by male participants"),     side = 3, adj = 0.5, line = -1.5, outer = TRUE)
  }
  
  if(sex == 'F'){
    mtext("Partner's age (Male)", side = 1, adj = 0.5, line = -1.5, outer = TRUE)
    mtext("Participant's age (Female)",     side = 2, adj = 0.5, line = -1.5, outer = TRUE)
    mtext(paste0("Reported by female participants"),     side = 3, adj = 0.5, line = -1.5, outer = TRUE)
  }
}
