# B splines
bspline = function(x, k, order, intervals)
{
  if (order == 1)
  {
    return(x >= intervals[k] & x < intervals[k + 1])
  }
  w1 <- 0
  w2 <- 0
  if (intervals[k] != intervals[k + order - 1])
    w1 <- (x - intervals[k]) / (intervals[k + order - 1] - intervals[k])
  if (intervals[k + 1] != intervals[k + order])
    w2 <- 1 - (x - intervals[k + 1]) / (intervals[k + order] - intervals[k + 1])

  spline <-  w1 * bspline(x, k, order - 1, intervals) + w2 * bspline(x, k + 1, order - 1, intervals)
  return(spline)
}

find_intervals = function(knots, degree, repeating = T)
{
  K <- length(knots)
  intervals <- vector(mode = 'double', length = 2*degree + K)
  # support of knots
  intervals[(degree + 1):(degree + K)] <- knots

  # extreme
  if (repeating)
  {
    intervals[1:degree] = min(knots)
    intervals[(degree + K + 1):(2*degree + K)] <- max(knots)
  } else {
    gamma <- 0.1
    intervals[1:degree] <- min(knots) - gamma*degree:1
    intervals[(degree + K + 1):(2*degree + K)] <- max(knots) + gamma*1:degree
  }
  return(intervals)
}

bsplines = function(data, knots, degree)
{
  K <- length(knots)
  num_basis <- K + degree - 1
  intervals <- find_intervals(knots, degree)
  m <- matrix(nrow = num_basis, ncol = length(data), 0)

  for (k in 1:num_basis)
  {
    m[k,] = bspline(data, k, degree + 1, intervals)
  }

  m[num_basis,length(data)] <- 1

  return(m)
}
