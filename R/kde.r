
#' tabulated.baseline
#' 
#' tabulate the baseline intensity function.
#' 
#' @param x integer. number of cylinder samples.
#' @param postcode.field numeric.
#' @importFrom KernSmooth bkde2D
#' @importFrom KernSmooth dpik
#' @importFrom tidyr expand_grid
tabulated.baseline<-function(case.df, date.time.field = 'week'){
  x.range = ceiling(range(case.df$x))
  y.range = ceiling(range(case.df$y))
  bdw.x = dpik(case.df$x) / 4
  bdw.y = dpik(case.df$y) / 4
  gridsize.x = ceiling(abs(diff(x.range)) / bdw.x)
  gridsize.y = ceiling(abs(diff(y.range)) / bdw.y)
  # gridsize.x = floor(abs(diff(x.range) / mean(diff(case.df$x))) )
  # gridsize.y = floor(abs(diff(y.range) / mean(diff(case.df$y))) ) 
  cat("Bandwidths are: ", bdw.x, bdw.y, '\n')
  cat("Grid size: ", gridsize.x,'x' ,gridsize.y, '\n')
  kde = bkde2D(case.df[,c('x','y')],
             bandwidth = c(bdw.x, bdw.y),
             gridsize = c(gridsize.x, gridsize.y),
            )

  r = expand.grid(kde$x1, kde$x2)
  r$z = c(unlist(kde$fhat))
  colnames(r) = c('x', 'y', 'z')
  times = min(case.df[,date.time.field], na.rm = T):max(case.df[,date.time.field], na.rm = T)
  r = as.data.frame(expand_grid(r, times))
  colnames(r) = c('x', 'y', 'z', 't')

  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    Parameters
  },
  error = function(e){
    Parameters = cmle(case.df[,date.time.field], n.cycles=20, start=NULL, save.on.dir=TRUE)
    return(Parameters)
  })

  delta2 = diff(head(kde$x2, 2)) * diff(head(kde$x1, 2))
  T = lambda(r$t, Parameters)
  r$z = r$z * T
  C = nrow(case.df) / sum(r$z)
  r$z = C * r$z
  attributes(r)$delta2 = delta2 # = diff(head(baseline.tab[[2]]$x2), 2)) * diff(head(baseline.tab[[2]]$x1), 2))
  # diff(x.range) / gridsize.x * diff(y.range) / gridsize.y
  # attributes(r)$kde2D = kde
  return(r)
}


#' postcode.in.england
#' 
#' Check if postcode is in England.
#' 
#' @param x integer. number of cylinder samples.
#' @param postcode.field numeric.
#' @importFrom jsonlite read_json
compute.from.tab.baseline<-function(cylinder, observation.matrix, tab.baseline, postcode.locations){
  t.low = as.numeric(cylinder['t.low'])
  t.upp = as.numeric(cylinder['t.upp'])
  x0 = as.numeric(cylinder['x'])
  y0 = as.numeric(cylinder['y'])
  rho = as.numeric(cylinder['rho'])
  t.range = as.character(as.integer(cylinder['t.low']):as.integer(cylinder['t.upp']))
#  cat("t.range",t.range,"\n")

  observations = observation.matrix[,t.range]


  d = sqrt(
    (tab.baseline$x - x0)^2 +
    (tab.baseline$y - y0)^2
  )
  in_circle = (d < rho) & !is.na(d)
  in_height = (tab.baseline$t >= t.low) & (tab.baseline$t <= t.upp)


  # print('XXX')
  # print(any(in_circle))
  # print(any(in_height))
  # print(tab.baseline[1,])
  # print(cylinder)
  # print(head(d))
  mu = sum(tab.baseline[in_circle & in_height,]$z) #* attributes(tab.baseline)$delta2 #multiply by delta2 as this is the bin size
  
#  in_of_square = sum(in_circle & in_height) * delta2
#  out_of_square = pi * rho* rho - in_of_square
#  correct by taking into account that outsise the square there is nothing. 
#  mu = mu * pi * rho *rho / in_of_square
  
  d = sqrt(
      (as.numeric(postcode.locations$longitude) - x0)^2 +
      (as.numeric(postcode.locations$latitude)  - y0)^2
  )
  in_circle = (d < rho) & (!is.na(d))

  n_cases_in_cylinder = sum(observations[in_circle, ], na.rm = T)
  
#  ci = qpois(c(0.25,0.95) , lambda=mu)
  p.val = ppois(n_cases_in_cylinder-1, lambda=mu, lower.tail=FALSE)
  return (c(n_cases_in_cylinder, mu, p.val))
}





# radia = sapply(1:5, function(h){sqrt(1 / mean.baseline / pi / h )} )
# radia_and_heights = cbind(1:5, radia)
# n.cylinders = 10000
# radia_and_heights = radia_and_heights[sample(1:nrow(radia_and_heights), n.cylinders, replace=T),]
# rho = radia_and_heights[,2]
# random_radia = runif(n.cylinders, 0, rho)
# theta = runif(n.cylinders, 0, 2* pi)
# idx = sample(1:nrow(observations), n.cylinders, replace=T)
# y0 = observations$y[idx] + sin(theta) * random_radia
# x0 = observations$x[idx] + cos(theta) * random_radia
# t = observations$t[idx] + round(runif(n.cylinders, -radia_and_heights[,1]/2, radia_and_heights[,1]/2))
# t.low = t - round(radia_and_heights[,1]/2)
# t.low = ifelse(t.low>0, t.low, 1)
# t.upp = t.low + round(radia_and_heights[,1]/2) 
# t.max = max(observations$t)
# t.upp = ifelse(t.upp <= t.max, t.upp, t.max)
# cylinders = data.frame(x=x0, y=y0, rho=rho, t.low=t.low, t.upp=t.upp, original.x=observations$x[idx],  original.y=observations$y[idx])


# cylinders[,c('n_obs', 'mu',  'p.val')] = t(apply(cylinders, 1, compute_cylinders,
#                                                   observations, tabulated.baseline))
# # cylinders$warning = apply(cylinders, 1, function(x){ifelse((x['p.val'] < 0.05) & (x['n_obs'] > 0), TRUE, FALSE)})
# cylinders$warning = (cylinders['p.val'] < 0.05) & (cylinders['n_obs'] > 0)



ppp=function(case.df, baseline.tab, x0=0){
  idx.case.df = (case.df$week==1) & (case.df$x < x0) & (case.df$y > 5900)
  print(sum(idx.case.df))
  idx.case.df = (case.df$week==1) & (case.df$x >= x0) & (case.df$y > 5900)
  print(sum(idx.case.df))
  print(sum(case.df$week==1))
  print(" ")
  idx = (baseline.tab$x < x0) & (baseline.tab$t == 1) & (baseline.tab$y > 5900)
  print(sum(baseline.tab[idx,]$z))
  idx = (baseline.tab$x >= x0) & (baseline.tab$t == 1) & (baseline.tab$y > 5900)
  print(sum(baseline.tab[idx,]$z))
  print(sum(baseline.tab[baseline.tab$t == 1,]$z))
}




#' tabulated.baseline
#' 
#' tabulate the baseline intensity function.
#' 
#' @param x integer. number of cylinder samples.
#' @param postcode.field numeric.
#' @importFrom KernSmooth dpik
#' @importFrom tidyr expand_grid
tabulated.baseline2<-function(case.df, date.time.field = 'week'){
  x.range = range(case.df$x)
  x.range = c(floor(x.range[1]), ceiling(x.range[2]))
  y.range = range(case.df$y)
  y.range = c(floor(y.range[1]), ceiling(y.range[2]))
  bdw.x = dpik(case.df$x)
  bdw.y = dpik(case.df$y)
  gridsize.x = ceiling(abs(diff(x.range)) / bdw.x)
  gridsize.y = ceiling(abs(diff(y.range)) / bdw.y)
  x.bin <- seq(x.range[1], x.range[2], length=gridsize.x+1)
  y.bin <- seq(y.range[1], y.range[2], length=gridsize.y+1)

  # gridsize.x = floor(abs(diff(x.range) / mean(diff(case.df$x))) )
  # gridsize.y = floor(abs(diff(y.range) / mean(diff(case.df$y))) ) 
  cat("Bandwidths are: ", bdw.x, bdw.y, '\n')
  cat("Grid size: ", gridsize.x,'x' ,gridsize.y, '\n')

  freq = as.data.frame(table(findInterval(case.df$x, x.bin),findInterval(case.df$y, y.bin)))
  freq[,1] = as.numeric(freq[,1])
  freq[,2] = as.numeric(freq[,2])

  x.bin.centers = x.bin[-1] - diff(x.bin[1:2]) / 2
  y.bin.centers = y.bin[-1] - diff(y.bin[1:2]) / 2 

  freq2D = matrix(data=0, ncol=gridsize.x, nrow=gridsize.y)
  freq2D[cbind(freq[,2], freq[,1])]=freq[,3]

  r1 = expand.grid(y.bin.centers, x.bin.centers)
  r1$z = c(unlist(freq2D))
  colnames(r1) = c('y', 'x', 'z')
  times = min(case.df[,date.time.field], na.rm = T):max(case.df[,date.time.field], na.rm = T)
  r1 = as.data.frame(expand_grid(r1, times))
  colnames(r1) = c('y', 'x', 'z', 't')

  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    Parameters
  },
  error = function(e){
    Parameters = cmle(case.df[,date.time.field], n.cycles=20, start=NULL, save.on.dir=TRUE)
    return(Parameters)
  })

  T = lambda(r1$t, Parameters)
  r1$z = r1$z * T
  C = nrow(case.df) / sum(r1$z)
  r1$z = C * r1$z

  return(list(r=r,r1=r1))
}
