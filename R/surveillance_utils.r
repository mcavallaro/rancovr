
distance.between.points<-function(x,y, n.neighbours=10){
  n.rows = length(x)
  n.cols = length(y)
  M = c()
  # M = matrix(NA, ncol = n.cols, nrow=n.rows)
  for (i in 1:n.rows){
    d = sqrt((x[i]-x)^2 + (y[i]-y)^2)
    # M[,i] = d
    d = sort(d)[1:n.neighbours]
    M = c(M, d)
  }
  return(quantile(M, 0.95))
}


#' Find optimal radia
#'
#' This function compute the optimal cylinder radia, such that the corresponding
#' cylinders contain one event in average 
#' for the chosen N heights (\code{heights}) and a given \code{baseline.matrix}.
#' It returns an Nx2 \code{Matrix} whose first column contains the heights and 
#' and the second column contains the corresponding radia.
#' @param baseline.matrix A \code{Matrix} enconding the baseline. 
#' @param heigths \code{integer}. 
#' @return A \code{Matrix}.
#' @examples
#' radia_and_heights = f_radia_and_heights(baseline.matrix, 1:10)
f_radia_and_heights<-function(baseline.matrix, heigths=1:100, target.n.event, coord.df){
  total.area = 130279 # total area of England in km2s
  # R = 3959 # earth radius in miles
  if(T){
    coordinates = vlatlong2km(coord.df[,c('latitude', 'longitude')])
  }
  else{
    coordinates = coord.df[,c('y', 'x')]
  }
  rho = distance.between.points(coordinates[,'x'], coordinates[,'y']) 
  print(sprintf("Reference distance is %f", rho))
  # a=diff(range(coordinates[,1], na.rm=T)) 
  # b=diff(range(coordinates[,2], na.rm=T)) 
  # total.area = a * b
  # total.n.events = sum(baseline.matrix, na.rm=T)
  # total.volume = total.area * ncol(baseline.matrix)
  # mean.baseline = mean(baseline.matrix, na.rm = T)
  # V = total.volume / total.n.events * target.n.event
  V = pi * rho^2
  radia = sapply(heigths, function(h){sqrt( V / pi / h)} )
  ret = cbind(heigths, radia)
  colnames(ret) = c('heights','radia')
  return(ret)
}



#' Find optimal radia
#'
#' This function compute the optimal cylinder radia, such that the corresponding
#' cylinders contain one event in average 
#' for the chosen N heights (\code{heights}) and a given \code{baseline.matrix}.
#' It returns an Nx2 \code{Matrix} whose first column contains the heights and 
#' and the second column contains the corresponding radia.
#' @param baseline.matrix A \code{Matrix} enconding the baseline. 
#' @param heigths \code{integer}. 
#' @return A \code{Matrix}.
#' @examples
#' radia_and_heights = f_radia_and_heights(baseline.matrix, 1:10)
f_radia_and_heights.<-function(baseline.matrix, heigths=1:100, target.n.event, coord.df){
  total.area = 130279 # total area of England in km2s
  # R = 3959 # earth radius in miles
  # if(T){
  #   coordinates = vlatlong2km(coord.df[,c('latitude', 'longitude')])
  # }
  # else{
  #   coordinates = coord.df[,c('y', 'x')]
  # }
  # a=diff(range(coordinates[,1], na.rm=T)) 
  # b=diff(range(coordinates[,2], na.rm=T)) 
  # total.area = a * b
  # total.n.events = sum(baseline.matrix, na.rm=T)
  total.volume = total.area * ncol(baseline.matrix)
  mean.baseline = mean(baseline.matrix, na.rm = T)
  # V = total.volume / total.n.events * target.n.event
  V = target.n.event / mean.baseline
  radia = sapply(heigths, function(h){sqrt( V / pi / h)} )
  ret = cbind(heigths, radia)
  colnames(ret) = c('heights','radia')
  return(ret)
}


#' Find optimal radia
#'
#' This function compute the optimal cylinder radia, such that the corresponding
#' cylinders contain one event in average 
#' for the chosen N heights (\code{heights}) and a given tabulated baseline.
#' It returns an Nx2 \code{Matrix} whose first column contains the heights and 
#' and the second column contains the corresponding radia.
#' @param baseline.matrix An \code{expand.grid} tab enconding the baseline. 
#' @param heigths \code{integer}. 
#' @return A \code{Matrix}.
#' @examples
#' radia_and_heights = f_radia_and_heights(baseline.matrix, 1:10)
f_radia_and_heights_<-function(baseline.tab, heights=1:100){
  n.heights = length(unique(baseline.tab$t))
  mean.baseline = mean(baseline.tab$z[baseline.tab$z>0]) # z strictly > 0 to exclude the grid points on unpopulated area.

  total.area = 130279 # total area of England in km2
  idx = baseline.tab$z > 0
  n.points = sum(idx) / n.heights
  delta2 = prod(attributes(baseline.tab)$delta2)
  mean.spatial.baseline = mean.baseline / (n.points * delta2) * total.area
  # we want A such that A * mean.spatial.baseline=1,  A * mean.spatial.baseline=1/2,  A * mean.spatial.baseline=1/3, etc
  radia = sapply(heights, function(x){sqrt(1 / mean.spatial.baseline / pi / x )} )
  ret = cbind(heigths, radia)
  colnames(ret) = c('heights','radia')
  return(ret)
}


#' Draw random cylinder coordinates
#' 
#' Find the coordinates (centers and height limits) of cylinders
#' that contain events defined in \code{observation.matrix}, with radia and heights
#' given in \code{radia_and_heights}.
#' 
#' @param n.cylinders An \code{integer}; the number of cylinders to draw.
#' @param observation.matrix A \code{sparseMatrix} object enconding the events.
#' @param time.range An \code{integer} vector.
#' @param radia_and_heights A \code{Matrix}.
#' @param postcode2coord A \code{data.frame} that maps the rows of \code{observation.matrix} to geographical coordinates.
#' @param only.last A \code{bool}; true if all cylinders must include the last time point (\code{time.range[2]}). For prospective analysis.
#' @importFrom truncnorm rtruncnorm
#' @import Matrix
#' @return A \code(data.frame).
#' @examples
#' cylinders=rcylinder(10, observation.matrix, time.range, radia_and_heights, postcode2coord)
rcylinder<-function(n.cylinders, observation.matrix, time.range, radia_and_heights, postcode2coord, only.last=F){
  if (any(rownames(observation.matrix)=='NA')){
    cat("WARNING: any(rownames(observation.matrix)=='NA'")
  }
  if (only.last){
    cols = as.character(time.range[2])
  }else{
    cols = as.character(time.range[1]:time.range[2])
  }
  cases = which(observation.matrix[, cols] > 0, arr.ind = T)
  # cases is of the form:
  #           row col
  # PL15 9NE 5851 231
  # SY11 3PN 6370 255
  if (sum(observation.matrix[, cols]) > 0){
    idx = sample(1:NROW(cases), n.cylinders, replace = T)
    if (is.matrix(cases)){
      idx_row = cases[idx, 1]
      idx_col = cases[idx, 2]
    }else{
      idx_row = cases[idx]
      idx_col = T
    }
    y = postcode2coord[idx_row, 'y']
    x = postcode2coord[idx_row, 'x']
    # t = cases[idx,2] + time.range[1] - 1
    # print(head(t))
    tt = as.integer(colnames(observation.matrix[,cols])[idx_col])

    # radia and heights are given as input in the matrix radia_and_heights

    radia_and_heights = radia_and_heights[sample(1:nrow(radia_and_heights), n.cylinders, replace=T),]
    # randomise whilst keeping same radius and height and avoiding negative t
    rho = radia_and_heights[, 2]
    random_radia = runif(n.cylinders, 0, rho)
    theta = runif(n.cylinders, 0, 2* pi)
    y = y + sin(theta) * random_radia
    x = x + cos(theta) * random_radia
    
    t.min = as.integer(time.range[1])
    t.max = as.integer(time.range[2])
    
    if(only.last){
      t.upp = rep(t.max, length(x))
      t.low = t.upp - as.integer(radia_and_heights[,1])
      t.low = ifelse(t.low < t.min, t.min, t.low)
    }else{
      rrr = v.sample.int(as.integer(radia_and_heights[,1]) + 1, 1) - 1
      t.low = tt - rrr
      t.upp = t.low + as.integer(radia_and_heights[,1])
      t.upp = ifelse(t.low > t.min, t.upp, t.upp + (t.min - t.low))    
      t.low = ifelse(t.low > t.min, t.low, t.min)
      t.low = ifelse(t.upp < t.max, t.low, t.low - (t.upp - t.max) )
      t.upp = ifelse(t.upp < t.max, t.upp, t.max)
      t.low = as.integer(ifelse( (t.upp==t.max) & (t.low < t.min), t.min, t.low))      
    }
    # tt = t + runif(n.cylinders, -radia_and_heights[,1]/2, radia_and_heights[,1]/2)
    # if (only.last){
    #   t.upp = as.integer(time.range[2])
    # }else{
    #   t.upp = ceiling(tt + radia_and_heights[,1] / 2)
    # }
    # t.max = as.integer(time.range[2])
    # if (only.last){
    #   t.low = floor(t.upp - radia_and_heights[,1])
    #   t.low = ifelse(t < t.low, t, t.low)
    # }else{
    #   t.low = floor(tt - radia_and_heights[,1] / 2)
    # }
    # t.min = as.integer(time.range[1])
    
    # t.low = ifelse(t.low >= t.min, t.low, t.min)
    # t.upp = ifelse(t.upp <= t.max, t.upp, t.max)
    
    # t.low = as.integer(ifelse(t.low == t.max, t.low - 1, t.low))
    # # isx = (t.low == t.upp)
    # # if (any(isx)){
    # #   print(c(t.low[isx], t.upp[isx]))
    # # }

    # # QUESTO FUNZIONA:
    # # t = cases[idx,2] + week.range[1] - 1
    # # # radia and heights are given as input in the matrix radia_and_heights

    # # radia_and_heights = radia_and_heights[sample(1:nrow(radia_and_heights), n.cylinders, replace=T),]
    # # #randomise wilst keeping same radius and height and avoiding negative t
    # # rho = radia_and_heights[,2]
    
    # # random_radia = runif(n.cylinders, 0, rho)
    # # theta = runif(n.cylinders, 0, 2* pi)

    # # y = yy + sin(theta) * random_radia
    # # x = xx + cos(theta) * random_radia
    # # tt = as.integer(t)
    # # # t = t + round(runif(n.cylinders, -radia_and_heights[,1]/2, radia_and_heights[,1]/2))
    # # v.sample.int<-Vectorize(sample.int, 'n')
    # # rrr = v.sample.int(as.integer(radia_and_heights[,1]) + 1, 1) - 1
    # # t.low = tt - rrr
    # # t.upp = t.low + as.integer(radia_and_heights[,1])
    # # # t.low = t - round(radia_and_heights[,1] / 2)
    # # t.min = as.integer(week.range[1])
    # # t.max = as.integer(week.range[2])
    
    # # t.upp = ifelse(t.low > t.min, t.upp, t.upp + (t.min - t.low))    
    # # t.low = ifelse(t.low > t.min, t.low, t.min)
    
    # # t.low = ifelse(t.upp < t.max, t.low, t.low - (t.upp - t.max) )
    # # t.upp = ifelse(t.upp < t.max, t.upp, t.max)
    
    # # # t.low = as.integer(ifelse( !(t.low == t.min) & (t.low == t.upp), t.low - 1, t.low))
    # # t.low = as.integer(ifelse( (t.upp==t.max) & (t.low < t.min), t.min, t.low))

    return(data.frame(x=x, y=y, rho=rho, t.low=t.low, t.upp=t.upp))
  }else{
    return(data.frame(x=double(), y=double(), rho=double(), t.low=integer(), t.upp=integer()))    
  }
}



#' Compute exceedance probabality in a cylinder.
#' 
#' A cylinder is defined by the circle coordinated (say, x,y, and radius) and lower and upper height limits (aay, t.low and t.upp, respectively).
#' For a given cylinder, this function computes the number of observed events (\code{n_cases}) in the cylinder according to
#' \code{observation.matrix}, the expected number \code{mu} of events according the Poisson point model (with intensity
#' defined in \code{baseline.matrix}), and the probability that .
#' The function returns \code{c(n_cases, mu, p.val)}.
#' 
#' @param cylinder 
#' @param observation.matrix A \code{sparseMatrix} object enconding the events.
#' @param baseline.matrix A \code{Matrix} object enconding the baseline.
#' @param postcode.locations A \code{data.frame} that maps the rows of \code{observation.matrix} to geographical coordinates.
#' @import Matrix
#' @return A \code{numeric} vector of dimension 3.
#' @examples
#' exceedance=compute(c(x,y,rho,t.low,t.upp), observation.matrix, baseline.matrix, postcode.locations)
compute<-function(cylinder, observation.matrix, baseline.matrix, postcode.locations){
  t.range = as.character(as.integer(cylinder['t.low']):as.integer(cylinder['t.upp']))
  observations = observation.matrix[,t.range]
  baselines = baseline.matrix[,t.range]

  d = sqrt(
    (as.numeric(postcode.locations$x) - as.numeric(cylinder['x']))^2 +
      (as.numeric(postcode.locations$y) - as.numeric(cylinder['y']))^2
  )
  in_circle = (d<as.numeric(cylinder['rho'])) & (!is.na(d))
  n_cases = sum(observations[in_circle, ], na.rm=T)
  # postcodes<-paste(rownames(observations[in_circle, ]), collapse = ",")
  # the sum of Poisson RVs is Poisson
  mu = sum(baselines[in_circle, ])
  #mu = sum(baselines[in_circle, ]) + 1
  
  #  ci = qpois(c(0.25, 0.95) , lambda=mu)
  p.val = ppois(n_cases-1, lambda=mu, lower.tail=FALSE)
  return (c(n_cases, mu, p.val))
}




#' Compute exceedance probabality in a cylinder.
#' 
#' A cylinder is defined by the circle coordinated (say, x,y, and radius) and lower and upper height limits (aay, t.low and t.upp, respectively).
#' For a given cylinder, this function computes the number of observed events (\code{n_cases}) in the cylinder according to
#' \code{observation.matrix}, the expected number \code{mu} of events according the Poisson point model (with intensity
#' defined in \code{tab.baseline}), and the probability that .
#' The function returns \code{c(n_cases, mu, p.val)}.
#' 
#' @param cylinder 
#' @param observation.matrix A \code{sparseMatrix} object enconding the events.
#' @param tab.baseline An \code{expand.grid} tab enconding the baseline.
#' @param postcode.locations A \code{data.frame} that maps the rows of \code{observation.matrix} to geographical coordinates.
#' @import Matrix
#' @return A \code{numeric} vector of dimension 3.
#' @examples
#' exceedance=compute(c(x,y,rho,t.low,t.upp), observation.matrix, tab.baseline, postcode.locations)
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





warning_ratio<-function(i, observation.matrix, cylinders, postcode.locations){
  # check if the location i
  x = as.numeric(postcode.locations[i, 'longitude'])
  y = as.numeric(postcode.locations[i, 'latitude'])
  
  in_circle = sqrt((as.numeric(cylinders['x']) - x)^2 + (as.numeric(cylinders['y']) - y)^2) < as.numeric(cylinders['rho'])
  
  times = which(observation.matrix[i,] > 0)

  if (length(times) > 0){
    in_cylinder_height = matrix(FALSE, nrow = nrow(cylinders), ncol=length(times))

    for (i in 1:length(times)){
      
      #da vettorizzare:
      in_cylinder_height[,i] = apply(cylinders, 1,  function(X){
        TT = as.numeric(names(times[i]))
        (as.numeric(X['t.low']) < TT) & (as.numeric(X['t.upp']) > TT)
      })
      
      
      # vettorizzato
      #TT = as.numeric(names(times[i]))
      #in_cylinder_height[,i] = (as.numeric(cylinders['t.low']) < TT) & (as.numeric(cylinders['t.upp']) > TT)

    }

    # number of cylinders that include locations `i`
    in_cylinder = sum(in_circle * in_cylinder_height)
    # number of cylinder with `warning` flag that include location `i`
    warning = sum(cylinders$warning * in_circle * in_cylinder_height)

    return(warning / in_cylinder)
  }else{
    return(NA)
  }
}
# 
# warning_ratio2<-function(i, observation.matrix, t, cylinders, postcode.locations){
#   times = which(observation.matrix > 0)
#   if (length(times) > 0){
#     # check if the location is in any circle
#     x = as.numeric(postcode.locations[i,'longitude'])
#     y = as.numeric(postcode.locations[i,'latitude'])
#     
#     
#     in_circle = apply(cylinders, 1, function(X){
#       ifelse(sqrt((as.numeric(X['x']) - x)^2 + (as.numeric(X['y']) - y)^2) < as.numeric(X['rho']), TRUE, FALSE)
#     })
#     in_cylinder_height = (cylinders['t.low'] < t) & (cylinders['t.upp'] > t)
#     # number of cylinders that include locations
#     in_cylinder = sum(in_circle * in_cylinder_height)
#     # number of cylinder with `warning` flag that include location `i`
#     warning = sum(cylinders$warning * in_circle * in_cylinder_height)
#     if (in_cylinder == 0){
#       return(0)
#     }else{
#       return(warning / in_cylinder)
#     }
#   }else{
#     return(0)
#   }
# }


warning.score<-function(case, cylinders, date.time.field = 'week'){
  x = as.numeric(case['x'])
  y = as.numeric(case['y'])
  TT = as.numeric(case[date.time.field])
  
  d = sqrt((cylinders$x - x)^2 + (cylinders$y - y)^2)
  in_circle = as.integer(d <= cylinders$rho)
  in_cylinder_height = as.integer((cylinders$t.low <= TT) & (cylinders$t.upp >= TT))
  
  # number of cylinders that include geo-coordinate of `case`
  in_cylinder = sum(in_circle * in_cylinder_height, na.rm=T)
  if (in_cylinder>0){
    # number of cylinder with `warning` flag that include location `i`
    warning = sum(cylinders$warning * in_circle * in_cylinder_height, na.rm=T)
    re = warning / in_cylinder
  }else{
    re = 0
  } 
  return(re)
}


warning.score2 <-function(case, TT, cylinders) {
  
  x = as.numeric(case['x'])
  y = as.numeric(case['y'])

  d = sqrt((cylinders$x - x)^2 + (cylinders$y - y)^2)
  in_circle = as.integer(d <= cylinders$rho)
  in_cylinder_height = as.integer((cylinders$t.low <= TT) & (cylinders$t.upp >= TT))
  
  # number of cylinders that include geo-coordinate of `case`
  in_cylinder = sum(in_circle * in_cylinder_height, na.rm=T)
  if (in_cylinder>0){
    # number of cylinder with `warning` flag that include location `i`
    warning = sum(cylinders$warning * in_circle * in_cylinder_height, na.rm=T)
    re = c(warning / in_cylinder,in_cylinder)
  }else{
    re = c(0, 0)
  } 
  return(re)
}




