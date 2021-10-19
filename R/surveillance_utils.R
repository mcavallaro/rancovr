library(truncnorm)

f_radia_and_heights<-function(baseline.matrix, weeks=1:100){
  # this function compute the optimal radia and no. of weeks, such that the
  # cylinders that these define contain one event in average, give the baseline.
  n.weeks = ncol(baseline.matrix)
  mean.baseline = mean(baseline.matrix[, 1:n.weeks])

  # R = 3959 # earth radius in miles
  # lat1 = min(postcode2coord$latitude)
  # lat2 = max(postcode2coord$latitude)
  # lon1 = min(postcode2coord$longitude)
  # lon2 = max(postcode2coord$longitude)
  # total.area = (pi / 180) * R * R  * abs(sin(lat1)-sin(lat2)) * abs(lon1-lon2)
  total.area = 130279  # total area of England in km2s
  mean.spatial.baseline = mean.baseline * nrow(baseline.matrix) / total.area
  
  # we want A such that A * mean.spatial.baseline=1,  A * mean.spatial.baseline=1/2,  A * mean.spatial.baseline=1/3, etc
  radia = sapply(weeks, function(x){sqrt(1 / mean.spatial.baseline / pi / x )} )
  return(cbind(weeks, radia))
}

rcylinder2<-function(n.cylinders, observation.matrix, week.range, radia_and_heights, postcode2coord){
  if (any(rownames(observation.matrix)=='NA')){
    cat("WARNING: any(rownames(observation.matrix)=='NA'")
  }
  cols = as.character(week.range[1]:week.range[2])
  cases = which(observation.matrix[, cols] > 0, arr.ind = T)
  # cases is of the form:
  #           row col
  # PL15 9NE 5851 231
  # SY11 3PN 6370 255
  if (sum(observation.matrix[, cols]) > 0){
    idx = sample(1:nrow(cases), n.cylinders, replace = T)
    y = postcode2coord[cases[idx,1], 'latitude']
    x = postcode2coord[cases[idx,1], 'longitude']
    t = cases[idx,2] + week.range[1] - 1
    # radia and heights are given as input in the matrix radia_and_heights

    radia_and_heights = radia_and_heights[sample(1:nrow(radia_and_heights), n.cylinders, replace=T),]
    #randomise wilst keeping same radius and height and avoiding negative t
    rho = radia_and_heights[,2]
    
    random_radia = runif(n.cylinders, 0, rho)
    theta = runif(n.cylinders, 0, 2* pi)

    y = y + sin(theta) * random_radia
    x = x + cos(theta) * random_radia
    tt = t
    t = t + round(runif(n.cylinders, -radia_and_heights[,1]/2, radia_and_heights[,1]/2))
    t.low = t - round(radia_and_heights[,1] / 2)
    t.min = as.integer(week.range[1])
    t.low = ifelse(t.low >= t.min, t.low, t.min)

    t.upp = tt + round(radia_and_heights[,1] / 2)
    t.max = as.integer(week.range[2])
    t.upp = ifelse(t.upp <= t.max, t.upp, t.max)

    t.low = as.integer(ifelse(t.low == t.max, t.low - 1, t.low))
    
    return(data.frame(x=x, y=y, rho=rho, t.low=t.low, t.upp=t.upp))
  }else{
    return(data.frame(x=double(), y=double(), rho=double(), t.low=integer(), t.upp=integer()))    
  }
}



rcylinder<-function(n,  X.range, Y.range, time.range, border=NULL, rs=0.5, a=0){
  # n is the number of samples
  # times
  # border is the distance from the boundary
  # rs is a reference size 
  
  X.min = X.range[1]
  X.max = X.range[2]
  Y.min = Y.range[1]
  Y.max = Y.range[2]
  
  if (is.null(border)){
    border = (X.max - X.min) / 100
  }
  
  # generate centers in a square:
  x = runif(n, X.min + border, X.max - border)
  y = runif(n, Y.min + border, Y.max - border)
  
  # the spatio-temporal portions und at time.range[1]
  t.upp = time.range[2]
  t.low = sample(time.range[1]:(time.range[2]-1), n, replace=T)
  
  
  # generate circle radia
  d_max = apply(data.frame(x=x-X.min, y=y-Y.min, xs=X.max-x, ys=Y.max-y), 1, min) # this is min vertical distance from the border.
  rho = rtruncnorm(n, mean=rs, sd=rs, a=a, b=d_max)
  return (data.frame(x=x, y=y, rho=rho, t.low=t.low, t.upp=t.upp))
}

is_in_circle<-function(Data, x, y, rho){
  d = sqrt((as.numeric(Data['longitude'])-x)^2 + (as.numeric(Data['latitude'])-y)^2 )
  res<-(d<rho) & (!is.na(d))
  return(res)
}

compute<-function(cylinder, observation.matrix, baseline.matrix, postcode.locations){
  # cylinder is a spatio-temporal portion
  # extract Postcodes and times contained in cyclinders
  t.range = as.character(as.integer(cylinder['t.low']):as.integer(cylinder['t.upp']))
  observations = observation.matrix[,t.range]
  baselines = baseline.matrix[,t.range]
  
  # # da vettorizzare:
  # in_circle<-apply(postcode.locations, 1, is_in_circle,
  #                  as.numeric(cylinder['x']), as.numeric(cylinder['y']), as.numeric(cylinder['rho']))
  # vettorizzato
  d = sqrt(
    (as.numeric(postcode.locations$longitude) - as.numeric(cylinder['x']))^2 +
      (as.numeric(postcode.locations$latitude)  - as.numeric(cylinder['y']))^2
  )
  in_circle = (d<as.numeric(cylinder['rho'])) & (!is.na(d))
  n_cases_in_cylinder = sum(observations[in_circle, ], na.rm=T)
  # postcodes<-paste(rownames(observations[in_circle, ]), collapse = ",")
  # the sum of Poisson RVs is Poisson
  mu = sum(baselines[in_circle, ])
  #mu = sum(baselines[in_circle, ]) + 1
  
  #  ci = qpois(c(0.25, 0.95) , lambda=mu)
  p.val = ppois(n_cases_in_cylinder-1, lambda=mu, lower.tail=FALSE)
  return (c(n_cases_in_cylinder, mu, p.val))
}

warning_ratio<-function(i, observation.matrix, cylinders, postcode.locations){
  # check if the location i
  x = as.numeric(postcode.locations[i, 'longitude'])
  y = as.numeric(postcode.locations[i, 'latitude'])
  
  ## da vettorizzare
  # in_circle = apply(cylinders, 1, function(X){
  #   ifelse(sqrt((as.numeric(X['x']) - x)^2 + (as.numeric(X['y']) - y)^2) < as.numeric(X['rho']), TRUE, FALSE)
  # })
  
  # vettorizzato
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

warning_ratio2<-function(i, observation.matrix, t, cylinders, postcode.locations){
  times = which(observation.matrix > 0)
  if (length(times) > 0){
    # check if the location is in any circle
    x = as.numeric(postcode.locations[i,'longitude'])
    y = as.numeric(postcode.locations[i,'latitude'])
    
    
    in_circle = apply(cylinders, 1, function(X){
      ifelse(sqrt((as.numeric(X['x']) - x)^2 + (as.numeric(X['y']) - y)^2) < as.numeric(X['rho']), TRUE, FALSE)
    })
    in_cylinder_height = (cylinders['t.low'] < t) & (cylinders['t.upp'] > t)
    # number of cylinders that include locations
    in_cylinder = sum(in_circle * in_cylinder_height)
    # number of cylinder with `warning` flag that include location `i`
    warning = sum(cylinders$warning * in_circle * in_cylinder_height)
    if (in_cylinder == 0){
      return(0)
    }else{
      return(warning / in_cylinder)
    }
  }else{
    return(0)
  }
}


warning.score<-function(case, cylinders, date.time.field = 'SAMPLE_DT_numeric'){
  x = as.numeric(case['x'])
  y = as.numeric(case['y'])
  TT = as.numeric(case[date.time.field])
  
  d = sqrt((cylinders$x - x)^2 + (cylinders$y - y)^2)
  in_circle = as.integer(d <= cylinders$rho)
  in_cylinder_height = as.integer((cylinders$t.low <= TT) & (cylinders$t.upp >= TT))
  
  # number of cylinders that include geo-coordinate of `case`
  in_cylinder = sum(in_circle * in_cylinder_height, na.rm=T)

  # number of cylinder with `warning` flag that include location `i`
  warning = sum(cylinders$warning * in_circle * in_cylinder_height, na.rm=T)
  if (in_cylinder>0){
    re = warning / in_cylinder
  }else{
    re = 0
  } 
  return(re)
}

