
#' Find cylinders' radii and heights
#'
#' This function returns cylinder radii and heights, such that the corresponding
#' cylinder volumes have determined characteristics.
#' For the chosen N heights (\code{heights}),
#' it returns an Nx2 \code{Matrix} whose first column contains the heights and 
#' and the second column contains the corresponding radii.
#' It works in two modes: 1) if \code{baseline.matrix} and \code{target.n.event} are provided (but \code{rho} is not),
#' the function computes the volume (such that each cylinder contains around \code{target.n.event}
#' events given the baseline); this is the mode used in Cavallaro, M. et al (2022) PLoS Comput Biol 18(11): e1010726
#' https://doi.org/10.1371/journal.pcbi.1010726 .
#' 2) if \code{rho} is provided (but  \code{baseline.matrix} and \code{target.n.event} are not)
#' the functions computes the largest volume as $pi * rho^2$. All other radii
#' are computed with heights from \code{heights} such the the cylinder volumes are constant.
#' @param heigths \code{integer}. 
#' @param baseline.matrix A \code{Matrix} encoding the baseline. 
#' @param target.n.event An \code{integer}, representing the desired expected number of events in each cylinder. It requires a baselineline matrix to be passed. 
#' @param rho A \code{numeric}. It is the maximum  
#' @return A \code{Matrix}.
#' @examples
#' radii_and_heights = f_radii_and_heights2(1:16, baseline.matrix, 100)
#' radii_and_heights = f_radii_and_heights2(1:16, rho=1.6446127)
f_radii_and_heights2<-function(heigths, baseline.matrix, target.n.event, rho){
  if (missing(baseline.matrix) & missing(target.n.event) & !missing(rho)){
    print(sprintf("Reference distance is %f", rho))
    V = pi * rho^2
    radii = sapply(heigths, function(h){sqrt( V / pi / h)} )
  }else if(missing(rho) & !missing(baseline.matrix) & !missing(target.n.event)){
    total.area = 130279
    total.volume = total.area * ncol(baseline.matrix)
    mean.baseline = mean(baseline.matrix, na.rm = T)
    V = target.n.event / mean.baseline
    radii = sapply(heigths, function(h){sqrt( V / pi / h)} )
  }else{
    print('suca')
  }
  ret = cbind(heigths, radii)
  colnames(ret) = c('heights', 'radii')
  return(ret)
}



# LTLA[,c('y', 'x')] = vlatlong2km(LTLA[,c('latitude', 'longitude')])
# if( ('x' %in% colnames(coord.df)) & ('y' %in% colnames(coord.df)) ){
#   coordinates = coord.df[,c('y', 'x')]
# }
# else{
#   coordinates = vlatlong2km(coord.df[,c('latitude', 'longitude')])
# }
# rho = distance.between.points(LTLA[,'x'], LTLA[,'y']) 


f_radii_and_heights<-function(baseline.matrix, heigths=1:100, target.n.event, coord.df){
  total.area = 130279 # total area of England in km2s
  # R = 3959 # earth radius in miles
  if( ('x' %in% colnames(coord.df)) & ('y' %in% colnames(coord.df)) ){
    coordinates = coord.df[,c('y', 'x')]
  }
  else{
    coordinates = vlatlong2km(coord.df[,c('latitude', 'longitude')])
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
  radii = sapply(heigths, function(h){sqrt( V / pi / h)} )
  ret = cbind(heigths, radii)
  colnames(ret) = c('heights', 'radii')
  return(ret)
}


f_radii_and_heights.<-function(baseline.matrix, heigths=1:100, target.n.event, coord.df){
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
  radii = sapply(heigths, function(h){sqrt( V / pi / h)} )
  ret = cbind(heigths, radii)
  colnames(ret) = c('heights','radii')
  return(ret)
}

#' This funcion is analolgous to f_radii_and_heights2, but is takes a tabular baseline
#' instead of a baseline matrix.
#' @inheritParams f_radii_and_heights2
f_radii_and_heights_<-function(baseline.tab, heights=1:100){
  n.heights = length(unique(baseline.tab$t))
  mean.baseline = mean(baseline.tab$z[baseline.tab$z>0]) # z strictly > 0 to exclude the grid points on unpopulated area.
  
  total.area = 130279 # total area of England in km2
  idx = baseline.tab$z > 0
  n.points = sum(idx) / n.heights
  delta2 = prod(attributes(baseline.tab)$delta2)
  mean.spatial.baseline = mean.baseline / (n.points * delta2) * total.area
  # we want A such that A * mean.spatial.baseline=1,  A * mean.spatial.baseline=1/2,  A * mean.spatial.baseline=1/3, etc
  radii = sapply(heights, function(x){sqrt(1 / mean.spatial.baseline / pi / x )} )
  ret = cbind(heigths, radii)
  colnames(ret) = c('heights','radii')
  return(ret)
}

