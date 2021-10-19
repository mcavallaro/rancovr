library(jsonlite)

save.and.tell<-function(variable.string, file=file.name){
  save(list=variable.string, file=file, envir = parent.frame())
  writeLines(sprintf("The variable `%s` has been saved on disk in file `%s`.", variable.string, file))
  writeLines(sprintf('Load on memory with `load("%s", verbose=1)`.', file))  
}


split<-function(x){
  y = strsplit(x, ' ')[[1]]
  l = length(y)
  return(tolower(y[l]))
}

postcode.to.location<-function(x, postcodes){
  A = gsub(" ", "", toupper(postcodes['Patient.Postcode']),  fixed = TRUE)
  B = gsub(" ", "", toupper(x['Patient Postcode']),  fixed = TRUE)
  idx = A == B 
  s = sum(idx)
  if (is.na(s)){
    res = c(NA, NA)
  }
  else if(s == 1){
    res = unlist(postcodes[idx, c("latitude", "longitude")])
  }else{
    res = c(NA, NA)
  }
  return(res)
}

postcode.to.location.and.population<-function(x, postcodes){
  library(jsonlite)
  # A = gsub(" ", "", toupper(postcodes['Patient.Postcode']),  fixed = TRUE)
  B = gsub(" ", "", toupper(x['Patient Postcode']),  fixed = TRUE)
  idx = B == postcodes['Patient.Postcode']
  s = sum(idx)
  if (is.na(s)){
    res=c(NA,NA,NA)
  }
  else if(s == 1){
    res = unlist(postcodes[idx, c("latitude", "longitude", "Total")])
  }else{
    lat=read_json(paste0("http://api.getthedata.com/postcode/", B))$data$latitude
    if (class(lat) == 'character'){
      long=read_json(paste0("http://api.getthedata.com/postcode/", B))$data$longitude
      res = as.numeric(c(lat, long, 0))
    }else{
      res=c(NA,NA,NA)
    }
  }
  return(res)
}

postcode.in.england<-function(x, column='Patient Postcode'){
  jsn = read_json(paste0("http://api.getthedata.com/postcode/", gsub(" ", "", x[column],  fixed = TRUE)))
  country = jsn$data$country
  if (class(country) == 'character'){
    return(country == 'England')
  }else if (!is.null(jsn$error)){
    if (strsplit(jsn$error, ' ')[[1]][1] == "Northern"){
      return(FALSE)
    }else{
      return(NA)
    }
  }else{
    return(NA)
  }
}

postcode.to.location2<-function(x){
  B = gsub(" ", "", x['Patient Postcode'],  fixed = TRUE)
  B = paste0("http://api.getthedata.com/postcode/", B)
  longitude = read_json(B)$data$longitude
  latitude = read_json(B)$data$latitude
  if (class(longitude) == 'character'){
    return(c(as.numeric(latitude), as.numeric(longitude)))
  }else{
    return(c(NA, NA))
  }
}

postcode.to.location3<-function(x){
  B = gsub(" ", "", x,  fixed = TRUE)
  B = paste0("http://api.getthedata.com/postcode/", B)
  longitude = read_json(B)$data$longitude
  latitude = read_json(B)$data$latitude
  if (class(longitude) == 'character'){
    return(c(as.numeric(latitude), as.numeric(longitude)))
  }else{
    return(c(NA, NA))
  }
}

postcode.to.region<-function(x, Area2Region_list){
  # B = gsub(" ", "", x,  fixed = TRUE)
  # B = paste0("http://api.getthedata.com/postcode/", B)
  # postcode_area = read_json(B)$data$postcode_area
  # print(postcode_area)
  postcode_area = strsplit(x, ' ')[[1]][1]
  if ((class(postcode_area) == 'character')  & !(postcode_area == 'NA')){
    ret = Area2Region_list[[postcode_area]]
    if (is.null(ret)){
      B = gsub(" ", "", x,  fixed = TRUE)
      B = paste0("http://api.getthedata.com/postcode/", B)
      return(Area2Region_list[[read_json(B)$data$postcode_district]])
    }else{
      return(ret)
    }
  }else{
    return(c("NA"))
  }
}

week2Date<-function(x){
  x = as.numeric(x)
  date = as.Date("2015-01-01") + as.difftime(x, unit='weeks')
  date = unlist(lapply(date, as.character))
  date = unlist(lapply(date, function(x){
      tmp = strsplit(x, '-')[[1]]
      return(paste( tmp[3],tmp[2],substring(tmp[1],3,4), sep='/') )
    }
    ))
  return(date)
}

week2Year<-function(x){
  x = as.numeric(x)
  date = as.Date("2015-01-01") + as.difftime(x, unit='weeks')
  date = unlist(lapply(date, as.character))
  date = unlist(lapply(date, function(x){
    tmp = strsplit(x, '-')[[1]]
    return(tmp[1])
  }
  ))
  return(date)
}

week2MonthYear<-function(x){
  x = as.numeric(x)
  date = as.Date("2015-01-01") + as.difftime(x, unit='weeks')
  date = unlist(lapply(date, as.character))
  date = unlist(lapply(date, function(x){
    tmp = strsplit(x, '-')[[1]]
    return(paste(tmp[2],substring(tmp[1],3,4), sep='/') )
  }
  ))
  return(date)
}

latlong2km<-function(latitude, longitude){
  y = 111.2 * latitude
  x = 111.2 * longitude * cos(pi / 130* latitude)
  return(c(y, x))
}

km2latlong<-function(y, x){
  latitude = y / 111.2
  longitude = x / 111.2 / cos(pi / 130 * latitude)
  return(c(latitude, longitude))
}

vlatlong2km<-function(latlon){
  latitude = latlon[,1]
  longitude = latlon[,2]
  y = 111.2 * latitude
  x = 111.2 * longitude * cos(pi / 130 * latitude)
  res = cbind(y,x)
  colnames(res)= c('y','x')
  return(res)
}

vkm2latlong<-function(yx){
  y = yx[,1]
  x = yx[,2]
  latitude = y / 111.2
  longitude = x / 111.2 / cos(pi / 130 * latitude)
  res = cbind(latitude,longitude)
  colnames(res) = c('latitude', 'longitude')
  return(res)
}


Simulate<-function(population, time.factor, total.average){
  #
  #' @param population
  #' @param time.factor
  #' @param total.average
  baseline.matrix = population %o% time.factor 
  baseline.matrix = baseline.matrix / sum(baseline.matrix) * total.average 
  
  flatten.baseline.matrix = c(baseline.matrix)
  n.col = ncol(baseline.matrix)
  simulation = rpois(length(flatten.baseline.matrix), lambda=flatten.baseline.matrix)
  simulation = matrix(simulation, ncol=n.col)
  rownames(simulation) = as.character(names(population))
  colnames(simulation) = as.character(names(time.factor))
  return(simulation)
}



ranScanSimulate<-Simulate

