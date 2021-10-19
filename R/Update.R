
default.postcode.file = "Data/Postcode_Estimates_Table_with_coordinates2.csv"
require(readxl)
source("utils.R")

ranScanUpdate<-function(case.file, postcode.file=default.postcode.file){
  source("utils.R")
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){return(NULL)}
  )
  if(is.null(case.df)){
    writeLines("There is no `case.df` to update, run first `ranScanInit`")
    return(0)
  }
  
  # UPDATE case.df
  case.df.new = read_excel(case.file)
  case.df.new = as.data.frame(case.df.new)
  nomi = c("FULLNO", "Patient Postcode", "SAMPLE_DT", "RECEPT_DT", "Isolation Site Decoded", "Sterile Site Y N",
           "emm gene type")
  case.df.new = case.df.new[, nomi]
  
  idx = !duplicated(case.df$FULLNO, case.df.new$FULLNO, fromLast = T)
  
  case.df.new = case.df.new[idx, ]
  
  case.df.new['emmtypes2'] = apply(case.df.new['emm gene type'], 1, FUN=split)
  
  idx = which(case.df.new['emmtypes2']  == "c74a.0\r")
  case.df.new$emmtypes2[idx[1]] = 'c74a.0'
  case.df.new$emmtypes2[idx[2]] = 'c74a.0'
  
  idx = case.df.new$`Sterile Site Y N` == "#s"
  idx2 = !is.na(case.df.new$`Sterile Site Y N`)
  case.df.new = case.df.new[idx & idx2, ]
  
  # convert dates to integers
  case.df.new$RECEPT_DT_numeric = as.integer(difftime(case.df.new$RECEPT_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
  case.df.new$SAMPLE_DT_numeric = as.integer(difftime(case.df.new$SAMPLE_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
  
  # clean data  
  idx = (case.df.new$SAMPLE_DT_numeric > -1) & !is.na(case.df.new$SAMPLE_DT_numeric)
  case.df.new<-case.df.new[idx,]
  idx = is.na(case.df.new$emmtypes2)
  case.df.new$emmtypes2[idx] = 'NA'
  
  idx = is.na(case.df.new$`Patient Postcode`)
  case.df.new$`Patient Postcode`[idx] = 'NA'
  
  case.df.new$lag = case.df.new$SAMPLE_DT_numeric - case.df.new$RECEPT_DT_numeric

  
  # Insert coordinates
  if(!is.null(postcode.file)){
    postcode.data = read.csv(postcode.file, stringsAsFactors = F, header = T)
    case.df.new[,c("latitude", "longitude")] = t(apply(case.df.new, 1, postcode.to.location, postcode.data))
  }else{
    case.df.new[,c("latitude", "longitude")] = t(apply(case.df.new, 1, postcode.to.location2))      
  }

  case.df.updated = rbind(case.df, case.df.new)
#  save("case.df", file=paste0(case.file,".Rdata"))
  
  # UPDATE observation.matrix
  observation.matrix<-tryCatch({
    load(paste0(case.file, "obs.Rdata"))
    observation.matrix
  },
  error = function(e){
    return(NULL)
  }
  )
  if(is.null(observation.matrix)){
    writeLines("There is no `observation.matrix` to update, run first `ranScanCreateObservationMatrix`.")
    return(0)
  }
  
  postcodes = unique(case.df$`Patient Postcode`)
  emmtypes = unique(case.df$emmtypes2)
  n.postcodes= length(postcodes)
  n.emmtypes = length(emmtypes)
  n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
  
  new.postcodes = setdiff(unique(case.df.new$`Patient Postcode`) , postcodes)
  new.emmtypes = setdiff(unique(case.df.new$emmtypes2) , emmtypes)
  #
  new.n.postcodes= length(new.postcodes)
  new.n.emmtypes = length(new.emmtypes)
  new.n.weeks = max(case.df.updated$SAMPLE_DT_numeric[!is.na(case.df.updated$SAMPLE_DT_numeric)]) - min(case.df.updated$SAMPLE_DT_numeric[!is.na(case.df.updated$SAMPLE_DT_numeric)]) + 1 - n.weeks
  
  observation.matrix.new = array(data=0,
                             dim = c(n.postcodes + new.n.postcodes, n.emmtypes + new.n.emmtypes, n.weeks + 2 + new.n.weeks),
                             dimnames=list(c(postcodes, new.postcodes), c(emmtypes, new.emmtypes), c('NA', as.character(0:(n.weeks+new.n.weeks)))))
  
  # the dimensions correspond to postcode, emmtype, time  
  
  
  # copy the data from the old matrix
  observation.matrix.new[1:n.postcodes, 1:n.emmtypes, 1:(n.weeks+2)] = observation.matrix
  observation.matrix = observation.matrix.new
  
  for (i in 1:nrow(case.df.new)){
    emmtype = case.df[i,'emmtypes2']
    postcode = case.df[i,'Patient Postcode']
    week = case.df[i, 'SAMPLE_DT_numeric']
    if (is.na(week)){
      week = 'NA'
    }
    observation.matrix[postcode, emmtype, as.character(week)] = observation.matrix[postcode, emmtype, as.character(week)] + 1
  }
  save("observation.matrix", file=paste0(case.file, 'obs.Rata'))
  return(observation.matrix)

  
  
  # UPDADE THE baseline.matrix
  baseline.matrix<-tryCatch({
    load(paste0(case.file, "bas.Rdata"))
    baseline.matrix
  },
  error = function(e){
    return(NULL)
  }
  )
  if(is.null(baseline.matrix)){
    writeLines("There is no `baseline.matrix` to update, run first `ranScanCreateBaselineMatrix`.")
    return(0)
  }
  source("time_utils.R")
  
  
  
  time.factor = apply(observation.matrix, 3, sum)
  aggregated.process = expand.histogram(time.factor[-1]) # exclude the first entry (-1) with undated data
  Parameters =  cmle(aggregated.process, 50)
  x <- as.numeric(dimnames(observation.matrix)[[3]][-1])
  prediction.cmle = predict.cmle(x, Parameters)
  
  baseline.matrix = array(data=0,
                          dim = c(n.postcodes, n.emmtypes, n.weeks + 2),
                          dimnames=list(postcodes, emmtypes, c('NA', as.character(0:n.weeks)))
  )
  ###
  
  time.factor = c(time.factor[1], prediction.cmle)
  emmtype.factor = matrix(0, ncol = length(observation.matrix[1,,1]), nrow=length(time.factor))
  location.factor = matrix(0, ncol = length(observation.matrix[,1,1]), nrow=length(time.factor))
  
  for(i in 1:length(time.factor)){
    emmtype.factor[i,] = apply(observation.matrix[,,1:i], 2, sum)
    location.factor[i,] = apply(observation.matrix[,,1:i], 1, sum)
  }
  for(i in 1:dim(observation.matrix)[1]){
    for(j in 1:dim(observation.matrix)[2]){
      for(k in 1:dim(observation.matrix)[3]){ #time
        baseline.matrix[i,j,k] <- time.factor[k] * emmtype.factor[k, j] * location.factor[k, i]
      }
    }
  }
  
  baseline.matrix = baseline.matrix * sum(observation.matrix) / sum(baseline.matrix.3)
  save("baseline.matrix", file=paste0(case.file, "bas.Rdata"))  
  
  
}
  


