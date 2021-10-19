
postcode.file = "Data/Postcode_Estimates_Table_with_coordinates.csv"
require(readxl)
source("utils.R")
case.file="Data/Full MOLIS dataset minus PII 20200918.xlsx"

ranScanInit<-function(case.file){ #}, postcode.file=default.postcode.file){
    case.df = read_excel(case.file)
    case.df = as.data.frame(case.df)
    nomi = c("FULLNO", "Patient Postcode", "SAMPLE_DT", "RECEPT_DT", "Isolation Site Decoded", "Sterile Site Y N",
             "emm gene type")
    case.df = case.df[, nomi]
    
    case.df['emmtype'] = apply(case.df['emm gene type'], 1, FUN=split)
    case.df['emmtype'] = apply(case.df['emmtype'], 1, FUN=function(x){strsplit(x, '\r')[[1]][1]})

    idx = case.df$`Sterile Site Y N` == "#s"
    idx2 = !is.na(case.df$`Sterile Site Y N`)
    case.df = case.df[idx & idx2, ]
    
    # Insert coordinates
    # if(!is.null(postcode.file)){
    #   postcode.data = read.csv(postcode.file, stringsAsFactors = F, header = T)
    #   case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location, postcode.data))
    # }else{
    #   case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location2))
    # }

    # convert dates to integers
    case.df$RECEPT_DT_numeric = as.integer(difftime(case.df$RECEPT_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
    case.df$SAMPLE_DT_numeric = as.integer(difftime(case.df$SAMPLE_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
  
    # clean data  
    idx = (case.df$SAMPLE_DT_numeric >= 0) & (!is.na(case.df$SAMPLE_DT_numeric)) & (case.df$RECEPT_DT_numeric >= 0)
    case.df<-case.df[idx,]
    
    idx = is.na(case.df$emmtype)
    case.df$emmtype[idx] = 'NA'
    
    idx = is.na(case.df$`Patient Postcode`)
    case.df$`Patient Postcode`[idx] = 'NA'
    
    case.df$lag = case.df$SAMPLE_DT_numeric - case.df$RECEPT_DT_numeric
    
    case.df$is.england = apply(case.df, 1, postcode.in.england)
    
    idx = is.na(case.df$is.england)
    case.df[idx,]$`Patient Postcode`= "NA"
    
    case.df[idx,]$is.england = TRUE
    
    case.df = case.df[case.df$is.england, ]
    
    case.df$is.england =NULL
    
    mm=mean(case.df$lag)
    sdt=sd(case.df$lag)
    idx = abs(case.df$lag) < mm + 3 * sdt
    case.df = case.df[idx,]
    save("case.df", file=paste0(case.file,".Rdata"))
    return(list(case.df=case.df, emmtypes=unique(case.df$emmtype)))
}

ranScanClean<-function(pattern){
  list_files = list.files('Data',pattern = paste0(pattern, '.Rdata'))
  for (file in list_files){
    file.remove(file.path("Data", file))
  }
}

ranScanCreateObservationMatrices<-function(case.file, emmtypes){
  # emmtype is an array of strings, e.g., 
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
    },
    error = function(e){
      case.df = ranScanInit(case.file)
      return(case.df)
    }
  )
  postcodes = unique(case.df$`Patient Postcode`)
  # emmtypes = unique(case.df$emmtype)
  n.postcodes= length(postcodes)
  n.emmtypes = length(emmtypes)
  n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
  for (e in 1:n.emmtypes){
    emmtype=emmtypes[e]
    cat("Step ",e, " of ", n.emmtypes, ", creating matrix for", emmtype, "emm type\n")
    idx = case.df$emmtype == emmtypes[e]
    case.df.idx = case.df[idx,]
    observation.matrix = array(data=0,
                             dim = c(n.postcodes, n.weeks + 2),
                             dimnames=list(postcodes, c('NA', as.character(0:n.weeks)) ))
    # the dimensions correspond to postcode, emmtype, time  
    for (i in 1:nrow(case.df.idx)){
#      emmtype = case.df[i,'emmtype']
      postcode = case.df.idx[i,'Patient Postcode']
      week = case.df.idx[i, 'SAMPLE_DT_numeric']
      if (is.na(week)){
        week = 'NA'
      }
      observation.matrix[postcode, as.character(week)] = observation.matrix[postcode, as.character(week)] + 1
    }
    save("observation.matrix", file=file.path(dirname(case.file), paste0(emmtype, '_obs.Rdata')))
  }
}

ranScanPostcodeMap<-function(observation.matrix, postcode.file.name=postcode.file){
  # Create a dataframe that maps postcodes to coordinates
  cat("Compiling the table that maps the postcodes to geo-coordinates and population.\n")
  source("utils.R")
  
  ret<-tryCatch({
    load("Data/postcode2coord.Rdata")
    postcode2coord
  },
  error = function(e){
    postcode2coord = data.frame("Patient Postcode" = rownames(observation.matrix))
    names(postcode2coord) = "Patient Postcode"
    # Insert coordinates and population density
    postcode.data = read.table(postcode.file.name, header=T, stringsAsFactors = F, sep=',')
    postcode.data['Patient.Postcode'] = c(sapply(postcode.data['Patient.Postcode'],
                                                 function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))
    postcode2coord[,c("latitude", "longitude", "Total")] = t(apply(postcode2coord, 1,
                                                                   postcode.to.location.and.population, postcode.data))
    save('postcode2coord', file="Data/postcode2coord.Rdata")
    return(postcode2coord)
  }
  )
  return(ret)
}

ranScanTimeFactor<-function(case.file, parameters=NULL){
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){
    case.df = ranScanInit(case.file)
    return(case.df)
  }
  )
  time.factor<-tryCatch({
    load(paste0(case.file, "_timefactor.Rdata"))
    cat("Temporal baseline loaded.\n")
    time.factor
  },
  error = function(e){
    cat("Computing the temporal baseline.\n")
    source("time_utils.R")
    Parameters = cmle(case.df$RECEPT_DT_numeric, 50, parameters)
    n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
    x = 0:n.weeks
    prediction.cmle = predict.cmle(x, Parameters)
    
    na = sum(is.na(case.df$RECEPT_DT_numeric))
    
    time.factor = c(na, prediction.cmle)
    names(time.factor) = c('NA', x)
    save('time.factor', file=paste0(case.file, "_timefactor.Rdata"))
    
    return(time.factor)
  }
  )
  return(time.factor)
}




ranScanEmmtypeFactor<-function(case.file){
  load(paste0(case.file, ".Rdata"))
  emmtype.factor = c(table(case.df$emmtype))
  emmtype.factor = emmtype.factor/ sum(emmtype.factor)
  xy.list=list()
  for (i in 1:length(emmtype.factor)) {
    xy.list[[i]]=emmtype.factor[i]
  }
  xy.list = setNames(xy.list, names(emmtype.factor))
  return(xy.list)
}

ranScanEmmtypeFactor.tau<-function(case.file, emmtypes){
  load(paste0(case.file, ".Rdata"))
  n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1

  xy.list=list()
  
  for (emmtype in emmtypes) {
    tmp = vector(mode = 'numeric', length=n.weeks)
    for (tau in 1:n.weeks){
      idx = (case.df$SAMPLE_DT_numeric <= tau) & (case.df$emmtype == emmtype)
      idx.tau = case.df$SAMPLE_DT_numeric <= tau
      tmp[tau] = sum(idx) / sum(idx.tau)
    }
    xy.list[[emmtype]]=tmp
  }
  return(xy.list)
}


ranScanCreateBaselineMatrix<-function(case.file){
  
  # emmtype is an array of strings, e.g., 
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){
    case.df = ranScanInit(case.file)
    return(case.df)
  }
  )
  
  postcodes = unique(case.df$`Patient Postcode`)
  n.postcodes= length(postcodes)
  n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
  
  baseline.matrix = array(data=0,
                          dim = c(n.postcodes, n.weeks + 2),
                          dimnames=list(postcodes, c('NA', as.character(0:n.weeks))))
  
  time.factor = ranScanTimeFactor(case.file)
  spatial.factor = ranScanPostcodeMap(array(data=0,
                                         dim = c(n.postcodes, n.weeks + 2),
                                         dimnames=list(postcodes, c('NA', as.character(0:n.weeks)) )))$Total
  for(i in 1:nrow(observation.matrix)){
    for(j in 1:ncol(observation.matrix)){
        baseline.matrix[i,j] = time.factor[j] * spatial.factor[i]
    }
  }
  baseline.matrix = baseline.matrix / sum(spatial.factor, na.rm = T)
  save("baseline.matrix", file=paste0(case.file, "_bas.Rdata"))
  return(baseline.matrix)
}


