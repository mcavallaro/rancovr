




# 
# Init<-function(case.file){ #}, postcode.file=default.postcode.file){
#     source("utils.R")
#     case.df<-tryCatch({
#       ""
#       load(paste0(case.file, ".Rdata"))
#       case.df
#     },
#     error = function(e){
#       ""
#       case.df = read_excel(case.file)
#       case.df = as.data.frame(case.df)
#       nomi = c("FULLNO", "Patient Postcode", "SAMPLE_DT", "RECEPT_DT", "Isolation Site Decoded", "Sterile Site Y N",
#                "emm gene type")
#       case.df = case.df[, nomi]
#       
#       case.df['emmtype'] = apply(case.df['emm gene type'], 1, FUN=split)
#       case.df['emmtype'] = apply(case.df['emmtype'], 1, FUN=function(x){strsplit(x, '\r')[[1]][1]})
#       
#       idx = case.df$`Sterile Site Y N` == "#s"
#       idx2 = !is.na(case.df$`Sterile Site Y N`)
#       if (any(idx & idx2)){
#         case.df = case.df[idx & idx2, ]      
#       }
#       
#       # Insert coordinates
#       # if(!is.null(postcode.file)){
#       #   postcode.data = read.csv(postcode.file, stringsAsFactors = F, header = T)
#       #   case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location, postcode.data))
#       # }else{
#       #   case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location2))
#       # }
#       # convert dates to integers
#       case.df$RECEPT_DT_numeric = as.integer(difftime(case.df$RECEPT_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
#       case.df$SAMPLE_DT_numeric = as.integer(difftime(case.df$SAMPLE_DT, as.POSIXct("2015-01-01 UTC", tz = "UTC"), units = "weeks"))
#       
#       # clean data  
#       idx = (!is.na(case.df$SAMPLE_DT_numeric)) & (!is.na(case.df$RECEPT_DT_numeric))
#       if (any(idx)){
#         case.df = case.df[idx,]      
#       }
#       
#       idx = (case.df$SAMPLE_DT_numeric >= 0) & (case.df$RECEPT_DT_numeric >= 0)
#       if (any(idx)){
#         case.df = case.df[idx,]      
#       }
#       
#       idx = case.df$SAMPLE_DT_numeric <= case.df$RECEPT_DT_numeric
#       if (any(idx)){
#         case.df = case.df[idx,]
#       }
#       
#       idx = is.na(case.df$emmtype)
#       case.df$emmtype[idx] = 'NA'      
#       
#       idx = is.na(case.df$`Patient Postcode`)
#       case.df$`Patient Postcode`[idx] = 'NA'      
#       
#       case.df$lag = case.df$SAMPLE_DT_numeric - case.df$RECEPT_DT_numeric
#       
#       case.df$is.england = apply(case.df, 1, postcode.in.england)
#       
#       idx = is.na(case.df$is.england)
#       if (any(idx)){
#         case.df[idx,]$`Patient Postcode`= "NA"
#         case.df[idx,]$is.england = TRUE
#       }
#       
#       case.df = case.df[case.df$is.england, ]
#       case.df$is.england = NULL
#       
#       mm=mean(case.df$lag)
#       sdt=sd(case.df$lag)
#       
#       idx = abs(case.df$lag) < mm + 3 * sdt
#       if(any(idx)){
#         case.df = case.df[idx,]      
#       }
#       
#       attribute_list = attributes(case.df)
#       attribute_list$emmtypes = unique(case.df$emmtype)
#       attributes(case.df) <- attribute_list
#       
#       save.and.tell("case.df", file=paste0(case.file,".Rdata"))
#       return(case.df)
#     }
#     )
# }
# 

Clean<-function(pattern){
  list_files = list.files('Data',pattern = paste0(pattern, '.Rdata'))
  for (file in list_files){
    file.remove(file.path("Data", file))
  }
}

CreateObservationMatrices<-function(case.df, types=NA, date.time.field = 'week', postcode.field = 'postcode'){
  # emmtype is an vector of strings, e.g., 
  postcodes = unique(case.df[,postcode.field])
  n.postcodes= length(postcodes)
  if (is.na(types)){
    n.types = 1
  }else{
    n.types = length(emmtypes)
  }
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
  for (e in 1:n.types){
    if (!is.na(types)){
      type=types[e]
      cat("Step", e, "of", n.types, ", creating observation matrix for", type, "emm type\n")
      idx = case.df$type == type
      case.df = case.df[idx,]
    }
    observation.matrix = as(matrix(data=0,
                                   nrow=n.postcodes,
                                   ncol=n.weeks + 2,
                                   dimnames=list(postcodes,
                                                 c('NA', as.character(0:n.weeks)) )), "sparseMatrix")
    # the dimensions correspond to postcode, emmtype, time  
    for (i in 1:nrow(case.df)){
      postcode = case.df[i, postcode.field]
      week = case.df[i, date.time.field]
      if (is.na(week)){
        week = 'NA'
      }
      observation.matrix[postcode, as.character(week)] = observation.matrix[postcode, as.character(week)] + 1
    }
    attribute_list = attributes(observation.matrix)
    if (is.na(types)){
      attribute_list$type = NA
    }else{
      attribute_list$type = type
    }
    attribute_list$date.time.field = date.time.field
    attributes(observation.matrix) <- attribute_list
    
    save.and.tell("observation.matrix",
                  file=file.path(getwd(), paste0('observation_matrix_tmp.Rdata')))
  }
}

#' CreateObservationMatrices_<-function(case.file, emmtypes, starting.week, n.weeks){
#'   #' This function rearranges the observation of `case.file` into observation Matrices
#'   #' similarly to `CreateObservationMatrices` with the difference that a) it
#'   #' creates matrices for each week starting from `starting.week` until the last week in the dataset
#'   #'  and b) it counts the cases that are not typed at "SAMPLE_DT" but only are at "RECEIPT_DT" datetimes.
#'   #' It creates and saves in two `.Rdata` files:
#'   #'  1) a list of matrices that contain the observations of untyped case.
#'   #'  2) a list of list of matrices that contain the observations of the typed cases.
#'   #'  
#'   #' The inner list index indicates the week at which the observations are been taken,
#'   #' while the outer index in list 2) indicates an emmtype in `emmtypes`.
#'   #'  
#'   #' The matrices are spaved in sparse format requiring `library(Matrix)`.
#'   #'    
#'   #' In a practical situation, one sets `starting.week` to the current week and extract
#'   #' matrices the observations of the last `n.weeks`. The value of `n.weeks`
#'   #' cannot be larger than `starting.week`, otherwise one would need
#'   #' events from before the beginning to the study.
#'   #' 
#'   #' In testing phase, e.g., to evaluate timeliness one can set an earlier `starting.week`'.
#'   #'  
#'   #' @param case.file A string
#'   #' @param emmtypes An vector of emmtypes.
#'   #' @param starting.week integer. This is the week corresponding to the last columns of observation and baseline matrices.
#'   #' @param n.weeks no. columns of each matrix
#'   #' @examples
#'   #' CreateObservationMatrices_(case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx", c('1.0', '12.0'), , )
#'   # emmtype is an vector of strings, e.g., 
#'   case.df<-tryCatch({
#'     load(paste0(case.file, ".Rdata"))
#'     case.df
#'   },
#'   error = function(e){
#'     case.df = Init(case.file)
#'     return(case.df)
#'   }
#'   )
#' 
#'   postcodes = unique(case.df$`Patient Postcode`)
#'   # emmtypes = unique(case.df$emmtype)
#'   n.postcodes= length(postcodes)
#'   n.emmtypes = length(emmtypes)
#'   MAX = max(c(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)], case.df$RECEPT_DT_numeric[!is.na(case.df$RECEPT_DT_numeric)]))
#'   MIN = min(c(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)], case.df$RECEPT_DT_numeric[!is.na(case.df$RECEPT_DT_numeric)]))
#'   cat("Max allowed week range is", MIN, MAX, "\n")
#'   if (starting.week < n.weeks){
#'     starting.week = n.weeks
#'     cat(sprintf("We enforced `starting.week=n.week=%d`", n.weeks), ".\n")
#'   }
#'   cat("We are creating", MAX - starting.week + 1, "matrices per emmtype, each with", n.weeks, "columns.\n")
#'   
#'   observation.matrices.untyped = list()
#'   for (week in starting.week:MAX){
#'     # this represents the week at week the scores are calculated
#'     dimnames = list(postcodes, c(as.character( (week - n.weeks+1):(week)) ))
#'     observation.matrices.untyped[[as.character(week)]] = as(matrix(data=0,
#'                                                                    nrow=n.postcodes,
#'                                                                    ncol=n.weeks, 
#'                                                                    dimnames=dimnames), "sparseMatrix")
#'   }
#'   observation.matrices.typed = list()
#'   for (emmtype in emmtypes){
#'     observation.matrices.typed[[emmtype]] = list()
#'     for (week in starting.week:MAX){
#'       dimnames = list(postcodes, c(as.character( (week-n.weeks+1):(week) ) ))
#'       observation.matrices.typed[[emmtype]][[as.character(week)]] = as(matrix(data=0,
#'                                                                               nrow=n.postcodes,
#'                                                                               ncol=n.weeks,
#'                                                                               dimnames=dimnames), "sparseMatrix")
#'     }
#'   }
#'   
#'   for (week in starting.week:MAX){
#'     # this represents the week at week the scores are calculated
#'     # fill the untyped observation.matrix with the cases that were detected by week "week"
#'     # (SAMPLE_DT<week) but whose typing occurred after the same "week" (RECEIPT_DT>week)
#'     idx = (case.df$SAMPLE_DT_numeric > week - n.weeks + 1) &
#'       (case.df$SAMPLE_DT_numeric <= week) &
#'       (case.df$RECEPT_DT_numeric > week)
#'     if (any(idx)){
#'       case.df.idx = case.df[idx,]
#'       for (i in 1:nrow(case.df.idx)){
#'         postcode = case.df.idx[i, 'Patient Postcode']
#'         week_ = case.df.idx[i, 'SAMPLE_DT_numeric']
#'         observation.matrices.untyped[[as.character(week)]][postcode, as.character(week_)] = 
#'           observation.matrices.untyped[[as.character(week)]][postcode, as.character(week_)] + 1
#'       }
#'     }
#'     
#'     for (emmtype in emmtypes){
#'       # fill the typed observation.matrix with the cases that were detected and typed by week "week"
#'       # (RECEIPT_DT<week)
#'       idx = (case.df$emmtype == emmtype) &
#'         (case.df$SAMPLE_DT_numeric > week - n.weeks + 1) &
#'         (case.df$RECEPT_DT_numeric <= week) &
#'         (case.df$SAMPLE_DT_numeric <= week)
#'       if (any(idx)){
#'         case.df.idx = case.df[idx,]
#'         for (i in 1:nrow(case.df.idx)){
#'           postcode = case.df.idx[i, 'Patient Postcode']
#'           week_ = case.df.idx[i, 'SAMPLE_DT_numeric']
#'           if (is.na(week_)){
#'             week_ = 'NA'
#'           }
#'           observation.matrices.typed[[emmtype]][[as.character(week)]][postcode, as.character(week_)] = 
#'             observation.matrices.typed[[emmtype]][[as.character(week)]][postcode, as.character(week_)] + 1
#'         }
#'       }
#'     }
#'   }
#'   attribute_list = list(names = names(observation.matrices.untyped), starting.week=starting.week, n.weeks=n.weeks, week.max=MAX)
#'   attributes(observation.matrices.untyped) <- attribute_list
#'   attribute_list = list(names = names(observation.matrices.typed),starting.week=starting.week, n.weeks=n.weeks, week.max=MAX)
#'   attributes(observation.matrices.typed) <- attribute_list
#' 
#'   save.and.tell("observation.matrices.untyped",
#'                 file=file.path(dirname(case.file), paste0('untyped', '_obs.Rdata')))
#'   save.and.tell("observation.matrices.typed",
#'                 file=file.path(dirname(case.file),paste0(c(emmtypes, '_typed_obs.Rdata'), collapse='')))
#' }


#' CreateObservationMatrices.prospective<-function(case.file, emmtypes, n.weeks){
#'   #' This is a simplified version of `CreateObservationMatrices.delay` that
#'   #' create a sparse observation matrix only for the last (current) week.
#'   #' This function rearranges the observation of `case.file` into observation Matrices
#'   #' similarly to `CreateObservationMatrices` with the difference that it counts the cases that are not typed at "SAMPLE_DT" but only are at "RECEIPT_DT" datetimes.
#'   #' It creates and saves in two `.Rdata` files:
#'   #'  1) a matrix that contains the observations of untyped case.
#'   #'  2) a list of matrices that contain the observations of the typed cases (one, for each emm type).
#'   #'  
#'   #' The matrices are spaved in sparse format requiring `library(Matrix)`.
#'   #'  
#'   #' @param case.file A string
#'   #' @param emmtypes An vector of emmtypes.
#'   #' @param n.weeks no. columns
#'   #' @examples
#'   #' CreateObservationMatrices.prospective(case.file = "Data/Full MOLIS dataset minus PII 20200918.xlsx", c('1.0', '12.0'), 100)
#'   library(Matrix)
#'   # emmtype is an vector of strings, e.g., 
#'   case.df<-tryCatch({
#'     load(paste0(case.file, ".Rdata"))
#'     case.df
#'   },
#'   error = function(e){
#'     case.df = Init(case.file)
#'     return(case.df)
#'   }
#'   )
#'   postcodes = unique(case.df$`Patient Postcode`)
#'   n.postcodes= length(postcodes)
#'   n.emmtypes = length(emmtypes)
#'   MAX = max(c(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)], case.df$RECEPT_DT_numeric[!is.na(case.df$RECEPT_DT_numeric)]))
#' 
#'   dimnames = list(postcodes, c(as.character( (week - n.weeks+1):(week) ) ))
#'   observation.matrix.untyped = as(matrix(data=0,
#'                                          nrow=n.postcodes,
#'                                          ncol=n.weeks, 
#'                                          dimnames=dimnames), "sparseMatrix")
#'   ### continua da qui
#'   observation.matrices.typed = list()
#'   for (emmtype in emmtypes){
#'     dimnames = list(postcodes, c(as.character( (week-n.weeks+1):(week) ) ))
#'     observation.matrices.typed[[emmtype]] = as(matrix(data=0,
#'                                                       nrow=n.postcodes,
#'                                                       ncol=n.weeks,
#'                                                       dimnames=dimnames), "sparseMatrix")
#'   }
#'   week = ???
#'   # this is the week at week the scores are calculated
#'   # fill the untyped observation.matrix with the cases that were detected by week "week"
#'   # (SAMPLE_DT<week) but whose typing occurred after the same "week" (RECEIPT_DT>weel)
#'   idx = (case.df$SAMPLE_DT_numeric > week - n.weeks + 1) &
#'       (case.df$SAMPLE_DT_numeric <= week) &
#'       (case.df$RECEPT_DT_numeric > week)
#'   if (any(idx)){
#'       case.df.idx = case.df[idx,]
#'       for (i in 1:nrow(case.df.idx)){
#'         postcode = case.df.idx[i, 'Patient Postcode']
#'         week_ = case.df.idx[i, 'SAMPLE_DT_numeric']
#'         observation.matrix.untyped[postcode, as.character(week_)] = 
#'           observation.matrices.untyped[postcode, as.character(week_)] + 1
#'       }
#'   }
#'     
#'   for (emmtype in emmtypes){
#'       # fill the typed observation.matrix with the cases that were detected and typed by week "week"
#'       # (RECEIPT_DT<week)
#'       idx = (case.df$emmtype == emmtype) &
#'         (case.df$SAMPLE_DT_numeric > week - n.weeks + 1) &
#'         (case.df$RECEPT_DT_numeric <= week) &
#'         (case.df$SAMPLE_DT_numeric <= week)
#'       if (any(idx)){
#'         case.df.idx = case.df[idx,]
#'         for (i in 1:nrow(case.df.idx)){
#'           postcode = case.df.idx[i, 'Patient Postcode']
#'           week_ = case.df.idx[i, 'SAMPLE_DT_numeric']
#'           if (is.na(week_)){
#'             week_ = 'NA'
#'           }
#'           observation.matrices.typed[[emmtype]][postcode, as.character(week_)] = 
#'             observation.matrices.typed[[emmtype]][postcode, as.character(week_)] + 1
#'         }
#'       }
#'   }
#'   attribute_list = list(n.weeks = n.weeks, week.max = MAX)
#'   attributes(observation.matrix.untyped) <- attribute_list
#'   attribute_list = list(n.weeks = n.weeks, week.max = MAX)
#'   attributes(observation.matrices.typed) <- attribute_list
#'   
#'   save.and.tell("observation.matrix.untyped",
#'                 file=file.path(dirname(case.file), paste0('untyped_prospective', '_obs.Rdata')))
#'   save.and.tell("observation.matrices.typed",
#'                 file=file.path(dirname(case.file),paste0(c(emmtypes, '_typed_prospective_obs.Rdata'), collapse='')))
#' }
#' 
#' 

PostcodeMap<-function(matrix){
  # Create a dataframe that maps postcodes to coordinates
  # `matrix` can be either an baseline.matrix or an baseline.matrix
  postcode.field = 'postcode'
  writeLines("Compiling the table that maps the rows of the observation/baseline matrix to geo-coordinates and population.")
  ret<-tryCatch({
    load("postcode2coord.Rdata", verbose = 1)
    if (all(as.character(postcode2coord[, postcode.field]) == rownames(observation.matrix))){
      writeLines("Using data loaded from `postcode2coord.Rdata`")
      postcode2coord
    }else{
      writeLines("Data loaded from `postcode2coord.Rdata` is for a different matrix and will be overwritten by the map for the current matrix.")
      stop() #raise error
    }
  },
  error = function(e){
    postcode2coord = data.frame(postcode.field = rownames(matrix))
    names(postcode2coord) = postcode.field
    # Insert coordinates and population density
    
    # postcode.data = 
    
    postcode2coord[,c("latitude", "longitude", "Total")] = t(apply(postcode2coord, 1,
                                                                   postcode.to.location.and.population,
                                                                   postcode.data))
    save.and.tell('postcode2coord', file=file.path(getwd(),
                                                   paste0("postcode2coord_tmp.Rdata")))

    return(postcode2coord)
  }
  )
  return(ret)
}

TimeFactor<-function(case.df, save.on.dir = TRUE, get.from.dir = FALSE, date.time.field = "week", parameters = NULL){
  time.factor<-tryCatch({
    if(get.from.dir == FALSE){
      stop("Not loading from local directory.")
    }
    load(file.path(getwd(), paste0('timefactor_tmp.Rdata')))
    cat("Temporal baseline loaded.\n")
    time.factor
  },
  error = function(e){
    cat("Computing the temporal baseline.\n")
    Parameters = cmle(case.df[,date.time.field], 10, parameters)
    n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
    x = 0:n.weeks
    prediction.cmle = predict.cmle(x, Parameters)
    
    na = sum(is.na(case.df[,date.time.field]))
    
    time.factor = c(na, prediction.cmle)
    names(time.factor) = c('NA', x)
    attr=attributes(time.factor)
    attr$Parameters = Parameters
    attributes(time.factor) = attr
    if(save.on.dir){
      save.and.tell('time.factor', file=file.path(getwd(), paste0('timefactor_tmp.Rdata')))
    }
    return(time.factor)
  })
  return(time.factor)
}


EmmtypeFactor<-function(case.file){
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


EmmtypeFactor.tau<-function(case.file, emmtypes, date.time.field = 'SAMPLE_DT_numeric'){
  load(paste0(case.file, ".Rdata"))
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1

  xy.list=list()
  
  for (emmtype in emmtypes) {
    tmp = vector(mode = 'numeric', length=n.weeks)
    for (tau in 1:n.weeks){
      idx = (case.df$SAMPLE_DT_numeric <= tau) & (case.df$emmtype == emmtype)
      idx.tau = case.df$SAMPLE_DT_numeric <= tau
      tmp[tau] = sum(idx) / sum(idx.tau)
    }
    names(tmp) = as.character(1:n.weeks)
    xy.list[[emmtype]]=tmp
  }

  return(xy.list)
}


#' At a given week, a fraction lambda_untyped \approx 0.6 of all cases are not typed.
#' the baselines e.g. are as follows:
#' - for the emmtype 33.0: (1-lambda_untyped) * lambda_33.0 * lambda_t * lambda_geo
#' - for the entyped: lambda_untyped * lambda_t * lambda_geo
#' 
EmmtypeFactor.delay_<-function(case.file, starting.week, n.weeks){
  if (starting.week < n.weeks){
    starting.week = n.weeks
    cat(sprintf("We enforced `starting.week=n.week=%d`", n.weeks), ".\n")
  }
  load(paste0(case.file, ".Rdata"))
  #n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
  MAX = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)])
  # MIN = min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)])
  # n.weeks=MAX-MIN+1
  # tmp = vector(mode = 'numeric', length=n.weeks)
  tmp = rep(0, n.weeks)
  # len = NROW(case.df)
  # case.df = case.df[case.df$SAMPLE_DT_numeric <= case.df$RECEPT_DT_numeric,]
  for(w in starting.week:MAX){
    for (i in 1:n.weeks){
      week = w - n.weeks + i
      idx1 = case.df$SAMPLE_DT_numeric == week
      idx2 = case.df$RECEPT_DT_numeric > w
      tmp[i] = tmp[i] + sum(idx1 & idx2) / sum(idx1)
    }
  }
  tmp = tmp / length(starting.week:MAX) # (MAX - starting.week + 1)
  # cum.tmp = cumsum(tmp) / 1:length(tmp)
  # for (i in 1:n.weeks){
  #   if (is.na(cum.tmp[i])){
  #     cum.tmp[i] = cum.tmp[i-1]
  #   }
  # }
  # names(cum.tmp) = as.character(1:n.weeks)
  return(tmp)
}


CreateBaselineMatrix<-function(case.df, save.on.dir=FALSE,
                               date.time.field='week', postcode.field='postcode'){
  
  postcodes = unique(case.df[,postcode.field])
  n.postcodes= length(postcodes)
  
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
  baseline.matrix = matrix(data=0,
                              nrow = n.postcodes,
                              ncol = n.weeks + 2,
                              dimnames=list(postcodes, c('NA', as.character(0:n.weeks))))

  time.factor = TimeFactor(case.df, save.on.dir, date.time.field)
  spatial.factor = PostcodeMap(matrix(data=0,
                                      nrow = n.postcodes,
                                      ncol = n.weeks + 2,
                                      dimnames=list(postcodes, c('NA', as.character(0:n.weeks)) )))$Total
  
  for(i in 1:nrow(baseline.matrix)){
    for(j in 1:ncol(baseline.matrix)){
        baseline.matrix[i,j] = time.factor[j] * spatial.factor[i]
    }
  }
  baseline.matrix = baseline.matrix / sum(spatial.factor, na.rm = T)
  
  attribute_list = attributes(baseline.matrix)
  attribute_list$date.time.field = date.time.field
  attributes(baseline.matrix) <- attribute_list
  
  if(save.on.dir){
    save.and.tell("baseline.matrix", file=file.path(getwd(), paste0('baseline_matrix_tmp.Rdata')))
  }
  return(baseline.matrix)
}
