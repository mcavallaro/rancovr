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
#' 
#' 
#' 


# warning_ratio<-function(i, observation.matrix, cylinders, postcode.locations){
#   # check if the location i
#   x = as.numeric(postcode.locations[i, 'longitude'])
#   y = as.numeric(postcode.locations[i, 'latitude'])
#   
#   in_circle = sqrt((as.numeric(cylinders['x']) - x)^2 + (as.numeric(cylinders['y']) - y)^2) < as.numeric(cylinders['rho'])
#   
#   times = which(observation.matrix[i,] > 0)
# 
#   if (length(times) > 0){
#     in_cylinder_height = matrix(FALSE, nrow = nrow(cylinders), ncol=length(times))
# 
#     for (i in 1:length(times)){
#       
#       #da vettorizzare:
#       in_cylinder_height[,i] = apply(cylinders, 1,  function(X){
#         TT = as.numeric(names(times[i]))
#         (as.numeric(X['t.low']) < TT) & (as.numeric(X['t.upp']) > TT)
#       })
#       
#       
#       # vettorizzato
#       #TT = as.numeric(names(times[i]))
#       #in_cylinder_height[,i] = (as.numeric(cylinders['t.low']) < TT) & (as.numeric(cylinders['t.upp']) > TT)
# 
#     }
# 
#     # number of cylinders that include locations `i`
#     in_cylinder = sum(in_circle * in_cylinder_height)
#     # number of cylinder with `warning` flag that include location `i`
#     warning = sum(cylinders$warning * in_circle * in_cylinder_height)
# 
#     return(warning / in_cylinder)
#   }else{
#     return(NA)
#   }
# }
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
#' 
#' 
#' 
#' 