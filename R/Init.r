
#' Clean
#'
#' Delete all '*.RData' files with given basename.
#' Can be used to clean the saved observation and baseline matrices.
#' Use with caution.
#' 
#' @param basename A \code{character} string.
#' @return None
#' @examples
#' Clean("33.0_observation_matrix_tmp")
#' Clean("observation_matrix_tmp")
Clean<-function(basename){
  list_files = list.files('Data',pattern = paste0(basename, '.RData'))
  for (file in list_files){
    file.remove(file.path("Data", file))
  }
}

#' CreateObservationMatrices
#'
#' Compute 2D Sparse matrices enconding observations recorded in \code{case.df}. 
#' Rows are locations (indexed by postcodes) and columns are time steps.
#' If a \code{character} vector of length N is passed as argument \code{types}, then N
#' matrices are created only for the events with matching \code{case.df$types}.
#' The matrices are saved on disk as '*observation_matrix.RData' files.
#' 
#' @param case.df A \code{data.frame} of events.
#' @param types \code{NULL} or a \code{character} vector.
#' @param date.time.field A \code{character} string.
#' @param postcode.field A \code{character} string.
#' @param more.postcodes A \code{character} vector
#' @param more.weeks A \code{integer} vector
#' @return None
#' @importFrom Matrix sparseMatrix
#' @examples
#' CreateObservationMatrices(case.df)
#' CreateObservationMatrices(case.df, types=c("1.0", "33.0"), date.time.field = "SAMPLE_DT_numeric", postcode.field = "Patient Postcode")
CreateObservationMatrices<-function(case.df, types=NULL, date.time.field = 'week', postcode.field = 'postcode',
  more.postcodes=c(), more.weeks=c()){

  if (length(more.postcodes) > 0){
    postcodes = sort(unique(c(unlist(case.df[,postcode.field]), unlist(more.postcodes))))
    postcodes = sort(postcodes)
  }else{
    postcodes = sort(unique(case.df[,postcode.field]))
  }
  n.postcodes = length(postcodes)

  if (is.null(types)){
    n.types = 1
  }else{
    n.types = length(types)
  }

  maxim = max(c(case.df[,date.time.field], unlist(as.integer(more.weeks))), na.rm = T)
  minim = min(c(case.df[,date.time.field], unlist(as.integer(more.weeks))), na.rm = T)

  n.weeks = maxim - minim + 1
  for (e in 1:n.types){
    if (!is.null(types)){
      type=types[e]
      cat("Step", e, "of", n.types, ", creating observation matrix for", type, "type\n")
      idx = case.df$type == type
      case.df = case.df[idx,] ## BUG here do not overwrite case.df
    }
    observation.matrix = as(matrix(data=0,
                                   nrow=n.postcodes,
                                   ncol=n.weeks + 1,
                                   dimnames=list(postcodes,
                                                 c('NA', as.character(minim:maxim)) )), "sparseMatrix")
    # the dimensions correspond to postcode, type, time  
    for (i in 1:nrow(case.df)){
      postcode = case.df[i, postcode.field]
      week = case.df[i, date.time.field]
      if (is.na(week)){
        week = 'NA'
      }
      observation.matrix[postcode, as.character(week)] = observation.matrix[postcode, as.character(week)] + 1
    }
    attribute_list = attributes(observation.matrix)
    if (is.null(types)){
      attribute_list$type = NA
    }else{
      attribute_list$type = type
    }
    attribute_list$date.time.field = date.time.field
    attributes(observation.matrix)<-attribute_list
    
    if (!is.null(types)){
      save.and.tell("observation.matrix",
                  file = file.path(getwd(), paste0(type, '_observation_matrix.RData')))
    }else{
      save.and.tell("observation.matrix",
                  file = file.path(getwd(), 'observation_matrix.RData'))
    }
  }
}


#' Map postcodes to coordinates
#
#' Returns and save 'postcode2coord.RData' a data frame that maps the postcodes included in \code{rownames(matrix)} to geographical coordinates.
#' \code{matrix} can be a \code{Matrix} or a \code{sparseMatrix} object storing baseline or observation data.
#' Requires the all postcode data tabulated in a \code{data.frame} called \code{postcode.data} available in the workspace.
#' 
#' @param matrix A \code{Matrix} or \code{sparseMatrix}.
#' @return A \code{data.frame}.
#' @examples
#' postcode2coord = PostcodeMap(observation.matrix)
PostcodeMap<-function(matrix, postcode.field = 'postcode'){
  writeLines("Compiling the table that maps the rows of the observation/baseline matrix to geo-coordinates and population.")
  ret<-tryCatch({
    load("postcode2coord.RData", verbose = 1)
    if (all(as.character(postcode2coord[, postcode.field]) == rownames(matrix))){
      writeLines("Using data loaded from `postcode2coord.RData`")
      postcode2coord
    }else{
      writeLines("Data loaded from `postcode2coord.RData` is for a different matrix and will be overwritten by the map for the current matrix.")
      stop() #raise error
    }
  },
  error = function(e){
    postcode2coord = data.frame(postcode.field = rownames(matrix))
    names(postcode2coord) = postcode.field
    
    postcode2coord$index = 1:NROW(matrix)
    postcode2coord['key'] = c(sapply(postcode2coord[postcode.field],
                                     function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))
    
    # Insert coordinates and population density
    # postcode.data is a data.frame to load.
    data("UK_population_per_postcode_with_coordinates")
    postcode.data['key'] = c(sapply(postcode.data['postcode'],
                                    function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))
    
    postcode.data[postcode.field] = NULL
    
    postcode2coord = merge(postcode2coord, postcode.data, by='key', all.x=T, all.y=F, sort=F)
    postcode2coord = postcode2coord[order(postcode2coord$index),]
    
    idx = is.na(postcode2coord$latitude)
    if(any(idx)){
      cat("Retrieving coordinates from api.getthedata for", sum(idx), "postcodes.\n")
      TTT=Sys.time()
      postcode2coord[idx, c("latitude", "longitude", "Total")] = t(apply(postcode2coord[idx,], 1,
                                                                         postcode.to.location.and.population,
                                                                         postcode.data, postcode.field = 'key'))
      print(Sys.time() - TTT)      
    }    
    
    postcode2coord['key'] = NULL
    save.and.tell('postcode2coord', file=file.path(getwd(),
                                                   paste0("postcode2coord.RData")))

    return(postcode2coord)
  }
  )
  return(ret)
}


#' Find temporal component of the baseline
#'
#' Returns a vector representing the temporal component of the baseline. This is obtained by calling
#' \code{cmle} and \code{predict.cmle}, which fit a seasonal trend model to aggregated data and return its best estimate, respectively. 
#' 
#' @param case.df A \code{data.frame} containing the events.
#' @param save.on.dir \code{logical}. If TRUE then the vector is saved in `timefactor.RData` file.
#' @param get.from.dir \code{logical}. If TRUE then the vector is obtained from the `timefactor.RData` file.
#' @param date.time.field A \code{character} string.
#' @param start \code{numeric}. Starting parameters (not yet implemented).
#' @param n.iterations An \code{integer}.
#' @return A \code{numeric} vector.
#' @examples
#' time.factor = TimeFactor(case.df)
TimeFactor<-function(case.df, save.on.dir = TRUE, get.from.dir = FALSE,
                     date.time.field = "week", start = NULL, n.iterations=20){
  time.factor<-tryCatch({
    if(get.from.dir == FALSE){
      stop("Not loading from local directory.")
    }
    load(file.path(getwd(), paste0('timefactor.RData')))
    cat("Temporal baseline loaded.\n")
    time.factor
  },
  error = function(e){
    cat("Computing the temporal baseline.\n")
    Parameters = cmle(case.df[,date.time.field], n.iterations, start)
    maxim = max(case.df[,date.time.field], na.rm = T)
    minim = min(case.df[,date.time.field], na.rm = T)
    n.weeks = maxim - minim + 1
    x = minim:maxim
    
    prediction.cmle = predict.cmle(x, Parameters)
    
    na = sum(is.na(case.df[,date.time.field]))
    
    time.factor = c(na, prediction.cmle)
    names(time.factor) = c('NA', x)
    attr=attributes(time.factor)
    attr$Parameters = Parameters
    attributes(time.factor) = attr
    if(save.on.dir){
      save.and.tell('time.factor', file=file.path(getwd(), paste0('timefactor.RData')))
    }
    return(time.factor)
  })
  return(time.factor)
}


EmmtypeFactor<-function(case.file){
  load(paste0(case.file, ".RData"))
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
  load(paste0(case.file, ".RData"))
  n.weeks = max(case.df[,date.time.field], na.rm = T) - min(case.df[,date.time.field], na.rm = T) + 1

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
  load(paste0(case.file, ".RData"))
  #n.weeks = max(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) - min(case.df$SAMPLE_DT_numeric[!is.na(case.df$SAMPLE_DT_numeric)]) + 1
  MAX = max(case.df$SAMPLE_DT_numeric, na.rm = T)
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

#' Create and save a baseline matrix
#'
#' Compute and save on disk a 2D dense matrix representing the
#' baseline, given a \code{data.frame} of events (\code{case.df}), where  rows represent the spatial component of the baseline
#' (indexed, e.g., by postcodes) and columns represent the temporal component.
#' It is possible to insert null rows corresponding to postcodes not included in \code{case.df} using \code{more.postcodes} field
#' for increased accuracy.
#' 
#' @param case.df A \code{data.frame} of events.
#' @param save.on.dir \code{logical}. If TRUE then the vector is saved in `baseline_matrix.RData` file.
#' @param date.time.field A \code{character} string.
#' @param postcode.field A \code{character} string.
#' @param more.postcodes A \code{character} vector of postcodes.
#' @param more.weeks A \code{integer} vector of postcodes.
#' @return A \code{Matrix}.
#' @examples
#' baseline.matrix = CreateBaselineMatrix(case.df, more.postcodes=c('CV31 1LS', 'E1 3BS'))
#' baseline.matrix = CreateBaselineMatrix(case.df, date.time.field = 'SAMPLE_DT_numeric', postcode.field = 'Patient Postcode')
CreateBaselineMatrix<-function(case.df, save.on.dir=FALSE,
                               date.time.field='week', postcode.field='postcode', 
                               more.postcodes=c(), more.weeks=c()){
  
  if (length(more.postcodes) > 0){
    postcodes = sort(unique(c(unlist(case.df[,postcode.field]), unlist(more.postcodes))))
  }else{
    postcodes = sort(unique(case.df[,postcode.field]))
  }

  n.postcodes= length(postcodes)
  maxim = max(c(case.df[,date.time.field], unlist(as.integer(more.weeks))), na.rm = T)
  minim = min(c(case.df[,date.time.field], unlist(as.integer(more.weeks))), na.rm = T)
  n.weeks = maxim - minim + 1

  baseline.matrix = matrix(data=0,
                           nrow = n.postcodes,
                           ncol = n.weeks + 1,
                           dimnames=list(postcodes, c('NA', as.character(minim:maxim))))
  
    
  time.factor = TimeFactor(case.df, save.on.dir, get.from.dir = TRUE, date.time.field = date.time.field)
  spatial.factor = PostcodeMap(matrix(data=0,
                                      nrow = n.postcodes,
                                      ncol = n.weeks + 1,
                                      dimnames=list(postcodes, c('NA', as.character(minim:maxim)) )))$Total
  
  for(j in 1:ncol(baseline.matrix)){
    for(i in 1:nrow(baseline.matrix)){
      baseline.matrix[i,j] = time.factor[j] * spatial.factor[i]
    }
  }
  baseline.matrix = baseline.matrix / sum(spatial.factor, na.rm = T)
  
  attribute_list = attributes(baseline.matrix)
  attribute_list$date.time.field = date.time.field
  attributes(baseline.matrix)<-attribute_list
  
  if(save.on.dir){
    save.and.tell("baseline.matrix", file=file.path(getwd(), paste0('baseline_matrix.RData')))
  }
  return(baseline.matrix)
}


#' Simulate an observation matrix
#'
#' Given  vectors of length M and N, representing spatial 
#' and temporal factors (e.g., population density per location and sesonal trends, respectively),
#' this function computes a MxN baseline matrix and simulate a Poisson point process.
#' 
#' @param population A \code{numeric} vector.
#' @param time.factor \code{numeric} vector.
#' @param total.average A \code{numeric}.
#' @param  save.baseline.matrix A \code{logical}.
#' @param n.unknown.time A \code{numeric}. Represents the baseline for cases that have unknown date report.
#' @return A \code{sparseMatrix} observation matrix.
#' @importFrom Matrix sparseMatrix
#' @examples
#' sim = Simulate(spatial.factor, time.factor, total.average)
Simulate<-function(population, time.factor, total.average, save.baseline.matrix=F, n.unknown.time=0){
  baseline.matrix = population %o% time.factor
  baseline.matrix = baseline.matrix / sum(baseline.matrix) * total.average 
  # colnames(baseline.matrix) = c("NA", as.character(names(time.factor)))
  
  flatten.baseline.matrix = c(baseline.matrix)
  n.col = ncol(baseline.matrix)
  simulation = rpois(length(flatten.baseline.matrix), lambda=flatten.baseline.matrix)
  simulation = matrix(simulation, ncol=n.col)
  rownames(simulation) = as.character(names(population))
  colnames(simulation) = as.character(names(time.factor))
  simulation = as(simulation, 'sparseMatrix')
  if(save.baseline.matrix){
    attribute_list = attributes(simulation)
    attribute_list$baseline.matrix = baseline.matrix
    attributes(simulation)<-attribute_list
  }
  return(simulation)
}

