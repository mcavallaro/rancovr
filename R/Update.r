#' Upgrade observation matrix
#' 
#' Create a new observation matrix and overwrite the old one
#'
#' @inheritParams CreateObservationMatrices
UpdateObservationMatrices<-function(case.df.old, case.df.new, types=NULL,
                                    date.time.field='week', postcode.field='postcode'){

  tmp = rbind(case.df.old, case.df.new)
  idx_tmp = !duplicated(tmp)
  tmp = NULL
  idx_tmp = tail(idx_tmp, nrow(case.df.new))
  case.df.toadd = case.df.new[idx_tmp,]

  if (is.null(types)){
    n.types = 1
  }else{
    n.types = length(types)
  }

  all.weeks = c(case.df.old[,date.time.field], case.df.new[,date.time.field])
  week.max = max(all.weeks, na.rm = T)
  new.postcodes = unique(case.df.new[,postcode.field])

  if (is.null(types)){
    observation.matrix<-tryCatch({
      load(paste0(case.file, "observation_matrix_tmp.Rdata"))
      observation.matrix
      },
      error = function(e){return(NULL)}
    )
    if(is.null(observation.matrix)){
      writeLines("There is no `observation.matrix` for type to update, run first `CreateObservationMatrix`.")
      return(0)
    }
    # create a new sparse matrix object
    previous.postcodes = rownames(observation.matrix)
    postcodes = c(previous.postcodes, setdiff(new.postcodes, previous.postcodes))
    n.postcodes = length(postcodes)
    previous.weeks = colnames(observation.matrix)
    new.weeks = (as.integer(previous.weeks[length(previous.weeks)]) + 1):week.max
    weeks = c(previous.weeks, as.character(new.weeks))
    n.weeks = length(weeks) - 2
    observation.matrix.new = as(matrix(data = 0,
                                   nrow = n.postcodes,
                                   ncol = n.weeks + 2,
                                   dimnames = list(postcodes,
                                                 weeks)), "sparseMatrix")

    # copy the entries of the old observation matrix into the new one
    observation.matrix.new[1:nrow(observation.matrix), 1:ncol(observation.matrix)] = observation.matrix
    observation.matrix<-observation.matrix.new
    observation.matrix.new<-NULL

    # add the events in `case.df.toadd` in the new observation matrix
    for (i in 1:nrow(case.df.toadd)){
      postcode = case.df.toadd[i, postcode.field]
      week = case.df.toadd[i, date.time.field]
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

    # save on disk
    save.and.tell("observation.matrix",
                  file = file.path(getwd(), 'observation_matrix_tmp.Rdata'))
  }else{
    for (e in 1:n.types){
          type=types[e]
          cat("Step", e, "of", n.types, ", creating observation matrix for", type, "emm type\n")
      }
  }
}


#' Upgrade baseline matrix
#' 
#' Create a new baseline matrix to match the new observation matrix and makes new fit
#' for the temporal factor. 
#'
#' @inheritParams CreateBaselineMatrix
#' @return A 2D matrix
UpgradeBaselineMatrix<-function(previous.baseline.matrix, case.df.old, case.df.new,
                                save.on.dir=FALSE, date.time.field='week', postcode.field='postcode', n.iterations=10){

  case.df = rbind(case.df.old, case.df.new)
  # remove duplicates
  case.df = case.df[!duplicated(case.df)]

  week.max = max(case.df[,date.time.field], na.rm = T)
  new.postcodes = unique(case.df.new[,postcode.field])

  previous.postcodes = rownames(previous.baseline.matrix)
  postcodes = c(previous.postcodes, setdiff(new.postcodes, previous.postcodes))
  n.postcodes = length(postcodes)

  previous.weeks = colnames(previous.baseline.matrix)
  new.weeks = (as.integer(previous.weeks[length(previous.weeks)])+1):week.max
  weeks = c(previous.weeks, as.character(new.weeks))
  n.weeks = length(weeks) - 2
  baseline.matrix = matrix(data=0,
    nrow = n.postcodes,
    ncol = n.weeks + 2,
    dimnames=list(postcodes, c('NA', as.character(0:n.weeks))))

  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    cat("Old parameters for temporal baseline loaded;\n")
    cat("This will be used as starting point for the new MLE estimation.\n")
    Parameters
  },
  error = function(e){
    cat("Parameters for temporal baseline loaded not avaiable.")    
    cat("Recomputing old the paramaters of temporal baseline.\n")
    Parameters = NULL
    })

  time.factor = TimeFactor(case.df, save.on.dir, date.time.field, parameters = Parameters, n.iterations = n.iterations)
  spatial.factor = PostcodeMap(matrix(data = 0,
                                      nrow = n.postcodes,
                                      ncol = n.weeks + 2,
                                      dimnames = list(postcodes, c('NA', as.character(0:n.weeks)) )))$Total
  
  for(j in 1:ncol(baseline.matrix)){
    for(i in 1:nrow(baseline.matrix)){    
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


#' Update baseline matrix
#' 
#' Create a new baseline matrix to match the new observation matrix.
#' This function is very similar to \code{UpgradeBaselineMatrix}
#' does not perform fit if time.factor parameters were saved,
#' only create a baseline matrix to match the new observation matrix.
#'
#' @inheritParams CreateBaselineMatrix
#' @return A 2D matrix
UpdateBaselineMatrix<-function(previous.baseline.matrix, recent.case.df,
                               save.on.dir=FALSE,
                               date.time.field = 'week', postcode.field = 'postcode'){

  week.max = max(recent.case.df[,date.time.field], na.rm = T)
  new.postcodes = unique(recent.case.df[,postcode.field])

  previous.postcodes = rownames(previous.baseline.matrix)
  postcodes = c(previous.postcodes, setdiff(new.postcodes, previous.postcodes))
  n.postcodes = length(postcodes)

  previous.weeks = colnames(previous.baseline.matrix)
  new.weeks = (as.integer(previous.weeks[length(previous.weeks)])+1):week.max
  weeks = c(previous.weeks, as.character(new.weeks))
  n.weeks = length(weeks) - 2
  
  baseline.matrix = matrix(data=0,
                           nrow = n.postcodes,
                           ncol = n.weeks + 2,
                           dimnames=list(postcodes, weeks))


  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    cat("Parameters for temporal baseline loaded.\n")
    Parameters
    },
    error = function(e){
      cat("Parameters for temporal baseline loaded not avaiable.")
      cat("Run bla bla to fit the seasonal model and save the parameters.")
    })

  time.factor = predict.cmle(as.integer(weeks), Parameters)
  spatial.factor = PostcodeMap(baseline.matrix)$Total

  baseline.matrix[1:nrow(previous.baseline.matrix),1:nrow(previous.baseline.matrix)] = previous.baseline.matrix * sum(spatial.factor[1:nrow(previous.baseline.matrix)])
  
  # block updates:
  for(j in (ncol(previous.baseline.matrix)+1):ncol(baseline.matrix)){
    for(i in 1:nrow(previous.baseline.matrix)){    
        baseline.matrix[i,j] = time.factor[j] * spatial.factor[i]
    }
  }
  for(j in 1:ncol(previous.baseline.matrix)){
    for(i in (nrow(previous.baseline.matrix)+1):nrow(baseline.matrix)){    
        baseline.matrix[i,j] = time.factor[j] * spatial.factor[i]
    }
  }
  for(j in (ncol(previous.baseline.matrix)+1):ncol(baseline.matrix)){
    for(i in (nrow(previous.baseline.matrix)+1):nrow(baseline.matrix)){
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

#' Update Time factor
#'
#' Returns a vector representing the temporal component of the baseline.
#' 
#' @inheritParams TimeFactor
#' @return A vector
#' @examples
#' time.factor = TimeFactor(case.df)
UpdateTimeFactor<-function(case.df, save.on.dir = TRUE, get.from.dir = FALSE,
                     date.time.field = "week", parameters = NULL, n.iterations=20){
  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    cat("Parameters for temporal baseline loaded.\n")
    Parameters
  },
  error = function(e){
    cat("The parameteres were not save as Rdata. Recomputing it with `cmle`.\n")
    Parameters = cmle(case.df[,date.time.field], n.iterations, parameters)
    })

  n.weeks = max(case.df[,date.time.field], na.rm = T) - min(case.df[,date.time.field], na.rm = T) + 1
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
}

