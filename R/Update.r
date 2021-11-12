#' Upgrade observation matrix
#' 
#' Create a new observation matrix and overwrite the old one
#'
#' @inheritParams CreateObservationMatrices
UpdateObservationMatrices<-function(case.df.new, case.df.old,
                                    types=NULL,
                                    date.time.field='week', postcode.field='postcode'){

  case.df = rbind(case.df.old, case.df.new)
  
  if (is.null(types)){
    n.types = 1
  }else{
    n.types = length(types)
  }
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
  for (e in 1:n.types){
    if (!is.null(types)){
      type=types[e]

      # UPDATE observation.matrix
      observation.matrix<-tryCatch({
        load(paste0(case.file, "observation_matrix_tmp.Rdata"))
        observation.matrix
      },
      error = function(e){
        return(NULL)
      }
      )
      if(is.null(observation.matrix)){
        writeLines("There is no `observation.matrix` for type to update, run first `CreateObservationMatrix`.")
        return(0)
      }

      cat("Step", e, "of", n.types, ", creating observation matrix for", type, "emm type\n")
      idx = case.df$type == type
      case.df = case.df[idx,]
    }
  }
}


#' Upgrade baseline matrix
#' 
#' Create a new baseline matrix to match the new observation matrix and makes new fit
#'
#' @inheritParams CreateBaselineMatrices
#' @return A 2D matrix
UpgradeBaselineMatrix<-function(case.df.new, case.df.old,
                               save.on.dir=FALSE,
                               date.time.field = 'week', postcode.field = 'postcode'){
  case.df = rbind(case.df.old, case.df.new)

  postcodes = unique(case.df[,postcode.field])
  n.postcodes= length(postcodes)
  
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
  baseline.matrix = matrix(data=0,
                              nrow = n.postcodes,
                              ncol = n.weeks + 2,
                              dimnames=list(postcodes, c('NA', as.character(0:n.weeks))))

  Parameters<-tryCatch({
    load(file.path(getwd(), paste0('timefactor_parameters_tmp.Rdata')))
    cat("Parameters for temporal baseline loaded.\n")
    Parameters
  },
  error = function(e){
    cat("Parameters for temporal baseline loaded not avaiable.")    
    cat("Recomputing old the paramaters of temporal baseline.\n")
    Parameters = cmle(case.df.old[,date.time.field], n.iterations)
    })

  time.factor = TimeFactor(case.df, save.on.dir, date.time.field, parameters = Parameters)
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


#' Upgrade baseline matrix
#' 
#' Create a new baseline matrix to match the new observation matrix and makes new a fit.
#' This function is very similar to \code{UpgradeBaselineMatrix}
#' does not perform fit if time.factor parameters were saved,
#' only create a baseline matrix to match the new observation matrix.
#'
#' @inheritParams CreateBaselineMatrices
#' @return A 2D matrix
UpdateBaselineMatrix<-function(case.df.new, case.df.old,
                               save.on.dir=FALSE,
                               date.time.field = 'week', postcode.field = 'postcode'){
  case.df = rbind(case.df.old, case.df.new)

  postcodes = unique(case.df[,postcode.field])
  n.postcodes= length(postcodes)
  
  n.weeks = max(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) - min(case.df[,date.time.field][!is.na(case.df[,date.time.field])]) + 1
  baseline.matrix = matrix(data=0,
                              nrow = n.postcodes,
                              ncol = n.weeks + 2,
                              dimnames=list(postcodes, c('NA', as.character(0:n.weeks))))

  time.factor = UpdateTimeFactor(case.df.new, save.on.dir, date.time.field)
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
    cat("Computing the temporal baseline.\n")
    Parameters = cmle(case.df[,date.time.field], n.iterations, parameters)
    })

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
}

