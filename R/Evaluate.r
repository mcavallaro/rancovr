#' Create cylinders to cover events.
#'
#' n.cylinders=10000 takes around 3 hours for the whole dataset, to end up with 300 non-empty cylinders
#' observation.matrix  and baseline.matrix have dimension
#' This function scans the matrices
#' -1 in observation.matrix index means that we are excluding from week NA
#' @param observation.matrix A 2D matrix.
#' @param baseline A 2D matrix or a Nx4 expand.grind data frame
#' @param week.range week range.
#' @param n.cylinders (integer) number of proposed cylinders per week interval.
#' @param p.val.threshold A numeric.
#' @param rs A numeric.
#' @param size_factor An integer.
CreateCylinders<-function(observation.matrix, baseline,
                          week.range, n.cylinders=1000, rs=0.1,
                          p.val.threshold=0.05,
                          size_factor=1){
  
  # load("~/Documents/Rancovr/Data/postcode2coord.Rdata")  
  coord.df = PostcodeMap(observation.matrix)
  #postcode2coord
  test = all(dim(baseline) == dim(observation.matrix))
  if (test){
    baseline = baseline[!(rownames(baseline) == 'NA'),]
    baseline = baseline[,!(colnames(baseline) == 'NA')]
    radia_and_heights = f_radia_and_heights(baseline, 1:24) * size_factor    
    if (sum(baseline) > 2 * sum(observation.matrix)){
      print(sum(baseline))
      print(sum(observation.matrix))    
      warning(
        sprintf(
          "Warning: the baseline might be too high. If you are working with typed data, please check that it is scaled by type factor."
        )
      )
    }
  }else{
    radia_and_heights = f_radia_and_heights_(baseline, 1:24) * size_factor    
  }

  observation.matrix = observation.matrix[!(rownames(observation.matrix) == 'NA'),]
  observation.matrix = observation.matrix[,!(colnames(observation.matrix) == 'NA')]

  init = Sys.time()
  coord.df = coord.df[!is.na(coord.df$latitude),]
  coord.km.df = coord.df
  coord.km.df[,2:3] = vlatlong2km(coord.df[,2:3])
  week.range = range(as.integer(week.range))
  if ((week.range[1] < min(as.integer(colnames(observation.matrix)))) | (week.range[2] > max(as.integer(colnames(observation.matrix))))){
    A = sprintf("%d-%d", week.range[1], week.range[2])
    B = range(as.integer(colnames(observation.matrix)))
    B = sprintf("%d-%d", B[1], B[2])
    stop(paste0("`week.range`` is ", A, ", while the range of `observation.matrix` is ", B, "." ))
  }
  # week.range = range(
  #   as.integer(week.range)) + 1 - min(as.integer(colnames(observation.matrix))
  #   )


  cat("Evaluating cylinder exceedances from ",
        as.character(week2Date(week.range[1])),   # this interactive output doesnt work in the prospective mode
        " to ",
        as.character(week2Date(week.range[2])), ".\n")
  
  # # generate cylinders
  # print(week.range)
  # week.range = week.range + 2 - min(as.integer(colnames(observation.matrix)), na.rm=T)

  cylinders = rcylinder(n.cylinders, observation.matrix, week.range, radia_and_heights, coord.km.df)

#  print(any(cylinders$t.low == cylinders$t.upp))
  if (NROW(cylinders) > 0){
    if (test){
      tmp = t(apply(cylinders, 1, compute, observation.matrix, baseline, coord.km.df))
    }else{
      tmp = t(apply(cylinders, 1, compute.from.tab.baseline, observation.matrix, baseline, coord.km.df))
    }

    cylinders[,c('n_obs', 'mu', 'p.val')] = tmp
    cylinders$warning = cylinders$p.val < p.val.threshold
  }else{
    cat("No cases in the selected week range.\n")
  }
  print(Sys.time() - init)
  return(cylinders)
}


#' 
#' n.cylinders=10000 takes around 3 hours for the whole dataset, to end up with 300 non-empty cylinders
#' observation.matrix  and baseline.matrix have dimension 
#' This function scans the matrices
#' -1 in observation.matrix index means that we are excluding from week NA
#' @inheritParams CreateCylinders
#' @param observation.matrix.untyped A 2D matrix.
#' @param baseline.matrix.untyped A 2dD matrix.
CreateCylinders.delay<-function(observation.matrix.typed, baseline.matrix.typed,
                                observation.matrix.untyped, baseline.matrix.untyped,
                                emmtype,
                                week.range, n.cylinders=1000, rs=0.1,
                                p.val.threshold=0.05, coord.df=postcode2coord,  size_factor=1){
  observation.matrix.typed = as.matrix(observation.matrix.typed[!(rownames(observation.matrix.typed) == 'NA'),])
  baseline.matrix.typed = as.matrix(baseline.matrix.typed[!(rownames(baseline.matrix.typed) == 'NA'),])
  observation.matrix.untyped = as.matrix(observation.matrix.untyped[!(rownames(observation.matrix.untyped) == 'NA'),])
  baseline.matrix.untyped = as.matrix(baseline.matrix.untyped[!(rownames(baseline.matrix.untyped) == 'NA'),])

  week.range = range(as.integer(week.range))
  c1 = sum(baseline.matrix.typed) > 2 * sum(observation.matrix.typed)
  c2 = sum(baseline.matrix.untyped) > 2 * sum(observation.matrix.untyped)
  if (c1){
    warning("Warning: the typed baseline is very high.")
  }
  if (c2){
    warning("Warning: the untyped baseline seems is very high.")
  }
  init = Sys.time()
  coord.df = coord.df[!is.na(coord.df$latitude),]
  coord.km.df = coord.df
  coord.km.df[,2:3] = vlatlong2km(coord.df[,2:3])
  if ((week.range[1] < min(as.integer(colnames(observation.matrix.typed)))) | (week.range[2] > max(as.integer(colnames(observation.matrix.typed))))){
    # line = sprintf("Try with `starting.week` < %d", ncol(observation.matrix)-2)
    # writeLines(c("`starting.week` is bigger than the matrix length.", line))
    A = sprintf("%d-%d", week.range[1], week.range[2])
    B = range(as.integer(colnames(observation.matrix.typed)))
    B = sprintf("%d-%d", B[1], B[2])
    writeLines(paste0("Error: `week.range` is ", A, ", while the range of `observation.matrix` is ", B, "." ))
    return(NA)
  }
  # weeks = as.integer(colnames(observation.matrix)[-seq(1:(starting.week+2))])
  # weeks = weeks[seq(1, length(weeks), 20)]
  # cylinders0 = data.frame(x=double(), y=double(), rho=double(), t.low=integer(), t.upp=integer(), n_obs=integer(), mu=double(), lower=integer(), upper=integer(), p.val=double(), warning=logical())
  radia_and_heights = f_radia_and_heights(baseline.matrix.typed, 1:24) * size_factor
  cat("Evaluating cylinder exceedances from ",
      as.character(week2Date(week.range[1])),   # this interactive output doesnt work in the prospective mode
      " to ",
      as.character(week2Date(week.range[2])),        
      " for emm type ", emmtype, ".\n")

  # generate cylinders
  cylinders = rcylinder(n.cylinders, observation.matrix.typed + observation.matrix.untyped, week.range, radia_and_heights, coord.km.df)
  if (NROW(cylinders) > 0){
    cylinders[,c('n_obs.typed', 'mu.typed', 'p.val.typed')] = t(apply(cylinders, 1, compute,
                                                                observation.matrix.typed,
                                                                baseline.matrix.typed,
                                                                coord.km.df))
    cylinders[,c('n_obs.untyped', 'mu.untyped', 'p.val.untyped')] = t(apply(cylinders, 1, compute,
                                                                            observation.matrix.untyped,
                                                                            baseline.matrix.untyped,
                                                                            coord.km.df))
  }else{
    cat("No (typed) cases in the selected week range.\n")
  }
  cylinders$warning = (cylinders$p.val.typed < p.val.threshold) | (cylinders$p.val.untyped < p.val.threshold)
  print(Sys.time() - init)
  return(cylinders)
}


SaveCylinders<-function(cylinders, file.basename){
  write.csv(cylinders, paste0(file.basename,'.csv'), quote = F, row.names = F)
  save.and.tell("cylinders", file = paste0(file.basename,'.Rdata'))
}

Evaluate<-function(case.file, cylinders, emmtype, p.val.threshold = 0.05,
                          warning.score.name = 'warning.score', date.time.field = 'SAMPLE_DT_numeric'){
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){
    case.df = Init(case.file)
    return(case.df$case.df)
  })
  if (!('x' %in% names(case.df) & ('y' %in% names(case.df)))){
    writeLines("Inserting coordinates...")
    if (!('latitude' %in% names(case.df) & ('longitude' %in% names(case.df)))){
      case.df[,c("latitude", "longitude")] = t(apply(case.df, 1, postcode.to.location2))      
    }
    case.df[,c("y", "x")] = vlatlong2km(case.df[,c("latitude", "longitude")])
    writeLines("...Done.")
    save.and.tell('case.df', file=paste0(case.file, ".Rdata"))
  }
  idx = case.df$emmtype == emmtype
  writeLines(paste0("Computing warning scores for emmtype ", emmtype, "..."))
  case.df[idx, warning.score.name] = apply(case.df[idx,], 1, FUN=warning.score, cylinders, date.time.field)
  writeLines("...Done.")
  return(case.df)
}


Cluster<-function(case.df, emmtype, makeplot=FALSE, warning.score='warning.score'){
  if(requireNamespace("tsne")){
    # idx = rownames(observation.matrix) != 'NA'
    # observation.matrix = observation.matrix[idx,]
    # idx = rownames(postcode2coord.km) != 'NA'
    # postcode2coord.km = postcode2coord.km[idx,]
    # cases = which(observation.matrix>0, arr.ind=TRUE)
    # cases = as.data.frame(cases)
    # c3 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 3]})
    # c2 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 2]})
    # cases2 = cbind(cases,c2,c3)
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    X = tsne(case.df[,c('x', 'y', 'SAMPLE_DT_numeric')])
    # palette = colorRampPalette(c('blue', 'red'))(max(case.df$SAMPLE_DT_numeric))
    if (makeplot == TRUE){
      if(requireNamespace(viridisLite)){
        palette = viridis(max(case.df$SAMPLE_DT_numeric), begin = 0, end = 1)
      }else{
        
      }
      warning.marker = 20 #ifelse(case.df[,warning.score] > 0.9, 20, 1)
      plot(X[,1],X[,2], xlab='C1', ylab='C2', col=palette[case.df$SAMPLE_DT_numeric], pch=warning.marker,
           main=emmtype)
      
      palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 1, begin = 0, end = 1)
      idx = ifelse(case.df[,warning.score] > 0.9, T, F)
      points(X[idx,1],X[idx,2], xlab='C1', ylab='C2', col=palette[case.df[idx, "SAMPLE_DT_numeric"]], pch=1, cex=2)
    }
    return(X)
  }else{
    writeLines("This function requires tsne package.")  
  }
}
