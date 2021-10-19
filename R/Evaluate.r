
source("utils.R")
#library(sf)
library(KernSmooth)
#library(spatstat)
#library(truncnorm)
source("surveillance_utils.R")
source("plot_utils.R")

#load('Data/postcode2coord.Rdata')
# UK range (latitude and longitude)
#Y.range = range(postcode2coord[!is.na(postcode2coord$latitude),]$latitude) + c(-2,2) * var(postcode2coord[!is.na(postcode2coord$latitude),]$latitude)
# [1] 50.83559 55.64878
#X.range = range(postcode2coord[!is.na(postcode2coord$longitude),]$longitude) +  c(-2,2) * var(postcode2coord[!is.na(postcode2coord$longitude),]$longitude)
# [1] -4.717408  0.348140


#' n.cylinders=10000 takes around 3 hours for the whole dataset, to end up with 300 non-empty cylinders
#' observation.matrix  and baseline.matrix have dimension 
#' This function scans the matrices
#' -1 in observation.matrix index means that we are excluding from week NA
#' @param emmtype
#' @param week.range
#' @param n.cylinders (integer) number of proposed cylinders per week interval.
#' @param rs
#' @param p.val.threshold
CreateCylinders<-function(observation.matrix, baseline.matrix, emmtype,
                                 week.range, n.cylinders=1000, rs=0.1,
                                 p.val.threshold=0.05, coord.df=postcode2coord, size_factor=1){
  observation.matrix = observation.matrix[!(rownames(observation.matrix) == 'NA'),]
  baseline.matrix = baseline.matrix[!(rownames(baseline.matrix) == 'NA'),]
  
  observation.matrix = observation.matrix[,!(colnames(observation.matrix) == 'NA')]
  baseline.matrix = baseline.matrix[,!(colnames(baseline.matrix) == 'NA')]
  
  week.range = range(as.integer(week.range))
  
  if (sum(baseline.matrix) > 2 * sum(observation.matrix)){
    print(sum(baseline.matrix))
    print(sum(observation.matrix))    
    warning(
      sprintf(
        "Warning: the baseline is too high. Have you multiplied for the emmtype factor? e.g., `CreateCylinders(observation.matrix, baseline.matrix*emmtype.factor[['%s']], '%s', %d:%d, n.cylinders=100, rs=0.1)`", emmtype, emmtype, as.integer(week.range[1]), as.integer(week.range[2])
      )
    )
  }
  init=Sys.time()
  coord.df = coord.df[!is.na(coord.df$latitude),]
  coord.km.df = coord.df
  coord.km.df[,2:3] = vlatlong2km(coord.df[,2:3])
  if ((week.range[1] < min(as.integer(colnames(observation.matrix)))) | (week.range[2] > max(as.integer(colnames(observation.matrix))))){
    A = sprintf("%d-%d", week.range[1], week.range[2])
    B = range(as.integer(colnames(observation.matrix)))
    B = sprintf("%d-%d", B[1], B[2])
    stop(paste0("`week.range`` is ", A, ", while the range of `observation.matrix` is ", B, "." ))
  }
  # weeks = as.integer(colnames(observation.matrix)[-seq(1:(starting.week+2))])
  # weeks = weeks[seq(1, length(weeks), 20)]
  # cylinders0 = data.frame(x=double(), y=double(), rho=double(), t.low=integer(), t.upp=integer(), n_obs=integer(), mu=double(), lower=integer(), upper=integer(), p.val=double(), warning=logical())
  radia_and_heights = f_radia_and_heights(baseline.matrix, 1:24) * size_factor
  cat("Evaluating cylinder exceedances from ",
        as.character(week2Date(week.range[1])),   # this interactive output doesnt work in the prospective mode
        " to ",
        as.character(week2Date(week.range[2])),        
        " for emm type ", emmtype, ".\n")
  
  # generate cylinders
  cylinders = rcylinder2(n.cylinders, observation.matrix, week.range, radia_and_heights, coord.km.df)
  if (NROW(cylinders) > 0){
    cylinders[,c('n_obs', 'mu', 'p.val')] = t(apply(cylinders, 1, compute,
                                                    observation.matrix, baseline.matrix, coord.km.df))
    ## da vettorizzare:
    # cylinders$warning = apply(cylinders, 1, function(x){ifelse((x['p.val'] < p.val.threshold) & (x['n_obs'] > 0), TRUE, FALSE)})
    # vettorizzato:
    cylinders$warning = cylinders$p.val < p.val.threshold
  }else{
    cat("No cases in the selected week range.\n")
  }
  print(Sys.time() - init)
  return(cylinders)
}

ranScanCreateCylinders<-CreateCylinders


#' 
#' n.cylinders=10000 takes around 3 hours for the whole dataset, to end up with 300 non-empty cylinders
#' observation.matrix  and baseline.matrix have dimension 
#' This function scans the matrices
#' -1 in observation.matrix index means that we are excluding from week NA
#' @inheritParams aFunction
#' @param observation.matrix.untyped
#' @param baseline.matrix.untyped
#' 
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
  init=Sys.time()
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
  cylinders = rcylinder2(n.cylinders, observation.matrix.typed + observation.matrix.untyped, week.range, radia_and_heights, coord.km.df)
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

ranScanCreateCylinders.delay<-CreateCylinders.delay


PlotCylindersCI<-function(cylinders, confidence.level=0.68, title=NULL){
  plot(n_obs.typed ~ mu.typed, cylinders, title=title, xlab = 'Expected cases', ylab = 'Observed cases', pch=20, cex=0.8)
  mu=max(cylinders$mu.typed)
  lines(0: mu, 0:mu)
  x = seq(0, mu, 0.1)
  upper=qpois(confidence.level, x)
  lower=qpois(1-confidence.level, x)
  polygon(c(x, x[length(x):1]),c(upper, lower[length(lower):1]), col="#a6a6a666", border = "#a6a6a600")
  points(cylinders[cylinders$warning,]$mu.typed, cylinders[cylinders$warning,]$n_obs.typed, col='red', pch=20, cex=0.9)
}

ranScanPlotCylindersCI<-PlotCylindersCI

SaveCylinders<-function(cylinders, file.basename){
  write.csv(cylinders, paste0(file.basename,'.csv'), quote = F, row.names = F)
  save.and.tell("cylinders", file = paste0(file.basename,'.Rdata'))
}

ranScanSaveCylinders<-SaveCylinders

Evaluate<-function(case.file, cylinders, emmtype, p.val.threshold = 0.05,
                          warning.score.name = 'warning.score', date.time.field = 'SAMPLE_DT_numeric'){
  case.df<-tryCatch({
    load(paste0(case.file, ".Rdata"))
    case.df
  },
  error = function(e){
    case.df = ranScanInit(case.file)
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

ranScanEvaluate<-Evaluate

plotBaseline<-function(week, emmtype, z, add=F, tf=factor, emmf=emmtype.factor){
  xmax=max(z$fhat) * tf[week] * emmf[emmtype][[1]]
  ColorRamp=colorRampPalette(c("white", blues9))(100)
  ColorRamp_ex=ColorRamp[1:(xmax / 3 * 100)]
  image(z$x1, z$x2,
        z$fhat * tf[week] *emmf[emmtype][[1]],
        col=ColorRamp_ex, add = add, axes=F, xlab=NA, ylab=NA)
}

# ranScanMovie<-function(cylinders, observation.matrix, postcode2coord, output.basename, emmtype){
#   idx = rownames(observation.matrix) != 'NA'
#   observation.matrix = observation.matrix[idx,]
#   idx = rownames(postcode2coord) != 'NA'
#   postcode2coord = postcode2coord[idx,]
#   z = bkde2D(cbind(postcode2coord$longitude, postcode2coord$latitude), bandwidth = 0.1, gridsize = c(400,400))
#   if (dim(postcode2coord)[1] != dim(observation.matrix)[1]){
#     writeLines("`postcode2coord` and observation matrices must have the same size")
#     return(0)
#   }
#   
#   cylinders[,c(2,1)] = vkm2latlong(cylinders[,c(1,2)])
#   
#   cRamp= colorRamp(c('Blue','Red'))
#   for(i in 2:ncol(observation.matrix)){
#     cat("week",i,"\r")
#     A=paste0(output.basename, '_', i, '_' , emmtype, '_.png')
#     png(A, width = 3.25, height = 3.25, units = 'in', res=600, pointsize = 4)
#     par(mar=c(2,2,4,2))
#     main=paste0("week ", i-2, ", emmtype ", emmtype)
#     plotBaseline(i, emmtype, z, add=F)
#     plotBaseMap(main='', add=T)
#     text(x = -3.5, y = 52.5, main) 
#     idx = observation.matrix[,i] > 0
#     warning.color = rgb2hex(cRamp(
#       sapply(1:nrow(postcode2coord[idx,]),
#              FUN=warning_ratio2, observation.matrix[idx,i],i, cylinders, postcode2coord[idx,])))
#     points(postcode2coord[idx,]$longitude, postcode2coord[idx,]$latitude, pch=20, cex=1.5, col=warning.color)
#     title(week2Date(i))
#     dev.off()
#   }
# #  shell("magick Fig/*.png Fig/movie.gif")
# }


# ranScanPlotCluster<-function(observation.matrix, postcode2coord.km, main=''){
#   library(tsne)
#   idx = rownames(observation.matrix) != 'NA'
#   observation.matrix = observation.matrix[idx,]
#   idx = rownames(postcode2coord.km) != 'NA'
#   postcode2coord.km = postcode2coord.km[idx,]
#   cases = which(observation.matrix>0, arr.ind=TRUE)
#   cases = as.data.frame(cases)
#   c3 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 3]})
#   c2 = apply(cases, 1, function(x){postcode2coord.km[x['row'], 2]})
#   cases2 = cbind(cases,c2,c3) 
#   cases2 = cases2[!is.na(cases2[,3]), ]
#   X = tsne(cases2[,2:4])
#   palette = colorRampPalette(c('blue', 'red'))(max(cases2[,2]))
#   warning.marker = ifelse(8, 20, sapply(1:nrow(postcode2coord[idx,]),
#                                         FUN=warning_ratio2, observation.matrix[idx,i],i, cylinders, postcode2coord[idx,])>0.9 )
#   
#   plot(X[,1],X[,2], xlab='C1', ylab='C2', col=colormap[cases2[,2]], pch=warning.marker, main=main)
# }


Cluster<-function(case.df, emmtype, makeplot=FALSE, warning.score='warning.score'){
  library(tsne)
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
    library(viridisLite)
    palette = viridis(max(case.df$SAMPLE_DT_numeric), begin = 0, end = 1)
    warning.marker = 20 #ifelse(case.df[,warning.score] > 0.9, 20, 1)
    plot(X[,1],X[,2], xlab='C1', ylab='C2', col=palette[case.df$SAMPLE_DT_numeric], pch=warning.marker,
         main=emmtype)
    
    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 1, begin = 0, end = 1)
    idx = ifelse(case.df[,warning.score] > 0.9, T, F)
    points(X[idx,1],X[idx,2], xlab='C1', ylab='C2', col=palette[case.df[idx, "SAMPLE_DT_numeric"]], pch=1, cex=2)
  }
  return(X)
}

ranScanCluster<-Cluster

PlotCluster0<-function(X, case.df, emmtype, warning.score='warning.score', legend.pos="bottomleft", ...){
  threshold = 0.95
  idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
  case.df = case.df[idx, ]
  palette = colorRampPalette(c('blue', 'red'))(101)
  par(fig=c(0, 1, 0, 1))
  par(mar=c(4, 4, 0.9, 4))
  #main=paste0('Embedding for emmtype ',emmtype)
  plot(X[,1], X[,2], xlab='C1', ylab='C2', col=palette[round(case.df[,warning.score] * 100+1)], pch=19, ...)
  idx2 = (case.df[,warning.score]>0.95)
  points(X[idx2,1], X[idx2,2], cex = 2, col=palette[round(case.df[,warning.score] * 100+1)][idx2])
  
  pos=legend(legend.pos, c("cases", as.expression(bquote(italic(w) ~ ">" ~ .(threshold)))),
         col = c(palette[as.integer(length(palette)/2)], 'red'), pch=c(19,1), pt.cex = c(1, 2))

  points((pos$text$x[2] + pos$rect$left)/2, pos$text$y[2], pch=19, col='red')
  
  Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
  par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
  par(mar=c(0,0,0,0))
  plot(c(0, 1), c(0, 1), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
  rasterImage(Cb[seq(length(palette), 1, -1),], 0, 0, 1, 1)
  rect(0, 0, 1, 1)
  axis(side = 4)
  mtext("warning score", side = 4, line = 2)
}

ranScanPlotCluster0<-PlotCluster0

PlotCluster<-function(X, case.df, emmtype, warning.score='warning.score', threshold=NULL, legend.pos="bottomleft", ...){
  if (is.null(threshold)){
    PlotCluster0(X, case.df, emmtype, warning.score, legend.pos, ...)
  }else{
    library(viridisLite)
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    palette = viridis(max(case.df$SAMPLE_DT_numeric)+1, alpha = 1, begin = 0, end = 1)
    par(fig=c(0, 1, 0, 1))
    par(mar=c(4, 4, 0.9, 4))
    plot(X[,1],X[,2], xlab='C1', ylab='C2', col=palette[case.df$SAMPLE_DT_numeric+1], pch=19, ...)
    # mtext(paste0('Embedding for emmtype ',emmtype), 3)
    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 1, begin = 0, end = 1)
    idx = ifelse(case.df[,warning.score] > threshold, T, F)
    points(X[idx,1], X[idx,2], xlab='C1', ylab='C2', col='black', pch=1, cex=2)
    legend(legend.pos, c("cases", as.expression(bquote(italic(w)~ ">" ~ .(threshold)))),
           col = c(palette[as.integer(length(palette)/2)], 'black'), pch=c(19,1), pt.cex = c(1, 2))
    
    Cb = matrix(rep(palette, 5), nrow = length(palette), ncol = 5)
    par(fig=c(0.86, 0.88, 0.2, 0.8), new=T)
    par(mar=c(0,0,0,0))
    plot(c(0, 1), c(1, length(palette)), yaxt='n', xaxt='n', col=NULL, frame.plot=F, ylab='time')
    rasterImage(Cb[seq(length(palette), 1, -1),], 0, 0, 1, length(palette))
    rect(0, 0, 1, length(palette))
    rr = max(case.df$SAMPLE_DT_numeric)
#    ticks = as.integer(seq(1, rr, (rr - 1)/ 10))
    ticks = c(0, 53, 53*2-1, 53*3-2, 53*4-3, 53*5-4)
    axis(side = 4, at=ticks, labels=week2Year(ticks))
#    mtext("time [week]", side = 4, line = 2)
  }
}

ranScanPlotCluster<-PlotCluster


PlotCluster3<-function(X, case.df, emmtype,
                              warning.score='warning.score',
                              threshold=0.95, legend.pos="bottomleft", ...){
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    ll = c("East Midlands", "East of England", "London", "North East",
           "North West", "South East", "South West", "West Midlands","Yorkshire&Humber")
    
    if(is.null(case.df$Postcode.Region)){
      load("Data/Area2Region_list.RData")
      
      case.df$Postcode.Region = ordered(unlist(sapply(case.df$`Patient Postcode`,
                                                     postcode.to.region,
                                                     Area2Region_list)), levels = ll)
    }
    palette = c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#0000FF")
    palette.alpha = c("#1B9E7755", "#D95F0255", "#7570B355",
                      "#E7298A55", "#66A61E55", "#E6AB0255",
                      "#A6761D55", "#66666655", "#0000FF55")
    
    case.df$color.Region = palette[case.df$Postcode.Region]
    case.df$color.Region.alpha = palette.alpha[case.df$Postcode.Region]
    
    par(fig=c(0, 1, 0, 1))
    par(mar=c(4, 4, 0.9, 4))
    
#    paste0('Embedding for emmtype ',emmtype)
#    palette = viridis(max(case.df$SAMPLE_DT_numeric), alpha = 0.1, begin = 0, end = 1)
    plot(X[,1],X[,2], xlab='C1', ylab='C2',
         col=case.df$color.Region.alpha, pch=19, ...)


    idx2 = case.df[,warning.score] > threshold
    while((sum(idx2) == 0) & (threshold > 0)){
      threshold = threshold - 0.1
      print(threshold)
      idx2 = ifelse(case.df[,warning.score] > threshold, T, F)
    }
    
    region = ordered(case.df[idx2,]$Postcode.Region)
    as.int = as.integer(region)

    points(X[idx2,1], X[idx2,2], xlab='C1', ylab='C2', col= case.df[idx2,]$color.Region,
           pch=1, cex=1.8)
    
    # if ((length(levels(region)) > 0) & (Legend == T)){
    #   legend("bottomleft", levels(region),
    #          col = palette[1:length(levels(region))],
    #          pch = 1, cex=0.7, pt.cex = 1.4)
    # }
    
    legend(legend.pos,
           c(as.character(sort(unique(case.df$Postcode.Region))),
             as.expression(bquote(italic(w) ~ '>' ~ .(threshold)))),
           col=c(palette.alpha[sort(unique(case.df$Postcode.Region))], 'black'),
           pch=c(rep(19, length(palette[sort(unique(case.df$Postcode.Region))])), 1),
           pt.cex = c(rep(1, length(palette[sort(unique(case.df$Postcode.Region))])), 1.4),
           cex=0.7
           )
}

ranScanPlotCluster3<-PlotCluster3

