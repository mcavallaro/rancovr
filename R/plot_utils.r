tab.blue = "#1f77b4"
tab.orange = "#ff7f0e"
tab.green = "#2ca02c"
tab.red = "#d62728"
tab.purple = "#9467bd"
tab.brown = "#8c564b"
tab.pink = "#e377c2"
tab.gray = "#7f7f7f"
tab.olive = "#bcbd22"
tab.cyan = "#17becf"

plotCases<-function(observation.matrix, week.range, postcode.locations, col='#ffa500', cex=0.8, pch=NULL){
  weeks<-week.range[1]:week.range[2]

  case_postcodes = apply(observation.matrix[,weeks], 1 , function(x)any(as.logical(x)))
  
  cex = apply(observation.matrix[case_postcodes,weeks], 1 , function(x){sum(x)*cex})
  
  if (is.null(pch)){
    pch=16
  }
  case_locations = postcode.locations[case_postcodes,2:3]
  points(case_locations$longitude, case_locations$latitude, col=col, pch=pch, cex=cex)
}

plotCylinders<-function(cylinders, color="#3182bd99",lwd=0.1){ #'#ffebeb88'
  for(i in 1:nrow(cylinders)){
    plot.owin(disc(cylinders$rho[i], c(cylinders$x[i],cylinders$y[i])), add=TRUE, col=NA, lwd=lwd, border=color)
  }
}


plotTypeFraction<-function(case.df, emmtypes, legend.position='topright', cols=NULL, ...){
  case.df.ordered = case.df[order(case.df$SAMPLE_DT_numeric),]
  case.df.ordered$cum.cases = 1:nrow(case.df)
  ylab = expression(lambda[m])
  y = 0
  lwd = 2
  l = length(emmtypes)
  if (is.null(cols)){
    cols = c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf")
  }
  for (i in 1:l){
    emmtype = emmtypes[i]
    y = c(y, range(cumsum(case.df.ordered$emmtype == emmtype)  / case.df.ordered$cum.cases))
  }
  x = case.df.ordered$SAMPLE_DT
  plot(range(x), range(y), type = 'n', ylab=ylab, ...)
  for (i in 1:l){
    emmtype = emmtypes[i]
    y = cumsum(case.df.ordered$emmtype == emmtype)  / case.df.ordered$cum.cases
    lines(x, y, col=cols[i], lwd=lwd)
  }
  legend(legend.position, lwd=lwd, col=cols[1:l], emmtypes)
}


plotThread<-function(prospective.w.df, ...){
  cols = 16:NCOL(prospective.w.df)
  idx = prospective.w.df[,NCOL(prospective.w.df)] > 0.95
  R = range( prospective.w.df[idx, 'SAMPLE_DT_numeric'])
  palette=inferno(length(R[1]:R[2]), begin=0.2, end=0.75, direction=-1)
  x = 1:len(cols)
  plot(range(x), c(0,1), type='n', xaxt='n', yaxt='n', axes = FALSE, ...)
  ret = list()
  abline(h=0.95, col="#7f7f7f")
  for (i in 1:2300){
    if( !is.na(prospective.w.df[i,'x'])){
      y = unlist(prospective.w.df[i, cols])
      idx = which(y>0)
      if(length(idx) > 0){
        if(idx[length(idx)] < length(y)){
          idx. = idx[length(idx)] + 1
          idx = c(idx,idx.)
        }        
      }
      if( y[length(y)] > 0.95){
        I = prospective.w.df[i, "SAMPLE_DT_numeric"] - R[1] + 1
        color = palette[I]
        print(c(I, color, as.character(prospective.w.df[i, "FULLNO"])))
        ret[[as.character(prospective.w.df[i, "FULLNO"])]] = color
      }else{
        color = 'black'
      }
      lines(x[idx], y[idx], col=color)
      
    }
  }
  by=1
  axis(1, at=seq(max(x),1, by=-by), labels=NA)
  by=20
  axis(1, at=seq(max(x),1, by=-by), labels=colnames(prospective.w.df)[cols][seq(max(x),1, by=-by)], tck=-0.05)
  return(ret)
}

plot.strategy.compare<-function(n.rows, prospective.w.df1, prospective.w.df2,
                                col1, col2, 
                                label1, label2, threshold=0.95, first.row=1, ...){
  ncol1 = NCOL(prospective.w.df1)
  ncol2 = NCOL(prospective.w.df2)
  cols = 16:ncol1
  
  ROW.NAMES = intersect(row.names(prospective.w.df1), row.names(prospective.w.df2))
  prospective.w.df1 = prospective.w.df1[ROW.NAMES,]
  prospective.w.df2 = prospective.w.df2[ROW.NAMES,]
  Idx = which(prospective.w.df1[, NCOL(prospective.w.df1)] > threshold)  
  x.labels = intersect(colnames(prospective.w.df1), colnames(prospective.w.df2))

  # remove the first 17 common columns
  x.labels = x.labels[-(1:17)]
  Len = length(x.labels)
  writeLines(sprintf("Max time steps is %d.", Len))
  x = 1:Len
  
  for (i in first.row:(first.row+n.rows)){
    if( !is.na(prospective.w.df1[i,'x']) & i %in% Idx){
      title = sprintf("%s %s - sampled on %d - typed on %d", prospective.w.df1[i,'Patient Postcode'], prospective.w.df1[i,'FULLNO'], prospective.w.df1[i,'SAMPLE_DT_numeric'], prospective.w.df1[i,'RECEPT_DT_numeric'] )
      plot(range(x), c(0,1), type = 'n', xaxt = 'n', main = title, ylab = expression(italic(w)),...)
      
      y = unlist(prospective.w.df1[i, x.labels])
      idx = which(y>0)
      lines(x[idx], y[idx], col = col1, lty=1)
      points(x[idx], y[idx], col = col1, pch=1)
      
      y = unlist(prospective.w.df2[i, x.labels])
      idx = which(y>0)
      lines(x[idx], y[idx], col = col2, lty=3)
      points(x[idx], y[idx], col = col2, pch=4)
      
      abline(h=0.95, col=tab.gray)

      legend("bottomright", legend = c(label1, label2), col = c(col1, col2), lty = c(1,3), pch=c(1,4))
      axis(side = 1, at = x, labels = x.labels)
    }
  }
}

plot.delay.compare<-function(n.rows, prospective.w.df1, delay.w.df2,
                                col1, col2, 
                                label1, label2, threshold=0.95, first.row=1, ...){
  ncol1 = NCOL(prospective.w.df1)
  ncol2 = NCOL(delay.w.df2)
  cols = 16:ncol1
  
  ROW.NAMES = intersect(row.names(prospective.w.df1), row.names(delay.w.df2))
  prospective.w.df1 = prospective.w.df1[ROW.NAMES,]
  delay.w.df2 = delay.w.df2[ROW.NAMES,]
  Idx = which(prospective.w.df1[, NCOL(prospective.w.df1)] > threshold)  
  x.labels = intersect(colnames(prospective.w.df1), colnames(delay.w.df2))
  
  # remove the first 17 common columns
  x.labels = x.labels[-(1:17)]
  Len = length(x.labels)
  writeLines(sprintf("Max time steps is %d.", Len))
  x = 1:Len
  
  for (i in first.row:(first.row+n.rows)){
    if( !is.na(prospective.w.df1[i,'x']) & i %in% Idx){
      title = sprintf("sampled %d - typed on %d", prospective.w.df1[i,'SAMPLE_DT_numeric'], prospective.w.df1[i,'RECEPT_DT_numeric'] )
      plot(range(x), c(0,1), type = 'n', xaxt = 'Diagnosis date',
           main = title, ylab = expression(italic(w)),...)
      
      y = unlist(prospective.w.df1[i, x.labels])
      idx = which(y>0)
      lines(x[idx], y[idx], col = col1, lty=1)
      points(x[idx], y[idx], col = col1, pch=1)
      
      y = unlist(delay.w.df2[i, x.labels])
      idx = which(y>0)
      lines(x[idx], y[idx], col = col2, lty=3)
      points(x[idx], y[idx], col = col2, pch=4)
      
      abline(h=0.95, col=tab.gray)
      
      legend("bottomright", legend = c(label1, label2), col = c(col1, col2), lty = c(1,3), pch=c(1,4))
      axis(side = 1, at = x, labels = x.labels)
    }
  }
}

#' @importFrom sf st_geometry
plotBaseMap<-function(main=NULL, add=T, boundaries=NULL, col=NA, ...){
  if (is.null(boundaries)){
    boundaries = GB_region_boundaries
  }
  plot(st_geometry(boundaries), col=col, lwd=0.8,  add=add,...)
}



plotCylindersCI<-function(cylinders, confidence.level=0.68, title=NULL){
  plot(n_obs ~ mu, cylinders, title=title, xlab = 'Expected cases', ylab = 'Observed cases', pch=20, cex=0.8)
  mu=max(cylinders$mu)
  lines(0: mu, 0:mu)
  x = seq(0, mu, 0.1)
  upper=qpois(confidence.level, x)
  lower=qpois(1-confidence.level, x)
  polygon(c(x, x[length(x):1]),c(upper, lower[length(lower):1]), col="#a6a6a666", border = "#a6a6a600")
  points(cylinders[cylinders$warning,]$mu, cylinders[cylinders$warning,]$n_obs, col='red', pch=20, cex=0.9)
}


plotBaseline<-function(week, emmtype, z, add=F, tf=factor, emmf=emmtype.factor){
  xmax=max(z$fhat) * tf[week] * emmf[emmtype][[1]]
  ColorRamp=colorRampPalette(c("white", blues9))(100)
  ColorRamp_ex=ColorRamp[1:(xmax / 3 * 100)]
  image(z$x1, z$x2,
        z$fhat * tf[week] *emmf[emmtype][[1]],
        col=ColorRamp_ex, add = add, axes=F, xlab=NA, ylab=NA)
}

# Movie<-function(cylinders, observation.matrix, postcode2coord, output.basename, emmtype){
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


# PlotCluster<-function(observation.matrix, postcode2coord.km, main=''){
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


plotCluster0<-function(X, case.df, emmtype, warning.score='warning.score', legend.pos="bottomleft", ...){
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

plotCluster<-function(X, case.df, emmtype, warning.score='warning.score', threshold=NULL, legend.pos="bottomleft", ...){
  if (is.null(threshold)){
    PlotCluster0(X, case.df, emmtype, warning.score, legend.pos, ...)
  }else{
    idx = ((case.df$emmtype == emmtype) & !is.na(case.df$x))
    case.df = case.df[idx, ]
    if(requireNamespace(viridisLite)){
      palette = viridis(max(case.df$SAMPLE_DT_numeric)+1, alpha = 1, begin = 0, end = 1)      
    }else{
      palette()
    }
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


#' centers: the data frame of centers with ID
#' radius: radius measured in kilometer
#'
make_circle<-function(x0, y0, radius, n.points = 100){
    # length per longitude changes with lattitude, so need correction
    # radiusLon <- radius /111 / cos(meanLat/57.3) 
    # radiusLat <- radius / 111
    angle = seq(0, 2*pi, length.out = n.points)

    circle.df = data.frame(
      x = x0 + radius * cos(angle),
      y = y0 + radius * sin(angle)
    )
    return(circle.df)
}

#' centers: the data frame of centers with ID
#' radius: radius measured in kilometer
#'
plotCircle<-function(cylinders, n.points = 100, add=T, ...){
  n = NROW(cylinders)
  for (i in 1:n){
    c = make_circle(cylinders[i,1], cylinders[i,2], cylinders[i,3], n.points)
    if (add==T){
      lines(c$x, c$y, ...)
    }else{
      plot(c$x, c$y, ...)
    }
  }
}