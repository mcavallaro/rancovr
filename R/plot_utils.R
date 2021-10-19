library(sf)
library(spatstat)


plot_cases<-function(observation.matrix, week.range, postcode.locations, col='#ffa500', cex=0.8, pch=NULL){
  weeks<-week.range[1]:week.range[2]

  case_postcodes = apply(observation.matrix[,weeks], 1 , function(x)any(as.logical(x)))
  
  cex = apply(observation.matrix[case_postcodes,weeks], 1 , function(x){sum(x)*cex})
  
  if (is.null(pch)){
    pch=16
  }
  case_locations = postcode.locations[case_postcodes,2:3]
  points(case_locations$longitude, case_locations$latitude, col=col, pch=pch, cex=cex)
}




plot_cylinders<-function(cylinders, color="#3182bd99",lwd=0.1){ #'#ffebeb88'
  for(i in 1:nrow(cylinders)){
    plot.owin(disc(cylinders$rho[i], c(cylinders$x[i],cylinders$y[i])), add=TRUE, col=NA, lwd=lwd, border=color)
  }
}


scale2one<-function(x){
  x = order(x)
  range = range(x)
  return( (x - range[1]) / (range[2] - range[1]))
  # x = (x - mean(reference)) / sqrt(sd(reference))  + 0.5
  # x = ifelse(x < 0, 0, x)
  # x = ifelse(x > 1, 1, x)
  # return(x)
}


plot.emm.fraction<-function(case.df, emmtypes, legend.position='topright', cols=NULL, ...){
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


plot.spaghetti<-function(prospective.w.df, ...){
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


