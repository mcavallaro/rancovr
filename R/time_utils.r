#' Lambda
#'  
#' Seasonal model for the intensity function, \eqn{\lambda(t) = a + b + b \sin(c\,2 \pi/365 + x 2 \,\pi / 365)}.
#' 
#' @param x \code{numeric} vector of time points.
#' @param params \code{numeric} vector of parameters of length = 4.
#' @param  names A vector. Can be used to label 
#' @return a \code{numeric} vector of the same length as \code{x}.
#' @examples
#' lambda(c(1,2,3,4,5,6), c(1,1,1,1) )
lambda<-function(x, params, names=F){    
  p1<-params[1]
  p2<-params[2]
  p3<-params[3]
  p4<-params[4]
   # period = 365 # (days)
  period = 52 # (weeks)
  
  ret = p1 + p2 + (p2 + p4 * x)* sin(2 * pi / period * p3 + x * 2 * pi / period) + p4* x
  
  if(names){
    ret = c(0, ret)
    names(ret) = c('NA', 0:(length(ret)-2))
    return(ret)
  }else{
    return(ret)    
  }
}


#' Negative log-likelihood
#'
#' \eqn{l(\lambda) = \log L(\lambda) = - \int_0^T\lambda(x) dx + \sum_{i=1}^n\log \lambda(x_i)}
#'
#' @param params A \code{numeric} vector of parameters of length = 4.
#' @param data A \code{numeric} vector of event times.
#' @return \code{numeric}.
#' @examples
#' neg.log.like(c(1,1,1,1),c(1,2,3,4,5,6))
neg.log.like<-function(params, data){
  lower<-min(data)
  upper<-max(data)
  integrate(lambda,lower,upper,params,subdivisions = 5000)$value - sum(log(lambda(data,params)))
}



#' Optimize function \code{neg.log.like} conditioned on having the last parameter fixed.
#' 
#' @param data A \code{numeric} vector of event times.
#' @param cpar A \code{numeric}.
#' @param iteration This is the maximum number of interations allowed for the function \code{optimx}.
#' @return 2D 4x2 matrix.
#' @importFrom optimx optimx
#' @examples
#' cmle1(data, 1, iterations=20)
cmle1<-function(data, cpar, iterations=10000){
  fn<-function(par){
    if (all(par>c(0,2,0)) & all(par<c(100, 100, 60)) ){
      res = neg.log.like(c(par, cpar), data)
    }
    else{
      res = Inf
    }
    return(res)
  }
  n = 10
  start<-data.frame(p1=runif(n,0, 100), p2=runif(n, 2, 100), p3=runif(n, 1, 60))
  opt<-function(init){
    res<-optimx(as.vector(init),
                fn,
                method='Nelder-Mead',
                control=list(maxit=iterations, dowarn=FALSE)) # all.methods=TRUE
    # if (res$convcode != 0){
    #   cat("optimisation did not converge with code:", res$convergence, "\n")
    # }
    return(c(res$value, res$convcode, unlist(res[,1:3])))
  }
  res = apply(start, 1, opt)
  id = which.min(res[1,])
  return(res[3:5, id])
}

#' Optimize function \code{neg.log.like} conditioned on having the first three parameters fixed.
#' 
#' @param data A \code{numeric} vector of event times.
#' @param cpar \code{numeric}.
#' @return \code{numeric}.
#' @examples
#' cmle2(data, c(1,1,1))
cmle2<-function(data, cpar){
  fn<-function(par){
    res = neg.log.like(c(cpar, par), data)
  }
  # else if ((par[4] < -0.1) | (par[4] > 0.1)){
  #   res = 1000000 * par[4] ^ 2
  # }
  res = optimize(fn, interval = c(0,10), tol = 0.00001)
  return(res$minimum)
}

#' Estimate the parameters for the seasonal model given the aggregated vector of event times.
#'
#' Optimize the function \code{neg.log.like} recursively calling \code{clme1} and \code{cmle2}.
#' 
#' @param data A \code{numeric} vector of event times.
#' @param n.cycles \code{integer}. Number of times the conditional optimizers are called.#
#' @param start \code{numeric}. Starting parameters (not yet implemented).
#' @param save.on.dir \code{logical}. If TRUE, will save the inferred parameter in `timefactor_parameters.RData`.
#' @return \code{numeric}. A 2D 2x4 matrix with the estimated parameters.
#' @examples
#' cmle(data)
cmle<-function(data, n.cycles=10, start=NULL, save.on.dir=TRUE){
  par = cmle1(data, 0)
  trace = matrix(data=0, nrow = n.cycles, ncol = 4)
  for (i in 1:n.cycles){
    cat("Estimating parameters for temporal trend, step ", i, " of ", n.cycles, ".\r")
    trace[i,1:3] = par
    par1 = cmle2(data, par)
    trace[i,4] = par1
    par = cmle1(data, par1)
  }
  par1 = cmle2(data, par)
  Parameters = c(par, par1)
  trace[i,] = Parameters
  attributes(Parameters) = list(trace=trace)
  if(save.on.dir){
    save.and.tell('Parameters', file=file.path(getwd(), paste0('timefactor_parameters.RData')))
  }
  return(Parameters)
}


#' Lambda
#'  
#' Seasonal model for the intensity function, \eqn{\lambda(t) = a + b + b \sin(c\,2 \pi/365 + x 2 \,\pi / 365)}.
#'
#' @inheritParams lambda
#' @return a \code{numeric} vector of the same length as \code{x}.
#' @examples
#' predict.cmle(c(1,2,3,4,5,6), c(1,1,1,1) )
predict.cmle<-lambda

predict.cmle.ci<-function(x, parameters){
  data.frame(upper=qpois(0.975,lambda(x,parameters)),
             lower=qpois(0.025,lambda(x,parameters)))
}

expand.histogram<-function(histo){
  Histo = data.frame(Names = as.integer(names(histo)), freq = unlist(histo))
  return(unlist(do.call(c, apply(Histo, 1, function(x){rep(x[1], x[2]) }) )))
}

Loess<-function(counts, span = 0.2){
  df = data.frame(time = 1:length(counts), counts = counts)
  loess(counts ~ time, df, span=span)
}

predict.Loess<-function(x, fit){
  df = data.frame(time = x)
  predict(fit, df)
}

predict.Loess.ci<-function(x, fit){
  df = data.frame(time = x)
  data.frame(upper= unlist(-predict(fit, df, se = T)[['se.fit']] + predict(fit, df, se = T)[['fit']]),
             lower= unlist(predict(fit, df, se = T)[['se.fit']] + predict(fit, df, se = T)[['fit']]))
}
