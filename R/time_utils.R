library(optimx)

lambda<-function(x, params){
  # $$ 
  #   \lambda(t) = a + b + b \sin(c\,2 \pi/365 + x 2 \,\pi / 365)
  # $$
  # $$
  #   l(\lambda) = \log L(\lambda) = - \int_0^T\lambda(x) dx + \sum_{i=1}^n\log \lambda(x_i)
  # $$
  #       
    
  p1<-params[1]
  p2<-params[2]
  p3<-params[3]
  p4<-params[4]
   # period = 365 # (days)
  period = 52 # (weeks)
  
  p1 + p2 + (p2 + p4 * x)* sin(2 * pi / period * p3 + x * 2 * pi / period) + p4* x
}

neg.log.like<-function(params, data){
  lower<-min(data)
  upper<-max(data)
  integrate(lambda,lower,upper,params,subdivisions = 5000)$value - sum(log(lambda(data,params)))
}

cmle1<-function(data, cpar, iterations=10000){
  fn<-function(par){
    if (all(par>c(0,2,0)) & all(par<c(100, 100, 60)) ){
      res = neg.log.like(c(par, cpar), data)
    }
    # else if ((par[4] < -0.1) | (par[4] > 0.1)){
    #   res = 1000000 * par[4] ^ 2
    # }
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

cmle2<-function(data, cpar){
  fn<-function(par){
    res = neg.log.like(c(cpar, par), data)
  }
  # else if ((par[4] < -0.1) | (par[4] > 0.1)){
  #   res = 1000000 * par[4] ^ 2
  # }
  res = optimize( fn, interval = c(0,10), tol = 0.00001)
  return(res$minimum)
}


cmle<-function(data, n.cycles=10, start=NULL){
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
  res = c(par, par1)
  trace[i,] = res
  attributes(res) = list(trace=trace)
  return(res)
}


predict.cmle<-lambda

predict.cmle.ci<-function(x, parameters){
  data.frame(upper=qpois(0.975,lambda(x,parameters)),
             lower=qpois(0.025,lambda(x,parameters)))
}

expand.histogram= function(histo){
  Histo = data.frame(Names = as.integer(names(histo)), freq = unlist(histo))
  return(unlist(do.call(c, apply(Histo, 1, function(x){rep(x[1], x[2]) }) )))
}

Loess = function(counts, span = 0.2){
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
