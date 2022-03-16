
devtools::load_all()
library(magrittr)
source("/mnt/hgfs/Documents/covid-exeedance/read_data.r")
source("/mnt/hgfs/Documents/covid-exeedance/smoother.R")
source("/mnt/hgfs/Documents/covid-exeedance/create_matrices.r")

ncols = length(300:(NCOL(baseline.matrix.f)-1))

LTLA[,c('y','x')] = vlatlong2km(LTLA[, c("latitude", "longitude")])

ws = matrix(NA, nrow = nrow(LTLA), ncol = ncols)

t = 300:(NCOL(baseline.matrix.f)-1)

for (TT in t){
  
  cylinders = CreateCylinders(Positives, baseline.matrix.p, week.range = c(TT-200, TT),
                              n.cylinders = 10000, coord.df=LTLA, only.last=T)
  I = TT -300 + 1
  ws[,I] = apply(LTLA, 1, FUN=warning.score2, TT, cylinders)
  
}


t = 300:(NCOL(baseline.matrix.p)-1)
plot(t,ws[1,], type='l', col=tab.gray, xlim=c(500,600), ylab = 'w (warning score)', main='Positives')
lines(t,ws[10,], type='l', col=tab.red)
# lines(t,ws[100,], type='l', col=tab.green)
# lines(t,ws[200,], type='l', col=tab.blue)
# lines(t,ws[300,], type='l', col=tab.orange)


norm<-function(x){
  return(x / max(x))
}
t=1:NCOL(baseline.matrix.f)
lines(t, norm(Positives[1,]), type='l', lty=3, col=tab.gray)
lines(t, norm(Positives[10,]), type='l', lty=3,col=tab.red)
# lines(t, norm(Positives[100,]), type='l', lty=3,col=tab.green)
# lines(t, norm(Positives[200,]), type='l', lty=3,col=tab.blue)
# lines(t, norm(Positives[300,]), type='l', lty=3, col=tab.orange)

lines(t, norm(colSums(baseline.matrix.p)), col='blue')

legend('topleft', c('baseline (not to scale)', 'w Hartepool', 'w Kingston', '# cases Hartep. (not to scale)', '# cases Kingst. (not to scale)'),
       col=c('blue', tab.gray,tab.red,tab.gray,tab.red), lty=c(1,1,1,3,3), cex=0.8)


idx = which(ws[, 483] > 0.95)
# LTLA[idx,1]
# [1] "South Bucks"           "Forest of Dean"        "Basingstoke and Deane" "Broxbourne"            "Hertsmere"            
# [6] "Watford"               "Ipswich"               "Stevenage" 

t = 300:(NCOL(baseline.matrix.p)-1)
plot(t,ws[idx[1],], type='l', col=tab.gray, xlim=c(700,785), ylab = 'ws (warning score)', main='Positives')
lines(t,ws[idx[2],], type='l', col=tab.red)
lines(t,ws[idx[7],], type='l', col=tab.blue)

t=1:NCOL(baseline.matrix.f)
lines(t, norm(Positives[idx[1],]), type='l', lty=3, col=tab.gray)
lines(t, norm(Positives[idx[2],]), type='l', lty=3,col=tab.red)
lines(t, norm(Positives[idx[7],]), type='l', lty=3,col=tab.blue)

lines(t, norm(colSums(baseline.matrix.p)), col='blue')
legend('bottomleft', c('baseline (not to scale)',  "ws South Bucks", "ws Forest of Dean", "ws Ipswich",
                    '# cases S.Bucks (not to scale)', '# cases FDean (not to scale)', '# Ipswich (not to scale)'),
       col=c('blue', tab.gray,tab.red,tab.blue,
             tab.gray,tab.red, tab.blue), lty=c(1,1,1,1,1,3,3,3,3), cex=0.7)


