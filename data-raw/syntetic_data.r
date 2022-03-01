# simulate endemic component
set.seed(1)
data("UK_population_per_postcode_with_coordinates")
idx = sample(NROW(postcode.data), 10000)
postcode.data. = postcode.data[idx,]

idx = order(postcode.data.$postcode)
sample.population = postcode.data.[idx,'Total']
names(sample.population)<-postcode.data.$postcode[idx]
rm(postcode.data.)

time.factor = lambda(0:99, c(34, 16, 55,  0), names=T)
total.average = 5000
simulation = Simulate(sample.population, time.factor, total.average, save.baseline.matrix=T)
baseline_for_sim=attributes(simulation)$baseline.matrix
save(baseline_for_sim, file = "data/baseline_for_sim.RData")


simulation_data = as.data.frame(which(simulation == 1, arr.ind = TRUE))
simulation_data$postcode = rownames(simulation)[simulation_data$row]
for (i in 2:max(unlist(simulation))){
  simulation_data.tmp = as.data.frame(which(simulation == i, arr.ind = TRUE))
  simulation_data.tmp$postcode = rownames(simulation)[simulation_data.tmp$row]
  simulation_data = rbind(simulation_data, simulation_data.tmp[rep(seq_len(NROW(simulation_data.tmp)), i),])
}
simulation_data$col = simulation_data$col - 2
# subtract 2 because the first two columns of baseline.matrix correspond c(NA, 0).

simulation_data = cbind(simulation_data, sample.population[simulation_data[,'row']])
names(simulation_data) = c('row','week', 'postcode', 'population')
simulation_data$sim='end.'

# simulate epidemic component
set.seed(1)
idx1 = grepl('AL1', names(sample.population), fixed=T)
idx1 = idx1 | grepl('AL2', names(sample.population), fixed=T)
epi = rpois(n = sum(idx1) * 20, lambda=rep(dnorm(-9:10, sd = 4) * 4, rep(sum(idx1), 20)))
epi.matrix = matrix(epi, ncol = 20)
rownames(epi.matrix) = row.names(simulation)[1:sum(idx1)]
colnames(epi.matrix) = 40:59 

new.simulation_data = as.data.frame(which(epi.matrix == 1, arr.ind = TRUE))
new.simulation_data$postcode = rownames(epi.matrix)[new.simulation_data$row]
for (i in 2:max(unlist(epi.matrix))){
  simulation_data.tmp = as.data.frame(which(epi.matrix == i, arr.ind = TRUE))
  simulation_data.tmp$postcode = rownames(epi.matrix)[simulation_data.tmp$row]
  new.simulation_data = rbind(new.simulation_data, simulation_data.tmp[rep(seq_len(NROW(simulation_data.tmp)), i),])
}
rm(simulation_data.tmp)
new.simulation_data$col = new.simulation_data$col + 39
new.simulation_data = cbind(new.simulation_data, sample.population[new.simulation_data[,'row']])
names(new.simulation_data) = c('row','week', 'postcode', 'population')
new.simulation_data$sim = 'epi.'

simulation_data = rbind(simulation_data, new.simulation_data)
rm(new.simulation_data)

# simulate serotypes
simulation_data$type = '1'
simulation_data$type[sample(1:NROW(simulation_data),  round(0.2 * NROW(simulation_data)))] = '2'

simulation_data = merge(simulation_data[,-1], postcode.data, by = 'postcode')
simulation_data[, c('y','x')] = vlatlong2km(simulation_data[,c('latitude', 'longitude')])
simulation_data = simulation_data[, -c(6,7,8,9)]

#write.table(simulation_data, file = 'data/simulation_data.csv', sep=';', row.names = FALSE)
save(simulation_data, file = "data/simulation_data.RData")
