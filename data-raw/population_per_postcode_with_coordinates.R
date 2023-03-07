

population.file.name = 'data-raw/Postcode_Estimates_Table_1.csv'
population.data = read.table(population.file.name, header=T, stringsAsFactors = F, sep=',', check.names=F)
population.data[,'postcode'] = c(sapply(population.data[,'Postcode'],
                                      function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))
postcode.file.name = 'data-raw/ukpostcodes.csv'
postcode.data = read.table(postcode.file.name, header=T, stringsAsFactors = F, sep=',', check.names=F)
postcode.data[,'postcode'] = c(sapply(postcode.data[,'postcode'],
                                      function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))
postcode.data$id=NULL

postcode.data = merge(postcode.data, population.data, by='postcode')
idx = order(postcode.data$postcode)
postcode.data = postcode.data[idx,]
rm(idx)
usethis::use_data(postcode.data, overwrite = T)
#save(postcode.data, file = "data/population_per_postcode_with_coordinates.RData")

# postcode.data[,c('latitude', 'longitude')] = t(apply(postcode.data, 1, postcode.to.location2))


