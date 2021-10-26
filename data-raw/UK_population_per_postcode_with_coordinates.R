

postcode.file.name = 'data-raw/UK_population_per_postcode_with_coordinates.csv'
postcode.data = read.table(postcode.file.name, header=T, stringsAsFactors = F, sep=',', check.names=F)
postcode.data[,'postcode'] = c(sapply(postcode.data[,'postcode'],
                                          function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))

save(postcode.data, file = "data/UK_population_per_postcode_with_coordinates.Rdata")
