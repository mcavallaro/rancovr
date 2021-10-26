

postcode.file.name = 'data-raw/UK_population_per_postcode_with_coordinates.csv'
postcode.data = read.table(postcode.file.name, header=T, stringsAsFactors = F, sep=',', check.names=F)
postcode.data[,'postcode'] = c(sapply(postcode.data[,'postcode'],
                                          function(x){gsub(" ", "", toupper(x),  fixed = TRUE)}))

diamonds <- read_csv("data-raw/diamonds.csv", col_types =
                       list(
                         cut = col_factor(c("Fair", "Good", "Very Good", "Premium", "Ideal"), TRUE),
                         color = col_factor(c("D", "E", "F", "G", "H", "I", "J"), TRUE),
                         clarity = col_factor(c("I1", "SI2", "SI1", "VS2", "VS1", "VVS2", "VVS1", "IF"), TRUE)
                       )
)

devtools::use_data(diamonds, overwrite = TRUE)

save(postcode.data, )