# Rancovr: Cluster detection in R with Random Neighbourhood Covering

`rancovr` is a statistical software package written in R for disease
cluster and anomaly detection. It implements the Random Neighbourhood
Covering (RaNCover) approach of reference \[1\]. RaNCover assigns a
score *w* ∈ \[0, 1\] to each records. A high score suggests that the
record is likely to be part of a cluster (e.g., it is an infection case
caused by a local outbreak), while a low score suggests that the record
is consistent with a baseline of sporadic cases.

``` r
install.packages("devtools")
devtools::install_github("mcavallaro/rancovr")
```

As a demonstration, we consider the spatio-temporal coordinates in
`Data/synthetic_dataset.csv`. The entries of this data set represent
(simulated) locations and detection dates of infected patients generated
by an endemic disease (`end.`) and cases due to a local outbreak
(`epi.`) in England. See also reference \[1\] for the simulation
details.

``` r
data("simulation_data")
head(simulation_data)
```

    ##   postcode week population  sim type latitude  longitude Postcode Total
    ## 1  AL100AZ   59         67 epi.    1 51.76421 -0.2309368  AL100AZ    67
    ## 2  AL100AZ   41         67 epi.    2 51.76421 -0.2309368  AL100AZ    67
    ## 3  AL100AZ   51         67 epi.    2 51.76421 -0.2309368  AL100AZ    67
    ## 4  AL100DR   50         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ## 5  AL100DR   47         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ## 6  AL100DR   51         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ##          y         x
    ## 1 5756.180 -8.074648
    ## 2 5756.180 -8.074648
    ## 3 5756.180 -8.074648
    ## 4 5756.123 -8.254003
    ## 5 5756.123 -8.254003
    ## 6 5756.123 -8.254003

``` r
data("GB_region_boundaries")
plotBaseMap(add=F, xlim=range(simulation_data$longitude), ylim=range(simulation_data$latitude))
points(simulation_data$longitude, simulation_data$latitude,
       col=ifelse(simulation_data$sim=='epi.', tab.red, tab.blue),
       pch=ifelse(simulation_data$sim=='epi.', 1, 20),
       cex=ifelse(simulation_data$sim=='epi.', 0.6, 0.2))
legend('topright',c('end.','epi.'), pch=c(20,1), col=c(tab.blue, tab.red))
```

![](README_files/figure-markdown_github/unnamed-chunk-4-1.png) For
convenience, all observations are arranged in a `sparseMatrix` object
named `observation.matrix` and saved on disk.

``` r
CreateObservationMatrices(simulation_data)
```

    ## The variable `observation.matrix` has been saved on disk in file `/home/massimo/Documents/rancovr/observation_matrix.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/observation_matrix.RData", verbose=1)`.

Observations must be compared with an appropriate baseline model. If the
numbers of these observations significantly exceeded the model
prediction, an outbreak might be occurring. Estimating the baseline
involves finding a temporal trend (using the function `TimeFactor`) and
a spatial trend based on the spatial population distribution.

``` r
load(file.path(getwd(), "observation_matrix.RData"), verbose=1)
```

    ## Loading objects:
    ##   observation.matrix

``` r
time.factor = TimeFactor(simulation_data, n.iterations=5)
```

    ## Computing the temporal baseline.
    ## Estimating parameters for temporal trend, step  1  of  5 .Estimating parameters for temporal trend, step  2  of  5 .Estimating parameters for temporal trend, step  3  of  5 .Estimating parameters for temporal trend, step  4  of  5 .Estimating parameters for temporal trend, step  5  of  5 .The variable `Parameters` has been saved on disk in file `/home/massimo/Documents/rancovr/timefactor_parameters.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/timefactor_parameters.RData", verbose=1)`.
    ## The variable `time.factor` has been saved on disk in file `/home/massimo/Documents/rancovr/timefactor.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/timefactor.RData", verbose=1)`.

``` r
baseline.matrix = CreateBaselineMatrix(simulation_data, save.on.dir = T)
```

    ## Temporal baseline loaded.
    ## Compiling the table that maps the rows of the observation/baseline matrix to geo-coordinates and population.
    ## Loading objects:
    ##   postcode2coord

    ## Warning in as.character(postcode2coord[, postcode.field]) == rownames(matrix):
    ## longer object length is not a multiple of shorter object length

    ## Data loaded from `postcode2coord.RData` is for a different matrix and will be overwritten by the map for the current matrix.
    ## The variable `postcode2coord` has been saved on disk in file `/home/massimo/Documents/rancovr/postcode2coord.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/postcode2coord.RData", verbose=1)`.
    ## The variable `baseline.matrix` has been saved on disk in file `/home/massimo/Documents/rancovr/baseline_matrix.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/baseline_matrix.RData", verbose=1)`.

``` r
load(file.path(getwd(), "observation_matrix.RData"), verbose=1)
```

    ## Loading objects:
    ##   observation.matrix

``` r
plot(time.factor, xlab = 'Week', ylab='Number of cases', xaxt='n')
# lines(colSums(baseline.matrix))
points(Matrix::colSums(observation.matrix), pch='+')
axis(side=1, at=1:length(time.factor), labels = names(time.factor))
legend('bottomright',legend=c('Baseline', 'Observations'), pch=c('o', '+'))
```

![](README_files/figure-markdown_github/unnamed-chunk-7-1.png)

Create 100,000 cylinders to cover the observed cases using the estimated
baseline.

``` r
cylinders = CreateCylinders(observation.matrix, baseline.matrix, week.range = c(0,99), n.cylinders = 100000)
```

    ## Compiling the table that maps the rows of the observation/baseline matrix to geo-coordinates and population.
    ## Loading objects:
    ##   postcode2coord
    ## Using data loaded from `postcode2coord.RData`
    ## Time difference of 1.860501 mins

``` r
head(cylinders)
```

    ##            x        y      rho t.low t.upp n_obs         mu     p.val warning
    ## 1  -4.029426 5725.620 7.358866    52    68    48 68.1715332 0.9957064   FALSE
    ## 2   3.555001 5730.511 8.163928    69    82    30 37.8589987 0.9170454   FALSE
    ## 3 -93.948032 5728.697 8.875126    86    97     6  4.4450999 0.2877209   FALSE
    ## 4 -81.656709 5934.709 7.358866    16    32     7  6.7911243 0.5187256   FALSE
    ## 5 -83.854503 6023.165 6.938005     5    23     2  0.7122419 0.1600714   FALSE
    ## 6 -75.820322 5833.160 8.163928    40    53     1  3.3078229 0.9634042   FALSE

Some cylinders contain much more cases than the baseline predicts. These
cylinders cover epidemic (outbreak) events.

``` r
plotCylindersCI(cylinders, confidence.level = 0.95)
```

![](README_files/figure-markdown_github/unnamed-chunk-9-1.png)

The “true” baseline matrix used to generate the endemic events is
available as `data()`. Let’s use it in place of the estimated baseline
matrix. Notice that the true baseline matrix has higher dimensionality
than the estimated baseline matrix (it includes entries for more
postcodes and times) and requires a matching observation matrix.

``` r
print(dim(baseline.matrix))
```

    ## [1] 3446  101

``` r
print(dim(observation.matrix))
```

    ## [1] 3446  101

``` r
data(baseline_for_sim)
print(dim(baseline_for_sim))
```

    ## [1] 10000   101

``` r
CreateObservationMatrices(simulation_data,
                          more.postcodes=rownames(baseline_for_sim),
                          more.weeks=colnames(baseline_for_sim))
```

    ## Warning in unlist(as.integer(more.weeks)): NAs introduced by coercion

    ## Warning in unlist(as.integer(more.weeks)): NAs introduced by coercion

    ## The variable `observation.matrix` has been saved on disk in file `/home/massimo/Documents/rancovr/observation_matrix.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/observation_matrix.RData", verbose=1)`.

``` r
load("/home/massimo/Documents/rancovr/observation_matrix.RData", verbose=1)
```

    ## Loading objects:
    ##   observation.matrix

``` r
print(dim(observation.matrix))
```

    ## [1] 10000   101

``` r
cylinders.2 = CreateCylinders(observation.matrix, baseline_for_sim, week.range = c(0,99), n.cylinders = 100000)
```

    ## Compiling the table that maps the rows of the observation/baseline matrix to geo-coordinates and population.
    ## Loading objects:
    ##   postcode2coord

    ## Warning in as.character(postcode2coord[, postcode.field]) == rownames(matrix):
    ## longer object length is not a multiple of shorter object length

    ## Data loaded from `postcode2coord.RData` is for a different matrix and will be overwritten by the map for the current matrix.
    ## The variable `postcode2coord` has been saved on disk in file `/home/massimo/Documents/rancovr/postcode2coord.RData`.
    ## Load on memory with `load("/home/massimo/Documents/rancovr/postcode2coord.RData", verbose=1)`.
    ## Time difference of 3.149846 mins

``` r
head(cylinders.2)
```

    ##            x        y      rho t.low t.upp n_obs        mu      p.val warning
    ## 1  -4.565732 5782.775 14.56173    66    78     7  5.230448 0.27222953   FALSE
    ## 2 -39.431829 6112.945 13.99046    12    25    27 17.817712 0.02536801    TRUE
    ## 3 -89.611223 5979.242 16.81444    56    65    12  9.115878 0.20837666   FALSE
    ## 4 -67.432738 5837.926 16.81444    17    26    19 23.455035 0.84773404   FALSE
    ## 5  -2.764860 5740.000 12.61083    10    26    69 79.836310 0.89999698   FALSE
    ## 6 -45.162555 5864.023 11.57249    17    36    13  8.069146 0.06718268   FALSE

``` r
plotCylindersCI(cylinders.2, confidence.level = 0.95)
```

![](README_files/figure-markdown_github/unnamed-chunk-11-1.png)

Compute the warning scores for each case:

``` r
simulation_data[,'warning.score'] = apply(simulation_data, 1, FUN=warning.score, cylinders)
simulation_data[,'warning.score.2'] = apply(simulation_data, 1, FUN=warning.score, cylinders.2)
head(simulation_data)
```

    ##   postcode week population  sim type latitude  longitude Postcode Total
    ## 1  AL100AZ   59         67 epi.    1 51.76421 -0.2309368  AL100AZ    67
    ## 2  AL100AZ   41         67 epi.    2 51.76421 -0.2309368  AL100AZ    67
    ## 3  AL100AZ   51         67 epi.    2 51.76421 -0.2309368  AL100AZ    67
    ## 4  AL100DR   50         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ## 5  AL100DR   47         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ## 6  AL100DR   51         64 epi.    1 51.76370 -0.2360576  AL100DR    64
    ##          y         x warning.score warning.score.2
    ## 1 5756.180 -8.074648     0.8916479       0.8492063
    ## 2 5756.180 -8.074648     0.8717340       0.7913043
    ## 3 5756.180 -8.074648     0.9716216       0.9369658
    ## 4 5756.123 -8.254003     0.9852744       0.9680511
    ## 5 5756.123 -8.254003     0.9869754       0.9783105
    ## 6 5756.123 -8.254003     0.9744624       0.9368984

Assess concordance with ROC-AUC:

``` r
library(pROC)
```

    ## Type 'citation("pROC")' for a citation.

    ## 
    ## Attaching package: 'pROC'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     cov, smooth, var

``` r
ROC = roc(ifelse(simulation_data$sim == 'end.', FALSE, TRUE), simulation_data$warning.score)
```

    ## Setting levels: control = FALSE, case = TRUE

    ## Setting direction: controls < cases

``` r
plot(ROC)
print(ROC$auc)
```

    ## Area under the curve: 0.9997

``` r
ROC = roc(ifelse(simulation_data$sim == 'end.', FALSE, TRUE), simulation_data$warning.score.2)
```

    ## Setting levels: control = FALSE, case = TRUE
    ## Setting direction: controls < cases

``` r
plot(ROC, add=T, col='red')
print(ROC$auc)
```

    ## Area under the curve: 0.9993

``` r
legend('bottomright', legend =  c('Using estimated baseline', 'Using true baseline'), lty=1, col=c('black','red'))
```

![](README_files/figure-markdown_github/unnamed-chunk-13-1.png)

With mean squared error:

``` r
simulation_data$Y = ifelse(simulation_data$sim == 'epi.',1,0)
simulation_data$sqerr = (simulation_data$Y - simulation_data$warning.score)^2
cat("MSE using estimated baseline:", mean(simulation_data$sqerr), '\n') 
```

    ## MSE using estimated baseline: 0.02540727

``` r
simulation_data$sqerr.2 = (simulation_data$Y - simulation_data$warning.score.2)^2
cat("MSE using true baseline:", mean(simulation_data$sqerr.2), '\n') 
```

    ## MSE using true baseline: 0.03199637

And with a map:

``` r
data("GB_region_boundaries")
#plotBaseMap(add=F, xlim=c(-0.6,0.6), ylim=c(51.648,51.65))
plotBaseMap(add=F, xlim=c(-1,1), ylim=c(50.648,52.65))
points(simulation_data$longitude, simulation_data$latitude,
       col=ifelse(simulation_data$sim=='epi.', tab.red, tab.blue),
       pch=ifelse(simulation_data$sim=='epi.', 4, 20),
       cex=ifelse(simulation_data$sim=='epi.', 1, 0.5))
points(simulation_data[simulation_data$warning.score>0.95,]$longitude,
       simulation_data[simulation_data$warning.score>0.95,]$latitude,
       col=tab.orange,
       pch=1,
       cex=1)
# case.df$color = rgb(colorRamp(c("blue", "red"))(case.df$warning.score) / 255)
# plot(case.df$longitude, case.df$latitude, col=case.df$color)

legend('topright',c('end.','true epi.', 'w>0.95'), pch=c(20,4,1), col=c(tab.blue, tab.red, tab.orange))
```

![](README_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
plotBaseMap(add=F, xlim=c(-0.6,0.6), ylim=c(51.648,51.65))
#plotBaseMap(add=F, xlim=c(-1,1), ylim=c(50.648,52.65))
points(simulation_data$longitude, simulation_data$latitude,
       col=ifelse(simulation_data$sim=='epi.', tab.red, tab.blue),
       pch=ifelse(simulation_data$sim=='epi.', 4, 20),
       cex=ifelse(simulation_data$sim=='epi.', 1, 0.5))
points(simulation_data[simulation_data$warning.score.2>0.95,]$longitude,
       simulation_data[simulation_data$warning.score.2>0.95,]$latitude,
       col=tab.orange,
       pch=1,
       cex=1)
# case.df$color = rgb(colorRamp(c("blue", "red"))(case.df$warning.score) / 255)
# plot(case.df$longitude, case.df$latitude, col=case.df$color)

legend('topright',c('end.','true epi.', 'w>0.95'), pch=c(20,4,1), col=c(tab.blue, tab.red, tab.orange))
```

![](README_files/figure-markdown_github/unnamed-chunk-16-1.png)

\[1\] M. Cavallaro, J. Coelho, D. Ready, V. Decraene, T. Lamagni, N. D.
McCarthy, D. Todkill, M. J. Keeling (2022) Cluster detection with random
neighbourhood covering: Application to invasive Group A Streptococcal
disease. PLoS Comput Biol 18(11): e1010726.
<https://doi.org/10.1371/journal.pbci.1010726>
