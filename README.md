`rancovr`. Cluster detection with Random Neighbourhood Covering
===============================================================

`rancovr` is a statistical software package written in R for the detection of disease clusters based on the Random Neighbourhood Covering (RaNCover) of reference \[1\]. `rancovr` assess whether a single recorded infection is part of a disease cluster (such that caused by a local outbreak) or is consistent with a baseline of sporadic cases.

As a demonstration, we consider the spatio-temporal coordinates stored in `Data/synthetic_dataset.csv`, which represent records of infection cases and is obtained aggregating data simulated from an endemic component (`end.`) and from an outbreak (`epi.`) in UK. See reference \[1\] for the details.

``` r
case.df = read.csv(file='Data/synthetic_dataset.csv', sep = ',')
head(case.df)
```

    ##   week postcode latitude longitude population        y         x
    ## 1    0  B14 6TN 52.42452 -1.906414         80 5829.607 -63.43771
    ## 2    0  B15 2BQ 52.47047 -1.908348         59 5834.716 -63.27718
    ## 3    0  B42 2RZ 52.53757 -1.902311        177 5842.178 -62.74949
    ## 4    0  B61 0DB 52.34426 -2.052947         18 5820.682 -68.73608
    ## 5    0  B91 3GX 52.40486 -1.775618         17 5827.421 -59.17484
    ## 6    0  BH178AN 50.75126 -1.961944        128 5643.541 -73.64748
    ##   warning.score  sim
    ## 1     0.1538462 end.
    ## 2     0.1339286 end.
    ## 3     0.1220339 end.
    ## 4     0.1150442 end.
    ## 5     0.1717557 end.
    ## 6     0.1230769 end.

    > # write.table(case.df[,1:8], file='Data/Simulation_experiment_case_df.csv', quote=T, sep=',', row.names = F)
    > case.df = read.table(file='Data/Simulation_experiment_case_df.csv',  sep=',', header=T) 
    > #
    > # saveRDS(b.matrix, 'Data/Simulation_experiment_baseline_matrix.rds')
    > baseline.matrix = readRDS('Data/Simulation_experiment_baseline_matrix.rds')
    > #
    > # saveRDS(sim2, 'Data/simulation_experiment_simulation_matrix.rds')
    > simulation.matrix = readRDS('Data/simulation_experiment_simulation_matrix.rds')
    > #
    > # write.table(df.cases2, file='Data/Simulation_experiment_coordinates_df.csv', quote=T, sep=',', row.names = F)
    > coordinate.df = read.table(file='Data/Simulation_experiment_coordinates_df.csv',  sep=',', header=T)

    > head(case.df)
      row SAMPLE_DT_numeric postcode latitude longitude population        y           x
    1  483                 0  BB120EZ 53.79724 -2.264622        157 5982.253  -67.346355
    2  488                 0  BB2 1HN 53.74896 -2.496710        113 5976.884  -74.560339
    3  555                 0  BB8 7AR 53.86487 -2.163336        109 5989.774  -63.955308
    4  973                 0  BN274EW 50.88265  0.263503        234 5658.151    9.803753
    5 1079                 0  BS106DD 51.50526 -2.603483         37 5727.385  -92.748101
    6 1142                 0  BS247EQ 51.34984 -2.917601        164 5710.103 -105.092027

    > baseline.matrix[1:6,1:6]
                     0           1           2           3
    AL1 1TA 0.008866186 0.009163252 0.009446268 0.009711113
    AL1 1UB 0.005344003 0.005523056 0.005693641 0.005853273
    AL1 2JT 0.002793456 0.002887052 0.002976221 0.003059666

    > simulation.matrix[1:6,1:6]
            0 1 2 3 4 5
    AL1 1TA 0 0 0 0 0 0
    AL1 1UB 0 0 0 0 0 0
    AL1 2JT 0 0 0 0 0 0
    AL1 4XG 0 0 0 0 0 0
    AL1 5DF 0 0 0 0 0 0
    AL1 5JQ 0 0 0 0 0 0

    > head(coordinate.df[1:3,])
            Postcode latitude longitude
    AL1 1TA  AL1 1TA 51.74220 -0.320578
    AL1 1UB  AL1 1UB 51.73746 -0.316577
    AL1 2JT  AL1 2JT 51.73965 -0.340915
    AL1 4XG  AL1 4XG 51.76646 -0.317093
    AL1 5DF  AL1 5DF 51.75159 -0.306273
    AL1 5JQ  AL1 5JQ 51.74941 -0.310961

Create 100,000 cyclinders to cover the detected cases:

    > source('R/Init2.R')
    > cylinders = CreateCylinders(observation.matrix = simulation.matrix, baseline.matrix = baseline.matrix,
    emmtype = 'sim', week.range = c(0,99), n.cylinders = 100000, coord.df=coordinate.df)

Compute the warning score for each case:

    > source('R/surveillance_utilis.R)
    > case.df[,'warning.score'] = apply(case.df, 1, FUN=warning.score, cylinders)
    > head(case.df)
       row SAMPLE_DT_numeric postcode latitude longitude population        y           x warning.score
    1  483                 0  BB120EZ 53.79724 -2.264622        157 5982.253  -67.346355   0.007518797
    2  488                 0  BB2 1HN 53.74896 -2.496710        113 5976.884  -74.560339   0.018292683
    3  555                 0  BB8 7AR 53.86487 -2.163336        109 5989.774  -63.955308   0.011764706
    4  973                 0  BN274EW 50.88265  0.263503        234 5658.151    9.803753   0.000000000
    5 1079                 0  BS106DD 51.50526 -2.603483         37 5727.385  -92.748101   0.014184397
    6 1142                 0  BS247EQ 51.34984 -2.917601        164 5710.103 -105.092027   0.014285714

\[1\] M. Cavallaro, J. Coelho, D. Ready, V. Decraene, T. Lamagni, N. D. McCarthy, D. Todkill, M. J. Keeling, Cluster detection with random neighbourhood covering: application to invasive Group A Streptococcal disease, 2021
