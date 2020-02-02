# STUDIES

This repository holds a collection of finalized studies in Julia

## Index

|         Name         |        Content           |
|:--------------------:|:------------------------:|
| AJD Algos Benchmark  | benchmark for speed performance of Approximate Joint Diagnalization Algorithms |

## Content

**AJD Algos Benchmark**

The *AJD Algos Benchmark* [folder](https://github.com/Marco-Congedo/STUDIES/tree/master/AJD%20Algos%20Benchmark) holds a script
that allow to benchmark the speed performance of *approximate joint diagonalization* (AJD) algorithms implemented in the
[Diagonalizations.jl](https://github.com/Marco-Congedo/Diagonalizations.jl) package. To know more about those AJD algorithms
see [here](https://marco-congedo.github.io/Diagonalizations.jl/dev/algorithms/).

The benchmark runs the AJD of *Fourier cospectra* estimated on a database of *84 eyes-closed resting state EEG recordings*. 
For each set of cospectra the algorithms are run several times and the median execution time is retained.

As an example, the script compares two AJD algorithms optimizing the log-likelihood criterion: the original *Pham's algorithm* (2000)
and the *quasi-Newton algorithm* of Ablin et al. (1999). 

| Figure 'AJD Benchmark'  |  Legend                |
|:-----------------------:|:-----------------------|
| ![](/AJD-Algos-Benchmark/Figure.png) | On the average the quasi-Newton algorithm execute in about 50ms on these real data. It is about one order of magnitude faster as compared to Pham's algorithm. Note that Pham's algorithm ialready runs several time faster in Julia as compared to Matlab and Python |



