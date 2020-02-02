# STUDIES

This repository holds a collection of finalized studies in Julia

## Index

|         Name         |         Content          |
|:--------------------:|:------------------------:|
| [AJD Algos Benchmark](AJD Algos Benchmark) | benchmark for speed performance of Approximate Joint Diagnalization Algorithms |


### AJD Algos Benchmark

The *AJD Algos Benchmark* [folder](https://github.com/Marco-Congedo/STUDIES/tree/master/AJD%20Algos%20Benchmark) holds a script
that allow to benchmark the speed performance of approximate joint diagonalization (AJD) algorithms implemented in the
[Diagonalizations.jl](https://github.com/Marco-Congedo/Diagonalizations.jl) package. The benchmark runs the AJD of Fourier
cospectra estimated on a database of 84 eyes-closed resting state EEG recordings. For each set of cospectra the algorithms
are run several times and the median time in milliseconds is retained.
