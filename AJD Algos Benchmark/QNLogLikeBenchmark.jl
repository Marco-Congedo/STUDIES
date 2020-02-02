#   Script "QNLogLikeBenchmark.jl" for Julia language
#
#   MIT License
#   Copyright (c) 2020
#   Marco Congedo, CNRS, UGA, Grenoble-INP, Grenoble, France
#   https://sites.google.com/site/marcocongedo/
#
# DIRECTIONS
#
# If you have never used Julia before, read first the accompanying file
# "WhatsJulia_ReadME.txt".
#
# Before running the script, install the required packages. For that,
# execute the following line in Julia's REPL and wait until
# all packages are installed:
#    ]add Statistics, LinearAlgebra, Dates, Random, MAT, Diagonalizations
#
# In order to run this script in ATOM, hit CTRL-A and then ENTER.
#
# ? CONTENTS
#
# This scipt benchmarks the quasi-Newton Log-Likelihood Approximate
# Joint Diagonalization (AJD) algorithm of Ablin, Cardoso and
# Gramfort (2019) and the log-likelihood algorithm of Pham (2000)
# on a real EEG database using the julia language.
#
# The accompanying file `Cospectra.mat` holds Fourier cospectra
# computed on 84 eyes-closed resting-state EEG recordings obtained
# on healthy individuals (sampling rate: 128, 19 electrodes placed
# according to the 10-20 int. system, manual artefact rejection of
# corrupted epochs).
#
# For each recording, data has been pre-whitened reducing the dimension,
# if needed, so as to keep the smallest dimension explaining at least
# 99.9% of the total variance. Then, Fourier cospectra have been estimated
# by means of the Welch methods, using 256-sample long 50% overlapping
# epochs and a Harris tapering window (see the documentation of the
# "FourierAnalysis.jl" package), yielding copsectra at 0.5Hz resolution.
#
# The 47 cospectra in the frequency band-pass region 1-24 Hz are submitted
# to the two AJD algorithms. This is a real-case-scenarion usage of AJD
# algorithms as the diagonalizations of such matrices yields a blind source
# separation method specifically tailored for EEG data (Congedo et al., 2008).
#
# For both algorithms the maximum number of iterations is set to 1000
# and the tolerance for convergence is set to 1e-6.
# In order to estimate execution time, for each one of the 84 sets of
# cospectra the two algorithms are run in randomized order `ntrials`
# times and the minimum completion time across trials is retained.
#
# The plot at the end of the script display the minimum completion time
# for the dataset sorted in ascending order for each algorithm.
#
# All computations for the benckmark are carried out on a single
# logical processor. The AJD algorithms are implemented in the
# "Diagonalizations.jl" package.
#
# REFERENCES
# P. Ablin, J.F. Cardoso, A. Gramfort (2019) Beyond Pham's algorithm
# for joint diagonalization, Proc. ESANN Conference.
# https://hal.archives-ouvertes.fr/hal-01936887v1
#
# Congedo M, Gouy-Pailler C, Jutten C (2018) On the blind source
# separation of human electroencephalogram by approximate joint
# diagonalization of second order statistics.
# Clinical Neurophysiology 119, 2677-2686.
# https://hal.archives-ouvertes.fr/hal-00343628/document
#
# Diagonalizations.jl
# https://github.com/Marco-Congedo/Diagonalizations.jl
#
# FourierAnalysis.jl.
# https://github.com/Marco-Congedo/FourierAnalysis.jl
#
# D.-T. Pham (2000) Joint approximate diagonalization of positive definite
# matrices, SIAM Journal on Matrix Analysis and Applications, 22(4), 1136–1152.
# https://pdfs.semanticscholar.org/0cb5/ca9de76b8893a2549ec278942bb6a5a37a35.pdf?_ga=2.131902607.124228321.1577610632-183244851.1563047228

using Statistics, LinearAlgebra, Dates, Random, MAT, Diagonalizations

const ntrials = 15
const ajdAlgoritms = [:QNLogLike, :LogLike]
# here above you can add other AJD algorithms supported in Diagonalizations.jl


# Given input matrices in array C, number of trials and algorithms provided
# in array of symbols `arg`, return the minimum execution time in milliseconds
# across # trials. For each trial the order ot algorithms is randomized
function bench1(C, trials, alg)
   t = zeros(Int64, trials, length(alg))
   order = randperm(length(alg))
   for i ∈ 1:trials, j ∈ 1:length(alg)
      ⌚=now()
      ajd(C; algorithm=alg[order[j]], sort=false, tol=1e-6, maxiter=1000)
      t[i, order[j]]=(now()-⌚).value # get the time in ms
   end
   return [minimum(t[:, i]) for i ∈ 1:size(t, 2)]
end

# read cospectra
vars = matread((@__DIR__)*"/Cospectra.mat")
nEEG = length(vars["S"]) # cospectra sets
R = zeros(Int64, length(ajdAlgoritms), nEEG)

# function to get the cospectra from 1 to 24 Hz in 0.5 Hz step as
# a vector f `Hermitian` matrices for subject 's'
getCospectra(data, s) =
   Vector{Hermitian}([Hermitian(data[s][i]) for i ∈ 1:length(vars["S"][s])])

# run all benchmarks
println("\nBenchmarking ", [string(a) for a ∈ ajdAlgoritms], " algorithms")
for s ∈ 1:nEEG
   print(rpad("EEG file $s of $nEEG,...", 21))
   𝗦 = getCospectra(vars["S"], s)
   # compute and store benchmark's results
   R[:, s]=bench1(𝗦, ntrials, ajdAlgoritms)
   println(" times(ms): ", R[:, s])
end

# plot the time in milliseconds sorted in ascending order separatedly for
# each AJD algorithm
using Plots
plotly()
# pyplot(), GR,...
📉 =  plot(sort(R', dims=1),
           label=["Ablin et al. (2018)" "Pham (2000)"],
           xtickfontsize=13, ytickfontsize=13,
           xguidefontsize=16, yguidefontsize=16,
           linewidth=2, titlefontsize=15,
           #title="Median execution time of AJD algorithms based on log-likelihood criterion",
           xlabel="EEG data sets", yscale=:log10, ylabel="time (ms)",
           legend=:top, legendfontsize=15)

R # print results in Julia's REPL
