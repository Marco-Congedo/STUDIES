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
# ? CONTENT
#
# This scipt benchmarks the quasi-Newton Log-Likelihood Approximate
# Joint Diagonalization (AJD) algorithm of Ablin, Cardoso and
# Gramfort (2019) and the log-likelihood algorithm of Pham (2000)
# on a real EEG database using the julia language.
# Read the accompanying `READ.ME` file for details.

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
