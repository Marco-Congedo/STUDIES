#   Script "choleskyBenchmark.jl" for Julia language
#
#   MIT License
#   Copyright (c) 2020
#   Marco Congedo, CNRS, UGA, Grenoble-INP, Grenoble, France
#   https://sites.google.com/site/marcocongedo/
#
# DIRECTIONS
#
# If you have never used Julia before, read first the accompanying file
# "Install-Julia-README.md" located at
#  https://github.com/Marco-Congedo/STUDIES.
#
# Before running the script, install the required packages. For that,
# execute the following line in Julia's REPL and wait until
# all packages are installed:
#    ]add Test, BenchmarkTools, PosDefManifold
#
# In order to run this script in ATOM, hit CTRL-A and then ENTER.
#
# ? CONTENT
#
# This scipt benchmarks Julia's `cholesky` method based on LAPACK
# against a simple modified symmetric Gaussian elimination procedure
# based on the seminal work of Alan Turing (1948).
# Read the accompanying `README.md` file for details.

using LinearAlgebra, Test, BenchmarkTools, PosDefManifold


"""
	trifact(P::Union{Hermitian, Symmetric, Matrix, LowerTriangular};
		   	check::Bool = true,
		   	tol::Real = √eps(real(eltype(P))))

 Compute the Cholesky factorization of a dense positive definite
 matrix `P`. Return the `LowerTriangular` factor ``L``,
 such that ``LL^H=P``.

 Input matrix `P` may be of type `Matrix`, `Symmetric` or `Hermitian`.
 Since only the lower triangle is used, `P` may also be a `LowerTriangular`
 view of a positive definite `Matrix`.

 The algorithm is based on the symmetric version of
 Alan Turing's *Gaussian elimination*.

 ## Examples
 	# test
	using Test, LinearAlgebra, PosDefManifold

	function testLLt(type::Type, n, etol)
		P = type <: Real ? randP(n) : randP(type, n)
		L = trifact(P)
		@test(norm(L*L'-P)/√n < etol)
	end

	testLLt(Float64, 10, 1e-9)		# real matrix
	testLLt(ComplexF64, 10, 1e-9)	# complex matrix

	# Benchmark
	using LinearAlgebra, BenchmarkTools
	n = 100
	P = randP(n)
	@benchmark(trifact(P))
	@benchmark(cholesky(P))
	Pc = randP(ComplexF64, n)
	@benchmark(trifact(Pc))
	@benchmark(cholesky(Pc))
"""
trifact( P::Union{Hermitian, Symmetric, Matrix, LowerTriangular};
	   	 check::Bool = true,
	   	 tol::Real = √eps(real(eltype(P)))) =
    trifact!(copy(Matrix(P)); check=check, tol=tol)


"""
    trifact!(P::AbstractMatrix{T};
			  tol::Real = √eps(T)) where T<:RealOrComplex
 The same thing as [`utrifact!`](@ref), but destroys the input matrix.
 `P` must be a general `Matrix` comprised of real or complex elements.
"""
function trifact!(	P::Matrix{T};
			  	 	check::Bool = true,
					tol::Real = √eps(real(T))) where T<:Union{Real, Complex}
	LinearAlgebra.require_one_based_indexing(P)
	n = LinearAlgebra.checksquare(P)

	@inbounds for j=1:n-1
		check && abs2(P[j, j])<tol && throw(LinearAlgebra.PosDefException(1))
		f = P[j, j]
		P[j, j] = sqrt(f)
		g = P[j, j]
		for i=j+1:n
			θ = P[i, j] / f
			cθ = conj(θ)
			for k=i:n P[k, i] -= cθ * P[k, j] end # update P and write D
			P[i, j] = θ*g # write L factor
		end
	end
	P[n, n] = sqrt(P[n, n]) # write last element of L factor

	return LowerTriangular(P)
end


# test
n=5
P=randP(n)
L=trifact(P)
@test(L*L'≈P)
P=randP(ComplexF64, n)
L=trifact(P)
@test(L*L'≈P)


# benchmark
N=[5, 10, 20, 30, 40, 50, 75, 100, 150, 200, 300]
a=Matrix{Float64}(undef, length(N), 2)
for (i, n) ∈ enumerate(N)
	println("running benchmark ", i, "(n=", n, ") of ", length(N))
	global P=randP(n)
	a[i, 1]=@belapsed(trifact(P))
	a[i, 2]=@belapsed(cholesky(P))
end

b=(a.*1_000_000)

using Plots
plotly()
plot(N, b; labels=["Gaussian elimination" "cholesky Julia-LAPACK"],
	xlabel="matrix size", ylabel="execution time in μs",
	xguidefontsize=12, yguidefontsize=12,
	xtickfontsize=10, ytickfontsize=10,
	yscale=:log10,
	xscale=:log10,
	legend=:top, legendfontsize=11)
