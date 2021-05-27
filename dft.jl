### A Pluto.jl notebook ###
# v0.14.7

using Markdown
using InteractiveUtils

# ╔═╡ 1d5ea4e7-dd3c-4583-a7dd-6bed050ada7a
begin
	# Install packages in a local sandbox:
    import Pkg
    Pkg.activate(mktempdir())
	ENV["PYTHON"] = "" # if empty, create own Conda environment for Julia
    Pkg.add([
        Pkg.PackageSpec(name="PyCall", version="1.92"),
        Pkg.PackageSpec(name="Conda", version="1.5"),
        Pkg.PackageSpec(name="LaTeXStrings", version="1"),
    	Pkg.PackageSpec(name="OffsetArrays", version="1.9"),
    ])
    using PlutoUI, LaTeXStrings
end

# ╔═╡ cb7a88cb-5f29-4feb-9560-e77013967dcd
using OffsetArrays

# ╔═╡ 3be5aa2a-92b2-4dbe-a9ff-bfafa2f6ecb0
using PyCall

# ╔═╡ c14506f3-c3f6-42a5-86f8-93812bd5e789
begin
	# Install numpy
	using Conda
	Conda.add(["numpy"])
end

# ╔═╡ 372dc7ad-5d86-4e09-a058-dc3c243bbb82
md"""
### Installation

If running locally we will first install Python and Conda and a couple Julia packages in a sandbox environment...just hang tight for 20 seconds or so...
"""

# ╔═╡ 8daa7390-bf0a-11eb-2f09-cb503508dca6
md"""
# DFT Challenge

This is a challeng to write the Discrete Fourier Tranform in any language to see how well it is suitable to electrical engineering tasks.  

There are many ways to write a faster DFT but this challenge is to just use the basic algorithm and to try to write it by hand like it is described in a textbook.  There are many issues that an electrical engineer will run into where they will have to write their own solution to and not depend on code written by someone else.

## DFT Description from a textbook:

$$H_k[k] = \sum_{n=0}^{N-1}x[n] e^{\frac{-j2\pi k n}{N}} \text{, where } k = 0 \ldots N-1$$

So for N samples of vector `x` there will be `N` results in the DFT vector `Hₖ`.  For the first results where `k=0` we have:

$$H_0[0] = \sum_{n=0}^{N-1}x[n] e^{0} \text{,    where } k = 0$$

And the second when `k=1` we have:

$$H_1[1] = \sum_{n=0}^{N-1}x[n] e^{\frac{-j2\pi 1 n}{N}} \text{,    where } k = 1$$

And so on until `k = N-1`.

## Julia Implementation
"""

# ╔═╡ ffb8cd4a-21f8-4507-945d-99b6967c9788
function DFT(x)
	N = length(x)
	Hₖ = [sum(x[n+1]*ℯ^(-im*2π*n*k/N) for n in 0:N-1) for k in 0:N-1]
end

# ╔═╡ 4a8624aa-2dad-4ce1-8442-32ec454930d6
md"""> **Note!** Julia uses `im` for the imaginary number and by default uses `1-based` indexing.  
> To write `π` type `\pi<tab>`, for Euler's number type `\euler<tab>` or use `exp()`, and for `Hₖ` type `H\_k<tab>`."""

# ╔═╡ 9f4a2c5f-a8b7-4fbf-93a4-2bc129ae7207
md"### Create some test data

The test data is a simle 2000 point sinusoild...nothing special."

# ╔═╡ 38c57971-eabc-4ecf-948b-ffbcffdf517d
t = range(0, stop=1, length=2000)

# ╔═╡ a0fd711d-2c86-4b30-b5b8-b1222c60486e
vsin = sin.(2pi .* t)

# ╔═╡ 883216a0-a527-4b7a-b266-b6b08e321c6c
md"""
> **Note:** Use "dot" broadcasting in Julia where any function can be executed element-by-element of the arguments inside.  Above `sin.` and `.*` operate element by element so it is like having a compact `for` loop to take the `sin` of each element of the vector passed to it.
"""

# ╔═╡ ac8f9841-8bd3-4913-8773-95923761c969
md"### Run Julia and record the execution time"

# ╔═╡ bd9842bb-dcca-4f75-a916-dbfe51a8469b
dft1 = DFT(vsin)

# ╔═╡ 4394b038-a5e5-4054-8bbb-b78d49329f00
t_julia = @elapsed DFT(vsin)

# ╔═╡ 3e57dd3b-c179-4d1e-aeec-42ef7fc7b9ca
md"""
The basic Julia version finishes in **$(round(t_julia, sigdigits=4))** seconds
"""

# ╔═╡ c6c22775-9d3a-4cf8-ba84-df3e7ea51b03
md"""
>**Note:** the **`2000`** point DFT must execute the loop `2000*2000` or **4 million** times.  Each loop has about 6 multiplies, 1 divide, and 1 `exp` of an imaginary number.  `exp` is a large function so it is hard to know how many operations are needed.  Fortunately, Julia can tell us by using the [GFlops.jl](https://github.com/triscale-innov/GFlops.jl) package.  It reports about 90 FLOPs per iteration so that is about **400 million** FLOPs for the complete DFT.
"""

# ╔═╡ e22aded1-cf5b-43ae-a923-fedca30f0fec
md"""
### Improved Julia version

Julia has optimized functions that are common in engineering but are not found in most programming languages.

In this case we can use the `cispi` function which is the same as:

$$\mathrm{cispi}(x) = \cos(\pi x) + \texttt{i} \sin(\pi x) = e^{-\texttt{i} \pi x}$$

where `i` is imaginary.  `cispi(x)` is a faster version of `exp(j*pi*x)` (see [Wikipedia](https://en.wikipedia.org/wiki/Cis_(mathematics)) for more info).

So let's rewrite it with `cispi`:
"""

# ╔═╡ 2ed07614-1715-4646-a09c-c7653a0102db
function DFT2(x)
	N = length(x)
	Hₖ = [sum(x[n+1]*cispi(-2*n*k/N) for n in 0:N-1) for k in 0:N-1]
end

# ╔═╡ 7be5b003-e8b4-4795-ad1e-350699b360fc
dft2 = DFT2(vsin)

# ╔═╡ 8350eed7-0e88-4055-8060-49bf8931c84e
t_julia2 = @elapsed DFT2(vsin)

# ╔═╡ bce85b30-df7f-4697-a9d6-2c2c281f5a4c
dft1 ≈ dft2  # are they the same (within floating point error)?

# ╔═╡ 418d36e3-0d7a-4b8a-aef6-ea0704e72320
md"""So the `cispi` version took $(round(t_julia2, sigdigits=3)) seconds.  It gets the same result and is $(round(t_julia/t_julia2, sigdigits=3))x faster."""

# ╔═╡ 485c624a-d950-46b5-8ed6-b0a100ae5008
md"""
### Improved Julia with 0-based indexing

The first index of a DFT is typically `0` for DC so users would like to have the first index be `0`.  Also, some users would like to use 0-based indexing and Julia can support this through a custom defined type.  The `OffsetArrays` package defines an array with user-defined indexing.  Let's use this to see how it works.
"""

# ╔═╡ 0f771326-6994-40a3-86c8-7dc49256d634
function DFT3(x::OffsetVector)
	N = length(x)
	# Main loop:
	Hₖ = [sum(x[n]*cispi(-2*n*k/N) for n in 0:N-1) for k in 0:N-1]
	OffsetVector(Hₖ, 0:N-1)
end

# ╔═╡ 7034ecce-8ac5-45cc-9786-9294b9c99fd9
vsin_0 = OffsetVector(vsin, 0:length(vsin)-1)

# ╔═╡ 3f6c22d2-aefd-4b3e-bdc5-09b15a02a516
dft3 = DFT3(vsin_0)

# ╔═╡ 2e9845bc-a6ad-4b12-962f-35d1def78e64
t_julia3 = @elapsed DFT3(vsin_0)

# ╔═╡ c35ce815-f1b9-4367-9f81-b2bd51850b12
md"""The 0-based index version took $(round(t_julia3, sigdigits=3)) seconds."""

# ╔═╡ b4b92f64-276f-4bdb-8b68-23406efc774a
md"""
## Python Implementation

We will use the `PyCall` package to call Python from Julia

"""

# ╔═╡ 470dfa78-66c8-40a7-90dd-749667244207
begin
	py"""
	import cmath
	import math
	def DFT_py(x):
		j = complex(0, 1)
		N = len(x)
		return [sum(x[n]*cmath.exp(-2*j*math.pi*n*k/N) for n in range(0,N)) for k in range(0,N)]
	"""
	DFT_py = py"DFT_py"  # copy python function over so it exists on the Julia side
end

# ╔═╡ 74abbc3f-784e-474c-bfdc-1e19718ca6cb
t_python = @elapsed DFT_py(vsin)  # this takes a long time to run (eg > 10 seconds)

# ╔═╡ bb6d5b1c-9382-43ed-8725-c0ab1b513479
md"""
The Python version took $(round(t_python, sigdigits=3)) seconds.
"""


# ╔═╡ 22de6bec-5774-40e4-9685-9347de93b21b
md""" ### Improved Numpy Version

We will write the same algorithm using Numpy.
"""

# ╔═╡ 55ed8230-7ab1-41dd-9337-01b7720791ac
begin
	py"""
	import numpy
	def DFT_numpy(x):
		j = complex(0, 1)
		N = len(x)
		return [sum(x[n]*numpy.exp(-2*j*numpy.pi*n*k/N) for n in range(0,N)) for k in range(0,N)]
	"""
	DFT_numpy = py"DFT_numpy"  # copy python function over so it exists on the Julia side
end

# ╔═╡ f4223f65-7c80-46ab-99ab-4736a2909517
md"Create a numpy array from `vsin`:"

# ╔═╡ 8ae5c333-0ab5-4c13-bebc-e5d7e023989d
vsin_numpy = PyObject(vsin)

# ╔═╡ 01b86d4e-d9f2-42cb-9c8c-f52cab1f4d31
t_numpy = @elapsed DFT_numpy(vsin_numpy) 

# ╔═╡ f0a5520a-c48d-424d-a12a-7d80f51a8071
md"""
The Numpy version took $(round(t_numpy, sigdigits=3)) seconds which is $(round(t_numpy/t_python, sigdigits=3))x slower than the regular Python version.
"""

# ╔═╡ 2a5cdcfc-3d72-458d-8a1a-471d7e2196a9
md"""
## Comparison between Python and Julia

### Pros for both:

1. Both Python and Julia have similar syntax and allow support array comprehensions which allow the user to write a complex double-nested `for` loop in one line.

### Pros for Python:

1. Some users may like `0-based` array indexing better if coming from a `C` background.

### Pros for Julia:

1. Some users may like `1-based` array indexing better if coming from an engineering background (e.g. MATLAB, Fortran, R, Maple).

2. The basic math functions/constants like `exp`, `pi`, `im` are available by default.

3. The `exp` function is generic in that it can accept reals, complex numbers and matrices.  In Python, for complex arguments `cmath.exp` is needed while for real arguments `math.exp` is needed (or `numpy.exp` for numpy).  Therefore it is harder for the user to write generic functions built on top of `exp`.  In Julia functions are generic with no loss in performance or additional complexity for the user.

4. Julia is about **$(round(t_python/t_julia3, sigdigits=3))x** faster than Python.  Often an interative loop is the easiest to write an algorithm but if is too slow then users will try to rewrite it in a vectorized form to speed things up.  Performance isn't an issue...until it is and then requires extra time and expertise to use other plug-in libraries.

5. Julia supports broadcasting generically so a function like `exp` can be distributed over a vector `x` like so `exp.(x)`.  This applies to all functions and makes for much easier development and usability.  To take the `exp` of a matrix use `exp(matrix)` (which is not the same as taking the `exp` of each element of the matrix with `exp.(matrix)`.

6. Julia supports many engineering math functions that are not available in other eccosystems. It is designed for engineers doing technical computing.

7. Julia syntax is more like an engineering textbook, supporting syntax like `f(x) = 2x^2 - 5x + 7` for function definitions, `0.0:0.001:0.5` for quick `start:step:stop` range definitions, and `≈` (`\approx<tab>`) for approximately equal.

### Notable differences

1. Julia uses 1-based indexing; Python uses 0-based indexing

2. Julia uses inclusive ranges (`1:10` is 1 to 10); Python uses exclusive ranges (`range(1,11)` is 1 to 10)
"""

# ╔═╡ Cell order:
# ╟─372dc7ad-5d86-4e09-a058-dc3c243bbb82
# ╠═1d5ea4e7-dd3c-4583-a7dd-6bed050ada7a
# ╟─8daa7390-bf0a-11eb-2f09-cb503508dca6
# ╠═ffb8cd4a-21f8-4507-945d-99b6967c9788
# ╟─4a8624aa-2dad-4ce1-8442-32ec454930d6
# ╟─9f4a2c5f-a8b7-4fbf-93a4-2bc129ae7207
# ╠═38c57971-eabc-4ecf-948b-ffbcffdf517d
# ╠═a0fd711d-2c86-4b30-b5b8-b1222c60486e
# ╟─883216a0-a527-4b7a-b266-b6b08e321c6c
# ╟─ac8f9841-8bd3-4913-8773-95923761c969
# ╠═bd9842bb-dcca-4f75-a916-dbfe51a8469b
# ╠═4394b038-a5e5-4054-8bbb-b78d49329f00
# ╟─3e57dd3b-c179-4d1e-aeec-42ef7fc7b9ca
# ╟─c6c22775-9d3a-4cf8-ba84-df3e7ea51b03
# ╟─e22aded1-cf5b-43ae-a923-fedca30f0fec
# ╠═2ed07614-1715-4646-a09c-c7653a0102db
# ╠═7be5b003-e8b4-4795-ad1e-350699b360fc
# ╠═8350eed7-0e88-4055-8060-49bf8931c84e
# ╠═bce85b30-df7f-4697-a9d6-2c2c281f5a4c
# ╠═418d36e3-0d7a-4b8a-aef6-ea0704e72320
# ╟─485c624a-d950-46b5-8ed6-b0a100ae5008
# ╠═cb7a88cb-5f29-4feb-9560-e77013967dcd
# ╠═0f771326-6994-40a3-86c8-7dc49256d634
# ╠═7034ecce-8ac5-45cc-9786-9294b9c99fd9
# ╠═3f6c22d2-aefd-4b3e-bdc5-09b15a02a516
# ╠═2e9845bc-a6ad-4b12-962f-35d1def78e64
# ╟─c35ce815-f1b9-4367-9f81-b2bd51850b12
# ╟─b4b92f64-276f-4bdb-8b68-23406efc774a
# ╠═3be5aa2a-92b2-4dbe-a9ff-bfafa2f6ecb0
# ╠═470dfa78-66c8-40a7-90dd-749667244207
# ╠═74abbc3f-784e-474c-bfdc-1e19718ca6cb
# ╟─bb6d5b1c-9382-43ed-8725-c0ab1b513479
# ╟─22de6bec-5774-40e4-9685-9347de93b21b
# ╠═c14506f3-c3f6-42a5-86f8-93812bd5e789
# ╠═55ed8230-7ab1-41dd-9337-01b7720791ac
# ╠═f4223f65-7c80-46ab-99ab-4736a2909517
# ╠═8ae5c333-0ab5-4c13-bebc-e5d7e023989d
# ╠═01b86d4e-d9f2-42cb-9c8c-f52cab1f4d31
# ╟─f0a5520a-c48d-424d-a12a-7d80f51a8071
# ╟─2a5cdcfc-3d72-458d-8a1a-471d7e2196a9
