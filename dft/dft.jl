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
		Pkg.PackageSpec(name="PlutoUI", version="0.7"),
		Pkg.PackageSpec(name="BenchmarkTools", version="1"),
    ])
    using LaTeXStrings, PlutoUI, BenchmarkTools
	TableOfContents(title="📚 Table of Contents", indent=true, depth=4, aside=true)
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

If running locally we will first install Python and Conda and a couple Julia packages in a sandbox environment...just hang tight for 10 seconds or so...and about 30 seconds to run the complete notebook the first time.
"""

# ╔═╡ 8daa7390-bf0a-11eb-2f09-cb503508dca6
md"""
# DFT Challenge

## Motivation

This is a challenge to write the Discrete Fourier Tranform in any language to see how well it is suitable to electrical engineering tasks.  Electrical engineers typically deal with data in the time and frequency domain so dealing with real and complex numbers is quite common along with processing lots of data.

There are many ways to write a faster DFT but this challenge is to just use the basic algorithm and to try to write it by hand like it is described in a textbook.  There are many issues that an electrical engineer will run into where they will have to write their own solution to and not depend on code written by someone else.

## DFT description from a textbook

$$H_k[k] = \sum_{n=0}^{N-1}x[n] e^{\frac{-j2\pi k n}{N}} \text{, where } k = 0 \ldots N-1$$

So for N samples of vector `x` there will be `N` results in the DFT vector `Hₖ`.  For the first results where `k=0` we have:

$$H_0[0] = \sum_{n=0}^{N-1}x[n] e^{0} \text{,    where } k = 0$$

And the second when `k=1` we have:

$$H_1[1] = \sum_{n=0}^{N-1}x[n] e^{\frac{-j2\pi 1 n}{N}} \text{,    where } k = 1$$

And so on until `k = N-1`.

## Julia
"""

# ╔═╡ f9295ff5-5a1c-4217-a455-bb12d37880fc
md"""
The first implementation is in [Julia](https://julialang.org).  A somewhat newer language built to be higher-level than Python and as fast as C or Fortran, yet easy for engineers to code in.
"""

# ╔═╡ ffb8cd4a-21f8-4507-945d-99b6967c9788
function DFT(x)
	N = length(x)
	Hₖ = [sum(x[n+1]*ℯ^(-im*2π*n*k/N) for n in 0:N-1) for k in 0:N-1]
end

# ╔═╡ 4a8624aa-2dad-4ce1-8442-32ec454930d6
md"""
!!! note
	Julia uses `im` for the imaginary number and by default uses `1-based` indexing.  

	To write `π` type `\pi<tab>`, for Euler's number type `\euler<tab>` or use `exp()`, and for `Hₖ` type `H\_k<tab>`."""

# ╔═╡ 9f4a2c5f-a8b7-4fbf-93a4-2bc129ae7207
md"### Create some test data

The test data is a simle 2000 point sinusoid...nothing special."

# ╔═╡ 38c57971-eabc-4ecf-948b-ffbcffdf517d
t = range(0, stop=1, length=2000)

# ╔═╡ a0fd711d-2c86-4b30-b5b8-b1222c60486e
vsin = sin.(2pi .* t)

# ╔═╡ 883216a0-a527-4b7a-b266-b6b08e321c6c
md"""
!!! note
	Use "dot" broadcasting in Julia where any function can be executed element-by-element of the arguments inside.  Above `sin.` and `.*` operate element by element so it is like having a compact `for` loop to take the `sin` of each element of the vector passed to it.
"""

# ╔═╡ ac8f9841-8bd3-4913-8773-95923761c969
md"### Run Julia and record the execution time"

# ╔═╡ bd9842bb-dcca-4f75-a916-dbfe51a8469b
dft1 = DFT(vsin)

# ╔═╡ adc93d58-0dce-4c6e-b100-e93ca481f53c
md"""
!!! note
	For measuring the run time accurately we will use the `@belapsed` macro from the `BenchmarkTools` package which re-runs the code over and over to measure the elapsed time with reduced noise.
"""

# ╔═╡ 4394b038-a5e5-4054-8bbb-b78d49329f00
t_julia = @belapsed DFT(vsin)

# ╔═╡ 3e57dd3b-c179-4d1e-aeec-42ef7fc7b9ca
md"""
The basic Julia version finishes in **$(round(t_julia, sigdigits=4))** seconds
"""

# ╔═╡ c6c22775-9d3a-4cf8-ba84-df3e7ea51b03
md"""
!!! warning
	The **`2000`** point DFT must execute the loop `2000*2000` or **4 million** times.  Each loop has about 6 multiplies, 1 divide, and 1 `exp` of an imaginary number.  `exp` is a large function so it is hard to know how many operations are needed.  Fortunately, Julia can tell us by using the [GFlops.jl](https://github.com/triscale-innov/GFlops.jl) package.  It reports about 90 FLOPs per iteration so that is about **400 million** FLOPs for the complete DFT.
"""

# ╔═╡ 2ed07614-1715-4646-a09c-c7653a0102db
function DFT2(x)
	N = length(x)
	Hₖ = [sum(x[n+1]*cispi(-2n*k/N) for n in 0:N-1) for k in 0:N-1]
end

# ╔═╡ db8c90c4-dba7-44b1-a2d5-4ee5756667c2
md"""
Let's check that the improved function returns the same result:
"""

# ╔═╡ 485c624a-d950-46b5-8ed6-b0a100ae5008
md"""
### Improved Julia version (with 0-based indexing)

The first index of a DFT is typically `0` for DC so users would like to have the first index be `0`.  Also, some users would like to use 0-based indexing and Julia can support this through a custom defined type.  The `OffsetArrays` package defines an array with user-defined indexing.  Let's use this to see how it works.
"""

# ╔═╡ 0f771326-6994-40a3-86c8-7dc49256d634
function DFT2(x::OffsetVector)
	N = length(x)
	Hₖ = [sum(x[n]*cispi(-2n*k/N) for n in 0:N-1) for k in 0:N-1]
	OffsetVector(Hₖ, 0:N-1)
end

# ╔═╡ 7be5b003-e8b4-4795-ad1e-350699b360fc
dft2 = DFT2(vsin)

# ╔═╡ bce85b30-df7f-4697-a9d6-2c2c281f5a4c
dft1 ≈ dft2

# ╔═╡ 8350eed7-0e88-4055-8060-49bf8931c84e
t_julia2 = @belapsed DFT2(vsin)

# ╔═╡ e22aded1-cf5b-43ae-a923-fedca30f0fec
md"""
### Improved Julia version (`cos` and `sin` = $(round(t_julia/t_julia2, sigdigits=2))x faster)

Julia has optimized functions that are common in engineering but are not found in most programming languages.

In this case we can use the `cispi` function which is the same as:

$$\mathrm{cispi}(x) = \cos(\pi x) + \texttt{i} \sin(\pi x) = e^{\texttt{i} \pi x}$$

where `i` is imaginary.  `cispi(x)` is a faster version of `exp(im*pi*x)` and also more accurate for large `x` since it doesn't depend on an inaccurate number for π 	(see [Wikipedia](https://en.wikipedia.org/wiki/Cis_(mathematics)) for more info).

So let's rewrite the DFT with `cispi`:
"""

# ╔═╡ 418d36e3-0d7a-4b8a-aef6-ea0704e72320
md"""
!!! note
	Above the `≈` operater checks for approximately equal (within floating point round-off error).  It can be typed with `\approx<tab>`.

So the `cispi` version took $(round(t_julia2, sigdigits=3)) seconds.  It gets the same result and is $(round(t_julia/t_julia2, sigdigits=2))x faster.
"""

# ╔═╡ fd7807ff-3a89-499e-937e-7550979145bd
md"""
!!! note
	This `DFT2` function definition specifies that the input must be an `OffsetVector`.  This allows users to write generic functions that can take different types of inputs yet run different code.  So by default `DFT2` will assume it is a regular vector (or any type since it does no type checking) but if an `OffsetVector` is passed to `DFT2` then the above method will be called instead.  Note how above the code Julia, it returns `DFT2 (generic function with 2 methods)`.


Now let's convert `vsin` into a 0-based index `OffsetVector`:
"""

# ╔═╡ 7034ecce-8ac5-45cc-9786-9294b9c99fd9
vsin_0 = OffsetVector(vsin, 0:length(vsin)-1)

# ╔═╡ 3f6c22d2-aefd-4b3e-bdc5-09b15a02a516
dft3 = DFT2(vsin_0)

# ╔═╡ 2e9845bc-a6ad-4b12-962f-35d1def78e64
t_julia3 = @belapsed DFT2(vsin_0)

# ╔═╡ c35ce815-f1b9-4367-9f81-b2bd51850b12
md"""The 0-based index version took $(round(t_julia3, sigdigits=3)) seconds which is $(round(t_julia2/t_julia3, sigdigits=2))x faster than the previous version."""

# ╔═╡ 40a76500-fb57-408f-baf2-032be784e7f0
function DFT2_threads(x)
	N = length(x)
    Hₖ = Vector{Complex{Float64}}(undef, N) # pre-allocate for thread-safety
    Threads.@threads for k in 0:N-1
        hₖ = 0.0 + im*0.0
        for n in 0:N-1
            hₖ += x[n+1]*cispi(-2n*k/N)
        end
        Hₖ[k+1] = hₖ
    end
    Hₖ
end

# ╔═╡ ba63aff4-f0f0-4f0d-acc7-7ce7d3ab60a4
dft2_threads = DFT2_threads(vsin)

# ╔═╡ ee637fb0-dba2-47fb-8cc2-117ad23cd5a6
t_julia2_threads = @belapsed DFT2_threads(vsin)

# ╔═╡ 4a4e3b04-6e00-46fb-a9b8-8a3ca06d9cf4
md"""
### Improved multi-threaded Julia version ($(Threads.nthreads()) CPUs = $(round(t_julia2/t_julia2_threads, sigdigits=2))x)

Julia is a bit unique in that it has built-in support for multi-threading.  To do so Julia must be started with `julia --threads N` and then put `Threads.@threads` before the `for` loop.  We will write it in a similar fashion to the final form for the other languages for comparision: 
"""

# ╔═╡ bb67d192-27e9-49e7-9853-202f498e6f46
dft2_threads ≈ dft2

# ╔═╡ c662aae3-fcd5-48dc-9e7c-fbbcbaf65e51
md"""
With using $(Threads.nthreads()) threads it finished in $(round(t_julia2_threads, sigdigits=3)) seconds which is $(round(t_julia2/t_julia2_threads, sigdigits=2))x faster than the `DFT2` Julia version without threads.

So with Julia, not only can you get really fast and easy to write code but it is also easy to make it multi-threaded and achieve close to linear CPU scaling with minimal effort.
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

# ╔═╡ 13196703-d411-4ad8-aa92-694812ab5508
md"""
!!! note

	It took a while to figure out how to handle complex numbers and `pi`.  The `math.exp` function doesn't work with complex numbers so `cmath.exp` had to be used instead.  It is clear that Python doesn't have good support for generics.
"""

# ╔═╡ 74abbc3f-784e-474c-bfdc-1e19718ca6cb
t_python = @elapsed dftpy = DFT_py(vsin) # this takes a long time to run (eg > 10 seconds)


# ╔═╡ bb6d5b1c-9382-43ed-8725-c0ab1b513479
md"""
The Python version takse too long to run so `@belapsed` isn't used to run it over and over.  It took $(round(t_python, sigdigits=3)) seconds which is $(round(Int, t_python/t_julia2_threads))x slower than the multi-threaded Julia version.  Let's check if they are equal:
"""


# ╔═╡ 1f084504-4f3a-4eb9-a4b9-b12b3c154497
dftpy ≈ dft2_threads

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
md"Create a numpy array from `vsin` (PyCall by default converts Julia arrays into Numpy arrays):"

# ╔═╡ 8ae5c333-0ab5-4c13-bebc-e5d7e023989d
vsin_numpy = PyObject(vsin)

# ╔═╡ 01b86d4e-d9f2-42cb-9c8c-f52cab1f4d31
t_numpy = @elapsed dftnumpy = DFT_numpy(vsin_numpy) 

# ╔═╡ 22de6bec-5774-40e4-9685-9347de93b21b
md""" ### Improved Numpy version ($(round(Int, t_numpy/t_julia2_threads))x slower)

We will write the same algorithm using Numpy.
"""

# ╔═╡ f0a5520a-c48d-424d-a12a-7d80f51a8071
md"""
The Numpy version took $(round(t_numpy, sigdigits=3)) seconds which is $(round(t_numpy/t_python, sigdigits=2))x slower than the regular Python version.  So the speed-up with using Numpy is only for vectorized operations.

Let's check if they are equal:
"""

# ╔═╡ 4bdb81af-7520-46ae-a151-90446ac4be61
dftnumpy ≈ dftpy

# ╔═╡ b7e7511a-0d0a-4b63-afdc-0c8de4352ff8
begin
	py"""
	import math
	def DFT_sincos_py(x):
		N = len(x)
		Hk = []
		for k in range(0,N):
			sum_real = 0.0
			sum_imag = 0.0
			for n in range(0,N):
				sum_real += x[n]*math.cos(-2*math.pi*n*k/N)
				sum_imag += x[n]*math.sin(-2*math.pi*n*k/N)
			Hk.append(complex(sum_real, sum_imag))
		return Hk
	"""
	DFT_sincos_py = py"DFT_sincos_py"  # copy python function over so it exists on the Julia side
	
end

# ╔═╡ cac67f07-f3ca-49de-b3df-d4f7c8b71a4c
t_python2 = @elapsed dft2py = DFT_sincos_py(vsin)

# ╔═╡ b4b92f64-276f-4bdb-8b68-23406efc774a
md"""
## Python ($(round(Int, t_python2/t_julia2_threads))x slower)

### Python `exp` version ($(round(Int, t_python/t_julia2_threads))x slower)

We will use the `PyCall` package to call Python from Julia

"""

# ╔═╡ 869b8945-3cc9-457f-a3cf-5cca441db802
md"""
### Improved Python (`cos` and `sin` = $(round(Int, t_python2/t_julia2_threads))x slower)

Because Julia can use `cis` instead of `exp`, let's create a similar Python version:

"""

# ╔═╡ 498a0aff-8879-447d-85f2-35dc10e5aae7
dft2py ≈ dft1

# ╔═╡ 5f078b8e-4068-44c8-932f-db769cbee742
md"""
This Python version finished in $(round(t_python2, sigdigits=2)) seconds which is $(round(t_python/t_python2, sigdigits=2))x faster than the `exp` version of Python and $(round(Int, t_python2/t_julia2_threads))x slower than comparable (multi-threaded) Julia version.
"""

# ╔═╡ 2a5cdcfc-3d72-458d-8a1a-471d7e2196a9
md"""
### Comparison between Python and Julia

#### Pros for both

1. Both Python and Julia have similar syntax and allow support array comprehensions which allow the user to write a complex double-nested `for` loop in one line.

#### Pros for Python

1. Some users may like `0-based` array indexing better if coming from a C, TCL, SKILL, or Perl background.

#### Pros for Julia

1. Some users may like `1-based` array indexing better if coming from an engineering background (e.g. MATLAB, Fortran, R, Maple).

2. The basic math functions/constants like `exp`, `pi`, `im` are available by default.

3. The `exp` function is generic in that it can accept reals, complex numbers and matrices.  In Python, for complex arguments `cmath.exp` is needed while for real arguments `math.exp` is needed (or `numpy.exp` for numpy).  Therefore it is harder for the user to write generic functions built on top of `exp`.  In Julia functions are generic with no loss in performance or additional complexity for the user.

4. For the simple DFT, Julia is **$(round(Int, t_python2/t_julia2_threads))x** faster than pure Python (both using `cos` and `sin`).  Often an interative loop is the easiest to write an algorithm but if the code is too slow then Python users will try to rewrite the algorithm in a vectorized form to speed things up with Numpy.  Performance isn't an issue...until it is...and then requires extra time and expertise to use other plug-in libraries.

4. The multi-threaded version of Julia with $(Threads.nthreads()) threads was $(round(t_julia2/t_julia2_threads, sigdigits=2))x faster than without threads.  It was very painless to do this.  Python has a dozen or so libraries which support a subset of the Python language (with different symantics) to speed things up in various ways but it is a big sacrifice in usablity and generality.  Looking at numerical computing code in Python they are usually a mix of C, Fortran and Cython.  Most engineers are not going to want to write a build system to do this and learn 1 or 2 more languages.

5. Julia supports broadcasting generically so a function like `exp` can be distributed over each element of vector `x` like so `exp.(x)`.  This applies to all functions and makes for much easier development and usability.  To take the `exp` of a matrix use `exp(matrix)` (which is not the same as taking the `exp` of each element of the matrix with `exp.(matrix)`.

6. Julia supports many engineering math functions that are not available in other ecosystems. It is designed for engineers doing technical computing.

7. Julia syntax is more like an engineering textbook, supporting syntax like `f(x) = 2x^2 - 5x + 7` for function definitions, `0.0:0.001:0.5` for quick `start:step:stop` range definitions, and `≈` (`\approx<tab>`) for approximately equal.

8. Julia has great Python interopability and can `import` Python packages and call them with ease.

#### Notable differences

1. Julia uses 1-based indexing; Python uses 0-based indexing

2. Julia uses inclusive ranges (`1:10` is 1 to 10); Python uses exclusive ranges (`range(1,11)` is 1 to 10)
"""

# ╔═╡ 527585d1-e9c7-42ad-aa1c-e443772f522c
t_tcl = 3.858

# ╔═╡ e62a71a3-7c54-4390-9263-5e6e55e70718
md"""
## TCL ($(round(Int, t_tcl/t_julia2_threads))x slower)

TCL is very popular in EDA so let's write the DFT function in TCL:

```tcl
proc dft {x} {
    # TCL doesn't support complex numbers so will keep them separate and use Euler's formula
    set N [llength $x]
    set pi [expr "acos(-1)"]
    set Hk {}
    for {set k 0} {$k < $N} {incr k} {
        set sum_real 0.0
        set sum_imag 0.0
        for {set n 0} {$n < $N} {incr n} {
            set xn [lindex $x $n]
            set real [expr {$xn*cos(-2*$pi*$n*$k/double($N))}]
            set imag [expr {$xn*sin(-2*$pi*$n*$k/double($N))}]
            set sum_real [expr {$sum_real + $real}]
            set sum_imag [expr {$sum_imag + $imag}]
        }
        lappend Hk $sum_real $sum_imag
    }
    return $Hk
}

proc listFromFile {filename} {
    set f [open $filename r]
    set data [split [string trim [read $f]]]
    close $f
    return $data
}

set x [listFromFile data.csv]
set N [llength $x]
puts "Read in $N lines for vsin"

puts "Calculating DFT"
set t0 [clock milliseconds]
set Hk [dft $x]
set t1 [clock milliseconds]
set dt [expr {($t1-$t0)/1000.0}]
puts "$N-point DFT took $dt seconds."

proc print_Hk {Hk} {
    foreach {real imag} $Hk {
        puts "$real $imag"
    }
}
```

And let's run it (on my local computer):
```bash
% tclsh dft.tcl                                                                       Read in 2000 lines for vsin
Calculating DFT
2000-point DFT took 3.858 seconds.
```


"""

# ╔═╡ f1dd33f2-744f-4288-92ab-7ff54108dd05
md"""
### Comparision between TCL, Julia and Python

1. It was very difficult to write the DFT in TCL.  TCL doesn't support complex numbers and the code in TCL is unlike any math textbook.  Instead of a 2 line function it is 14 lines and is much less usable.

2. In this case at least, TCL is a bit faster than pure Python.  But being that TCL doesn't support complex numbers it wouldn't make sense to use TCL.

2. I made many errors trying to get the TCL version to work and couldn't understand the error messages.  It turned out I had missed putting the `$` in front of a variable name and it gave weird error messages that referenced lines far away from the actual issue.  

3. TCL does integer division so `1/4` equals `0` (not `0.25`).  I had to make sure to do `/double($N)` to get regular division.  These are the sorts of things that can easily cause code to break (say a user passed in an integer) in non-obvious ways and lead to long debugging sessions.  I've seen this bug in production code where the Verilog-A compiler (written in C which also does integer division) was doing `V/R` but the resistance was an integer and the programmer forgot to convert it to a real and the bug was in the field many years before it could be tracked down. 

4. Basic math constants like `pi` are not available.  Had to use `acos(-1)` to get the value of `pi`.

5. Without complex number support using the output of the DFT is going to be very painful as functions like `abs` would need to be rewritten for complex numbers.  There are TCL libraries out there for complex numbers but they are typically not shipped/used with EDA software.

6. TCL doesn't have math as part of its syntax.  Math syntax is handled specifically by the `expr` function.  I don't know of a way to use my newly created `dft` function along with the `expr` function so I can use it in other math expressions.  User defined math functions in TCL are not composable.

7. The run time of $t_tcl seconds was $(round(Int, t_tcl/t_julia2_threads))x slower than Julia and $(round(t_python2/t_tcl, sigdigits=2))x faster than Python.  A bit surprising to that TCL was faster than Python.  But TCL doesn't support complex numbers so it isn't a fair comparision as someone who needed complex numbers wouldn't want to use TCL.
"""

# ╔═╡ c20504f8-56a6-47f6-9bd8-29c127c860f6
md"""
## Synopsys ACE/TCL

Synopsys has a math calculator language (called ACE/TCL) that is accessed from TCL using the `sx_equation` command.  It is only for expressions (no control flow or looping) so the DFT below is a combination of TCL and `sx_equation` but `sx_equation` is used as much as possible where it can be used.


```tcl
proc DFT {x} {
    sx_equation "N=[llength $x]"
    sx_equation "pi=acos(-1)"
    set Hk {}
    for {set k 0} {$k < $N} {incr k} {
        sx_equation "sum_real=0.0"
        sx_equation "sum_imag=0.0"
        for {set n 0} {$n < $N} {incr n} {
            sx_equation "xn=[lindex $x $n]}
            sx_equation "n=$n"
            sx_equation "k=$k"
            sx_equation "real=xn*cos(-2*pi*n*k/N)"
            sx_equation "imag=xn*sin(-2*pi*n*k/N)"
            sx_equation "sum_real = sum_real + real"
            sx_equation "sum_imag = sum_imag + imag"
        }
        lappend Hk $sum_real $sum_imag
    }
    return $Hk
}

proc listFromFile {filename} {
    set f [open $filename r]
    set data [split [string trim [read $f]]]
    close $f
    return $data
}

set x [listFromFile data.csv]
set N [llength $x]
puts "Read in $N lines for vsin"

puts "Calculating DFT"
set t0 [clock milliseconds]
set Hk [dft $x]
set t1 [clock milliseconds]
set dt [expr {($t1-$t0)/1000.0}]
puts "$N-point DFT took $dt seconds."
```

The above assumes the data passed to the `DFT` proc is a TCL list.  If it was a waveform type then accessing the values from the waveform would be more complex as well as if the result was returned as a complex waveform (I suspect it would be slower too).  The above code is completely untested and hasn't been run but if someone can run it and send the results, that would be appreciated -- my guess is it will be at least an order of magnitude slower than anything else.

### Comparision between ACE/TCL and Julia

1. The speed will be prohibitively slow.  The calculation may take longer than running the whole simulation.

2. The ACE/TCL version is very difficult to write because the user must know how to write both TCL and ACE/TCL.  `sx_equation` is needed because TCL's math processing is too limited but then `sx_equation` has no `for` loops or control flow so the user must also know TCL.  

3. The new user-defined `DFT` function is not available in `sx_equation` because `sx_equation` doesn't know about it.  There are other ACE/TCL functions for registering a user-defined function into `sx_equation` but the are buggy and don't work realiably.

4. The error messages from ACE/TCL are typically unhelpful, often referencing the line of the top-most `proc` and not enough information to do debugging.

5. If the user wanted to write the formula in `exp(-2*pi*j*n*k/N)` form does the `exp` in `sx_equation` handle complex numbers?
"""

# ╔═╡ 3e48d6b3-3bf6-4d2e-a6d1-8ad3db80fc1a
md"""
## Cadence SKILL

!!! note
	If you have access to SKILL and create an implementation, please let me know (see below).
"""

# ╔═╡ c12c124c-62e8-4a57-a46a-0b1a44f14571
md"""
## MATLAB

!!! note

	If you have access to MATLAB and create an implementation, please let me know (see below).
"""

# ╔═╡ 08bc918d-3c42-4d34-b4ea-7927fea3a2dc
begin
	slow(x) = string(round(t_python2/x, sigdigits=2), "x")
	sec(x) = round(x, sigdigits=3)
md"""
# Summary

## Performance
	
To summarize these are the run times (using the Python `sin` `cos` version as a baseline):

| Language | Implementation | Run time (s)             |  Speed-up (x)  |
|:---------|:---------------|:------------------------:|:----------:|
| Julia  |`cispi` 4CPU | $(sec(t_julia2_threads)) | $(slow(t_julia2_threads)) |
| Julia  |`cispi` (1CPU) | $(sec(t_julia2)) | $(slow(t_julia2)) |
| Julia  |`cispi` (0-based) | $(sec(t_julia3)) | $(slow(t_julia3)) |
| Python | `sin` `cos` | $(sec(t_python2)) | $(slow(t_python2)) |
| TCL    | `sin` `cos` | $(sec(t_tcl)) | $(slow(t_tcl)) |
| Julia  |`exp` | $(sec(t_julia)) | $(slow(t_julia)) |
| Python |`exp` | $(sec(t_python)) | $(slow(t_python)) |
| Python | `exp` Numpy | $(sec(t_numpy)) | $(slow(t_numpy)) |
| TCL    | `exp`      | NA | NA |
	
## Ease of use
	
Looking at comparison of other factors the ratings are from `0` (not available) to `10` (great support) for how closely the language resembles an engineering textbook (as a proxy for ease of use).
	
| Ease of use aspect | Julia | Python | TCL |
|:-------|:-----:|:------:|:---:|
|Math multiplication syntax (`3x + 5`) | 10 | 5 | 3 |
|One line `for` loops | 8 | 8 | 0 |
|Element-by-element operations (`abs.(vec)`) | 7 | 2 | 0 |
|Imaginary numbers (`3 + j4`) | 9 | 7 | 0 |
|Constants like `pi`, `π`, `ℯ`       | 10 | 7  | 1 |
|Math notation for names (`Hₖ`, `ĝ`, etc.)  | 10 | ? | ? |
|Floating point approximately equal (`≈`) | 10 | 0 | 0 |
|Real division (`1/4 = 0.25`) | 10 | 10 | 0 |
|Generic functions (`exp` of complex, reals)     | 10 | 2  | 0 |
|User extendable generic functions | 10 | 2  | 0 |
|Error message quality | 8 | 9 | 2 |
	
"""	

end

# ╔═╡ 0532bbb3-428e-4d32-876e-af92c1c7bb01
md"""
# Getting in contact

If you have any comments or suggestions, please let me know at `hertz` at `juliacomputing.com`.  I'd really appreciate examples from other programming languages common in EDA.  


If you are interested in using Julia in your company then reach out for technical support, professional services or engineering/EDA design software.  See more info at [juliacomputing.com](https://juliacomputing.com/).  We enable engineers to be more productive.

# Appendix

1. This site was written in [Pluto](https://github.com/fonsp/Pluto.jl) which is an interactive Julia environment in a web browser.

2. MIT has many good university level courses using Julia and Pluto on their [Computational Thinking](https://computationalthinking.mit.edu/Spring21/) course.  It is really a game changer for someone learning technical computing in that it gives students an interactive notebook where they can interact with the mathematics and get a good understanding of how the algorithms work.  



"""

# ╔═╡ Cell order:
# ╟─372dc7ad-5d86-4e09-a058-dc3c243bbb82
# ╠═1d5ea4e7-dd3c-4583-a7dd-6bed050ada7a
# ╟─8daa7390-bf0a-11eb-2f09-cb503508dca6
# ╟─f9295ff5-5a1c-4217-a455-bb12d37880fc
# ╠═ffb8cd4a-21f8-4507-945d-99b6967c9788
# ╟─4a8624aa-2dad-4ce1-8442-32ec454930d6
# ╟─9f4a2c5f-a8b7-4fbf-93a4-2bc129ae7207
# ╠═38c57971-eabc-4ecf-948b-ffbcffdf517d
# ╠═a0fd711d-2c86-4b30-b5b8-b1222c60486e
# ╟─883216a0-a527-4b7a-b266-b6b08e321c6c
# ╟─ac8f9841-8bd3-4913-8773-95923761c969
# ╠═bd9842bb-dcca-4f75-a916-dbfe51a8469b
# ╟─adc93d58-0dce-4c6e-b100-e93ca481f53c
# ╠═4394b038-a5e5-4054-8bbb-b78d49329f00
# ╟─3e57dd3b-c179-4d1e-aeec-42ef7fc7b9ca
# ╟─c6c22775-9d3a-4cf8-ba84-df3e7ea51b03
# ╟─e22aded1-cf5b-43ae-a923-fedca30f0fec
# ╠═2ed07614-1715-4646-a09c-c7653a0102db
# ╠═7be5b003-e8b4-4795-ad1e-350699b360fc
# ╠═8350eed7-0e88-4055-8060-49bf8931c84e
# ╟─db8c90c4-dba7-44b1-a2d5-4ee5756667c2
# ╠═bce85b30-df7f-4697-a9d6-2c2c281f5a4c
# ╟─418d36e3-0d7a-4b8a-aef6-ea0704e72320
# ╟─485c624a-d950-46b5-8ed6-b0a100ae5008
# ╠═cb7a88cb-5f29-4feb-9560-e77013967dcd
# ╠═0f771326-6994-40a3-86c8-7dc49256d634
# ╟─fd7807ff-3a89-499e-937e-7550979145bd
# ╠═7034ecce-8ac5-45cc-9786-9294b9c99fd9
# ╠═3f6c22d2-aefd-4b3e-bdc5-09b15a02a516
# ╠═2e9845bc-a6ad-4b12-962f-35d1def78e64
# ╟─c35ce815-f1b9-4367-9f81-b2bd51850b12
# ╟─4a4e3b04-6e00-46fb-a9b8-8a3ca06d9cf4
# ╠═40a76500-fb57-408f-baf2-032be784e7f0
# ╠═ba63aff4-f0f0-4f0d-acc7-7ce7d3ab60a4
# ╠═ee637fb0-dba2-47fb-8cc2-117ad23cd5a6
# ╠═bb67d192-27e9-49e7-9853-202f498e6f46
# ╟─c662aae3-fcd5-48dc-9e7c-fbbcbaf65e51
# ╟─b4b92f64-276f-4bdb-8b68-23406efc774a
# ╠═3be5aa2a-92b2-4dbe-a9ff-bfafa2f6ecb0
# ╠═470dfa78-66c8-40a7-90dd-749667244207
# ╟─13196703-d411-4ad8-aa92-694812ab5508
# ╠═74abbc3f-784e-474c-bfdc-1e19718ca6cb
# ╟─bb6d5b1c-9382-43ed-8725-c0ab1b513479
# ╠═1f084504-4f3a-4eb9-a4b9-b12b3c154497
# ╟─22de6bec-5774-40e4-9685-9347de93b21b
# ╠═c14506f3-c3f6-42a5-86f8-93812bd5e789
# ╠═55ed8230-7ab1-41dd-9337-01b7720791ac
# ╟─f4223f65-7c80-46ab-99ab-4736a2909517
# ╠═8ae5c333-0ab5-4c13-bebc-e5d7e023989d
# ╠═01b86d4e-d9f2-42cb-9c8c-f52cab1f4d31
# ╟─f0a5520a-c48d-424d-a12a-7d80f51a8071
# ╠═4bdb81af-7520-46ae-a151-90446ac4be61
# ╟─869b8945-3cc9-457f-a3cf-5cca441db802
# ╠═b7e7511a-0d0a-4b63-afdc-0c8de4352ff8
# ╠═cac67f07-f3ca-49de-b3df-d4f7c8b71a4c
# ╠═498a0aff-8879-447d-85f2-35dc10e5aae7
# ╟─5f078b8e-4068-44c8-932f-db769cbee742
# ╟─2a5cdcfc-3d72-458d-8a1a-471d7e2196a9
# ╟─e62a71a3-7c54-4390-9263-5e6e55e70718
# ╠═527585d1-e9c7-42ad-aa1c-e443772f522c
# ╟─f1dd33f2-744f-4288-92ab-7ff54108dd05
# ╟─c20504f8-56a6-47f6-9bd8-29c127c860f6
# ╟─3e48d6b3-3bf6-4d2e-a6d1-8ad3db80fc1a
# ╟─c12c124c-62e8-4a57-a46a-0b1a44f14571
# ╟─08bc918d-3c42-4d34-b4ea-7927fea3a2dc
# ╟─0532bbb3-428e-4d32-876e-af92c1c7bb01
