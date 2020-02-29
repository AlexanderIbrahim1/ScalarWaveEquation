# An example of how the solution for the wave equation
# using a gaussian initial condition converges

using ScalarWaveEquation
using PyPlot

function get_gaussian(N::Int)
	# returns a normalized Gaussian distribution with (mean, std) = (0.5, 0.1)
	mean = 0.5
	std  = 0.1
	
	xdata      = range(0.0, stop = 1.0, length = N)
	exponent   = -0.5 * (xdata.-mean) .* (xdata.-mean) / (std*std)
	norm_coeff = std*sqrt(2*pi)
	
	return exp.(exponent)/norm_coeff
end

function get_gaussian_solution(Nx::Int, Nt::Int, Ttot::T) where {T}
	# returns a FunArr2D instance for the solution to the gaussian
	# initial condition, performed at difference resolutions
	
	# initial conditions are accepted in the form of FunArr instances
	# the x-range is assumed to be [0, 1]
	init_positions  = get_gaussian(Nx)
	init_velocities = zeros(Nx)
	ufun  = FunArr{Float64, Float64}(init_positions)
	utfun = FunArr{Float64, Float64}(init_velocities)
	
	# use the method of lines to calculate u(x, t) over [0, 1]x[0, Ttot]
	u2fun, energies = scalarwave_and_energy(ufun, utfun, Nt, Ttot)
	
	return u2fun
end

# get the solution at different resolutions for the same duration
u2fun1 = get_gaussian_solution(100, 150, 0.2)
u2fun2 = get_gaussian_solution(200, 300, 0.2)
u2fun3 = get_gaussian_solution(300, 450, 0.2)
u2fun4 = get_gaussian_solution(400, 600, 0.2)
u2fun5 = get_gaussian_solution(500, 750, 0.2)
u2fun6 = get_gaussian_solution(600, 900, 0.2)

# get the integrated squared difference between solutions
# for successive resolutions, for a 2000x2000 spacing
area12 = integrate_difference(u2fun1, u2fun2, 2000, 2000)
area23 = integrate_difference(u2fun2, u2fun3, 2000, 2000)
area34 = integrate_difference(u2fun3, u2fun4, 2000, 2000)
area45 = integrate_difference(u2fun4, u2fun5, 2000, 2000)
area56 = integrate_difference(u2fun5, u2fun6, 2000, 2000)

# the differences between increasing resolutions decreases
println("The integrated difference between 100x150 and 200x300 is $area12")
println("The integrated difference between 200x300 and 300x450 is $area23")
println("The integrated difference between 300x450 and 400x600 is $area34")
println("The integrated difference between 400x600 and 500x750 is $area45")
println("The integrated difference between 500x750 and 600x900 is $area56")
