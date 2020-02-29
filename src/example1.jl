# An example of the ScalarWaveEquation code for the propagation
# of a gaussian wave in time

using ScalarWaveEquation
using Printf
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

Ttot = 0.3  # total duration
Nt   = 1000  # number of time divisions
Nx   = 400  # number of divisions on position axis

# initial conditions are accepted in the form of FunArr instances
# we create them from 1D Arrays
# the x-range is assumed to be [0, 1]
# both initial conditions must have the same size
init_positions  = get_gaussian(Nx)
init_velocities = zeros(Nx)
ufun  = FunArr{Float64, Float64}(init_positions)
utfun = FunArr{Float64, Float64}(init_velocities)

# use the method of lines to calculate u(x, t) over [0, 1]x[0, Ttot]
# returns a FunArr2D instance with the 2D Array
# also returns the energy at each of the Nt time steps
u2fun, energies = scalarwave_and_energy(ufun, utfun, Nt, Ttot)

# the energy tends to increase over time
# a small energy increase indicates good convergence
@printf("The lowest energy calculated is %.8f\n", minimum(energies))
@printf("The highest energy calculated is %.8f\n", maximum(energies))

# u(x, t) can be plotted as a colour map
xlabel("t")
ylabel("x")
pcolormesh(u2fun.tdat, u2fun.xdat, u2fun.u2dat)
colorbar()
