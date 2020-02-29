module ScalarWaveEquation

export FunArr
export FunArr2D
export scalarwave_and_energy
export integrate_difference

mutable struct FunArr{T, U}
	xmin::T            # x-range minimum
	xmax::T            # x-range maximum
	dx::T              # step size of z-data
	ydat::Vector{U}    # range data
	
	function FunArr{T, U}(ydat_) where {T, U}
		size = length(ydat_)
		if size < 3
			error("The y-data Vector must have a length 3 or greater.")
		end
		
		xmin_ = 0.0
		xmax_ = 1.0
		dx_   = (xmax_ - xmin_)/(size - 1)
		new(xmin_, xmax_, dx_, ydat_)
	end
end

mutable struct FunArr2D{T, U, V}
	xmin::T
	xmax::T
	dx::T
	Nx::Int
	tmin::U
	tmax::U
	dt::U
	Nt::Int
	xdat::Array{T}
	tdat::Array{U}
	u2dat::Array{V, 2}
	
	function FunArr2D{T, U, V}(u2dat_, t_duration) where {T, U, V}
		# test for appropriate sizes
		Nx_ = size(u2dat_, 1)
		if Nx_ < 3
			error("The first dimension of the array must have 3 or more elements.")
		end
		Nt_ = size(u2dat_, 2)
		if Nt_ < 3
			error("The second dimension of the array must have 3 or more elements.")
		end
		
		# the duration of elapsed time must be positive
		if t_duration <= 0.0
			error("The duration of time that passes must be positive.")
		end
		
		xmin_ = 0.0
		xmax_ = 1.0
		tmin_ = 0.0
		tmax_ = t_duration
		dx_ = (xmax_ - xmin_)/(Nx_ - 1)
		dt_ = (tmax_ - tmin_)/(Nt_ - 1)
		xdat_ = range(0.0, stop = 1.0, length = Nx_)
		tdat_ = range(0.0, stop = t_duration, length = Nt_)
		new(xmin_, xmax_, dx_, Nx_, tmin_, tmax_, dt_, Nt_, xdat_, tdat_, u2dat_)
	end
end

function first_deriv(f::FunArr{T, U})::FunArr{T, U} where {T, U}
	N  = length(f.ydat)
	dx = U(f.dx)
	dydat = similar(f.ydat)
	
	# take second derivative of end-points
	dydat[1] = (f.ydat[2] - f.ydat[1]  ) / dx
	dydat[N] = (f.ydat[N] - f.ydat[N-1]) / dx
	
	# take the second derivative of all points in between
	for i in 2:(N-1)
		dydat[i] = (f.ydat[i+1] - f.ydat[i-1]) / (2*dx)
	end
	
	FunArr{T, U}(dydat)
end

function second_deriv(f::FunArr{T, U})::FunArr{T, U} where {T, U}
	N  = length(f.ydat)
	dx = U(f.dx)
	d2ydat = similar(f.ydat)
	
	# take second derivative of end-points
	d2ydat[1] = (f.ydat[3] - 2.0*f.ydat[2] + f.ydat[1]) / dx^2
	d2ydat[N] = (f.ydat[N] - 2.0*f.ydat[N-1] + f.ydat[N-2]) / dx^2
	
	# take the second derivative of all points in between
	for i in 2:(N-1)
		d2ydat[i] = (f.ydat[i+1] - 2*f.ydat[i] + f.ydat[i-1]) / dx^2
	end
	
	FunArr{T, U}(d2ydat)
end

function forward_euler2_disc(i0::Int, yfun::FunArr{T, U}, vfun::FunArr{T, U}, afun::FunArr{T, U}, h::T) where {T, U}
	# i    : index value between 1 and N inclusive
	# yfun : the discretized position function
	# vfun : the discretized velocity function
	# afun : the discretized acceleration function
	y0 = yfun.ydat[i0]
	v0 = vfun.ydat[i0]
	a0 = afun.ydat[i0]
	dt = U(h)
	v1 = v0 + a0*dt
	y1 = y0 + v0*dt + 0.5*a0*dt*dt
	
	return y1, v1
end

function get_energy(utfun::FunArr{T, U}, uxfun::FunArr{T, U}) where {T, U}
	etotal = 0.0
	N = length(utfun.ydat)
	h = utfun.dx
	
	for i in 1:N
		etotal += utfun.ydat[i]^2 + uxfun.ydat[i]^2
	end
	
	etotal *= 0.5*h
	
	return etotal
end

function nearest_lower_index(N::I, xmin::T, xmax::T, x::T) where {T, I<:Int}
	# takes care of the rand() = 0.0 corner case
	if x == 0.0
		x += eps(1.0)
	end
	@assert xmax > x > xmin
	ratio = (x - xmin)/(xmax - xmin)
	lower_index = ceil(typeof(N), (N-1)*ratio)
	return lower_index
end

function bilinear_interp_regular(u2fun::FunArr2D{T, U, V}, x::T, t::U) where {T, U, V}
	# performs bilinear interpolation on a function discretized on a regular grid
	# find the x's and y's around the given (x, y) point
	# they will be (x1, y1), (x1, y2), (x2, y1), (x2, y2)
	
	xsize = u2fun.Nx; xmin = u2fun.xmin; xmax = u2fun.xmax; dx = u2fun.dx
	tsize = u2fun.Nt; tmin = u2fun.tmin; tmax = u2fun.tmax; dt = u2fun.dt
	
	Nx1 = nearest_lower_index(xsize, xmin, xmax, x)
	Nt1 = nearest_lower_index(tsize, tmin, tmax, t)
	
	x1 = xmin + (Nx1-1)*dx
	x2 = xmin + Nx1*dx
	t1 = tmin + (Nt1-1)*dt
	t2 = tmin + Nt1*dt
	
	u11 = u2fun.u2dat[Nx1  , Nt1  ]
	u12 = u2fun.u2dat[Nx1  , Nt1+1]
	u21 = u2fun.u2dat[Nx1+1, Nt1  ]
	u22 = u2fun.u2dat[Nx1+1, Nt1+1]
	
	coeff = 1.0/((x2 - x1)*(t2 - t1))
	xV    = [x2 - x, x - x1]
	tV    = [t2 - t, t - t1]
	uM    = [u11 u12; u21 u22]
	
	interp_value = (xV'*uM*tV)*coeff
	return interp_value
end

function integrate_difference(u2fun1::FunArr2D{T, U, V}, u2fun2::FunArr2D{T, U, V}, Nxdiv::I, Ntdiv::I) where {T, U, V, I<:Int}
	@assert u2fun1.xmin == u2fun2.xmin
	@assert u2fun1.xmax == u2fun2.xmax
	@assert u2fun1.tmin == u2fun2.tmin
	@assert u2fun1.tmax == u2fun2.tmax
	
	xmin = u2fun1.xmin; xmax = u2fun1.xmax
	tmin = u2fun1.tmin; tmax = u2fun1.tmax
	
	dx = (xmax - xmin)/T(Nxdiv - 1)
	dt = (tmax - tmin)/U(Ntdiv - 1)
	
	total = 0.0
	for ix in 1:(Nxdiv-1)
		x = (ix - 0.5)*dx
		for it in 1:(Ntdiv-1)
			t = (it - 0.5)*dt
			term_fun1 = bilinear_interp_regular(u2fun1, x, t)
			term_fun2 = bilinear_interp_regular(u2fun2, x, t)
			total += dx * dt * (term_fun1 - term_fun2)^2
		end
	end
	
	return total
end

function scalarwave_and_energy(ufun::FunArr{T, U}, utfun::FunArr{T, U}, Nt::I, duration::T) where {T, U, I<:Int}
	# ufun  : initial positions at t = 0
	# utfun : initial first time derivatives at t = 0
	# Nt    : total number of time steps
	# duration : total amount of time elapsed
	# returns : FunArr2D instance for the u(x, t) grid, and an array of energies
	
	Nx    = length(ufun.ydat)
	Nx_ut = length(utfun.ydat)
	if (Nx != Nx_ut)
		error("The initial positions and initial velocities Arrays must be of the same size.")
	end
	
	h = duration/(Nt-1)  # time step
	u2fun    = FunArr2D{Float64, Float64, Float64}(zeros(Nx, Nt), duration)
	energies = zeros(U, Nt)
	
	for it in 1:Nt
		# propagate forward in time
		uxfun = first_deriv(ufun)
		afun  = second_deriv(ufun)
		for ix in 1:Nx
			u_1, ut_1 = forward_euler2_disc(ix, ufun, utfun, afun, h)
			ufun.ydat[ix] = u_1
			utfun.ydat[ix] = ut_1
		end
		
		# set the values in the 2D matrix
		for ix in 1:Nx
			u2fun.u2dat[ix, it] = ufun.ydat[ix]
		end
		
		# collect the energy
		energies[it] = get_energy(utfun, uxfun)
	end
	
	return u2fun, energies
end

end # module
