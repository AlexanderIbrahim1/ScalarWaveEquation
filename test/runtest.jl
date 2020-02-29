using ScalarWaveEquation
using Test

function check_similar1D(ufun1::FunArr{T, U}, ufun2::FunArr{T, U}, tol::T) where {T, U}
	# checks if two FunArr instances are equal to within a tolerance
	N  = length(ufun1.ydat)
	N2 = length(ufun2.ydat)
	@assert N == N2
	
	dx = ufun1.dx
	total = 0.0
	for ix in 1:N
		diff = ufun1.ydat[ix] - ufun2.ydat[ix]
		total += dx*diff*diff
	end
	
	return total < tol
end

@testset "FunArr" begin
	# FunArr cannot hold an Array smaller than 3
	@test_throws ErrorException ufun = FunArr{Float64, Float64}([1 2])
	
	# Neither of FunArr2D's dimensions can be smaller than 3
	A = [1.0 2.0 3.0; 4.0 5.0 6.0]
	B = [1.0 2.0; 3.0 4.0; 5.0 6.0]
	C = ones(3, 3)
	@test_throws ErrorException u2fun = FunArr2D{Float64, Float64, Float64}(A, 1.0)
	@test_throws ErrorException u2fun = FunArr2D{Float64, Float64, Float64}(B, 1.0)
	
	# the duration of elapsed time must be positive
	@test_throws ErrorException u2fun = FunArr2D{Float64, Float64, Float64}(C, 0.0)
end

# test for 1st deriv and 2nd deriv
@testset "Derivatives" begin
	import ScalarWaveEquation
	first_deriv  = ScalarWaveEquation.first_deriv
	second_deriv  = ScalarWaveEquation.second_deriv
	
	size  = 100
	xdata = range(0.0, stop = 1.0, length = size)
	
	# create the constant array of value 1.0
	ufun1 = FunArr{Float64, Float64}(ones(Float64, size))
	
	# create the array of y = x, then take the 1st derivative
	ufunx = FunArr{Float64, Float64}(xdata)
	ufun2 = first_deriv(ufunx)
	
	# create the array of y = 0.5*x^x, then take the 2nd derivative
	x2data = 0.5 * xdata .* xdata
	ufun2x = FunArr{Float64, Float64}(x2data)
	ufun3  = second_deriv(ufun2x)
	
	tolerance = 1e-12
	@test check_similar1D(ufun1, ufun2, tolerance)
	@test check_similar1D(ufun1, ufun3, tolerance)
end

@testset "Forward Integration" begin
	import ScalarWaveEquation
	forward_euler2_disc = ScalarWaveEquation.forward_euler2_disc
	
	# test forward integration using euler2
	# y'' = 1, y'(0) = 0, y(0) = 0, should get y = 0.5*x^2
	
	size = 100
	afun = FunArr{Float64, Float64}(ones(Float64, size))
	vfun = FunArr{Float64, Float64}(zeros(Float64, size))
	yfun = FunArr{Float64, Float64}(zeros(Float64, size))
	
	dx = yfun.dx
	for it in 2:size
		y1, v1 = forward_euler2_disc(it-1, yfun, vfun, afun, dx)
		yfun.ydat[it] = y1
		vfun.ydat[it] = v1
	end
	
	# compare yfun and vfun to their exact versions
	yexact = [0.5*x^2 for x in range(0.0, stop = 1.0, length = size)]
	yfun_exact = FunArr{Float64, Float64}(yexact)
	vexact = range(0.0, stop = 1.0, length = size)
	vfun_exact = FunArr{Float64, Float64}(vexact)
	
	tolerance = 1e-12
	@test check_similar1D(yfun, yfun_exact, tolerance)
	@test check_similar1D(vfun, vfun_exact, tolerance)
end
