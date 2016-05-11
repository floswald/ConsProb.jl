module test_accelerator


using ConsProb, FactCheck


facts("testing searchsorted1step") do

	n = 11
	y = collect(linspace(1,11.0,n))
	x = 8.9
	i = 1
	@fact ConsProb.searchsorted1step(y,x,i,n) --> 8

	n = 13
	x = 0.9
	i = 4
	@fact ConsProb.searchsorted1step(y,x,i,n) --> 3

	n = 15
	y = collect(linspace(1,15.0,n))
	x = 1.999
	i = 1
	@fact ConsProb.searchsorted1step(y,x,i,n) --> 1

	x = 2.999
	@fact ConsProb.searchsorted1step(y,x,i,n) --> 2
	x = 8.999
	@fact ConsProb.searchsorted1step(y,x,i,n) --> 8

end


facts("testing linearapprox with accelerator and 1 step bracket search") do

	acc = ConsProb.Accelerator(1)
	n = 2000
	x = collect(linspace(0.0,n,n))
	y = 2.29 .* x

	@fact ConsProb.linearapprox(x,y,0.0,n,acc) --> 0.0
	@fact acc.i --> 1
	@fact ConsProb.linearapprox(x,y,1.1,n,acc) --> roughly(1.1*2.29)
	@fact acc.i --> 2
	@fact ConsProb.linearapprox(x,y,1500.0,n,acc) --> roughly(1500*2.29)
	@fact acc.i --> 1500

	ConsProb.resetAccel!(acc)

	@fact ConsProb.linearapprox2(x,y,0.0,n,acc) --> 0.0
	@fact acc.i --> 1
	@fact ConsProb.linearapprox2(x,y,1.1,n,acc) --> roughly(1.1*2.29)
	@fact acc.i --> 2
	@fact ConsProb.linearapprox(x,y,154.0,n,acc) --> roughly(154*2.29)
	@fact acc.i --> 154

end

facts("benchmarking against searchsortedlast") do


	function f(n)
		acc = ConsProb.Accelerator(1)
		x = collect(linspace(0.0,n,n))
		y = 2.29 .* x
		y0 = 0.0
		for i=1:n
			y0 = ConsProb.linearapprox(x,y,x[i],n,acc)
		end
	end

	function f2(n)
		acc = ConsProb.Accelerator(1)
		x = collect(linspace(0.0,n,n))
		y = 2.29 .* x
		y0 = 0.0
		for i=1:n
			y0 = ConsProb.linearapprox2(x,y,x[i],n,acc)
		end
	end

	# warm up
	f(10)
	f2(10)

	# run
	# n = 200
	# t1 = @elapsed f(n)
	# t2 = @elapsed f2(n)
	# println("benchmarking with n=$n:")
	# println("	linearapprox:  $t1")
	# println("	linearapprox2: $t2")

	# n = 2000
	# t1 = @elapsed f(n)
	# t2 = @elapsed f2(n)
	# println("benchmarking with n=$n:")
	# println("	linearapprox:  $t1")
	# println("	linearapprox2: $t2")

	n = 20000
	t1 = @elapsed f(n)
	t2 = @elapsed f2(n)
	println("benchmarking with n=$n:")
	println("	linearapprox:  $t1")
	println("	linearapprox2: $t2")

	n = 2000000
	t1 = @elapsed f(n)
	t2 = @elapsed f2(n)
	println("benchmarking with n=$n:")
	println("	linearapprox:  $t1")
	println("	linearapprox2: $t2")

end


end