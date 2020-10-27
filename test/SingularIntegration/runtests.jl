using SafeTestsets

@safetestset "Singularity handler" begin include("singularityhandler_test.jl") end

#     @testset "Quadrature 1d" begin
#         using SingularIntegration: Fejer1d, Generic1d, Gauss1d, integrate

#         N     = 40
#         kress = Kress()
#         q     = Fejer1d(N)
#         xs    = 0.1 # singularity location
#         qnew  = kress(q;xs=xs)

#         # for smooth integrands, both quadratures should work
#         f = (x) -> cos(x)
#         val, _ = quadgk(f,-1,1)
#         @test integrate(f,q) ≈ val
#         @test integrate(f,qnew) ≈ val

#         # for singular integrand the regular quadrature should fail, while the transformed quad shoudl work
#         f = (x) -> cos(x)*log(abs(x-xs))
#         val, _ = quadgk(f,-1,xs,1)
#         @test integrate(f,q) ≉ val
#         @test integrate(f,qnew) ≈ val

#     end

#     @testset "Quadrature 2d" begin
#         using SingularIntegration: Fejer2d, Generic2d

#         M,N   = 100,100
#         kress = Kress()
#         q     = Fejer2d(M,N)
#         xs     = (0.1,0.) # singularity location
#         qnew  = kress(q,xs=xs)

#         # for smooth integrands, both quadratures should work
#         using HCubature: hcubature
#         f = (x) -> cos(x[1]*x[2])
#         val, _ = hcubature(f,(-1,-1),(1,1))
#         @test integrate(f,q) ≈ val
#         @test integrate(f,qnew) ≈ val

#         # for singular integrand the regular quadrature should fail, while the transformed quad should work
#         f = (x) -> cos(x[1]*x[2])*log(sum((x.-xs).^2))
#         val, _ = hcubature(f,(-1,-1),(1,1))
#         @test integrate(f,q) ≉ val
#         @test integrate(f,qnew) ≈ val

#         f = (x) -> cos(x[1]*x[2])/sqrt(sum((x.-xs).^2))
#         val, _ = hcubature(f,(-1,-1),(1,1))
#         @test integrate(f,q) ≉ val
#         @test integrate(f,qnew) ≈ val
#         @btime integrate($f,$qnew)

#         # using SingularIntegration:  RectangularPolar2d
#         # sing_handler = RectangularPolar2d{2}((0.,0.))

#         # f = (x) -> cos(x[1]^2 + x[2]^2)
#         # kress = Kress{10}()
#         # hcubature(f,(-1,-1),(1,1))
#         # hcubature(kress,f,(-1,-1),(1,1),xs=(0,0))

#         # f = (x) -> log(x[1]^2 + x[2]^2)
#         # kress = Kress{2}()
#         # hcubature(f,(-1,-1),(0,0))
#         # hcubature(kress,f,(-1,-1),(0,0),xs=(0,0))

#         # @btime hcubature($f,(-1,-1),(0,0))
#         # @btime hcubature($kress,$f,(-1,-1),(0,0),xs=(0,0))

#         # f = (x) -> log(x[1]^2 + x[2]^2)
#         # imt = IMT{1,1}()
#         # hcubature(f,(-1,-1),(0,0))
#         # hcubature(imt,f,(-1,-1),(0,0),xs=(0,0))

#         # @btime hcubature($f,(-0.1,-0.1),(0,0))
#         # @btime hcubature($kress,$f,(-0.1,-0.1),(0,0),xs=(0,0))

#         # f = (x) -> (x[1] == 0 && x[2] == 0 ) ? 0.0 : log(sqrt(x[1]^2 + x[2]^2))
#         # val,_  = hcubature(f,(-1,-1),(1,1))
#         # val2,_ = hcubature(sing_handler,f,(-1,-1),(1,1))
#         # @test val ≈ val2

#         # @btime hcubature($f,(-1,-1),(1,1))
#         # @btime hcubature($sing_handler,$f,(-1,-1),(1,1))

#         # f = (x) -> (x[1] == 0 && x[2] == 0 ) ? 0.0 : 1/(sqrt(x[1]^2 + x[2]^2))
#         # val,_  = hcubature(f,(-1,-1),(0,0))
#         # val2,_ = hcubature(sing_handler,f,(-1,-1),(0,0))
#         # @test val ≈ val2

#         # @btime hcubature($f,(-1,-1),(0,0))
#         # @btime hcubature($sing_handler,$f,(-1,-1),(0,0))
#         # hcubature(f,(-1,-1),(1,1))

#     end
#     # end

#     # using SingularIntegration: Fejer1d, fejer1, ReferenceSegment, integrate
#     # @testset "Reference element" begin
#     #     N  = 10
#     #     f  = (x) -> cos(x)
#     #     qrule = Fejer1d(N)
#     #     val = integrate(f,qrule)
#     #     @test val ≈ sin(1) - sin(-1)
#     #     val = integrate(f,()->fejer1(10))
#     #     @test val ≈ sin(1) - sin(-1)

#     #     using SingularIntegration: ReferenceSquare, Fejer2d
#     #     M,N  = 10,12
#     #     f  = (x) -> cos(x[1])*cos(x[2])
#     #     qrule = Fejer2d(M,N)
#     #     val = integrate(f,qrule)
#     #     @test val ≈ (sin(1) - sin(-1))*(sin(1) - sin(-1))
#     # end

#     # @testset "Element" begin
#     #     using SingularIntegration: GenericElement
#     #     N  = 10
#     #     f  = (x) -> 1
#     #     ref = ReferenceSegment()
#     #     el  = GenericElement(t->(cos(π*t),sin(π*t)),ref)
#     #     SingularIntegration.jacobian(::typeof(el),x) = π
#     #     qrule = Fejer1d(N)
#     #     val = integrate(f,el,qrule)
#     #     @test val ≈ 2π
#     #     using SingularIntegration: precompute_nodes_and_weights
#     #     x,w = precompute_nodes_and_weights(el,qrule)
#     #     @test mapreduce(+,zip(x,w)) do (x,w)
#     #         f(x)*w
#     #     end ≈ 2π
#     #     # using SingularIntegration: ReferenceSquare, Fejer2d
#     #     # M,N  = 10,12
#     #     # f  = (x) -> cos(x[1])*cos(x[2])
#     #     # element = ReferenceSquare()
#     #     # qrule = Fejer2d(M,N)
#     #     # val = integrate(f,element,qrule)
#     #     # @test val ≈ (sin(1) - sin(-1))*(sin(1) - sin(-1))
#     # end


# #     x,w = precompute_nodes_and_weights(element,quadrature)

# #     x,w = precompute_nodes_and_weights(kernel,element,quadrature)

# #     w = precompute_weights(kernel,element,quadrature,nodes)

# #     N  = 100
# #     x,w = fejer1(N)
# #     @test sum(w) ≈ 2
# #     @test mapreduce(+,zip(x,w)) do (x,w)
# #         cos(x)*w
# #     end ≈ sin(1) - sin(-1)

# #     Q = Fejer1d(N)
# #     @test integrate(cos,Q) ≈ sin(1) - sin(-1)
# #     alloc = @ballocated integrate($cos,$Q)
# #     @test alloc == 0
# #     @test Q(cos) ≈ sin(1) - sin(-1)

#     # xs = 0.1
#     # cov = SingularIntegration.RectangularPolar{7}(xs)
#     # @test maximum(abs(ForwardDiff.derivative(cov,x) - SingularIntegration.derivative(cov,x)) for x in -1:0.1:1) < 1e-10

# #     f  = (x) -> log(abs(x-xs))
# #     fc = x -> f(cov(x))*ForwardDiff.derivative(cov,x)
# #     @test quadgk(f,-1,xs,1)[1] ≈ quadgk(fc,-1,xs,1)[1]
# #     val = quadgk(f,-1,xs,1,rtol=1e-14)[1]
# #     @test Q(f) ≉ val
# #     @test integrate(f,Q,cov) ≈ val

# #     # @btime quadgk(fc,-1,xs,1)
# #     # @btime integrate($f,$Q,$cov)
# #     # @btime integrate($f,$Q)

# #     # integrate(f,Q,singularity_handler=RecangularPolar{2})

# end



# # # using SingularIntegration: derivative
# # # using ForwardDiff

# # # x   = -1:0.1:1 |> collect
# # # plot(x,cov.(x))

# # # y = [ForwardDiff.derivative(cov,x) for x in x]
# # # plot(x,y)

# # # using Plots
# # # fig = plot()
# # # for p=2:5
# # #     cov = RectangularPolar{p}(0.6)
# # #     x   = -1:0.1:1 |> collect
# # #     xt  = cov.(x)
# # #     using Plots
# # #     plot!(fig,x,xt,m=:x,label="p=$p")
# # # end
# # # display(fig)
