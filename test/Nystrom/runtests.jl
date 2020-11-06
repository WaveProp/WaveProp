using SafeTestsets

@safetestset "Pde tests" begin include("pde_test.jl") end

@safetestset "Kernels tests" begin include("kernels_test.jl") end

@safetestset "Density tests" begin include("density_test.jl") end

@safetestset "Potential tests" begin include("potential_test.jl") end

@safetestset "Operator tests" begin include("integraloperator_test.jl") end

@safetestset "Density interpolation tests" begin include("densityinterpolation_test.jl") end
