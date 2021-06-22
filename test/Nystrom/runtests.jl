using SafeTestsets

@safetestset "Nystrom mesh" begin include("nystrommesh_test.jl") end

@safetestset "Pde" begin include("pde_test.jl") end

@safetestset "Kernels" begin include("kernels_test.jl") end

@safetestset "Density" begin include("density_test.jl") end

# @safetestset "Potential tests" begin include("potential_test.jl") end

# @safetestset "Operator tests" begin include("integraloperator_test.jl") end

# @safetestset "Default assemble" begin include("assemble_test.jl") end

# @safetestset "Density interpolation" begin include("dim_test.jl") end

# @safetestset "Local density interpolation" begin include("ldim_test.jl") end

# @safetestset "Martensen-Kussmaul tests" begin include("martensen_kussmaul_test.jl") end
