using WaveProp
using WaveProp.Geometry
using WaveProp.Interpolation

d = ReferenceTriangle()
basis = MonomialBasis(d,1)
for el in basis
    @show el
end    
xi = [(0,0),(1,0),(0,1)]
vandermond(xi,basis)

C = lagrange_basis(xi,basis)
