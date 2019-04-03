srcroot = "$(dirname(@__FILE__))"

using MapperMDS
using Test
import LinearAlgebra: rand, norm
import MultivariateStats: classical_mds

Points_X = rand(100,3);
X = zeros(100, 100)

for i in 1:size(Points_X)[1]
    for j in i+1:size(Points_X)[1]
        dij = norm(Points_X[i,:] - Points_X[j,:])
        X[i,j] = dij
        X[j,i] = dij
    end
end

X_projected = classical_mds(X, 1)'[:,1]

mpr = MapperMDS.mapper(X, X_projected)

@test typeof(mpr.patches) == Vector{Vector{<:Integer}}
