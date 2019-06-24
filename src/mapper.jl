import Clustering
import SparseArrays
using RecipesBase

export mapper

"""
    Mapper result type
"""
struct Mapper
    adj::AbstractMatrix{<:Integer}
    patches::Vector{Vector{<:Integer}}
    filter::Vector{<:Float64}
    covering::Vector{Vector{<:Integer}}
end

Base.show(io::IO, mpr::Mapper) = show(io, "Mapper[$(length(mpr.patches))]")

"""
    mapper(X, filter; kwargs...)

X is the dissimilary matrix and filter is a vector in Rn, where
n is the number of points in the dataset and each entry corresponds
to a point, it is the value of the filter function.
"""
function mapper(X::Matrix{<:Float64}, filter::Array{<:Float64}; kwargs...)

    # setup parameters
    coverfn = balanced
    clusterfn = Clustering.hclust
    linkage = :single
    clusterselectionfn = noselection
    k=2
    for (p,v) in kwargs
        if p == :cover
            @assert isa(v, Function) "`$p` parameter must be function"
            coverfn = v
        elseif p == :clustering
            @assert isa(v, Function) "`$p` parameter must be function"
            clusterfn = v
        elseif p == :clustselection
            @assert isa(v, Function) "`$p` parameter must be function"
            clusterselectionfn = v
        elseif p == :seed
            Random.seed!(v)
        elseif p == :n_clusters
            k = v
        elseif p == :linkage
            linkage = p
        end
    end

    # construct cover of the filter range
    covering = coverfn(filter; kwargs...)

    # using clustering algorithm create patches from cover elements
    patches = Vector{Int}[]

    for c in covering
        if length(c) == 1
            push!(patches, c)
        else
            if clusterselectionfn == silhouette
                k_and_silh = []
                lbls, best_k, silhouette_score = clusterselectionfn(clusterfn, X[c,c]; linkage = linkage)
                push!(k_and_silh, (best_k, silhouette_score))
            else
                lbls = clusterselectionfn(clusterfn, X[c,c]; kwargs...)
            end
            for i in unique(lbls)
                push!(patches, c[findall(isequal(i), lbls)])
            end
        end
    end

    # combine all patches & determine which has overlaps
    P = length(patches)
    adj = SparseArrays.spzeros(UInt8, P, P)
    for i in 1:P
        for j in i+1:P
            overlap = intersect(patches[i], patches[j])
            if length(overlap) > 0
                adj[i,j] = 0x01
                adj[j,i] = 0x01
            end
        end
    end

    Mapper(adj, patches, filter, covering)
end
