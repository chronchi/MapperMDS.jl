import Clustering
import SparseArrays
using RecipesBase

"""
    Mapper result type
"""
struct Mapper
    adj::AbstractMatrix{<:Integer}
    patches::Vector{Vector{<:Integer}}
    filter
end

Base.show(io::IO, mpr::Mapper) = show(io, "Mapper[$(length(mpr.patches))]")

"""
    mapper(X, X_embedding; kwargs...)

X is the dissimilary matrix and X_embedding is the embedding coordinates
in ℜ obtained by MDS.
"""
function mapper(X::AbstractMatrix{<:Real}, X_embedding::Array{Float64,1}; kwargs...)

    # setup parameters
    coverfn = balanced
    clusterfn = Clustering.hclust
    linkage = :single
    clusterselectionfn = silhouette
    k=2
    filter = X_embedding
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
    k_and_silh = []
    for c in covering
        #println(c)
        if length(c) == 1
            push!(patches, c)
        end
        # return the cluster assingments
        if clusterselectionfn == silhouette
            lbls, best_k, silhouette_score = clusterselectionfn(clusterfn, X[c,c]; linkage = linkage)
            push!(k_and_silh, (best_k, silhouette_score))
        else
            lbls = clusterselectionfn(clusterfn, X[c,c]; n_clusters=k)
        end
        plen = length(unique(lbls))
        #println(plen)
        for i in 1:plen
            push!(patches, c[findall(isequal(i), lbls)])
        end
    end

    # combine all patches & determine which has overlaps
    P = length(patches)
    adj = SparseArrays.spzeros(UInt8, P, P)
    for i in 1:P
        #println("$i => ", patches[i])
        for j in i+1:P
            overlap = intersect(patches[i], patches[j])
            if length(overlap) > 0
                adj[i,j] = 0x01
            end
        end
    end

    Mapper(adj, patches, filter)
end

@recipe function f(mpr::Mapper; complex_layout = circular_layout,
                   minvsize = 15, maxvsize = 35)

    xpos, ypos = complex_layout(mpr)

    # set image limits
    xlims --> extrema(xpos) .* 1.2
    ylims --> extrema(ypos) .* 1.2

    # show 1-skeleton
    for i in 1:size(mpr.adj, 1)
        idxs = findall(e->e>0, view(mpr.adj, i, :))
        push!(idxs, i)
        @series begin
            seriestype := :path
            linewidth --> 2
            linecolor --> :black
            label --> ""
            xpos[idxs], ypos[idxs]
        end
    end

    # calculate vertex attribues
    n = length(mpr.patches)
    xcrd = zeros(n)
    ycrd = zeros(n)
    zcol = zeros(n)
    msize = fill(1,n)
    for (i, p) in enumerate(mpr.patches)
        zcol[i] = mean(mpr.filter[p])
        msize[i] = length(p)
        xcrd[i] = xpos[i]
        ycrd[i] = ypos[i]
    end
    manno = map(string, msize)
    smin, smax = extrema(msize)
    srng = smax-smin
    msize = (maxvsize-minvsize).*(msize .- smin)./srng .+ minvsize

    # show nodes
    @series begin
        seriestype := :scatter
        markersize := msize
        label --> ""
        zcolor := zcol
        series_annotations := manno
        xcrd, ycrd
    end
end
