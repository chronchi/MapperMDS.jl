import Clustering
import LinearAlgebra: dot
import Clustering: cutree

"""
    noselection(algo, data; kwargs...)

Perform no clustering selection.
"""
function noselection(algo, data; kwargs...)
    k = 2
    linkage = :single
    for (p,v) in kwargs
        if p == :n_clusters
            k = v
        elseif p == :linkage
            linkage = p
        end
    end
    cl = algo(data, linkage=p)
    # returns the cluster assignment vector
    z = cutree(cl, k=k)
end

"""
    silhouette(algo, data)

Perform automatic selection of the best clustering of `data` by `algo` using silhouette method.
Clustering algorithm `algo` must be a function with one argument which determines number of clusters.
"""
function silhouette(algo, data; kwargs...)
    linkage = :single
    for (p, v) in kwargs
        if p == :linkage
            linkage = v
        end
    end

    # maximum number of clusters
    Kmax = size(data,2)
    # cluster
    cl = algo(data)
    # get the silhouette scorings for each k up to Kmax
    # since we are only using hierarchical clustering.
    S = [let C=cutree(cl, k=k)
            mean(Clustering.silhouettes(C, data)) => C
        end
    for k in 2:Kmax]
    best_k = findmax(map(first, S)) |> last
    # return a clustering corresponding to a maximum average silhouette
    # the number of clusters
    # the maximum average silhouette
    return S[best_k] |> last, best_k, S[best_k] |> first
end
