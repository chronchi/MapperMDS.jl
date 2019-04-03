MapperMDS
===

Adaptation from [TDA.jl](https://github.com/wildart/TDA.jl) mapper
tool. Here the input is a distance matrix and it uses a hierarchical
clustering. It was inspired by the problems when one only has the
distance matrix of a dataset and not the point cloud itself.

Installation
---
In the julia REPL
```julia
using Pkg
Pkg.add("https://github.com/chronchi/MapperMDS.jl.git")
```

Usage
---
The usage here is pretty similar to wildart's module.
We have to pass a distance matrix and filter values of the
distance matrix, e.g. MDS using the distance matrix and projecting
to the real line.

The package performs a hierarchical clustering, since it only needs
a distance matrix (or a dissimilarity matrix).

The linkage criteria can be given as a parameter to the function ```mapper``` via ```linkage = :criteria```, where criteria is one of
the options [here](http://juliastats.github.io/Clustering.jl/stable/hclust.html#Hierarchical-Clustering-1).

You can choose the method to evaluate the quality of the clustering.
The silhouette method is the default. You can choose no method passing
as argument ```clustselection = noselection```. You can then specify the number
of clusters to the hierarchical clustering with ```n_clusters = k```.
The default value is 2. 

Example
---

```julia
using MapperMDS
import LinearAlgebra: rand, norm
import MultivariateStats: classical_mds

# generate a random symmetric matrix with non negative entries
X = rand(100, 100)
X -= X'
X = map(abs, X)

# embedding of the data represented by the distance matrix in ℜ.
X_embedding = classical_mds(X, 1)'[:,1]

# call mapper. note we don't have to specify the filter function,
# as X_embedding is the filtered data already.
mpr = MapperMDS.mapper(X, X_projected, intervals=5, overlap=0.2)

# plot the corresponding graph
using Plots
plot(mpr, c=:viridis)
```
