MapperMDS
===

[![Build Status](https://travis-ci.org/chronchi/MapperMDS.jl.svg?branch=master)](https://travis-ci.org/chronchi/MapperMDS.jl)


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
The usage here is pretty similar to [wildart's module](https://github.com/wildart/TDA.jl).
We have to pass a distance matrix and filter values of the
distance matrix, e.g., MDS using the distance matrix and projecting
to the real line.

The package performs a hierarchical clustering, since it only needs
a distance matrix (or a dissimilarity matrix).

The linkage criteria can be given as a parameter to the function ```mapper```
via ```linkage = :criteria```, where criteria is one of the options [here](http://juliastats.github.io/Clustering.jl/stable/hclust.html#Hierarchical-Clustering-1).

You can choose the method to evaluate the quality of the clustering.
The default is no method. You can choose the silhouette
using ```clustselection = MapperMDS.silhouette```. You can then specify the number
of clusters to the hierarchical clustering with ```n_clusters = k```.
The default value is 2. Notice we can use DBScan from the package Clustering. When using DBScan it is important to use the default value of ```clustselection```.

To plot we use the package [GraphPlot](https://github.com/JuliaGraphs/GraphPlot.jl)
together with [LightGraphs](https://github.com/JuliaGraphs/LightGraphs.jl).


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

# embedding of the data represented by the distance matrix in R.
filter = classical_mds(X, 1)'[:,1]

# call mapper. note we don't have to specify the filter function,
# as filter is the filtered data already.
mpr = MapperMDS.mapper(X, filter, intervals=5, overlap=0.2)

# plot the corresponding graph
plot_graphmpr(mpr)
```
