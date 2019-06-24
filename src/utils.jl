import Statistics: cov, mean, normalize
import LinearAlgebra: inv, det, norm
import MultivariateStats
import LightGraphs
import GraphPlot
using Plots

export plot_graphmpr

"""
    circlepoints(r, n)
Generate `n` evenly spaced points on circle of radius `r`
"""
function circlepoints(n, r; noise = 0.0)
    t = LinRange(0.0,2π,n)
    x = r .* cos.(t) .+ rand(n).*noise
    y = r .* sin.(t) .+ rand(n).*noise
    return x, y
end


"""
    balanced(x[; intervals=10, overlap=0.2])
Balanced cover with `intervals` and `overlap` fraction.
The interval boundaries are distributed so that each patch covers the same fraction of the data set.
"""
function balanced(x::AbstractVector{<:Real}; kwargs...)

    # setup parameters
    intervals=10
    overlap=0.5
    for (p,v) in kwargs
        if p == :intervals
            intervals = v
        elseif p == :overlap
            overlap = v
        end
    end

    xmin, xmax = extrema(x)
    xrng = xmax - xmin
    ilen = 1.0 / (intervals - (intervals-1)*overlap)
    istep = ilen*(1-overlap)

    # compute ranges of intervals with an overlap
    irngs = [let rmin = (i-1)*istep;
                (rmin, (rmin+ilen)) .* xrng .+ xmin
            end
            for i in 1:intervals]

    # generate indexes for cover elements
    patches = map(i->findall(e->i[1]<=e<=i[2], x), irngs)
    return patches
end

"""
    plot_graphmpr(mpr[; layout=spring_layout])
Plot the mapper output using the structure mpr given by the function mapper.
"""
function plot_graphmpr(mpr; kwargs...)
    # get the adjacency matrix form the mpr struct
    A = mpr.adj
    # turn it into a graph
    G = LightGraphs.Graph(A)
    # for each patch calculate its mean
    means = [mean(mpr.filter[p]) for p in mpr.patches]
    # normalize the means
    mns = normalize(means);
    # define the colors to plot the mapper output
    C(g::ColorGradient) = RGB[g[z] for z = mns]
    g = :algae
    d = C(cgrad(g))
    # plot the mapper graph
    GraphPlot.gplot(G, nodelabel=[size(p,1) for p in mpr.patches], nodefillc=d; kwargs...)
end
