# Problem5.jl
# timeseries problem:
#import Pkg
#Pkg.add("Plots")
#Pkg.add("PyPlot")
using Plots

function ARForwards(x,α,ϵ)
    return α*x+ϵ
end

function Iterate(X,α)
    for i in 1:(200-1)
        ϵ = randn()
        X[i+1] = ARForwards(X[i],α,ϵ)
    end
    return X
end

function MyMain()
    α = [0.0, 0.5, 0.9]
    y = zeros(200,3)
    for j in 1:3
        y[:,j] = Iterate(y[:,j],α[j])
    end
    #display(plot(y, label = α'))
    return 0
end

@time MyMain()