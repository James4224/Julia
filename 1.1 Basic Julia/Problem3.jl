# Problem3.jl
using Statistics

function MyBinomial(n,p)
    # n is the number of draws, p is the probability
    # the output is the number of successes
    s = 0
    for i in 1:n
        s += p<rand()
    end
    return s
end

# test:
MyBinomial(100,0.5)