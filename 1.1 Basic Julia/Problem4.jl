# Problem4.jl
# Monte carlo for pi:

function MonteCarloPi(n)
    c = 0
    for i in 1:n
        x = rand()
        y = rand()
        if x*x + y*y <= 1
            c += 1
        end
    end
    return 4*c/n
end

MonteCarloPi(1000000)