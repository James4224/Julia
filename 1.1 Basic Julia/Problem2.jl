# solve problem 2
function MyFactorial(n)
    m = one(n)
    for i in 1:n
        m = m*i;
    end
    return m
end

MyFactorial(6)