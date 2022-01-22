# Operators
# define type StrangMatix


function StrangMatrixMultiply!(y,x)
# function to apply the strang matrix to a vector in an operator stype form
    y[1] = 2*x[1] + x[2]
    for i in 2:(length(x)-1)
        y[i] = x[i-1] + 2*x[i] + x[i+1]
    end
    i = length(x)
    y[i] = x[i-1] + 2*x[i]
end

function MyMain()
    x = ones(4,1)
    y = zeros(4,1)

    StrangMatrixMultiply!(y,x)
    print(y)
end

MyMain()
