# Problem1.jl
# my range and my linspace functionality
struct MyRange
    start::Number
    step::Number
    stop::Number
end

A = MyRange(0.0, 1.0, 10.0)

function _MyRange(a::MyRange, i::Integer)
    nSteps = (a.stop - a.start)/a.step
    if i > nSteps
       print("Out of bounds ")
        i = nSteps
    elseif i < 0
        print("Out of bounds ")
        i = 0
    end
    a = a.start + i * a.step
    return a
end

for i in 0:11
    print( _MyRange(A, i) )
    print(" ")
end

print(range(0.0,10.0,2))

print(0.0:1.0:10.0)

struct MyLinspace
    start::Number
    step::Number
    stop::Number
end

B = MyLinspace()