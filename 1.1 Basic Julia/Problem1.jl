# script to solve problem 1:
# use the array and control flow syntax to define the NxN strang matrix:

using LinearAlgebra
N = 4
M = zeros(N,N)
for i in 1:N
    if i == 1
        M[i,i]      = -2.0
        M[i,i+1]    = 1.0
    elseif i == N
        M[i,i]      = -2.0
        M[i,i-1]    = 1.0
    else
        M[i,i-1]    = 1.0
        M[i,i]      = -2.0
        M[i,i+1]    = 1.0
    end
end

print(M)