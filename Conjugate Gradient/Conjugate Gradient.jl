using LinearAlgebra
using BenchmarkTools

function ConjugateGradient!(x,A,b,tolerance)
    # conjugate gradient method implemented in Julia
    # Solve b/A = x where x'Ax > 0 for all real x
    r = b - A*x
    err = dot(r,r)
    p = r
    k = 0
    while err > tolerance
        Ap = A*p
        alpha = err/dot(p,Ap)
        # need to fill x for proper mutation
        x[:] = x + alpha * p
        rk = r - alpha * Ap
        err = rk'*rk
        beta = err/dot(r,r)
        p = rk + beta*p

        # update the loop
        r = rk
        k = k+1
    end
end

function ConjugateGradient(x,A,b,tolerance)
    y = copy(x)
    ConjugateGradient!(y,A,b,tolerance)
    return y 
end

function PreconditionedConjugateGradient!(x, Minv, A, b, tolerance)
    # conjugate gradient method implemented in Julia
    # Solve b/A = x where x'Ax > 0 for all real x
    r = b - A*x
    err = dot(r,r)
    z = Minv*r
    p = z
    k = 0
    while err > tolerance
        Ap = A*p
        alpha = err/dot(p,Ap)
        # need to fill x for proper mutation
        x[:] = x + alpha * p
        rk = r - alpha * Ap
        err = dot(rk,rk)
        zk = Minv*rk
        beta = dot(rk,zk)/dot(r,z)
        p = zk + beta*p

        # update the loop
        r = rk
        k = k+1
        z = zk
    end
end

function PreconditionedConjugateGradient(x, Minv, A, b, tolerance)
    y = copy(x)
    PreconditionedConjugateGradient!(y, Minv, A, b, tolerance)
    return y
end

function IncompleteCholeskyPreconditioner!(L::Matrix,A::Matrix)
    for i in 1:1:size(A,1)
        s = 0.0
        for k in 1:1:(i-1)
            s = s + L[i,k]*L[i,k]
        end
        L[i,i] = sqrt(A[i,i] - s)

        for j in (i+1):1:(size(A,1))
            s = 0.0
            for k in 1:1:(i-1)
                s = s + L[i,k]*L[j,k]
            end
            if (L[i,i]>0.0)
                L[j,i] = (A[j,i] - s)/L[i,i] 
            end
        end
    end
end

function IncompleteCholeskyPreconditioner(A)
    L = 0.0.*A
    IncompleteCholeskyPreconditioner!(L,A)
    return L
end

# wikipedia test case:
# A = [[4.0 1.0];[1.0 3.0]]
# x = [2.0, 1.0]
# b = [1.0, 2.0]
tolerance = 0.000000001

#@btime ConjugateGradient!(x, A, b, tolerance);
#L = [[1.0 0.0];[0.0 1.0]]
#@btime PreconditionedConjugateGradient!(x, Minv, A, b, tolerance)
#y = PreconditionedConjugateGradient(x, L, A, b, tolerance)
#print(y)

## Test the factorisation:
A = [[3.0, 0.0, -1.0, -1.0, 0.0, -1.0] [0.0, 2.0, 0.0, -1.0, 0.0, 0.0] [-1.0, 0.0, 3.0, 0.0, -1.0, 0.0] [-1.0, -1.0, 0.0, 2.0, 0.0, -1.0] [0.0, 0.0, -1.0, 0.0, 3.0, -1.0] [-1.0, 0.0, 0.0, -1.0, -1.0, 4.0]]
L = zeros(6,6)
IncompleteCholeskyPreconditioner!(L, A)
Minv = I/L
x = ones(6)::Vector
b = ones(6)::Vector

#y = PreconditionedConjugateGradient(x, Minv, A, b, tolerance)
ConjugateGradient!(x,A,b,tolerance)
print(x)