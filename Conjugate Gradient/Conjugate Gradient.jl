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
    n = size(A,1)
    for k in 1:1:n
        @inbounds L[k,k] = sqrt(A[k,k])
        for i in (k+1):1:n
        @inbounds if A[i,k] != 0.0
            @inbounds L[i,k] = A[i,k]/L[k,k]
            end
        end

        for j in (k+1):1:n
            for i in j:1:n
                @inbounds if A[i,j] != 0.0
                    @inbounds L[i,j] = L[i,j]-L[i,k]*L[j,k]
                end
            end
        end
    end
end

function IncompleteCholeskyPreconditioner(A)
    L = copy(A)
    IncompleteCholeskyPreconditioner!(L,A)
    return L
end

# wikipedia test case:
 A = [[4.0 1.0];[1.0 3.0]]
 x = [2.0, 1.0]
 b = [1.0, 2.0]
tolerance = 1e-8

@btime ConjugateGradient!(x, A, b, tolerance);
L = [[1.0 0.0];[0.0 1.0]]
@btime PreconditionedConjugateGradient!(x, L, A, b, tolerance)

## Test the factorisation:
A = [[3.0, 0.0, -1.0, -1.0, 0.0, -1.0] [0.0, 2.0, 0.0, -1.0, 0.0, 0.0] [-1.0, 0.0, 3.0, 0.0, -1.0, 0.0] [-1.0, -1.0, 0.0, 2.0, 0.0, -1.0] [0.0, 0.0, -1.0, 0.0, 3.0, -1.0] [-1.0, 0.0, 0.0, -1.0, -1.0, 4.0]]
L = zeros(6,6)
@btime IncompleteCholeskyPreconditioner!(L, A)
x = ones(6)::Vector
b = ones(6)::Vector
@btime PreconditionedConjugateGradient!(x, L, A, b, tolerance)

print(A*x)