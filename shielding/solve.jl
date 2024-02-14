#Diagonalizing Function 

mutable struct sol
    R
    Ef
    vec
    val
end

function Hpot(H::moleculeInteraction, R::T, Edc::T)  where (T<:Number)
    Hermitian(Matrix(R == 0.0 ? 0 : 1/R^3*H.Hdip + Edc*H.Hdc + H.Hrot + H.Hcent/R^2))
end
function diagonalize(H::moleculeInteraction, R::T, Edc::T) where {T<:Number}
    J2uK = 1/(1.38e-23)*1e6
    vals, vecs = eigen(Hpot(H, R, Edc))
    sorted_indices = sortperm(vals)
    vals = vals[sorted_indices]
    vecs = vecs[:, sorted_indices]
    sol(R, Edc, vecs, vals*J2uK)
end
function scan(H, R::Array, Edc)
    solutions = Vector{sol}(undef, length(R))
    for (i, R_i) in enumerate(R)
        solutions[i] = diagonalize(H, R_i , Edc)
    end
    solutions
end
function scan(H, R, Edc::Array)
    solutions = Vector{sol}(undef, length(Edc))
    
    for (i, E_i) in enumerate(Edc)
        
        solutions[i] = diagonalize(H, R, E_i)
    end
    solutions
end