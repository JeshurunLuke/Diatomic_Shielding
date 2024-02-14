
State(comp::Array{ComplexF64, 2}, QN::Array{Float64, 2}, mol::moleculeInteraction) = State(normalize!(comp), QN, getBasisUC(mol.M.basisTree))
State(QN::Vector{<:Float64}, mol::moleculeInteraction) = State([1.0], [QN], getBasisUC(mol.M.basisTree) )
State(comp::Vector{<:ComplexF64}, mol::moleculeInteraction) = State(comp, getBasisUC(mol.M.basisTree) , getBasisUC(mol.M.basisTree) )
KetName(state::Vector{<:ComplexF64}, mol::moleculeInteraction; QMorder = [6, 5, 4, 3, 2, 1]) = KetName(state, getBasisUC(mol.M.basisTree), QMorder = QMorder )
KetName(state::State, mol::moleculeInteraction; QMorder = [6, 5, 4, 3, 2, 1]) = KetName(Ket(state), getBasisUC(mol.M.basisTree), QMorder = QMorder )
KetName(state::State) = KetName(Ket(state), state.basis)


function findMaxOverlap(basisState::Vector{<:ComplexF64}, eigvec)
    argmax(abs.(transpose(eigvec)*basisState))
end
function findMaxComponent(basisState::Vector{<:ComplexF64}, eigvec)
    B = zeros(ComplexF64, size(eigvec))
    for j in 1:size(eigvec, 2)
        max_index = argmax(abs.(eigvec[:, j]))
        B[max_index, j] = eigvec[max_index, j]
    end
    argmax(abs.(transpose(B)*basisState))

end

findState(stateOI::State, eigsol::Array{<:ComplexF64, 2}) = findState(Ket(stateOI), eigsol)
findState(stateOI::State, eigsol::sol) = findState(Ket(stateOI), eigsol.vec)
findState(stateOI::Vector{<:Complex}, eigsol::sol) = findState(stateOI, eigsol.vec)
function findState(stateOI::Vector{<:Complex}, eigsol::Array{ComplexF64, 2})
    findMaxOverlap(stateOI, eigsol)
end
findStateMax(stateOI::State, eigsol::Array{<:ComplexF64, 2}) = findStateMax(Ket(stateOI), eigsol)
findStateMax(stateOI::State, eigsol::sol) = findStateMax(Ket(stateOI), eigsol.vec)
findStateMax(stateOI::Vector{<:Complex}, eigsol::sol) = findStateMax(stateOI, eigsol.vec)
function findStateMax(stateOI::Vector{<:Complex}, eigsol::Array{ComplexF64, 2})
    findMaxComponent(stateOI, eigsol)
end

adiabats( sols::Vector{sol}, state::State) = adiabats( sols::Vector{sol}, Ket(state))
function adiabats(sols::Vector{sol},state::Vector{<:ComplexF64},  offset = true)
    startingInd = argmax([sol_i.R for sol_i in sols])
    adiabatInd = findStateMax(state, sols[startingInd])
    println("state: $adiabatInd")
    Upot = offset ? [real(sol_i.val[adiabatInd]) for sol_i in sols] .- sols[startingInd].val[adiabatInd] : [real(sol_i.val[adiabatInd]) for sol_i in sols] 
    stateVecs = [sol_i.vec[adiabatInd] for sol_i in sols]
    stateVecs, Upot
end
function energies(sols::Vector{sol}; units = "K")
    yConv = Dict("K"=> 1/1e6, "uK" => 1)
    hcat([i.val*yConv[units] for i in sols]...)
end

