
using LinearAlgebra
using SparseArrays

function generateSpinOperator(F::Node)
    EndNodes = endNode(F)
    Operators = []
    for (i, node) in enumerate(EndNodes)
        dimsA = i == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(i+ 1):end]])
        dimsB = i == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(i-1):-1:1]])

        I_a = sparse(1.0I, dimsA, dimsA)
        I_b = sparse(1.0I, dimsB, dimsB)

        push!(Operators, [kron(I_b, kron(x_operator(node), I_a)), kron(I_b, kron(y_operator(node), I_a)), kron(I_b, kron(z_operator(node), I_a)), kron(I_b, kron(raising_operator(node), I_a))])
    end
    Operators
end

function generateMolecule(MoleculeInfo, Nmax, Lmax)

    N1 = angularMomentum("N1", [0:Nmax...])
    N2 = angularMomentum("N2", [0:Nmax...])
    L = nothing
    if MoleculeInfo["Sym"] == 0
        L = angularMomentum("L", [0:2:Lmax...])
    else MoleculeInfo["Sym"] == 1
        L = angularMomentum("L", [1:2:Lmax...])
    end


    N_c = couple("N_c", N2, N1)
    F = couple("F", L, N_c)
    LOp, N2Op, N1Op = generateSpinOperator(F)
    moleculeProperties(MoleculeInfo,F, LOp, N2Op, N1Op)
end