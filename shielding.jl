module shielding_jl
using Reexport
include("QuantumToolbox/QuantumWrapper.jl")
using .QuantumWrapper
import .QuantumWrapper: AM_Toolbox
import .QuantumWrapper.State_Toolbox.State
import .QuantumWrapper.State_Toolbox.KetName
import .QuantumWrapper.State_Toolbox.Ket
import .QuantumWrapper.AM_Toolbox: getBasisUC, getBasis


export AM_Toolbox, getBasisUC, getBasis, State, KetName, Ket


@reexport module MoleculeTypes
include("shielding/MoleculeConstants.jl")
export Rb87Cs133, K41Cs133, K40Rb87 ,Na23K40, Na23Rb87, K40Rb87T, Na23Cs133
end

@reexport module Molecule

import ..QuantumWrapper.AM_Toolbox: Basis, Node, angularMomentum, couple, endNode, AngularNode, getBasis, getBasisUC
using ..QuantumWrapper.SpinOperator
include("shielding/Molecule.jl")
struct moleculeProperties{T<:Node}
    Constants
    basisTree::T
    L
    N2
    N1
end
export moleculeProperties, generateMolecule

end


@reexport module Hamiltonian
import ..Molecule: moleculeProperties
import ..AM_Toolbox: getBasisUC, endNode, getBasis
using LinearAlgebra
using SparseArrays
using Base.Threads
using WignerSymbols

struct moleculeInteraction{T}
    M::moleculeProperties
    Hrot::T
    Hdc::T
    Hdip::T
    Hcent::T
end

include("shielding/Hamiltonian.jl")

export generateHamiltonian, moleculeInteraction
getBasisUC(mol::moleculeInteraction) = getBasisUC(mol.M.basisTree)

end

@reexport module solve
import ..Hamiltonian: moleculeInteraction
using LinearAlgebra
include("shielding/solve.jl")
export sol, Hpot, diagonalize, scan
end

@reexport module calculate
import ..QuantumWrapper.State_Toolbox.State
import ..QuantumWrapper.State_Toolbox.KetName
import ..QuantumWrapper.State_Toolbox.Ket
import ..AM_Toolbox: getBasisUC, endNode, getBasis

import ..Hamiltonian: moleculeInteraction
import ..solve: sol

include("shielding/calculate.jl")
export findMaxComponent, findState,findStateMax,  adiabats, energies

end

@reexport module plotting
using Reexport
import ..Hamiltonian: moleculeInteraction
import ..solve: sol
import ..AM_Toolbox: getBasisUC, endNode, getBasis
import ..QuantumWrapper.State_Toolbox.State
import ..QuantumWrapper.State_Toolbox.KetName
import ..QuantumWrapper.State_Toolbox.Ket

using PhysicalConstants.CODATA2018
include("shielding/plotting.jl")

@reexport using PlotlyJS
export plotScan



    

end




end