import Base.*
using PhysicalConstants.CODATA2018

const eps0 = VacuumElectricPermittivity.val

const hbar = PlanckConstant.val/(2*pi)

function compute_SHC_elements!(operator, q, N; order = 2)
    
    for (j, state) in enumerate(getBasis(N.spin)), (k, state2) in enumerate(getBasis(N.spin))
        if true #abs(state[2] - state2[2]) <= order && -1*state[1] + state2[1] + q == 0
            mN, N1 = state
            mNp, Np = state2 

            operator[j, k] = (-1)^(mN) * sqrt((2*N1 + 1)*(2*Np + 1)) * wigner3j(N1, order, Np, -mN, q, mNp) * wigner3j(N1, order, Np, 0, 0, 0)
        end
    end
end
function DipoleMatrix(mol::moleculeProperties, basisOI)
    F = mol.basisTree
    EndNodes = endNode(F)
    dimsA = basisOI == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(basisOI+ 1):end]])
    dimsB = basisOI == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(basisOI-1):-1:1]])

    I_a = sparse(1.0I, dimsA, dimsA)
    I_b = sparse(1.0I, dimsB, dimsB)
    nodeOI = EndNodes[basisOI]
    dipoleOperator = [spzeros(ComplexF64, length(nodeOI.spin), length(nodeOI.spin)) for _ in -1:1] 
    @threads for i in 1:3
        q = i - 2
        compute_SHC_elements!(dipoleOperator[i], q, nodeOI, order = 1)
    end
    [kron(I_b, kron(sparse(dip), I_a)) for dip in dipoleOperator]

end
function getDipoleMatrix(mol::moleculeProperties, basisOI)
    dipm, dip0, dipp =  DipoleMatrix(mol, basisOI)
    dip = -1*mol.Constants["d0"]*[(dipm + adjoint(dipm))/sqrt(2), 1im*(dipm - adjoint(dipm))/sqrt(2), dip0]
end

function DCStark(M::moleculeProperties, dirE, basisOI)
    sum(dirE.* getDipoleMatrix(M, basisOI))
end

*(a::Vector{T}, b::Vector{T}) where {T<:AbstractSparseMatrix} = sum([a[i]*b[i] for i = 1:size(a)[1]])
function generateRotational(Nrot, Brot)
    Nsq = Nrot*Nrot
    Brot*Nsq
end
function generateCent(Lrot, redMass)
    Lsq = Lrot*Lrot
    hbar^2/(2*redMass)*Lsq
end

function T2_C(mol::moleculeProperties, basisOI; order = 2)
    F = mol.basisTree
    EndNodes = endNode(F)
    dimsA = basisOI == length(EndNodes) ? 1 : prod([length(node.spin) for node in EndNodes[(basisOI+ 1):end]])
    dimsB = basisOI == 1 ? 1 : prod([length(node.spin) for node in EndNodes[(basisOI-1):-1:1]])

    I_a = sparse(1.0I, dimsA, dimsA)
    I_b = sparse(1.0I, dimsB, dimsB)
    
    nodeOI = EndNodes[basisOI]
    T2C = [spzeros(ComplexF64, length(nodeOI.spin), length(nodeOI.spin)) for _ in -order:order] 
    @threads for i in 1:(2*order + 1)
        q = i - (order + 1)
        
        compute_SHC_elements!(T2C[i], q, nodeOI, order = order)
    end

    [kron(I_b, kron(sparse(q), I_a)) for q in T2C]
end



function makeT2_n(I1, I2)


    T2m2 = I1[1]*I2[1]
    T2p2 = I1[3]*I2[3]
    T2p1 = 1/sqrt(2)*(I1[3]*I2[2] + I1[2]*I2[3])
    T2m1 = 1/sqrt(2)*(I1[1]*I2[2] + I1[2]*I2[1])
    T20 = 1/sqrt(6)*(I1[3]*I2[1] + I1[1]*I2[3] + 2*I1[2]*I2[2])
   
    T = [T2m2, T2m1, T20, T2p1, T2p2]

    return T 
end


function makeT2(I1, I2)
    T2m2 = 0.5 * (I1[1] * I2[1] - 1im*I1[1] * I2[2] - 1im*I1[2] * I2[1] - I1[2] * I2[2])
    T2p2 = 0.5 * (I1[1] * I2[1] + 1im*I1[1] * I2[2] + 1im*I1[2] * I2[1] - I1[2] * I2[2])

    T2m1 = 0.5 * (I1[1] * I2[3] - 1im*I1[2] * I2[3] + I1[3] * I2[1] - 1im*I1[3] * I2[2])
    T2p1 = -0.5 * (I1[1] * I2[3] + 1im*I1[2] * I2[3] + I1[3] * I2[1] + 1im*I1[3] * I2[2])

    T20 = -sqrt(1/6) * (I1[1] * I2[1] + I1[2] * I2[2]) + sqrt(2/3) * I1[3] * I2[3]

    T = [T2m2, T2m1, T20, T2p1, T2p2]

    return T
end


function tensor_dot(T1, T2; order = 2)
    Tprod = spzeros(ComplexF64, size(T1[1])...)
    T2 = reverse(T2)
    for (i, q) in enumerate(-order:order)
        Tprod += ((-1)^q)*T1[i]*T2[i]
    end
    Tprod
end


function DDInteraction(M::moleculeProperties)
    T2 = T2_C(M, 1, order = 2)
    
    T1 = makeT2_n(T2_C(M, 2, order = 1), T2_C(M, 3, order = 1))


    -(1/(4*pi*eps0))*sqrt(6)*(M.Constants["d0"]^2)*tensor_dot(T1, T2, order = 2)
end


function generateHamiltonian(M::moleculeProperties; dirE = [0, 0, 1.0])
    consts = M.Constants
    Hrot = generateRotational(M.N1[1:3], consts["Brot"]) + generateRotational(M.N2[1:3], consts["Brot"])
    Hdc = DCStark(M, normalize!(dirE), 3) + DCStark(M, normalize!(dirE), 2)
    Hdip = DDInteraction(M)
    HCent = generateCent(M.L[1:3], consts["redMass"])
    moleculeInteraction(M,Matrix(Hrot), Matrix(Hdc), Matrix(Hdip), Matrix(HCent))
end

