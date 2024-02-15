
using PhysicalConstants.CODATA2018
# Constants
h = PlanckConstant.val
muN = 5.0507837461e-27
bohr = BohrRadius.val
eps0 = VacuumElectricPermittivity.val
amu = AtomicMassConstant.val
c_l = SpeedOfLightInVacuum.val

DebyeSI = 3.33564e-30

# Rb87Cs133 Constants

Rb87Cs133 = Dict(
    "Name" => "Rb87Cs133",
    "redMass" => (132.9054 + 86.909184)*amu/2,
    "d0" => 1.225 * DebyeSI,
    "Sym" => 0,
    "Brot" => 490.173994326310e6 * h,
    "Drot" => 207.3 * h,
    "a0" => 2020*4*pi*eps0*bohr^3, #1064nm
    "a2" => 1997*4*pi*eps0*bohr^3, #1064nm
)
# K41Cs133 Constants
K41Cs133 = Dict(
    "Name" => "K41Cs133",
    "redMass" => (40.96182526 + 132.9054)*amu/2,
    "Sym" => 0,
    "d0" => 1.84 * DebyeSI,
    "Brot" => 880.326e6 * h,
    "Drot" => 0 * h,
    "a0" => 7.783e6 * h,
    "a2" => 0,
)

# K40Rb87 Constants
K40Rb87 = Dict(
    "Name" => "K40Rb87",
    "redMass" => (39.96399848 + 86.909184)*amu/2,
    "Sym" => 1,
    "d0" => 0.566 * DebyeSI,
    "Brot" => 1113.950e6 * h,
    "Drot" => 0 * h,
    "a0" => 55.3*1e-4*h*eps0*c_l*2, #0,  #1/3*(100 +2*33) *h *1e-4 , #5.53e-5 * 1e6 * h,
    "a2" => 44.7*1e-4*h*eps0*c_l*2, #2/3*(100 - 33)* h *1e-4*eps0*c_l/2, #4.47e-5 * 1e6 * h,
)


K40Rb87T = Dict(
    "Name" => "K40Rb87_Till",
    "redMass" => (39.96399848 + 86.909184)*amu/2,
    "Sym" => 1,
    "d0" => 0.573999 * DebyeSI,
    "Brot" => 1113.9514e6 * h,
    "Drot" => 0 * h,
    "a0" => 55.3*1e-4*h*eps0*c_l*2, #0,  #1/3*(100 +2*33) *h *1e-4 , #5.53e-5 * 1e6 * h,
    "a2" => 44.7*1e-4*h*eps0*c_l*2, #2/3*(100 - 33)* h *1e-4*eps0*c_l/2, #4.47e-5 * 1e6 * h,
)


Na23Cs133 = Dict(
    "Name" => "Na23Cs133",
    "redMass" => (22.98977 + 132.9054)*amu/2,
    "Sym" => 0,
    "d0" => 4.75 * DebyeSI,
    "Brot" => 1736e6 * h,
    "Drot" => 0 * h,
    "a0" => 0 * h,
    "a2" => 0 * h,
)


# Na23K40 Constants
Na23K40 = Dict(
    "Name" => "Na23K40",
    "Sym" => 1,
    "redMass" => (22.98977 + 39.96399848)*amu/2,
    "d0" => 2.72 * DebyeSI,
    "Brot" => 2.8217297e9 * h,
    "Drot" => 0 * h,
    "a0" => 0 * h,
    "a2" => 0 * h,
)

# Na23Rb87 Constants
Na23Rb87 = Dict(
    "Name" => "Na23Rb87",
    "Sym" => 0,
    "redMass" => (22.98977 + 86.909184)*amu/2,
    "d0" => 3.2 * DebyeSI,
    "Brot" => 2.0896628e9 * h,
    "Drot" => 0 * h,
    "a0" => 0 * h,
    "a2" => 0 * h,
)



