using Plots, LaTeXStrings, TightBindingToolkit

"""
This script sets up a simple triangular lattice with first, second, and third neighbour hoppings.
It then calculates the magnetic susceptibility at zero energy, along a path in the brillouin zone through the Γ, K, and M points.
"""

"""
The lattice primitive vector for a triangular lattice : a1 and a2.
"""
a1 = [ 1/2 , sqrt(3)/2 ]
a2 = [-1/2 , sqrt(3)/2 ]
UC = UnitCell( [a1 , a2] , 2) ##### localDim=2 since we are working with spin-1/2 particles now

"""
Unit cell has only one sub-lattice.
"""
b1 = [ 0.0 , 0.0 ]
AddBasisSite!( UC , b1 )

"""
Adding structure to the lattice now, through the bond objects.
"""
SpinVec     =   SpinMats(1//2) ##### Working with spin-1/2
################ Nearest Neighbour #################
const t1  =   -1.0
const NNdistance  =   1.0
AddIsotropicBonds!(UC , NNdistance, t1 * SpinVec[4] , "t1")
############### 2nd Nearest Neighbour #################
const t2  =   0.0
const secondNNdistance  =   sqrt(3)
AddIsotropicBonds!(UC , secondNNdistance, t2 * SpinVec[4] , "t2" ; checkOffsetRange = 2)
############### 3rd Nearest Neighbour #################
const t3  =   -0.0
const thirdNNdistance  =   2.0
AddIsotropicBonds!(UC , thirdNNdistance, t3 * SpinVec[4] , "t3", checkOffsetRange = 3)

""" 
Interaction on the lattice
"""
JMatrix =   zeros(Float64, 4, 4)
JMatrix[1, 1]   =   1.0
JMatrix[2, 2]   =   1.0
JMatrix[3, 3]   =   1.0

JParam  =   Param(1.0, 2)
AddIsotropicBonds!(JParam, UC, NNdistance, JMatrix, "NN Heisenberg")



##### BZ with size kSize
const n       =   40
const kSize   =   6 * n + 3
bz            =   BZ(kSize)
FillBZ!(bz, UC)

"""
Some model parameters : temperature and chemical potential
"""
const T     =   0.001
const mu    =   1.0   
H       =   Hamiltonian(UC, bz)
DiagonalizeHamiltonian!(H)
p = Plot_UnitCell(UC ; range=2)
display(p)
"""
Filling the model at this chemical potential with spin-1/2 fermions.
"""
TriangleModel   =   Model(UC, bz, H; T=T, mu=mu)
SolveModel!(TriangleModel)
println("The filling at μ = 0 is ", TriangleModel.filling)

kPoints     =   [bz.HighSymPoints["G"], bz.HighSymPoints["K1"], bz.HighSymPoints["M2"]]
"""
A path in the discretized BZ which goes through Γ -> K -> M -> Γ.
"""
path = CombinedBZPath(bz, kPoints ; nearest=true)

# band_structure  =   Plot_Band_Structure!(TriangleModel, kPoints ; labels = [L"$\Gamma$", L"$K_1$", L"$M_2$"])
# FS              =   Plot_FS!(TriangleModel.Ham , TriangleModel.bz , [TriangleModel.mu] , [1])

Omegas  =   [0.0]
const spread  = 5e-3  
# dos = DOS(Omegas, H)    ##### the single particle density of states

suscep  =   Susceptibility(path, Omegas, TriangleModel ; eta=spread)
# FillChis!(suscep, TriangleModel)
FillRPAChis!(suscep, TriangleModel, [JParam] ; Generators = [[1, 1], [2, 2], [3, 3]])
##### Plotting the real part of the magnetic susceptibility at 0 energy, and momentum points along the path.
# plot(real.(suscep.chis["zz"][1, :]), labels="χ_zz(Ω=0, Q)", lw=2.0, markershape=:circle, markercolor="red")



