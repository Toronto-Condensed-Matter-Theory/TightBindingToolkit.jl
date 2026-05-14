using Plots, LaTeXStrings, TightBindingToolkit

"""
This script sets up the famous boxaev chain in 1d.
"""

"""
The lattice primitive vector for a honeycomb lattice : a1 and a2.
"""
a1 = [ 1.0 ]
UC_hop  = UnitCell( [a1] , 1 , 2) ##### localDim = 1 cause its a spinless Hamiltonian
UC_pair = UnitCell( [a1] , 1 , 2) ##### localDim = 1 cause its a spinless Hamiltonian

"""
The single sub-lattice of the chain
"""
b1 = [ 0.0 ]
AddBasisSite!( UC_hop  , b1 )
AddBasisSite!( UC_pair , b1 )

"""
Adding structure to the lattice now, through the bond objects.
"""
################# Nearest Neighbour hopping #################
const t1    =   1.0
t1Param     =   Param(t1, 2)   
const NNdistance  =   1.0
AddIsotropicBonds!(t1Param, UC_hop , NNdistance, 1.0 , "t1")
################# Nearest Neighbour pairing #################
const p1    =   0.5
p1Param     =   Param(p1, 2)   
AddIsotropicBonds!(p1Param, UC_pair , NNdistance, 1.0 , "p1")

CreateUnitCell!(UC_hop , [t1Param])
CreateUnitCell!(UC_pair, [p1Param])

##### BZ with size kSize
const n       =   40
const kSize   =   6 * n 
bz            =   BZ(kSize, 1)
FillBZ!(bz, UC_hop)

##### Some model parameters
const T     =   0.001
const mu    =   0.0
ModifyIsotropicFields!(UC_hop, mu, length(UC_hop.fields[begin]))
H       =   Hamiltonian(UC_hop , UC_pair , bz)
DiagonalizeHamiltonian!(H)

boxaevModel   =   BdGModel(UC_hop , UC_pair, bz, H; T=T, mu=mu)
SolveModel!(boxaevModel)

kPoints = [bz.HighSymPoints["G"], bz.HighSymPoints["M1"]]

Plot_Band_Structure!(boxaevModel, kPoints ; labels = [L"$\Gamma$", L"$M_1$"])



