using LaTeXStrings, Plots, LinearAlgebra, TightBindingToolkit

"""
The lattice primitive vector for a square lattice : a1 and a2.
"""
a1 = [ 1.0 , 0.0 ]
a2 = [ 0.0 , 1.0 ]
UC = UnitCell( [a1 , a2] , 2, 2) ##### localDim=2 since we are working with spin-1/2 particles now

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

"""
Now add a second unit cell for pairing
"""
UCp = UnitCell( [a1 , a2] , 2)
AddBasisSite!( UCp , b1 )
J = 2.0
AddIsotropicBonds!(UCp , NNdistance, J * SpinVec[2] , "p1" ; checkOffsetRange = 1)


##### BZ with size kSize
const n       =   40
const kSize   =   6 * n + 3
bz            =   BZ(kSize, 2)
FillBZ!(bz, UC)


"""
Some model parameters : temperature and chemical potential
"""
const T     =   0.001
const mu    =   -1.0
# ModifyIsotropicFields!(UC, mu, 4)
H       =   Hamiltonian(UC , UCp , bz)
DiagonalizeHamiltonian!(H)

"""
Filling the model at this chemical potential with spin-1/2 fermions.
"""
SquareModel   =   BdGModel(UC , UCp , bz , H; T=T, mu = mu)
SolveModel!(SquareModel)

println("The filling at μ = ",SquareModel.mu ," is ", SquareModel.filling)
Plot_UnitCell!(UC)



# plot_band_structure!(SquareModel, [bz.HighSymPoints["G"], bz.HighSymPoints["M1"]] ; labels = [L"$\Gamma$", L"$M_1$"])




