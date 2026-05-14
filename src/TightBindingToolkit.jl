module TightBindingToolkit

####### Write your package code here.

##### Module to get SU(2) representation for any half-integer / integer spin. Also have a function to return traceless matrix basis for Hermitian NxN matrices.
include("SpinMats.jl")
using .SpinMatrices
export SpinMats, HermitianBasis

include("Useful.jl")
using .Useful
export GetAllOffsets, VecAngle, Meshgrid, BinarySearch, DistFunction, DeriDistFunction , GetIndexPath, FFTArrayofMatrix, GetPhase #, Central_Diff, Arrayfy, DeArrayfy

##### Module for defining Unit Cell structure and some funcitons to manipulate it.
include("UnitCell.jl")
using .UCell
export Bond , BondRank, IsSameBond , FlipBond , UnitCell , IsSameUnitCell , AddBasisSite! , GetDistance , GetRealSpacePositions, IsSameUnitCell

include("DesignUnitCell.jl")
using .DesignUCell
export AddAnisotropicBond! , AddIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!, Lookup

include("ExpandUnitCell.jl")
using .ExpandUCell
export ChangePrimitives!, ExpandUnitCell, ExpandBonds!

##### Module to define Param structure which represents a general hopping parameter, and can be used to construct a UnitCell.
include("Params.jl")
using .Parameters
export Param, AddAnisotropicBond!, AddIsotropicBonds!, AddSimilarBonds! , CreateUnitCell!, ModifyUnitCell!, GetParams, Lookup

##### Module to make a lattice in real-space using a UnitCell.
include("Lattice.jl")
using .LatticeStruct
export Lattice, FillSites!, FillBonds!, FillLattice!, GetBCPhase, ApplyBCToSite

include("DesignLattice.jl")
using .DesignLattice
export CreateLattice, ModifyLattice!, ScaleLatticeBonds!, ScaleLattice!, RemoveLatticeBonds!

include("Plaquettes.jl")
using .Plaqs
export Plaquette, PlaquetteParam, AddIsotropicPlaquettes!, GetPlaquetteFlux, GetPlaquetteSites

include("Flux.jl")
using .Flux
export GetBondPhase, CheckGaugeValidity, ModifyGauge!, GetStringPhase, InsertMonopolePair!

include("Gauge.jl")
using .Gauge
export LandauGauge, RadialGauge

##### Module to define Brillouin Zone structure and some functions to manipulate it.
include("BZ.jl")
using .BZone
export GetRLVs , BZ , Monkhorst , FillBZ! , ReduceQ , GetQIndex , CombinedIndexPath , GetBZPath , CombinedBZPath , MomentumPhaseFFT

##### module to define a Hamiltonian structure and solving for some of its properties.
include("Hamiltonian.jl")
using .Hams
export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , DiagonalizeHamiltonian! , DOS, ModifyHamiltonianField!, IsBandGapped, GetVelocity!

include("LatticeHamiltonian.jl")
using .LatHam
export FillHamiltonian, LatticeHamiltonian, DiagonalizeHamiltonian!, Slater, SingleParticleFidelity, SlaterOverlap, GaugeFix!

##### Module to define a Tight-Binding Model structure which takes into account thermodynamical parameters such as temperature and filling etc.
include("Model.jl")
using .TBModel
export Model , FindFilling , GetMu! , GetFilling! , GetCount , GetGk! , GetGr!, SolveModel!, GetGap!, FreeEnergy, GetOrderParameter, FindEntropy, GetBondCoorelation, GetRSEnergy

include("LatticeModel.jl")
using .LatModel
export LatticeModel , GetMu! , GetFilling! , GetGr!, SolveModel!, GetGap!

################# Symmetries acting on lattice #################
include("LatticeSymmetries.jl")
using .LatticeSymmetries
export Translations, Degeneracies, FindQuantumNumbers

##### Module to define the equivalent but for bdG systems with pairing.
include("BdGModel.jl")
using .BdG
export BdGModel, FindFilling , GetMu! , GetFilling! , GetGk! , GetGr!, SolveModel!, GetGap!, FreeEnergy, FindEntropy, GetGroundEnergy, GetRSEnergy

##### Module to calculate generalized Chern numbers
include("Chern.jl")
using .Chern
export FindLinks , FieldStrength , ChernNumber, CheckValidity, PartialChernNumber, FilledChernNumber

##### Module to calculate generalized magnetic susceptibility
include("Susceptibility.jl")
using .suscep
export Susceptibility , FindChi , FillChis!, FillRPAChis!, PerformRPA!, FindReducedChi, nF_factor, ResolvedQInteraction, ReduceChis!

##### Module to calculate generalized electrical conductivity
include("Conductivity.jl")
using .conduct
export Conductivity, SpectralFunction, SpectralContribution, GetConductivity!

##### Module containing some plotting functions
include("Plot.jl")
using .PlotTB
export Plot_UnitCell , Plot_Band_Contour, Plot_Band_Structure, Plot_FS, Plot_Fields, Plot_Lattice

##### Module for quantum geometric tensor and related metric quantities
include("Metric.jl")
using .MetricTB
export GeoTensor, Curvature, KuboChern, get_metric_measures, hexagon, plot_metric_data

##### Module for spectral response functions
include("SpectralResponse.jl")
using .SpecResponse
export spectral, response

end
