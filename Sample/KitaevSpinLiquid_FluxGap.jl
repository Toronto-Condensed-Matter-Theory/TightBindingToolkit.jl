using Plots, LinearAlgebra, LaTeXStrings#, TightBindingToolkit
include("../src/TightBindingToolkit.jl")
using .TightBindingToolkit

const a1 = [1/2, sqrt(3)/2]
const a2 = [-1/2, sqrt(3)/2]
##### Complex Fermions made from Majorana Fermions live on the triangular lattice #####
const b1 = [0.0, 0.0]
const NNdistance  = 1.0
#### Lattice hopping and pairing of the complexified Majorana fermions ##### 
function KitaevLattice(N::Int64 ; 
    Kx::Float64 = 1.0, Ky::Float64 = 1.0, Kz::Float64 = 1.0,)
    HoppingUC = UnitCell([a1 , a2], 1, 2)
    AddBasisSite!( HoppingUC , b1 )

    PairingUC = UnitCell([a1 , a2], 1, 2)
    AddBasisSite!( PairingUC , b1 )

    ################# Nearest Neighbour hopping #################
    t1Param = Param(1.0, 2)   
    AddAnisotropicBond!( t1Param, HoppingUC, 1, 1, [-1, 0], Kx, NNdistance, "Hopping")
    AddAnisotropicBond!( t1Param, HoppingUC, 1, 1, [0, -1], Ky, NNdistance, "Hopping")

    t0Param = Param(1.0, 2)
    AddAnisotropicBond!( t0Param, HoppingUC, 1, 1, [0, 0], 2.0*Kz , 0.0, "Hopping on-site")

    p1Param = Param(1.0, 2)   
    AddAnisotropicBond!( p1Param, PairingUC, 1, 1, [-1, 0], 1.0*Kx, NNdistance, "Pairing")
    AddAnisotropicBond!( p1Param, PairingUC, 1, 1, [0, -1], 1.0*Ky, NNdistance, "Pairing")

    HoppingLattice = CreateLattice(HoppingUC, [t1Param, t0Param], [N, N])
    PairingLattice = CreateLattice(PairingUC, [p1Param], [N, N])

    return (HoppingLattice, PairingLattice)
end

Ns = 6:6:48
E0s = zeros(length(Ns))
gaps = zeros(length(Ns))
for (i, N) in enumerate(Ns)
    HoppingLattice, PairingLattice = KitaevLattice(N)
    ##### zero-flux sector #####
    bdgHamiltonian = LatticeHamiltonian(HoppingLattice, PairingLattice)
    DiagonalizeHamiltonian!(bdgHamiltonian)
    E0 = sum(bdgHamiltonian.bands[1:N*N])
    E0s[i] = E0/(N*N)
    ##### two-flux sector #####
    #####! Flipping one of the Z2 gauge fields to create a pair of Z2 fluxes #####
    HoppingLattice.bondMats[1, 3] = -HoppingLattice.bondMats[1, 3]
    bdgHamiltonian = LatticeHamiltonian(HoppingLattice, PairingLattice)
    DiagonalizeHamiltonian!(bdgHamiltonian)
    E1 = sum(bdgHamiltonian.bands[1:N*N])
    gaps[i] = E1 - E0
end

p_ground = plot(
    framestyle = :box,
    xlabel = L"1/L", ylabel = L"E_0",
    guidefont = font(14, "Computer Modern"),
    tickfont = font(14, "Computer Modern"),
)
plot!(p_ground, (1 ./ (Ns .^ 2)), (E0s/8), lw = 4.0, label = "", marker=true, markersize = 6.0)
ylims!(p_ground, -0.401, -0.39)

p_gap = plot(
    framestyle = :box,
    xlabel = L"1/L", ylabel = L"\Delta E_{\mathrm{2-flux}}",
    guidefont = font(14, "Computer Modern"),
    tickfont = font(14, "Computer Modern"),
)
plot!(p_gap, (1 ./ (Ns .^ 2))[1:end], (gaps/8)[1:end], lw = 4.0, label = "", marker=true, markersize = 6.0)
ylims!(p_gap, 0.0, ylims(p_gap)[2]+0.01)