module TBModel

    export Model , FindFilling , GetMu! , GetFilling! , GetCount , GetGk! , GetGr!, SolveModel!, GetGap!, FreeEnergy, GetOrderParameter, FindEntropy, GetRSEnergy, GetBondCoorelation

    using ..TightBindingToolkit.Useful: Meshgrid, DistFunction, DeriDistFunction, FFTArrayofMatrix, BinarySearch
    using ..TightBindingToolkit.UCell: UnitCell, Bond
    using ..TightBindingToolkit.BZone: BZ, MomentumPhaseFFT
    using ..TightBindingToolkit.Parameters: Param
    using ..TightBindingToolkit.Hams:Hamiltonian

    using LinearAlgebra, Logging, Statistics


@doc """
`Model` is a data type representing a general Tight Binding system.

# Attributes
- `uc      ::  UnitCell`: the Unit cell of the lattice.
- `bz      ::  BZ`: The discretized Brillouin Zone.
- `Ham     ::  Hamiltonian`: the Hamiltonian at all momentum-points.
- `T       ::  Float64`: the temperature of the system.
- `filling ::  Float64`: The filling of the system.
- `mu      ::  Float64`: The chemical potential of the system.
- `stat    ::  Int64` : ±1 for bosons and fermions.
- `gap     ::  Float64` : the energy gap of excitations at the given filling.
- `Gk      ::  Array{Matrix{ComplexF64}}` : An array (corresponding to the grid of k-points in `BZ`) of Greens functions.

Initialize this structure using 
```julia
Model(uc::UnitCell, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)
```
You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
"""
    mutable struct Model
        uc      ::  UnitCell{2}
        bz      ::  BZ
        Ham     ::  Hamiltonian
        """
        Thermodynamic properties
        """
        T       ::  Float64         ##### Temperature
        filling ::  Float64         ##### Filling fraction
        mu      ::  Float64         ##### Chemical potential
        stat    ::  Int64           ##### +1 for bosons, -1 for fermions
        gap     ::  Float64
        """
        Correlations
        """
        Gk      ::  Array{Matrix{ComplexF64}}
        Gr      ::  Array{Matrix{ComplexF64}}
        
        Model(uc::UnitCell{2}, bz::BZ, Ham::Hamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1) = new{}(uc, bz, Ham, T, filling, mu, stat, -999.0, Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc.primitives))...))
        ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
    end


@doc """
```julia
FindFilling(bands::Vector{Float64}, mu::Float64, T::Float64, stat::Int64=-1) --> Float64
```
Find filling at given `T`=temperature and `mu`=chemical potential, for a given `bands`.

"""
    function FindFilling( mu::Float64, bands::Vector{Float64}, T::Float64, stat::Int64=-1) :: Float64
        return sum(DistFunction(bands; T=T, mu=mu, stat=stat)) / length(bands)
    end

    function FindFilling( mu::Float64, bands::Array{Vector{Float64}, S}, T::Float64, stat::Int64) :: Float64 where {S}
        return sum(sum(DistFunction.(bands; T=T, mu=mu, stat=stat))) / sum(length.(bands))
    end


@doc """
```julia
GetMu!(M::Model; tol::Float64=0.001)
```
Function to get chemical potential for a given `Model`, within a tolerance.

"""
    function GetMu!(M::Model ;  guess::Float64 = 0.0, mu_tol::Float64 = 0.001, filling_tol::Float64 = 1e-6)
        ##### TODO Get an initial guess and pass to BinarySearch
        if M.T≈0.0
            energies  =   sort(reduce(vcat, M.Ham.bands))
            M.mu      =   (energies[floor(Int64, length(energies)*M.filling)] + energies[floor(Int64, length(energies)*M.filling)+1])/2
        else
            M.mu      =   BinarySearch(M.filling, 10 .* M.Ham.bandwidth, FindFilling, (M.Ham.bands, M.T, M.stat) ; initial = guess, x_tol=mu_tol, target_tol = filling_tol)
            @info "Found chemical potential μ = $(M.mu) for given filling = $(M.filling)."
        end
    end


@doc """
```julia
GetFilling!(M::Model)
```
Find filling for a given `Model`.

"""
    function GetFilling!(M::Model)
        if M.T≈0.0
            energies    =   sort(reduce(vcat, M.Ham.bands))
            M.filling   =   searchsortedlast(energies, M.mu) / length(energies)
        else
            M.filling   =   FindFilling( M.mu, M.Ham.bands, M.T, M.stat)
        end
    end


@doc """
```julia
GetCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) --> Matrix{Float64}
```
Function to return a diagonal matrix whose entries are M[i, i] = θ(-(E^i(k)-μ)) ----> 1 if the energy is below the chemical potential, otherwise 0.

"""
    function GetCount(Es::Vector{Float64}, mu::Float64, T::Float64, stat::Int64) :: Matrix{Float64}
        return diagm(DistFunction(Es; T=T, mu=mu, stat=stat))
    end


@doc """
```julia
GetGk!(M::Model)
```
Finding the equal-time Greens functions in momentum space of a `Model`.

"""
    function GetGk!(M::Model)
        quasiCount 	=	GetCount.(M.Ham.bands, Ref(M.mu), Ref(M.T), Ref(M.stat))   ##### Matrix with 1s and 0s along the diagonal. The count of the quasiparticles at each k point determined by the bands and the chemical potential
        M.Gk        =   transpose.(M.Ham.states .* quasiCount .* adjoint.(M.Ham.states))
    end


@doc """
```julia
GetGr!(M::Model)
```
Finding the equal-time Greens functions in real space of a `Model`.

"""
    function GetGr!(M::Model)
        
        Gr           =  FFTArrayofMatrix(M.Gk)  
        phaseShift   =  MomentumPhaseFFT(M.bz, M.uc)
        M.Gr         =  Gr .* phaseShift   
    end


@doc """
```julia
SolveModel!(M::Model)
```
one-step function to find all the attributes in Model after it has been initialized.

"""
    function SolveModel!(M::Model ; mu_guess::Float64 = M.Ham.bandwidth[1] + M.filling * (M.Ham.bandwidth[2] - M.Ham.bandwidth[1]), get_correlations::Bool = true, get_gap::Bool = false, verbose::Bool = true, mu_tol::Float64 = 1e-3, filling_tol::Float64 = 1e-6)
        @assert M.Ham.is_BdG==false "BdG Hamiltonians should be solved using a BdGModel"
        ##### TODO Pass initial guess of chemical potential to GetMu! if provided by user
        if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
            GetFilling!(M)
        else
            GetMu!(M ; guess = mu_guess, mu_tol = mu_tol, filling_tol = filling_tol)
        end

        if get_gap

            energies  =   sort(reduce(vcat, M.Ham.bands))
            M.gap     =   energies[min(floor(Int64, length(energies)*M.filling) + 1, length(energies))] - energies[floor(Int64, length(energies)*M.filling)]
        end
    
        if get_correlations
            GetGk!(M)
            GetGr!(M)
        end
        
        if verbose
            @info "System Filled!"
        end
    end


@doc """
```julia
GetGap!(M::BdGModel)
```
Calculate the BdG gap of the system.
"""
    function GetGap!(M::Model)
        energies    =   sort(reduce(vcat, M.Ham.bands))
        M.gap       =   energies[min(floor(Int64, length(energies)*M.filling) + 1, length(energies))] - energies[floor(Int64, length(energies)*M.filling)]

    end


@doc """
```julia
FreeEnergy(M::Model; F0::Float64 = 0.0) --> Float64
```
Calculate the free energy of the given `Model`.

"""
    function FreeEnergy(M::Model; F0::Float64 = 0.0) :: Float64

        Es      =   reduce(vcat, M.Ham.bands)
        F       =   log.(1 .+ exp.(-(Es .- M.mu) / (M.T)))
        F       =   -M.T * (sum(F) / length(F))

        return F - F0
    end


    function GetOrderParameter(M::Model, param::Param)

        order   =   Float64[]
        for bond in param.unitBonds

            index       =   mod.((-bond.offset) , M.bz.gridSize) .+ ones(Int64, length(bond.offset)) 
            ##### TODO : the extra - sign in offset is because right now G[r] = <f^{dagger}_0 . f_{-r}> ===> NEED TO FIX THIS
            b1          =   M.uc.localDim * (bond.base   - 1) + 1
            b2          =   M.uc.localDim * (bond.target - 1) + 1
            G           =   M.Gr[index...][b1 : b1 + M.uc.localDim - 1, b2 : b2 + M.uc.localDim - 1]

            decomposition   =   (tr( adjoint(bond.mat) * G) / (tr(adjoint(bond.mat) * bond.mat)))
            push!(order, real(decomposition))
        end

        order   =   sum(order) / length(order)
    end

@doc """
```julia
FindEntropy(M::Model;) --> Float64
```
Calculate the entropy of given Model.

"""
function FindEntropy(M::Model,tol=1e-10)
    sum_entropy = 0.0
    for n in 1:length(M.Ham.bands[1])
        fk = DistFunction.(getindex.(M.Ham.bands,n),T=M.T,mu=M.mu)
        one_m_fk = 1.0 .- fk
        #Check elements of one_m_fk are non-zero and positive and change them to 1.0 in this case
        #Check elements of fk are non-zero and positive and change them to 1.0 in this case
        for i in eachindex(one_m_fk)
            if one_m_fk[i] <= tol
                one_m_fk[i] = 1.0
            end
            if fk[i] <= tol
                fk[i] = 1.0
            end
        end
        partial_entropy = -sum(one_m_fk .* log.(one_m_fk)) - sum(fk .* log.(fk))
        # println("Partial entropy for band $n: ", partial_entropy)
        sum_entropy += partial_entropy
    end
    return sum_entropy/(prod(M.bz.gridSize)*length(M.uc.basis))
end

function GetBondCoorelation(Gr::Array{Matrix{ComplexF64}, T}, base::Int64, target::Int64, offset::Vector{Int64}, uc::UnitCell, bz::BZ) :: Matrix{ComplexF64} where {T}
    index       =   mod.((-offset) , bz.gridSize) .+ ones(Int64, length(offset))
    ##### TODO : the extra - sign in offset is because right now G[r] = <f^{dagger}_0 . f_{-r}> ===> NEED TO FIX THIS
    b1          =   uc.localDim * (base   - 1) + 1
    b2          =   uc.localDim * (target - 1) + 1

    G           =   Gr[index...][b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]
    return G
end

function GetBondCoorelation(Gr::Array{Matrix{ComplexF64}, T}, bond::Bond, uc::UnitCell, bz::BZ) :: Matrix{ComplexF64} where {T}
    index       =   mod.((-bond.offset) , bz.gridSize) .+ ones(Int64, length(bond.offset))
    ##### TODO : the extra - sign in offset is because right now G[r] = <f^{dagger}_0 . f_{-r}> ===> NEED TO FIX THIS
    b1          =   uc.localDim * (bond.base   - 1) + 1
    b2          =   uc.localDim * (bond.target - 1) + 1

    G           =   Gr[index...][b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]
    return G
end

@doc """
```julia
GetMFTEnergy(bdgMFT::BdGMFT{T, R}) --> Float64
```
Returns the total mean-field energy of the BdG model including decomposed interactions.

"""
function GetRSEnergy(model::Model) :: Float64

    Energy             =   0.0
    HoppingLookup      =   Lookup(model.uc_hop)

    for BondKey in keys(HoppingLookup)

        G_ij        =   GetBondCoorelation(model.Gr, BondKey..., model.uc_hop, model.bz)
        t_ij        =   HoppingLookup[BondKey]

        if base == target && offset == [0, 0]
                Energy      +=  sum((t_ij .* G_ij))
                #Here only sum the conjuguate if the pairing/hopping involves bond between different sites, as the UC does not keep track of bond and its conjuguate.
                #If it is between different d.o.f, we dont sum the conjugate because the conjuguate part is encoded in the matrix itself.
        else
                Energy      +=  sum((t_ij .* G_ij) + conj.((t_ij .* G_ij)))
        end
    end
    return real(Energy) / length(model.uc_hop.basis)
end

end