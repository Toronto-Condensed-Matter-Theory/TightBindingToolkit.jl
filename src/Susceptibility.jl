module suscep
    export Susceptibility , FindChi , FillChis!, FillRPAChis!, PerformRPA!, FindReducedChi, nF_factor, ResolvedQInteraction, ReduceChis!

    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.Useful: DistFunction, SwitchKroneckerBasis, ImpIndices
    using ..TightBindingToolkit.BZone:BZ, GetQIndex
    using ..TightBindingToolkit.TBModel:Model
    using ..TightBindingToolkit.Parameters:Param
    using ..TightBindingToolkit.UCell:UnitCell


    using LinearAlgebra, TensorCast, Tullio, Logging

    directions = ["x" , "y" , "z"]


    function ResolvedQInteraction(Q::Vector{Float64}, interaction::Param{T, Float64}, uc::UnitCell{T2})::Matrix{Array{ComplexF64, T}} where {T, T2}

        resolvedInt     =   repeat([zeros(ComplexF64, size(interaction.unitBonds[begin].mat))], length(uc.basis), length(uc.basis))     ##### Initializing the sublattice resolved interaction matrix
        flippedIndices 	=	collect(T:-1:1)


        for bond in interaction.unitBonds

            contribution    =   interaction.value[end] * bond.mat * exp.(im .* dot(Q, sum(bond.offset .* uc.primitives)))     ##### Adding the bond contribution to the resolved interaction matrix
            contribution    +=  interaction.value[end] * collect(conj.(permutedims(bond.mat, flippedIndices))) * exp.(-im .* dot(Q, sum(bond.offset .* uc.primitives)))    ##### Adding the hermitian conjugate bond contribution to the resolved interaction matrix

            resolvedInt[bond.base, bond.target]     +=  contribution / 2
            resolvedInt[bond.target, bond.base]     +=  contribution / 2
        end

        return resolvedInt
    end





@doc """
`Susceptibility` is a data type representing the magnetic response, ``χ^{ab}(Q , Ω)`` for a general tight-binding `Model`.

# Attributes
- `M       ::  Model `: the given model.
- `Qs      ::  Vector{Vector{Float64}}`: the set of momentum points over which ``χ^{ab}(Q , Ω)`` is calculated.
- `Omegas  ::  Vector{Float64}`: the set of energies over which ``χ^{ab}(Q , Ω)`` is calculated.
- `Spread  ::  Float64` : the finite spread when summing over delta functions.
- `chis    ::  Dict`: a dictionary containing ``χ^{ab}(Q , Ω)`` for the different directions e.g. `chis["xx"]` etc.

Initialize this structure using 
```julia
Susceptibility(M::Model , Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, [], Omegas, eta, Dict())
Susceptibility(M::Model , Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64} ;  eta::Float64 = 1e-2) = new{}(M, Qs, Omegas, eta, Dict())
```
"""
    mutable struct Susceptibility
        
        Qs      ::  Vector{Vector{Float64}}
        Omegas  ::  Vector{Float64}
        Spread  ::  Float64
        Bare        ::  Matrix{Matrix{ComplexF64}}  ##### The bare susceptibility
        BareReduced ::  Matrix{Matrix{ComplexF64}}  ##### The reduced bare susceptibility
        RPA         ::  Matrix{Matrix{ComplexF64}}  ##### The RPA susceptibility
        RPAReduced  ::  Matrix{Matrix{ComplexF64}}  ##### The reduced RPA susceptibility

        function Susceptibility(Omegas::Vector{Float64} , M::Model;  fill_bz::Bool = true , eta::Float64 = 1e-2)

            if fill_bz
                Qs      =   reshape(M.bz.ks, length(M.bz.ks))
            else
                Qs      =   [zeros(Float64, length(M.uc.primitives))]
            end

            Bare    =   [zeros(ComplexF64, length(M.uc.basis) * length(M.uc.OnSiteMats), length(M.uc.basis) * length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]
            BareReduced     =   [zeros(ComplexF64, length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]

            RPA    =   [zeros(ComplexF64, length(M.uc.basis) * length(M.uc.OnSiteMats), length(M.uc.basis) * length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]
            RPAReduced  =   [zeros(ComplexF64, length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]

            return new{}(Qs, Omegas, eta, Bare, BareReduced, RPA, RPAReduced)
        end

        function Susceptibility(Qs::Vector{Vector{Float64}}, Omegas::Vector{Float64},  M::Model;  eta::Float64 = 1e-2) 
            
            Bare    =   [zeros(ComplexF64, length(M.uc.basis) * length(M.uc.OnSiteMats), length(M.uc.basis) * length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]
            BareReduced     =   [zeros(ComplexF64, length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]

            RPA    =   [zeros(ComplexF64, length(M.uc.basis) * length(M.uc.OnSiteMats), length(M.uc.basis) * length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]
            RPAReduced  =   [zeros(ComplexF64, length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)) for i in 1:length(Omegas), j in 1:length(Qs)]

            return new{}(Qs, Omegas, eta, Bare, BareReduced, RPA, RPAReduced)
        end
    end


@doc raw"""
```julia
nF_factor(Q_index::Vector{Int64}, Omega::Float64 , M::Model; eta::Float64=1e-2) --> Matrix{Matrix{ComplexF64}}
```
function to calculate ``(nF(E(k)) - nF(E(k+Q))) / (Ω + i * η - (E(k+Q) - E(k)))`` for a general multi-band system, at all k-points.

"""
    function nF_factor(Q_index::Vector{Int64}, Omega::Float64 , M::Model; eta::Float64=1e-2) :: Vector{Matrix{ComplexF64}}
        
        Es_k    =   M.Ham.bands
        Es_kpQ  =   circshift(Es_k, -Q_index)  ##### The bands with momentum indices rolled to k+Q instead of k.
        dEs     =   reshape(map((x, y) -> x .- y', Es_kpQ, Es_k), length(M.Ham.bands))  ##### The energy differences between the bands at k+Q and k. Flattened momentum index.

        ##### The fermi distribution functions for each band at each k-point
        nf_k    =   DistFunction.(Es_k ; T=M.T, mu=M.mu, stat=M.stat) 
        nf_kpQ  =   circshift(nf_k, -Q_index)
        dnFs    =   reshape(map((x, y) -> x .- y', nf_k, nf_kpQ), length(M.Ham.bands))  ##### The fermi distribution function differences between the bands at k+Q and k. Flatened momentum index.

        @cast nF[k][i, j] |= dnFs[k][i, j] / (Omega + im * eta - dEs[k][j, i]);
        return nF
    end 


##### Grid index of the Γ point. The Monkhorst grid puts Γ at N/2 + 1 for even N and (N + 1)/2 for odd N,
##### and nowhere at all for a shifted or twisted grid, so it is looked up rather than assumed.
    function GammaIndex(bz::BZ) :: Vector{Int64}
        return GetQIndex(zeros(Float64, length(bz.gridSize)), bz ; nearest = false)
    end


@doc """
```julia
FindReducedChi(Omega::Float64, Q::Vector{Float64}, M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2, Gamma_index::Vector{Int64}=GammaIndex(M.bz)) --> ComplexF64
```
function to calculate susceptibility at a fixed Ω=`Omega`, and `Q`, and along a fixed direction given by `a` and `b`.
`Q` must lie on the grid of `M.bz`, and `Gamma_index` is the grid index of the Γ point, found from the grid by default.

"""
    function FindReducedChi(Omega::Float64, Q::Vector{Float64}, M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2, Gamma_index::Vector{Int64} = GammaIndex(M.bz)) ::ComplexF64

        expQ    =   diagm(exp.(-im .* dot.(Ref(Q) , M.uc.basis)))
        Qa      =   kron(expQ , M.uc.OnSiteMats[a])
        Qb      =   kron(expQ , M.uc.OnSiteMats[b])

        Q_index =   GetQIndex(Q, M.bz) .- Gamma_index
        nF      =   nF_factor(Q_index, Omega , M; eta=eta)

        Uk      =   M.Ham.states
        UkQ     =   circshift(M.Ham.states, -Q_index)

        if a!=b
            Ma      =   reshape(adjoint.(UkQ) .* Ref(Qa) .* Uk, length(M.Ham.bands))
            Mb      =   reshape(adjoint.(UkQ) .* Ref(Qb) .* Uk, length(M.Ham.bands))
            @tullio chi     :=   - Ma[k][i, j] * conj(Mb[k][i, j]) * nF[k][j, i]
        else
            Ma      =   reshape(adjoint.(UkQ) .* Ref(Qa) .* Uk, length(M.Ham.bands))
            @tullio chi     :=   - Ma[k][i, j] * conj(Ma[k][i, j]) * nF[k][j, i]
        end

        return chi 
    end


    function FindChi(Omega::Float64, Q::Vector{Float64}, M::Model; a::Int64=3, b::Int64=3, eta::Float64=1e-2, Gamma_index::Vector{Int64} = GammaIndex(M.bz)) ::Matrix{ComplexF64}

        Q_index =   GetQIndex(Q, M.bz) .- Gamma_index
        nF      =   nF_factor(Q_index, Omega , M; eta=eta)

        Uk      =   M.Ham.states
        UkQ     =   circshift(M.Ham.states, -Q_index)

        chiMat     =   zeros(ComplexF64, length(M.uc.basis), length(M.uc.basis))

        for (i, base) in enumerate(M.uc.basis)
            for (j, target) in enumerate(M.uc.basis)

                Uis     =   selectdim.(Uk, 1, Ref((i - 1) * M.uc.localDim + 1 : i * M.uc.localDim))
                Ujs     =   selectdim.(Uk, 1, Ref((j - 1) * M.uc.localDim + 1 : j * M.uc.localDim))

                UisQ    =   selectdim.(adjoint.(UkQ), 2, Ref((i - 1) * M.uc.localDim + 1 : i * M.uc.localDim))
                UjsQ    =   selectdim.(adjoint.(UkQ), 2, Ref((j - 1) * M.uc.localDim + 1 : j * M.uc.localDim))

                Ma      =   reshape(UisQ .* Ref(M.uc.OnSiteMats[a]) .* Uis, length(M.Ham.bands))
                Mb      =   reshape(UjsQ .* Ref(M.uc.OnSiteMats[b]) .* Ujs, length(M.Ham.bands))

                @tullio chi     :=   - Ma[k][i, j] * conj(Mb[k][i, j]) * nF[k][j, i]
                chiMat[i, j]    +=  chi
            end
        end

        return chiMat 
    end


    function PerformRPA!(chi::Susceptibility, GeneratorInteractions::Vector{Param{2, Float64}}, uc::UnitCell{T}) where {T}
        ##### GeneratorInteraction refers to the interaction written in the basis of all the generators of the local Hilbert space (including the identity matrix). Eg. Id + SU(2) spins -> GeneratorInteraction is a 4x4 matrix.
        resolvedInts    =   sum([ResolvedQInteraction.(chi.Qs, Ref(interaction), Ref(uc)) for interaction in GeneratorInteractions])
        resolvedInts    =   map(x -> hvcat(length(uc.basis), permutedims(x, [2, 1])...), resolvedInts)  ##### Reshaping the resolved interactions to be a vector of matrices, with the first index being the momentum index, and the second being the sublattice x Generator index.
        
        permutation     =   SwitchKroneckerBasis((length(uc.OnSiteMats), length(uc.basis)))
        map!(x -> x[permutation, permutation], resolvedInts, resolvedInts)  ##### Switching the basis of the resolved interactions to be a vector of matrices into the generator x sublattice index
        
        resolvedInts    =   permutedims(reshape(repeat(resolvedInts, length(chi.Omegas)), (length(chi.Qs), length(chi.Omegas))), [2, 1])    ##### Repeating the resolved interactions to be a matrix of matrices, with the outermatrix having the (energy, momentum) index, and the inner being sublattice x Generator index.
        
        corrections     =   inv.(Ref(I) .+ (resolvedInts .* chi.Bare))  ##### The RPA corrections to the bare susceptibility coming from resummation
        chi.RPA         =   chi.Bare .* corrections  ##### The RPA susceptibility
    end

@doc """
```julia
FillChis!(chi::Susceptibility; fill_BZ::Bool=false, a::Int64=3, b::Int64=3)
```
function to calculate susceptibility at a all given Ω=`Omegas`, but for all `Q` present in the given path, and along a fixed direction given by `a` and `b`.

"""
    function FillChis!(chi::Susceptibility, M::Model; a::Int64=3, b::Int64=3, reduce::Bool = true)
        GammaIndex(M.bz)    ##### throws early if the Γ point is not on the BZ grid

        if reduce   ##### Calculate only the reduced susceptibility
            chis    =   FindReducedChi.(chi.Omegas, reshape(chi.Qs, 1, length(chi.Qs)), Ref(M), ; a=a, b=b, eta=chi.Spread)
            chis    =   chis ./ length(M.bz.ks)

            setindex!.(chi.BareReduced, chis, Ref(CartesianIndex(a, b)))
            @info "sublattice reduced bare Chis filled along $(a)-$(b)."

        else    ##### Calculate the sublattice resolved susceptibility
            chis    =   FindChi.(chi.Omegas, reshape(chi.Qs, 1, length(chi.Qs)), Ref(M), ; a=a, b=b, eta=chi.Spread)
            chis    =   chis ./ length(M.bz.ks)

            inds    =   Ref(length(M.uc.basis) * (a - 1)) .+ collect(1:length(M.uc.basis))
            ##### Setting the appropriate indices of the susceptibility matrix which is generator x sublattice.
            for i in 1:length(chi.Omegas), j in 1:length(chi.Qs)
                chi.Bare[i, j][inds, inds]    =   chis[i, j]
            end

            @info "bare Chis filled along $(a)-$(b)."
        end
    
    end


    function FillRPAChis!(chi::Susceptibility, M::Model, GeneratorInteractions::Vector{Param{2, Float64}} ; Generators::Vector{Vector{Int64}} = [[i, i] for i in 1:length(M.uc.OnSiteMats)])

        for generator in Generators
            FillChis!(chi, M, ; a=generator[1], b=generator[2], reduce=false)
        end

        PerformRPA!(chi, GeneratorInteractions, M.uc)

        @info "RPA performed!."
    end


    function ReduceChis!(chi::Susceptibility, M::Model)

        subPhases   =   [zeros(ComplexF64, length(M.uc.basis) * length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)) for i in 1:length(chi.Omegas), j in 1:length(chi.Qs)]
        
        for i in 1:length(chi.Omegas), j in 1:length(chi.Qs)
            subPhases[i, j]     =   kron(Matrix{ComplexF64}(I, length(M.uc.OnSiteMats), length(M.uc.OnSiteMats)), exp.(im .* dot.(Ref(chi.Qs[j]) , M.uc.basis)))
        end

        chi.BareReduced   =   adjoint.(subPhases) .* chi.Bare .* subPhases
        chi.RPAReduced    =   adjoint.(subPhases) .* chi.RPA .* subPhases
    end


end


