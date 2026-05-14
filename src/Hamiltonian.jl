module Hams
    export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , DiagonalizeHamiltonian! , DOS, ModifyHamiltonianField!, IsBandGapped, GetVelocity!

    using ..TightBindingToolkit.Useful: Central_Diff, Arrayfy, DeArrayfy
    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.UCell:Bond, UnitCell, IsSameUnitCell
    using ..TightBindingToolkit.DesignUCell: ModifyFields!, ModifyIsotropicFields!
    using ..TightBindingToolkit.BZone:BZ

    using LinearAlgebra, TensorCast, Logging

@doc """
```julia
FillHoppingHamiltonian(uc::UnitCell, k::Vector{Float64} ; OnSiteMatrices::Vector{Matrix{ComplexF64}}) --> Matrix{ComplexF64}
```
Returns the hopping Hamiltonian at momentum point `k`, corresponding to the bonds present in `UnitCell`. `OnSiteMatrices` are used for the fields, with the convention that the last matrix is the one corresponding to the chemimcal potential.

"""
    function FillHoppingHamiltonian(uc::UnitCell{2}, k::Vector{Float64})
        dims    =   uc.localDim * length(uc.basis)
        H       =   zeros(ComplexF64, dims, dims)
        
        ##### On-site terms in the Hamiltonian
        for site in 1:length(uc.basis)
            b1  =   uc.localDim * (site - 1) + 1
            b2  =   uc.localDim * (site - 1) + 1
            H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .-=   sum(uc.fields[site] .* uc.OnSiteMats) 
        end
        ##### Inter-site terms which depend on the momentum
        for bond in uc.bonds
            b1  =   uc.localDim * (bond.base - 1) + 1
            b2  =   uc.localDim * (bond.target - 1) + 1
            
            if b1==b2 && bond.offset==zeros(length(uc.primitives))
                H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   (bond.mat + bond.mat') / 2

            else
                H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=    exp( im .* dot(k, sum(bond.offset .* uc.primitives))) .* bond.mat
                H[b2 : b2 + uc.localDim - 1, b1 : b1 + uc.localDim - 1]  .+=    exp(-im .* dot(k, sum(bond.offset .* uc.primitives))) .* bond.mat' 
            end 
        end
    
        return H
    end


@doc """
```julia
FillPairingHamiltonian(uc::UnitCell, k::Vector{Float64}) --> Matrix{ComplexF64}
```
Returns the pairing Hamiltonian at momentum point `k`, corresponding to the bonds present in `UnitCell`.

"""
    function FillPairingHamiltonian(uc::UnitCell{2}, k::Vector{Float64}) :: Matrix{ComplexF64}
        dims    =   uc.localDim * length(uc.basis)
        H       =   zeros(ComplexF64, dims, dims)
    
        for bond in uc.bonds
            b1  =   uc.localDim * (bond.base - 1) + 1
            b2  =   uc.localDim * (bond.target - 1) + 1
            
            H[b1 : b1 + uc.localDim - 1, b2 : b2 + uc.localDim - 1]  .+=   exp( im .* dot(k, sum(bond.offset .* uc.primitives))) .* bond.mat

        end
    
        return H
    end


@doc raw"""
```julia
FullHamiltonian(uc::UnitCell, bz::BZ) --> Matrix{Matrix{ComplexF64}}
```
Returns the full Hamiltonian at all momentum points in `BZ`, corresponding to the bonds present in `UnitCell`.

"""
    function FillHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, k::Vector{Float64}) :: Matrix{ComplexF64}

        @assert IsSameUnitCell(uc_hop, uc_pair) "Inconsistent unit cells for hopping and pairing!"

        Tk      =   FillHoppingHamiltonian( uc_hop,  k)
        Tmk     =   FillHoppingHamiltonian( uc_hop, -k)
        Δk      =   FillPairingHamiltonian(uc_pair,  k ) - transpose(FillPairingHamiltonian(uc_pair, -k ))

        ##### The full BdG Hamiltonian in the nambu basis.
        Hk      =   hcat(Tk , Δk)
        Hk      =   vcat(Hk , hcat(Δk' , -transpose(Tmk)))
        return (Hk)
    end

    function FillHamiltonian(uc::UnitCell{2}, bz::BZ)
        
        return FillHoppingHamiltonian.(Ref(uc), bz.ks)
    end
    
    function FillHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, bz::BZ)

        return FillHamiltonian.(Ref(uc_hop), Ref(uc_pair), bz.ks)
    end


@doc """
`Hamiltonian` is a data type representing a general momentum-space Hamiltonian corresponding to the given `UnitCell` and `BZ` (or 2 Unit Cells if it is a BdG Hamiltonian).

# Attributes
- `H           :: Array{Matrix{ComplexF64}}`: A Array (corresponding to the grid of k-points in `BZ`) of Hamiltonian matrices.
- `bands       :: Array{Vector{Float64}}`: A Array (corresponding to the grid of k-points in `BZ`) of band spectrums.
- `states      :: Array{Matrix{ComplexF64}}`: A Array (corresponding to the grid of k-points in `BZ`) of band wavefunctions.
- `bandwidth   :: Tuple{Float64, Float64}` : the tuple of minimum and maximum energies in the band structure.
- `is_BdG      :: Bool` : is the Hamiltonian a bdG hamiltonian or a pure hopping hamiltonian.

Initialize this structure using
```julia
Hamiltonian(uc::UnitCell, bz::BZ) --> Hopping Hamiltonian
Hamiltonian(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ) --> BdG Hamiltonian
```
"""
    mutable struct Hamiltonian
        H           ::  Array{Matrix{ComplexF64}}
        bands       ::  Array{Vector{Float64}}
        states      ::  Array{Matrix{ComplexF64}}
        bandwidth   ::  Tuple{Float64, Float64}
        is_BdG      ::  Bool
        velocity    ::  Vector{Array{Matrix{ComplexF64}}}
    
        Hamiltonian(uc::UnitCell{2}, bz::BZ) = new{}(FillHamiltonian(uc, bz), Array{Vector{Float64}}(undef, zeros(Int64, length(uc.primitives))...), Array{Matrix{ComplexF64}}(undef, zeros(Int64, length(uc.primitives))...), (0.0, 0.0), false, Array{Matrix{ComplexF64}, length(uc.primitives)}[])

        Hamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, bz::BZ) = new{}(FillHamiltonian(uc_hop, uc_pair, bz), Array{Vector{Float64}}(undef, zeros(Int64, length(uc_hop.primitives))...), Array{Matrix{ComplexF64}}(undef,zeros(Int64, length(uc_hop.primitives))...), (0.0, 0.0), true, Array{Matrix{ComplexF64}, length(uc_hop.primitives)}[])
    end


@doc """
```julia
DiagonalizeHamiltonian!(Ham::Hamiltonian)
```
Diagonalize the `Hamiltonian` at all momentum points in the `BZ`. `verbose` is an optional argument to print when the Hamiltonian is diagonalized.

"""
    function DiagonalizeHamiltonian!(Ham::Hamiltonian ; verbose::Bool=true)
        sols            =   eigen.(Hermitian.(Ham.H))
        Ham.bands       =   getfield.(sols, :values)
        Ham.states      =   getfield.(sols, :vectors)
        Ham.bandwidth   =   (minimum(first.(Ham.bands)) , maximum(last.(Ham.bands)))
        if verbose
            @info "Hamiltonian Diagonalized!"
        end
    end


@doc """
```julia
ModifyHamiltonianField!(Ham::Hamiltonian, uc::UnitCell, newFields::Vector{Float64} ; dim::Int64=4, verbose::Bool=false)
```
Faster implementation of modifying ONLY the on-site field part of a `Hamiltonian`. `newFields` must be a vector of the same length as `uc.basis`.
`dim` is an optional argument which determines which element of on-site field is being replaced.

"""
    function ModifyHamiltonianField!(Ham::Hamiltonian, uc::UnitCell, newFields::Vector{Float64} ; dim::Int64=length(uc.fields[begin]),  verbose::Bool=false)

        @assert length(newFields)==length(uc.basis) "Inconsistent number of basis sites and fields given"
    
        if Ham.is_BdG==false 
            Ham.H   .+=   Ref(kron(diagm(getindex.(uc.fields, dim) .- newFields) , uc.OnSiteMats[dim]))
            DiagonalizeHamiltonian!(Ham ; verbose=verbose)
        else
            Ham.H   .+=   Ref(kron(SpinMats(1//2)[3] * 2, kron(diagm(getindex.(uc.fields, dim) .- newFields) , uc.OnSiteMats[dim])))  ##### The extra Sz corresponds to the nambu basis.
            DiagonalizeHamiltonian!(Ham ; verbose=verbose)
        end
    
        ModifyFields!(uc, newFields, dim)
    end


@doc """
```julia
isBandGapped(H::Hamiltonian ; tol::Float64 = 1e-3) --> BitMatrix
```
Returns a matrix of booleans marked as `true` if the band corresponding to the row and column of the matrix are gapped (greater than the tolerance), and false otherwise. 

"""
    function IsBandGapped(H::Hamiltonian ; tol::Float64 = 1e-3) :: BitMatrix
        bands   =   reshape(H.bands, length(H.bands))
        @cast bandDiff[i, j][k] |= abs(bands[k][i] - bands[k][j])
        @cast gapped[i, j] |= minimum(bandDiff[i, j]) > tol

        return gapped
    end


@doc """
```julia
DOS(Omega::Float64, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Float64
DOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-3) --> Vector{Float64}
```
Calculate the Density of State correspondingto the given energies in `Omegas`, for the lowest bands upto `till_band`.
The calculation is done at a finite `spread` of the delta-function sum. 
"""
    function DOS(Omega::Float64, energies::Vector{Float64}; spread::Float64=1e-2) :: Float64

        dos         =   @. 1 / ((Omega - energies) + im * spread)
        return sum(imag(dos))
    end

    function DOS(Omegas::Vector{Float64}, Ham::Hamiltonian; till_band::Int64=length(Ham.bands[1, 1]), spread::Float64=1e-2) :: Vector{Float64}

        energies    =   reduce(vcat, Ham.bands)
        n_bands     =   length(Ham.bands[1, 1])
        energies    =   energies[filter(i -> (i-1) % n_bands + 1 <= till_band, 1:length(energies))]

        dos     =   DOS.(Omegas, Ref(energies); spread=spread)

        dOmega  =   Omegas[2:end] .- Omegas[1:end-1]
        norm    =   sum(dos[1:end-1] .* dOmega)

        dos     .=  dos ./ norm
        return dos
    end


@doc """
```julia
GetVelocity!(H::Hamiltonian, bz::BZ) --> Vector{Array{Matrix{ComplexF64}, T}}
```
returns a vector of velocity matrices at each k point defined as ∂H(k)/∂k^{μ}.
"""
    function GetVelocity!(H::Hamiltonian, bz::BZ) :: Vector{typeof(H.H)}

        dH      =   Central_Diff(H.H ; delta = (1 ./ bz.gridSize), PBC = repeat([true], length(bz.gridSize)))
        dH      =   Arrayfy.(dH)

        bMatrix =   inv(transpose(hcat(bz.basis...)))   ##### Changing velocity from reciprocal basis to cartesian basis.
        
        v       =   bMatrix * dH
        H.velocity       =   DeArrayfy.(v, Ref(bz.gridSize))

    end






end