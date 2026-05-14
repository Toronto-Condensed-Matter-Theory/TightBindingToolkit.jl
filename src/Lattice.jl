module LatticeStruct

    export Lattice, FillSites!, FillBonds!, FillLattice!, GetBCPhase, ApplyBCToSite, GetMaxCoordinationNumber

    using LinearAlgebra, Bijections

    using ..TightBindingToolkit.Useful: Meshgrid
    using ..TightBindingToolkit.UCell: Bond, UnitCell


@doc """
```julia
GetMaxCoordinationNumber(uc::UnitCell) --> Int64
```
Returns the maximum coordination number amongs all the sublattices of the `UnitCell`.

"""
    function GetMaxCoordinationNumber(uc::UnitCell) :: Int64

        coords  =   zeros(Int64, length(uc.basis))
        for bond in uc.bonds

            coords[bond.base]    +=  1           
        end

        return max(coords...)
    end

    
@doc """
`Lattice{T}` is a data type representing a general real-space lattice constructed out of a `UnitCell{T}`.

# Attributes
- `uc          ::  UnitCell{T}`: Unit Cell of the lattice.
- `size        ::  Vector{Int64}`: size of the lattice in units of the UnitCell `primitives`.
- `length      ::  Int64`: total number of sites.
- `sites       ::  Dict{ Tuple{Int64, Vector{Int64}}, Int64}` : a dictionary with key = (sublattice, offset), and the corresponding value being the site number on the lattice.
- `positions   ::  Vector{ Vector{Float64}}`: real-space positions of all the lattice sites.
- `fields      ::  Vector{ Vector{Float64}}`: the fields on all the lattice sites.
- `bondSites   ::  Matrix{Int64}`: a matrix with `lattice.length` rows such that `bondSites[s, i]` is the site number of the ith neighbour of site-s on the lattice.
- `bondDists   ::  Matrix{Float64}`: a matrix with `lattice.length` rows such that `bondDists[s, i]` is the distance to the ith neighbour of site-s on the lattice.
- `bondLabels  ::  Matrix{String}`: a matrix with `lattice.length` rows such that `bondLabels[s, i]` is the label of the bond to the ith neighbour of site-s on the lattice.
- `bondMats    ::  Matrix{Array{ComplexF64, T}}`: a matrix with `lattice.length` rows such that `bondMats[s, i]` is the `Array{ComplexF64, T}` of the bond connecting to the ith neighbour of site-s on the lattice.

Initialize this structure using 
```julia
Lattice( uc::UnitCell{T}, size::Vector{Int64} ; null_dist::Float64 = -1.0, null_label::String = "-")
```

where `null_dist` and `null_label` are used for bonds which are not allowed due to given boundary conditions, but still tracked in the code.
"""
    mutable struct Lattice{T}

        uc          ::  UnitCell{T}
        size        ::  Vector{Int64}
        length      ::  Int64
        """
        All sites and fields
        """
        sites       ::  Bijection{Int64, Tuple{Int64, Vector{Int64}}}
        positions   ::  Bijection{Int64, Vector{Float64}}
        fields      ::  Vector{ Vector{Float64}}
        """
        Bond information
        """
        bondSites   ::  Matrix{Int64}
        bondDists   ::  Matrix{Float64}
        bondLabels  ::  Matrix{String}
        bondMats    ::  Matrix{Array{ComplexF64, T}}
        bondShifts  ::  Matrix{Vector{Int64}}

        function Lattice( uc::UnitCell{T}, size::Vector{Int64} ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}

            Ncoords         =   GetMaxCoordinationNumber(uc)
            L               =   length(uc.basis) * prod(size)

            bondSites       =   zeros( Int64, L, Ncoords)
            bondDists       =   null_dist .* ones( Float64, L, Ncoords)
            bondLabels      =   repeat([null_label], L, Ncoords)
            bondMats        =   repeat([zeros(ComplexF64, repeat([uc.localDim], T)...)], L, Ncoords)
            bondShifts      =   repeat([zeros(Int64, length(uc.primitives))], L, Ncoords)

            return new{T}( uc, size, L, Bijection{ Int64, Tuple{Int64, Vector{Int64}}}(), Bijection{ Int64, Vector{Float64}}(), Vector{Float64}[], bondSites, bondDists, bondLabels, bondMats, bondShifts)
        end
    end


@doc """
```julia
FillSites!(lattice::Lattice{T})
```
Fills all the sites, positions, and fields of the lattice using information in the `UnitCell`.

"""
    function FillSites!(lattice::Lattice{T} ; precision::Int64 = 8) where {T}

        offsets     =   collect.(Meshgrid(lattice.size .- Ref(1) ; starts = zeros(Int64, length(lattice.size))))
        site        =   1

        for offset in offsets
            for (b, basis) in enumerate(lattice.uc.basis)

                position    =   sum(offset .* lattice.uc.primitives) + basis

                lattice.sites[site]  =   (b, offset)
                lattice.positions[site]     =   round.(position, digits = precision)
                site        =   site + 1

                push!(lattice.fields, lattice.uc.fields[b])

            end
        end

        lattice.positions[0]    =   round.(999 * ones(Float64, length(lattice.uc.primitives)), digits = precision)
        lattice.sites[0]        =   (0, zeros(Int64, length(lattice.uc.primitives)))
    
    end

@doc """
```julia
ApplyBCToSite(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) -->  Vector{Int64}
```
Returns the effective real-space offset of `site` given a lattice size `L`, and boundary conditions `BC`.

"""
    function ApplyBCToSite(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) :: Vector{Int64}

        BCStrength  =   Bool.(abs.(BC))     ##### strength of the boundary condition refers to it being 1 (for a periodic system with arbitrary phase), or 0 for open systems.
        BCPhase     =   angle.(BC)

        LEffective  =   L .* BCStrength + .!(BCStrength)    ##### Leff_i = L if system is closed, or 1 if system is open along that dimension
        return mod.(site, LEffective) + (.!(BCStrength)) .* (site)  ##### returns site_i mod L_i if system is closed, or site_i if system is open along that dimension.
    end


@doc """
```julia
GetBCPhase(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) --> ComplexF64
```
Returns the effective phase of `site` given a lattice size `L`, and boundary conditions `BC`.

"""
    function GetBCPhase(site::Vector{Int64}, L::Vector{Int64}, BC::Vector{ComplexF64}) :: ComplexF64

        BCPhase     =   angle.(BC)  ##### the phase, ϕ_i corresponding to the boundary condition along each dimension
        phase       =   sum((div.(site, L, Ref(RoundDown))) .* BCPhase)     ##### a vector of [θ_i] s.t. θ_i = n_i * ϕ_i, where n_i  = site_i ÷ L_i
        return exp(im * phase)

    end


@doc """
```julia
FillBonds!(lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-" ) where {T}
```
Fills all the bond information in the lattice using the bonds in `UnitCell`.

"""
    function FillBonds!(lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-" ) where {T}

        bases           =   getproperty.(lattice.uc.bonds, :base)
        coord           =   ones(Int64, lattice.length)

        lookup          =   inv(lattice.sites)

        flippedIndices 	=	collect(T:-1:1)

        for (s, site) in lattice.sites

            sub, offset     =   site

            OutgoingBonds   =   lattice.uc.bonds[findall(==(sub), bases)]

            for bond in OutgoingBonds

                targetOffset    =   ApplyBCToSite(offset + bond.offset, lattice.size, lattice.uc.BC)
                targetPhase     =   GetBCPhase(   offset + bond.offset, lattice.size, lattice.uc.BC)

                target          =   get(lookup, (bond.target, targetOffset), 0)

                lattice.bondSites[ s, coord[s]]    =   target
                lattice.bondDists[ s, coord[s]]    =   !iszero(target) * bond.dist + iszero(target) * null_dist
                lattice.bondLabels[s, coord[s]]    =   iszero(target) ? null_label : bond.label
                lattice.bondMats[  s, coord[s]]    =   !iszero(target) * targetPhase * bond.mat
                lattice.bondShifts[s, coord[s]]    =   div.((offset + bond.offset), lattice.size, Ref(RoundDown))
                
                coord[s]   =   coord[s] + 1

            end

        end
    end


@doc """
```julia
FillLattice!( lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-") where {T}
```
Wrapper function to fill both site and bond information in the lattice.

"""
    function FillLattice!( lattice::Lattice{T} ; null_dist::Float64 = -1.0, null_label::String = "-", precision::Int64 = 8) where {T}

        FillSites!(lattice ; precision = precision)
        FillBonds!(lattice ; null_dist = null_dist, null_label = null_label)

    end





end