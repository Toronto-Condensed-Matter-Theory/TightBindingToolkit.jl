module Hams
    export Hamiltonian , FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian , FillHamiltonian!, FourierHamiltonian, DiagonalizeHamiltonian! , DOS, ModifyHamiltonianField!, IsBandGapped, GetVelocity!

    using ..TightBindingToolkit.Useful: Central_Diff, Arrayfy, DeArrayfy
    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.UCell:Bond, UnitCell, IsSameUnitCell
    using ..TightBindingToolkit.DesignUCell: ModifyFields!, ModifyIsotropicFields!
    using ..TightBindingToolkit.BZone:BZ

    using LinearAlgebra, TensorCast, Logging

    ##### k-independent parts of the Hamiltonian, computed once instead of at every k-point.
    BondDisplacements(uc::UnitCell{2}) = Vector{Float64}[sum(bond.offset .* uc.primitives) for bond in uc.bonds]
    OnSiteBlocks(uc::UnitCell{2}) = Matrix{ComplexF64}[sum(uc.fields[site] .* uc.OnSiteMats) for site in 1:length(uc.basis)]

    ##### Phase arguments of the bonds. Small dot products are written out: `dot` on two Vector{Float64} calls BLAS, which costs more than the arithmetic.
    @inline function SmallDot(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) :: Float64
        length(a) == length(b) || throw(DimensionMismatch("vectors of lengths $(length(a)) and $(length(b))"))
        s   =   0.0
        @inbounds for i in eachindex(a, b)
            s   +=  a[i] * b[i]
        end
        return s
    end

    ##### k . aⱼ for every primitive vector, computed once per momentum.
    PrimitivePhases(k::AbstractVector{<:Real}, uc::UnitCell) = Float64[SmallDot(k, a) for a in uc.primitives]

    ##### k . (offset . primitives) from the precomputed k . aⱼ, or k . deltas[n] when explicit displacements are given.
    @inline function BondPhaseArgument(ka::Vector{Float64}, bond::Bond, deltas::Nothing, k, n::Int) :: Float64
        length(bond.offset) == length(ka) || throw(DimensionMismatch("bond offset of length $(length(bond.offset)) in $(length(ka)) dimensions"))
        arg     =   0.0
        @inbounds for j in eachindex(ka)
            arg +=  bond.offset[j] * ka[j]
        end
        return arg
    end
    @inline BondPhaseArgument(ka, bond::Bond, deltas::AbstractVector, k, n::Int) :: Float64 = SmallDot(k, deltas[n])

    ##### Validates what the @inbounds loops below index: sites, bond matrices, the output matrix and explicit displacements.
    function CheckBonds(uc::UnitCell{2}, H = nothing, deltas = nothing)
        nSites  =   length(uc.basis)
        dims    =   uc.localDim * nSites
        isnothing(H) || size(H) == (dims, dims) || throw(DimensionMismatch("H has size $(size(H)), the unit cell needs ($dims, $dims)"))
        isnothing(deltas) || length(deltas) == length(uc.bonds) || throw(DimensionMismatch("$(length(deltas)) displacements for $(length(uc.bonds)) bonds"))
        for bond in uc.bonds
            (1 <= bond.base <= nSites && 1 <= bond.target <= nSites) || throw(ArgumentError("bond $(bond.base) => $(bond.target) refers to a site outside the $(nSites)-site basis"))
            size(bond.mat) == (uc.localDim, uc.localDim) || throw(DimensionMismatch("bond matrix of size $(size(bond.mat)) in a unit cell with localDim = $(uc.localDim)"))
        end
    end

    ##### H .-= on-site fields, block by block.
    function AddOnSite!(H::AbstractMatrix{ComplexF64}, uc::UnitCell{2}, onsite::AbstractVector{<:AbstractMatrix})
        L   =   uc.localDim
        for site in 1:length(uc.basis)
            o   =   L * (site - 1)
            @views H[o + 1 : o + L, o + 1 : o + L]  .-=   onsite[site]
        end
        return H
    end

    ##### H += hopping bonds and their hermitian conjugates. `deltas = nothing` computes the phases from the bond offsets.
    function AddHoppings!(H::AbstractMatrix{ComplexF64}, uc::UnitCell{2}, k::AbstractVector{<:Real}, deltas)
        CheckBonds(uc, H, deltas)
        L   =   uc.localDim
        ka  =   PrimitivePhases(k, uc)
        @inbounds for (n, bond) in enumerate(uc.bonds)
            o1  =   L * (bond.base - 1)
            o2  =   L * (bond.target - 1)
            M   =   bond.mat

            if bond.base == bond.target && iszero(bond.offset)
                for j in 1:L, i in 1:L
                    H[o1 + i, o1 + j]   +=  (M[i, j] + conj(M[j, i])) / 2
                end
            else
                phase   =   cis(BondPhaseArgument(ka, bond, deltas, k, n))
                for j in 1:L, i in 1:L
                    v   =   phase * M[i, j]
                    H[o1 + i, o2 + j]   +=  v
                    H[o2 + j, o1 + i]   +=  conj(v)
                end
            end
        end
        return H
    end

    ##### H += pairing bonds (no hermitian conjugate).
    function AddPairings!(H::AbstractMatrix{ComplexF64}, uc::UnitCell{2}, k::AbstractVector{<:Real}, deltas)
        CheckBonds(uc, H, deltas)
        L   =   uc.localDim
        ka  =   PrimitivePhases(k, uc)
        @inbounds for (n, bond) in enumerate(uc.bonds)
            o1  =   L * (bond.base - 1)
            o2  =   L * (bond.target - 1)
            M   =   bond.mat
            phase   =   cis(BondPhaseArgument(ka, bond, deltas, k, n))
            for j in 1:L, i in 1:L
                H[o1 + i, o2 + j]   +=  phase * M[i, j]
            end
        end
        return H
    end

@doc """
```julia
FillHoppingHamiltonian(uc::UnitCell, k::Vector{Float64}) --> Matrix{ComplexF64}
```
Returns the hopping Hamiltonian at momentum point `k`, corresponding to the bonds present in `UnitCell`. `uc.OnSiteMats` are used for the fields, with the convention that the last matrix is the one corresponding to the chemical potential.
To evaluate the same unit cell at many momenta, build a [`FourierHamiltonian`](@ref) once instead.

"""
    function FillHoppingHamiltonian(uc::UnitCell{2}, k::AbstractVector{<:Real}, deltas = nothing, onsite = OnSiteBlocks(uc)) :: Matrix{ComplexF64}
        dims    =   uc.localDim * length(uc.basis)
        H       =   zeros(ComplexF64, dims, dims)
        AddOnSite!(H, uc, onsite)
        AddHoppings!(H, uc, k, deltas)
        return H
    end


@doc """
```julia
FillPairingHamiltonian(uc::UnitCell, k::Vector{Float64}) --> Matrix{ComplexF64}
```
Returns the pairing Hamiltonian at momentum point `k`, corresponding to the bonds present in `UnitCell`.

"""
    function FillPairingHamiltonian(uc::UnitCell{2}, k::AbstractVector{<:Real}, deltas = nothing) :: Matrix{ComplexF64}
        dims    =   uc.localDim * length(uc.basis)
        H       =   zeros(ComplexF64, dims, dims)
        AddPairings!(H, uc, k, deltas)
        return H
    end


@doc """
`FourierHamiltonian` is the k-independent form of a Bloch Hamiltonian,
```math
H(k) = H_0 + \\sum_n e^{i k \\cdot d_n} h_n \\, |r_n\\rangle\\langle c_n| ,
```
with every bond of a `UnitCell` (and its hermitian conjugate) merged into one entry per (row, column, displacement).
Build it once and evaluate it at many momenta with `FillHamiltonian(F, k)` or `FillHamiltonian!(H, F, k)`: each phase is computed once per distinct displacement instead of once per bond.
It is a snapshot: rebuild it after the bonds or fields of the unit cell change.

# Attributes
- `H0            :: Matrix{ComplexF64}`: the k-independent part (on-site fields and bonds within a site).
- `displacements :: Matrix{Float64}`: the distinct bond displacements, one per column.
- `rows, cols    :: Vector{Int64}`: matrix entry of each term.
- `disp          :: Vector{Int64}`: column of `displacements` of each term.
- `vals          :: Vector{ComplexF64}`: amplitude of each term.

Initialize this structure using
```julia
FourierHamiltonian(uc::UnitCell) --> hopping Hamiltonian, as FillHoppingHamiltonian
FourierHamiltonian(uc_hop::UnitCell, uc_pair::UnitCell) --> BdG Hamiltonian, as FillHamiltonian(uc_hop, uc_pair, k)
```
"""
    struct FourierHamiltonian
        H0              ::  Matrix{ComplexF64}
        displacements   ::  Matrix{Float64}
        rows            ::  Vector{Int64}
        cols            ::  Vector{Int64}
        disp            ::  Vector{Int64}
        vals            ::  Vector{ComplexF64}
    end

    ##### Accumulates terms keyed by (row, col, integer offset), so repeated bonds collapse into one entry.
    struct FourierTerms
        H0      ::  Matrix{ComplexF64}
        offsets ::  Dict{Vector{Int64}, Int64}
        terms   ::  Dict{Tuple{Int64, Int64, Int64}, ComplexF64}
    end
    FourierTerms(dims::Int) = FourierTerms(zeros(ComplexF64, dims, dims), Dict{Vector{Int64}, Int64}(), Dict{Tuple{Int64, Int64, Int64}, ComplexF64}())

    function AddTerm!(F::FourierTerms, row::Int, col::Int, offset::Vector{Int64}, val::ComplexF64)
        iszero(val) && return F
        if iszero(offset)
            F.H0[row, col]  +=  val
        else
            d   =   get!(F.offsets, offset, length(F.offsets) + 1)
            key =   (row, col, d)
            F.terms[key]    =   get(F.terms, key, zero(ComplexF64)) + val
        end
        return F
    end

    ##### Iterates over the hopping terms of `uc` exactly as FillHoppingHamiltonian adds them: f(row, col, offset, value).
    function ForEachHopping(f::Function, uc::UnitCell{2})
        CheckBonds(uc)
        L   =   uc.localDim
        for (site, block) in enumerate(OnSiteBlocks(uc))
            o   =   L * (site - 1)
            for j in 1:L, i in 1:L
                f(o + i, o + j, zeros(Int64, length(uc.primitives)), -block[i, j])
            end
        end
        for bond in uc.bonds
            o1  =   L * (bond.base - 1)
            o2  =   L * (bond.target - 1)
            M   =   bond.mat
            if bond.base == bond.target && iszero(bond.offset)
                for j in 1:L, i in 1:L
                    f(o1 + i, o1 + j, bond.offset, (M[i, j] + conj(M[j, i])) / 2)
                end
            else
                minus   =   -bond.offset
                for j in 1:L, i in 1:L
                    f(o1 + i, o2 + j, bond.offset, M[i, j])
                    f(o2 + j, o1 + i, minus, conj(M[i, j]))
                end
            end
        end
    end

    ##### Iterates over the pairing terms of `uc` as FillPairingHamiltonian adds them.
    function ForEachPairing(f::Function, uc::UnitCell{2})
        CheckBonds(uc)
        L   =   uc.localDim
        for bond in uc.bonds
            o1  =   L * (bond.base - 1)
            o2  =   L * (bond.target - 1)
            for j in 1:L, i in 1:L
                f(o1 + i, o2 + j, bond.offset, bond.mat[i, j])
            end
        end
    end

    function FourierHamiltonian(F::FourierTerms, primitives::Vector{Vector{Float64}})
        nDisp   =   length(F.offsets)
        D       =   length(primitives)
        displacements   =   zeros(Float64, D, nDisp)
        for (offset, d) in F.offsets
            displacements[:, d]     .=  sum(offset .* primitives)
        end
        ##### Sorted by displacement, then column, so that evaluation walks memory in order.
        keys_   =   sort!(collect(keys(F.terms)), by = t -> (t[3], t[2], t[1]))
        return FourierHamiltonian(F.H0, displacements, getindex.(keys_, 1), getindex.(keys_, 2), getindex.(keys_, 3), [F.terms[t] for t in keys_])
    end

    function FourierHamiltonian(uc::UnitCell{2})
        F   =   FourierTerms(uc.localDim * length(uc.basis))
        ForEachHopping((r, c, o, v) -> AddTerm!(F, r, c, o, ComplexF64(v)), uc)
        return FourierHamiltonian(F, uc.primitives)
    end

    ##### BdG in the nambu basis: [T(k)  Δ(k) - Δ(-k)ᵀ ; h.c.  -T(-k)ᵀ], matching FillHamiltonian(uc_hop, uc_pair, k).
    function FourierHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2})
        @assert IsSameUnitCell(uc_hop, uc_pair) "Inconsistent unit cells for hopping and pairing!"
        dims    =   uc_hop.localDim * length(uc_hop.basis)
        F       =   FourierTerms(2 * dims)
        ForEachHopping(uc_hop) do r, c, o, v
            AddTerm!(F, r, c, o, ComplexF64(v))                     ##### T(k)
            AddTerm!(F, dims + c, dims + r, -o, -ComplexF64(v))     ##### -T(-k)ᵀ
        end
        ForEachPairing(uc_pair) do r, c, o, v
            v   =   ComplexF64(v)
            AddTerm!(F, r, dims + c, o, v)                          ##### Δ(k)
            AddTerm!(F, c, dims + r, -o, -v)                        ##### -Δ(-k)ᵀ
            AddTerm!(F, dims + c, r, -o, conj(v))                   ##### hermitian conjugates of the two above
            AddTerm!(F, dims + r, c, o, -conj(v))
        end
        return FourierHamiltonian(F, uc_hop.primitives)
    end

    Base.size(F::FourierHamiltonian) = size(F.H0)
    Base.size(F::FourierHamiltonian, d::Int) = size(F.H0, d)

@doc """
```julia
FillHamiltonian!(H::Matrix{ComplexF64}, F::FourierHamiltonian, k::Vector{Float64}) --> H
FillHamiltonian(F::FourierHamiltonian, k::Vector{Float64}) --> Matrix{ComplexF64}
```
Evaluates the Bloch Hamiltonian `F` at momentum `k`, in place or into a new matrix.

"""
    function FillHamiltonian!(H::AbstractMatrix{ComplexF64}, F::FourierHamiltonian, k::AbstractVector{<:Real}, phases::Vector{ComplexF64} = Vector{ComplexF64}(undef, size(F.displacements, 2)))
        size(H) == size(F.H0) || throw(DimensionMismatch("H has size $(size(H)), the Hamiltonian $(size(F.H0))"))
        D, nDisp    =   size(F.displacements)
        (length(k) == D && length(phases) == nDisp) || throw(DimensionMismatch("k has length $(length(k)) in $D dimensions, or phases has length $(length(phases)) for $nDisp displacements"))

        @inbounds for d in 1:nDisp
            arg     =   0.0
            for i in 1:D
                arg +=  k[i] * F.displacements[i, d]
            end
            phases[d]   =   cis(arg)
        end

        copyto!(H, F.H0)
        rows, cols, disp, vals  =   F.rows, F.cols, F.disp, F.vals
        @inbounds for n in eachindex(vals)
            H[rows[n], cols[n]]     +=  phases[disp[n]] * vals[n]
        end
        return H
    end

    FillHamiltonian(F::FourierHamiltonian, k::AbstractVector{<:Real}) :: Matrix{ComplexF64} = FillHamiltonian!(similar(F.H0), F, k)


@doc raw"""
```julia
FillHamiltonian(uc::UnitCell, bz::BZ) --> Matrix{Matrix{ComplexF64}}
FillHamiltonian(uc_hop::UnitCell, uc_pair::UnitCell, bz::BZ) --> Matrix{Matrix{ComplexF64}}
FillHamiltonian(uc_hop::UnitCell, uc_pair::UnitCell, k::Vector{Float64}) --> Matrix{ComplexF64}
```
Returns the full Hamiltonian at all momentum points in `BZ` (or at the single momentum `k`), corresponding to the bonds present in `UnitCell`.

"""
    function FillHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, k::AbstractVector{<:Real}) :: Matrix{ComplexF64}

        @assert IsSameUnitCell(uc_hop, uc_pair) "Inconsistent unit cells for hopping and pairing!"

        return FillHamiltonian(uc_hop, uc_pair, k, nothing, nothing, OnSiteBlocks(uc_hop))
    end

    function FillHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, k::AbstractVector{<:Real}, dHop, dPair, onsite) :: Matrix{ComplexF64}
        dims    =   uc_hop.localDim * length(uc_hop.basis)
        Tk      =   FillHoppingHamiltonian( uc_hop,  k, dHop, onsite)
        Tmk     =   FillHoppingHamiltonian( uc_hop, -k, dHop, onsite)
        Δk      =   FillPairingHamiltonian(uc_pair,  k, dPair)
        Δmk     =   FillPairingHamiltonian(uc_pair, -k, dPair)

        ##### The full BdG Hamiltonian in the nambu basis.
        Hk      =   zeros(ComplexF64, 2 * dims, 2 * dims)
        @inbounds for j in 1:dims, i in 1:dims
            Δ                       =   Δk[i, j] - Δmk[j, i]
            Hk[i, j]                =   Tk[i, j]
            Hk[i, dims + j]         =   Δ
            Hk[dims + j, i]         =   conj(Δ)
            Hk[dims + i, dims + j]  =   -Tmk[j, i]
        end
        return Hk
    end

    function FillHamiltonian(F::FourierHamiltonian, bz::BZ)
        phases  =   Vector{ComplexF64}(undef, size(F.displacements, 2))
        return map(k -> FillHamiltonian!(similar(F.H0), F, k, phases), bz.ks)
    end

    FillHamiltonian(uc::UnitCell{2}, bz::BZ) = FillHamiltonian(FourierHamiltonian(uc), bz)

    function FillHamiltonian(uc_hop::UnitCell{2}, uc_pair::UnitCell{2}, bz::BZ)

        @assert IsSameUnitCell(uc_hop, uc_pair) "Inconsistent unit cells for hopping and pairing!"

        return FillHamiltonian(FourierHamiltonian(uc_hop, uc_pair), bz)
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