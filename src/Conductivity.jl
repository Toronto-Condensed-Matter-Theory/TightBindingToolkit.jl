module conduct
    export Conductivity, SpectralFunction, SpectralContribution, GetConductivity!, GetSpectralFunction!

    using LinearAlgebra, Logging

    using ..TightBindingToolkit.SpinMatrices:SpinMats
    using ..TightBindingToolkit.Useful: DistFunction, DeriDistFunction
    using ..TightBindingToolkit.Hams: GetVelocity!
    using ..TightBindingToolkit.BZone:BZ, GetQIndex
    using ..TightBindingToolkit.TBModel:Model

    directions = ["x" , "y" , "z"]


@doc """
`Conductivity` is a data type to track electrical conductivity of a tight-binding model.

# Attributes
- `M       ::  Model`: the tight-binding model whose conductivity is to be calculated.
- `omegas  ::  Vector{Float64}`: The range of energies to integrate over to get the DC conductivity.
- `spread  ::  Float64`: the disorder/spread in the real-frequency Greens function.
- `spectral::  Vector{Array{Matrix{ComplexF64}}}`: the matrix Spectral function at the given frequencies, and all momentum points.
- `sigma   ::  Dict{String, Float64}`: The dictionary containing the conductivity along the different directios.

Initialize this structure using 
```julia
Conductivity(M::Model ; spread::Float64 = 1e-3)
Conductivity(M::Model , omegas::Vector{Float64} ; spread::Float64 = 1e-3)
```
You can either input a filling, or a chemical potential. The corresponding μ for a given filling, or filling for a given μ is automatically calculated.
"""
    mutable struct Conductivity

        M           ::  Model
        omegas      ::  Vector{Float64}
        spread      ::  Float64
        spectral    ::  Vector{Array{Matrix{ComplexF64}}}
        sigma       ::  Dict{String, Float64}

        function Conductivity(M::Model ; spread::Float64 = 1e-3)

            return new{}(M, [0.0], spread, Array{Matrix{ComplexF64}}[], Dict{String, Float64}())
        end

        function Conductivity(M::Model, omegas::Vector{Float64} ; spread::Float64 = 1e-3)

            return new{}(M, omegas, spread, Array{Matrix{ComplexF64}}[], Dict{String, Float64}())
        end
    end


@doc """
```julia
SpectralFunction(H::Array{Matrix{ComplexF64}, N}, omega::Float64 ;  spread::Float64 = 1e-3):: Array{Matrix{ComplexF64}} where {N}
```
Returns the matrix spectral function given the Hamiltonian, the frequency, and the spread.

"""
    function SpectralFunction(H::Array{Matrix{ComplexF64}, N}, omega::Float64 ;  spread::Float64 = 1e-3):: Array{Matrix{ComplexF64}} where {N}

        if isapprox(abs(omega + im * spread), 0.0, atol=1e-6, rtol=1e-6)
            A       =   repeat([zeros(ComplexF64, size(H[begin])...)], size(H)...)

        else
            G       =   inv.(Ref((omega + im * spread) * I) .- H)
            A       =   G - adjoint.(G)

        end
        
        return A ./ (2 * pi * im)
    end


@doc """
```julia
GetSpectralFunction!(cond::Conductivity)
```
Fills the spectral function in the conductivity class.

"""
    function GetSpectralFunction!(cond::Conductivity)

        cond.spectral   =   SpectralFunction.(Ref(cond.M.Ham.H), cond.omegas;  spread = cond.spread)
    end


@doc """
```julia
SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, H::Array{Matrix{ComplexF64}, N}, omega::Float64;  spread::Float64 = 1e-3, a::Int64 = 1):: ComplexF64 where {N}
SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, A::Array{Matrix{ComplexF64}, N};  a::Int64 = 1):: ComplexF64 where {N}
SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, H::Array{Matrix{ComplexF64}, N}, omega1::Float64, omega2::Float64;  spread::Float64 = 1e-3, a::Int64 = 1, b::Int64 = 2):: ComplexF64 where {N}
SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, A1::Array{Matrix{ComplexF64}, N}, A2::Array{Matrix{ComplexF64}, N} ;  a::Int64 = 1, b::Int64 = 2):: ComplexF64 where {N}
```
The contribution to the conductivity coming from the spectral function for directions `a` (and `b`). Optionally can also pass the frequencies `omega` (and `omega2`).

"""
    function SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, H::Array{Matrix{ComplexF64}, N}, omega::Float64;  spread::Float64 = 1e-3, a::Int64 = 1):: ComplexF64 where {N}

        A       =   SpectralFunction(H, omega;  spread = spread)
        Tab     =   sum(tr.(v[a] .* A .* v[a] .* A)) / length(H)

        return Tab
    end

    function SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, A::Array{Matrix{ComplexF64}, N};  a::Int64 = 1):: ComplexF64 where {N}

        Tab     =   sum(tr.(v[a] .* A .* v[a] .* A)) / length(A)

        return Tab
    end

    function SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, H::Array{Matrix{ComplexF64}, N}, omega1::Float64, omega2::Float64;  spread::Float64 = 1e-3, a::Int64 = 1, b::Int64 = 2):: ComplexF64 where {N}

        A1      =   SpectralFunction(H, omega1;  spread = spread)
        A2      =   SpectralFunction(H, omega2;  spread = spread)

        Tab     =   sum(tr.(v[a] .* A1 .* v[b] .* A2)) / length(H)

        return Tab
    end

    function SpectralContribution(v::Vector{Array{Matrix{ComplexF64}}}, A1::Array{Matrix{ComplexF64}, N}, A2::Array{Matrix{ComplexF64}, N} ;  a::Int64 = 1, b::Int64 = 2):: ComplexF64 where {N}

        Tab     =   sum(tr.(v[a] .* A1 .* v[b] .* A2)) / length(A1)

        return Tab
    end


@doc """
```julia
GetConductivity!(Cond::Conductivity, as::Vector{Int64}  ; get_velocity::Bool = false, get_spectral::Bool = false)
```
Calculates the full DC conductivity. Optional Booleans to calculate the velocity matrix of the Hamiltonian, calculate spectral functions also.

"""
    function GetConductivity!(Cond::Conductivity, as::Vector{Int64}  ; get_velocity::Bool = false, get_spectral::Bool = true)

        if get_velocity
            GetVelocity!(Cond.M.Ham, Cond.M.bz)
        end

        if get_spectral
            GetSpectralFunction!(Cond)
        end

        local DistContribution

        if length(Cond.omegas)>1

            dw = (Cond.omegas[2] - Cond.omegas[1]) 
            if dw >= 5 * Cond.M.T
                @warn "Frequency density, $(round(dw, digits = 4)), should be of the order of temperature = $(Cond.M.T)."
            end
            
            DistContribution    =   DeriDistFunction.(Cond.omegas ; T = Cond.M.T, mu = 0.0, stat = Cond.M.stat)
            @assert !(any(isnan, DistContribution) || any(isinf, DistContribution)) "Temperature might be too low : Derivative of nF is underflowing!"
            
        else
            dw = 1.0
            DistContribution    =   -1.0
        end

        for a in as

            spectralContribution=   real.(SpectralContribution.(Ref(Cond.M.Ham.velocity), Cond.spectral; a = a))
            Cond.sigma["$(directions[a])$(directions[a])"]  =   -dw * sum(DistContribution .* spectralContribution) * pi / length(Cond.M.uc.basis)
        end

    end


end