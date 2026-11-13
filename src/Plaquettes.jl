module Plaqs

    export Plaquette, PlaquetteParam, AddIsotropicPlaquettes!, GetPlaquetteFlux, GetPlaquetteSites

    using LinearAlgebra, Bijections

    using ..TightBindingToolkit.LatticeStruct: Lattice
    using ..TightBindingToolkit.UCell:UnitCell
    using ..TightBindingToolkit.DesignUCell: Lookup


    mutable struct Plaquette{T}

        loopVectors     ::  Vector{Vector{Float64}}
        sites           ::  Vector{Tuple{Int64, Vector{Int64}}}
        mat             ::  Array{ComplexF64, T}
        label           ::  String

        function Plaquette( loopVectors::Vector{Vector{Float64}}, sites::Vector{Tuple{Int64, Vector{Int64}}}, mat::Array{ComplexF64, T}, label::String ) where {T}

            return new{T}( loopVectors, sites, mat, label )
        end

    end


    mutable struct PlaquetteParam{T, R}

        value        ::  R
        plaquettes   ::  Vector{Plaquette{T}}
        label        ::  String

        function PlaquetteParam( value::R, plaquettes::Vector{Plaquette{T}}, label::String ) where {T, R<:Number}

            return new{T, R}( value, plaquettes, label )
        end

        function PlaquetteParam( value::R, rank::Int64) where {R <: Number}

            return new{rank, R}( value, Plaquette{rank}[], "" )
        end

    end


    function AddAnisotropicPlaquette(param::PlaquetteParam{T, R}, plaq::Plaquette{T}) where {T, R}

        if length(param.plaquettes)==0
            param.label =   plaq.label
        end

        push!(param.plaquettes, plaq)

    end


    function AddIsotropicPlaquettes!(param::PlaquetteParam{T, R}, lattice::Lattice{T2}, loopVectors::Vector{Vector{Float64}}, mat::Array{ComplexF64, T}, label::String ; precision::Int64 = 8) where {T, R<:Number, T2}

        positions   =   collect(image(lattice.positions))

        PrimitiveBasis  =   (reduce(hcat, lattice.uc.primitives))
        InvBasis        =   inv(PrimitiveBasis)

        BCStrength  =   abs.(lattice.uc.BC) .> 1e-8     ##### strength of the boundary condition refers to it being 1 (for a periodic system with arbitrary phase), or 0 for open systems. Tolerance since e.g. abs(cis(2π/3)) = 0.9999999999999999.
        LEffective  =   lattice.size .* BCStrength + .!(BCStrength)    ##### Leff_i = L if system is closed, or 1 if system is open along that dimension

        for (site, basis) in lattice.sites

            sub, offset =   basis

            participants    =   [ (sub, offset) ]
            position        =   lattice.positions[site]

            for vector in loopVectors

                position    =   position + vector

                displacements   =   Ref(InvBasis) .* (positions .- Ref(position))   ##### displacements = (r - r_i) in the basis of the UnitCell primitives
                displacements   =   map.(x -> rem.(x, LEffective, Ref(RoundNearest)), displacements) ####### apply the boundary conditions to the displacements
                displacements   =   Ref(PrimitiveBasis) .* displacements ######## convert the displacements back to the original basis 
                
                val, ind    =   findmin(norm.(displacements))   ##### find the site closest to the original site

                if val < 5 * 10.0^(-precision)
                    newSite =   lattice.positions(positions[ind])
                    push!(participants, lattice.sites[newSite])

                else
                    break
                end

            end

            if length(participants) == length(loopVectors) + 1 && participants[begin] == participants[end]
                AddAnisotropicPlaquette(param, Plaquette(loopVectors, participants[begin:end-1], mat, label))
            end

        end

    end


    function GetPlaquetteSites(param::PlaquetteParam{T, R}, lattice::Lattice{T2}) where {T, R, T2}

        sites = getproperty.(param.plaquettes, :sites)
        return [lattice.sites.(site) for site in sites]

    end


    function GetPlaquetteFlux(plaq::Plaquette, uc::UnitCell{T} ; mat::Matrix{ComplexF64} = Matrix{ComplexF64}(I, uc.localDim, uc.localDim)) :: Float64  where {T}

        lookup      =   Lookup(uc)
        plaqSize    =   length(plaq.sites)

        FluxMat     =   Matrix{ComplexF64}(I, uc.localDim, uc.localDim)

        for (s, site) in enumerate(plaq.sites)

            nextSite    =   plaq.sites[(s % plaqSize) + 1]
            
            base    =   site[begin]
            target  =   nextSite[begin]
            offset  =   nextSite[end] - site[end]

            bond    =   get(lookup, (base, target, offset), zeros(ComplexF64, uc.localDim, uc.localDim))
            FluxMat =   FluxMat * bond

        end

        if norm(FluxMat) > 1e-8
            flux    =   tr(adjoint(mat) * FluxMat) / tr(adjoint(mat) * mat)
            flux    =   angle(flux) / (2 * pi)
            return flux

        else
            return 0.0
        end

    end


    function GetPlaquetteFlux(plaq::Plaquette, lattice::Lattice{T} ; mat::Matrix{ComplexF64} = Matrix{ComplexF64}(I, lattice.uc.localDim, lattice.uc.localDim)) :: Float64  where {T}

        plaqSize    =   length(plaq.sites)

        FluxMat     =   Matrix{ComplexF64}(I, lattice.uc.localDim, lattice.uc.localDim)

        for (s, site) in enumerate(plaq.sites)

            nextSite    =   plaq.sites[(s % plaqSize) + 1]
            
            base    =   lattice.sites(site)
            target  =   lattice.sites(nextSite)

            indsOut =   findall(==(target), lattice.bondSites[base, :])
            indsIn  =   findall(==(base), lattice.bondSites[target, :])

            bondMat     =   zeros(ComplexF64, lattice.uc.localDim, lattice.uc.localDim)
            if !isempty(indsOut)
                bondMat     +=   sum(lattice.bondMats[base, indsOut])
                
            elseif !isempty(indsIn)
                bondMat     +=   adjoint(sum(lattice.bondMats[target, indsIn]))
            end

            FluxMat =   FluxMat * bondMat

        end

        if norm(FluxMat) > 1e-8
            flux    =   tr(adjoint(mat) * FluxMat) / tr(adjoint(mat) * mat)
            flux    =   angle(flux) / (2 * pi)
            return flux

        else
            return 0.0
        end

    end
























end