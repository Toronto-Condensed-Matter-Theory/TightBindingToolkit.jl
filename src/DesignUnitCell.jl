module DesignUCell
    export AddAnisotropicBond! , AddIsotropicBonds! , ModifyBonds! , ScaleBonds! , RemoveBonds! , ModifyFields!, ModifyIsotropicFields!, Lookup

    using LinearAlgebra, Logging

    using ..TightBindingToolkit.Useful: GetAllOffsets
    using ..TightBindingToolkit.UCell: Bond, BondRank, UnitCell, IsSameBond


@doc """
```julia
AddAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String )
AddAnisotropicBond!( uc::UnitCell , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Matrix{<:Number} , dist::Float64, label::String )
```
Add a bond with the given attributes to `UnitCell`.
If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.

"""	
	function AddAnisotropicBond!( uc::UnitCell{T} , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Array{<:Number, T} , dist::Float64, label::String ) where {T}

		dims 	=	repeat([uc.localDim], T)
		@assert size(mat) == Tuple(dims) "Given Interaction matrix has the inconsistent dimensions as compared to UnitCell!"
		@assert size(offset) == length(uc.primitives) "Inconsistent offset vector dimension as compared to the UnitCell basis dimension!"
		if base <= length(uc.basis) && target <= length(uc.basis)
			if norm( sum(offset .* uc.primitives) .+ (uc.basis[target] .- uc.basis[base] ) ) ≈ dist
				push!( uc.bonds , Bond( base , target , offset , ComplexF64.(mat) , dist, label ) )
			else
				@warn "Inconsistency in bond with label $label" 
			end
		else
			@warn "One or both of basis sites ($base, $target) have not been added to the UnitCell object."
		end
	end

	function AddAnisotropicBond!( uc::UnitCell{T} , base::Int64 , target::Int64 , offset::Vector{Int64} , mat::Number , dist::Float64, label::String ) where {T}
	
		@assert uc.localDim == 1 "Passing a scalar to a bond is only possible if localDim of UnitCell is 1"
		@assert size(offset) == length(uc.primitives) "Inconsistent offset vector dimension as compared to the UnitCell basis dimension!"

		dims 	=	repeat([uc.localDim], T)
	
		AddAnisotropicBond!( uc, base, target, offset, ComplexF64.(reshape([mat], dims...)), dist, label)
	end


@doc raw"""
```julia
AddIsotropicBonds!( uc::UnitCell{T} , dist::Float64 , mats::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)))
AddIsotropicBonds!( uc::UnitCell{T} , dist::Float64 , mats::Array{<:Number, T} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) )
```
Add a set of "isotropic" bonds, which are the same for each pair of sites at the given distance. 
If given `mat` attribute is a number, it is converted into a 1x1 matrix when entered into the bond.
The input `checkOffsetRange` must be adjusted depending on the input distance. 
The optional input `subs` is meant for isotropic bonds when only a subset of sublattices are involved.

"""	
	function AddIsotropicBonds!( uc::UnitCell{T} , dist::Float64 , mat::Array{<:Number, T} , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis)) ) where {T}

		dims 	=	repeat([uc.localDim], T)
		@assert size(mat) == Tuple(dims) "Interaction matrix has the inconsistent dimensions as compared to UnitCell!"
		@assert size(offset) == length(uc.primitives) "Inconsistent offset vector dimension as compared to the UnitCell basis dimension!"

		offsets 		=	GetAllOffsets(checkOffsetRange, length(uc.primitives))    
	
		for i in subs
			for j in subs
				for offset in offsets
					if norm( sum( offset.*uc.primitives ) + (uc.basis[j] - uc.basis[i] ) ) ≈ dist
						proposal 	=	Bond(i, j, offset, ComplexF64.(mat), dist, label)
						if sum(IsSameBond.( Ref(proposal) , uc.bonds ))==0
							push!( uc.bonds , proposal )
						end
					end
				end
	
			end
		end
	end

	function AddIsotropicBonds!( uc::UnitCell{T} , dist::Float64 , mat::Number , label::String; checkOffsetRange::Int64=1 , subs::Vector{Int64}=collect(1:length(uc.basis))) where {T}

		@assert uc.localDim == 1 "Passing a scalar to a bond is only possible if localDim of UnitCell is 1"
		@assert size(offset) == length(uc.primitives) "Inconsistent offset vector dimension as compared to the UnitCell basis dimension!"

		dims 	        =	repeat([uc.localDim], T)
		
		AddIsotropicBonds!( uc, dist, ComplexF64.(reshape([mat], dims...)), label; checkOffsetRange = checkOffsetRange, subs = subs )
	end


@doc """
```julia
ModifyBonds!(uc::UnitCell, dist::Float64, newMat::Matrix{<:Number})
ModifyBonds!(uc::UnitCell, label::String, newMat::Matrix{<:Number})
```
Modify an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`, to the given bond matrix.

"""
	function ModifyBonds!(uc::UnitCell{T}, dist::Float64, newMat::Array{<:Number, T}) where {T}

        dims 	=	repeat([uc.localDim], T)
        @assert size(newMat) == Tuple(dims) "New entry is incompatible with existing bonds!"

        distances 	=	getfield.(uc.bonds, :dist)
        map(x -> x.mat = ComplexF64.(newMat), uc.bonds[findall(≈(dist), distances)])
    end

    function ModifyBonds!(uc::UnitCell{T}, label::String, newMat::Array{<:Number, T}) where {T}

        dims 	=	repeat([uc.localDim], T)
        @assert size(newMat) == Tuple(dims) "New entry is incompatible with existing bonds!"

        labels 	=	getfield.(uc.bonds, :label)
        map(x -> x.mat = ComplexF64.(newMat), uc.bonds[findall(==(label), labels)])
    end


@doc """
```julia
ScaleBonds!(uc::UnitCell, dist::Float64, scale::Number)
ScaleBonds!(uc::UnitCell, label::String, scale::Number)
```
Scale the matrix of an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`, by the given scaling factor.

"""
	function ScaleBonds!(uc::UnitCell, dist::Float64, scale::Number)
		distances 	=	getfield.(uc.bonds, :dist)
		map(x -> x.mat = scale * x.mat, uc.bonds[findall(≈(dist), distances)])
	end

	function ScaleBonds!(uc::UnitCell, label::String, scale::Number)
		labels 	=	getfield.(uc.bonds, :label)
		map(x -> x.mat = scale * x.mat, uc.bonds[findall(==(label), labels)])
	end


@doc """
```julia
RemoveBonds!(uc::UnitCell, dist::Float64)
ScaleBonds!(uc::UnitCell, label::String)
```
Remove an existing bond in the `UnitCell` with the given `label`, or at a given distance=`dist`.

"""
	function RemoveBonds!(uc::UnitCell, label::String )
		labels 	=	getfield.(uc.bonds, :label)
		deleteat!( uc.bonds , findall(==(label), labels) )
	end

	function RemoveBonds!(uc::UnitCell, dist::Float64 )
		labels 	=	getfield.(uc.bonds, :dist)
		deleteat!( uc.bonds , findall(≈(dist), labels) )
	end


@doc """
```julia
ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})
ModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})
```
Modify the on-site fields in the `UnitCell`, either one at a time, or all of them.

"""
	function ModifyFields!(uc::UnitCell, site::Int64, newField::Vector{Float64})
		uc.fields[site] 	=	newField
	end

	function ModifyFields!(uc::UnitCell, newField::Vector{Float64}, dim::Int64)
		@assert length(newField) == length(uc.basis)
		setindex!.(uc.fields, newField, Ref(dim))
	end

	function ModifyFields!(uc::UnitCell, newField::Vector{Vector{Float64}})
		uc.fields 	=	newField
	end


@doc """
```julia
ModifyIsotropicFields!(uc::UnitCell, newField::Vector{Float64})
ModifyIsotropicFields!(uc::UnitCell, newField::Float64, dim::Int64)
```
Modify the on site field uniformly, on all sublattices. The optional argument `dim` is if you want to only modify one of the 4 elements of on-site fields (3 Zeeman and 1 chemical potential).

"""
	function ModifyIsotropicFields!(uc::UnitCell, newField::Vector{Float64})
		uc.fields 	= 	repeat(newField, length(uc.basis))				
	end
	
	function ModifyIsotropicFields!(uc::UnitCell, newField::Float64, dim::Int64)
		map(x -> x[dim] = newField, uc.fields )			
	end


@doc """
```julia
Lookup(uc::UnitCell) --> Dict
```
Returns a dictionary with keys = (base, target, offset) for bond ∈ `UnitCell` bond list, and the entry being the bond matrix. 
If there are multiple bonds with the same identifier, it adds them up.

"""
	function Lookup(uc::UnitCell{T}) :: Dict where {T}
		lookupTable 	=	Dict()

		for bond in uc.bonds
			identifier 	=	(bond.base, bond.target, bond.offset)

			if haskey(lookupTable, identifier)
				lookupTable[identifier] 	=	lookupTable[identifier] + bond.mat

			elseif haskey(lookupTable, (bond.target, bond.base, -bond.offset))
				##### ///TODO: Fix for arrays, when bond rank >2
				##### TODO : TEST
				flippedIndices 	=	collect(T:-1:1)
				lookupTable[(bond.target, bond.base, -bond.offset)] 	=	lookupTable[(bond.target, bond.base, -bond.offset)] + collect(conj.(permutedims(bond.mat, flippedIndices)))

			else
				lookupTable[identifier] 	=	bond.mat
			end
		end

		return lookupTable
	end

end