module BZone
    export GetRLVs , BZ , Monkhorst , FillBZ! , ReduceQ , GetQIndex , CombinedIndexPath , GetBZPath , CombinedBZPath, MomentumPhaseFFT

    using LinearAlgebra, Logging

    using ..TightBindingToolkit.Useful: Meshgrid, VecAngle, GetIndexPath
    using ..TightBindingToolkit.UCell: UnitCell


@doc """
```julia
GetRLVs( uc::UnitCell ) --> Vector{Vector{Float64}}
```
Returns the reciprocal lattice vectors corresponding to the given Unit Cell.

"""
    function GetRLVs( uc::UnitCell ) :: Vector{Vector{Float64}}
        if length(uc.primitives) == 1
            a1      =   uc.primitives[1]
            return [(2 *pi / dot(a1, a1)) * a1]

        elseif length(uc.primitives) == 2
            d   =   length(uc.primitives[begin])

            if d==2
                a1  = vcat( uc.primitives[1] , 0.0 )
                a2  = vcat( uc.primitives[2] , 0.0 )
            end

            a3  = cross(a1, a2)
            V   = dot( a1 , cross( a2 , a3 ) )
            b1  = 2 * pi / V * ( cross( a2 , a3 ) )
            b2  = 2 * pi / V * ( cross( a3 , a1 ) )
            return [ b1[1:d] , b2[1:d] ]

        elseif length(uc.primitives) == 3
            V   = dot( uc.primitives[1] , cross( uc.primitives[2] , uc.primitives[3] ) )
            b1  = 2 * pi / V * ( cross( uc.primitives[2] , uc.primitives[3] ) )
            b2  = 2 * pi / V * ( cross( uc.primitives[3] , uc.primitives[1] ) )
            b3  = 2 * pi / V * ( cross( uc.primitives[1] , uc.primitives[2] ) )
            return [ b1 , b2 , b3 ]

        else
            @warn "Getting reciprocal lattice vectors only works for up to d=3 lattices right now."
        end
    end

@doc """
`BZ` is a data type representing a discretized Brillouin Zone in momentum space.

# Attributes
- `basis           :: Vector{ Vector{ Float64 } }`: reciprocal lattice vectors of the Brillouin Zone.
- `gridSize        :: Vector{Int64}`: The number of points along each dimension of the grid.
- `kInds           :: Vector{Array{Float64}}`: the [`Monkhorst`](@ref) grid corresponding to k along bs.
- `ks   	       :: Array{Vector{Float64}}`: The grid of momentum vectors k.
- `HighSymPoints   :: Dict`: A dictionary containing the HIgh-Symmetry points Γ, K(2), and M(3).
- `shift           :: Vector{Int64}` : how shifted the grid is from its centre point at the Γ point, in units of `1/gridSize`.

Initialize this structure using 
```julia
BZ(gridSize::Int64) --> defaults to dims = 2
BZ(gridSize::Int64, dims) --> Uniform grid in dimension = dims
BZ(gridSize::Vector{Int64})
```
"""
    mutable struct BZ
        basis	        ::  Vector{Vector{Float64}}
        gridSize        ::  Vector{Int64}
        kInds           ::  Vector{Array{Float64}}
        ks   	        ::	Array{Vector{Float64}}
        HighSymPoints   ::  Dict
        shift           ::  Vector{Int64}
    
        BZ(gridSize::Vector{Int64})         =   new{}(Vector{Float64}[],  gridSize, Array{Float64, length(gridSize)}[], Array{Vector{Float64}}(undef, zeros(Int64, length(gridSize))...), Dict(), Int64[])
        BZ(gridSize::Int64, dims::Int64) 	=	new{}(Vector{Float64}[], repeat([gridSize], dims), Array{Float64, dims}[], Array{Vector{Float64}}(undef, zeros(Int64, dims)...), Dict(), Int64[])
        
        function BZ(gridSize::Int64)
            #@warn "Positional argument `dims' not passed in BZ. Resorting to its default value of 2."
            dims 	=	2
            return   new{}(Vector{Float64}[], repeat([gridSize], dims), Array{Float64, dims}[], Array{Vector{Float64}}(undef, zeros(Int64, dims)...), Dict(), Int64[])
        end
    end


@doc raw"""
```julia
Monkhorst(ind::Int64, N::Int64) --> Float64
Monkhorst(ind::Int64, N::Int64, shift::Int64, BC::Float64) --> Float64
```
The usual Monkhorst grid is defined as follows
`` [\frac{2i - (N+1)}{2N}, i∈[1, N]] ``
The modified Monkhorst grid takes into account the desired boundary condition `BC`, and an integer `shift` (if required, to change the starting point), and shifts the momentum grid accordingly.

"""
    function Monkhorst(ind::Int64, N::Int64) :: Float64
        return (2 * ind - (N + 1)) / (2 * N)
    end

    function Monkhorst(ind::Int64, N::Int64, shift::Int64, BC::Float64) :: Float64
        return (1 / N) * (ind + shift - ((N + 2 - (N % 2)) / 2) + (BC / (2 * pi)))
    end

    
    ##### Adds high symmetry points to the bz dictionary
    function AddHighSymPoints!(bz::BZ)
        dims    =   length(bz.basis)
        
        ##### Gamma point at the origin.
        bz.HighSymPoints["G"]   =   zeros(Float64, dims)
        
        if dims == 1    ##### 1d BZ
            bz.HighSymPoints["M1"]      =   @.  0.5 * bz.basis[1]
            bz.HighSymPoints["-M1"]     =   @. -0.5 * bz.basis[1]    ##### M point at the midpoint of the BZ (pi)
            bz.HighSymPoints["K1"]      =   @. (1/3) * bz.basis[1]  ##### The two K points at pi/3 and 2pi/3
            bz.HighSymPoints["K2"]      =   @. (2/3) * bz.basis[1]  #####
            bz.HighSymPoints["K1"]      =   @. (1/3) * bz.basis[1]  ##### The two K points at pi/3 and 2pi/3
            bz.HighSymPoints["K2"]      =   @. (2/3) * bz.basis[1]

        elseif dims == 2    ##### 2d BZ
    
            bz.HighSymPoints["M1"]      =    @. 0.5 * bz.basis[1] + 0.0 * bz.basis[2]     ##### The 3 possible M points
            bz.HighSymPoints["X"]      =    @. 0.0 * bz.basis[1] + 0.5 * bz.basis[2]
            bz.HighSymPoints["M2"]      =    @. 0.5 * bz.basis[1] + 0.5 * bz.basis[2]
            bz.HighSymPoints["-M1"]     =  -(@. 0.5 * bz.basis[1] + 0.0 * bz.basis[2])     ##### The 3 possible M points
            bz.HighSymPoints["-X"]     =  -(@. 0.0 * bz.basis[1] + 0.5 * bz.basis[2])
            bz.HighSymPoints["-M2"]     =  -(@. 0.5 * bz.basis[1] + 0.5 * bz.basis[2])
    
            if isapprox(VecAngle(bz.basis[1], bz.basis[2]), 2*pi/3, atol=1e-4, rtol=1e-4)   ##### The K points depend on the relative angle of the reciprocal basis.
                bz.HighSymPoints["K1"]      =   @. (2/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2"]      =   @. (1/3) * bz.basis[1] + (2/3) * bz.basis[2]
                bz.HighSymPoints["K1'"]     =   @. (-1/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2'"]     =   @. (-2/3) * bz.basis[1] + (-1/3) * bz.basis[2]
                bz.HighSymPoints["K1''"]    =   @. (-1/3) * bz.basis[1] + (-2/3) * bz.basis[2]
                bz.HighSymPoints["K2''"]    =   @. (1/3) * bz.basis[1] + (-1/3) * bz.basis[2]
    
            elseif isapprox(VecAngle(bz.basis[1], bz.basis[2]), pi/3, atol=1e-4, rtol=1e-4)
                bz.HighSymPoints["K1"]      =   @. (1/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2"]      =   @. (-1/3) * bz.basis[1] + (-1/3) * bz.basis[2]
                bz.HighSymPoints["K1'"]     =   @. (-2/3) * bz.basis[1] + (1/3) * bz.basis[2]
                bz.HighSymPoints["K2'"]     =   @. (2/3) * bz.basis[1] + (-1/3) * bz.basis[2]
                bz.HighSymPoints["K1''"]    =   @. (1/3) * bz.basis[1] + (-2/3) * bz.basis[2]
                bz.HighSymPoints["K2''"]    =   @. (-1/3) * bz.basis[1] + (2/3) * bz.basis[2]
            end
    
        end
    
    end

@doc """
```julia
FillBZ!(bz::BZ, uc::UnitCell ; shift::Vector{Float64}=zeros(Float64, length(uc.primitives)))
```
Fills the `BZ` object with the relevant attributes, after it has been initialized with `BZ(gridSize)`.

"""
    function FillBZ!(bz::BZ, uc::UnitCell; shift::Vector{Int64}=zeros(Int64, length(uc.primitives)))

        bz.basis    =   GetRLVs(uc)   ##### Get the reciprocal basis
        bz.shift    =   shift
        dims 		=	length(uc.primitives)

        @assert length(bz.gridSize) == dims "Inconsistent dimensions of UnitCell and BZ"
        @assert length(bz.basis)<=3 "Sorry, the code is only written for upto 3 dimensions right now!"

        inds    =   Meshgrid(bz.gridSize)

        ##### Getting the Meshgrid of Monkhorst coefficients of the BZ grid
        for dim in 1:dims
            index   =   getindex.(inds, dim)
            push!(bz.kInds, Monkhorst.(index, Ref(bz.gridSize[dim]), Ref(shift[dim]), Ref(angle(uc.BC[dim]))))
        end

        bz.ks   =   (bz.kInds[1] .* Ref(bz.basis[1])) 
        for dim in 2:dims
            bz.ks   .+=     (bz.kInds[dim] .* Ref(bz.basis[dim]))
        end 

        AddHighSymPoints!(bz)

    end


@doc """
```julia
ReduceQ(Q::Vector{Float64}, bz::BZ) --> Vector{Float64}
```
Reduces a given momentum back to the range covered by the discretized Brillouin Zone.

"""
    function ReduceQ(Q::Vector{Float64}, bz::BZ) :: Vector{Float64}

        U           =   reduce(hcat, bz.basis) ##### basis transformation from x, y -> b1, b2
        Q_reduced   =   inv(U) * (Q)  ##### Q in the basis of b1 and b2
        Q_reduced   =   Q_reduced .- round.(Q_reduced .- (bz.shift ./ bz.gridSize))  ##### Q in the basis of b1 and b2, and shifted back to the first BZ
        Q_reduced   =   sum((Q_reduced) .* bz.basis)   ##### Q back in the x ,y basis, but shifted to the first BZ
        return Q_reduced
    end


@doc """
```julia
GetQIndex(Q::Vector{Float64}, bz::BZ ; nearest::Bool = false) --> Vector{Int64}
```
Returns the index in the discretized `BZ` of the momentum point corresponding to the fiven momentum `Q`. 
If the input `nearest` is set to `true`, will return the index of the momentum point on the grid closest to `Q`, if `Q` does not exist on the grid. 

"""
    function GetQIndex(Q::Vector{Float64}, bz::BZ ; nearest::Bool = false) :: Vector{Int64}

        Q_reduced   =   ReduceQ(Q, bz)

        if !nearest
            inds    =   findfirst(isapprox(Q_reduced, rtol=1e-5, atol=1e-5), bz.ks)
            @assert !isnothing(inds)  Q_reduced, "Given momentum does not exist on the BZ grid" 
        else
            inds    =   findmin(norm.(bz.ks .- Ref(Q_reduced)))[2]
        end

        return collect(Tuple(inds))
    end


@doc """
```julia
CombinedIndexPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) --> Vector{Vector{Int64}}
```
Returns a path in index-space joins the given points present in `points` as point[1] --> point[2] --> ... --> point[end] --> point[1].
The optional input `nearest` is the same as in [`GetQIndex`](@ref), and `closed` determines if the path is a closed loop or not.

"""
    function CombinedIndexPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) :: Vector{Vector{Int64}}

        if closed
            starts  =   copy(points)
            endings =   circshift(starts, -1)
        else
            starts  =   copy(points[1:end-1])
            endings =   copy(points[2:end])
        end
    
        paths   =   GetIndexPath.(GetQIndex.(starts, Ref(bz) ; nearest = nearest), GetQIndex.(endings, Ref(bz), nearest=nearest) ; exclusive = true)
        paths   =   reduce(vcat, paths)
    
        return paths
    end


@doc """
```julia
GetBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) --> Vector{Vector{Float64}}
```
Returns the actual path in momentum-space of the discretized `BZ` which joins the two momentums `start` and `ending`. 
The optional input `nearest` is the same as in [`GetQIndex`](@ref), and `exclusive` is the same as in [`GetIndexPath`](@ref).

"""
    function GetBZPath(bz::BZ, start::Vector{Float64}, ending::Vector{Float64} ; nearest::Bool = false, exclusive::Bool = true) :: Vector{Vector{Float64}}

        start_ind   =   GetQIndex(start, bz ; nearest=nearest)
        end_ind     =   GetQIndex(ending, bz ; nearest=nearest)
        path_ind    =   GetIndexPath(start_ind, end_ind ; exclusive = exclusive)

        return getindex.(Ref(bz.ks), CartesianIndex.(Tuple.(path_ind)))
    end


@doc """
```julia
CombinedBZPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) --> Vector{Vector{Float64}}
```
Returns a path in momentum-space of the discretized `BZ` which joins the given momentum points present in `points` as point[1] --> point[2] --> ... --> point[end] --> point[1].
The optional input `nearest` is the same as in [`GetQIndex`](@ref), and `closed` determines if the path is a closed loop or not.

"""
    function CombinedBZPath(bz::BZ, points::Vector{Vector{Float64}} ; nearest::Bool = false, closed::Bool = true) :: Vector{Vector{Float64}}

        if closed
            starts  =   copy(points)
            endings =   circshift(starts, -1)
        else
            starts  =   copy(points[1:end-1])
            endings =   copy(points[2:end])
        end

        paths   =   GetBZPath.(Ref(bz), starts, endings ; nearest = nearest , exclusive = true)
        paths   =   reduce(vcat, paths)

        return paths
    end


    ##### Momentum pahse factors needed when doing FFT
    function MomentumPhaseFFT(bz::BZ, uc::UnitCell)

        momentumShift  =   (2 * pi) .* (bz.shift .+ Ref(1) .- ((bz.gridSize .+ 2 .- (bz.gridSize .% 2)) / 2) + (angle.(uc.BC) / (2 * pi))) ./ (bz.gridSize)
        grid           =    collect.(Meshgrid(bz.gridSize)) .- Ref(ones(Int64, length(bz.gridSize)))
        phaseShift     =    exp.(-im .* dot.(grid, Ref(momentumShift)))

        return  phaseShift
    end

end