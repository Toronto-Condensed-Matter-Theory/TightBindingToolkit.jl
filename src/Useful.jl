module Useful
    export GetAllOffsets, VecAngle, Meshgrid, BinarySearch, DistFunction, DeriDistFunction , GetIndexPath, FFTArrayofMatrix, Central_Diff, Arrayfy, DeArrayfy, GetPhase, SegmentIntersection, SwitchKroneckerBasis, SecantSearch, ImpIndices

    using LinearAlgebra, Statistics, FFTW, TensorCast

@doc """
```julia
GetAllOffsets(OffsetRange::Int64, dim::Int64) --> Vector{Vector{Int64}}
```
Given a range, returns the set of all possible Bond offsets such that each element of the offset vector lies in [-`OffsetRange`, `OffsetRange`].

"""
    function GetAllOffsets(OffsetRange::Int64, dim::Int64) :: Vector{Vector{Int64}}
        if dim==1
            offsets 	=	collect([i] for i in OffsetRange:-1:-OffsetRange)
        elseif  dim==2
            offsets 	=	reshape([[i, j] for i in OffsetRange:-1:-OffsetRange, j in OffsetRange:-1:-OffsetRange], (2*OffsetRange+1)^2)
        elseif  dim==3
            offsets 	=	reshape([[i, j, k] for i in OffsetRange:-1:-OffsetRange, j in OffsetRange:-1:-OffsetRange, k in OffsetRange:-1:-OffsetRange], (2*OffsetRange+1)^3)
        else
            println("Does not work for dimensions = ", dim)
        end
        return offsets
    end


@doc """
```julia
VecAngle(a::Vector{Float64}, b::Vector{Float64})
```
returns the angle b/w two given vectors.

"""
    function VecAngle(a::Vector{Float64}, b::Vector{Float64})
        return acos(clamp(a⋅b/(norm(a)*norm(b)), -1, 1))
    end


@doc """
```julia
Meshgrid(grid::Vector{Int64} ; starts::Vector{Int64} = zeros(Int64, length(grid)))
```
returns a meshgrid of (i_1, i_2, ..., i_n) such that i_j ∈ [starts[j], starts[j] + 1, ..., grid[j]] ∀ j in [1, 2, ..., n = length(grid)]

"""
    function Meshgrid(grid::Vector{Int64} ; starts::Vector{Int64} = ones(Int64, length(grid)))
        inds    =   UnitRange.( starts, grid)
        inds    =   collect(Base.product(inds...))
        return inds
    end


@doc """
```julia
BinarySearch(target::Float64, xRange::Tuple{Float64, Float64}, f::T, args::Tuple ; tol::Float64=0.001)
```
General function implementing Binary search on a monotonic function f(x, args...)=target, in the range x ∈ xRange, with tolerance tol. 

"""
    function BinarySearch(target::Float64, xRange::Tuple{Float64, Float64}, f::T, args::Tuple ; initial::Float64 = mean(xRange), x_tol::Float64 = 1e-3, target_tol::Float64 = 1e-6)::Float64 where T<:Function
        ##### ///TODO pass an initial guess 
        xExt = collect(xRange)
        ##### Choosing initial point to start Binary search with. 
        if initial < xExt[begin]
            current = xExt[begin]
        elseif xExt[end] < initial
            current = xExt[end]
        else
            current = initial
        end 

        ##### ///TODO : Fix bug when steps = 0
        if xExt[end]!= xExt[begin]
            steps       =   max(Int(ceil(log2((xExt[end]-xExt[begin])/(x_tol)))), 1)
        end
        
        for i in 1:steps
            check = f(current, args...)

            if check - target > target_tol
                xExt[end] = current
            elseif check - target < -target_tol
                xExt[begin] = current
            else
                break
            end

            current     =   mean(xExt)
        end
    
        return current
    end


    function SecantSearch(target::Float64, xInitials::Vector{Float64}, yInitials::Vector{Float64}, f::T, args::Tuple ; x_tol::Float64 = 1e-3, target_tol::Float64 = 1e-4, max_iter::Int64 = 10)::Float64 where T<:Function

        x0, x1  =   xInitials
        y0, y1  =   yInitials
        
        for i in 1:max_iter
            x  =   x1 - (x1 - x0) * (y1 - target) / (y1 - y0)
            check   =   f(x, args...)

            if abs(check - target) < target_tol || abs(x - x1) < x_tol
                break
            end

            x0, x1  =   x1, x
            y0, y1  =   y1, check
        end
    
        return x
    end


@doc """
```julia
DistFunction(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
```
Distribution function. `stat`=1 --> Bose-Einstein distribution, and `stat`=-1 --> Fermi distribution.

"""
    function DistFunction(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
        return @. 1 / (exp((E - mu) / T) - stat)
    end


@doc """
```julia
DeriDistFunction(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
```
derivative of [`DistFunction`](@ref) w.r.t the energy. `stat`=1 --> Bose-Einstein distribution, and `stat`=-1 --> Fermi distribution.

"""
    function DeriDistFunction(E ; T::Float64, mu::Float64=0.0, stat::Int64=-1)
        df  =   @. - (1/T) * exp((E - mu) / T) / ((exp((E - mu) / T) - stat) ^ 2)
        if isnan(df)
            df  =   0.0
        end

        return df 
    end


@doc """
```julia
GetIndexPath(start::Vector{Int64}, ending::Vector{Int64} ; exclusive::Bool=true) --> Vector{Vector{Int64}}
```
Returns a path in index-space of a discretized lattice which joins the two sets of indices `start` and `ending`. 
If the input `exclusive` is set to `true`, the returned path will NOT contain the `ending` point itself.
Works in any dimensions.

"""
    function GetIndexPath(start::Vector{Int64}, ending::Vector{Int64} ; exclusive::Bool=true) :: Vector{Vector{Int64}}
    
        dx  = ending - start
        i   = findfirst(!=(0), dx)    ##### find the first dimension where the two given points differ.
        ts  = (1/abs(dx[i])) .* collect(0:abs(dx[i]) - exclusive)    ##### choosing the paramtrization of a straight line in generalized dimensiosn s.t. the increment is ±1 along this dimension.
    
        coords = []
        for j in 1:length(start)
            xj = start[j] .+ (dx[j] .* ts)  ##### a straight line in generalized dimensions b/w two given points : (x[j] - start[j]) / (ending[j] - start[j]) = t ∈ [0, 1] 
            push!(coords, round.(Int64, xj))
        end
    
        coords = reduce(hcat, coords)   ##### combining all sets of indices into sets of points
        path = Vector{eltype(coords)}[eachrow(coords)...]
    
        return path
    end


    function Arrayfy(A::Array{Matrix{T}, N}) :: Array{T, 3} where {T<:Number, N}

        grid    =   size(A)

        G       =   reshape(A, prod(grid))       ##### Flattening all super-indices for casting 
        @cast G[k, i, j] |= G[k][i, j] 
        
        return G
    end


    function DeArrayfy(A::Array{T, 3}, grid::Vector{Int64}) :: Array{Matrix{T}, length(grid)} where {T<:Number}

        @cast G[k][i, j] |= A[k, i, j]
        G   =   reshape(G, grid...)
        
        return G
    end

    function FFTArrayofMatrix(Gs::Array{Matrix{ComplexF64}, T}) where {T}

        N       =   size(Gs[begin])      ##### size of the matrix being FFTed
        grid    =   size(Gs)        ##### size of the array of matrices

        G   =   reshape(Gs, prod(grid))       ##### Flattening all super-indices for casting 
        @cast G[k, i, j] |= G[k][i, j]        ##### Casting into an array for FFTW
        G   =   reshape(G, grid..., N...)     ##### Unflattenng all super-indices back to their original ranks

        Gr  =   fft(G, collect(1 : length(grid)))   ##### FFT on the momentum indices
        Gr  =   reshape(Gr, prod(grid), N...)
        @cast Gr[k][i, j] |= Gr[k, i, j]

        Gr    =   reshape(Gr, grid...) / prod(grid)

        return Gr
    end


@doc """
f : Can be an array of anything (even static arrays)
delta : vector of real numbers with length = spatial dimensions = rank of f ---> the displacement vector when calculating the gradient
PBC   : vector of boolean with length = spatial dimensions = rank of f ---> if the ith dimension has PBC or not.

"""
    function Central_Diff(f::Array{Any} ;  delta::Vector{Float64} = ones(Float64, length(size(f))), PBC::Vector{Bool} = repeat([true], length(size(f))))
        dims    =   length(size(f))
        IdMat   =   Matrix(1I, dims, dims)  ##### Identity matrix of size dims x dims where dims = rank of f = "spatial "dimensions. Eg in 3d f=f(x, y, z) = f[ix, iy, iz] => dims=3
        f_plus  =   circshift.(Ref(f), eachrow(IdMat))  
        """ circshift(f, v) shifts the array indices of f by a vector v. 
            eachrow(IdMat) basically gives [1, 0, ..., 0], [0, 1, 0, ..., 0] and so on.
            f_plus is a vector of arrays of same dimensions as f s.t. f_plus[i] is f shifted by one index to the right along the ith dimension""" 
        f_minus =   circshift.(Ref(f), eachrow(-IdMat))
        """ f_minus is a vector of arrays of same dimensions as f s.t. f_minus[i] is f shifted by one index to the left along the ith dimension""" 
        diff    =   (f_minus .- f_plus) ./ (2*delta)
        """ diff is another vector of arrays of same dimensions as f.
            diff[i][ix1,..., ixn] = d/dx^i(f)(ix1, ..., ixn)"""
        
    
        for i in 1:dims
            if !PBC[i]
                """
                selectdim(x, i, ind) = x[:, :, ..., ind, :, ..., :] where ind is in the ith dimension
                the relevant boundary for d/dx^i is ind=1, end for the ith dimension
                For such cases, instead of f_plus-f_minus, we need f_plus-f, and f-f_minus
                """
                left_edge    =   selectdim(diff[i], i, 1)
                left_edge   .=   (selectdim(f_minus[i], i, 1) .- selectdim(f, i, 1)) ./ delta[i]   ##### Pointer black magic???
                right_edge   =   selectdim(diff[i], i, size(f)[i])
                right_edge  .=   (selectdim(f, i, size(f)[i]) .- selectdim(f_plus[i], i, size(f)[i])) ./ delta[i] 
            end
        end
    
        return diff
    end

    function Central_Diff(f::Array{T, N} ;  delta::Vector{Float64} = ones(Float64, length(size(f))), PBC::Vector{Bool} = repeat([true], length(size(f)))) where {T<:Union{Number, Vector{<:Number}, Matrix{<:Number}}, N}
        dims    =   length(size(f))
        IdMat   =   Matrix(1I, dims, dims)  ##### Identity matrix of size dims x dims where dims = rank of f = "spatial "dimensions. Eg in 3d f=f(x, y, z) = f[ix, iy, iz] => dims=3
        f_plus  =   circshift.(Ref(f), eachrow(IdMat))  
        """ circshift(f, v) shifts the array indices of f by a vector v. 
            eachrow(IdMat) basically gives [1, 0, ..., 0], [0, 1, 0, ..., 0] and so on.
            f_plus is a vector of arrays of same dimensions as f s.t. f_plus[i] is f shifted by one index to the right along the ith dimension""" 
        f_minus =   circshift.(Ref(f), eachrow(-IdMat))
        """ f_minus is a vector of arrays of same dimensions as f s.t. f_minus[i] is f shifted by one index to the left along the ith dimension""" 
        diff    =   (f_minus .- f_plus) ./ (2*delta)
        """ diff is another vector of arrays of same dimensions as f.
            diff[i][ix1,..., ixn] = d/dx^i(f)(ix1, ..., ixn)"""
        
    
        for i in 1:dims
            if !PBC[i]
                """
                selectdim(x, i, ind) = x[:, :, ..., ind, :, ..., :] where ind is in the ith dimension
                the relevant boundary for d/dx^i is ind=1, end for the ith dimension
                For such cases, instead of f_plus-f_minus, we need f_plus-f, and f-f_minus
                """
                left_edge    =   selectdim(diff[i], i, 1)
                left_edge   .=   (selectdim(f_minus[i], i, 1) .- selectdim(f, i, 1)) ./ delta[i]   ##### Pointer black magic???
                right_edge   =   selectdim(diff[i], i, size(f)[i])
                right_edge  .=   (selectdim(f, i, size(f)[i]) .- selectdim(f_plus[i], i, size(f)[i])) ./ delta[i] 
            end
        end
    
        return diff
    end


    function GetPhase(v::Vector{ComplexF64}) :: Vector{Float64}

        return mod.(angle.(v) / (2 * pi), Ref(1.0))
    end


@doc """
```julia
SegmentIntersection(SegmentPoints::Vector{Vector{Float64}}, LinePoints::Vector{Vector{Float64}})-->Float64
```
Given two points on a line and two points on a segment, returns the parameter t such that the intersection point is given by (1-t) * SegmentPoints[1] + t * SegmentPoints[2]. 
If the intersection point is not on the segment, returns a value outside the range [0, 1].
"""
    function SegmentIntersection(SegmentPoints::Vector{Vector{Float64}}, LinePoints::Vector{Vector{Float64}})::Float64

        M1, M2  =   SegmentPoints
        r1, r2  =   LinePoints

        dx  =   (r2[1] - r1[1])
        dy  =   (r2[2] - r1[2])
        dMx =   (M2[1] - M1[1])
        dMy =   (M2[2] - M1[2])

        t       =   (dx * (r1[2] - M1[2]) - dy * (r1[1] - M1[1])) / (dx * dMy - dy * dMx)

        return t
    end

    function Linearize(Ns::Tuple{Int64, Int64}, index::Tuple{Int64, Int64})::Int64

        return Ns[2] * (index[1] - 1) + index[2]
    end

    function SwitchKroneckerBasis(Ns::Tuple{Int64, Int64})::Vector{Int64}

        indices     =   Meshgrid(collect(Ns))
        indices     =   vcat(indices...)
        order       =   Linearize.(Ref(Ns), indices)

        return order
    end


    function ImpIndices(A::Matrix{Float64}, sizeB::Tuple{Int64, Int64})::Vector{Tuple{Int64, Int64}}
        inds = Tuple{Int64, Int64}[]

        for (i, j) in Meshgrid(collect(sizeB))

            B = zeros(Float64, sizeB...)
            B[i, j] = 1.0

            C = B * A

            if sum(.!(isapprox.(C, zeros(Float64, size(C))))) > 0   ##### if the matrix product is not zero => this element of B can lead to non-zero result
                push!(inds, (i, j))
            end
        end

        return inds
    end



end