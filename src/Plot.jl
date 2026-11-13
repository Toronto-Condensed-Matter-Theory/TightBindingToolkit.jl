module PlotTB
    export Plot_UnitCell, Plot_Band_Contour, Plot_Band_Structure, Plot_FS, Plot_FS!, FirstBZBoundary, Plot_Fields, Plot_Lattice

    using LinearAlgebra, LaTeXStrings, Plots

    using ..TightBindingToolkit.Useful: GetAllOffsets
    using ..TightBindingToolkit.SpinMatrices: SpinMats
    using ..TightBindingToolkit.UCell:UnitCell
    using ..TightBindingToolkit.DesignUCell:Lookup
    using ..TightBindingToolkit.LatticeStruct: Lattice
    using ..TightBindingToolkit.BZone:BZ, GetRLVs, GetQIndex, GetBZPath, CombinedBZPath, ReduceQ
    using ..TightBindingToolkit.Hams:Hamiltonian, FillHoppingHamiltonian, BondDisplacements, OnSiteBlocks
    using ..TightBindingToolkit.TBModel:Model
    using ..TightBindingToolkit.BdG:BdGModel


@doc """
```julia
Plot_UnitCell(uc::UnitCell ; range::Int64 = 1, bond_cmp::Symbol = :matter, sub_cmp::Symbol = :rainbow, plot_conjugate::Bool=false, plot_labels::Vector{String} = unique(getproperty.(uc.bonds, :label)), plot_arrows::Bool = true, bond_opacity::Float64 = 0.6, site_size::Float64 = 16.0, bond_thickness::Tuple{Float64, Float64} = (4.0, 2.0), bond_rev::Bool=false, plot_lattice::Bool=false)
```
Function to plot the `UnitCell`. 
- `range` determines the range of UnitCells plotted in real-space. default is ±1.
- `bond_cmp` determines the colorScheme chosen to differentiate b.w different hopping types.
- `plot_conjugate` switches on whether the conjugate of the bonds in `UnitCell` are also plotted or not.
- `sub_cmp` is the choden colorScheme for the different sublattices.
- `plot_labels` is the vector of bond labels to be plotted.
- `plot_arrows` : whether to plot arrows on bonds.
- `bond_opacity` : alpha value of bonds being plotted.
- `bond_thickness` : linewidth of bonds being plotted.
- `site_size`: size of sublattice points.
- `plot_lattice`: plot bonds on all sites or only on one unit cell.

"""
    function Plot_UnitCell(uc::UnitCell ; range::Int64 = 1, bond_cmp::Symbol = :matter, sub_cmp::Symbol = :rainbow, plot_conjugate::Bool=false, plot_labels::Vector{String} = unique(getproperty.(uc.bonds, :label)), plot_arrows::Bool = true, bond_opacity::Float64 = 0.6, site_size::Float64 = 16.0, bond_thickness::Tuple{Float64, Float64} = (4.0, 2.0), bond_rev::Bool=false, plot_lattice::Bool=false)
        
        dim         =   length(uc.primitives)
        @assert dim == 2 "Unit Cell plotting only works for 2d right now!"
        offsets     =   GetAllOffsets(range, dim)

        subCmp     =   cgrad(sub_cmp, max(2 , length(uc.basis)), categorical=true, rev=true)

        p       =   plot(aspect_ratio=:equal, grid=false)
        ##### Plotting sites
        for offset in offsets
            shift   =   sum(offset .* uc.primitives)

            for (b, basis) in enumerate(uc.basis)

                site   =   basis + shift
                ##### Lattice sites
                scatter!(Tuple(site), label = "", markercolor=:darkorange1, markersize=site_size, markeralpha = 0.7, markerstrokecolor = subCmp[b], markerstrokewidth = site_size / 2)
                annotate!(site..., Plots.text(L"\mathbf{%$(b)}", :hcenter, :vcenter, round(Int64, site_size), :beige)) ##### Sublattice labels


                if offset == zeros(Int64, dim) && b==1
                    ##### Lattice primitive vectors
                    for (p, primitive) in enumerate(uc.primitives)
                        base    =   uc.basis[begin]
                        target  =   uc.basis[begin] .+ primitive

                        if p==1
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=true , color=:seagreen , linewidth=8 , label = "primitives", linealpha=0.75)
                        else
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=true , color=:seagreen , linewidth=8 , label = "", linealpha=0.75)
                        end
                    end

                end
            end

        end

        ##### Different labels, colors, and thickness of each bond types.
        labels      =   unique(getproperty.(uc.bonds, :label))
        cmp         =   cgrad(bond_cmp, max(2 , length(labels)), categorical=true, rev=bond_rev)
        thickness   =   collect(bond_thickness[1]:-(bond_thickness[1] - bond_thickness[2])/(max(2 , length(labels)) -1):bond_thickness[2])   
        counts      =   zeros(Int64, length(labels))
        
        ##### Plotting bonds
        if !plot_lattice

            for bond in uc.bonds
                base    =   uc.basis[bond.base]
                target  =   uc.basis[bond.target] .+ sum(bond.offset .* uc.primitives)
                base_conj       =   uc.basis[bond.target]
                target_conj     =   uc.basis[bond.base] .- sum(bond.offset .* uc.primitives)

                index   =   findfirst(==(bond.label), labels)
                counts[index]   +=  1

                if counts[index] == 1 && bond.label in plot_labels
                    plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "$(bond.label)", linealpha=bond_opacity)
                elseif counts[index] > 1 && bond.label in plot_labels
                    plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=bond_opacity)
                end
                ##### Plotting the conjugate of each bond
                if plot_conjugate && bond.label in plot_labels
                    plot!([base_conj[begin] , target_conj[begin]] , [base_conj[end], target_conj[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=0.75 * bond_opacity)
                end
            end

        else
            for offset in offsets
                counts      =   zeros(Int64, length(labels))
                shift       =   sum(offset .* uc.primitives)

                for bond in uc.bonds
                    base    =   shift + uc.basis[bond.base]
                    target  =   shift + uc.basis[bond.target] .+ sum(bond.offset .* uc.primitives)

                    base_conj       =   shift + uc.basis[bond.target]
                    target_conj     =   shift + uc.basis[bond.base] .- sum(bond.offset .* uc.primitives)
    
                    index   =   findfirst(==(bond.label), labels)
                    counts[index]   +=  1
                    
                    if counts[index] == 1 && bond.label in plot_labels 
                        if offset == zeros(Int64, dim)
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "$(bond.label)", linealpha=bond_opacity)
                        else
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=bond_opacity)
                        end
                    elseif counts[index] > 1 && bond.label in plot_labels 
                        plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=bond_opacity)
                    end

                    if plot_conjugate && bond.label in plot_labels
                        plot!([base_conj[begin] , target_conj[begin]] , [base_conj[end], target_conj[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=0.75 * bond_opacity)
                    end

                end
            end
        end

        xlabel!(L"x")
        ylabel!(L"y")
        title!("Unit Cell")

        return p

    end


    function Plot_Lattice(lat::Lattice{T} ; bond_cmp::Symbol = :matter, sub_cmp::Symbol = :rainbow, plot_labels::Vector{String} = unique(getproperty.(lat.uc.bonds, :label)), plot_arrows::Bool = true, bond_opacity::Float64 = 0.6, site_size::Float64 = 16.0, bond_thickness::Tuple{Float64, Float64} = (4.0, 2.0), bond_rev::Bool=false) where {T}

        dim         =   length(lat.uc.primitives)
        @assert dim == 2 "Unit Cell plotting only works for 2d right now!"

        subCmp     =   cgrad(sub_cmp, max(2 , length(lat.uc.basis)), categorical=true, rev=true)

        p       =   plot(aspect_ratio=:equal, grid=false)

        ##### Plotting sites
        for (site, info) in lat.positions

            position    =   info[1]
            sub, offset =   info[2]

            ##### Lattice sites
            scatter!(Tuple(position), label = "", markercolor=:darkorange1, markersize=site_size, markeralpha = 0.7, markerstrokecolor = subCmp[sub], markerstrokewidth = 0.4 * site_size )
            annotate!(position..., Plots.text(L"\mathbf{%$(site)}", :hcenter, :vcenter, round(Int64, site_size), :beige)) 

        end

        ##### Different labels, colors, and thickness of each bond types.
        labels      =   unique(getproperty.(lat.uc.bonds, :label))
        cmp         =   cgrad(bond_cmp, max(2 , length(labels)), categorical=true, rev=bond_rev)
        thickness   =   collect(bond_thickness[1]:-(bond_thickness[1] - bond_thickness[2])/(max(2 , length(labels)) -1):bond_thickness[2])   
        counts      =   zeros(Int64, length(labels))

        for site in 1:lat.length
            for (coord, neighbour) in enumerate(lat.BondSites[site, :])

                base    =   lat.positions[site][1]
                if neighbour != 0
                    target  =   lat.positions[neighbour][1]
                    label   =   lat.BondLabels[site, coord]
                    dist    =   lat.BondDists[site, coord]

                    if norm(target - base) ≈ dist

                        index   =   findfirst(==(label), labels)
                        counts[index]   +=  1

                        if counts[index] == 1 && label in plot_labels
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "$(label)", linealpha=bond_opacity)
                        elseif counts[index] > 1 && label in plot_labels
                            plot!([base[begin] , target[begin]] , [base[end], target[end]] , arrow=plot_arrows , color=cmp[index] , linewidth=thickness[index] , label = "", linealpha=bond_opacity)
                        end
                    end
                end

            end
        end

        xlabel!(L"x")
        ylabel!(L"y")
        title!("Lattice")

        return p

    end

@doc """
```julia
Plot_Fields!(uc::UnitCell ; OnSiteMatrices::Vector{Matrix{ComplexF64}}=SpinMats((uc.localDim - 1)//2), scale::Float64 = 1.0, range::Int64 = 1, cmp::Symbol = :thermal, field_thickness::Float64 = 1.0, field_opacity::Float64 = 0.6, use_lookup::Bool = false, site_size::Float64 = 12.0)
```
Function to plot the Fields of the UnitCell.
``
- `range` determines the range of UnitCells plotted in real-space. default is ±1.
- `cmp` determines the colorScheme of chosen to differentiate b.w different fields.
- `field_thickness` linewidth of the fields being plotted.
- `field_opacity` : alpha value of the fields being plotted.
- `use_lookup` : if you have implemented on-site terms using Bonds and not fields. 
- `OnSiteMatrices` : the matrices w.r.t which the lookup table on each site will be decomposed to get a vector field on each site.
- `site_size` : size of sublattice points being plotted.


"""
    function Plot_Fields(uc::UnitCell ; OnSiteMatrices::Vector{Matrix{ComplexF64}}=SpinMats((uc.localDim - 1)//2), scale::Float64 = 1.0, range::Int64 = 1, cmp::Symbol = :thermal, field_thickness::Float64 = 1.0, field_opacity::Float64 = 0.6, use_lookup::Bool = false, site_size::Float64 = 12.0)
        
        dim         =   length(uc.primitives)
        @assert dim == 2 "Unit Cell plotting only works for 2d right now!"
        offsets     =   GetAllOffsets(range, dim)

        fields  =   Vector{Float64}[]
        if use_lookup
            lookup  =   Lookup(uc)

            for (b, basis) in enumerate(uc.basis)
                fieldMat   =   get(lookup, (b, b, zeros(Int64, 2)), zeros(Float64, 3))
                field      =   real.(tr.( adjoint.(OnSiteMatrices) .* Ref(fieldMat)) ./ (tr.(adjoint.(OnSiteMatrices) .* OnSiteMatrices)))[1:3]
                push!(fields, scale .* field)
            end

        else
            for (b, basis) in enumerate(uc.basis)
                field   =   uc.fields[b][1:3]
                push!(fields, scale .* field)
            end
        end

        fieldZs     =   (getindex.(fields, 3) ./ norm.(fields))
        cmp         =   cgrad(cmp)

        p       =   plot(aspect_ratio=:equal, grid=false, colorbar = true, bg_legend = :transparent)
        ##### Plotting sites
        for offset in offsets
            shift   =   sum(offset .* uc.primitives)

            for (b, basis) in enumerate(uc.basis)

                site   =   basis + shift
                field  =   fields[b][1:2]
                ##### Lattice sites
                if offset == zeros(Int64, dim)
                    scatter!(Tuple(site), label = "", markercolor=:darkorange1, markersize=site_size, markerstrokecolor = :darkred, markerstrokewidth = 4) ##### Differentiationg the unit cell
                    theta   =   round(acos(fieldZs[b]) / pi, digits = 3)
                    plot!([site[begin] - (field[begin] / 2) , site[begin] + (field[begin] / 2)] , [site[end] - (field[end] / 2), site[end] + (field[end] / 2)] , arrow=true , color=cmp[(fieldZs[b] + 1)/2] , linewidth=field_thickness , label = L"\theta = %$(theta) \pi", linealpha=field_opacity)

                else
                    scatter!(Tuple(site), label = "", markercolor=:darkorange1, markersize=site_size)
                    plot!([site[begin] - (field[begin] / 2) , site[begin] + (field[begin] / 2)] , [site[end] - (field[end] / 2), site[end] + (field[end] / 2)] , arrow=true , color=cmp[(fieldZs[b] + 1)/2] , linewidth=field_thickness , label = "" , linealpha=field_opacity)
                end                

                if offset == zeros(Int64, dim)
                    annotate!(site..., Plots.text(L"\mathbf{%$(b)}", :hcenter, :vcenter, Int(site_size), :beige)) ##### Sublattice labels
                end
            end

        end

        xlabel!(L"x")
        ylabel!(L"y")
        title!("Fields")

        return p

    end


@doc """
```julia
plot_band_contour!(Ham::Hamiltonian , bz::BZ , band_index::Int64) --> Plots.plot()
```
Function to draw equal energy contours of the bands in `Hamiltonian`, specifically for the band with the given `band_index`.

"""
    function Plot_Band_Contour(Ham::Hamiltonian , bz::BZ , band_index::Int64 ; cmp::Symbol = :turbo)
        gr()
        @assert band_index <= length(Ham.bands[begin]) "Given band does not exist in the given Hamiltonian!"
        @assert length(size(Ham.bands)) == 2 "Contour plots only work for 2d Hamiltonians"

        xgrid = getindex.(bz.ks , 1)
        ygrid = getindex.(bz.ks , 2)

        plt   = contourf(xgrid , ygrid , getindex.(Ham.bands , band_index), levels=20, color=cmp, size=(600,600))
        xlabel!(L"k_x")
        ylabel!(L"k_y")
        title!(L"Band :  %$band_index")
        return plt
    end 


@doc """
```julia
plot_band_structure(M<:Union{Model, BdGModel}, path::Vector{Vector{Float64}},  band_index::Vector{Int64} = collect(1:length(M.Ham.bands[begin])) ; labels::Vector{} = repeat([""], length(path)), closed::Bool=true, nearest::Bool=true) --> Plots.plot()
```
Function to plot band structures of the `Model` along a `path` in the BZ determined by the given critical points.
Can take in multiple bands into account ∈ `band_index`.
`labels` are the Plot labels of the critical points.

"""
    function Plot_Band_Structure(M::T, path::Vector{Vector{Float64}}, band_index::Vector{Int64} = collect(1:length(M.Ham.bands[begin])) ; labels::Vector{} = repeat([""], length(path)), closed::Bool=true, nearest::Bool=true, plot_legend::Bool = true, framestyle::Symbol = :box, guidefontsize::Int64 = 14, tickfontsize::Int64 = 12, font::String = "Helvetica", plot_title::Bool=true) where {T<:Union{Model, BdGModel}}
        
        ##### the k-points taken along the path joining the given points
        bzpath     = CombinedBZPath(M.bz, path ; nearest = nearest, closed = closed)
        path_index = GetQIndex.(bzpath, Ref(M.bz) ; nearest = nearest)
        ##### only plotting the given bands
        bands_from_index = getindex.(Ref(M.Ham.bands), CartesianIndex.(Tuple.(path_index)))

        label_indices    = getindex.(findmin.([norm.(Ref(ReduceQ(x,M.bz)).-bzpath) for x in path]) , 2)

        plt = plot(grid=false, legend = plot_legend, bg_legend = :transparent, framestyle = framestyle, guidefontsize = guidefontsize, tickfontsize = tickfontsize)
        for j in band_index
            plot!(getindex.(bands_from_index , j), labels=L"Band : %$j", lw = 2.0)
        end
        if typeof(M)==BdGModel
            hline!([0.0], linestyle = :dash, label="", linecolor=:black) ##### Effective chemical potential in BdG
        else
            hline!([M.mu], linestyle = :dash, label=L"\mu", lw=0.5, linecolor=:black) 
        end

        
        if !closed 
            xticks!(label_indices, labels)
            vline!(label_indices, linestyle=:dash, linecolor=:indigo, label="", lw=0.5)
        else
            xticks!(vcat(label_indices , [length(bzpath)]), vcat(labels , [labels[begin]]))
            vline!(label_indices, linestyle=:dash, linecolor=:indigo, label="", lw=0.5)
        end
        ylabel!("Energy", guidefontsize = guidefontsize, guidefont=font)
        if plot_title
            xlabel!("Path", guidefontsize = guidefontsize, guidefont=font)
            title!("Band Structure along path", titlefontsize = 12, guidefont=font)
        end

        return plt
    end


@doc """
```julia
plot_FS(Ham::Hamiltonian , bz::BZ , Efermi::Vector{Float64} , band_index::Vector{Int64})--> Plots.plot()
```
Function to draw the fermi surface at `Efermi` for the given `Hamiltonian` on the given `BZ`. 

- `cmp` determines the colorMap used for the contours.
- `cbar` determines whether to plot the colorbar or not.    

"""
    function Plot_FS(Ham::Hamiltonian , bz::BZ , Efermi::Vector{Float64} , band_index::Vector{Int64} ; cmp::Symbol = :turbo, cbar::Bool=false)
        gr()
        @assert length(size(Ham.bands)) == 2 "Fermi surface plots only work for 2d Hamiltonians"
        offsets     =   GetAllOffsets(1, 2)

        plt = Plots.plot(aspect_ratio=:equal)

        max_range   =   maximum(norm.(bz.basis))

        for j in band_index
            energies    =   getindex.(Ham.bands , j)

            for offset in offsets
                shift   =   sum(offset .* bz.basis)
                ks      =   bz.ks .+ Ref(shift) 

                x       =   getindex.(ks, 1)
                y       =   getindex.(ks, 2)
                ##### a single level has no colour range to spread over
                scale   =   length(Efermi) > 1 ? (clims = extrema(Efermi), colorbar_ticks = round.(Efermi[1:max(1, length(Efermi) ÷ 10):end], digits = 3)) : NamedTuple()
                contour!(x, y, energies; levels=Efermi, color=cmp, size=(600,600), cbar = cbar, lw=2, fill = true, scale...)

                for key in keys(bz.HighSymPoints)
                    p = Tuple(shift .+ bz.HighSymPoints[key])

                    if prod(-max_range .< p .< max_range) && norm(shift)<1e-6

                        scatter!(p, label = "", markercolor=:yellow, markersize = 15.0)
                        annotate!(p..., Plots.text(L"%$key", :hcenter, :vcenter, 8, :red4))
                    end
                end

                if offset == [0, 0]
                    boundary_1  =   vcat(x[begin, :], x[:, end], x[end, :], x[:, begin])
                    boundary_2  =   vcat(y[begin, :], y[:, end], y[end, :], y[:, begin])
                    scatter!(boundary_1, boundary_2, label = "", alpha=0.5, markercolor="orange", markersize=0.5)
                end

            end
        end

        for (i, b) in enumerate(bz.basis)
            Plots.plot!([0 , getindex(b , 1)] , [0 , getindex(b , 2)] , arrow=true , color=:black , lw=2 , label = L"b_{%$i}", linestyle=:dash, linealpha=0.75)
        end

        xlabel!(L"k_x", guidefontsize = 9)
        ylabel!(L"k_y", guidefontsize = 9)
        Ef  =   round.(Efermi, digits=2)
        title!("Fermi Surface for bands : $(band_index) at " * L"E_f \in %$(extrema(Efermi))", titlefontsize = 12)

        max_range   =   maximum(norm.(bz.basis))
        xlims!((-max_range , max_range))
        ylims!((-max_range , max_range))

        return plt
    end


    ##### Band energies of the hopping Hamiltonian at an arbitrary momentum, with the k-independent parts computed once.
    function BandFunction(uc::UnitCell{2})
        deltas  =   BondDisplacements(uc)
        onsite  =   OnSiteBlocks(uc)
        return k -> eigvals(Hermitian(FillHoppingHamiltonian(uc, Float64.(collect(k)), deltas, onsite)))
    end

    ##### Distance from Γ to the first-BZ (Wigner-Seitz) boundary along the unit vector u.
    function ZoneEdge(u::Vector{Float64}, Gs::Vector{Vector{Float64}}) :: Float64
        return minimum(dot(G, G) / (2 * dot(u, G)) for G in Gs if dot(u, G) > 1e-9)
    end

    function ZoneVectors(uc::UnitCell{2}) :: Vector{Vector{Float64}}
        b1, b2  =   GetRLVs(uc)
        return [i * b1 + j * b2 for i in -2:2 for j in -2:2 if (i, j) != (0, 0)]
    end

    InZone(k::Vector{Float64}, Gs::Vector{Vector{Float64}} ; tol::Float64 = 1e-9) = all(dot(k, G) <= dot(G, G) / 2 + tol for G in Gs)


@doc """
```julia
FirstBZBoundary(uc::UnitCell{2} ; n::Int64 = 721) --> Vector{Vector{Float64}}
```
Points on the boundary of the first Brillouin zone (the Wigner-Seitz cell of the reciprocal lattice) of `uc`, going once around anticlockwise, closed.

"""
    function FirstBZBoundary(uc::UnitCell{2} ; n::Int64 = 721) :: Vector{Vector{Float64}}
        Gs  =   ZoneVectors(uc)
        return [ZoneEdge([cos(θ), sin(θ)], Gs) .* [cos(θ), sin(θ)] for θ in range(0, 2π, length = n)]
    end


@doc """
```julia
Plot_FS(uc::UnitCell{2}, mu::Real ; kwargs...) --> Plots.plot()
Plot_FS!(plt, uc::UnitCell{2}, mu::Real ; kwargs...) --> Plots.plot()
```
Draws the Fermi surface ``ε_n(k) = μ`` of the hopping Hamiltonian of `uc` in the first Brillouin zone, with the zone boundary.
The bands are evaluated directly at every point of a Cartesian grid, so no `BZ` or `Hamiltonian` is needed and any 2d lattice works.

- `Q` : a nesting vector (Cartesian). Also draws the Fermi surface translated by `Q` (folded back into the zone), dashed, and an arrow `k₀ → k₀ + Q` from a point of the Fermi surface, chosen so that `k₀ + Q` lies as close to the Fermi surface as possible. Where the two surfaces overlap, `Q` nests them.
- `bands` : the band indices to draw (default: all). Bands degenerate with the previous one (e.g. spin) are drawn once.
- `color`, `lw`, `ls`, `label` : style of the Fermi surface and its legend entry.
- `nest_color` : colour of the translated surface and of the arrow.
- `zone` : whether to draw the first-BZ boundary. `clip` : whether to hide the contours outside it.
- `npix` : grid points per side.

"""
    function Plot_FS(uc::UnitCell{2}, mu::Real ; kwargs...)
        plt     =   Plots.plot(aspect_ratio = 1, grid = false, framestyle = :box)
        return Plot_FS!(plt, uc, mu ; kwargs...)
    end

    function Plot_FS!(plt, uc::UnitCell{2}, mu::Real ; Q::Union{Nothing, AbstractVector{<:Real}} = nothing, bands = nothing,
                      color = :royalblue, lw::Real = 2, ls::Symbol = :solid, label = "",
                      nest_color = :crimson, zone::Bool = true, clip::Bool = true, npix::Int64 = 241)

        @assert length(uc.primitives) == 2 && length(uc.primitives[1]) == 2 "Fermi surface plots only work for 2d unit cells"
        energy      =   BandFunction(uc)
        Gs          =   ZoneVectors(uc)
        boundary    =   FirstBZBoundary(uc)
        box         =   1.08 * maximum(norm, boundary)
        ks          =   collect(range(-box, box, length = npix))
        dk          =   step(range(-box, box, length = npix))
        inside      =   [InZone([kx, ky], Gs ; tol = dk) for ky in ks, kx in ks]

        ##### one matrix per band, rows = ky (the convention of contour)
        band_index  =   something(bands, 1:length(energy(zeros(2))))
        function band_grids(shift)
            grid    =   [energy([kx, ky] .- shift) for ky in ks, kx in ks]
            per     =   [getindex.(grid, b) for b in band_index]
            keep    =   [i == 1 || maximum(abs.(per[i] .- per[i - 1])) > 1e-8 for i in eachindex(per)]
            return per[keep]
        end
        masked(E)   =   clip ? ifelse.(inside, E, NaN) : E

        ##### estimated distance from k to the Fermi surface, |ε - μ| / |∇ε|, over the drawn bands
        function distance(k)
            h   =   dk / 2
            e0, ex1, ex2, ey1, ey2  =   (energy(k .+ s)[band_index] for s in ([0, 0], [h, 0], [-h, 0], [0, h], [0, -h]))
            return minimum(@. abs(e0 - mu) / max(hypot(ex1 - ex2, ey1 - ey2) / (2h), 1e-12))
        end

        Es  =   band_grids(zeros(2))
        for E in Es
            contour!(plt, ks, ks, masked(E) ; levels = [Float64(mu)], color = color, lw = lw, ls = ls, colorbar = false, label = "")
        end
        label == "" || Plots.plot!(plt, [NaN], [NaN] ; color = color, lw = lw, ls = ls, label = label)

        if !isnothing(Q)
            Q   =   Float64.(collect(Q))
            for E in band_grids(Q)
                contour!(plt, ks, ks, masked(E) ; levels = [Float64(mu)], color = nest_color, lw = lw, ls = :dash, colorbar = false, label = "")
            end
            Plots.plot!(plt, [NaN], [NaN] ; color = nest_color, lw = lw, ls = :dash, label = L"FS $+\, Q$")

            ##### The arrow k0 → k0 + S, with S = Q + G for a small reciprocal lattice vector G (equivalent to Q), so that
            ##### both ends lie in the zone: nesting often joins opposite sides of the zone only up to a G.
            b1, b2  =   GetRLVs(uc)
            shifts  =   [(Q .+ i * b1 .+ j * b2, i == j == 0) for i in -1:1 for j in -1:1]
            ##### candidates: k0 within a grid step of the Fermi surface, and the distance of k0 + S from it
            candidates  =   Tuple{Vector{Float64}, Vector{Float64}, Bool, Float64}[]
            for ix in 2:npix-1, iy in 2:npix-1
                inside[iy, ix] || continue
                d   =   minimum(abs(E[iy, ix] - mu) / max(hypot(E[iy, ix+1] - E[iy, ix-1], E[iy+1, ix] - E[iy-1, ix]) / (2 * dk), 1e-12) for E in Es)
                d < dk || continue
                k   =   [ks[ix], ks[iy]]
                for (S, bare) in shifts
                    InZone(k .+ S, Gs) && push!(candidates, (k, S, bare, distance(k .+ S)))
                end
            end
            ##### k0 + S on the Fermi surface too if possible, taking the arrow whose midpoint is farthest from the surface
            ##### (so it does not run along it) and S = Q on a tie; else the closest miss. With no Fermi surface to start
            ##### from (e.g. Dirac points), Q is drawn from Γ.
            nested  =   filter(c -> c[4] < dk, candidates)
            k0, S, bare =   if !isempty(nested)
                                argmax(c -> (distance(c[1] .+ c[2] / 2), c[3]), nested)[1:3]
                            elseif !isempty(candidates)
                                argmin(c -> (c[4], !c[3]), candidates)[1:3]
                            else
                                zeros(2), Q, true
                            end
            ##### Q = 0 (up to a reciprocal lattice vector) has no arrow to draw
            if minimum(norm(first(s)) for s in shifts) > dk
                Plots.plot!(plt, [k0[1], k0[1] + S[1]], [k0[2], k0[2] + S[2]] ; arrow = true, color = nest_color, lw = lw,
                            label = bare ? L"Q" : L"Q + G")
            end
        end

        zone && Plots.plot!(plt, first.(boundary), last.(boundary) ; color = :black, lw = 1, label = "")
        Plots.plot!(plt ; aspect_ratio = 1, xlims = (-box, box), ylims = (-box, box), xlabel = L"k_x", ylabel = L"k_y")

        return plt
    end

end