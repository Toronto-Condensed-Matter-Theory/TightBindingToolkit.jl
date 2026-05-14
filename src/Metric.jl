module MetricTB

using LinearAlgebra
using LaTeXStrings
using Plots

using ..TightBindingToolkit.Hams: Hamiltonian, DiagonalizeHamiltonian!, GetVelocity!
using ..TightBindingToolkit.BZone: BZ
using ..TightBindingToolkit.Chern: FindLinks, FieldStrength, ChernNumber
using ..TightBindingToolkit.UCell: UnitCell
using ..TightBindingToolkit.Parameters: Param, CreateUnitCell!
using ..TightBindingToolkit.Useful: GetAllOffsets

export GeoTensor, Curvature, KuboChern, get_metric_measures, hexagon, plot_metric_data

@doc """
```julia
GeoTensor(Ham::Hamiltonian, subset::Vector{Int64}) --> Array{Matrix{ComplexF64}}
```
Computes the quantum geometric tensor for the given `subset` of bands in the `Hamiltonian`.
The velocity must have been computed first via `GetVelocity!`.

"""
function GeoTensor(Ham::Hamiltonian, subset::Vector{Int64})

    Vx = conj.(permutedims.(Ham.states)) .* Ham.velocity[1] .* Ham.states
    Vy = conj.(permutedims.(Ham.states)) .* Ham.velocity[2] .* Ham.states

    geotensor = similar(Ham.states, Matrix{ComplexF64})

    for k in eachindex(Ham.bands)
        Es = Ham.bands[k]
        vx = Vx[k]
        vy = Vy[k]

        geotensor[k] = zeros(ComplexF64, 2, 2)

        for i in subset
            for j in setdiff(1:length(Es), subset)
                denom = ((Es[j] - Es[i])^2)
                geotensor[k][1, 1] += (vx[i, j] * vx[j, i] ) / denom
                geotensor[k][1, 2] += (vx[i, j] * vy[j, i] ) / denom
                geotensor[k][2, 1] += (vy[i, j] * vx[j, i] ) / denom
                geotensor[k][2, 2] += (vy[i, j] * vy[j, i] ) / denom
            end
        end

    end

    return geotensor / length(Ham.bands)

end

@doc """
```julia
Curvature(Ham::Hamiltonian, subset::Vector{Int64}) --> Matrix{Float64}
```
Computes the Berry curvature for the given `subset` of bands using wavefunction link variables.

"""
function Curvature(Ham::Hamiltonian, subset::Vector{Int64})::Matrix{Float64}

    Links   =   FindLinks(Ham, subset)
    Field   =   FieldStrength(Links)
    curvature = angle.(Field)
    return curvature
end

@doc """
```julia
KuboChern(Ham::Hamiltonian, bz::BZ, mu::Float64) --> Float64
```
Computes the Chern number using the Kubo formula from velocity matrix elements.
All bands below `mu` are included in the calculation.

"""
function KuboChern(Ham::Hamiltonian, bz::BZ, mu::Float64)

    Vx = conj.(permutedims.(Ham.states)) .* Ham.velocity[1] .* Ham.states
    Vy = conj.(permutedims.(Ham.states)) .* Ham.velocity[2] .* Ham.states

    chern = 0.0 + im * 0.0
    for k in eachindex(Ham.bands)
        Es = Ham.bands[k]
        vx = Vx[k]
        vy = Vy[k]

        ind = searchsortedfirst(Es, mu)
        if ind == 1 || ind == length(Es)
            continue
        else
            for i in 1:ind-1
                for j in ind:length(Es)
                    chern += (vx[i, j] * vy[j, i] - vx[j, i] * vy[i, j]) / ((Es[j] - Es[i])^2)
                end
            end
        end

    end

    b1 = [bz.basis[1]; 0.0]
    b2 = [bz.basis[2]; 0.0]
    bzUnitArea = cross(b1, b2)[3] / (4 * pi^2)

    return imag(chern) * bzUnitArea * 2 * pi / length(Ham.bands)
end


@doc """
```julia
get_metric_measures(Js, UC, bz, HoppingParams, jhParam, band ; measures=true) --> Dict
```
Sweeps over values of `Js`, computing Chern numbers, quantum metric weight, and volume
for the specified `band` at each value.

"""
function get_metric_measures(Js::Vector{Float64}, UC::UnitCell, bz::BZ,
    HoppingParams::Vector{T}, jhParam::Param, band::Int64 ;
    measures::Bool = true) where {T}

    # Compute BZ unit area once
    b1 = [bz.basis[1]; 0.0]
    b2 = [bz.basis[2]; 0.0]
    bzUnitArea = abs(cross(b1, b2)[3])

    gap_belows    =   Float64[]
    gap_aboves    =   Float64[]
    mus     =   Float64[]
    Cherns_wfns  =   Float64[]
    Cherns_metric =   Float64[]
    weights =   Float64[]
    volumes =   Float64[]
    bandwidths = Float64[]

    for J in Js
        jhParam.value = [J]
        UC.bonds = []
        CreateUnitCell!(UC, HoppingParams)

        H = Hamiltonian(UC, bz)
        DiagonalizeHamiltonian!(H)

        energies = getindex.(H.bands, band)
        energy_below = getindex.(H.bands, band-1)
        energy_above = getindex.(H.bands, band+1)
        gap_below =  minimum((energies-energy_below))
        gap_above =  minimum((energy_above-energies))
        bandwidth = maximum(energies) - minimum(energies)

        push!(gap_belows, gap_below)
        push!(gap_aboves, gap_above)
        push!(bandwidths, bandwidth)

        if measures
            GetVelocity!(H, bz)
            if gap_above>1e-3 && gap_below>1e-3
                chern = ChernNumber(H, [band])
                push!(Cherns_wfns, chern)

                geo = GeoTensor(H, [band])
                metric = [real(mat) for mat in geo]
                berry = [-2*imag(mat) for mat in geo]
                metric_traces = [tr(mat) for mat in metric]
                berry_curvature = [mat[1, 2] for mat in berry]

                weight = sum(metric_traces) * bzUnitArea / (2*pi)
                push!(weights, weight)
                push!(Cherns_metric, sum(berry_curvature) * bzUnitArea / (2*pi))

                push!(volumes, sum(sqrt.(abs.(det.(metric))) * bzUnitArea / (pi)))
                println("W = $(round(weight, digits=3)), and C = $(round(chern, digits=3)) at Jh = $(J)")
            else
                push!(Cherns_wfns, 0.0)
                push!(weights, 0.0)
                push!(Cherns_metric, 0.0)
                push!(volumes, 0.0)
                println("W = 0.0, and C = 0.0 at Jh = $(J)")
            end
        end
    end

    data = Dict("Js" => Js, "Cherns_wfns" => Cherns_wfns, "Cherns_metric" => Cherns_metric,
        "weights" => weights, "volumes" => volumes, "bandwidths" => bandwidths, "gap_above"=>gap_aboves, "gap_below"=>gap_belows)
    return data
end

@doc """
```julia
hexagon(corner::Vector{Float64}) --> Vector{Tuple}
```
Returns the edges of a hexagon starting from `corner`, each edge as a tuple of (xs, ys) arrays.
Useful for drawing hexagonal Brillouin zone boundaries.

"""
function hexagon(corner::Vector{Float64})
    RotMat = [cos(pi/3) sin(pi/3); -sin(pi/3) cos(pi/3)]
    lines = []
    for i in 1:6
        p1 = corner
        p2 = RotMat * corner

        if abs(p2[1]-p1[1])>1e-6
            slope = (p2[2] - p1[2]) / (p2[1] - p1[1])
            xs = LinRange(p1[1], p2[1], 100)
            ys = slope .* (xs .- p1[1]) .+ p1[2]
        else
            ys = LinRange(p1[2], p2[2], 100)
            xs = p1[1] .* ones(100)
        end
        push!(lines, (xs, ys))
        corner = RotMat * corner
    end

    return lines
end

@doc """
```julia
plot_metric_data(data, kxs, kys, bz ; kwargs...) --> Plots.Plot
```
Plots a heatmap of `data` on a hexagonal Brillouin zone using Plots.jl.

Optional keyword arguments:
- `colorbar_title`: label for the colorbar (default: Berry curvature symbol)
- `to_save`: tuple `(Bool, String)` — if true, saves to the given filename
- `labels`: whether to show axis labels
- `clims`: color limits tuple
- `cmap`: colormap symbol
- `annotation`: annotation text

"""
function plot_metric_data(data::Matrix{T},
    kxs::Vector{Float64}, kys::Vector{Float64},
    bz::BZ;
    colorbar_title::AbstractString = L"\\Omega_n(\\mathbf{k})",
    to_save::Tuple{Bool, String} = (false, ""),
    labels::Bool=true,
    clims::Tuple{T, T} = extrema(data),
    cmap::Symbol=:blues,
    annotation::AbstractString = "") where {T}

    font = "Computer Modern"

    ssf_plot = plot(framestyle=:box, grid=false, aspect_ratio=:equal,
        guidefontstyle = font, tickfont = font, legendfont = font,
        guidefontsize = 12, tickfontsize = 12, legendfontsize = 10)

    tick_start = round(Int, clims[1])
    tick_end = round(Int, clims[2])

    heatmap!(kxs, kys, data', c=cmap,
        clim=clims, colorbar_title=colorbar_title,
        colorbar_ticks = collect(range(tick_start, step=2, stop=tick_end)))

    xlims!(-2.025, 2.025)
    ylims!(-2.025, 2.025)
    if labels
        xlabel!(L"k_x")
        ylabel!(L"k_y")
        xticks!([0.0], [""])
        yticks!([0.0], [""])
    else
        xticks!([0.0], [""])
        yticks!([0.0], [""])
    end

    if hasproperty(bz, :HighSymPoints) && haskey(bz.HighSymPoints, "K1")
        inner_hex = hexagon(bz.HighSymPoints["K1"])
        for edge in inner_hex
            xs, ys = edge
            plot!(xs, ys, label = "", lw=2.0, lc=:orange)
        end
    end

    annotate!([(1.8, 1.8, text(annotation, "Computer Modern", 12, :black))])

    if to_save[1]
        savefig(to_save[2])
    end

    return ssf_plot
end

end # module MetricTB
