#=
Regenerates the README figures of TightBindingToolkit.jl into docs/src/assets/readme/.

    julia --project=docs/readme -t 8 docs/readme/figures.jl

Every figure is computed with the package; Plots only draws the results.
=#
include("style.jl")
using TightBindingToolkit, LinearAlgebra

##### -------------------------------------------------------------------------------------------
##### Helpers
##### -------------------------------------------------------------------------------------------

"Reciprocal vectors bounding the first Brillouin zone of a 2d unit cell."
zone_vectors(uc) = (b = GetRLVs(uc); [i * b[1] + j * b[2] for i in -1:1 for j in -1:1 if (i, j) != (0, 0)])
in_zone(k, Gs) = all(dot(k, G) <= dot(G, G) / 2 + 1e-9 for G in Gs)

"Value of a quantity stored on the `bz` grid at the Cartesian momentum `k` (periodic lookup)."
function grid_lookup(values::AbstractMatrix, bz::BZ, k)
    B   =   hcat(bz.basis...)
    f0  =   B \ bz.ks[1, 1]
    f   =   B \ k .- f0
    i   =   mod(round(Int, f[1] * bz.gridSize[1]), bz.gridSize[1]) + 1
    j   =   mod(round(Int, f[2] * bz.gridSize[2]), bz.gridSize[2]) + 1
    return values[i, j]
end

"Cartesian grid covering the first BZ of `uc`, with `NaN` outside it."
function zone_map(f, uc ; n = 401)
    Gs      =   zone_vectors(uc)
    kmax    =   maximum(norm.(FirstBZBoundary(uc)))
    ks      =   range(-kmax, kmax, length = n)
    Z       =   [in_zone([kx, ky], Gs) ? f([kx, ky]) : NaN for ky in ks, kx in ks]
    return ks, Z
end

function draw_zone!(plt, uc ; kwargs...)
    edge    =   FirstBZBoundary(uc)
    plot!(plt, first.(edge), last.(edge) ; color = COLORS.ink, lw = 1.2, label = "", kwargs...)
end

bands_at(uc, k) = eigvals(Hermitian(FillHoppingHamiltonian(uc, k)))

"Bands of `uc` along straight segments joining `points`, and the path coordinate of each point."
function bands_along(uc, points ; n = 200)
    xs, Es, ticks = Float64[], Vector{Float64}[], [0.0]
    x = 0.0
    for (p, q) in zip(points[1:end-1], points[2:end])
        for t in range(0, 1, length = n)[1:end-1]
            k = p .+ t .* (q .- p)
            push!(xs, x + t * norm(q - p)); push!(Es, bands_at(uc, k))
        end
        x += norm(q - p); push!(ticks, x)
    end
    push!(xs, x); push!(Es, bands_at(uc, points[end]))
    return xs, reduce(hcat, Es)', ticks
end

##### -------------------------------------------------------------------------------------------
##### 1. Unit cell and bonds: honeycomb lattice with 1st, 2nd and 3rd neighbour hoppings
##### -------------------------------------------------------------------------------------------
function figure_unitcell()
    a1, a2  =   [1/2, sqrt(3)/2], [-1/2, sqrt(3)/2]
    uc      =   UnitCell([a1, a2], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddBasisSite!(uc, [0.0, 1/sqrt(3)])
    AddIsotropicBonds!(uc, 1/sqrt(3), -1.0, "t1")
    AddIsotropicBonds!(uc, 1.0,       -0.2, "t2")
    AddIsotropicBonds!(uc, 2/sqrt(3), -0.1, "t3" ; checkOffsetRange = 3)

    pos(b, o) = uc.basis[b] .+ o[1] .* a1 .+ o[2] .* a2
    plt = plot(aspect_ratio = :equal, axis = false, ticks = false, framestyle = :none,
               legend = :outerright, size = (560, 420), xlims = (-1.6, 1.6), ylims = (-1.25, 1.85))

    ##### the lattice: nearest-neighbour bonds of a patch of unit cells
    cells   =   [[i, j] for i in -3:3 for j in -3:3]
    for o in cells, bond in uc.bonds
        bond.label == "t1" || continue
        r1, r2 = pos(bond.base, o), pos(bond.target, o .+ bond.offset)
        plot!(plt, [r1[1], r2[1]], [r1[2], r2[2]] ; color = COLORS.muted, lw = 1, alpha = 0.6, label = "")
    end
    for o in cells, b in 1:2
        r = pos(b, o)
        scatter!(plt, [r[1]], [r[2]] ; color = (b == 1 ? COLORS.ink : :white), markerstrokecolor = COLORS.ink,
                 markersize = 5, label = "")
    end

    ##### the bonds leaving one site of sublattice A, coloured by neighbour shell
    for (label, color, name) in (("t1", COLORS.blue, L"t_1"), ("t2", COLORS.orange, L"t_2"), ("t3", COLORS.aqua, L"t_3"))
        first_bond = true
        for bond in uc.bonds
            (bond.label == label) || continue
            for (b, t, o) in ((bond.base, bond.target, bond.offset), (bond.target, bond.base, -bond.offset))
                b == 1 || continue
                r2 = pos(t, o) .- uc.basis[1]
                plot!(plt, [0, 0.93 * r2[1]], [0, 0.93 * r2[2]] ; color = color, lw = 2.5, arrow = arrow(:closed, :head, 0.6, 0.6),
                      label = first_bond ? name : "")
                first_bond = false
            end
        end
    end
    scatter!(plt, [0.0], [0.0] ; color = COLORS.red, markerstrokecolor = COLORS.ink, markersize = 8, label = "")
    savefig_readme(plt, "honeycomb_bonds")
end

##### -------------------------------------------------------------------------------------------
##### 2-3. Kagome lattice: band structure and density of states
##### -------------------------------------------------------------------------------------------
function kagome_unitcell()
    a1, a2  =   [1.0, 0.0], [1/2, sqrt(3)/2]
    uc      =   UnitCell([a1, a2], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddBasisSite!(uc, a1 / 2)
    AddBasisSite!(uc, a2 / 2)
    AddIsotropicBonds!(uc, 1/2, -1.0, "t1")
    return uc
end

const KAGOME_YLIMS = (-4.3, 2.4)

function figure_kagome_bands()
    uc      =   kagome_unitcell()
    b1, b2  =   GetRLVs(uc)
    Γ, K, M =   [0.0, 0.0], (2b1 + b2) / 3, b1 / 2
    xs, Es, ticks = bands_along(uc, [Γ, K, M, Γ])
    plt = plot(xs, Es ; color = COLORS.blue, label = "",
               xticks = (ticks, [L"\Gamma", L"K", L"M", L"\Gamma"]), xlims = extrema(xs), ylims = KAGOME_YLIMS,
               ylabel = L"E / |t|", minorticks = false)
    vline!(plt, ticks[2:end-1] ; color = COLORS.muted, lw = 0.8, label = "")
    savefig_readme(plt, "kagome_bands")
end

function figure_kagome_dos()
    uc      =   kagome_unitcell()
    bz      =   BZ([240, 240])
    FillBZ!(bz, uc)
    H       =   Hamiltonian(uc, bz)
    DiagonalizeHamiltonian!(H ; verbose = false)
    Ωs      =   collect(range(KAGOME_YLIMS..., length = 1201))
    ρ       =   DOS(Ωs, H ; spread = 0.03)
    plt = plot(Shape(vcat(0.0, ρ, 0.0), vcat(Ωs[1], Ωs, Ωs[end])) ; color = COLORS.blue, fillalpha = 0.15, linewidth = 0, label = "")
    plot!(plt, ρ, Ωs ; color = COLORS.blue, label = "",
          xlabel = L"\rho(E)", ylabel = L"E / |t|", ylims = KAGOME_YLIMS, xlims = (0, 0.4), minorticks = false)
    savefig_readme(plt, "kagome_dos")
end

##### -------------------------------------------------------------------------------------------
##### 4. Fermi surfaces of the triangular lattice at several fillings
##### -------------------------------------------------------------------------------------------
function figure_fermi_surfaces()
    uc      =   UnitCell([[1.0, 0.0], [1/2, sqrt(3)/2]], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddIsotropicBonds!(uc, 1.0, -1.0, "t1")

    bz      =   BZ([120, 120])
    FillBZ!(bz, uc)
    H       =   Hamiltonian(uc, bz)
    DiagonalizeHamiltonian!(H ; verbose = false)

    fillings=   [0.1, 0.3, 0.5, 0.65, 0.75]
    colors  =   [SEQUENTIAL[x] for x in range(0.35, 1.0, length = length(fillings))]
    ks, E   =   zone_map(k -> bands_at(uc, k)[1], uc ; n = 501)

    plt = plot(aspect_ratio = :equal, xlabel = L"k_x", ylabel = L"k_y", legend = :outerright,
               xlims = extrema(ks), ylims = extrema(ks), size = (560, 420))
    for (n, c) in zip(fillings, colors)
        model = Model(uc, bz, H ; filling = n, T = 1e-3)
        SolveModel!(model ; get_correlations = false, verbose = false)
        contour!(plt, ks, ks, E ; levels = [model.mu], color = c, lw = 2.2, colorbar = false)
        plot!(plt, [NaN], [NaN] ; color = c, label = L"n = %$n")   ##### legend entry for the contour
    end
    draw_zone!(plt, uc)
    savefig_readme(plt, "triangular_fermi_surfaces")
end

##### -------------------------------------------------------------------------------------------
##### 5-6. Haldane model: Chern-number phase diagram and Berry curvature
##### -------------------------------------------------------------------------------------------
##### Spinless Haldane model (localDim = 1): one orbital per site, so the lower band has C = 0 or ±1.
function haldane_unitcell(M, t2, φ ; t1 = -1.0)
    uc  =   UnitCell([[1/2, sqrt(3)/2], [-1/2, sqrt(3)/2]], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddBasisSite!(uc, [0.0, 1/sqrt(3)])
    AddAnisotropicBond!(uc, 1, 1, [0, 0],  M, 0.0, "M")
    AddAnisotropicBond!(uc, 2, 2, [0, 0], -M, 0.0, "M")
    AddIsotropicBonds!(uc, 1/sqrt(3), t1, "t1")
    for o in ([1, 0], [0, -1], [-1, 1])     ##### complex second-neighbour hopping, opposite chirality on A and B
        AddAnisotropicBond!(uc, 1, 1,  o, t2 * cis(φ), 1.0, "t2")
        AddAnisotropicBond!(uc, 2, 2, -o, t2 * cis(φ), 1.0, "t2")
    end
    return uc
end

function figure_haldane_phase_diagram()
    t2      =   0.2
    φs      =   range(-π, π, length = 181)
    Ms      =   range(-1.4, 1.4, length = 121)
    bz      =   BZ([36, 36])
    FillBZ!(bz, haldane_unitcell(0.0, t2, 0.0))
    C       =   zeros(length(Ms), length(φs))
    Threads.@threads for i in eachindex(φs)
        for (j, M) in enumerate(Ms)
            H = Hamiltonian(haldane_unitcell(M, t2, φs[i]), bz)
            DiagonalizeHamiltonian!(H ; verbose = false)
            C[j, i] = round(ChernNumber(H, [1]))
        end
    end
    phases  =   ((-1, COLORS.blue, "C = −1"), (0, colorant"#e4e2dc", "C = 0"), (1, COLORS.red, "C = +1"))
    plt = heatmap(φs ./ π, Ms ./ t2, C ; color = cgrad([p[2] for p in phases], 3, categorical = true), clims = (-1.5, 1.5),
                  colorbar = false, xlabel = L"\phi / \pi", ylabel = L"M / t_2", size = (560, 400), minorticks = false,
                  xlims = (-1, 1), ylims = extrema(Ms ./ t2),
                  legend = :outerright)
    for (_, c, name) in phases                                  ##### legend swatches for the three phases
        scatter!(plt, [NaN], [NaN] ; marker = :square, markersize = 8, color = c, markerstrokewidth = 0.5, label = name)
    end
    savefig_readme(plt, "haldane_chern_phase_diagram")
end

function figure_haldane_curvature()
    uc      =   haldane_unitcell(0.7, 0.2, π/2)            ##### C = 1, close to the transition at M = 3√3 t2
    bz      =   BZ([180, 180])
    FillBZ!(bz, uc)
    H       =   Hamiltonian(uc, bz)
    DiagonalizeHamiltonian!(H ; verbose = false)
    F       =   Curvature(H, [1])                          ##### Berry flux through each plaquette of the k-grid
    cellarea=   abs(det(hcat(bz.basis...))) / prod(bz.gridSize)
    Ω       =   F ./ cellarea
    ks, Z   =   zone_map(k -> grid_lookup(Ω, bz, k), uc ; n = 501)
    lim     =   maximum(abs, filter(!isnan, Z))
    plt = heatmap(ks, ks, Z ; color = DIVERGING, clims = (-lim, lim), aspect_ratio = :equal, colorbar = false,
                  xlabel = L"k_x", ylabel = L"k_y", xlims = extrema(ks), ylims = extrema(ks))
    draw_zone!(plt, uc)
    println("Haldane curvature: Chern number of band 1 = ", round(ChernNumber(H, [1]), digits = 4))
    savefig_readme(with_colorbar(plt, (-lim, lim), L"\Omega(\mathbf{k})" ; cmap = DIVERGING, size = (540, 420)), "haldane_berry_curvature")
end

##### -------------------------------------------------------------------------------------------
##### 7. Static susceptibility of the square lattice with next-nearest-neighbour hopping
##### -------------------------------------------------------------------------------------------
function figure_susceptibility()
    uc      =   UnitCell([[1.0, 0.0], [0.0, 1.0]], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddIsotropicBonds!(uc, 1.0,       -1.0,  "t1")
    AddIsotropicBonds!(uc, sqrt(2.0),  0.3,  "t2")
    bz      =   BZ([80, 80])
    FillBZ!(bz, uc)
    H       =   Hamiltonian(uc, bz)
    DiagonalizeHamiltonian!(H ; verbose = false)
    model   =   Model(uc, bz, H ; filling = 0.42, T = 0.02)
    SolveModel!(model ; get_correlations = false, verbose = false)

    χ       =   zeros(size(bz.ks))
    Threads.@threads for I in CartesianIndices(bz.ks)
        χ[I] = real(FindReducedChi(0.0, bz.ks[I], model ; a = 1, b = 1, eta = 1e-3)) / prod(bz.gridSize)
    end
    ##### At q = 0 exactly, the static intraband term (f(ε) - f(ε)) / 0 is dropped; use the q -> 0 limit instead.
    Γ       =   CartesianIndex(Tuple(GetQIndex([0.0, 0.0], bz)))
    χ[Γ]    =   sum(χ[mod1(Γ[1] + di, end), mod1(Γ[2] + dj, end)] for (di, dj) in ((1, 0), (-1, 0), (0, 1), (0, -1))) / 4
    ks, Z   =   zone_map(k -> grid_lookup(χ, bz, k), uc ; n = 401)
    clims   =   extrema(filter(!isnan, Z))
    plt = heatmap(ks ./ π, ks ./ π, Z ; color = SEQUENTIAL, clims = clims, aspect_ratio = :equal, colorbar = false,
                  xlabel = L"q_x / \pi", ylabel = L"q_y / \pi", xlims = (-1, 1), ylims = (-1, 1))
    savefig_readme(with_colorbar(plt, clims, L"|t|\,\chi_0(\mathbf{q}, 0)" ; size = (540, 420)), "square_susceptibility")
end

##### -------------------------------------------------------------------------------------------
##### 8. Hofstadter butterfly on a real-space square lattice in the Landau gauge
##### -------------------------------------------------------------------------------------------
function figure_hofstadter()
    Lx, Ly  =   6, 200
    αs, Es  =   Float64[], Float64[]
    for n in 0:Ly
        uc  =   UnitCell([[1.0, 0.0], [0.0, 1.0]], 1)
        AddBasisSite!(uc, [0.0, 0.0])
        AddIsotropicBonds!(uc, 1.0, -1.0, "t1")
        lat =   Lattice(uc, [Lx, Ly])
        FillLattice!(lat)
        B   =   2π * n / Ly                                ##### flux per plaquette in units of the flux quantum: n / Ly
        ModifyGauge!(lat, LandauGauge ; _kwargs = Dict{Symbol, Any}(:B => B))
        E   =   eigvals(Hermitian(FillHamiltonian(lat)))
        append!(Es, E); append!(αs, fill(n / Ly, length(E)))
    end
    plt = scatter(αs, Es ; markersize = 0.6, markerstrokewidth = 0, color = COLORS.ink, alpha = 0.6, label = "",
                  xlabel = L"\Phi / \Phi_0", ylabel = L"E / |t|", xlims = (0, 1), ylims = (-4.2, 4.2), size = (540, 420))
    savefig_readme(plt, "hofstadter_butterfly")
end

##### -------------------------------------------------------------------------------------------
const FIGURES = Dict(
    "unitcell"     => figure_unitcell,
    "bands"        => figure_kagome_bands,
    "dos"          => figure_kagome_dos,
    "fermi"        => figure_fermi_surfaces,
    "chern"        => figure_haldane_phase_diagram,
    "curvature"    => figure_haldane_curvature,
    "chi"          => figure_susceptibility,
    "hofstadter"   => figure_hofstadter,
)

##### run all figures, or only those named on the command line
if abspath(PROGRAM_FILE) == @__FILE__
for name in (isempty(ARGS) ? sort(collect(keys(FIGURES))) : ARGS)
    @time FIGURES[name]()
end
end
