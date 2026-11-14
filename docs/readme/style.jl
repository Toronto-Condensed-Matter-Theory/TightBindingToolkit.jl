# Shared plot style for the README figures of the Toronto-Condensed-Matter-Theory toolkits.
# The same file is used in every package so the figures look alike.
ENV["GKSwstype"] = "100"   # headless GR
using Plots, LaTeXStrings
gr()

# Categorical colors, in a fixed order that stays distinguishable under colour-vision deficiency.
const COLORS = (
    blue    = colorant"#2a78d6",
    orange  = colorant"#eb6834",
    aqua    = colorant"#1baf7a",
    violet  = colorant"#6250d6",
    red     = colorant"#e34948",
    yellow  = colorant"#eda100",
    ink     = colorant"#1b1b1b",
    muted   = colorant"#8a8984",
)
const SERIES = [COLORS.blue, COLORS.orange, COLORS.aqua, COLORS.violet, COLORS.red, COLORS.yellow]

# Sequential map (magnitude, one hue light -> dark) and diverging map (sign, neutral midpoint).
const SEQUENTIAL = cgrad([colorant"#f7f9fc", colorant"#9cc3ee", colorant"#2a78d6", colorant"#123a6b"])
const DIVERGING  = cgrad([colorant"#2a78d6", colorant"#f2f1ee", colorant"#e34948"])

default(
    fontfamily          = "Computer Modern",
    framestyle          = :box,
    grid                = false,
    minorticks          = true,
    tick_direction      = :in,
    linewidth           = 2,
    palette             = SERIES,
    size                = (520, 380),
    dpi                 = 200,
    guidefontsize       = 13,
    tickfontsize        = 10,
    legendfontsize      = 10,
    colorbar_titlefontsize = 12,
    foreground_color_legend = nothing,
    background_color    = :white,
    left_margin         = 2Plots.mm,
    bottom_margin       = 2Plots.mm,
)

const FIGDIR = normpath(joinpath(@__DIR__, "..", "src", "assets", "readme"))
mkpath(FIGDIR)

"Save `plt` as `name.png` in the README figure folder and report the path."
function savefig_readme(plt, name::AbstractString)
    path = joinpath(FIGDIR, name * ".png")
    savefig(plt, path)
    println("wrote ", path)
    return path
end

"""
    with_colorbar(p, clims, label; cmap = SEQUENTIAL, size = (560, 380), widths = [0.9, 0.1])

`p` (drawn with `colorbar = false`) beside a colorbar for `clims` with `label` on its right-hand
axis. GR's own colorbar title collides with long tick labels, so the bar is a thin heatmap panel.
"""
function with_colorbar(p, clims, label; cmap = SEQUENTIAL, size = (560, 380), widths = [0.9, 0.1])
    v = collect(range(clims..., length = 256))
    bar = heatmap([0.0, 1.0], v, repeat(v, 1, 2); color = cmap, colorbar = false, clims,
                  xticks = false, ymirror = true, ylabel = label, minorticks = false, ylims = clims,
                  left_margin = -2Plots.mm, right_margin = 3Plots.mm)
    return plot(p, bar; layout = grid(1, 2; widths), size)
end
