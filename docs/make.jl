using Documenter
using TightBindingToolkit

makedocs(
    build       =   "build" ,
    sitename    =   "TightBindingToolkit.jl"    ,
    modules     =   [TightBindingToolkit, TightBindingToolkit.UCell, TightBindingToolkit.DesignUCell, TightBindingToolkit.ExpandUCell, TightBindingToolkit.Parameters, TightBindingToolkit.LatticeStruct, TightBindingToolkit.DesignLattice, TightBindingToolkit.BZone, TightBindingToolkit.Hams, TightBindingToolkit.TBModel, TightBindingToolkit.Chern, TightBindingToolkit.suscep, TightBindingToolkit.conduct]   ,
    pages = [
        "Introduction"              =>  "index.md",
        "Utilities"                 =>  "Utilities.md",
        "Unit Cell"                 =>  "UnitCell.md",
        "Parameters"                =>  "Params.md",
        "Lattice"                   =>  "Lattice.md",
        "Brillouin Zone"            =>  "BZ.md",
        "Hamiltonian"               =>  "Hamiltonian.md",
        "Tight Binding Model"       =>  "Model.md",
        "BdG Model"                 =>  "BdGModel.md",
        "Chern Numbers"             =>  "Chern.md",
        "Magnetic susceptibility"   =>  "susceptibility.md"  ,
        "Electrical Conductivity"   =>  "conductivity.md",
        "Plotting"                  =>  "Plot.md"
    ]
)

deploydocs(
    repo = "github.com/andrewkhardy/TightBindingToolkit.jl",
    devbranch = "main"
)
