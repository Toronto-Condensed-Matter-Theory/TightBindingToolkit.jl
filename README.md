# TightBindingToolkit

[![Build Status](https://https://github.com/Toronto-Condensed-Matter-Theory/TightBindingToolkit.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://https://github.com/Toronto-Condensed-Matter-Theory/TightBindingToolkit.jl/actions/workflows/CI.yml?query=branch%3Amain)

TightBindingToolkit.jl is a Julia package meant for constructing, and obtaining useful properties of generic tight-binding models. It supports any lattice structure, with any user-defined bonds on that lattice. It also has support for any spin of the particle hopping on the lattice. 
In the Julia REPL, Push "]" to enter the package mode.
```julia
add TightBindingToolkit
```
Or equivalently
```julia
using Pkg
Pkg.add("TightBindingToolkit")
```

Until the package is registered in the General registry, install via URL:
```julia
using Pkg
Pkg.add(url="https://https://github.com/Toronto-Condensed-Matter-TheoryTightBindingToolkit.jl")
```

## Features
Currently supported :
* Custom Unit Cell Construction. (v1.3) : now also supports construction in dimensions up to $d=3$. Can change primitives, or expand unit cell on the go.
    |![Alt text](./docs/src/Figures/Honeycomb_123NN_UC.png)|
    |:--:| 
    |*Honeycomb Model with 1st, 2nd and 3rd neighbour hoppings.*|
* Corresponding Brillouin Zone Construction. (v1.3) : now also supports construction in dimensions up to $d=3$.
* Hamiltonian, given a Unit Cell and a Brillouin Zone. (v1.3) : now also support BdG Hamiltonians.
* Diagonalizing the Hamiltonian in momentum space to get band structures and wavefunctions.
    |![Alt text](./docs/src/Figures/Triangle_123NN_bandStructure.png)|
    |:--:| 
    |*Band structure for a triangular lattice with 1st, 2nd and 3rd neighbour hopping.*|
* Density of State 
    |![Alt text](./docs/src/Figures/Graphene_DOS.png)|
    |:--:| 
    |*Density of state of graphene.*|
* Filling the model at given chemical potential, and calculating gaps. (v1.3) : also supported for BdG systems.
    |![Alt text](./docs/src/Figures/KitaevChain.png)|
    |:--:| 
    |*Band structure of the Kitaev Chain in 1-d.*|
* Fermi surfaces of systems
    |![Alt text](./docs/src/Figures/Triangle_123NN_FS.png)|
    |:--:| 
    |*Fermi surface for a triangular lattice with 1st, 2nd and 3rd neighbour hopping.*|
* Getting correlation functionsin momentum space and real space. (v1.3) : also supported for BdG systems.
* Getting Berry curvature and Chern numbers. Also workjs for topological superconductors.
    |![Alt text](./docs/src/Figures/Haldane_Chern.png)|
    |:--:| 
    |*Chern numbers in a spin-ful Haldane model as a function of t2.*|
* Getting magnetic susceptibility in any direction, at any momentum, and energy.
    |![Alt text](./docs/src/Figures/Triangle_123NN_chi.png)|
    |:--:| 
    |*imaginary part of zero-energy susceptibility for a triangular lattice with 1st, 2nd and 3rd neighbour hopping.*|
* Real-sapce lattice construction with arbitrary boundary condition. (v2.1)

## Documentation
For more details, please see [Documentation](https://Toronto-Condensed-Matter-Theory.github.io/TightBindingToolkit.jl/dev/)
