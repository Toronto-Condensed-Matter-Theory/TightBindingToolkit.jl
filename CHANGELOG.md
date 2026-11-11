# Changelog

## 3.2.0

### Changed

- **Real-space Hamiltonian (`FillHamiltonian(::Lattice{2})`, `LatticeHamiltonian`, `LatticeModel`): on-site fields now enter with a minus sign**, `H_ii -= ∑ fields[i] .* OnSiteMats`, matching the momentum-space `Hamiltonian`/`FillHoppingHamiltonian` and the BdG field code. Previously the real-space code added them, so a lattice with nonzero `fields` (Zeeman fields or the chemical-potential entry) had the opposite on-site energy from the same `UnitCell` in momentum space. On-site bonds (zero distance) are now symmetrized as `(M + M')/2`, as in momentum space, instead of adding `M` alone. `LatticeModel` users with nonzero fields or non-Hermitian on-site bond matrices will see different spectra; flip the sign of the fields to recover the old behaviour.
- Field vectors are no longer shared by reference. `Lattice` sites, `AddBasisSite!`, `ModifyFields!` and `ModifyIsotropicFields!` now store copies, so changing one site's field no longer changes every site of that sublattice or the `UnitCell`. Changing `uc.fields` after `FillLattice!` no longer changes `lattice.fields`; refill or edit `lattice.fields` directly.
- `FillBonds!` applies the boundary-condition phase only to rank-2 bonds. Rank-4 (interaction) lattices with twisted or antiperiodic boundary conditions no longer get the phase multiplied onto bonds that wrap around the boundary.

### Fixed

- `ModifyIsotropicFields!(uc, ::Vector{Float64})` threw a `MethodError`; it now sets the same field (copied) on every sublattice.
