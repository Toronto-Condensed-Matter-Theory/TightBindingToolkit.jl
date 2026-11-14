# Changelog

## 3.3.0

### Added

- `FourierHamiltonian(uc)` and `FourierHamiltonian(uc_hop, uc_pair)`: the k-independent form of the hopping or BdG Bloch Hamiltonian. Every bond and its hermitian conjugate are merged into one entry per (row, column, displacement), so the phases are computed once per distinct displacement instead of once per bond. Evaluate it with `FillHamiltonian(F, k)` or in place with `FillHamiltonian!(H, F, k)`. It is a snapshot of the unit cell: rebuild it after the bonds or fields change.

### Changed

- `Hamiltonian(uc, bz)` and `Hamiltonian(uc_hop, uc_pair, bz)` build a `FourierHamiltonian` once and evaluate it at every momentum. For the twisted-WSe2 Wannier model (6 × 6, 1317 bonds) a 30 × 30 grid takes about 20 ms instead of 120 ms (BdG: 50–70 ms instead of 430 ms); for its √3 × √3 supercell (18 × 18, 3951 bonds), about 30 ms instead of 400 ms. MeanFieldToolkit rebuilds the Hamiltonian this way at every iteration.
- `FillHoppingHamiltonian`, `FillPairingHamiltonian` and the single-momentum BdG `FillHamiltonian(uc_hop, uc_pair, k)` write the bonds element by element and compute the bond phases from k · aⱼ, without allocating displacement vectors or calling BLAS for two-component dot products: 60–70 µs instead of about 460 µs per call on the same model. Explicit `deltas` and `onsite` arguments are still accepted. Results agree with 3.2 to rounding (≲ 1e-13 relative).
- Momenta may be any `AbstractVector{<:Real}`.
- A bond to a site outside the basis, or a bond matrix that is not `localDim × localDim`, now throws an error instead of indexing out of range.

## 3.2.0

### Changed

- **Real-space Hamiltonian (`FillHamiltonian(::Lattice{2})`, `LatticeHamiltonian`, `LatticeModel`): on-site fields now enter with a minus sign**, `H_ii -= ∑ fields[i] .* OnSiteMats`, matching the momentum-space `Hamiltonian`/`FillHoppingHamiltonian` and the BdG field code. Previously the real-space code added them, so a lattice with nonzero `fields` (Zeeman fields or the chemical-potential entry) had the opposite on-site energy from the same `UnitCell` in momentum space. On-site bonds (zero distance) are now symmetrized as `(M + M')/2`, as in momentum space, instead of adding `M` alone. `LatticeModel` users with nonzero fields or non-Hermitian on-site bond matrices will see different spectra; flip the sign of the fields to recover the old behaviour.
- Field vectors are no longer shared by reference. `Lattice` sites, `AddBasisSite!`, `ModifyFields!` and `ModifyIsotropicFields!` now store copies, so changing one site's field no longer changes every site of that sublattice or the `UnitCell`. Changing `uc.fields` after `FillLattice!` no longer changes `lattice.fields`; refill or edit `lattice.fields` directly.
- `FillBonds!` applies the boundary-condition phase only to rank-2 bonds. Rank-4 (interaction) lattices with twisted or antiperiodic boundary conditions no longer get the phase multiplied onto bonds that wrap around the boundary.

### Fixed

- `ModifyIsotropicFields!(uc, ::Vector{Float64})` threw a `MethodError`; it now sets the same field (copied) on every sublattice.
