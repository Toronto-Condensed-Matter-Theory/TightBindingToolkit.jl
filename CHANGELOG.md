# Changelog

## Unreleased

### Changed

- **Real-space Hamiltonian (`FillHamiltonian(::Lattice{2})`, `LatticeHamiltonian`, `LatticeModel`): on-site fields now enter with a minus sign**, `H_ii -= ∑ fields[i] .* OnSiteMats`, matching the momentum-space `Hamiltonian`/`FillHoppingHamiltonian` and the BdG field code. Previously the real-space code added them, so a lattice with nonzero `fields` (Zeeman fields or the chemical-potential entry) had the opposite on-site energy from the same `UnitCell` in momentum space. On-site bonds (zero distance) are now symmetrized as `(M + M')/2`, as in momentum space, instead of adding `M` alone. `LatticeModel` users with nonzero fields or non-Hermitian on-site bond matrices will see different spectra; flip the sign of the fields to recover the old behaviour.
