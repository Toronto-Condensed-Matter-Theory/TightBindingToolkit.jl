##### The Bloch Hamiltonian: FillHoppingHamiltonian, FillPairingHamiltonian, FillHamiltonian and FourierHamiltonian
##### must all agree with the block-by-block construction below (the v3.2 implementation).

using Random

function reference_hopping(uc, k)
    L       =   uc.localDim
    H       =   zeros(ComplexF64, L * length(uc.basis), L * length(uc.basis))
    for site in 1:length(uc.basis)
        b   =   L * (site - 1) + 1
        H[b : b + L - 1, b : b + L - 1]    .-=  sum(uc.fields[site] .* uc.OnSiteMats)
    end
    for bond in uc.bonds
        b1  =   L * (bond.base - 1) + 1
        b2  =   L * (bond.target - 1) + 1
        if b1 == b2 && iszero(bond.offset)
            H[b1 : b1 + L - 1, b1 : b1 + L - 1]    .+=  (bond.mat + bond.mat') / 2
        else
            phase   =   cis(dot(k, sum(bond.offset .* uc.primitives)))
            H[b1 : b1 + L - 1, b2 : b2 + L - 1]    .+=  phase .* bond.mat
            H[b2 : b2 + L - 1, b1 : b1 + L - 1]    .+=  conj(phase) .* bond.mat'
        end
    end
    return H
end

function reference_pairing(uc, k)
    L       =   uc.localDim
    H       =   zeros(ComplexF64, L * length(uc.basis), L * length(uc.basis))
    for bond in uc.bonds
        b1  =   L * (bond.base - 1) + 1
        b2  =   L * (bond.target - 1) + 1
        H[b1 : b1 + L - 1, b2 : b2 + L - 1]    .+=  cis(dot(k, sum(bond.offset .* uc.primitives))) .* bond.mat
    end
    return H
end

function reference_bdg(uc_hop, uc_pair, k)
    T   =   reference_hopping(uc_hop, k)
    Tm  =   reference_hopping(uc_hop, -k)
    Δ   =   reference_pairing(uc_pair, k) - transpose(reference_pairing(uc_pair, -k))
    return [T Δ; Δ' -transpose(Tm)]
end

##### A honeycomb cell with every kind of bond: on-site, intra-cell, same site in another cell, repeated bonds,
##### and complex matrices that are not hermitian.
function messy_cells(rng)
    uc      =   UnitCell([[1.0, 0.0], [0.5, sqrt(3) / 2]], 2, 2)
    AddBasisSite!(uc, [0.0, 0.0])
    AddBasisSite!(uc, [0.5, 1 / (2 * sqrt(3))])
    uc.fields   =   [randn(rng, 4) for _ in uc.basis]
    bonds   =   [(1, 1, [0, 0]), (2, 2, [0, 0]), (1, 2, [0, 0]), (1, 2, [-1, 0]), (1, 2, [0, -1]),
                 (1, 1, [1, 0]), (2, 2, [0, 1]), (1, 2, [-1, 0]), (2, 1, [1, -1]), (1, 1, [2, -1])]
    up      =   deepcopy(uc)
    for (base, target, offset) in bonds
        push!(uc.bonds, Bond(base, target, offset, randn(rng, ComplexF64, 2, 2), 1.0, "t"))
        push!(up.bonds, Bond(base, target, offset, randn(rng, ComplexF64, 2, 2), 1.0, "Δ"))
    end
    return uc, up
end

@testset "Bloch Hamiltonian" begin
    rng     =   MersenneTwister(7)
    uc, up  =   messy_cells(rng)
    ks      =   [randn(rng, 2) .* 4 for _ in 1:10]

    @testset "single momentum" begin
        for k in ks
            H   =   FillHoppingHamiltonian(uc, k)
            @test H ≈ reference_hopping(uc, k) atol = 1e-12
            @test H ≈ H' atol = 1e-12
            @test FillPairingHamiltonian(up, k) ≈ reference_pairing(up, k) atol = 1e-12
            HB  =   FillHamiltonian(uc, up, k)
            @test HB ≈ reference_bdg(uc, up, k) atol = 1e-12
            @test HB ≈ HB' atol = 1e-12
        end
    end

    @testset "explicit displacements" begin
        deltas  =   TightBindingToolkit.Hams.BondDisplacements(uc)
        onsite  =   TightBindingToolkit.Hams.OnSiteBlocks(uc)
        @test FillHoppingHamiltonian(uc, ks[1], deltas, onsite) ≈ reference_hopping(uc, ks[1]) atol = 1e-12
        @test_throws DimensionMismatch FillHoppingHamiltonian(uc, ks[1], deltas[1:end-1], onsite)
    end

    @testset "FourierHamiltonian" begin
        F   =   FourierHamiltonian(uc)
        FB  =   FourierHamiltonian(uc, up)
        @test size(F) == (4, 4)
        @test size(FB) == (8, 8)
        ##### repeated bonds such as (1, 2, [-1, 0]) merge; the displacements are the distinct ±offsets
        offsets =   [b.offset for b in uc.bonds if !iszero(b.offset)]
        @test size(F.displacements, 2) == length(unique(vcat(offsets, .-offsets))) == 8
        @test length(F.vals) < 2 * 4 * count(b -> !iszero(b.offset), uc.bonds)
        H   =   zeros(ComplexF64, 4, 4)
        for k in ks
            @test FillHamiltonian(F, k) ≈ reference_hopping(uc, k) atol = 1e-12
            @test FillHamiltonian!(H, F, k) === H
            @test H ≈ reference_hopping(uc, k) atol = 1e-12
            @test FillHamiltonian(FB, k) ≈ reference_bdg(uc, up, k) atol = 1e-12
        end
        @test_throws DimensionMismatch FillHamiltonian!(zeros(ComplexF64, 3, 3), F, ks[1])
        @test_throws DimensionMismatch FillHamiltonian(F, [0.1, 0.2, 0.3])
    end

    @testset "Hamiltonian on a BZ grid" begin
        bz  =   BZ([6, 9])
        FillBZ!(bz, uc)
        H   =   Hamiltonian(uc, bz)
        @test size(H.H) == size(bz.ks)
        @test all(H.H[i] ≈ reference_hopping(uc, bz.ks[i]) for i in eachindex(bz.ks))
        HB  =   Hamiltonian(uc, up, bz)
        @test HB.is_BdG
        @test all(HB.H[i] ≈ reference_bdg(uc, up, bz.ks[i]) for i in eachindex(bz.ks))
    end

    @testset "invalid bonds" begin
        bad     =   deepcopy(uc)
        push!(bad.bonds, Bond(1, 3, [0, 1], zeros(ComplexF64, 2, 2), 1.0, "bad"))
        @test_throws ArgumentError FillHoppingHamiltonian(bad, ks[1])
        @test_throws ArgumentError FourierHamiltonian(bad)
        bad     =   deepcopy(uc)
        push!(bad.bonds, Bond(1, 2, [0, 1], zeros(ComplexF64, 3, 3), 1.0, "bad"))
        @test_throws DimensionMismatch FillHoppingHamiltonian(bad, ks[1])
    end
end
