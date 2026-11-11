##### Real-space Lattice tests: the real-space spectrum must equal the union of the
##### momentum-space spectra on the commensurate (and boundary-twisted) k-grid.

const SpinVec2 = SpinMats(1//2)

function SquareUnitCell(BC::Vector{<:Number})
    uc = UnitCell([[1.0, 0.0], [0.0, 1.0]], 2)
    AddBasisSite!(uc, [0.0, 0.0], [0.0, 0.0, 0.3, 0.7], HermitianBasis(2))
    ##### nearest-neighbour hopping, with a spin-dependent imaginary part
    AddIsotropicBonds!(uc, 1.0, -1.0 * SpinVec2[4] + 0.2im * SpinVec2[3], "t1")
    ##### on-site bond with a non-Hermitian matrix, so the (M + M')/2 symmetrization matters
    AddAnisotropicBond!(uc, 1, 1, [0, 0], 0.4 * SpinVec2[1] + 0.3im * SpinVec2[2], 0.0, "onsite")
    uc.BC = ComplexF64.(BC)
    return uc
end

function HoneycombUnitCell(BC::Vector{<:Number})
    uc = UnitCell([[1/2, sqrt(3)/2], [-1/2, sqrt(3)/2]], 2)
    AddBasisSite!(uc, [0.0, 0.0], [0.1, 0.0, 0.2, 0.5], HermitianBasis(2))
    AddBasisSite!(uc, [0.0, 1/sqrt(3)], [0.0, -0.1, -0.2, 0.4], HermitianBasis(2))
    AddIsotropicBonds!(uc, 1/sqrt(3), -1.0 * SpinVec2[4], "t1")
    AddIsotropicBonds!(uc, 1.0, 0.1im * SpinVec2[4] + 0.05 * SpinVec2[1], "t2")
    uc.BC = ComplexF64.(BC)
    return uc
end

function RealSpaceSpectrum(uc::UnitCell, L::Vector{Int64})
    lat = Lattice(uc, L)
    FillLattice!(lat)
    H = FillHamiltonian(lat)
    return sort(eigvals(Hermitian(H))), H
end

function MomentumSpaceSpectrum(uc::UnitCell, L::Vector{Int64})
    bz = BZ(L)
    FillBZ!(bz, uc)     ##### shifts the k-grid by the twist angle(uc.BC) / L
    H = Hamiltonian(uc, bz)
    DiagonalizeHamiltonian!(H ; verbose = false)
    return sort(reduce(vcat, vec(H.bands)))
end

@testset "Real-space Lattice" begin

    @testset "Square lattice: real space == k space, $name" for (name, BC) in (("PBC", [1, 1]), ("APBC in y", [1, -1]), ("twist 2π/3 in x", [cis(2π/3), 1]))
        uc = SquareUnitCell(BC)
        Er, H = RealSpaceSpectrum(uc, [4, 4])
        Ek = MomentumSpaceSpectrum(uc, [4, 4])
        @test ishermitian(H)
        @test length(Er) == length(Ek) == 4 * 4 * 2
        @test Er ≈ Ek atol = 1e-10
    end

    @testset "Honeycomb (two sublattices): real space == k space, $name" for (name, BC) in (("PBC", [1, 1]), ("APBC in y", [1, -1]), ("twist 2π/3 in y", [1, cis(2π/3)]))
        uc = HoneycombUnitCell(BC)
        Er, H = RealSpaceSpectrum(uc, [4, 4])
        Ek = MomentumSpaceSpectrum(uc, [4, 4])
        @test ishermitian(H)
        @test length(Er) == length(Ek) == 4 * 4 * 2 * 2
        @test Er ≈ Ek atol = 1e-10
    end

    @testset "On-site field sign matches k space" begin
        uc = UnitCell([[1.0, 0.0], [0.0, 1.0]], 2)
        AddBasisSite!(uc, [0.0, 0.0], [0.0, 0.0, 0.0, 0.7], HermitianBasis(2))
        lat = Lattice(uc, [2, 2])
        FillLattice!(lat)
        ##### last on-site matrix is the identity; fields enter as -field * mat
        @test FillHamiltonian(lat) ≈ -0.7 * Matrix{ComplexF64}(I, 8, 8)
    end

    @testset "Twisted boundary condition does not throw" begin
        uc = SquareUnitCell([cis(2π/3), 1])
        @test abs(cis(2π/3)) != 1.0      ##### the case that used to throw InexactError in Bool.(abs.(BC))
        lat = Lattice(uc, [3, 3])
        @test (FillLattice!(lat); true)
        @test ApplyBCToSite([3, 0], [3, 3], uc.BC) == [0, 0]
        @test ApplyBCToSite([3, 0], [3, 3], ComplexF64[0, 1]) == [3, 0]     ##### open along x
    end

    @testset "Lattice fields are not aliased" begin
        uc = SquareUnitCell([1, 1])
        lat = Lattice(uc, [2, 2])
        FillLattice!(lat)
        @test lat.fields[1] == lat.fields[2] == uc.fields[1]
        @test lat.fields[1] !== lat.fields[2]
        @test lat.fields[1] !== uc.fields[1]
        lat.fields[1][3] = 5.0
        @test lat.fields[2][3] == 0.3
        @test uc.fields[1][3] == 0.3

        field = [0.0, 0.0, 0.1, 0.2]
        uc2 = UnitCell([[1.0, 0.0], [0.0, 1.0]], 2)
        AddBasisSite!(uc2, [0.0, 0.0], field, HermitianBasis(2))
        AddBasisSite!(uc2, [0.5, 0.5], field, HermitianBasis(2))
        @test uc2.fields[1] !== field && uc2.fields[1] !== uc2.fields[2]

        ModifyIsotropicFields!(uc2, [0.0, 0.0, 0.0, 1.0])
        @test uc2.fields == [[0.0, 0.0, 0.0, 1.0], [0.0, 0.0, 0.0, 1.0]]
        @test uc2.fields[1] !== uc2.fields[2]
        ModifyIsotropicFields!(uc2, 0.5, 3)
        @test uc2.fields == [[0.0, 0.0, 0.5, 1.0], [0.0, 0.0, 0.5, 1.0]]
    end

    @testset "Boundary phase only on rank-2 bonds" begin
        uc4 = UnitCell([[1.0, 0.0], [0.0, 1.0]], 2, 4)
        AddBasisSite!(uc4, [0.0, 0.0])
        J = zeros(ComplexF64, 2, 2, 2, 2)
        J[1, 1, 2, 2] = 1.0
        AddIsotropicBonds!(uc4, 1.0, J, "J")
        uc4.BC = ComplexF64[1, -1]
        lat4 = Lattice(uc4, [3, 3])
        FillLattice!(lat4)
        @test all(m -> m ≈ J, lat4.bondMats[lat4.bondSites .!= 0])

        uc2 = SquareUnitCell([1, -1])
        lat2 = Lattice(uc2, [3, 3])
        FillLattice!(lat2)
        t1 = lat2.bondLabels .== "t1"
        ##### only the bonds wrapping around y (3 of the 18 "t1" bonds) pick up the -1
        @test count(m -> real(m[1, 1]) > 0, lat2.bondMats[t1]) == 3
    end

    @testset "DesignLattice functions" begin
        uc = UnitCell([[1.0, 0.0], [0.0, 1.0]], 2)
        AddBasisSite!(uc, [0.0, 0.0])
        t1 = Param(-1.0, 2)
        AddIsotropicBonds!(t1, uc, 1.0, SpinVec2[4], "t1")
        lat = CreateLattice(uc, [t1], [3, 3])
        @test all(m -> m ≈ -SpinVec2[4], lat.bondMats[lat.bondLabels .== "t1"])

        push!(t1.value, -2.0)
        ModifyLattice!(lat, [t1])
        @test all(m -> m ≈ -2 * SpinVec2[4], lat.bondMats[lat.bondLabels .== "t1"])

        ScaleLatticeBonds!(lat, "t1", 0.5)
        @test all(m -> m ≈ -SpinVec2[4], lat.bondMats[lat.bondLabels .== "t1"])

        RemoveLatticeBonds!(lat, "t1")
        @test !any(==("t1"), lat.bondLabels)
        @test all(iszero, lat.bondMats)
        @test all(==(-1.0), lat.bondDists)
    end

end
