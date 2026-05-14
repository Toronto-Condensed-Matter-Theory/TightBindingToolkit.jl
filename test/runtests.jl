using TightBindingToolkit
using LinearAlgebra
using Test

@testset "TightBindingToolkit.jl" begin

    @testset "SpinMats" begin
        S = SpinMats(1//2)
        @test length(S) == 4                      # Sx, Sy, Sz, I
        @test size(S[1]) == (2, 2)
        @test S[4] ≈ [1.0 0.0; 0.0 1.0]          # identity
    end

    @testset "UnitCell construction" begin
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        UC = UnitCell([a1, a2], 2, 2)
        AddBasisSite!(UC, [0.0, 0.0])
        @test length(UC.basis) == 1
        @test UC.localDim == 2
    end

    @testset "Square lattice Hamiltonian" begin
        # Build a simple square lattice tight-binding model
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        UC = UnitCell([a1, a2], 2, 2)
        AddBasisSite!(UC, [0.0, 0.0])

        SpinVec = SpinMats(1//2)
        t1 = -1.0
        NNdist = 1.0
        AddIsotropicBonds!(UC, NNdist, t1 * SpinVec[4], "t1")

        kSize = 6 * 6 + 3  # small grid for testing
        bz = BZ([kSize, kSize])
        FillBZ!(bz, UC)

        H = Hamiltonian(UC, bz)
        DiagonalizeHamiltonian!(H)

        # Bands should exist and have the right shape
        @test length(H.bands) == kSize * kSize
        # Each k-point should have 2 bands (localDim=2, 1 site)
        @test length(H.bands[1]) == 2
        # Bandwidth should be finite
        all_energies = vcat(H.bands...)
        @test maximum(all_energies) - minimum(all_energies) > 0.0
    end

    @testset "BZ construction" begin
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        UC = UnitCell([a1, a2], 2, 2)
        AddBasisSite!(UC, [0.0, 0.0])

        bz = BZ([9, 9])
        FillBZ!(bz, UC)
        @test length(bz.ks) == 81
        @test length(bz.basis) == 2
    end

    @testset "Model at half-filling" begin
        a1 = [1.0, 0.0]
        a2 = [0.0, 1.0]
        UC = UnitCell([a1, a2], 2, 2)
        AddBasisSite!(UC, [0.0, 0.0])

        SpinVec = SpinMats(1//2)
        t1 = -1.0
        NNdist = 1.0
        AddIsotropicBonds!(UC, NNdist, t1 * SpinVec[4], "t1")

        kSize = 6 * 6 + 3
        bz = BZ([kSize, kSize])
        FillBZ!(bz, UC)

        H = Hamiltonian(UC, bz)
        DiagonalizeHamiltonian!(H)

        M = Model(UC, bz, H; filling=0.5)
        SolveModel!(M)
        # Chemical potential should be found
        @test isfinite(M.mu)
        # Filling should be close to 0.5
        @test abs(M.filling - 0.5) < 0.05
    end

    @testset "Param-based construction" begin
        a1 = [1/2, sqrt(3)/2]
        a2 = [-1/2, sqrt(3)/2]
        UC = UnitCell([a1, a2], 2, 2)
        AddBasisSite!(UC, [0.0, 0.0])
        AddBasisSite!(UC, [0.0, 1/sqrt(3)])

        SpinVec = SpinMats(1//2)
        t1Param = Param(-1.0, 2)
        NNdist = 1/sqrt(3)
        AddIsotropicBonds!(t1Param, UC, NNdist, SpinVec[4], "t1")
        CreateUnitCell!(UC, [t1Param])

        @test length(UC.bonds) > 0
    end

end
# Test change
