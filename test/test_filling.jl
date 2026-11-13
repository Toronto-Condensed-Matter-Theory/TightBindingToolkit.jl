##### The chemical potential must be found even when every band lies on one side of zero
##### (e.g. a large uniform on-site energy, as from a Hartree shift). Before v3.2 GetMu! searched 10 .* bandwidth
##### (BdG: 2 .* bandwidth), which then excludes μ, and every level was filled.

@testset "Chemical potential with all bands above zero" begin
    SpinVec =   SpinMats(1 // 2)
    shift   =   10.0

    function shifted_square()
        uc  =   UnitCell([[1.0, 0.0], [0.0, 1.0]], 2, 2)
        AddBasisSite!(uc, [0.0, 0.0])
        t   =   Param(-1.0, 2)
        AddIsotropicBonds!(t, uc, 1.0, SpinVec[4], "t")
        ε   =   Param(shift, 2)
        AddAnisotropicBond!(ε, uc, 1, 1, [0, 0], SpinVec[4], 0.0, "on-site energy")
        CreateUnitCell!(uc, [t, ε])
        return uc
    end

    @testset "Model" begin
        uc  =   shifted_square()
        bz  =   BZ([24, 24])
        FillBZ!(bz, uc)
        H   =   Hamiltonian(uc, bz)
        DiagonalizeHamiltonian!(H ; verbose = false)
        @test minimum(reduce(vcat, H.bands)) > 0            ##### the situation that used to fail
        for filling in (0.2, 0.5, 0.8)
            M   =   Model(uc, bz, H ; T = 0.1, filling = filling)
            SolveModel!(M ; mu_tol = 1e-10, filling_tol = 1e-10, verbose = false)
            @test FindFilling(M.mu, M.Ham.bands, M.T, M.stat) ≈ filling atol = 1e-8
            @test shift - 4 < M.mu < shift + 4
        end
    end

    @testset "BdGModel" begin
        hop     =   shifted_square()
        pair    =   UnitCell([[1.0, 0.0], [0.0, 1.0]], 2, 2)
        AddBasisSite!(pair, [0.0, 0.0])
        Δ       =   Param(0.2, 2)
        AddIsotropicBonds!(Δ, pair, 1.0, ComplexF64[0 1; -1 0], "s-wave singlet")
        CreateUnitCell!(pair, [Δ])
        bz      =   BZ([24, 24])
        FillBZ!(bz, hop)
        H       =   Hamiltonian(hop, pair, bz)
        DiagonalizeHamiltonian!(H ; verbose = false)
        for filling in (0.3, 0.5)
            M   =   BdGModel(hop, pair, bz, H ; T = 0.05, filling = filling)
            SolveModel!(M ; mu_tol = 1e-10, filling_tol = 1e-10, verbose = false)
            @test FindFilling(M) ≈ filling atol = 1e-6
            @test shift - 4 < M.mu < shift + 4

            ##### Solving again (as a self-consistency loop does): the Hamiltonian now carries μ ≈ 10 as its on-site field,
            ##### so the BdG bands ±|ε - μ| span only a few t around zero, and 2 .* bandwidth no longer contains μ.
            H2  =   Hamiltonian(hop, pair, bz)          ##### built with the fields written by the previous solve
            DiagonalizeHamiltonian!(H2 ; verbose = false)
            M2  =   BdGModel(hop, pair, bz, H2 ; T = 0.05, filling = filling)
            SolveModel!(M2 ; mu_tol = 1e-10, filling_tol = 1e-10, verbose = false)
            @test FindFilling(M2) ≈ filling atol = 1e-6
            @test M2.mu ≈ M.mu atol = 1e-6
        end
    end
end
