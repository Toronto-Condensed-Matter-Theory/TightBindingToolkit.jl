##### FindReducedChi must equal the Lindhard sum on any grid. Its default Γ index, (N + 1) ÷ 2,
##### used to be right only for odd N; on even grids every momentum was shifted by one grid point.

@testset "Bare susceptibility on even and odd grids" begin
    uc  =   UnitCell([[1.0, 0.0], [0.0, 1.0]], 1)
    AddBasisSite!(uc, [0.0, 0.0])
    AddIsotropicBonds!(uc, 1.0, -1.0, "t1")
    AddIsotropicBonds!(uc, sqrt(2.0), 0.3, "t2")
    ε(k)    =   real(FillHoppingHamiltonian(uc, k)[1, 1])

    @testset "grid $(Tuple(L))" for L in ([12, 12], [9, 9], [10, 7])
        bz  =   BZ(L)
        FillBZ!(bz, uc)
        H   =   Hamiltonian(uc, bz)
        DiagonalizeHamiltonian!(H ; verbose = false)
        M   =   Model(uc, bz, H ; filling = 0.4, T = 0.1)
        SolveModel!(M ; get_correlations = false, verbose = false)

        @test bz.ks[TightBindingToolkit.suscep.GammaIndex(bz)...] ≈ [0.0, 0.0] atol = 1e-12

        for Q in (bz.ks[3, 5], bz.ks[end, 2]), Ω in (0.0, 0.4)
            η       =   0.05
            f(e)    =   DistFunction(e ; T = M.T, mu = M.mu, stat = M.stat)
            lindhard=   -sum((f(ε(k)) - f(ε(k + Q))) / (Ω + im * η - (ε(k + Q) - ε(k))) for k in bz.ks)
            @test FindReducedChi(Ω, Q, M ; a = 1, b = 1, eta = η) ≈ lindhard rtol = 1e-8
        end
    end
end
