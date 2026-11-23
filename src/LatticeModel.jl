module LatModel

    export LatticeModel , GetMu! , GetFilling! , GetGr!, SolveModel!, GetGap!

    using ..TightBindingToolkit.Useful: Meshgrid, DistFunction, DeriDistFunction, FFTArrayofMatrix, BinarySearch
    using ..TightBindingToolkit.LatticeStruct: Lattice
    using ..TightBindingToolkit.LatHam: LatticeHamiltonian
    using ..TightBindingToolkit.TBModel: FindFilling, GetCount

    import ..TightBindingToolkit.TBModel: GetMu!, GetFilling!, GetGr!, SolveModel!, GetGap!, GetGroundStateEnergyDensity


    mutable struct LatticeModel
        lattice ::  Lattice{2}
        Ham     ::  LatticeHamiltonian
        """
        Thermodynamic properties
        """
        T       ::  Float64         ##### Temperature
        filling ::  Float64         ##### Filling fraction
        mu      ::  Float64         ##### Chemical potential
        stat    ::  Int64           ##### +1 for bosons, -1 for fermions
        gap     ::  Float64
        """
        Correlations
        """
        Gr      ::  Matrix{ComplexF64}

        ##### Chosen the default value of filling to be -1 which is unphysical so that the code knows when filling has not been provided and has to be calculated from mu instead!
        function LatticeModel(lattice::Lattice{2}, Ham::LatticeHamiltonian ; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)

            localDim    = lattice.uc.localDim
            Hamsize     = lattice.length * localDim
            return new{}(lattice, Ham, T, filling, mu, stat, -999.0, Matrix{ComplexF64}(undef, Hamsize, Hamsize))
        end
        ##### Automatically cosntruct the lattice hamiltonian if only the unit cell is passed #####
        function LatticeModel(lattice::Lattice{2}; T::Float64=1e-3, filling::Float64=-1.0, mu::Float64=0.0, stat::Int64=-1)

            localDim    = lattice.uc.localDim
            Hamsize     = lattice.length * localDim
            return new{}(lattice, LatticeHamiltonian(lattice), T, filling, mu, stat, -999.0, Matrix{ComplexF64}(undef, Hamsize, Hamsize))
        end
    end


    function GetMu!(M::LatticeModel ;  mu_tol::Float64 = 0.001, filling_tol::Float64 = 1e-6)
    
        if M.T≈0.0
            M.mu      =   (M.Ham.bands[floor(Int64, length(M.Ham.bands)*M.filling)] + M.Ham.bands[floor(Int64, length(M.Ham.bands)*M.filling)+1])/2
        else
            M.mu      =   BinarySearch(M.filling, M.Ham.bandwidth, FindFilling, (M.Ham.bands, M.T, M.stat) ; x_tol=mu_tol, target_tol = filling_tol)
            @info "Found chemical potential μ = $(M.mu) for given filling = $(M.filling)."
        end
    end


    function GetFilling!(M::LatticeModel)

        if M.T≈0.0
            M.filling   =   searchsortedlast(M.Ham.bands, M.mu) / length(M.Ham.bands)
        else
            M.filling   =   FindFilling( M.mu, M.Ham.bands, M.T, M.stat)
        end
    end


    function GetGr!(M::LatticeModel)

        quasiCount 	=	GetCount(M.Ham.bands, M.mu, M.T, M.stat)   ##### Matrix with 1s and 0s along the diagonal. The count of the quasiparticles at each k point determined by the bands and the chemical potential
        M.Gr        =   transpose(M.Ham.states * quasiCount * adjoint(M.Ham.states))
    end


    function GetGap!(M::LatticeModel)

        index     =   floor(Int64, length(M.Ham.bands)*M.filling)   
        M.gap     =   M.Ham.bands[clamp(index + 1, 1, length(M.Ham.bands))] - M.Ham.bands[clamp(index, 1, length(M.Ham.bands))]
    end

    function SolveModel!(M::LatticeModel ; get_correlations::Bool = true, get_gap::Bool = false, verbose::Bool = true, mu_tol::Float64 = 1e-3, filling_tol::Float64 = 1e-6)
        @assert M.Ham.is_BdG==false "BdG Hamiltonians should be solved using a BdGModel"

        ##### Check if Hamiltonian has already been diagonalized. If not diagonalize it.
        if isempty(M.Ham.bands)
            DiagonalizeHamiltonian!(M.Ham)
        end

        if M.filling<0    ##### Must imply that filling was not provided by user and hence needs to be calculated from given mu
            GetFilling!(M)
        else
            GetMu!(M ; mu_tol = mu_tol, filling_tol = filling_tol)
        end

        if get_gap
            GetGap!(M)
        end
    
        if get_correlations
            GetGr!(M)
        end
        
        if verbose
            @info "System Filled!"
        end
    end

    function GetGroundStateEnergyDensity(M::LatticeModel)
        @assert !isempty(M.Ham.bands) "Model hamiltonian has not been diagonalized!"
        ##### energy eigenstates are already sorted in real-space #####
        N = prod(M.lattice.size)
        ##### sum over states with energies below the chemical potential #####
        mu_ind = searchsortedfirst(M.Ham.bands, M.mu)
        if mu_ind > 1
            return @views sum(M.Ham.bands[1:mu_ind - 1]) / N
        else
            return 0.0
        end
    end


end