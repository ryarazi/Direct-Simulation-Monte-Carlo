include("gas.jl")

mutable struct System
    #spatial parameters
    const L::NTuple{3, Float64} #the size in each dimension
    const Vol::Float64 #volume of the system

    #spatial resolution
    const Ncell::NTuple{3, Int64} #the number of cells in each dimension
    const Ncell_tot::Float64 #total number of cells
    const Volcell::Float64 #Volume of cell
    const dL::NTuple{3, Float64} #the length of a cell in each dimensions

    #particles parameters
    const gas::Gas
    const Nsim::Int64 #number of simulation particles
    const Nphys::Float64 #number of physical particles
    const Neff::Float64 #number of physical particles for each simulation particles
    
    #particles locations+velocities
    r::Matrix{Float64}
    v::Matrix{Float64}

    function System(L::NTuple{3, Float64}, Ncell::NTuple{3, Int64},
                    gas::Gas, ρ::Real, Nsim::Integer, 
                    v_init::Real)
        n = ρ_to_n(ρ, gas)
        Vol = prod(L)
        Ncell_tot = prod(Ncell)
        dL = L./Ncell

        Nphys = n*Vol
        Neff = Nphys/Nsim

        r = rand(3, Nsim)
        for i in 1:3 r[i,:] *= L[i] end
        
        v = (rand(Bool, 3, Nsim) .* 2. .- 1.)*v_init
        new(L, Vol, Ncell, Ncell_tot, Vol/Ncell_tot, dL,
            gas, Nsim, Nphys, Neff, r, v)
    end
end


const kb = 1.380649e-23 #Boltzmann constant
struct SystemProperties
    E::Vector{Float64}
    Etot::Float64

    kT::Float64
    T::Float64

    function SystemProperties(sys::System)
        D = length(sys.L)

        E = vec(0.5*sys.gas.m*sum(sys.v.^2, dims=1) * sys.Neff)
        Etot = sum(E)
        kT = Etot/sys.Nphys * (2/D)
        T = kT/kb

        new(E, Etot, kT, T)
    end
end