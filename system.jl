include("gas.jl")

mutable struct System
    #spatial parameters
    const L::Tuple{Vararg{Float64, 3}} #the size in each dimension
    const Vol::Float64 #volume of the system

    #spatial resolution
    const Ncell::Tuple{Vararg{Int64, 3}} #the number of cells in each dimension
    const Ncell_tot::Float64 #total number of cells
    const Volcell::Float64 #Volume of cell
    const dL::Tuple{Vararg{Float64, 3}} #the length of a cell in each dimensions

    #particles parameters
    const gas::Gas
    const Nsim::Int64 #number of simulation particles
    const Nphys::Float64 #number of physical particles
    const Neff::Float64 #number of physical particles for each simulation particles
    
    #particles locations+velocities
    r::Matrix{Float64}
    v::Matrix{Float64}

    function System(L::Tuple{Vararg{Float64, 3}}, Ncell::Tuple{Vararg{Int64, 3}},
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

