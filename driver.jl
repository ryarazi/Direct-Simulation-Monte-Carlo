using LinearAlgebra

include("visual.jl")
include("sampler.jl")
include("gas.jl")
include("system.jl")
include("boundary.jl")

function simulate_dsmc(;sys::System,
               Nsteps, Nsample,
               boundary::T=SpecularWall()) where {T<:BoundaryCondition}
    #collision parameter
    vr_max = 5*maximum(sys.v)

    dt = 0.1*minimum(sys.dL)/vr_max

    samp = Sampler()

    #preallocate memory for arrays
    cc2cl = LinearIndices(sys.Ncell) #"cell cartesian to linear" - transform cartesian cell index to linear cell index
    p2cl = Vector{Float64}(undef, sys.Nsim) #"particle to cell linear" - index of particle to the index of it's cell
    sp2p = Vector{Int64}(undef, sys.Nsim) #"sorted particle to particle" - index of particles as they are sorted by their cell

    #main loop
    t=0.
    i_loop=0
    while i_loop<Nsteps
        sys.r += dt * sys.v

        apply_boundary(sys, boundary)

        #index particles to cells
        for ip in 1:sys.Nsim
            i = (floor(Int, sys.r[d,ip]/sys.dL[d]) + 1 for d in 1:length(sys.L))
            p2cl[ip] = cc2cl[i...]
        end
        sortperm!(sp2p, p2cl)

        #calculate collisions
        cell_start = 1
        cell_index = p2cl[sp2p[cell_start]]
        for i in 2:sys.Nsim+1
            if i!=sys.Nsim+1 && p2cl[sp2p[i]] == cell_index continue end

            cell_end = i-1
            cell_range = cell_start:cell_end
            
            if i!=sys.Nsim+1
                cell_start = i
                cell_index = p2cl[sp2p[cell_start]]
            end

            n_particles = length(cell_range)

            if n_particles<2 continue end #not enough particles for collisions
            
            candidates = n_particles^2 * sys.Neff * pi * sys.gas.d^2 * vr_max * dt / 2. / sys.Volcell
            for _ in 1:candidates

                #generate two different random particles in the cell
                i1 = rand(cell_range)
                i2 = rand(cell_range)
                while i1==i2 i2=rand(cell_range) end

                ip1 = sp2p[i1]
                ip2 = sp2p[i1]

                v1 = sys.v[:, ip1]
                v2 = sys.v[:, ip2]
                v_relative_norm = norm(v1.-v2)

                if v_relative_norm/vr_max < rand() continue end #accept/reject

                v_center_mass = 0.5.*(v1.+v2)
                phi = 2*pi*rand()
                cos_phi = cos(phi)
                sin_phi = sin(phi)
                cos_theta = 2*rand()-1
                sin_theta = sqrt(1-cos_theta^2)
                v_relative_new = v_relative_norm .* (sin_theta*cos_phi, sin_theta*sin_phi, cos_theta)

                sys.v[:, ip1] = v_center_mass .+ 0.5.*v_relative_new
                sys.v[:, ip2] = v_center_mass .- 0.5.*v_relative_new
            end
        end

        #sample
        if i_loop%Nsample == 0
            sample!(samp, t, sys)
        end

        t += dt
        i_loop += 1
    end

    return samp
end