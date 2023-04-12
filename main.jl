using LinearAlgebra

include("visual.jl")
include("sampler.jl")
include("gas.jl")
include("system.jl")
include("boundary.jl")

function main(;sys::System,
               Nsteps, Nsample,
               boundary::T=SpecularWall()) where {T<:BoundaryCondition}
    #collision parameter
    vr_max = 5*maximum(sys.v)

    dt = 0.1*minimum(sys.dL)/vr_max

    samp = Sampler()

    #main loop
    t=0.
    i_loop=0
    while i_loop<Nsteps
        sys.r += dt * sys.v

        apply_boundary(sys, boundary)

        #index particles to cells
        cell_indices = LinearIndices(sys.Ncell)
        particles_in_cells = [Int64[] for _ in 1:sys.Ncell_tot]

        particle_cell = floor.(Int, sys.r./sys.dL) .+ 1
        particle_cell_linear = [cell_indices[xyz_index...] for xyz_index in eachcol(particle_cell)]

        for ip in 1:sys.Nsim
            push!(particles_in_cells[particle_cell_linear[ip]], ip)
        end

        #calculate collisions
        for particle_list in particles_in_cells
            n_particles = length(particle_list)
            if n_particles<2 continue end #not enough particles for collisions
            
            candidates = n_particles^2 * sys.Neff * pi * sys.gas.d^2 * vr_max * dt / 2. / sys.Volcell
            for _ in 1:candidates
                p1 = rand(1:n_particles)
                p2 = rand(1:n_particles)
                while p1==p2 p2=rand(1:n_particles) end
                ip1 = particle_list[p1]
                ip2 = particle_list[p2]

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

sys = System((1.e-6, 1.e-6, 1.e-6), (20, 1, 1), argon, 1.78, 2000, 400.)
main(sys=sys,
     Nsteps=2000)