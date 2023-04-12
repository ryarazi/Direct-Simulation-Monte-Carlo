using LinearAlgebra

include("visual.jl")
include("sampler.jl")
include("gas.jl")
include("system.jl")

function main(;sys::System, Nsteps)
    #collision parameter
    v_init = abs(sys.v[1,1])
    vr_max = 5*v_init

    dt = 0.09*minimum(sys.dL)/v_init

    samp = Sampler()

    #main loop
    t=0.
    i_loop=0
    while i_loop<Nsteps
        sys.r += dt * sys.v

        #check for specular reflections
        for ip in 1:sys.Nsim #ip = i_particle
            for axis in 1:3
                if sys.r[axis, ip] >= sys.L[axis]
                    sys.r[axis, ip]=2*sys.L[axis]-sys.r[axis, ip]
                    sys.v[axis, ip]*=-1
                elseif sys.r[axis, ip] <= 0.
                    sys.r[axis, ip]=-sys.r[axis, ip]
                    sys.v[axis, ip]*=-1
                end
            end
        end

        #index particles to cells
        particles_in_cells = reshape([Int64[] for _ in 1:sys.Ncell_tot], sys.Ncell)

        particle_cell = similar(sys.r, Int64)
        for ip in 1:sys.Nsim #ip = i_particle
            particle_cell[:, ip] = floor.(Int, sys.r[:, ip]./sys.dL[:]) .+ 1
        end

        for ip in 1:sys.Nsim
            ix, iy, iz = particle_cell[:, ip]
            push!(particles_in_cells[ix,iy,iz], ip)
        end

        #calculate collisions
        for ix in 1:sys.Ncell[1], iy in 1:sys.Ncell[2], iz in 1:sys.Ncell[3]
            particle_list = particles_in_cells[ix,iy,iz]
            n_particles = length(particle_list)
            if n_particles<2 continue end #not enough cells for collisions
            
            n_candidates = n_particles^2 * sys.Neff * pi * sys.gas.d^2 * vr_max * dt / 2. / sys.Volcell
            for _ in 1:n_candidates
                p1 = rand(1:n_particles)
                p2 = rand(1:n_particles)
                while p1==p2 p2=rand(1:n_particles) end
                i1 = particle_list[p1]
                i2 = particle_list[p2]

                v1 = sys.v[:, i1]
                v2 = sys.v[:, i2]
                v_relative = v1.-v2
                v_relative_norm = norm(v_relative)

                if v_relative_norm/vr_max < rand() continue end #accept/reject

                v_center_mass = 0.5.*(v1.+v2)
                phi = 2*pi*rand()
                cos_phi = cos(phi)
                sin_phi = sin(phi)
                cos_theta = 2*rand()-1
                sin_theta = sqrt(1-cos_theta^2)
                v_relative_new = v_relative_norm .* (sin_theta*cos_phi, sin_theta*sin_phi, cos_theta)

                sys.v[:, i1] = v_center_mass .+ 0.5.*v_relative_new
                sys.v[:, i2] = v_center_mass .- 0.5.*v_relative_new
            end
        end

        #sample
        if i_loop%100 == 0
            sample!(samp, t, sys)
        end

        t += dt
        i_loop += 1
    end
    
    #plot
    #plot_system(x,y,z, xlim=(0,Lx),ylim=(0,Ly),zlim=(0,Lz))
    #animate_system(samp, xlim=(0,Lx),ylim=(0,Ly),zlim=(0,Lz))
    #plot_energy(samp)
    plot_thermalization(samp)
end

sys = System((1.e-6, 1.e-6, 1.e-6), (1, 1, 1), argon, 1.78, 2000, 400.)
main(sys=sys,
     Nsteps=10000)