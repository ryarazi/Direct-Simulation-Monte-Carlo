include("system.jl")

abstract type BoundaryCondition end
struct SpecularWall <: BoundaryCondition end
struct Periodic <: BoundaryCondition end

function apply_boundary(sys::System, ::SpecularWall)
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
end

function apply_boundary(sys::System, ::Periodic)
    sys.r = mod.(sys.r, sys.L)
end