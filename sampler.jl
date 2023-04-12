include("gas.jl")
include("system.jl")

struct Sampler
    #raw data
    t::Vector{Float64}
    sys::Vector{System}

    #calculated data
    E::Vector{Vector{Float64}}
end
Sampler() = Sampler([], [], [])

function sample!(sampler::Sampler, t, sys)
    push!(sampler.t, t)
    push!(sampler.sys, deepcopy(sys))
    push!(sampler.E, 0.5*sys.gas.m*(sys.v[1,:].^2 .+ sys.v[2,:].^2 .+ sys.v[3,:].^2))
end