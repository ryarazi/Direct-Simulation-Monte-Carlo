include("gas.jl")
include("system.jl")

struct Sampler
    #raw data
    t::Vector{Float64}
    sys::Vector{System}
    sysProp::Vector{SystemProperties}

end
Sampler() = Sampler([], [], [])

function sample!(sampler::Sampler, t, sys)
    push!(sampler.t, t)
    push!(sampler.sys, deepcopy(sys))
    push!(sampler.sysProp, SystemProperties(sys))
end