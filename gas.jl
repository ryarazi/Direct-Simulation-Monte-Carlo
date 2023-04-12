struct Gas
    d::Float64 #gas diameter
    m::Float64 #gas mass 
end

ρ_to_n(ρ, gas::Gas) = ρ/gas.m
n_to_ρ(n, gas::Gas) = n*gas.m

const argon = Gas(3.66e-10, 6.63e-26)