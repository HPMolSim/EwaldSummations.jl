struct Ewald2DInteraction{T} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    α::T
    r_c::T
    k_c::T
    L::NTuple{3,T}
    k_set::Vector{Tuple{T, T, T}}
end

struct Ewald3DInteraction{T} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    ϵ_inf::T # boundary condition at infinity, default is Inf (conductor)
    α::T
    r_c::T
    k_c::T
    L::NTuple{3,T}
    k_set::Vector{Tuple{T, T, T, T}}
end