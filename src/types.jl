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

struct IcmEwald2DInteraction{T} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    α::T
    r_c::T
    k_c::T
    γ::Tuple{T, T}
    L::NTuple{3,T}
    N_image::Int
    k_set::Vector{Tuple{T, T, T}}
end

struct IcmEwald3DInteraction{T} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    α::T
    r_c::T
    k_c::T
    γ::Tuple{T, T}
    L::NTuple{3,T}
    N_image::Int # layer of the image charges
    N_pad::Int # layer of the z padding, L_z_pad = (2 * N_pad + 1) * L_z
    k_set::Vector{Tuple{T, T, T, T}}
end