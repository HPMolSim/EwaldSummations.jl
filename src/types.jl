struct Ewald2DInteraction{T} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    α::T
    r_c::T
    k_c::T
    L::NTuple{3,T}
    k_set::Vector{Tuple{T, T, T}}
end

Base.show(io::IO, interaction::Ewald2DInteraction) = print(io, "Ewald2DInteraction(n_atoms = $(interaction.n_atoms), ϵ = $(interaction.ϵ), α = $(interaction.α), r_c = $(interaction.r_c), k_c = $(interaction.k_c), L = $(interaction.L))")

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

Base.show(io::IO, interaction::Ewald3DInteraction) = print(io, "Ewald3DInteraction(n_atoms = $(interaction.n_atoms), ϵ = $(interaction.ϵ), ϵ_inf = $(interaction.ϵ_inf), α = $(interaction.α), r_c = $(interaction.r_c), k_c = $(interaction.k_c), L = $(interaction.L))")

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

Base.show(io::IO, interaction::IcmEwald2DInteraction) = print(io, "IcmEwald2DInteraction(n_atoms = $(interaction.n_atoms), ϵ = $(interaction.ϵ), α = $(interaction.α), r_c = $(interaction.r_c), k_c = $(interaction.k_c), γ = $(interaction.γ), L = $(interaction.L), N_image = $(interaction.N_image))")

struct IcmEwald3DInteraction{T, TP} <: ExTinyMD.AbstractInteraction
    n_atoms::Int
    ϵ::T
    α::T
    r_c::T
    k_c::T
    γ::Tuple{T, T}
    L::NTuple{3,T}
    N_image::Int # layer of the image charges
    N_pad::TP # layer of the z padding, L_z_pad = (2 * N_pad + 1) * L_z
    k_set::Vector{Tuple{T, T, T, T}}
end

Base.show(io::IO, interaction::IcmEwald3DInteraction) = print(io, "IcmEwald3DInteraction(n_atoms = $(interaction.n_atoms), ϵ = $(interaction.ϵ), α = $(interaction.α), r_c = $(interaction.r_c), k_c = $(interaction.k_c), γ = $(interaction.γ), L = $(interaction.L), N_image = $(interaction.N_image), N_pad = $(interaction.N_pad))")