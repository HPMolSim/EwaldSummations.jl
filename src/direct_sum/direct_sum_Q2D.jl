struct SysQ2D{T}
    γ::NTuple{2, T} # (γ_up, γ_down)
    L::NTuple{3, T} # (Lx, Ly, Lz)
    N_real::Int
    N_img::Int
    ϵ::T
end

SysQ2D(γ::NTuple{2, T}, L::NTuple{3, T}, N_real::Int, N_img::Int; ϵ::T = one(T)) where{T} = SysQ2D{T}(γ, L, N_real, N_img, ϵ)

function ICM_reflect(γ::NTuple{2, T}, L::NTuple{3, T}, N_img::Int, position::Vector{Point{3, T}}, charge::Vector{T}) where{T}
    γ_up = γ[1]
    γ_down = γ[2]
    reflect_position = Vector{Point{3, T}}()
    reflect_charge = Vector{T}()

    for i in 1:length(charge)
        x, y, z = position[i]
        q = charge[i]
        push!(reflect_position, Point(x, y, z))
        push!(reflect_charge, q)
    end

    for i = 1:length(charge)
        x, y, z = position[i]
        q = charge[i]
        position_up = Vector{T}()
        position_down = Vector{T}()
        charge_up = Vector{T}()
        charge_down = Vector{T}()
        push!(position_up, 2.0 * L[3] - z)
        push!(position_down, - z)
        push!(charge_up, γ_up * q)
        push!(charge_down, γ_down * q)
        for m in 2:N_img
            push!(position_up, 2 * L[3] - position_down[m-1])
            push!(position_down, - position_up[m-1])
            push!(charge_up, γ_up * charge_down[m-1])
            push!(charge_down, γ_down * charge_up[m-1])
        end
        for m in 1:N_img
            push!(reflect_position, Point(x, y, position_up[m]))
            push!(reflect_charge, charge_up[m])
            push!(reflect_position, Point(x, y, position_down[m]))
            push!(reflect_charge, charge_down[m])
        end
    end

    return reflect_position, reflect_charge
end

function SysQ2DInit(sys::SysQ2D, position::Vector{Point{3, T}}, charge::Vector{T}) where T <: Number

    return ICM_reflect(sys.γ, sys.L, sys.N_img, position, charge)
end

function Energy_Q2D(sys::SysQ2D, position::Vector{Point{3, T}}, charge::Vector{T}, reflect_position::Vector{Point{3, T}}, reflect_charge::Vector{T}) where T<:Number
    energy = zero(T)
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(reflect_charge)
            q_j = reflect_charge[j]
            pos_j = reflect_position[j]
            for mx in -sys.N_real:sys.N_real
                for my in -sys.N_real:sys.N_real
                    energy += CoulumbEnergy(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], zero(T))) / T(2) / sys.ϵ
                end
            end
        end
    end
    return energy
end

function Force_Q2D(sys::SysQ2D, position::Vector{Point{3, T}}, charge::Vector{T}, reflect_position::Vector{Point{3, T}}, reflect_charge::Vector{T}) where T<:Number
    force = [Point(zero(T), zero(T), zero(T)) for _=1:length(charge)]
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(reflect_charge)
            q_j = reflect_charge[j]
            pos_j = reflect_position[j]
            for mx in -sys.N_real:sys.N_real
                for my in -sys.N_real:sys.N_real
                    force[i] += CoulumbForce(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], 0.0)) / sys.ϵ
                end
            end
        end
    end
    return force 
end