struct Sys3D{T}
    L::NTuple{3, T} # (Lx, Ly, Lz)
    N::NTuple{3, Int} # (Nx, Ny, Nz)
    ϵ::T
end

Sys3D(L::NTuple{3, T}, N::NTuple{3, Int}; ϵ::T = one(T)) where{T} = Sys3D{T}(L, N, ϵ)

function Energy_3D(sys::Sys3D, position::Vector{Point{3, T}}, charge::Vector{T}) where T<:Number
    @warn "This will be very slow!"
    energy = zero(T)
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(charge)
            q_j = charge[j]
            pos_j = position[j]
            for mx in -sys.N[1]:sys.N[1]
                for my in -sys.N[2]:sys.N[2]
                    for mz in -sys.N[3]:sys.N[3]
                        energy += CoulumbEnergy(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], mz * sys.L[3])) / T(2) / sys.ϵ
                    end
                end
            end
        end
    end
    return energy
end

function Force_3D(sys::Sys3D, position::Vector{Point{3, T}}, charge::Vector{T}) where T<:Number
    @warn "This will be very slow!"
    force = [Point(zero(T), zero(T), zero(T)) for _=1:length(charge)]
    for i in 1:length(charge)
        q_i = charge[i]
        pos_i = position[i]
        for j in 1:length(charge)
            q_j = charge[j]
            pos_j = position[j]
            for mx in -sys.N[1]:sys.N[1]
                for my in -sys.N[2]:sys.N[2]
                    for mz in -sys.N[3]:sys.N[3]
                        force[i] += CoulumbForce(q_i, q_j, pos_i, pos_j + Point(mx * sys.L[1], my * sys.L[2], mz * sys.L[3])) / sys.ϵ
                    end
                end
            end
        end
    end
    return force 
end
