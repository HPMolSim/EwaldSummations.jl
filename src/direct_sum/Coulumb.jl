function CoulumbEnergy(q_i::T, q_j::T, coo_i::Point{3, T}, coo_j::Point{3, T}) where T<:Number
    r = dist2(coo_i, coo_j)
    if iszero(r) == false
        return q_i * q_j / (4π * sqrt(r))
    else
        return 0.0
    end
end

function CoulumbForce(q_i::T, q_j::T, coo_i::Point{3, T}, coo_j::Point{3, T}) where T<:Number
    r = dist2(coo_i, coo_j)
    if iszero(r) == false
        rho = sqrt(r)
        F = q_i * q_j / (4π * r)
        angle = (coo_i - coo_j) * (1/rho)
        return F * angle
    else
        return Point(0.0, 0.0, 0.0)
    end
end