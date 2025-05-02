# refer to "Accurate Error Estimates and Optimal Parameter Selection in Ewald Summation for Dielectrically Confined Coulomb Systems" https://arxiv.org/abs/2503.18126

function icm_energy_error(Lx, Ly, H, M, ϵ, Cq, gu, gd)
    return (8 * π^2 / (Lx * Ly * ϵ)) * (Cq * (abs(gu * gd))^(floor((M + 1) / 2)) * exp(- (4 * π * H * floor((M + 1) / 2)) / max(Lx, Ly))) / (1 - abs(gu * gd) * exp(- (4 * π * H) / max(Lx, Ly)))
end

function elc_energy_error(Lx, Ly, Lz, H, ϵ, Cq)
    t = 4 * π^2 * Cq / (ϵ * Lx * Ly * (1 - exp(- (2 * π * Lz) / max(Lx, Ly))))
    return t * exp(- (2 * π * (Lz - H)) / max(Lx, Ly)) / (Lz - H)
end

function icm_elc_energy_error(Lx, Ly, Lz, H, M, ϵ, Cq, gu, gd)
    t = 4 * π^2 * Cq / (ϵ * Lx * Ly * (1 - exp(- (2 * π * Lz) / max(Lx, Ly))))
    e = t * exp(- (2 * π * (Lz - H)) / max(Lx, Ly)) / abs(Lz - H)
    for l in 1:M
        Cl = abs(gu^(ceil(l / 2)) * gd^(floor(l / 2))) + abs(gu^(floor(l / 2)) * gd^(ceil(l / 2)))
        e += t * Cl * exp(- (2 * π * (Lz - (l + 1) * H)) / max(Lx, Ly)) / abs(2 * (Lz - (l + 1) * H))
    end
    e += t * 4 * exp(- (2 * π * Lz) / max(Lx, Ly)) / (Lz * (1 - exp(- (4 * π * H) / max(Lx, Ly))))

    return e
end