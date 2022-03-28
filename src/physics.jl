using LinearAlgebra

const r0 = 6.371 * 1e6 # average earth radius
const mp = 1.6726219 * 1e-27 # proton mass
const q = 1.60217662 * 1e-19 # proton charge
const μ0 = 4.0 * pi * 1e-7 # magnetic vacuum permeability
const m0 = 8.22 * 1e22 # earth magnetic moment
const C = q * μ0 * m0 * 1e-5 / (2.5 * r0^2 * mp * 4.0 * pi) # time scale factor
const θ = deg2rad(11.7) # earth magnetic moment polar angle
const ϕ = 0.0 # earth magnetic moment azimuth angle
const m = [cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)] # magnetic moment vector

# constants used in for RK4-5 method
const b1 = 1.0 / 4.0
const c1 = 3.0 / 32.0
const c2 = 9.0 / 32.0
const d1 = 1932.0 / 2197.0
const d2 = -7200.0 / 2197.0
const d3 = 7296.0 / 2197.0
const e1 = 439.0 / 216.0
const e2 = -8.0
const e3 = 3680.0 / 513.0
const e4 = -845.0 / 4104.0
const f1 = -8.0 / 27.0
const f2 = 2.0
const f3 = -3544.0 / 2565.0
const f4 = 1859.0 / 4104.0
const f5 = -11.0 / 40.0
const λ1 = 25.0 / 216.0
const λ3 = 1408.0 / 2565.0
const λ4 = 2197.0 / 4101.0
const λ5 = -1.0 / 5.0
const Λ1 = 16.0 / 135.0
const Λ3 = 6656.0 / 12825.0
const Λ4 = 28561.0 / 56430.0
const Λ5 = -9.0 / 50.0
const Λ6 = 2.0 / 55.0

function f(λ::Array{Float64, 2})::Array{Float64, 2}
    @views r = λ[1, :]
    @views v = λ[2, :]
    rLength = norm(r)
    rHat = r ./ rLength
    B = C ./ rLength^3 .* (3.0 .* dot(m, rHat) .* rHat .- m)
    return vcat(v', cross(v, B)')
end

function RK45(n::Int64, hMin::Float64, hMax::Float64,
              errorMin::Float64, errorMax::Float64,
              nMax::Int64, init::Array{Float64, 2})::Array{Float64, 3}
    h = hMax
    Λ = zeros(Float64, (n + 1, 2, 3))
    @views Λ[1, :, :] = init
    for i in 1:n
        for j in 1:nMax
            Λi = Λ[i, :, :]
            k1 = h .* f(Λi)
            k2 = h .* f(Λi .+ b1 .* k1)
            k3 = h .* f(Λi .+ c1 .* k1 .+ c2 .* k2)
            k4 = h .* f(Λi .+ d1 .* k1 .+ d2 .* k2 .+ d3 .* k3)
            k5 = h .* f(Λi .+ e1 .* k1 .+ e2 .* k2 .+ e3 .* k3 .+ e4 .* k4)
            k6 = h .* f(Λi .+ f1 .* k1 .+ f2 .* k2 .+ f3 .* k3 .+ f4 .* k4 .+ f5 .* k5)

            λ = Λi .+ λ1 .* k1 .+ λ3 .* k3 .+ λ4 .* k4 .+ λ5 .* k5
            @views Λ[i+1, :, :] .= Λi .+ Λ1 .* k1 .+ Λ3 .* k3 .+ Λ4 .* k4 .+ Λ5 .* k5 .+ Λ6 .* k6

            @views error = abs(norm(Λ[i+1, 1, :]) .- norm(λ[1, :]))
            @views error += abs(norm(Λ[i+1, 2, :]) .- norm(λ[2, :]))
            if error < errorMin
                h = min(2.0 * h, hMax)
            elseif error > errorMax
                h = max(h / 2.0, hMin)
            else
                break
            end
        end
    end
    return Λ
end
