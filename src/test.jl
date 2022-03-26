using LinearAlgebra, Plots

const r0 = 6.371 * 1e6 # average earth radius
const mp = 1.6726219 * 1e-27 # proton mass
const q = 1.60217662 * 1e-19 # proton charge
const μ0 = 4.0 * pi * 1e-7 # magnetic vacuum permeability
const m0 = 8.22 * 1e22 # earth magnetic moment
const C = q * μ0 * m0 * 1e-5 / (2.5 * r0^2 * mp * 4.0 * pi) # time scale factor
const θ = deg2rad(11.7) # earth magnetic moment polar angle
const ϕ = 0.0 # earth magnetic moment azimuth angle
const m = [cos(ϕ) * sin(θ), sin(ϕ) * sin(θ), cos(θ)] # magnetic moment vector

function f(λ::Array{Float64, 2})::Array{Float64, 2}
    @views r = λ[1, :]; v = λ[2, :]
    rLength = norm(r)
    rHat = r ./ rLength
    B = C ./ rLength^3 .* (3.0 * dot(m, rHat) .* rHat .- m)
    return vcat(v', cross(v, B)')
end

function RK4(n::Int64, h::Float64, init::Array{Float64, 2})::Array{Float64, 3}
    Λ = zeros(Float64, (n + 1, 2, 3))
    @views Λ[1, :, :] = init
    for i in 1:n
        Λi = Λ[i, :, :]
        k1 = f(Λi)
        k2 = f(Λi .+ 0.5 .* h .* k1)
        k3 = f(Λi .+ 0.5 .* h .* k2)
        k4 = f(Λi .+ h .* k3)
        k = k1 .+ 2.0 .* k2 .+ 2.0 .* k3 .+ k4
        Λ[i+1, :, :] = Λi .+ h ./ 6.0 .* k 
    end
    return Λ
end

function toFile(Λ, filename)
    io = open(filename, "w")
    write(io, "x,y,z,vx,vy,vz\n")
    for i in 1:size(Λ, 1)
        for j in 1:2
            for k in 1:3 
                write(io, string(Λ[i, j, k], ","))
            end
        end
        write(io, "\n")
    end
    close(io)
end

init = [10.0 0.0 5.0; -1.0 0.0 0.0]
steps = 100000
h = 0.0005
Λ = RK4(steps, h, init)
toFile(Λ, "hei.txt")
plot(xlim=(-15.0, 15.0), ylim=(-15.0, 15.0))
plot!([0], [0], seriestype = :scatter, markersize = 10.0)
plot!(Λ[:, 1, 1], Λ[:, 1, 3])
