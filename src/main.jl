include("physics.jl")

function toFile(Λ, filename)
    io = open(filename, "w")
    write(io, "x,y,z,vx,vy,vz\n")
    for i in 1:size(Λ, 1)
        for j in 1:2
            for k in 1:3 
                write(io, string(Λ[i, j, k], " "))
            end
        end
        write(io, "\n")
    end
    close(io)
end

const steps = 5000000
const hMin = 1e-7
const hMax = 5e-1
const errorMin = 1e-9
const errorMax = 1e-6
const nMax = 500

# init format = [x, y, z; vx, vy, vz]
# add more rows for more trajectories
init = cat([-15. 0. 11.; 1. 0. 0.], 
           [15. 0. 10.; -1. 0. 0.], 
           [10. 0. 4.; -1. 0. 0.], 
           dims = 3)

for i in 1:length(init[1, 1, :])
    @time Λ = RK45(steps, hMin, hMax, errorMin, 
                   errorMax, nMax, init[:, :, i])
    toFile(Λ, string("../data/trajectory", i, ".txt"))
end
