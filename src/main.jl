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

steps = 500000
hMin = 1e-7
hMax = 5e-1
errorMin = 1e-9
errorMax = 1e-6
nMax = 500

initX = [20.0 -20.0]
initZ = [10.0 0.0]
initVx = [-1.0 1.0]
tot = 0
for i in 1:length(initX)
    for z in initZ
        init = [initX[i] 0.0 z; initVx[i] 0.0 0.0]
        @time Λ = RK45(steps, hMin, hMax, errorMin, 
                       errorMax, nMax, init)
        toFile(Λ, string("../data/trajectory", tot, ".txt"))
        global tot += 1
    end
end
