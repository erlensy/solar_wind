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
h = 0.00001
initX = [-10.0 10.0]
initZ = [-10.0 -5.0 -2.5 0.0 2.5 5.0 10.0]
initVx = [-1.0 1.0]
tot = 0
for i in 1:length(initX)
    for z in initZ
        init = [initX[i] 0.0 z; initVx[i] 0.0 0.0]
        Λ = RK4(steps, h, init)
        toFile(Λ, string("../data/trajectory", tot, ".txt"))
        global tot += 1
    end
end
