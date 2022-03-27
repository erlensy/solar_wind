import numpy as np
from matplotlib import pyplot as plt

def plotMagField(XZ = False, XY = False):
    x = y = np.linspace(-15, 15, 100)
    theta = np.deg2rad(11.7)
    phi = 0.0
    m = np.array([np.cos(phi) * np.sin(theta), np.sin(phi) * np.sin(theta), np.cos(theta)])
    X, Y = np.meshgrid(x, y)
    U, V = np.zeros(X.shape), np.zeros(Y.shape)
    for i in range(len(X)):
        for j in range(len(X[0])):
            if XZ: r = np.array([X[i, j], 0.0, Y[i, j]])
            elif XY: r = np.array([X[i, j], Y[i, j], 0.0])
            rLength = np.linalg.norm(r)
            rHat = r / rLength
            if rLength >= 1.0:
                if XZ: U[i, j], _, V[i, j] = -(3.0 * np.dot(m, rHat) * rHat - m)
                elif XY: U[i, j], V[i, j], _ = -(3.0 * np.dot(m, rHat) * rHat - m)
    plt.streamplot(X, Y, U, V, color = "lightblue", density = 0.9) 

    scale = 5.0
    if XZ: plt.plot([-scale * m[0], scale * m[0]], [-scale * m[2], scale * m[2]], color = "blue", linestyle = "--", alpha = 0.8)
    elif XY: 
        scale = 20.0
        plt.plot([-scale * m[0], scale * m[0]], [-scale * m[1], scale * m[1]], color = "blue", linestyle = "--", alpha = 0.8)

def plotEarth():
    im = plt.imread("../data/earth.jpeg")
    plt.imshow(im, extent = [-1.1, 1.1, -1.1, 1.1], alpha = 1.0)

def plotTrajectories(n, XZ = False, XY = False):
    plt.xlabel("x/r$_0$"); plt.ylabel("z/r$_0$")
    for i in range(0, n):
        data = np.loadtxt(f"../data/trajectory{i}.txt", skiprows = 1)
        if XZ: yData = data[:, 2]
        if XY: yData = data[:, 1]
        plt.plot(data[:, 0], yData, color = "red", linewidth = 0.6)

def setupFig(XZ = False, XY = False):
    fig = plt.figure()
    plt.xlim(-15, 15); plt.ylim(-15, 15)
    if XZ: 
        plt.xlabel("x/r$_0$") 
        plt.ylabel("z/r$_0$")
    if XY: 
        plt.xlabel("x/r$_0$") 
        plt.ylabel("y/r$_0$")

def plot1():
    setupFig(XZ = True)
    plotEarth()
    plotMagField(XZ = True)
    plt.savefig("../figures/magFieldXZ.pdf", dpi = 600)
    plt.close()

def plot2():
    setupFig(XY = True)
    plotEarth()
    plotMagField(XY = True)
    plt.savefig("../figures/magFieldXY.pdf", dpi = 600)
    plt.close()

def plot3():
    setupFig(XZ = True)
    plotEarth()
    plotMagField(XZ = True)
    plotTrajectories(14, XZ = True)
    plt.savefig("../figures/trajectoriesXZ.pdf", dpi = 600)
    plt.close()

def plot4():
    setupFig(XY = True)
    plotEarth()
    plotMagField(XY = True)
    plotTrajectories(14, XY = True)
    plt.savefig("../figures/trajectoriesXY.pdf", dpi = 600)
    plt.close()

if __name__ == "__main__":
    plot1()
    plot2()
    plot3()
    plot4()
