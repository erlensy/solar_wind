from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import numpy as np

# all constants are in SI units
earthRadius = 6371.0 * 1e3
protonMass = 1.6726219 * 1e-27
electronMass = 9.10938356 * 1e-31
electronCharge = -1.60217662 * 1e-19
vacuumPermeability = 4.0 * np.pi * 1e-7

earthMagneticMomentTilt = 11.7 * np.pi / 180.0
earthMagneticMoment = 8.22 * 1e22 * np.array([0.0, np.sin(earthMagneticMomentTilt), np.cos(earthMagneticMomentTilt)])

# charged particle
protonMass = 1.6726219 * 1e-27
mass = protonMass
charge = -1.0 * electronCharge

def plot_magneticField():
    scale = 4
    X, Y, Z = np.meshgrid(
            np.arange(-scale * earthRadius, scale * earthRadius, earthRadius),
            np.arange(-scale * earthRadius, scale * earthRadius, earthRadius),
            np.arange(-scale * earthRadius, scale * earthRadius, earthRadius))

    U, V, W = np.zeros(X.shape), np.zeros(Y.shape), np.zeros(Z.shape)

    for i in range(len(X)):
        for j in range(len(X[0])):
            for k in range(len(X[0, 0])):
                r = np.array([X[i, j, k], Y[i, j, k], Z[i, j, k]])
                if np.linalg.norm(r) >= 1.5 * earthRadius:
                    U[i, j, k], V[i, j, k], W[i, j, k] = magneticField(r)


    ax = plt.figure().add_subplot(projection='3d')
    ax.quiver(X / earthRadius, Y / earthRadius, Z / earthRadius, U, V, W, length = 1, normalize=True)

    # draw sphere
    phi, theta = np.mgrid[0:2*np.pi:20j, 
                          0:np.pi:10j]
    x = np.cos(phi) * np.sin(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(theta)
    ax.plot_wireframe(x, y, z, color="black")

    ax.set_xlim(-scale, scale)
    ax.set_ylim(-scale, scale)
    ax.set_zlim(-scale, scale)
    plt.show()

def magneticField(r):
    rLength = np.linalg.norm(r)
    return vacuumPermeability / (4. * np.pi) * (3. * r * np.dot(earthMagneticMoment, r) / (rLength**5) - earthMagneticMoment / (rLength**3))

def f(q, m, z):
    """z = [x, v], returns d/dt of z"""
    return np.array([z[1, :], q / m * np.cross(z[1, :], magneticField(z[0, :]))])

def RK4(n, h, q, m, r0, v0):
    z = np.zeros((n, 2, 3))
    z[0, 0, :] = r0
    z[0, 1, :] = v0

    for i in range(n-1):
        k1 = f(q, m, z[i, :, :])
        k2 = f(q, m, z[i, :, :] + 0.5 * h * k1)
        k3 = f(q, m, z[i, :, :] + 0.5 * h * k2)
        k4 = f(q, m, z[i, :, :] + h * k3)

        z[i+1, :, :] = z[i, :, :] + h / 6. * (k1 + 2.*k2 + 2.*k3 + k4)

    return z

if __name__ == "__main__":
    #plot_magneticField()
    r0 = np.array([2.*earthRadius, 0, 0])

    v0 = np.array([489. * 1e3, 489. * 1e3, 0.])
    h = 0.005
    z = RK4(10000, h, charge, mass, r0, v0)
    ax = plt.figure().add_subplot(projection='3d')
    ax.plot(z[:, 0, 0], z[:, 0, 1], z[:, 0, 2])

    # earth plot
    phi, theta = np.mgrid[0:2*np.pi:80j, 0:np.pi:20j]
    x = earthRadius * np.cos(phi) * np.sin(theta)
    y = earthRadius * np.sin(phi) * np.sin(theta)
    z = earthRadius * np.cos(theta)
    #ax.plot_wireframe(x, y, z, color="black")
    plt.show()
