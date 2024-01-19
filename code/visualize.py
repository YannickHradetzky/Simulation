import glob

import matplotlib.pyplot as plt
import numpy as np
import sys


def PlotParticlePositions(show: int, filename: str):
    # create interactive 3D plot
    data = np.loadtxt(filename)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    x = data[:, 0]
    y = data[:, 1]
    z = data[:, 2]
    ax.scatter(x, y, z, c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    # Set title
    ax.set_title('Particle positions')
    plt.savefig('out/positions.png')
    if show == 1:
        plt.show()


def PlotParticlePositions2D(show: int, filename: str):
    # create 2D plot of particle positions
    data = np.loadtxt(filename)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = data[:, 0]
    y = data[:, 1]
    ax.scatter(x, y, c='r', marker='o')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    # Set title
    ax.set_title('Particle positions')
    plt.savefig('out/positions.png')
    if show == 1:
        plt.show()


def PlotParticleVelocityDistribution(show: int, filename: str):
    # Plot Particle velocities distribution as histogram using black color
    data = np.loadtxt(filename)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.hist(data, bins=100, color='black', histtype='step')
    ax.set_xlabel('Velocity')
    ax.set_ylabel('Number of particles')
    ax.set_title('Particle velocities distribution')
    plt.savefig('out/velocity_distribution.png')

    if show == 1:
        plt.show()


def PlotRDF(show: int, filename: str):
    # Plot Particle velocities distribution as histogram using black color
    data = np.loadtxt(filename)
    x = data[:, 0]
    y = data[:, 1]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x, y, color='black')
    ax.set_xlabel('r')
    ax.set_ylabel('g(r)')
    # remove .txt and folder structure from filename
    filename = filename[:-4]
    filename = filename.split('/')[-1]
    ax.set_title("RDF for " + filename)

    plt.savefig('out/images/' + filename + '.png')

    if show == 1:
        plt.show()


if __name__ == '__main__':
    FunctionName = sys.argv[1]
    arg1 = int(sys.argv[2])
    arg2 = str(sys.argv[3])

    print(FunctionName)
    print(arg1)
    print(arg2)

    if FunctionName == "PlotParticlePositions":
        PlotParticlePositions(arg1, arg2)
    if FunctionName == "PlotParticleVelocityDistribution":
        PlotParticleVelocityDistribution(arg1, arg2)
    if FunctionName == "PlotParticlePositions2D":
        PlotParticlePositions2D(arg1, arg2)
    if FunctionName == "PlotRDF":
        PlotRDF(arg1, arg2)
