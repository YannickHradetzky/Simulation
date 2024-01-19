import numpy as np
import matplotlib.pyplot as plt
import glob

NSteps = 10000
NStepsEqui = 100
NStepsSample = 10



def PlotAbsVel(N : int):
    filenames = glob.glob("out/vicsek/N" + str(N) + "/*.txt")
    filenames.sort(key = lambda x: float(x.split("Noise")[1].split("txt")[0][:-1]))
    # create array of eta values between 0 and 6.2 with 0.1 stepsize
    eta = np.arange(0, 6.3, 0.1)
    # create array of mean values
    mean = []
    std = []
    std = np.zeros(len(eta))
    for f in filenames:
        data = np.loadtxt(f)
        absvel = data[1:]
        mean.append(np.mean(absvel))
        std = np.std(absvel)
    # plot values with errorbars
    print(mean)
    print(std)
    plt.errorbar(eta, mean, yerr=std, fmt='x')
    plt.xlabel('eta')
    plt.ylabel('Average velocity')
    plt.grid()
    plt.show()
        

        


def PlotAllAbsVel():
    files = glob.glob("out/vicsek/1*.txt")
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # chose 4 different symbols
    symbols = ['o', 's', 'D', '^']
    for f in files:
        data = np.loadtxt(f)
        eta = data[:, 0]
        absvel = data[:, 1]
        label = f.split('.'[0])
        ax.plot(eta, absvel, label=label, marker=symbols.pop())
    ax.set_xlabel('eta')
    ax.set_ylabel('Average velocity')
    ax.grid()
    ax.legend()
    plt.show()

PlotAbsVel(1000)

