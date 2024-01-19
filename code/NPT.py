import matplotlib.pyplot as plt
import pandas as pd
import glob

N = 1000

def PlotAllPressures(NParticels : int):
    n = str(NParticels)
    files = glob.glob("out/NPT/" + "*" + n + ".txt")
    # Sort...start with lowest pressure
    sortedFiles = sorted(files, key=lambda x: float(x.split('_of_P_')[1].split('_forN')[0]))
    print(sortedFiles)

    df = pd.read_csv(sortedFiles[0], delimiter=" ")
    print(df.columns)

    for file in sortedFiles:
        df = pd.read_csv(file, skipfooter=2, delimiter=" ", engine="python")
        steps = df["Step"].values
        density = df["Density"].values
        volume = df["Volume"].values
        pressure = str(df["Pressure"].values[0])

        plt.plot(steps, volume, label=pressure)


    plt.tight_layout()
    plt.legend(loc = "upper right", bbox_to_anchor = (1.20, 1))
    plt.show()


PlotAllPressures(N)




