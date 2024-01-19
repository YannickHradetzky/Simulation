#include <iostream>

#include "Simulation.h"
#include "Potentials.h"

// Compile with: g++ -std=c++11 -O3  code/main.cpp code/Simulation.cpp  -o main


using namespace std;
int main() {

    srand(12);
    Simulation sim;

    sim.NStepsSimulation = 1000;
    sim.NStepsEquilibration = 100;
    sim.NStepsSampling = 10;
    sim.TimeStepSize = 1;

    ofstream SettingsFile;
    SettingsFile.open("out/vicsek/Settings.txt");
    SettingsFile << "NStepsSimulation = " << sim.NStepsSimulation << endl;
    SettingsFile << "NStepsEquilibration = " << sim.NStepsEquilibration << endl;
    SettingsFile << "NStepsSampling = " << sim.NStepsSampling << endl;
    SettingsFile.close();

    sim.PrintInitInfo = false;
    sim.NParticles = 1000;
    int times = 10;
    double Noise;
    ofstream ResultsFile;
    system("rm -rf out/vicsek/*");
    for(int n = 10; n < 100000; n*=10)
    {
        system(("mkdir -p out/vicsek/" + to_string(n)).c_str());
        sim.NParticles = n;
        for(Noise = 0; Noise <= 2 * M_PI; Noise += 0.1)
        {
            string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Noise_" + to_string(Noise) + ".txt";
            ResultsFile.open(filenameNoise);
            for(int i = 1 ; i <= times; ++i)
            {
                // Seed the random number generator
                srand(i);
                // Run the simulation for a given noise
                vector<double> Results = sim.RunVicsekForNoise(Noise);
                // Write the results to a file
                for(int j = 0; j < Results.size(); ++j)
                {
                    ResultsFile << Results[j] << " ";
                }
                ResultsFile << endl;
            }
            ResultsFile.close();
        }
        ResultsFile.close();
    }
        


    return 0;
}

