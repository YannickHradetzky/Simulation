#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <mutex>

#include "Simulation.h"
#include "Potentials.h"

// Compile with: g++ -std=c++11 -O3 code/main.cpp code/Simulation.cpp -o main

// Settings
int NSim = 1000;
int NEquil = 100;
int NSamp = 10;
double DensityInit = 4;
double NoiseForDensity = 0.5;
double TimeStepSize = 1;
double VelInit = 0.3; // Initial velocity of particles
int times = 10;



using namespace std;


int main() {
    vector<thread> threads; // Store threads in a vector
    system("mkdir -p out/vicsek"); // Create output directory (if it doesn't exist)
    system("rm -rf out/vicsek/*"); // Remove previous results (if any)

    mutex globalMutex; // Mutex to synchronize access to shared variables

    for (int n = 10; n <= 1000; n *= 10) {
        // create a folder for each value of n
        system(("mkdir -p out/vicsek/" + to_string(n)).c_str());

        cout << "Running for n = " << n << endl;
        vector<thread> nThreads; // Store threads for a specific number of particles
        // Parallelize the loop for different values of Noise
        for (double Noise = 0; Noise <= 2 * M_PI; Noise += 0.1) {
            // create a thread for each value of Noise
            nThreads.emplace_back([n, Noise, &globalMutex]() {
                // Create a new Simulation object for each thread
                Simulation sim;
                sim.NParticles = n;
                sim.DensityInit = DensityInit;
                sim.NStepsSimulation = NSim;
                sim.NStepsEquilibration = NEquil;
                sim.NStepsSampling = NSamp;
                sim.TimeStepSize = TimeStepSize;
                sim.PrintInitInfo = false;
                sim.PrintProgress = false;

                string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Corr" + to_string(Noise) + ".txt";

                // Open the file in
                ofstream ResultsFile(filenameNoise);

                // For each Noise value, run the simulation multiple times
                for (int i = 1; i <= times; ++i) {
                    // Use unique seed for each thread
                    srand(i * static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id())));

                    // Initialize the system
                    sim.InitVicsek(VelInit);

                    // Equilibrate the system
                    for(int i = 0; i < sim.NStepsEquilibration; ++i) {
                        sim.VicsekUpdate(Noise);
                    }

                    vector<double> Temp;
                    vector<double> Results(200, 0);
                    
                    // Main simulation loop
                    for(int j = 0; j < sim.NStepsSimulation; ++j){
                        sim.VicsekUpdate(Noise);
                        if (j % sim.NStepsSampling == 0){
                            // Sample the system
                            // 1) Calculate Average Velocity
                            double AvgX = 0; 
                            double AvgY = 0;
                            for (const auto &p : sim.Particles) {
                                AvgX += p->vx;
                                AvgY += p->vy;
                            }
                            AvgX /= sim.NParticles;
                            AvgY /= sim.NParticles;
                            double AvgVel = sqrt(AvgX * AvgX + AvgY * AvgY);
                            // 2) Calculate Correlation
                            Temp = sim.ComputeOrientationCorrelationVicsek(AvgX, AvgY);

                            // 3) Add to Results
                            for(int k = 0; k < Temp.size(); ++k){
                                Results[k] += Temp[k];
                            }
                        }
                    // Normalize and Print Results
                    for(int k = 0; k < Results.size(); ++k){
                        Results[k] /= (sim.NStepsSimulation / sim.NStepsSampling);
                        ResultsFile << Results[k] << " ";
                    }
                    ResultsFile << endl;   
                    }
                }
                ResultsFile.close();
                cout << "Finished Noise = " << Noise << endl;
            });
        }
        // Join all threads for a specific number of particles
        for (auto &t : nThreads) {
            t.join();
        }
    }

    return 0;

}




