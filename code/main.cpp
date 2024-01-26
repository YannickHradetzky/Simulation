#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <mutex>
#include <algorithm> // for std::min

#include "Simulation.h"
#include "Potentials.h"

// Settings
int NSim = 1000;
int NEquil = 100;
int NSamp = 10;
double DensityInit = 4;
double NoiseForDensity = 0.5;
double TimeStepSize = 1;
int times = 10;
int nMax = 10000;

using namespace std;

int main() {
    mutex globalMutex; // Mutex to synchronize access to shared variables

    int num_cores = thread::hardware_concurrency() - 2;
    cout << "Number of CPU cores: " << num_cores << endl;

    vector<thread> threads;

    for (int n = 10; n <= nMax; n *= 10) {
        // create a folder for each value of n
        system(("mkdir -p out/vicsek/" + to_string(n)).c_str());
        cout << "Running for n = " << n << endl;
        
        for (double Noise = 0; Noise <= 2 * M_PI; Noise += 0.1) {
            while (true) {
                if (threads.size() <= num_cores) {
                    threads.emplace_back([n, Noise, &globalMutex]() {
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

                        string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Noise_" + to_string(Noise) + ".txt";

                        // Open the file
                        ofstream ResultsFile(filenameNoise);

                        // For each Noise value, run the simulation multiple times
                        for (int i = 1; i <= times; ++i) {
                            // Use unique seed for each thread
                            srand(i * static_cast<unsigned int>(hash<thread::id>{}(this_thread::get_id())));

                            // Run the simulation for a given noise
                            vector<double> Results = sim.RunVicsekForNoise(Noise);

                            // Use the global mutex to synchronize file writing
                            lock_guard<mutex> lock(globalMutex);

                            // Write the results to a file
                            for (int j = 0; j < Results.size(); ++j) {
                                ResultsFile << Results[j] << " ";
                            }
                            ResultsFile << endl;
                        }
                        ResultsFile.close();
                        cout << "Finished Noise = " << Noise << endl;
                    });
                    break;
                }
                else {
                    // Wait for any thread to finish, then join and erase it
                    for (auto it = threads.begin(); it != threads.end(); ++it) {
                        if (it->joinable()) {
                            it->join();
                            threads.erase(it);
                            break;
                        }
                    }
                }
            }
        }
    }

    // for(int n = 10; n <= nMax; n *= 10 ){
    //     cout << "Running for n = " << n << endl;
        
    //     for (double Density = 0.1; Density <= 10; Density += 0.1){
    //         while (true) {
    //             if (threads.size() <= num_cores) {
    //                 threads.emplace_back([n, Density, &globalMutex]() {
    //                     // Create a new Simulation object for each thread
    //                     Simulation sim;
    //                     sim.NParticles = n;
    //                     sim.DensityInit = Density;
    //                     sim.NStepsSimulation = NSim;
    //                     sim.NStepsEquilibration = NEquil;
    //                     sim.NStepsSampling = NSamp;
    //                     sim.TimeStepSize = TimeStepSize;
    //                     sim.PrintInitInfo = false;
    //                     sim.PrintProgress = false;

    //                     string filenameDensity = "out/vicsek/" + to_string(sim.NParticles) + "/Density_" + to_string(Density) + ".txt";

    //                     // Open the file
    //                     ofstream ResultsFile(filenameDensity);

    //                     // For each Density value, run the simulation multiple times
    //                     for (int i = 1; i <= times; ++i) {
    //                         // Use unique seed for each thread
    //                         srand(i * static_cast<unsigned int>(hash<thread::id>{}(this_thread::get_id())));

    //                         // Run the simulation for a given density
    //                         vector<double> Results = sim.RunVicsekForDensity(Density, NoiseForDensity);

    //                         // Use the global mutex to synchronize file writing
    //                         lock_guard<mutex> lock(globalMutex);

    //                         // Write the results to a file
    //                         for (int j = 0; j < Results.size(); ++j) {
    //                             ResultsFile << Results[j] << " ";
    //                         }
    //                         ResultsFile << endl;
    //                     }
    //                     ResultsFile.close();
    //                     cout << "Finished Density = " << Density << endl;
    //                 });
    //                 break;
    //             }
    //             else {
    //                 // Wait for any thread to finish, then join and erase it
    //                 for (auto it = threads.begin(); it != threads.end(); ++it) {
    //                     if (it->joinable()) {
    //                         it->join();
    //                         threads.erase(it);
    //                         break;
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // Join remaining threads
    for (auto &t : threads) {
        if (t.joinable()) {
            t.join();
        }
    }

    return 0;
}
