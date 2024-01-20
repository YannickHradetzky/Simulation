﻿#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include <cmath>
#include <string>
#include <mutex>

#include "Simulation.h"
#include "Potentials.h"

using namespace std;

int main() {
    int times = 10;
    double Noise;
    double Density;

    // Outer loop parameters
    for (int n = 10; n < 10000; n *= 10) {
        system(("mkdir -p out/vicsek/" + to_string(n)).c_str());

        vector<thread> threads;

        // Parallelize the loop for different values of Noise
        for (Noise = 0; Noise <= 2 * M_PI; Noise += 0.1) {
            threads.emplace_back([n, times, Noise]() {
                // Create a new Simulation object for each thread
                Simulation sim;
                sim.NParticles = n;
                sim.PrintInitInfo = false;
                sim.PrintProgress = false;
                sim.NStepsSampling = 10000;
                sim.NStepsEquilibration = 1000;
                sim.NStepsSampling = 100;

                string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Noise_" + to_string(Noise) + ".txt";
                ofstream ResultsFile(filenameNoise);

                // Use unique seed for each thread
                srand(static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id())));

                for (int i = 1; i <= times; ++i) {
                    // Run the simulation for a given noise
                    vector<double> Results = sim.RunVicsekForNoise(Noise);

                    // Use a lock to synchronize file writing
                    static mutex fileMutex; // Static ensures it's shared among all threads
                    lock_guard<mutex> lock(fileMutex);

                    // Write the results to a file
                    for (int j = 0; j < Results.size(); ++j) {
                        ResultsFile << Results[j] << " ";
                    }
                    ResultsFile << endl;
                }

                // Close the file after the loop
                ResultsFile.close();
            });
        }

        // Join threads after processing the Noise loop
        for (auto& thread : threads) {
            thread.join();
        }

        // Clear threads vector for reuse
        threads.clear();

        // Parallelize the loop for different values of Density
        for (Density = 0.1; Density <= 10; Density += 0.1) {
            threads.emplace_back([n, times, Density]() {
                // Create a new Simulation object for each thread
                Simulation sim;
                sim.NParticles = n;
                sim.PrintInitInfo = false;
                sim.PrintProgress = false;
                sim.NStepsSampling = 10000;
                sim.NStepsEquilibration = 1000;
                sim.NStepsSampling = 100;

                string filenameDensity = "out/vicsek/" + to_string(sim.NParticles) + "/Density_" + to_string(Density) + ".txt";
                ofstream ResultsFile(filenameDensity);

                // Use unique seed for each thread
                srand(static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id())));

                for (int i = 1; i <= times; ++i) {
                    // Run the simulation for a given density
                    vector<double> Results = sim.RunVicsekForDensity(Density);

                    // Use a lock to synchronize file writing
                    static mutex fileMutex; // Static ensures it's shared among all threads
                    lock_guard<mutex> lock(fileMutex);

                    // Write the results to a file
                    for (int j = 0; j < Results.size(); ++j) {
                        ResultsFile << Results[j] << " ";
                    }
                    ResultsFile << endl;
                }

                // Close the file after the loop
                ResultsFile.close();
            });
        }

        // Join threads after processing the Density loop
        for (auto& thread : threads) {
            thread.join();
        }
    }

    return 0;
}










// #include <iostream>

// #include "Simulation.h"
// #include "Potentials.h"

// // Compile with: g++ -std=c++11 -O3  code/main.cpp code/Simulation.cpp  -o main


// using namespace std;
// int main() {
//     // Create a simulation object
//     Simulation sim;


//     sim.NStepsSimulation = 10000;
//     sim.NStepsEquilibration = 1000;
//     sim.NStepsSampling = 100;
//     sim.TimeStepSize = 1;
//     sim.PrintInitInfo = false;


//     ofstream SettingsFile;
//     SettingsFile.open("out/vicsek/Settings.txt");
//     SettingsFile << "NStepsSimulation = " << sim.NStepsSimulation << endl;
//     SettingsFile << "NStepsEquilibration = " << sim.NStepsEquilibration << endl;
//     SettingsFile << "NStepsSampling = " << sim.NStepsSampling << endl;
//     SettingsFile.close();


//     int times = 10;
//     double Noise;
//     double Density;
//     ofstream ResultsFile;

//     for(int n = 10; n < 10000; n*=10)
//     {
//         system(("mkdir -p out/vicsek/" + to_string(n)).c_str());
//         sim.NParticles = n;
//         for(Noise = 0; Noise <= 2 * M_PI; Noise += 0.1)
//         {
//             string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Noise_" + to_string(Noise) + ".txt";
//             ResultsFile.open(filenameNoise);
//             for(int i = 1 ; i <= times; ++i)
//             {
//                 // Seed the random number generator
//                 srand(i);
//                 // Run the simulation for a given noise
//                 vector<double> Results = sim.RunVicsekForNoise(Noise);
//                 // Write the results to a file
//                 for(int j = 0; j < Results.size(); ++j)
//                 {
//                     ResultsFile << Results[j] << " ";
//                 }
//                 ResultsFile << endl;
//             }
//             ResultsFile.close();
//         }



//         for(Density = 0.1; Density <= 10; Density += 0.1)
//         {
//             string filenameDensity = "out/vicsek/" + to_string(sim.NParticles) + "/Density_" + to_string(Density) + ".txt";
//             ResultsFile.open(filenameDensity);
//             for(int i = 1 ; i <= times; ++i)
//             {
//                 // Seed the random number generator
//                 srand(i);
//                 // Run the simulation for a given density
//                 vector<double> Results = sim.RunVicsekForDensity(Density);
//                 // Write the results to a file
//                 for(int j = 0; j < Results.size(); ++j)
//                 {
//                     ResultsFile << Results[j] << " ";
//                 }
//                 ResultsFile << endl;
//             }
//             ResultsFile.close();
//         }
//     }
        


//     return 0;
// }





// #include <iostream>
// #include <thread>
// #include <vector>
// #include <fstream>
// #include <cmath>
// #include <string>
// #include <mutex>

// #include "Simulation.h"
// #include "Potentials.h"

// using namespace std;

// int main() {
//     int times = 10;
//     double Noise;
//     double Density;
//     ofstream ResultsFile;

//     // Outer loop parameters
//     for (int n = 10; n < 10000; n *= 10) {
//         system(("mkdir -p out/vicsek/" + to_string(n)).c_str());

//         vector<thread> threads;

//         // Parallelize the loop for different values of Noise
//         for (Noise = 0; Noise <= 2 * M_PI; Noise += 0.1) {
//             threads.emplace_back([n, times, Noise]() {
//                 // Create a new Simulation object for each thread
//                 Simulation sim;
//                 sim.NParticles = n;

//                 string filenameNoise = "out/vicsek/" + to_string(sim.NParticles) + "/Noise_" + to_string(Noise) + ".txt";
//                 ofstream ResultsFile(filenameNoise);

//                 for (int i = 1; i <= times; ++i) {
//                     // Use unique seed for each thread
//                     srand(i * static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id())));

//                     // Run the simulation for a given noise
//                     vector<double> Results = sim.RunVicsekForNoise(Noise);

//                     // Use a lock to synchronize file writing
//                     static mutex fileMutex; // Static ensures it's shared among all threads
//                     lock_guard<mutex> lock(fileMutex);

//                     // Write the results to a file
//                     for (int j = 0; j < Results.size(); ++j) {
//                         ResultsFile << Results[j] << " ";
//                     }
//                     ResultsFile << endl;
//                 }
//             });
//         }

//         // Join threads after processing the Noise loop
//         for (auto& thread : threads) {
//             thread.join();
//         }

//         // Clear threads vector for reuse
//         threads.clear();

//         // Parallelize the loop for different values of Density
//         for (Density = 0.1; Density <= 10; Density += 0.1) {
//             threads.emplace_back([n, times, Density]() {
//                 // Create a new Simulation object for each thread
//                 Simulation sim;
//                 sim.NParticles = n;

//                 string filenameDensity = "out/vicsek/" + to_string(sim.NParticles) + "/Density_" + to_string(Density) + ".txt";
//                 ofstream ResultsFile(filenameDensity);

//                 for (int i = 1; i <= times; ++i) {
//                     // Use unique seed for each thread
//                     srand(i * static_cast<unsigned int>(std::hash<std::thread::id>{}(std::this_thread::get_id())));

//                     // Run the simulation for a given density
//                     vector<double> Results = sim.RunVicsekForDensity(Density);

//                     // Use a lock to synchronize file writing
//                     static mutex fileMutex; // Static ensures it's shared among all threads
//                     lock_guard<mutex> lock(fileMutex);

//                     // Write the results to a file
//                     for (int j = 0; j < Results.size(); ++j) {
//                         ResultsFile << Results[j] << " ";
//                     }
//                     ResultsFile << endl;
//                 }
//             });
//         }

//         // Join threads after processing the Density loop
//         for (auto& thread : threads) {
//             thread.join();
//         }
//     }

//     return 0;
// }