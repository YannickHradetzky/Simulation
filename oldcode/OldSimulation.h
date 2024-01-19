//
// Created by Yannick Hradetzky on 14.12.23.
//

#ifndef SIMULATION_OLDSIMULATION_H
#define SIMULATION_OLDSIMULATION_H

#include "OldSystem.h"

using namespace std;

class OldSimulation: public System{
public:
    OldSimulation(System *System0){
        System = System0;
    }
    OldSimulation(){
        System = new class System;
    }
    ~OldSimulation(){
        delete System;
    }
    System *System;
    //----------------------------------------------------------------------------------------------
    // OldSimulation Variables
    int NParticles;
    int NSteps;
    int NStepsEquilibration;
    int NStepsSampling;
    int NStepsThermo;
    // Output Variables
    string OutputFolder;
    string OutputFileObservables;
    string OutputFileParticles;
    // Settings
    bool ShowProgressBar;
    bool ShowPlots;

    //----------------------------------------------------------------------------------------------
    // Update / Step Functions
    void VicsekUpdate(double InteractionRadius, double Noise, double Timestep);
    int MCStepHardSpheres(double StepSize);

    //----------------------------------------------------------------------------------------------
    // Run Functions
    void RunVicsekSimulationVariableNoise(double InteractionRadius, double NoiseBegin, double NoiseEnd,
                                          double Timestep, double Rho, double Radius, int NCells,
                                          double Velocity);
    void RunVicsekSimulationVariableDensity(double InteractionRadius, double Noise, double Timestep,
                                            double RhoBegin, double RhoEnd, double Radius, int NCells,
                                            double Velocity);
    void RunHardSpheresSimulation(double Radius, double Phi, double StepSize);

private:
    //----------------------------------------------------------------------------------------------
    // Helper Functions
    void PrintProgressBar(int i, int NSteps, int Width);
};

#endif //SIMULATION_OLDSIMULATION_H
