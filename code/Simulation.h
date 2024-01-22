//
// Created by Yannick Hradetzky on 20.12.23.
//

#ifndef SIMULATION_SIMULATION_H
#define SIMULATION_SIMULATION_H

//----------------------------------------------------------------------------------------------
// Include Files
#ifndef System
#include "System.h"
#endif

#ifndef iomanip
#include "iomanip"
#endif

#ifndef iostream
#include "iostream"
#endif

#ifndef fstream
#include "fstream"
#endif

#ifndef dirent
#include "dirent.h"
#endif

#ifndef unistd
#include "unistd.h"
#endif

//----------------------------------------------------------------------------------------------

class Simulation: public System {
public:
    Simulation() = default; // Default Constructor
    ~Simulation() = default; // Default Destructor
    //----------------------------------------------------------------------------------------------
    // Variables
    int NStepsSimulation, NStepsEquilibration, NStepsSampling;
    double TimeStepSize;
    // Export Variables
    std::string OutputFolder;
    std::string OutputFileObservables;
    std::string OutputFilePositions;
    std::string OutputFileVelocities;
    //----------------------------------------------------------------------------------------------
    // Boolean Variables
    bool PrintProgress = true;
    bool ShowPlots = true;
    //----------------------------------------------------------------------------------------------
    // Simulation Functions
    void VicsekUpdate(double Noise) {
        double ThetaSum, ThetaMean, dr;
        int Count;
        for (const auto &Current : Particles) {
            ThetaMean = ThetaSum = 0;
            Count = 0;
            for (const auto &Other : Current->mycell->ParticlesInCell) {
                dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < InteractionRadius) {
                    ThetaSum += Other->theta;
                    Count++;
                }
            }
            for (const auto &NeighborCell : Current->mycell->NeighbourCells) {
                for (const auto &Other : NeighborCell->ParticlesInCell) {
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    if (dr < InteractionRadius) {
                        ThetaSum += Other->theta;
                        Count++;
                    }
                }
            }
            if (Count > 0) ThetaMean = ThetaSum / Count;
            Current->theta = ThetaMean + Noise * (rand() / (double)RAND_MAX);
            
            // Update positions using the calculated velocities
            Current->vx = Current->velocity * cos(Current->theta);
            Current->vy = Current->velocity * sin(Current->theta);
            Current->x += Current->vx * TimeStepSize;
            Current->y += Current->vy * TimeStepSize;
            Current->ApplyCenteredPeriodicBoundaryConditions(Lx, Ly, Lz);
            UpdateCell(Current);
        }
    }
    int MCNParticleMove(double StepSize){
        double xold, yold, zold;
        int Naccepted = 0;
        int count = 0;
        for(const auto &Current : Particles){
            count++;
            xold = Current->x; yold = Current->y; zold = Current->z;
            Current->x += StepSize * (rand() / (double) RAND_MAX - 0.5);;
            Current->y += StepSize * (rand() / (double) RAND_MAX - 0.5);;
            Current->z += StepSize * (rand() / (double) RAND_MAX - 0.5);;
            Current->ApplyCenteredPeriodicBoundaryConditions(Lx, Ly, Lz);
            UpdateCell(Current); // Update Cell of Particle
            UpdateOverlap(); //
            if(!ParticlesOverlap){
                Naccepted ++;
            }else if(ParticlesOverlap){
                Current->x = xold;
                Current->y = yold;
                Current->z = zold;
                UpdateCell(Current);
            }
            UpdateOverlap();
            if(ParticlesOverlap){
                std::cout << "Overlap not resolved in Step: " << count  << std::endl;
                std::cout << "Particle " << Current->myindex << std::endl;
                std::cout << "  xold: " << xold << " current " << Current->x << std::endl;
                std::cout << "  yold: " << yold << " current " << Current->y << std::endl;
                std::cout << "  zold: " << zold << " current " << Current->z << std::endl;
                exit(1);
            }
        }
        return Naccepted;
    }
    bool MCParticleMove(double StepSize){
        bool Accepted = false;
        // Choose random particle
        std::shared_ptr<Particle> Chosen = Particles[(rand()/ (double) RAND_MAX) * Particles.size()];
        double xold = Chosen->x, yold = Chosen->y, zold = Chosen->z;
        // Move the particle
        Chosen->x += StepSize * (rand() / (double) RAND_MAX - 0.5);
        Chosen->y += StepSize * (rand() / (double) RAND_MAX - 0.5);
        Chosen->z += StepSize * (rand() / (double) RAND_MAX - 0.5);
        Chosen->ApplyCenteredPeriodicBoundaryConditions(Lx, Ly, Lz);
        UpdateCell(Chosen); // Update Cell of Particle
        UpdateOverlap(Chosen); // Update ParticlesOverlap from Container
        if (!ParticlesOverlap) {
            Accepted = true;
        } else if (ParticlesOverlap) {
            Chosen->x = xold;
            Chosen->y = yold;
            Chosen->z = zold;
            UpdateCell(Chosen);
        }
        return Accepted;
    }
    bool MCVolumeMove(double StepSize) {
        bool Accepted = false;
        double Vold = L * L * L;
        double Random1 = (rand() / (double) RAND_MAX - 0.5); // Random Number between -0.5 and 0.5
        double Random2 = (rand() / (double) RAND_MAX); // Random Number between 0 and 1

        // Calculate new Volume from Random Number
        double Vnew = exp(log(Vold) + Random1 * StepSize);
        double AcceptCoeff = exp(-Beta * Pressure * (Vnew - Vold) + (NParticles + 1) * log(Vnew / Vold));
        double Lnew = pow(Vnew, 1.0 / 3.0);
        double Lold = pow(Vold, 1.0 / 3.0);

        // Scale the Container (Coordinates) and check if particles overlap
        Lx = Ly = Lz = L = Lnew;
        for (const auto &Current: Particles) {
            Current->x *= Lnew / Lold;
            Current->y *= Lnew / Lold;
            Current->z *= Lnew / Lold;
        }
        if(Random1 < 0) UpdateOverlap(); // Only check for overlap if the volume is decreased
        if (ParticlesOverlap) {
            // Scale the Container (Coordinates) back
            Lx = Ly = Lz = L = Lold;
            for (const auto &Current: Particles) {
                Current->x *= Lold / Lnew;
                Current->y *= Lold / Lnew;
                Current->z *= Lold / Lnew;
            }
        } else{
            if (Random2 < AcceptCoeff) {
                // Accept the move with probability AcceptCoeff and update the system
                Accepted = true;
                Surface = 6 * L * L;
                Density = NParticles / (L*L*L);
                Volume = L*L*L;
                CalculatePackingFraction(true);
            } else {
                // Scale the Container (Coordinates) back
                L = Lx = Ly = Lz = Lold;
                for (const auto &Current: Particles) {
                    Current->x *= Lold / Lnew;
                    Current->y *= Lold / Lnew;
                    Current->z *= Lold / Lnew;
                }
            }
        }
        return Accepted;
    }
    int MCStepNPT(double StepSizeParticle, double StepSizeVolume){
        double Random = (rand() / (double) RAND_MAX); // Random Number between 0 and 1
        bool ParticleMove = false;
        bool VolumeMove = false;
        bool VolumeAccepted = false;
        bool ParticleAccepted = false;
        int Feedback;

        if(Random < 1./ NParticles){
            VolumeMove = true;
            VolumeAccepted = MCVolumeMove(StepSizeVolume);
        }else{
            ParticleMove = true;
            ParticleAccepted = MCParticleMove(StepSizeParticle);
        }
        //------------------------------------------------------------------------------------------
        // Check for Errors
        if(ParticleMove && VolumeMove){
            std::cout << "Error in MCStepNPT: Both Particle and Volume Move" << std::endl;
            exit(1);
        }
        //------------------------------------------------------------------------------------------

        if(ParticleMove){
            if(ParticleAccepted) Feedback =  1;
            else Feedback = 0;
        }
        if(VolumeMove){
            if(VolumeAccepted) Feedback = 3;
            else Feedback =  2;
        }
        return Feedback;
    }
    void VerletIntegrator(double Timestep){
        // First Half Step : Update Velocities
        for(const auto &Current : Particles){
            Current->vx += 0.5 * Current->fx * Timestep / Current->mymass;
            Current->vy += 0.5 * Current->fy * Timestep / Current->mymass;
            Current->vz += 0.5 * Current->fz * Timestep / Current->mymass;
        }
        // Full Step : Update Positions and Cells
        for(const auto &Current : Particles){
            Current->x += Current->vx * Timestep;
            Current->y += Current->vy * Timestep;
            Current->z += Current->vz * Timestep;
            Current->ApplyPeriodicBoundaryConditions(Lx, Ly, Lz);
        }
        // Update Forces and Cells
        InitCells(InteractionRadius);
        CalculateForces();
        // Second Half Step : Update Velocities
        for(const auto &Current : Particles){
            Current->vx += 0.5 * Current->fx * Timestep / Current->mymass;
            Current->vy += 0.5 * Current->fy * Timestep / Current->mymass;
            Current->vz += 0.5 * Current->fz * Timestep / Current->mymass;
        }
    }

    //----------------------------------------------------------------------------------------------
    // Helper Functions
    void PrintProgressBar(int Step, int TotalSteps, int Width) {
        float progress = static_cast<float>(Step) / TotalSteps;
        int progressBarWidth = static_cast<int>(progress * Width);

        std::cout << "\r[";
        for (int i = 0; i < progressBarWidth; ++i) {
            std::cout << "=";
        }
        for (int i = progressBarWidth; i < Width; ++i) {
            std::cout << " ";
        }

        std::cout << "] " << std::fixed << std::setprecision(1) << (progress * 100.0) << "%";
        std::cout.flush();

        if (Step == TotalSteps) {
            std::cout << std::endl; // Move to the next line after the loop completes
        }
    }

    //----------------------------------------------------------------------------------------------
    // PotentialFunction Functions


    //----------------------------------------------------------------------------------------------
    //Run Functions (defined in Simulation.cpp)
    std::vector<double> RunVicsekForNoise(double Noise);
    std::vector<double> RunVicsekForDensity(double Density);
    void RunHardSphere();
    void RunNPTPressure(double PressureStart, double PressureEnd, double PressureStep);

    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
};
#endif //SIMULATION_SIMULATION_H
