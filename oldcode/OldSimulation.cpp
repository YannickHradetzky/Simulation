//
// Created by Yannick Hradetzky on 14.12.23.
//

#include <fstream>
#include "OldSimulation.h"

//----------------------------------------------------------------------------------------------
// Update / Step Functions
void Simulation::VicsekUpdate(double InteractionRadius, double Noise, double Timestep) {
    double ThetaSum, ThetaMean, dr;
    int ThetaCount;
    shared_ptr<Particle> CurrentParticle;
    shared_ptr<Particle> OtherParticle;
    shared_ptr<Cell> CurrentCell;
    shared_ptr<Cell> NeighborCell;
    // Loop over all Particles
    for(int i = 0; i < NParticles; ++i){
        ThetaSum = ThetaMean = ThetaCount = 0;
        CurrentParticle = System->Container->Particles[i];
        CurrentCell = CurrentParticle->mycell;
        // Loop over all Particles in Current Cell
        for (int j = 0; j < CurrentCell->ParticlesInCell.size(); ++j){
            OtherParticle = CurrentCell->ParticlesInCell[j];

            dr = CurrentParticle->CalculateDistance(OtherParticle,
                                                       System->Container->Lx,
                                                       System->Container->Ly,
                                                       System->Container->Lz);
            if (dr < InteractionRadius){
                ThetaSum += OtherParticle->theta;
                ThetaCount++;
            }
        }
        // Loop over all Neighbor Cells
        for(int k = 0; k < System->CellListSystem->Nneighbours; ++k){
            NeighborCell = CurrentParticle->mycell->NeighbourCells[k];
            for(int j = 0; j < NeighborCell->ParticlesInCell.size(); ++j){
                OtherParticle = NeighborCell->ParticlesInCell[j];

                dr = CurrentParticle->CalculateDistance(OtherParticle,
                                                        System->Container->Lx,
                                                        System->Container->Ly,
                                                        System->Container->Lz);
                if (dr < InteractionRadius) {
                    ThetaSum += OtherParticle->theta;
                    ThetaCount++;
                }
            }
        }
        ThetaMean = ThetaSum / ThetaCount;
        CurrentParticle->theta = ThetaMean + (Noise * ((rand() / (RAND_MAX + 1.0)) - 0.5));
        CurrentParticle->vx = cos(CurrentParticle->theta);
        CurrentParticle->vy = sin(CurrentParticle->theta);
        CurrentParticle->x += Timestep * CurrentParticle->vx;
        CurrentParticle->y += Timestep + CurrentParticle->vy;
        CurrentParticle->ApplyCenteredPeriodicBoundaryConditions(System->Container->Lx,
                                                                System->Container->Ly,
                                                                System->Container->Lz);
        System->CellListSystem->SetCellForParticle(i);
    }
    System->Container->UpdateAbsoluteAverageVelocity();
}
//
int Simulation::MCStepHardSpheres(double StepSize) {
    //Perform a Monte Carlo step with step size StepSize for each particle
    //Returns the number of accepted moves

    int NAcceptedMoves = 0;
    double xold, yold, zold, dx, dy, dz, dr;
    shared_ptr<Particle> CurrentParticle;
    for(int i = 0; i < NParticles; ++i){
        CurrentParticle = System->Container->Particles[i];
        dx = StepSize * (rand() / (RAND_MAX + 1.0) - 0.5);
        dy = StepSize * (rand() / (RAND_MAX + 1.0) - 0.5);
        dz = StepSize * (rand() / (RAND_MAX + 1.0) - 0.5);
        xold = CurrentParticle->x;
        yold = CurrentParticle->y;
        zold = CurrentParticle->z;
        // Move Particle
        CurrentParticle->x += dx;
        CurrentParticle->y += dy;
        CurrentParticle->z += dz;
        // Apply Periodic Boundary Conditions
        CurrentParticle->ApplyCenteredPeriodicBoundaryConditions(System->Container->Lx,
                                                                System->Container->Ly,
                                                                System->Container->Lz);
        // Check if Particle overlaps with other Particles
        System->Container->UpdateOverlap();
        if (!System->Container->Overlap){
            //cout << "Accepted Move" << endl;
            System->CellListSystem->SetCellForParticle(i);
            NAcceptedMoves++;
        }
        else{
            // Move Particle back
            //cout << "Rejected Move" << endl;
            CurrentParticle->x = xold;
            CurrentParticle->y = yold;
            CurrentParticle->z = zold;
            Overlap = false;
        }
    }
    return NAcceptedMoves;
}
//



//----------------------------------------------------------------------------------------------
// Run Functions
void Simulation::RunVicsekSimulationVariableNoise(double InteractionRadius, double NoiseBegin, double NoiseEnd,
                                                  double Timestep, double Rho, double Radius, int NCells,
                                                  double Velocity) {
    double NoiseStep = 0.1;
    double Noise;
    double NoiseOfInterest = 2.0;
    // Open Output File
    ofstream OutputAbsVel;

    // Print Information to Screen
    cout << "Vicsek Simulation" << endl;
    cout << "   NParticles: " << NParticles << endl;
    cout << "   NSteps: " << NSteps << endl;
    cout << "   NStepsEquilibration: " << NStepsEquilibration << endl;
    cout << "   NStepsSampling: " << NStepsSampling << endl;

    // Open Output File
    OutputFileObservables = "AbsVelEtaN" + to_string(NParticles) + "Rho" + to_string(Rho) + ".txt";
    OutputFileParticles = "ParticlesEtaN" + to_string(NParticles) + "Rho" + to_string(Rho) + ".txt";
    OutputAbsVel.open(OutputFolder + OutputFileObservables);
    if (!OutputAbsVel.is_open()) {
        std::cerr << "Error opening file: " << OutputFolder + OutputFileObservables << std::endl;
        exit(1);
        return;
    }


    for(Noise = NoiseBegin; Noise <= NoiseEnd; Noise += NoiseStep)
    {

        // Calculate NCells from InteractionRadius and BoxSize from Rho
        double L0 = sqrt(NParticles / Rho);
        NCells = static_cast<int>(L0 / InteractionRadius);


        // Initialize System
        System->InitVicsek(NParticles, Radius, Velocity, Rho, NCells);
        cout << "Initializing System" << endl;
        cout << "   Box Dimensions: " << endl;
        cout << "      Lx: " << System->Container->Lx << endl;
        cout << "      Ly: " << System->Container->Ly << endl;
        cout << "   Rho: " << System->Container->CalculateDensity() << endl;
        cout << "   Noise: " << Noise << endl;
        cout << "   NCells: " << NCells << endl;



        // Equilibrate System
        cout << "Equilibrating System" << endl;
        int Times = 5;
        int NPrints = NStepsEquilibration / 5;
        for(int i = 0; i <= NStepsEquilibration; ++i){
            VicsekUpdate(InteractionRadius, Noise, Timestep);
            // print Absolut Average Velocity every 10 Steps to Console
            if (i % NPrints == 0){
                cout <<  "   Step: " << i <<"   AbsVel: " << System->Container->AbsolutAverageVelocity << endl;
            }

        }
        // Main Loop
        cout << "Sampling System" << endl;
        OutputAbsVel << Noise << " ";
        for(int i = 0; i <= NSteps; ++i)
        {
            VicsekUpdate(InteractionRadius, Noise, Timestep);
            if (i % NStepsSampling == 0)
            {
                if(ShowProgressBar) PrintProgressBar(i, NSteps, 20);
                OutputAbsVel << System->Container->AbsolutAverageVelocity << " ";
            }
        }
        OutputAbsVel << endl;
    }
    OutputAbsVel.close();
}
//
void
Simulation::RunVicsekSimulationVariableDensity(double InteractionRadius, double Noise, double Timestep, double RhoBegin,
                                               double RhoEnd, double Radius, int NCells, double Velocity) {
    double RhoStep = 0.1;
    double Rho;

    // Open Output File
    ofstream OutputPosition;
    ofstream OutputAbsVel;

    // Print Information to Screen
    cout << "Vicsek Simulation" << endl;
    cout << "   NParticles: " << NParticles << endl;
    cout << "   NSteps: " << NSteps << endl;
    cout << "   NStepsEquilibration: " << NStepsEquilibration << endl;
    cout << "   NStepsSampling: " << NStepsSampling << endl;

    // Open Output File
    OutputFileObservables = "AbsVelRhoN" + to_string(NParticles) + "Noise" + to_string(Noise) + ".txt";
    OutputAbsVel.open(OutputFolder + OutputFileObservables);

    if (!OutputAbsVel.is_open()) {
        std::cerr << "Error opening file: " << OutputFolder + OutputFileObservables << std::endl;
        exit(1);
        return;
    }

    for(Rho = RhoBegin; Rho <= RhoEnd; Rho += RhoStep)
    {
        // Calculate NCells from InteractionRadius and BoxSize from Rho
        double L0 = sqrt(NParticles / Rho); // Calculate BoxSize from Rho
        NCells = static_cast<int>(L0 / InteractionRadius); // Calculate NCells from InteractionRadius and BoxSize
        // Check if NCells is at least 3
        if (NCells < 3){
            cout << "NCells < 3" << endl;
            NCells = 3;
        }
        if (L0 / NCells < InteractionRadius){ // Check if LCutOffx is smaller than InteractionRadius
            cout << "Warning:" << endl;
            cout << "   LCutOffx < InteractionRadius" << endl;
            cout << "   LCutOffx: " << System->CellListSystem->LCutOffx << endl;
            cout << "   InteractionRadius: " << InteractionRadius << endl;
        }

        // Initialize System with NParticles, Radius, Velocity, Rho, NCells
        System->InitVicsek(NParticles, Radius, Velocity, Rho, NCells);
        cout << "Initializing System" << endl;
        cout << "   Box Dimensions: " << endl;
        cout << "      Lx: " << System->Container->Lx << endl;
        cout << "      Ly: " << System->Container->Ly << endl;
        cout << "   Rho: " << System->Container->CalculateDensity() << endl;
        cout << "   Noise: " << Noise << endl;
        cout << "   NCells: " << NCells << endl;

        // Equilibrate System
        cout << "Equilibrating System" << endl;
        int Times = 5;
        int NPrints = NStepsEquilibration / 5;
        for(int i = 0; i <= NStepsEquilibration; ++i){
            VicsekUpdate(InteractionRadius, Noise, Timestep);
            // print Absolut Average Velocity every 10 Steps to Console
            if (i % NPrints == 0){
                cout <<  "   Step: " << i <<"   AbsVel: " << System->Container->AbsolutAverageVelocity << endl;
            }
        }
        // Main Loop
        cout << "Sampling System" << endl;
        OutputAbsVel << Rho << " ";
        for(int i = 0; i <= NSteps; ++i)
        {
            VicsekUpdate(InteractionRadius, Noise, Timestep);
            if (i % NStepsSampling == 0)
            {
                if(ShowProgressBar) PrintProgressBar(i, NSteps, 20);
                OutputAbsVel << System->Container->AbsolutAverageVelocity << " ";
            }
        }
        OutputAbsVel << endl;
    }
    OutputAbsVel.close();
}
//
void Simulation::RunHardSpheresSimulation(double Radius, double Phi, double StepSize) {

    cout << "Hard Spheres Simulation" << endl;
    cout << "   NParticles: " << NParticles << endl;
    cout << "   NSteps: " << NSteps << endl;
    cout << "   NStepsEquilibration: " << NStepsEquilibration << endl;
    cout << "   NStepsSampling: " << NStepsSampling << endl;

    // Open Output File
    ofstream MeanDistance;
    MeanDistance.open(OutputFolder + OutputFileObservables);
    if (!MeanDistance.is_open()) {
        std::cerr << "Error opening file: " << OutputFolder + OutputFileObservables << std::endl;
        exit(1);
        return;
    }

    // Initialize System
    System->InitHardSpheres(NParticles, Radius, Phi);
    cout << "Initializing System" << endl;
    cout << "   Box Dimensions: " << endl;
    cout << "      Lx: " << System->Container->Lx << endl;
    cout << "      Ly: " << System->Container->Ly << endl;
    cout << "      Lz: " << System->Container->Lz << endl;
    cout << "      NCells: " << System->CellListSystem->Ncellsx << endl;
    cout << "   Rho: " << System->Container->CalculateDensity() << endl;
    cout << "   Phi: " << System->Container->CalculatePackingFraction() << endl;
    cout << "   Radius of Particle[0]: " << System->Container->Particles[0]->myradius << endl;

    cout << "Calculating Initial RDF" << endl;
    System->Container->ExportCalculateRDF(OutputFolder + "RDFInitial.txt", ShowPlots);


    // Equilibrate System to find the right StepSize
    cout << "Equilibrating System" << endl;
    int Times = 5;
    int NPrints = NStepsEquilibration / Times;
    int NAcceptedMoves = 0, NMonteCarloSteps = 0;
    double AcceptanceRatio;
    for(int i = 0; i <= NStepsEquilibration; ++i)
    {
        NAcceptedMoves += MCStepHardSpheres(StepSize);
        NMonteCarloSteps += NParticles;
        AcceptanceRatio = ((double )NAcceptedMoves / NMonteCarloSteps);
        // Change StepSize if AcceptanceRatio is too high or too low
        if (AcceptanceRatio > 0.55){
            StepSize *= 1.1;
        }
        else if (AcceptanceRatio < 0.45){
            StepSize *= 0.9;
        }
        if (i % NPrints == 0){
            cout << "   Step: " << i << endl;
            cout << "       AcceptanceRatio: " << AcceptanceRatio << "      StepSize: " << StepSize << endl;
        }
    }
    System->Container->ExportCalculateRDF(OutputFolder + "RDFEquil.txt", ShowPlots);

    // Main Loop
    cout << "Sampling System" << endl;
    for(int i = 0; i <= NSteps; ++i)
    {
        MCStepHardSpheres(StepSize);
        if (ShowProgressBar) PrintProgressBar(i, NSteps, 20);
        if (i % NStepsSampling == 0)
        {
            System->Container->UpdateMeanMinDistance();
            MeanDistance << System->Container->MeanMinDistance << endl;
        }
    }
    MeanDistance.close();

    // Calculate RDF
    System->Container->ExportCalculateRDF(OutputFolder + "RDFFinal.txt", ShowPlots);

}


//----------------------------------------------------------------------------------------------
// Helper Functions
void Simulation::PrintProgressBar(int Step, int TotalSteps, int Width) {
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







