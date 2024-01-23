//
// Created by Yannick Hradetzky on 20.12.23.
//

#include "Simulation.h"

std::vector<double> Simulation::RunVicsekForNoise(double Noise){
    double VelocityInit = 0.3;
    InitVicsek(VelocityInit);
    // Equilibrate the System
    for(int i = 0; i <= NStepsEquilibration; ++i){
        VicsekUpdate(Noise);
    }
    std::vector<double> Results;
    for (int i = 0; i <= NStepsSimulation; ++i){
        VicsekUpdate(Noise);
        // Calculate absolut average velocity
        if(i % NStepsSampling == 0){
            
            double sumx = 0,  sumy = 0;
            for(const auto &Current : Particles)
            {
                sumx += Current->vx;
                sumy += Current->vy;
            }
            double vabs = sqrt(sumx*sumx + sumy*sumy) / (NParticles * VelocityInit);
            Results.push_back(vabs);
        } 
    }
    return Results;
}

std::vector<double> Simulation::RunVicsekForDensity(double Density, double NoiseForDensity){
    // Initialize the System
    double VelocityInit = 0.3;
    double Noise = NoiseForDensity;
    DensityInit = Density;
    InitVicsek(VelocityInit);
    // Equilibrate the System
    for(int i = 0; i <= NStepsEquilibration; ++i){
        VicsekUpdate(Noise);
    }
    std::vector<double> Results;
    for (int i = 0; i <= NStepsSimulation; ++i){
        VicsekUpdate(Noise);
        if(i % NStepsSampling == 0){
            // Calculate absolut average velocity
            double sumx = 0,  sumy = 0;
            for(const auto &Current : Particles)
            {
                sumx += Current->vx;
                sumy += Current->vy;
            }
            double vabs = sqrt(sumx*sumx + sumy*sumy) / (NParticles * VelocityInit);
            Results.push_back(vabs);
        } 
    }
    return Results;
}


std::vector<double> Simulation::ComputeOrientationCorrelationVicsek(double AvgVelocityX, double AvgVelocityY){
    int NBins = 200; 
    double RStep = 0.5 * L / NBins;

    std::vector<double>Correlation (NBins, 0.0); 
    std :: vector<double> CorrelationCount (NBins, 0.0);
    for(const auto &Current : Particles){
        for(const auto &Other : Particles){
            if(Current != Other){
                double dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                int Bin = (int) (dr / RStep + 0.5);
                if(Bin < NBins ) {
                    Correlation[Bin] += Current->theta - Other->theta;
                    CorrelationCount[Bin] += 1;
                }
            }
        }
    }
    for(int i = 0; i < NBins; ++i){
        if(CorrelationCount[i] != 0){
            Correlation[i] /= CorrelationCount[i];
        }
    }
    return Correlation;
}


void Simulation::RunHardSphere() {
    std::cout << "Hard Sphere OldSimulation" << std::endl;
    std::cout << "  with " << NParticles << " particles" << std::endl;
    std::cout << "  for " << NStepsSimulation << " steps" << std::endl;
    std::cout << "  with " << NStepsEquilibration << " equilibration steps" << std::endl;
    std::cout << "  and " << NStepsSimulation / NStepsSampling << " sampling steps" << std::endl;

    OutputFolder = "out/hardspheres/";
    OutputFileObservables = "mean_distance.txt";

    // Initialize the System
    InitLattice(0.5);
    CalculateExportRDF(OutputFolder + "RDFInitialisation.txt", ShowPlots, true);

    std::cout << "Equilibrating System" << std::endl;
    int Times = 5;
    int NPrints = NStepsEquilibration / Times;
    int NAcceptedMoves = 0, NMonteCarloSteps = 0;
    double AcceptanceRatio;
    double StepSize = 0.1;
    for(int i = 0; i <= NStepsEquilibration; ++i)
    {
        NAcceptedMoves += MCNParticleMove(StepSize);
        NMonteCarloSteps += Particles.size();
        if(ParticlesOverlap) {
            std::cout << "Overlap" << std::endl;
            std::cout << "  Step: " << i << std::endl;
            exit(1);
        }
        AcceptanceRatio = ((double )NAcceptedMoves / NMonteCarloSteps);
        if(AcceptanceRatio > 0.55) StepSize *= 1.1;
        else if(AcceptanceRatio < 0.45) StepSize *= 0.9;
        if (i % NPrints == 0 && i != 0){
            // Change StepSize if AcceptanceRatio is too high or too low
            std::cout << "   Step: " << i << std::endl;
            std::cout << "       AcceptanceRatio: " << AcceptanceRatio << "      StepSize: " << StepSize << std::endl;
            NAcceptedMoves = NMonteCarloSteps = 0;
        }
    }
    CalculateExportRDF(OutputFolder + "RDFEquilibration.txt", ShowPlots, true);

    std::cout << "Sampling System" << std::endl;
    std::ofstream OutputFile;
    OutputFile.open(OutputFolder + OutputFileObservables);
    for(int i = 0; i <= NStepsSimulation; ++i)
    {
        if(PrintProgress) PrintProgressBar(i, NStepsSimulation, 30);
        MCNParticleMove(StepSize);
        if(ParticlesOverlap) {
            std::cout << "Overlap" << std::endl;
            std::cout << "  Step: " << i << std::endl;
            exit(1);
        }
        if(i % NStepsSampling == 0){
            CalculateExportRDF(OutputFolder + "RDF/" + "RDF" + std::to_string(i) + ".txt", false, false);
            OutputFile << CalculateMinMeanDistance(true) << std::endl;
        }
    }
    OutputFile.close();
    CalculateExportRDF(OutputFolder + "RDFinal.txt", ShowPlots, true);
}
void Simulation::RunNPTPressure(double PressureStart, double PressureEnd, double PressureStep) {
    std::cout << "NPT Pressure Simulation" << std::endl;
    std::cout << "  Pressure from " << PressureStart << " to " << PressureEnd << std::endl;
    std::cout << "  with " << NParticles << " particles" << std::endl;
    std::cout << "  for " << NStepsSimulation << " steps" << std::endl;
    std::cout << "  with " << NStepsEquilibration << " equilibration steps" << std::endl;
    std::cout << "  and " << NStepsSimulation / NStepsSampling << " sampling steps" << std::endl;

    OutputFolder = "out/NPT/";


    double ParticleRadius = 0.5;
    double StepSizeParticleMove = 0.1;
    double StepSizeVolumeMove = 0.0001;

    for(double PressureSim = PressureStart; PressureSim <= PressureEnd; PressureSim += PressureStep){
        Pressure = PressureSim;
        std::cout << "Simulating Pressure: " << Pressure << std::endl;
        std::cout << "and Temperature: " << Temperature << std::endl;
        InitLattice(ParticleRadius);
        OutputFileObservables = "NPT_of_P_" + std::to_string(Pressure) +
                                "_forN_" + std::to_string(NParticles) + ".txt";

        std::ofstream OutputFile;
        OutputFile.open(OutputFolder + OutputFileObservables);
        OutputFile << "Step PackingFraction Volume Density Pressure Temperature" << std::endl;

        int NaccParticle = 0, NaccVolume = 0, NmcParticle = 0, NmcVolume = 0;
        for(int i = 0; i <= NStepsSimulation; ++i)
        {

            int Feedback = MCStepNPT(StepSizeParticleMove, StepSizeVolumeMove);
            if(Feedback == 0) NmcParticle++;
            else if(Feedback == 1) NaccParticle++, NmcParticle++;
            else if(Feedback == 2) NmcVolume++;
            else if(Feedback == 3) NaccVolume++, NmcVolume++;


            if(i % NStepsSampling == 0) {
                // Write Observables to File
                OutputFile << i << " " << PackingFraction << " " << L * L * L << " "
                            << (double) NParticles / (L * L * L) << " " << Pressure << " " << Temperature << std::endl;




                if (i <= NStepsEquilibration) {
                    // Change StepSize if AcceptanceRatios are too high or too low
                    if (NmcParticle == 0) NmcParticle = 1;
                    if (NmcVolume == 0) NmcVolume = 1;

                    if ((NaccParticle / (double) NmcParticle) > 0.55) {
                        if (StepSizeParticleMove * 1.1 <= 0.5 * ParticleRadius) {
                            StepSizeParticleMove *= 1.1; // 0.5 times the radius is the maximum step size
                        } else {
                            StepSizeParticleMove = 0.5 * ParticleRadius;
                        }
                        NaccParticle = NmcParticle = 0;
                    }
                    if ((NaccParticle / (double) NmcParticle) < 0.45) {
                        StepSizeParticleMove *= 0.9;
                        NaccParticle = NmcParticle = 0;
                    }
                    if ((NaccVolume / (double) NmcVolume) > 0.35) {
                        if (StepSizeVolumeMove * 1.2 < 0.1) {
                            StepSizeVolumeMove *= 1.2; // 0.1 is the maximum step size

                        }else{
                            StepSizeVolumeMove = 0.1;
                        }
                        NaccVolume = NmcVolume = 0;
                    }
                    if((NaccVolume / (double) NmcVolume) < 0.25) {
                        StepSizeVolumeMove *= 0.9;
                        NaccVolume = NmcVolume = 0;
                    }
                }

                if(PrintProgress && i != 0 && i % (NStepsSampling * 20) == 0) {
                    std::cout<< std::endl;
                    std::cout << "  Step: " << i << std::endl;
                    std::cout << "      Length: " << L << std::endl;
                    std::cout << "      PackingFraction: " << CalculatePackingFraction(true) << std::endl;
                    std::cout << "      Volume: " << L * L * L << std::endl;
                    std::cout << "      Density: " << (double) NParticles / (L * L * L) << std::endl;
                    std::cout << "      AcceptanceRatioParticle: " << NaccParticle / (double)NmcParticle << std::endl;
                    std::cout << "          StepSizeParticleMove: " << StepSizeParticleMove << std::endl;
                    std::cout << "      AcceptanceRatioVolume: " << NaccVolume / (double)NmcVolume << std::endl;
                    std::cout << "          StepSizeVolumeMove: " << StepSizeVolumeMove << std::endl;
                    std::cout<< std::endl;
                    std::cout<< std::endl;
                }
            }
        }
        // Print Information about the Simulation into the OutputFile like Steps, EquilibrationSteps, SamplingSteps, ...
        OutputFile << "Steps EquilibrationSteps SamplingSteps" << std::endl;
        OutputFile << NStepsSimulation << " " << NStepsEquilibration << " " << NStepsSampling << std::endl;
        OutputFile.close();
    }
    std::cout << "Finished Simulation" << std::endl;
}

