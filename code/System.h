//
// Created by Yannick Hradetzky on 19.12.23.
//

#ifndef SIMULATION_SYSTEM_H
#define SIMULATION_SYSTEM_H

#ifndef Container
#include "Container.h"
#endif

#ifndef fstream
#include "fstream"
#endif




// define general PotentialFunction and ForceFunction Functions
using PotentialFunction1 = double (*)(double); // for example distance r between two particles
using ForceFunction1 = double (*)(double); // for example distance r between two particles




class System : public Container {
public:
    System() = default;
    ~System() = default;
    //----------------------------------------------------------------------------------------------
    // Variables
    double  Temperature = 0, PackingFraction = 0, Density = 0;
    double KineticEnergy = 0, PotentialEnergy = 0, TotalEnergy = 0;
    double MinMeanDistance = 0;
    double Beta = 0, Pressure = 0;
    double InteractionRadius = 1;
    //----------------------------------------------------------------------------------------------
    // Init Variables
    int NParticles;
    double  TemperatureInit, PackingFractionInit, DensityInit;
    double PressureInit;
    //----------------------------------------------------------------------------------------------
    // Boolean Variables
    bool PrintInitInfo = true;
    //----------------------------------------------------------------------------------------------
    // Dynamic Variables
    PotentialFunction1 PotentialFunction = nullptr;
    ForceFunction1 ForceFunction = nullptr;

    //------------------------------- ---------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
    // Init Models
    void InitLattice(double Radius) {
        std::cout << "Initialize Lattice" << std::endl;
        Particles.clear();

        if (NParticles == 0){
            std::cout << "  No NParticles given, set to 27" << std::endl;
            NParticles = 27;
        }else{
            std::cout << "  with " << NParticles << " particles" << " and Radius " << Radius << std::endl;
        }
        int ParticlesPerDim = (int) pow(NParticles + 0.00001, 1.0/3.0);
        if(pow(ParticlesPerDim, 3) != NParticles) {
            std::cout << "  NParticles is not a perfect cube! " << std::endl;
            exit(1);
        }
        double L0 = ParticlesPerDim * Radius * 2;

        // Initialize Particles
        double xsum = 0, ysum = 0, zsum = 0;
        for(int i = 0; i < ParticlesPerDim; ++i) {
            for(int j = 0; j < ParticlesPerDim; ++j) {
                for(int k = 0; k < ParticlesPerDim; ++k) {
                    int index = i * ParticlesPerDim * ParticlesPerDim + j * ParticlesPerDim + k;
                    Particle tempParticle;
                    tempParticle.myindex = index;
                    tempParticle.myradius = Radius;
                    tempParticle.x = i * Radius * 2;
                    tempParticle.y = j * Radius * 2;
                    tempParticle.z = k * Radius * 2;
                    xsum += tempParticle.x;
                    ysum += tempParticle.y;
                    zsum += tempParticle.z;
                    Particles.push_back(std::make_shared<Particle>(tempParticle));
                }
            }
        }
        // Set Center of Mass to Zero
        xsum /= NParticles;
        ysum /= NParticles;
        zsum /= NParticles;
        for(int i = 0; i < NParticles; ++i) {
            Particles[i]->x -= xsum;
            Particles[i]->y -= ysum;
            Particles[i]->z -= zsum;
        }
        // Rescale the Box to the desired Packing Fraction
        if(PackingFractionInit != 0.0){
            L = pow(NParticles * 4.0 / 3.0 * M_PI * pow(Radius, 3) / PackingFractionInit, 1.0/3.0);
            for(int i = 0; i < NParticles; ++i) {
                Particles[i]->x *= L / L0;
                Particles[i]->y *= L / L0;
                Particles[i]->z *= L / L0;
            }
            PackingFraction = NParticles * 4.0 / 3.0 * M_PI * pow(Radius, 3) / pow(L, 3);
        }
        else{
            L = L0;
            PackingFraction = NParticles * 4.0 / 3.0 * M_PI * pow(Radius, 3) / pow(L, 3) - 0.01;
            std::cout << "  No PackingFractionInit given " << std::endl;
            std::cout << "  Set to slightly below maximum where particles touch " << PackingFraction << std::endl;
            for(int i = 0; i < NParticles; ++i) {
                Particles[i]->x *= L / L0;
                Particles[i]->y *= L / L0;
                Particles[i]->z *= L / L0;
            }
        }
        // Compute Container Properties
        Volume = L*L*L;
        Surface = L*L*6;
        Area = nan("Not 2D");
        Lx = Ly = Lz = L;
        Dim = 3;
        // Compute System Properties
        CalculateTemperature(true);
        CalculateKineticEnergy(true);
        CalculatePackingFraction(true);
        CalculateTemperature(true);
        CalculateMinMeanDistance(true);
        Density = NParticles / Volume;
        // Initialize Cells
        InteractionRadius = 2.5 * Radius;
        InitCells(InteractionRadius);
        if(PrintInitInfo){
            // Print Container Properties
            std::cout << "  Container Properties: " << std::endl;
            std::cout << "      L = " << L << std::endl;
            std::cout << "      Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << std::endl;
            std::cout << "      Volume = " << Volume << std::endl;
            std::cout << "      Surface = " << Surface << std::endl;
            // Print Cells Properties
            std::cout << "  Cells Properties: " << std::endl;
            std::cout << "      Ncellsx = " << Ncellsx << std::endl;
            std::cout << "      Ncellsy = " << Ncellsy << std::endl;
            std::cout << "      Ncellsz = " << Ncellsz << std::endl;
            // Print System Properties
            std::cout << "  System Properties: " << std::endl;
            std::cout << "      Temperature = " << Temperature << std::endl;
            std::cout << "      Kinetic Energy = " << KineticEnergy << std::endl;
            std::cout << "      PotentialFunction Energy = " << PotentialEnergy << std::endl;
            std::cout << "      Total Energy = " << TotalEnergy << std::endl;
            std::cout << "      Packing Fraction = " << PackingFraction << std::endl;
            std::cout << "      Density = " << Density << std::endl;
        }
    }
    void InitVicsek(double Velocity)
    {
        //std::cout << "Initialize Vicsek" << std::endl;
        Particles.clear();

        if (NParticles == 0){
            std::cout << "  No NParticles given, set to 100" << std::endl;
            NParticles = 100;
        }else{
            //std::cout << "  with " << NParticles << " particles" << std::endl;
        }
        if (DensityInit == 0){
            std::cout << "  No DensityInit given, set to 1" << std::endl;
            DensityInit = 1;
            Density = DensityInit;
            L = pow(NParticles / DensityInit, 1.0/2.0);
            Lx = Ly = L;
            Lz = 0;
            Volume = nan("Not 3D");
            Surface = nan("Not 3D");
            Area = L*L;
            Dim = 2;
        }else{
            //std::cout << "  with Density " << DensityInit << std::endl;
            Density = DensityInit;
            L = pow(NParticles / DensityInit, 1.0/2.0);
            Lx = Ly = L;
            Lz = 0;
            Volume = nan("Not 3D");
            Surface = nan("Not 3D");
            Area = L*L;
            Dim = 2;
        }
        // Initialize Particles
        double sumx = 0, sumy = 0;
        for(int i = 0; i < NParticles; ++i) {
            Particle tempParticle;
            tempParticle.myindex = i;
            tempParticle.myradius = 0.0;
            tempParticle.x = (double) rand() / RAND_MAX * L;
            tempParticle.y = (double) rand() / RAND_MAX * L;
            sumx += tempParticle.x;
            sumy += tempParticle.y;
            tempParticle.theta = (double) rand() / RAND_MAX * 2 * M_PI;
            tempParticle.velocity = Velocity;
            tempParticle.vx = Velocity * cos(tempParticle.theta);
            tempParticle.vy = Velocity * sin(tempParticle.theta);
            Particles.push_back(std::make_shared<Particle>(tempParticle));
        }
        // Set Center of Mass to Zero
        sumx /= NParticles;
        sumy /= NParticles;
        for(int i = 0; i < NParticles; ++i) {
            Particles[i]->x -= sumx;
            Particles[i]->y -= sumy;
        }
        InteractionRadius = 1;
        InitCells(InteractionRadius);
        // Compute System Properties
        CalculateMinMeanDistance(true);
        CalculateTemperature(true);
        CalculateKineticEnergy(true);
        // Print Container Properties
        if(PrintInitInfo){
            std::cout << "  Container Properties: " << std::endl;
            std::cout << "      L = " << L << std::endl;
            std::cout << "      Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << std::endl;
            std::cout << "      Area = " << Area << std::endl;
            std::cout << "      Volume = " << Volume << std::endl;
            std::cout << "      Surface = " << Surface << std::endl;
            // Print Cells Properties
            std::cout << "  Cells Properties: " << std::endl;
            std::cout << "      Ncellsx = " << Ncellsx << std::endl;
            std::cout << "      Ncellsy = " << Ncellsy << std::endl;
            std::cout << "      Ncellsz = " << Ncellsz << std::endl;
            // Print System Properties
            std::cout << "  System Properties: " << std::endl;
            std::cout << "      Temperature = " << Temperature << std::endl;
            std::cout << "      Kinetic Energy = " << KineticEnergy << std::endl;
            std::cout << "      PotentialFunction Energy = " << PotentialEnergy << std::endl;
            std::cout << "      Total Energy = " << TotalEnergy << std::endl;
            std::cout << "      Packing Fraction = " << PackingFraction << std::endl;
            std::cout << "      Density = " << Density << std::endl;
            std::cout << "      MinMeanDistance = " << MinMeanDistance << std::endl;
        }

    }
    void InitRandom(double Mass, double Radius) {
        std::cout << "Initialize Random Particles with distributed Radius and Mass" << std::endl;
        std::cout << "and distributed positions and velocities" << std::endl;
        Particles.clear();

        if (NParticles == 0){
            std::cout << "  No NParticles given, set to 100" << std::endl;
            NParticles = 100;
        }
        else{
            std::cout << "  with " << NParticles << " particles" << std::endl;
        }
        if (DensityInit == 0){
            std::cout << "  No DensityInit given, set to 1" << std::endl;
            DensityInit = 1;
            Density = DensityInit;
            L = pow(NParticles / DensityInit, 1.0/3.0);
            Lx = Ly = Lz = L;
            Volume = L*L*L;
            Surface = L*L*6;
            Area = nan("Not 2D");
            Dim = 3;
        }
        else{
            std::cout << "  with Density " << DensityInit << std::endl;
            Density = DensityInit;
            L = pow(NParticles / DensityInit, 1.0/3.0);
            Lx = Ly = Lz = L;
            Volume = L*L*L;
            Surface = L*L*6;
            Area = nan("Not 2D");
            Dim = 3;
        }

        // Initialize Particles
        double sumx = 0, sumy = 0, sumz = 0;
        for(int i = 0; i < NParticles; ++i) {
            double Rand = ((double) rand() / RAND_MAX - 0.5); // Random Number between -0.5 and 0.5
            Particle tempParticle;
            tempParticle.myindex = i;
            tempParticle.myradius = Radius * (1 + Rand);
            tempParticle.mymass = Mass * (1 + Rand);
            tempParticle.x = (double) rand() / RAND_MAX * L;
            tempParticle.y = (double) rand() / RAND_MAX * L;
            tempParticle.z = (double) rand() / RAND_MAX * L;
            sumx += tempParticle.x;
            sumy += tempParticle.y;
            sumz += tempParticle.z;
            Particles.push_back(std::make_shared<Particle>(tempParticle));
        }
        // Set Center of Mass to Zero
        sumx /= NParticles;
        sumy /= NParticles;
        sumz /= NParticles;
        for(int i = 0; i < NParticles; ++i) {
            Particles[i]->x -= sumx;
            Particles[i]->y -= sumy;
            Particles[i]->z -= sumz;
        }
        InitCells(InteractionRadius);
        InitGausianDistributedVelocities();
        // Compute System Properties
        CalculateMinMeanDistance(true);
        CalculateTemperature(true);
        CalculateKineticEnergy(true);

        if(PrintInitInfo){
            std::cout << "  Container Properties: " << std::endl;
            std::cout << "      L = " << L << std::endl;
            std::cout << "      Lx = " << Lx << ", Ly = " << Ly << ", Lz = " << Lz << std::endl;
            std::cout << "      Area = " << Area << std::endl;
            std::cout << "      Volume = " << Volume << std::endl;
            std::cout << "      Surface = " << Surface << std::endl;
            // Print Cells Properties
            std::cout << "  Cells Properties: " << std::endl;
            std::cout << "      Ncellsx = " << Ncellsx << std::endl;
            std::cout << "      Ncellsy = " << Ncellsy << std::endl;
            std::cout << "      Ncellsz = " << Ncellsz << std::endl;
            // Print System Properties
            std::cout << "  System Properties: " << std::endl;
            std::cout << "      Temperature = " << Temperature << std::endl;
            std::cout << "      Kinetic Energy = " << KineticEnergy << std::endl;
            std::cout << "      PotentialFunction Energy = " << PotentialEnergy << std::endl;
            std::cout << "      Total Energy = " << TotalEnergy << std::endl;
            std::cout << "      Packing Fraction = " << PackingFraction << std::endl;
            std::cout << "      Density = " << Density << std::endl;
            std::cout << "      MinMeanDistance = " << MinMeanDistance << std::endl;
        }

    }
    void InitLennardJonesFluid(double RadiusInit, double MassInit){


    }
    //----------------------------------------------------------------------------------------------
    // Init Functions System
    void InitRandomDistributedVelocities(){
        double sumx = 0, sumy = 0, sumz = 0;
        double T0, factor;
        if (TemperatureInit == 0){
            std::cout << "  No TemperatureInit given, set to 1" << std::endl;
            TemperatureInit = 1;
        }
        switch (Dim) {
            case 3:
                for(const auto &Current : Particles){
                    Current->vx = (double) rand() / RAND_MAX  - 0.5;
                    Current->vy = (double) rand() / RAND_MAX  - 0.5;
                    Current->vz = (double) rand() / RAND_MAX  - 0.5;
                    sumx += Current->vx;
                    sumy += Current->vy;
                    sumz += Current->vz;
                }
                // Shift so that the total momentum is zero
                sumx /= NParticles;
                sumy /= NParticles;
                sumz /= NParticles;
                for(const auto &Current : Particles){
                    Current->vx -= sumx;
                    Current->vy -= sumy;
                    Current->vz -= sumz;
                }
                // Rescale velocities to the desired TemperatureInit
                T0 = CalculateTemperature(false);
                factor = sqrt(TemperatureInit / T0);
                for(const auto &Current : Particles){
                    Current->vx *= factor;
                    Current->vy *= factor;
                    Current->vz *= factor;
                }
                T0 = CalculateTemperature(true);
                break;
            case 2:
                for(const auto &Current : Particles){
                    Current->vx = (double) rand() / RAND_MAX  - 0.5;
                    Current->vy = (double) rand() / RAND_MAX  - 0.5;
                    sumx += Current->vx;
                    sumy += Current->vy;
                }
                // Shift so that the total momentum is zero
                sumx /= NParticles;
                sumy /= NParticles;
                for(const auto &Current : Particles){
                    Current->vx -= sumx;
                    Current->vy -= sumy;
                }
                // Rescale velocities to the desired TemperatureInit
                T0 = CalculateTemperature(false);
                factor = sqrt(TemperatureInit / T0);
                for(const auto &Current : Particles){
                    Current->vx *= factor;
                    Current->vy *= factor;
                }
                T0 = CalculateTemperature(true);
                break;
        }
    }
    void InitGausianDistributedVelocities(){
        // calculate mean velocity from TemperatureInit
        if(TemperatureInit == 0){
            std::cout << "  No TemperatureInit given, set to 1" << std::endl;
            TemperatureInit = 1;
        }
        double v_mean = sqrt(2 * TemperatureInit);
        // calculate standard deviation from TemperatureInit
        double sigma = sqrt(TemperatureInit);
        // Draw random velocities with Box-Muller-Algorithm
        double u1, u2, u3, v1, v2, v3, s;
        double sumx = 0, sumy = 0, sumz = 0;
        switch (Dim) {
            case 2:
                for(const auto &Current : Particles){
                    do {
                        u1 = (double) rand() / RAND_MAX;
                        u2 = (double) rand() / RAND_MAX;
                        v1 = 2 * u1 - 1;
                        v2 = 2 * u2 - 1;
                        s = v1 * v1 + v2 * v2;
                    } while (s >= 1 || s == 0);
                    Current->vx = v_mean * v1 * sqrt(-2 * log(s) / s);
                    Current->vy = v_mean * v2 * sqrt(-2 * log(s) / s);
                    sumx += Current->vx;
                    sumy += Current->vy;
                }
                // Shift so that the total momentum is zero
                sumx /= NParticles;
                sumy /= NParticles;
                for(const auto &Current : Particles){
                    Current->vx -= sumx;
                    Current->vy -= sumy;
                }

                break;
            case 3:
                for(const auto &Current : Particles){
                    do {
                        u1 = (double) rand() / RAND_MAX;
                        u2 = (double) rand() / RAND_MAX;
                        u3 = (double) rand() / RAND_MAX;
                        v1 = 2 * u1 - 1;
                        v2 = 2 * u2 - 1;
                        v3 = 2 * u3 - 1;
                        s = v1 * v1 + v2 * v2 + v3 * v3;
                    } while (s >= 1 || s == 0);
                    Current->vx = v_mean * v1 * sqrt(-2 * log(s) / s);
                    Current->vy = v_mean * v2 * sqrt(-2 * log(s) / s);
                    Current->vz = v_mean * v3 * sqrt(-2 * log(s) / s);
                    sumx += Current->vx;
                    sumy += Current->vy;
                    sumz += Current->vz;
                }
                // Shift so that the total momentum is zero
                sumx /= NParticles;
                sumy /= NParticles;
                sumz /= NParticles;
                for(const auto &Current : Particles){
                    Current->vx -= sumx;
                    Current->vy -= sumy;
                    Current->vz -= sumz;
                }
                break;
        }
        // Rescale velocities to the desired TemperatureInit
        double T0 = CalculateTemperature(false);
        double factor = sqrt(TemperatureInit / T0);
        for(const auto &Current : Particles){
            Current->vx *= factor;
            Current->vy *= factor;
            Current->vz *= factor;
        }
        T0 = CalculateTemperature(true);
        KineticEnergy = CalculateKineticEnergy(true);
        TotalEnergy = KineticEnergy + PotentialEnergy;

    }
    //----------------------------------------------------------------------------------------------
    // Calculate and Set Function
    double CalculateTemperature(bool Set)
    {
        double sum = 0.0, T0;
        for (const auto &particle : Particles) {
            sum += particle->vx * particle->vx + particle->vy * particle->vy + particle->vz * particle->vz;
        }
        T0 = sum / (3.0 * NParticles);
        if(Set)  Temperature = T0;
        return T0;
    }
    double CalculateKineticEnergy(bool Set){
        double sum = 0.0;
        for (const auto &particle : Particles) {
            sum += particle->vx * particle->vx + particle->vy * particle->vy + particle->vz * particle->vz;
        }
        if(Set)  {
            KineticEnergy = 0.5 * sum;
            TotalEnergy = KineticEnergy + PotentialEnergy;
        }
        return 0.5 * sum;
    }
    double CalculatePackingFraction(bool Set){
        double sum = 0.0;
        for (const auto &particle : Particles) {
            sum += 4.0/3 * M_PI * pow(particle->myradius, 3);
        }
        if(Set)  {
            PackingFraction = sum / Volume;
        }
        return sum / Volume;
    }
    double CalculateMinMeanDistance(bool Set)
    {
        double mean = 0.0, dr, min;

        for(const auto &Current : Particles){
            min = L * 2;
            for(const auto &Other : Particles){
                if(Current != Other){
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    if(dr < min) min = dr;
                }
            }
            mean += min;
        }
        if(Set) MinMeanDistance = mean / NParticles;
        return mean / NParticles;
    }
        // All using Cell Lists
    void CalculatePotentialEnergy(){
        PotentialEnergy = 0;
        for(const auto &Current : Particles){
            for(const auto &Other : Current->mycell->ParticlesInCell){
                if(Current != Other  && Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz) < InteractionRadius){
                    PotentialEnergy += PotentialFunction(Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz));
                }
            }
            for(const auto &NCell : Current->mycell->NeighbourCells){
                for(const auto &Other : NCell->ParticlesInCell){
                    if(Current != Other  && Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz) < InteractionRadius){
                        PotentialEnergy += PotentialFunction(Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz));
                    }
                }
            }
        }
        PotentialEnergy = PotentialEnergy / 2;
        TotalEnergy = KineticEnergy + PotentialEnergy;
    }
    void CalculatePotentialEnergyGlobal(){
        double Potential = 0;
        for(const auto &Current : Particles){
            for(const auto &Other : Particles){
                if(Current != Other  &&
                    Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz) < InteractionRadius)
                {
                    Potential += PotentialFunction(Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz));
                }
            }
        }
        PotentialEnergy = Potential / 2;
    }
    void CalculateForces() {
        double dr, F;
        // Reset Forces
        for (const auto &Current: Particles) {
            Current->fx = 0;
            Current->fy = 0;
            Current->fz = 0;
        }
        // Calculate Forces
        for (const auto &Current: Particles) {
            for (const auto &Other: Current->mycell->ParticlesInCell) {
                if (Current != Other  &&
                    Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz) < InteractionRadius)
                {
                    F = ForceFunction(Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz));
                    Current->fx += F * (Current->x - Other->x);
                    Current->fy += F * (Current->y - Other->y);
                    Current->fz += F * (Current->z - Other->z);
                }
            }
            for (const auto &NCell: Current->mycell->NeighbourCells) {
                for (const auto &Other: NCell->ParticlesInCell) {
                    if (Current != Other  &&
                        Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz) < InteractionRadius)
                    {
                        F = ForceFunction(Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz));
                        Current->fx += F * (Current->x - Other->x);
                        Current->fy += F * (Current->y - Other->y);
                        Current->fz += F * (Current->z - Other->z);
                    }
                }
            }
        }
    }
    void CalculateForcesGlobally(){
        double dr, Fx, Fy, Fz;
        // Reset Forces
        for (const auto &Current: Particles) {
            Current->fx = 0;
            Current->fy = 0;
            Current->fz = 0;
        }
        // Calculate Forces
        for (const auto &Current: Particles) {
            for (const auto &Other: Particles) {
                if (Current != Other) {
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    Fx = ForceFunction(dr) * (Current->x - Other->x) / dr;
                    Fy = ForceFunction(dr) * (Current->y - Other->y) / dr;
                    Fz = ForceFunction(dr) * (Current->z - Other->z) / dr;
                    Current->fx += Fx;
                    Current->fy += Fy;
                    Current->fz += Fz;
                }
            }
        }
    }

    //----------------------------------------------------------------------------------------------
    // Calculate and Export Functions
    void CalculateExportRDF(std::string Filename, bool show, bool makeplot) {
        double Phi0 = CalculatePackingFraction(false);
        int Nhisto = 200, bin; // Number of Histogram Bins and bin for placing particles
        double drhisto = (L / Nhisto) / 2; // Histogram Bin Size
        double dr; // Distance between two particles
        std::vector<double> Histogram(Nhisto, 0); // Histogram Vector
        for(const auto &Current : Particles){
            for(const auto &Other : Particles){
                if(Current != Other){
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    bin = (int) (dr / drhisto + 0.5);
                    if(bin < Nhisto) Histogram[bin] += 1.0;
                }
            }
        }
        // Normalize Histogram
        double Volume;
        for(int i = 0; i < Nhisto; ++i){
            Volume = 4.0/3 * M_PI * (pow((i+1)*drhisto, 3) - pow(i*drhisto, 3));
            Histogram[i] = Histogram[i] / (Volume * NParticles * Phi0 * 6/M_PI);
        }
        // Export Histogram
        std::ofstream OutputFile;
        OutputFile.open(Filename);
        for(int i = 0; i < Nhisto; ++i){
            OutputFile << i*drhisto << " " << Histogram[i] << std::endl;
        }
        OutputFile.close();
        // Plot Radial Distribution Function
        if(makeplot)
        {
            if (show) {
                std::string command = "python code/visualize.py PlotRDF 1 " + Filename;
                system(command.c_str());
            }
            else {
                std::string command = "python code/visualize.py PlotRDF 0 " + Filename;
                system(command.c_str());
            }
        }
    }
    void CalculateExportSSF(std::string Filename, bool show, bool makeplot){
    }




    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
};


#endif //SIMULATION_SYSTEM_H
