//
// Created by Yannick Hradetzky on 16.12.23.
//

#include "OldSystem.h"

void System::InitLattice(int N, double Phi, double Radius, int NCells) {
    // Initialize a Lattice with N particles and a packing fraction Phi with no Velocity
    // The particles have a radius of Radius
    // The lattice is initialized in a cubic box with periodic boundary conditions
    // Initialize the Cell List with 10 cells in each direction

    Container->InitializeLatticeParticlesPhi(N, Radius, Phi);
    Container->V = Container->CalculateVolume();
    Container->Rho = Container->CalculateDensity();
    Container->Phi = Container->CalculatePackingFraction();
    Container->T = Container->CalculateTemperature();
    Container->SetCenterOfMassToZero();

    CellListSystem->Nneighbours = 26;
    CellListSystem->LCutOffx = Container->Lx / NCells;
    CellListSystem->LCutOffy = Container->Ly / NCells;
    CellListSystem->LCutOffz = Container->Lz / NCells;
    CellListSystem->DimCellList = Container->Dim;
    CellListSystem->Ncellsx = NCells;
    CellListSystem->Ncellsy = NCells;
    CellListSystem->Ncellsz = NCells;

    CellListSystem->SetParticles(Container->Particles);
    CellListSystem->InitializeCells();



    Container->SetCellList(*CellListSystem);

}
//
void System::InitRandom2D(int N, double Rho, double Radius, int NCells) {
    // Initialize a random system with N particles in a box with density Rho
    // The particles have a radius of Radius
    double L0 = sqrt(N / Rho);
    Container->SetQuadraticBoxDimensions(L0);
    Container->InitializeRandomParticles(N);
    Container->SetCenterOfMassToZero();
    Container->SetParticleRadius(Radius);

    CellListSystem->Nneighbours = 8;
    CellListSystem->LCutOffx = Container->Lx / NCells;
    CellListSystem->LCutOffy = Container->Ly / NCells;
    CellListSystem->DimCellList = Container->Dim;
    CellListSystem->Ncellsx = NCells;
    CellListSystem->Ncellsy = NCells;

    CellListSystem->SetParticles(Container->Particles);
    CellListSystem->InitializeCells();

    Container->SetCellList(*CellListSystem);

    SystemParticles = Container->Particles;
    SystemCells = CellListSystem->Cells;

}
//
void System::InitRandom3D(int N, double Rho, double Radius, int NCells) {
    // Initialize a random system with N particles in a box with density Rho
    // The particles have a radius of Radius
    double L0 = pow(N / Rho, 1.0 / 3.0);
    Container->SetCubicBoxDimensions(L0);
    Container->InitializeRandomParticles(N);
    Container->SetCenterOfMassToZero();
    Container->SetParticleRadius(Radius);

    CellListSystem->Nneighbours = 26;
    CellListSystem->LCutOffx = Container->Lx / NCells;
    CellListSystem->LCutOffy = Container->Ly / NCells;
    CellListSystem->LCutOffz = Container->Lz / NCells;
    CellListSystem->DimCellList = Container->Dim;
    CellListSystem->Ncellsx = NCells;
    CellListSystem->Ncellsy = NCells;
    CellListSystem->Ncellsz = NCells;

    CellListSystem->SetParticles(Container->Particles);
    CellListSystem->InitializeCells();

    Container->SetCellList(*CellListSystem);

    SystemParticles = Container->Particles;
    SystemCells = CellListSystem->Cells;
}

void System::InitVicsek(int N, double Radius, double Velocity, double Rho, int NCells) {
    // Initialize a Vicsek system with N particles in a box with density Rho
    // The particles have a radius of Radius

    double L0 = sqrt(N/Rho);

    Container->SetQuadraticBoxDimensions(L0);
    Container->InitializeVicsekParticles(N, Radius, Velocity);
    Container->SetCenterOfMassVelocityToZero();
    Container->SetCenterOfMassToZero();
    Container->SetParticleRadius(Radius);

    CellListSystem->Nneighbours = 8;
    CellListSystem->LCutOffx = Container->Lx / NCells;
    CellListSystem->LCutOffy = Container->Ly / NCells;
    CellListSystem->DimCellList = Container->Dim;
    CellListSystem->Ncellsx = NCells;
    CellListSystem->Ncellsy = NCells;

    CellListSystem->SetParticles(Container->Particles);
    CellListSystem->InitializeCells();

    Container->SetCellList(*CellListSystem);

    Container->Ekin = Container->CalculateKineticEnergy();
    Container->Etot = Container->Ekin + Container->Epot;
    Container->T = Container->CalculateTemperature();
    Container->V = Container->CalculateVolume();
    Container->Rho = Container->CalculateDensity();


    SystemParticles = Container->Particles;
    SystemCells = CellListSystem->Cells;

}

void System::InitHardSpheres(int N, double Radius, double Phi) {
    Container->InitializeLatticeParticlesPhi(N, Radius, Phi);
    Container->V = Container->CalculateVolume();
    Container->Rho = Container->CalculateDensity();
    Container->Phi = Container->CalculatePackingFraction();
    Container->T = Container->CalculateTemperature();
    Container->UpdateMeanMinDistance();

    int NCells0 = floor(Container->L / (3 * Radius));
    if (NCells0 < 3) NCells0 = 3;

    CellListSystem->Nneighbours = 26;
    CellListSystem->LCutOffx = Container->Lx / NCells0;
    CellListSystem->LCutOffy = Container->Ly / NCells0;
    CellListSystem->LCutOffz = Container->Lz / NCells0;
    CellListSystem->DimCellList = Container->Dim;
    CellListSystem->Ncellsx = CellListSystem->Ncellsy = CellListSystem->Ncellsz = NCells0;
    CellListSystem->SetParticles(Container->Particles);
    CellListSystem->InitializeCells();

    Container->SetCellList(*CellListSystem);
    SystemParticles = Container->Particles;
    SystemCells = CellListSystem->Cells;
}



