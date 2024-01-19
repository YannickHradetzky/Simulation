//
// Created by Yannick Hradetzky on 11.12.23.
//

#ifndef SIMULATION_BOX_H
#define SIMULATION_BOX_H

#include "../code/Particle.h"
#include "Celllist.h"

using namespace std;
// Potential Functions
using PotentialFunctionPosition = function<double(double, double , double )>;
using PotentialFunctionDistance = function<double(double)>;
using ForceFunctionDistance = function<double(double)>;

class Box {
public:
    Box();
    // System Variables
    double Lx, Ly, Lz, L;
    double V, AbsolutAverageVelocity;
    double Rho, Phi;
    double Ekin, Epot, Etot;
    double MeanMinDistance;
    int N;
    double T;
    int Dim;
    bool Overlap;
    // Center of Mass Variables
    double CMVx, CMVy, CMVz;
    double CMx, CMy, CMz;
    // Particles and Cells
    vector<shared_ptr<Particle>> Particles;
    CellList *CellListBox;
    // System Funtions
    vector<double> BoxRDFFunction;
    vector<double> BoxSSFFunction;

    // Setter Functions
    void SetBoxDimensions(double Lx0, double Ly0, double Lz0);
    void SetBoxDimensions(double Lx0, double Ly0);
    void SetCubicBoxDimensions(double L0);
    void SetQuadraticBoxDimensions(double L0);
    void SetTemperature(double T0);
    void SetDensity(double Rho0);
    void SetCenterOfMassToZero();
    void SetCenterOfMassVelocityToZero();
    void SetParticleRadius(double Radius0);
    void SetCellList(CellList Cell_List0);


    // Initialization Functions
    void InitializeRandomParticles(int N0);
    void InitializeLatticeParticlesPhi(int N0, double Radius, double PackingFraction);
    void InitializeRandomVelocities();
    void InitializeVicsekParticles(int N0, double Radius, double velcity0);


    // Calculation Functions
    [[nodiscard]] double CalculateTemperature() const;
    [[nodiscard]] double CalculateKineticEnergy() const;
    [[nodiscard]] double CalculateDensity () const;
    [[nodiscard]] double CalculatePackingFraction() const;
    [[nodiscard]] double CalculateVolume() const;
    [[nodiscard]] double* CalculateCenterOfMass() const;
    [[nodiscard]] double* CalculateCenterOfMassVelocity() const;
    [[nodiscard]] vector<double> CalculateSSF(int axis) const;
    [[nodiscard]] double CalculatePotentialEnergy(PotentialFunctionPosition PotentialFunction);
    [[nodiscard]] double CalculatePotentialEnergy(PotentialFunctionDistance PotentialFunction);
    [[nodiscard]] double CalculateMeanMinDistance() const;
    [[nodiscard]] double CalculateAbsolutAverageVelocity() const;


    // Update Functions
    void UpdateKineticEnergy();
    void UpdatePotentialEnergy(PotentialFunctionPosition PotentialFunction);
    void UpdatePotentialEnergy(PotentialFunctionDistance PotentialFunction);
    void UpdateForces(ForceFunctionDistance ForceFunction);
    void UpdateTotalEnergy();
    void UpdateAbsoluteAverageVelocity();
    void UpdateMeanMinDistance();
    void UpdateOverlap();

    // Export Functions
    void ExportDisplayParticlePositions(string Filename, bool show) const;
    void ExportDisplayParticleVelocities(string Filename, bool show) const;
    void ExportCalculateRDF(string Filename, bool show) const;


    // Print Functions
    void PrintBoxInfo() const;
    void PrintParticleInfo(int i) const;
    void PrintParticleInfo() const;
    void PrintCellsOfParticle(int i) const;
    void PrintCellOfParticle() const;



private:
    shared_ptr<CellList> BoxCellList;
};


#endif //SIMULATION_BOX_H
