//
// Created by Yannick Hradetzky on 11.12.23.
//

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "Box.h"
#include "../code/Particle.h"
#include "Celllist.h"
#include "vector"
#include "fstream"


using namespace std;
//----------------------------------------------------------------------
// Constructors

Box::Box() {
    Box::Lx = 0.0;
    Box::Ly = 0.0;
    Box::Lz = 0.0;
    Box::L = 0.0;
    Box::Dim = 0;
    Box::N = 0;
    Box::Rho = 0.0;
    Box::T = 0.0;
    Box::V = 0.0;
    Box::Ekin = 0.0;
    Box::Epot = 0.0;
    Box::Etot = 0.0;
}

//----------------------------------------------------------------------
// Setter Functions
void Box::SetBoxDimensions(double Lx0, double Ly0, double Lz0) {
    Lx = Lx0;
    Ly = Ly0;
    Lz = Lz0;
    Dim = 3;
    V = Lx0 * Ly0 * Lz0;
}
//
void Box::SetBoxDimensions(double Lx0, double Ly0) {
    Lx = Lx0;
    Ly = Ly0;
    Lz = 0.0;
    Dim = 2;
    V = Lx0 * Ly0;
}
//
void Box::SetCubicBoxDimensions(double L0) {
    Lx = L0;
    Ly = L0;
    Lz = L0;
    L = L0;
    Dim = 3;
    V = L0 * L0 * L0;
}
//
void Box::SetQuadraticBoxDimensions(double L0) {
    Lx = L0;
    Ly = L0;
    L = L0;
    Lz = 0.0;
    Dim = 2;
    V = Ly * Lx;
}
//
void Box::SetTemperature(double T0) {
    double Tcurrent = CalculateTemperature();
    double factor = sqrt(T0 / Tcurrent);
    for (int i = 0; i < N; ++i) {
        Particles[i]->vx *= factor;
        Particles[i]->vy *= factor;
        Particles[i]->vz *= factor;
    }
    T = CalculateTemperature();
}
//
void Box::SetDensity(double Rho0) {
    double L0 = L;
    double factor;
    switch (Dim){
        case 2:
            L = sqrt(N / Rho0);
            Lx = L;
            Ly = L;
            V = Lx * Ly;
            factor = L / L0;
            for (int i = 0; i < N; ++i) {
                Particles[i]->x *= factor;
                Particles[i]->y *= factor;
            }
            break;
        case 3:
            L = pow(N / Rho0, 1. / 3.);
            Lx = L;
            Ly = L;
            Lz = L;
            V = Lx * Ly * Lz;
            factor = L / L0;
            for (int i = 0; i < N; ++i) {
                Particles[i]->x *= factor;
                Particles[i]->y *= factor;
                Particles[i]->z *= factor;
            }
            break;
    }
    Rho = CalculateDensity();

}
//
void Box::SetCenterOfMassToZero() {
    double * CM = CalculateCenterOfMass();
    for (int i = 0; i < N; ++i) {
        if (Dim == 3) {
            Particles[i]->x -= CM[0];
            Particles[i]->y -= CM[1];
            Particles[i]->z -= CM[2];
        } else if (Dim == 2) {
            Particles[i]->x -= CM[0];
            Particles[i]->y -= CM[1];
        } else {
            cout << "Error: Dimension is neither 2 nor 3!" << endl;
        }
    }
    CM = CalculateCenterOfMass();
    if (Dim == 3) {
        CMx = CM[0];
        CMy = CM[1];
        CMz = CM[2];
    } else if (Dim == 2) {
        CMx = CM[0];
        CMy = CM[1];
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
    }
}
//
void Box::SetCenterOfMassVelocityToZero() {
    double * CMV = CalculateCenterOfMassVelocity();
    for (int i = 0; i < N; ++i) {
        if (Dim == 3) {
            Particles[i]->vx -= CMV[0];
            Particles[i]->vy -= CMV[1];
            Particles[i]->vz -= CMV[2];
        } else if (Dim == 2) {
            Particles[i]->vx -= CMV[0];
            Particles[i]->vy -= CMV[1];
        } else {
            cout << "Error: Dimension is neither 2 nor 3!" << endl;
        }
    }
    CMV = CalculateCenterOfMassVelocity();
    if (Dim == 3) {
        CMVx = CMV[0];
        CMVy = CMV[1];
        CMVz = CMV[2];
    } else if (Dim == 2) {
        CMVx = CMV[0];
        CMVy = CMV[1];
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
    }
    SetTemperature(T);
}
//
void Box::SetParticleRadius(double Radius0) {
    for (int i = 0; i < N; ++i) {
        Particles[i]->myradius = Radius0;
    }
}
//
void Box::SetCellList(CellList Cell_List0) {
    CellListBox = &Cell_List0;
    // TODO Think about if this is efficient
}


//----------------------------------------------------------------------
// Initialization Functions
void Box::InitializeRandomParticles(int N0) {
    N = N0;
    Particles.clear();
    for (int i = 0; i < N; ++i) {
        Particle n;
        n.myindex = i;
        n.x = Lx * (double) rand() / (double) RAND_MAX;
        n.y = Ly * (double) rand() / (double) RAND_MAX;
        n.z = Lz * (double) rand() / (double) RAND_MAX;
        Particles.push_back(make_shared<Particle>(n));

    }
    Rho = CalculateDensity();
    double *CM = CalculateCenterOfMass();
    if (Dim == 3) {
        CMx = CM[0];
        CMy = CM[1];
        CMz = CM[2];
    } else if (Dim == 2) {
        CMx = CM[0];
        CMy = CM[1];
        CMz = 0.0;
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
    }
}
//
void Box::InitializeLatticeParticlesPhi(int N0, double Radius, double PackingFraction) {
    Dim = 3;
    N = N0;
    Particles.clear();
    int ParticlesPerDimension = (int) pow(N + 0.00002,1./3); // Number of particles per dimension
    N = ParticlesPerDimension * ParticlesPerDimension * ParticlesPerDimension; // Total Number of particles
    if (N != N0){
        cout << "Warning: Number of Particles is not a perfect cube!" << endl;
    }

    double L0 = ParticlesPerDimension * 2 * Radius; // Length of the box with Atoms touching each other
    double LatticeSpacing = 2 * Radius; // Lattice Spacing
    double MyX, MyY, MyZ; // Coordinates of the current particle
    double CMX, CMY, CMZ; // Center of Mass of the current Configuration
    int i, j, k, u; // Counters for the loops
    int Index; // Index of the current particle
    CMX = CMY = CMZ = 0.0; // Initialize the Center of Mass

    // Compute the Length of the Box with the given packing fraction
    double VParticle = 4. / 3. * M_PI * Radius * Radius * Radius;
    if (Radius == 0.0){
        cout << "Warning: Radius is 0.0!" << endl;
        cout << "Initialize Particles Spacing 1.0!" << endl;
        LatticeSpacing = 1.0;
        L0 = ParticlesPerDimension * LatticeSpacing;
        L = ParticlesPerDimension * LatticeSpacing;
        Lx = Ly = Lz = L;
        V = CalculateVolume();

    }else if (Radius > 0.0){
        L0 = ParticlesPerDimension * 2 * Radius; // Length of the box with Atoms touching each other
        L = pow(N * VParticle / PackingFraction, 1. / 3.);
        Lx = Ly = Lz = L;
        V = CalculateVolume();
    }
    if (L < L0)
    {
        printf("The given packing fraction is too large!\n");
        printf("The maximum packing fraction is %f\n", N * VParticle / (L0 * L0 * L0));
        exit(1);
    }else{
        Overlap = false;
    }

    // Create the particles and add them to the Box
    for (i = 0; i < ParticlesPerDimension; ++i)
    {
        MyX = i * LatticeSpacing;
        for (j = 0; j < ParticlesPerDimension; ++j)
        {
            MyY = j * LatticeSpacing;
            for (k = 0; k < ParticlesPerDimension; ++k)
            {
                MyZ = k * LatticeSpacing;
                Index = i * ParticlesPerDimension * ParticlesPerDimension + j * ParticlesPerDimension + k;
                Particle n;
                CMX += MyX;
                CMY += MyY;
                CMZ += MyZ;
                n.myindex = Index;
                n.x = MyX;
                n.y = MyY;
                n.z = MyZ;
                n.myradius = Radius;
                Particles.push_back(make_shared<Particle>(n));
            }
        }
    }
    CMX /= N;
    CMY /= N;
    CMZ /= N;

    // Center the Box around the Center of Mass
    for (i = 0; i < N; ++i)
    {
        Particles[i]->x -= CMX;
        Particles[i]->y -= CMY;
        Particles[i]->z -= CMZ;
    }
    // Rescale the Box to the desired length
    for (i = 0; i < N; ++i)
    {
        Particles[i]->x *= L / L0;
        Particles[i]->y *= L / L0;
        Particles[i]->z *= L / L0;
    }
    // Calculate the Center of Mass of the current Configuration
    double* CM = CalculateCenterOfMass();
    CMx = CM[0];
    CMy = CM[1];
    CMz = CM[2];
}
void Box::InitializeRandomVelocities() {
    double xsum = 0, ysum = 0, zsum = 0;
    for (int i = 0; i < N; ++i) {
        Particles[i]->vx = (double) rand() / (double) RAND_MAX - 0.5;
        Particles[i]->vy = (double) rand() / (double) RAND_MAX - 0.5;
        Particles[i]->vz = (double) rand() / (double) RAND_MAX - 0.5;
        xsum += Particles[i]->vx;
        ysum += Particles[i]->vy;
        zsum += Particles[i]->vz;
    }
    // shift velocities so total momentum is zero
    for (int i = 0; i < N; ++i) {
        Particles[i]->vx -= xsum / (double) N;
        Particles[i]->vy -= ysum / (double) N;
        Particles[i]->vz -= zsum / (double) N;
    }
    //calculate current temperature
    T = CalculateTemperature();
    Ekin = CalculateKineticEnergy();
    Etot = Ekin + Epot;
}
//
void Box::InitializeVicsekParticles(int N0, double Radius, double velcity0) {
    N = N0;
    Particles.clear();
    Dim = 2;
    double CMX = 0, CMY = 0;
    for (int i = 0; i < N; ++i) {
        Particle n;
        n.myindex = i;
        n.myradius = Radius;
        n.theta = 2.0 * M_PI * (double) rand() / (double) RAND_MAX;
        n.x = Lx * (double) rand() / (double) RAND_MAX ;
        n.y = Ly * (double) rand() / (double) RAND_MAX;
        n.z = 0;
        CMX += n.x;
        CMY += n.y;
        n.vx = velcity0 * cos(n.theta);
        n.vy = velcity0 * sin(n.theta);
        n.vz = 0.0;

        Particles.push_back(make_shared<Particle>(n));
    }

    CMX /= N;
    CMY /= N;
    // Center the Box around the Center of Mass
    for (int i = 0; i < N; ++i) {
        Particles[i]->x -= CMX;
        Particles[i]->y -= CMY;
    }

    Rho = CalculateDensity();
    AbsolutAverageVelocity = CalculateAbsolutAverageVelocity();
    double *CM = CalculateCenterOfMass();
    if (Dim == 3) {
        CMx = CM[0];
        CMy = CM[1];
        CMz = CM[2];
    } else if (Dim == 2) {
        CMx = CM[0];
        CMy = CM[1];
        CMz = 0.0;
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
    }
    CM = CalculateCenterOfMassVelocity();
    if (Dim == 3) {
        CMVx = CM[0];
        CMVy = CM[1];
        CMVz = CM[2];
    } else if (Dim == 2) {
        CMVx = CM[0];
        CMVy = CM[1];
        CMVz = 0.0;
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
    }
    T = CalculateTemperature();
}

//----------------------------------------------------------------------
// Calculation Functions (can be const)
double Box::CalculateTemperature() const {
    double sum = 0.0, T0;
    for (int i = 0; i < N; ++i) {
        sum += Particles[i]->vx * Particles[i]->vx + Particles[i]->vy * Particles[i]->vy + Particles[i]->vz * Particles[i]->vz;
    }
    T0 = sum / (3.0 * N);
    return T0;
}
//
double Box::CalculateKineticEnergy() const {
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += Particles[i]->vx * Particles[i]->vx + Particles[i]->vy * Particles[i]->vy + Particles[i]->vz * Particles[i]->vz;
    }
    return 0.5 * sum;
}
//
double Box::CalculateDensity() const {
    return N / V;
}
//
double Box::CalculateVolume() const {
    double V0;
    switch (Dim){
        case 2:
            V0 = Lx * Ly;
            break;
        case 3:
            V0 = Lx * Ly * Lz;
            break;
    }
    return V0;
}
//
double *Box::CalculateCenterOfMass() const {
    if (Dim == 3){
        double *CM = new double[3];
        CM[0] = 0.0;
        CM[1] = 0.0;
        CM[2] = 0.0;
        for (int i = 0; i < N; ++i) {
            CM[0] += Particles[i]->x;
            CM[1] += Particles[i]->y;
            CM[2] += Particles[i]->z;
        }
        CM[0] /= N;
        CM[1] /= N;
        CM[2] /= N;
        return CM;
    } else if (Dim == 2){
        double *CM = new double[2];
        CM[0] = 0.0;
        CM[1] = 0.0;
        for (int i = 0; i < N; ++i) {
            CM[0] += Particles[i]->x;
            CM[1] += Particles[i]->y;
        }
        CM[0] /= N;
        CM[1] /= N;
        return CM;
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
        return nullptr;
    }
}
//
double *Box::CalculateCenterOfMassVelocity() const {
    if (Dim == 3){
        double *CMV = new double[3];
        CMV[0] = 0.0;
        CMV[1] = 0.0;
        CMV[2] = 0.0;
        for (int i = 0; i < N; ++i) {
            CMV[0] += Particles[i]->vx;
            CMV[1] += Particles[i]->vy;
            CMV[2] += Particles[i]->vz;
        }
        CMV[0] /= N;
        CMV[1] /= N;
        CMV[2] /= N;
        return CMV;
    } else if (Dim == 2){
        double *CMV = new double[2];
        CMV[0] = 0.0;
        CMV[1] = 0.0;
        for (int i = 0; i < N; ++i) {
            CMV[0] += Particles[i]->vx;
            CMV[1] += Particles[i]->vy;
        }
        CMV[0] /= N;
        CMV[1] /= N;
        return CMV;
    } else {
        cout << "Error: Dimension is neither 2 nor 3!" << endl;
        return nullptr;
    }
}
//
vector<double>  Box::CalculateSSF(int axis) const {
    // Calculate Structure Factor along certain axis (x, y or z)
    int Nq = 100;
    vector<double> SSF (Nq, 0.0);
    double q, r;
    int counter = 0;
    switch (axis)
    {
        case 1:
            // x axis
            for (int i = 0; i < Nq; ++i)
            {
                q = 2.0 * M_PI * i / Lx;
                for (int j = 0; j < N; ++j)
                {
                    for (int k = 0; k < N; ++k)
                    {
                        r = Particles[j]->x - Particles[k]->x;
                        SSF[i] += cos(q * r);

                    }
                }
            }
            break;
        case 2:
            // y axis
            for (int i = 0; i < Nq; ++i)
            {
                q = 2.0 * M_PI * i / Ly;
                for (int j = 0; j < N; ++j)
                {
                    for (int k = 0; k < N; ++k)
                    {
                        r = Particles[j]->y - Particles[k]->y;
                        SSF[i] += cos(q * r);

                    }
                }
            }
            break;
        case 3:
            // z axis
            for (int i = 0; i < Nq; ++i)
            {
                q = 2.0 * M_PI * i / Lz;
                for (int j = 0; j < N; ++j)
                {
                    for (int k = 0; k < N; ++k)
                    {
                        r = Particles[j]->z - Particles[k]->z;
                        SSF[i] += cos(q * r);

                    }
                }
            }
            break;
    }
    return SSF;
}
//
double Box::CalculatePotentialEnergy(PotentialFunctionPosition PotentialFunction) {
    double PotentialEnergy = 0.0;
    for (int i = 0; i < N; ++i){
        PotentialEnergy += PotentialFunction(Particles[i]->x, Particles[i]->y, Particles[i]->z);
    }
    return PotentialEnergy;
}
double Box::CalculatePotentialEnergy(PotentialFunctionDistance PotentialFunction) {
    double PotentialEnergy = 0.0;
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            if (i != j){
                double dx = Particles[i]->x - Particles[j]->x;
                double dy = Particles[i]->y - Particles[j]->y;
                double dz = Particles[i]->z - Particles[j]->z;
                double dr = sqrt(dx * dx + dy * dy + dz * dz);
                PotentialEnergy += PotentialFunction(dr);
            }
        }
    }
    return PotentialEnergy/2;
}
//
double Box::CalculateMeanMinDistance() const {
    // calculate the average minimum distance between particles
    double mean_dist = 0, min_dist;
    int i,j;
    double r2;
    for(i=0; i<N; ++i)
    {
        min_dist = L + 10; // set to some high value
        for(j=0; j<N; j++)
        {
            if(j!=i)
            {
                double dx = Particles[i]->x - Particles[j]->x;
                double dy = Particles[i]->y - Particles[j]->y;
                double dz = Particles[i]->z - Particles[j]->z;
                // pbc
                if (dx > L/2.0) dx -= L;
                if (dx < -L/2.0) dx += L;
                if (dy > L/2.0) dy -= L;
                if (dy < -L/2.0) dy += L;
                if (dz > L/2.0) dz -= L;
                if (dz < -L/2.0) dz += L;
                r2 = dx*dx + dy*dy + dz*dz;
                // determine the next neighbour
                if (r2 < min_dist) min_dist = r2;
            }
        }
        min_dist = sqrt(min_dist);
        // avergage NN-distance over all particles
        mean_dist += min_dist;
    }
    return mean_dist/((double) N);
}
//
double Box::CalculatePackingFraction() const {
    double VParticle = 4. / 3. * M_PI * Particles[0]->myradius * Particles[0]->myradius * Particles[0]->myradius;
    return N * VParticle / V;
}
//
double Box::CalculateAbsolutAverageVelocity() const {
    double sumx = 0.0;
    double sumy = 0.0;
    for (int i = 0; i < N; ++i)
    {
        sumx += Particles[i]->vx;
        sumy += Particles[i]->vy;
    }
    return sqrt(sumx * sumx + sumy * sumy) / (N);
}

//----------------------------------------------------------------------
// Update Functions
void Box::UpdateKineticEnergy() {
    Ekin = CalculateKineticEnergy();
}
//
void Box::UpdatePotentialEnergy(PotentialFunctionPosition PotentialFunction) {
    Epot = CalculatePotentialEnergy(PotentialFunction);
}
//
void Box::UpdatePotentialEnergy(PotentialFunctionDistance PotentialFunction) {
    Epot = CalculatePotentialEnergy(PotentialFunction);
}
void Box::UpdateTotalEnergy() {
    Etot = Ekin + Epot;
}
//
void Box::UpdateForces(ForceFunctionDistance ForceFunction) {
    for (int i = 0; i < N; ++i) {
        Particles[i]->fx = 0.0;
        Particles[i]->fy = 0.0;
        Particles[i]->fz = 0.0;
        for (int j = 0; j < N; ++j) {
            if (i != j) {
                double dx = Particles[i]->x - Particles[j]->x;
                double dy = Particles[i]->y - Particles[j]->y;
                double dz = Particles[i]->z - Particles[j]->z;
                double dr = sqrt(dx * dx + dy * dy + dz * dz);
                double F = ForceFunction(dr);
                Particles[i]->fx += F * dx / dr;
                Particles[i]->fy += F * dy / dr;
                Particles[i]->fz += F * dz / dr;
            }
        }
    }
}
//
void Box::UpdateAbsoluteAverageVelocity() {
    AbsolutAverageVelocity = CalculateAbsolutAverageVelocity();
}
//
void Box::UpdateMeanMinDistance(){
    MeanMinDistance = CalculateMeanMinDistance();
}
//
void Box::UpdateOverlap(){
    // Check if there is overlap between particles using the Celllist
    shared_ptr<Particle> Current, Other;
    shared_ptr<Cell> OwnCell, NeighborCell;
    double dr;

    Overlap = false;
    for(int i = 0; i < N; ++i)
    {
        Current = Particles[i];
        OwnCell = Particles[i]->mycell;
        for(int j = 0; j < OwnCell->ParticlesInCell.size(); ++j)
        {
            Other = OwnCell->ParticlesInCell[j];
            if (Current == Other) continue; // skip self
            dr = Current->CalculateDistance(Other, Lx, Ly, Lz);
            if(dr < Current->myradius + Other->myradius){
                //cout << "Overlap between Particle " << Current->MyIndex << " and Particle " << Other->MyIndex << endl;
                Overlap = true;
                return;
            }
        }
        for(int k = 0; k < OwnCell->NeighbourCells.size(); ++k)
        {
            NeighborCell = OwnCell->NeighbourCells[k];
            for(int l = 0; l < NeighborCell->ParticlesInCell.size(); ++l)
            {
                Other = NeighborCell->ParticlesInCell[l];
                if (Current == Other) continue; // skip self
                dr = Current->CalculateDistance(Other, Lx, Ly, Lz);
                if(dr < Current->myradius + Other->myradius) {
                    //cout << "Overlap between Particle " << Current->MyIndex << " and Particle " << Other->MyIndex << endl;
                    Overlap = true;
                    return;
                }
            }
        }
    }
}


//----------------------------------------------------------------------
// Print Functions
void Box::PrintBoxInfo() const {
    cout << "Box Information:" << endl;
    cout << "N = " << N << endl;
    cout << "Lx = " << Lx << endl;
    cout << "Ly = " << Ly << endl;
    cout << "Lz = " << Lz << endl;
    cout << "V = " << V << endl;
    cout << "Rho = " << Rho << endl;
    cout << "Packing Fraction = " << Phi << endl;
    cout << "T = " << T << endl;
    cout << "Ekin = " << Ekin << endl;
    cout << "Epot = " << Epot << endl;
    cout << "Etot = " << Etot << endl;
    cout << "AbsolutAverageVelocity = " << AbsolutAverageVelocity << endl;
    cout << "Dim = " << Dim << endl;
    cout << "CMVx = " << CMVx << endl;
    cout << "CMVy = " << CMVy << endl;
    cout << "CMVz = " << CMVz << endl;
    cout << "CMx = " << CMx << endl;
    cout << "CMy = " << CMy << endl;
    cout << "CMz = " << CMz << endl;
    cout << endl;
}
//
void Box::PrintParticleInfo(int i) const {
    cout << "Particle " << i << ":" << endl;
    cout << "x = " << Particles[i]->x << endl;
    cout << "y = " << Particles[i]->y << endl;
    cout << "z = " << Particles[i]->z << endl;
    cout << "vx = " << Particles[i]->vx << endl;
    cout << "vy = " << Particles[i]->vy << endl;
    cout << "vz = " << Particles[i]->vz << endl;
    cout << "fx = " << Particles[i]->fx << endl;
    cout << "fy = " << Particles[i]->fy << endl;
    cout << "fz = " << Particles[i]->fz << endl;
    cout << "mycell = " << Particles[i]->mycell << endl;
    cout << endl;
}
//
void Box::PrintParticleInfo() const {
    for (int i = 0; i < N; ++i) {
        PrintParticleInfo(i);
    }
    cout << endl;

}
//
void Box::PrintCellsOfParticle(int i) const {
    cout << "Cell of Particle " << i << ":" << endl;
    cout << "mycell = " << Particles[i]->mycell << endl;
    cout << "with Index: " << Particles[i]->mycell->MyIndex << endl;
    cout << endl;
}
//
void Box::PrintCellOfParticle() const {
    for (int i = 0; i < N; ++i) {
        PrintCellsOfParticle(i);
    }
    cout << endl;
}

//----------------------------------------------------------------------
// Export Functions
void Box::ExportDisplayParticlePositions(string Filename, bool show) const {
    ofstream File;
    File.open(Filename);

    if (!File.is_open()) {
        std::cerr << "Error opening file: " << Filename << std::endl;
        exit(1);
        return;
    }

    switch (Dim) {
        case 3:
            for (int i = 0; i < N; ++i) {
                File << Particles[i]->x << " " << Particles[i]->y << " " << Particles[i]->z << endl;
            }
            File.close();

            if (show) {
                string command = "python code/visualize.py PlotParticlePositions 1";
                system(command.c_str());
            } else
            {
                string command = "python code/visualize.py PlotParticlePositions 0";
                system(command.c_str());
            }
            break;
        case 2:
            for (int i = 0; i < N; ++i) {
                File << Particles[i]->x << " " << Particles[i]->y << endl;
            }
            File.close();

            if (show) {
                string command = "python code/visualize.py PlotParticlePositions2D 1";
                system(command.c_str());
            } else
            {
                string command = "python code/visualize.py PlotParticlePositions2D 0";
                system(command.c_str());
            }
            break;
    }

}
//
void Box::ExportDisplayParticleVelocities(string Filename, bool show) const {
    ofstream File;
    File.open(Filename);

    for (int i = 0; i < N; ++i) {
        File << Particles[i]->vx << " " << Particles[i]->vy << " " << Particles[i]->vz
        << " " << Particles[i]->CalculateVelocityNorm() <<endl;
    }
    File.close();

    if (show) {
        string command = "python code/visualize.py PlotParticleVelocityDistribution 1 " + Filename;
        system(command.c_str());
    }else{
        string command = "python code/visualize.py PlotParticleVelocityDistribution 0 " + Filename;
        system(command.c_str());
    }
}
//
void Box::ExportCalculateRDF(string Filename, bool show) const {
    // compute the radial distribution function and
    // print it on a file with file name:  "gr_step"+mc_step+".dat"

    double Phi0 = CalculatePackingFraction();
    // constants for histogram
    int Nhisto = 200;
    double drhisto = L/((double) Nhisto)/2, dr;
    // declare and initialize histogram
    double histo[Nhisto];
    int i,j,drint;
    for(i=0; i<Nhisto; ++i)
    {
        histo[i] = 0.0;
    }
    // compute histogram
    shared_ptr<Particle> Current, Other;
    for(i=0; i<N; ++i)
    {
        Current = Particles[i];
        for(j=0; j<N; j++)
        {
            Other = Particles[j];
            if(Current != Other)
            {
                dr = Current->CalculateDistance(Other, Lx, Ly, Lz);
                drint = (int) ( dr/drhisto + 0.5);
                // bin the particle distance into a histogram
                if (drint < Nhisto) histo[drint]+=1.0;
            }
        }
    }
    // normalization
    double Volume;
    for(i=0; i<Nhisto; ++i)
    {
        Volume = 4.0/3.0*M_PI* ((i*drhisto + drhisto)*(i*drhisto + drhisto)*(i*drhisto + drhisto) - (i*drhisto)*(i*drhisto)*(i*drhisto) );
        histo[i] = histo[i]/((double) N  * Volume * 6.0/M_PI*Phi0);
    }

    // write histogram data
    ofstream File;
    File.open(Filename);
    if (!File.is_open()) {
        std::cerr << "Error opening file: " << Filename << std::endl;
        exit(1);
    }
    for(i=0; i<Nhisto; ++i)
    {
        File << i*drhisto << " " << histo[i] << endl;
    }
    File.close();

    if (show) {
        string command = "python code/visualize.py PlotRDF 1 " + Filename;
        system(command.c_str());
    }else{
        string command = "python code/visualize.py PlotRDF 0 " + Filename;
        system(command.c_str());
    }


}
//

//----------------------------------------------------------------------



















































