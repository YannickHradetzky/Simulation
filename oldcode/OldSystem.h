//
// Created by Yannick Hradetzky on 16.12.23.
//

#ifndef SIMULATION_SYSTEM_H
#define SIMULATION_SYSTEM_H

#include "Box.h"
#include "Celllist.h"
#include "../code/Particle.h"

#include <cstdlib>
#include "iostream"

class System: public Box, public CellList {
public:
    //----------------------------------------------------------------------------------------------
    // Constructor and Destructor
    System(){
         Container = new class Box;
        CellListSystem = new class CellList;
    }
    ~System(){
        delete Container;
        delete CellListSystem;
    }
    Box *Container;
    CellList *CellListSystem;
    // Make Pointers to the Particles and Cells
    vector<shared_ptr<Particle>> SystemParticles;
    vector<shared_ptr<Cell>> SystemCells;

    //----------------------------------------------------------------------------------------------
    // System Functions
    void InitLattice(int N, double Phi, double Radius, int NCells);
    void InitRandom2D(int N, double Rho, double Radius, int NCells);
    void InitRandom3D(int N, double Rho, double Radius, int NCells);
    void InitVicsek(int N, double Radius, double Velocity, double Rho, int NCells);
    void InitHardSpheres(int N, double Radius, double Phi);

};


#endif //SIMULATION_SYSTEM_H
