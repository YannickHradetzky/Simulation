//
// Created by Yannick Hradetzky on 12.12.23.
//

#ifndef SIMULATION_CELL_H
#define SIMULATION_CELL_H

#ifndef vector
#include "vector"
#endif

#ifndef Particle
#include "Particle.h"
#endif



class Cell {
public:
    Cell(){
        MyIndex = -1;
    };
    std::vector<std::shared_ptr<Particle>> ParticlesInCell; // contains pointers to all particles in the cell
    std::vector<std::shared_ptr<Cell>> NeighbourCells; // contains pointers to all neighbour cells
    int MyIndex; // index of the cell
};


#endif //SIMULATION_CELL_H
