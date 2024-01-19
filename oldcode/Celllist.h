//
// Created by Yannick Hradetzky on 11.12.23.
//

#ifndef SIMULATION_CELLLIST_H
#define SIMULATION_CELLLIST_H




#include "vector"
#include "../code/Particle.h"
#include "../code/Cell.h"
#include "memory"
using namespace std;

class CellList {
public:
    //----------------------------------------------------------------------------------------------
    // Constructor
    CellList();
    CellList(double LCutOffx0, double LCutOffy0, double LCutOffz0, int DimCellList0,
             int Nmax0, double Lx0, double Ly0, double Lz0);
    CellList(double LCutOff, int DimCellList0, double L0);

    int Ncellsx, Ncellsy, Ncellsz;
    int DimCellList;
    double LCutOffx;
    double LCutOffy;
    double LCutOffz;
    double Lx, Ly, Lz;
//----------------------------------------------------------------------------------------------
    // Cell List Variables
    int Nneighbours;
    vector<shared_ptr<Particle>> CellListParticles;
    vector<shared_ptr<Cell>> Cells;

    //----------------------------------------------------------------------------------------------
    // Setter Functions
    void SetParticles(vector<shared_ptr<Particle>> Particles0);
    void SetCellsForParticles();
    void SetCellForParticle(int index);

    //----------------------------------------------------------------------------------------------
    // Initialization Functions
    void InitializeCells();

    //----------------------------------------------------------------------------------------------
    // Print Functions
    void PrintCellListInfo() const;
    void PrintCellListNeighbourInfo() const;



};


#endif //SIMULATION_CELLLIST_H
