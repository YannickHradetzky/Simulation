//
// Created by Yannick Hradetzky on 11.12.23.
//

#include <iostream>
#include <utility>
#include "Celllist.h"
#include "../code/Particle.h"

using namespace std;

//----------------------------------------------------------------------
// Constructors
CellList::CellList() {
    LCutOffx = 0.0;
    LCutOffy = 0.0;
    LCutOffz = 0.0;
    Lx = 0.0;
    Ly = 0.0;
    Lz = 0.0;
    Ncellsx = Ncellsy = Ncellsz = 0;
    DimCellList = 0;
    Nneighbours = 0;
    Ncellsx = Ncellsy = Ncellsz = 0;
}
//
CellList::CellList(double LCutOffx0, double LCutOffy0, double LCutOffz0,
                   int DimCellList0, int Nmax0, double Lx0, double Ly0, double Lz0) {
    LCutOffx = LCutOffx0;
    LCutOffy = LCutOffy0;
    LCutOffz = LCutOffz0;
    Lx = Lx0;
    Ly = Ly0;
    Lz = Lz0;
    DimCellList = DimCellList0;
    Ncellsx = (int) (Lx0 / LCutOffx);
    Ncellsy = (int) (Ly0 / LCutOffy);
    Ncellsz = (int) (Lz0 / LCutOffz);
    if (Ncellsx < 3 || Ncellsy < 3 || Ncellsz < 3) {
        cout << "Error: Ncellsx, Ncellsy and Ncellsz must be at least 3 (Constructor CellListSystem)" << endl;
        exit(1);
    }
    // Set number of neighbours depending on dimension of cell list
    if (DimCellList == 2){
        Nneighbours = 8;
    } else if (DimCellList == 3) {
        Nneighbours = 26;
    } else {
        cout << "Error: DimCellList must be 2 or 3 (Constructor CellListSystem)" << endl;
        exit(1);
    }
}
//
CellList::CellList(double LCutOff, int DimCellList0, double L0) {
    LCutOffx = LCutOff;
    LCutOffy = LCutOff;
    LCutOffz = LCutOff;
    Lx = Ly = Lz = L0;
    DimCellList = DimCellList0;
    Ncellsx = (int) (L0 / LCutOff);
    Ncellsy = (int) (L0 / LCutOff);
    Ncellsz = (int) (L0 / LCutOff);
    if (Ncellsx < 3 || Ncellsy < 3 || Ncellsz < 3) {
        cout << "Error: Ncellsx, Ncellsy and Ncellsz must be at least 3 (Constructor CellListSystem)" << endl;
        exit(1);
    }
    // Set number of neighbours depending on dimension of cell list
    if (DimCellList == 2){
        Nneighbours = 8;
    } else if (DimCellList == 3) {
        Nneighbours = 26;
    } else {
        cout << "Error: DimCellList must be 2 or 3 (Constructor CellListSystem)" << endl;
        exit(1);
    }
}
//


//----------------------------------------------------------------------
// Setter Functions
void CellList::SetParticles(vector<shared_ptr<Particle>> Particles0) {
    CellListParticles = Particles0;
}
//
void CellList::SetCellsForParticles() {
    int xint, yint, zint;
    int CellIndex;
    switch (DimCellList) {
        case 3:
            for (const auto & CellListParticle : CellListParticles) {
                xint = (int) (Ncellsx * (CellListParticle->x / Lx + 0.5));
                yint = (int) (Ncellsy * (CellListParticle->y / Ly + 0.5));
                zint = (int) (Ncellsz * (CellListParticle->z / Lz + 0.5));
                if (DimCellList == 2) zint = 0;
                // Periodic Boundary Conditions
                if (xint < 0) xint = 0;
                if (yint < 0) yint = 0;
                if (zint < 0) zint = 0;
                if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                if (zint > (Ncellsz - 1)) zint = Ncellsz - 1;
                CellIndex = xint * Ncellsx * Ncellsy + yint * Ncellsz + zint;
                Cells[CellIndex]->ParticlesInCell.push_back(CellListParticle);
                Cells[CellIndex]->NParticlesInCell = Cells[CellIndex]->ParticlesInCell.size(); // Update NParticlesInCell
                CellListParticle->mycell = Cells[CellIndex];
            }
            break;
        case 2:
            for (const auto & CellListParticle : CellListParticles) {
                xint = (int) (Ncellsx * (CellListParticle->x / Lx + 0.5));
                yint = (int) (Ncellsy * (CellListParticle->y / Ly + 0.5));
                // Periodic Boundary Conditions
                if (xint < 0) xint = 0;
                if (yint < 0) yint = 0;
                if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                CellIndex = xint * Ncellsy + yint;
                Cells[CellIndex]->ParticlesInCell.push_back(CellListParticle);
                Cells[CellIndex]->NParticlesInCell = Cells[CellIndex]->ParticlesInCell.size(); // Update NParticlesInCell
                CellListParticle->mycell = Cells[CellIndex];
            }
            break;
    }
}
//
void CellList::SetCellForParticle(int index) {
    int xint, yint, zint;
    int CellIndex, IndexInCurrentCell = -1;
    shared_ptr<Particle> P = CellListParticles[index];
    shared_ptr<Cell> CurrentCell = P->mycell;
    shared_ptr<Cell> NewCell = nullptr;

    switch (DimCellList) {
        case 3:
            xint = (int) (Ncellsx * (CellListParticles[index]->x / Lx + 0.5));
            yint = (int) (Ncellsy * (CellListParticles[index]->y / Ly + 0.5));
            zint = (int) (Ncellsz * (CellListParticles[index]->z / Lz + 0.5));
            if (xint < 0) xint = 0;
            if (yint < 0) yint = 0;
            if (zint < 0) zint = 0;
            if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
            if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
            if (zint > (Ncellsz - 1)) zint = Ncellsz - 1;
            CellIndex = xint * Ncellsx * Ncellsy + yint * Ncellsz + zint;
            NewCell = Cells[CellIndex];
            if (NewCell != CurrentCell) {
                // find index of the Particle in Current Cell
                for (int i = 0; i < CurrentCell->ParticlesInCell.size(); ++i) {
                    if (CurrentCell->ParticlesInCell[i]->myindex == P->myindex) {
                        IndexInCurrentCell = i;
                        break;
                    }
                }
                // Delete from Current Cell
                CurrentCell->ParticlesInCell.erase(CurrentCell->ParticlesInCell.begin() + IndexInCurrentCell);
                // Add to new Cell
                NewCell->ParticlesInCell.push_back(P);
                // Update mycell
                P->mycell = NewCell;
            }
            break;
        case 2:
            xint = (int) (Ncellsx * (CellListParticles[index]->x / Lx + 0.5));
            yint = (int) (Ncellsy * (CellListParticles[index]->y / Ly + 0.5));
            // Periodic Boundary Conditions
            if (xint < 0) xint = 0;
            if (yint < 0) yint = 0;
            if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
            if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
            CellIndex = xint * Ncellsy + yint;
            NewCell = Cells[CellIndex];
            if (NewCell != CurrentCell)
            {
                // find index of the Particle in Current Cell
                for (int i = 0; i < CurrentCell->ParticlesInCell.size(); ++i)
                {
                    if (CurrentCell->ParticlesInCell[i]->myindex == P->myindex)
                    {
                        IndexInCurrentCell = i;
                        break;
                    }
                }
                // Delete from Current Cell
                CurrentCell->ParticlesInCell.erase(CurrentCell->ParticlesInCell.begin() + IndexInCurrentCell);
                // Add to new Cell
                NewCell->ParticlesInCell.push_back(P);
                // Update mycell
                P->mycell = NewCell;
            }
            break;
    }
}

//----------------------------------------------------------------------
// Initialize Cell List
void CellList::InitializeCells() {
    Cells.clear();
    // Fill Cells with Pointers to empty Cells and set MyIndex as well as the Neighbour Cells (via Pointer)
    switch(DimCellList){
        case 3:
            // Add Cells to Cell List
            for(int i = 0; i < Ncellsx; ++i) {
                for(int j = 0; j < Ncellsy; ++j) {
                    for(int k = 0; k < Ncellsz; ++k) {
                        Cell tempCell;
                        tempCell.MyIndex = i * Ncellsx * Ncellsy + j * Ncellsz + k;
                        Cells.push_back(make_shared<Cell>(tempCell));
                    }
                }
            }
            // Set Neighbour Cells
            for (int i = 0; i < Ncellsx; ++i) {
                for (int j = 0; j < Ncellsy; ++j) {
                    for (int k = 0; k < Ncellsz; ++k) {
                        int CellIndex = i * Ncellsx * Ncellsy + j * Ncellsz + k;
                        Cells[CellIndex]->NeighbourCells.clear();
                        for(int ii = i-1; ii <= i+1; ++ii) {
                            for(int jj = j-1; jj <= j+1; ++jj) {
                                for(int kk = k-1; kk <= k+1; ++kk) {
                                    if (ii == i && jj == j && kk == k) continue;
                                    int a = ii;
                                    int b = jj;
                                    int c = kk;
                                    // Periodic Boundary Conditions
                                    if(a < 0) a += Ncellsx;
                                    if(a >= Ncellsx) a -= Ncellsx;
                                    if(b < 0) b += Ncellsy;
                                    if(b >= Ncellsy) b -= Ncellsy;
                                    if(c < 0) c += Ncellsz;
                                    if(c >= Ncellsz) c -= Ncellsz;
                                    int NeighbourCellIndex = a * Ncellsy * Ncellsx + b * Ncellsz + c;
                                    Cells[CellIndex]->NeighbourCells.push_back(Cells[NeighbourCellIndex]);
                                }
                            }
                        }
                    }
                }
            }
            break;
        case 2:
            Ncellsz = 0;
            // Add Cells to Cell List
            for(int i = 0; i < Ncellsx; ++i) {
                for(int j = 0; j < Ncellsy; ++j) {
                    Cell tempCell;
                    tempCell.MyIndex = i * Ncellsy + j;
                    Cells.push_back(make_shared<Cell>(tempCell));
                }
            }
            // Set Neighbour Cells
            for (int i = 0; i < Ncellsx; ++i) {
                for (int j = 0; j < Ncellsy; ++j) {
                    int CellIndex = i * Ncellsy + j;
                    Cells[CellIndex]->NeighbourCells.clear();
                    for(int ii = i-1; ii <= i+1; ++ii) {
                        for(int jj = j-1; jj <= j+1; ++jj) {
                            if(!(ii == i && jj == j)) {
                                int a = ii;
                                int b = jj;
                                // Periodic Boundary Conditions
                                if(a < 0) a += Ncellsx;
                                if(a >= Ncellsx) a -= Ncellsx;
                                if(b < 0) b += Ncellsy;
                                if(b >= Ncellsy) b -= Ncellsy;
                                int NeighbourCellIndex = a * Ncellsy + b;
                                Cells[CellIndex]->NeighbourCells.push_back(Cells[NeighbourCellIndex]);
                            }
                        }
                    }
                }
            }
            break;
    }
    SetCellsForParticles();
}
//

//----------------------------------------------------------------------
// Print Functions
void CellList::PrintCellListInfo() const {
    cout << "CellListSystem Information:" << endl;
    cout << "LCutOffx = " << LCutOffx << endl;
    cout << "LCutOffy = " << LCutOffy << endl;
    cout << "LCutOffz = " << LCutOffz << endl;
    cout << "Ncellsx = " << Ncellsx << endl;
    cout << "Ncellsy = " << Ncellsy << endl;
    cout << "Ncellsz = " << Ncellsz << endl;
    cout << "DimCellList = " << DimCellList << endl;
    cout << "N = " << CellListParticles.size() << endl;
    cout << "Nneighbours = " << Nneighbours << endl;
    cout << "CellListParticles.size() = " << CellListParticles.size() << endl;
    cout << "Cells.size() = " << Cells.size() << endl;
    cout << endl;
}
//
void CellList::PrintCellListNeighbourInfo() const {
    cout << "CellListSystem Neighbour Information:" << endl;
    cout << "Cells.size() = " << Cells.size() << endl;
    for(int i = 0; i < Cells.size(); ++i) {
        cout << "Cell " << i << " has " << Cells[i]->NeighbourCells.size() << " neighbours:" << endl;
        for(int j = 0; j < Cells[i]->NeighbourCells.size(); ++j) {
            cout << "Cell " << Cells[i]->NeighbourCells[j]->MyIndex << " ";
        }
        cout << endl;
        cout << "and contains " << Cells[i]->ParticlesInCell.size() << " particles with Indexes:." << endl;
        for(int j = 0; j < Cells[i]->ParticlesInCell.size(); ++j) {
            cout << Cells[i]->ParticlesInCell[j]->myindex << " ";
        }

        cout << endl;
    }
    cout << endl;
}

//----------------------------------------------------------------------