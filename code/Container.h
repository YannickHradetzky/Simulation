//
// Created by Yannick Hradetzky on 19.12.23.
//

#ifndef SIMULATION_CONTAINER_H
#define SIMULATION_CONTAINER_H

#ifndef fstream
#include "fstream"
#endif

#ifndef iostream
#include "iostream"
#endif

#ifndef vector
#include "vector"
#endif

#ifndef Particle
#include "Particle.h"
#endif

#ifndef Cell
#include "Cell.h"
#endif


class Container {
public:
    //----------------------------------------------------------------------------------------------
    // Constructor and Destructor
    Container(){
        Lx = 0;
        Ly = 0;
        Lz = 0;
        L = 0;
        Dim = 0;
        Volume = 0;
        Surface = 0;
        Area = 0;
        Ncellsx = 0;
        Ncellsy = 0;
        Ncellsz = 0;
        Ncells = 0;
        //std::cout << "Default Container initialized" << std::endl;
    };
    Container(double lx0, double ly0, double lz0){
        Lx = lx0;
        Ly = ly0;
        Lz = lz0;
        if (Lx == Ly && Ly == Lz){
            L = Lx;
        }
        else{
            L = nan("Not Cubic");
        }
        Dim = 3;
        Volume = Lx*Ly*Lz;
        Surface = Lx*Ly*2 + Lx*Lz*2 + Ly*Lz*2;
        Area = nan("Not 2D");
        InitCells(1.0);
        std::cout << "Container initialized" << std::endl;
        std::cout << "  Lx = " << Lx << std::endl;
        std::cout << "  Ly = " << Ly << std::endl;
        std::cout << "  Lz = " << Lz << std::endl;
        std::cout << "  Ncellsx = " << Ncellsx << std::endl;
        std::cout << "  Ncellsy = " << Ncellsy << std::endl;
        std::cout << "  Ncellsz = " << Ncellsz << std::endl;

    };
    Container(double lx0, double ly0){
        Lx = lx0;
        Ly = ly0;
        Lz = 0;
        if (Lx == Ly){
            L = Lx;
        }
        else{
            L = nan("Not Quadratic");
        }
        Dim = 2;
        Volume = Lx*Ly*Lz;
        Surface = Lx*Ly;
        Area = Lx*Ly*2;
        InitCells(1.0);
        std::cout << "Container initialized" << std::endl;
        std::cout << "  Lx = " << Lx << std::endl;
        std::cout << "  Ly = " << Ly << std::endl;
        std::cout << "  Ncellsx = " << Ncellsx << std::endl;
        std::cout << "  Ncellsy = " << Ncellsy << std::endl;

    }
    ~Container(){};
    //----------------------------------------------------------------------------------------------
    // Container Variables
    double Lx, Ly, Lz, L;
    int Ncellsx, Ncellsy, Ncellsz, Ncells;
    double Volume, Surface, Area;
    int Dim; bool ParticlesOverlap;
    //----------------------------------------------------------------------------------------------
    // Container Functions
    void InitCells(double RInteraction)
    {
        Cells.clear();
        // Calculate Number of Cells in each Direction
        Ncellsx = floor(Lx / RInteraction);
        Ncellsy = floor(Ly / RInteraction);
        Ncellsz = floor(Lz / RInteraction);
        if (Ncellsy == Ncellsx && Ncellsx == Ncellsz) {
            Ncells = Ncellsy;
        } else {
            Ncells = nan("Not Cubic");
        }
        if(Ncellsx < 3 && Ncellsx > 0) {
            Ncellsx = 3;
            //std::cout << "  Ncellsx < 3, set to 3" << std::endl;
        }
        if(Ncellsy < 3 && Ncellsy > 0) {
            Ncellsy = 3;
            //std::cout << "  Ncellsy < 3, set to 3" << std::endl;
        }
        if(Ncellsz < 3 && Ncellsz > 0) {
            Ncellsz = 3;
            //std::cout << "  Ncellsz < 3, set to 3" << std::endl;
        }

        // Fill Cells with Pointers to empty Cells and set MyIndex as well as the Neighbour Cells (via Pointer)
        switch(Dim){
            case 3:
                // Add Cells to Cell List
                for(int i = 0; i < Ncellsx; ++i) {
                    for(int j = 0; j < Ncellsy; ++j) {
                        for(int k = 0; k < Ncellsz; ++k) {
                            Cell tempCell;
                            tempCell.MyIndex = i * Ncellsx * Ncellsy + j * Ncellsz + k;
                            Cells.push_back(std::make_shared<Cell>(tempCell));
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
                        Cells.push_back(std::make_shared<Cell>(tempCell));
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
        // Place Particles in Cells
        int xint, yint, zint;
        int CellIndex;
        switch (Dim) {
            case 3:
                for (const auto & CurrentParticle : Particles) {
                    xint = (int) (Ncellsx * (CurrentParticle->x / Lx + 0.5));
                    yint = (int) (Ncellsy * (CurrentParticle->y / Ly + 0.5));
                    zint = (int) (Ncellsz * (CurrentParticle->z / Lz + 0.5));
                    // Periodic Boundary Conditions
                    if (xint < 0) xint = 0;
                    if (yint < 0) yint = 0;
                    if (zint < 0) zint = 0;
                    if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                    if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                    if (zint > (Ncellsz - 1)) zint = Ncellsz - 1;
                    CellIndex = xint * Ncellsx * Ncellsy + yint * Ncellsz + zint;
                    Cells[CellIndex]->ParticlesInCell.push_back(CurrentParticle);
                    CurrentParticle->mycell = Cells[CellIndex];
                }
                break;
            case 2:
                for (const auto & CurrentParticle: Particles) {
                    xint = (int) (Ncellsx * (CurrentParticle->x / Lx + 0.5));
                    yint = (int) (Ncellsy * (CurrentParticle->y / Ly + 0.5));
                    // Periodic Boundary Conditions
                    if (xint < 0) xint = 0;
                    if (yint < 0) yint = 0;
                    if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                    if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                    CellIndex = xint * Ncellsy + yint;
                    Cells[CellIndex]->ParticlesInCell.push_back(CurrentParticle);
                    CurrentParticle->mycell = Cells[CellIndex];}
                break;
        }
    };
    void UpdateCell(std::shared_ptr<Particle> P) {
        int xint, yint, zint, CellIndex;
        std::shared_ptr<Cell> NewCell;
        switch (Dim)
        {
            case 3:
                xint = (int) (Ncellsx * (P->x / Lx + 0.5));
                yint = (int) (Ncellsy * (P->y / Ly + 0.5));
                zint = (int) (Ncellsz * (P->z / Lz + 0.5));
                if (xint < 0) xint = 0;
                if (yint < 0) yint = 0;
                if (zint < 0) zint = 0;
                if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                if (zint > (Ncellsz - 1)) zint = Ncellsz - 1;
                CellIndex = xint * Ncellsx * Ncellsy + yint * Ncellsz + zint;
                NewCell = Cells[CellIndex];
                if (NewCell != P->mycell) {
                    // find index of the Particle in Current Cell
                    int IndexInPCell = 0;
                    for (int i = 0; i < P->mycell->ParticlesInCell.size(); ++i) {
                        if (P->mycell->ParticlesInCell[i]->myindex == P->myindex) {
                            IndexInPCell = i;
                            break;
                        }
                    }
                    // Delete from Current Cell
                    P->mycell->ParticlesInCell.erase(P->mycell->ParticlesInCell.begin() + IndexInPCell);
                    // Add to new Cell
                    NewCell->ParticlesInCell.push_back(P);
                    // Update mycell
                    P->mycell = NewCell;
                }
                break;
            case 2:
                xint = (int) (Ncellsx * (P->x / Lx + 0.5));
                yint = (int) (Ncellsy * (P->y / Ly + 0.5));
                if (xint < 0) xint = 0;
                if (yint < 0) yint = 0;
                if (xint > (Ncellsx - 1)) xint = Ncellsx - 1;
                if (yint > (Ncellsy - 1)) yint = Ncellsx - 1;
                CellIndex = xint * Ncellsy + yint;
                NewCell = Cells[CellIndex];
                if (NewCell != P->mycell) {
                    // find index of the Particle in Current Cell
                    int IndexInPCell = 0;
                    for (int i = 0; i < P->mycell->ParticlesInCell.size(); ++i) {
                        if (P->mycell->ParticlesInCell[i]->myindex == P->myindex) {
                            IndexInPCell = i;
                            break;
                        }
                    }
                    // Delete from Current Cell
                    P->mycell->ParticlesInCell.erase(P->mycell->ParticlesInCell.begin() + IndexInPCell);
                    // Add to new Cell
                    NewCell->ParticlesInCell.push_back(P);
                    // Update mycell
                    P->mycell = NewCell;
                }
                break;
        }

    }
    void UpdateOverlap(){
        double dr;
        for(const auto& Current : Particles){
            for(const auto &Other : Current->mycell->ParticlesInCell){
                if(Other==Current) continue;
                dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < (Current->myradius + Other->myradius)){
                    ParticlesOverlap = true;
                    return;
                }else ParticlesOverlap = false;
            }
            for(const auto &NCell : Current->mycell->NeighbourCells){
                for(const auto &Other : NCell->ParticlesInCell){
                    if(Other==Current) continue;
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    if (dr < (Current->myradius + Other->myradius)){
                        ParticlesOverlap = true;
                        return;
                    }else ParticlesOverlap = false;
                }
            }
        }
    }
    void UpdateOverlap(const std::shared_ptr<Particle> P){
        // Check if Particle P overlaps with any other Particle at its current position
        // --> only check Particles in the same Cell and in the Neighbour Cells
        double dr;
        for(const auto &Other : P->mycell->ParticlesInCell){
            if(Other==P) continue;
            dr = P->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
            if (dr < (P->myradius + Other->myradius)){
                ParticlesOverlap = true;
                return;
            }else ParticlesOverlap = false;
        }
        for(const auto &NCell : P->mycell->NeighbourCells){
            for(const auto &Other : NCell->ParticlesInCell){
                if(Other==P) continue;
                dr = P->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < (P->myradius + Other->myradius)){
                    ParticlesOverlap = true;
                    return;
                }else ParticlesOverlap = false;
            }
        }

    }
    void UpdateOverlapLinear(){
        double dr;
        ParticlesOverlap = false;
        for(const auto& Current : Particles){
            for(const auto &Other : Current->mycell->ParticlesInCell){
                if(Other->myindex==Current->myindex) continue;
                dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < (Current->myradius + Other->myradius)){
                    ParticlesOverlap = true;
                    return;
                } else ParticlesOverlap = false;
            }
        }
    }
    bool CalculateOverlap(){
        double dr;
        for(const auto& Current : Particles){
            for(const auto &Other : Current->mycell->ParticlesInCell){
                if(Other==Current) continue;
                dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < (Current->myradius + Other->myradius)){
                    return true;
                }
            }
            for(const auto &NCell : Current->mycell->NeighbourCells){
                for(const auto &Other : NCell->ParticlesInCell){
                    if(Other==Current) continue;
                    dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                    if (dr < (Current->myradius + Other->myradius)){
                        return true;
                    }
                }
            }
        }
        return false;
    }
    bool CalculateOverlapLinear(){
        double dr;
        for(const auto& Current : Particles){
            for(const auto &Other : Current->mycell->ParticlesInCell){
                if(Other==Current) continue;
                dr = Current->CalculateDistanceToParticle(Other, Lx, Ly, Lz);
                if (dr < (Current->myradius + Other->myradius)){
                    return true;
                }
            }
        }
        return false;
    }
    void ScaleSystem(double Lnew, double Lold){
        double ScaleFactor = Lnew / Lold;
        for(const auto &Current : Particles){
            Current->x *= ScaleFactor;
            Current->y *= ScaleFactor;
            Current->z *= ScaleFactor;
        }
        Lx *= ScaleFactor;
        Ly *= ScaleFactor;
        Lz *= ScaleFactor;
        L = Lnew;
    }
    void SetParticleMass(double Mass){
            for(const auto &Current : Particles){
                Current->mymass = Mass;
            }
        }
    //----------------------------------------------------------------------------------------------
    std::vector<std::shared_ptr<Particle>> Particles;
    std::vector<std::shared_ptr<Cell>> Cells;
    //----------------------------------------------------------------------------------------------
    // Print Functions
    void PrintParticleInfo(int i) const {
        using namespace std;
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
        cout << "myindex = " << Particles[i]->myindex << endl;
        cout << "myradius = " << Particles[i]->myradius << endl;
        cout << endl;
    }
    void PrintParticleInfo() const {
        for(int i = 0; i < Particles.size(); ++i){
            PrintParticleInfo(i);
        }
    }
    void ExportVelocityDistribution(std::string Filename, bool show, bool makeplot) const {
        std::ofstream OutputFile;
        OutputFile.open(Filename);
        for(const auto &Current : Particles){
            OutputFile << Current->CalculateVelocityNorm() << std::endl;
        }
        OutputFile.close();

        if(makeplot){
            if(show){
                std::string command = "python code/visualize.py PlotParticleVelocityDistribution 1 " + Filename;
                system(command.c_str());
            }
            else{
                std::string command = "python code/visualize.py PlotParticleVelocityDistribution 0 " + Filename;
                system(command.c_str());
            }
        }
    }

    void ExportPositions(std::string Filename) const {
        std::ofstream OutputFile;
        OutputFile.open(Filename);
        for(const auto &Current : Particles){
            OutputFile << Current->x << " " << Current->y << " " << Current->z << std::endl;
        }
        OutputFile.close();
    }
    void ExportVelocities(std::string Filename) const {
        std::ofstream OutputFile;
        OutputFile.open(Filename);
        for(const auto &Current : Particles){
            OutputFile << Current->vx << " " << Current->vy << " " << Current->vz << std::endl;
        }
        OutputFile.close();
    }
    //----------------------------------------------------------------------------------------------
    //----------------------------------------------------------------------------------------------
};




#endif //SIMULATION_CONTAINER_H
