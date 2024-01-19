//
// Created by Yannick Hradetzky on 11.12.23.
//

#ifndef SIMULATION_PARTICLE_H
#define SIMULATION_PARTICLE_H


#include "memory"
#include "cmath"

class Cell;

class Particle {
    public:
        Particle(){
            x = 0.0;
            y = 0.0;
            z = 0.0;
            vx = 0.0;
            vy = 0.0;
            vz = 0.0;
            fx = 0.0;
            fy = 0.0;
            fz = 0.0;
            myindex = 0;
            myradius = 0.0;
            mymass = 0.0;
            mycell = nullptr;
            velocity = 0.0;
            theta = 0.0;
        }
        double x, y, z;
        double vx, vy, vz;
        double velocity;
        double theta;
        double fx, fy, fz;
        double mymass;
        double myradius;
        int myindex;
        // Create shared pointer which will point to Cell the particle is in
        std::shared_ptr<Cell> mycell;

    void ApplyPeriodicBoundaryConditions(double Lx, double Ly, double Lz) {
        if (x < 0.0) {
            x += Lx;
        }
        if (x > Lx) {
            x -= Lx;
        }
        if (y < 0.0) {
            y += Ly;
        }
        if (y > Ly) {
            y -= Ly;
        }
        if (z < 0.0) {
            z += Lz;
        }
        if (z > Lz) {
            z -= Lz;
        }
    }
    void ApplyCenteredPeriodicBoundaryConditions(double Lx, double Ly, double Lz) {
        if (x < -Lx/2.0) {
            x += Lx;
        }
        if (x > Lx/2.0) {
            x -= Lx;
        }
        if (y < -Ly/2.0) {
            y += Ly;
        }
        if (y > Ly/2.0) {
            y -= Ly;
        }
        if (z < -Lz/2.0 && Lz != 0.0) {
            z += Lz;
        }
        if (z > Lz/2.0 && Lz != 0.0) {
            z -= Lz;
        }
    }
    double CalculateVelocityNorm() const {
        return sqrt(vx*vx + vy*vy + vz*vz);
    }
    double CalculateDistanceToParticle(const std::shared_ptr<Particle>& Particle2, double Lx, double Ly, double Lz) const {
            double dx = x - Particle2->x;
            double dy = y - Particle2->y;
            double dz = z - Particle2->z;
            if (dx > Lx/2.0) {
                dx -= Lx;
            }
            if (dx < -Lx/2.0) {
                dx += Lx;
            }
            if (dy > Ly/2.0) {
                dy -= Ly;
            }
            if (dy < -Ly/2.0) {
                dy += Ly;
            }
            if (dz > Lz/2.0 && Lz != 0.0) {
                dz -= Lz;
            }
            if (dz < -Lz/2.0 && Lz != 0.0) {
                dz += Lz;
            }
            return sqrt(dx*dx + dy*dy + dz*dz);
        }
};



#endif //SIMULATION_PARTICLE_H
