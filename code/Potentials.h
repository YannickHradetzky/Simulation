//
// Created by Yannick Hradetzky on 28.12.23.
//

#ifndef SIMULATION_POTENTIALS_H
#define SIMULATION_POTENTIALS_H
// Here we will define functions that will be used to calculate the forces between particles

#include "cmath"
#include "memory"

double LennardJonesPotential(double r){
    // Calculate the Lennard-Jones potential between two particles
    // r is the distance between the particles
    // epsilon is the depth of the potential well
    // sigma is the distance at which the potential is zero
    double epsilon = 1.0;
    double sigma = 1.0;
    double r6 = pow(sigma / r, 6);
    double r12 = pow(r6, 2);
    return 4.0 * epsilon * (r12 - r6) + epsilon;
}
double LennardJonesForce(double r){
    // Calculate the Lennard-Jones force between two particles
    // r is the distance between the particles
    // epsilon is the depth of the potential well
    // sigma is the distance at which the potential is zero
    double epsilon = 1.0;
    double sigma = 1.0;
    double r2 = r * r;
    double r2inv = sigma / r2;
    double r6inv = r2inv * r2inv * r2inv;
    double r12inv = r6inv * r6inv;
    return 48.0 * epsilon * (r12inv - 0.5 * r6inv) * r2inv / sigma;
}

double GravityPotential(double r){
    // Calculate the gravitational potential between two particles
    // r is the distance between the particles
    // G is the gravitational constant
    // m is the mass of the particles
    double G = 1.0;
    double m = 1.0;
    return -G * m * m / r;
}
double GravityForce(double r){
// Calculate the gravitational force between two particles
    // r is the distance between the particles
    // G is the gravitational constant
    // m is the mass of the particles
    double G = 1.0;
    double m = 1.0;
    return G * m * m / (r * r);
}




#endif //SIMULATION_POTENTIALS_H
