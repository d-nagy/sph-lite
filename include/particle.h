// SPH Particle class
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>

class Particle
{
    public:
        bool isBoundary = false;
        bool isActive = true;
        double x[3], v[3], a[3];
        double pressure, density, energyPerMass, energyPerMassDerivative;
        std::vector<int> neighbours;
        std::vector<double> neighbourDist;
        unsigned int gridCellNo;
};

#endif
