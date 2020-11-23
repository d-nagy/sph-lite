// SPH Particle class
#ifndef PARTICLE_H
#define PARTICLE_H

#include <vector>
#include <array>

class Particle
{
    public:
        bool isBoundary = false;
        bool isActive = true;
        std::array<double, 3> x, v, a;
        double pressure, density, energyPerMass, energyPerMassDerivative;
        std::vector<int> neighbours;
        std::vector<double> neighbourDist;
        unsigned int gridCellNo;
};

#endif
