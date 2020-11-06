// SPH particle system class
#ifndef SPH_H
#define SPH_H

#include "particle.h"
#include "eos.h"
#include "kernels.h"

#include <vector>
#include <string>

class SPH
{
    public:
        double maxSpeedSquared = 0.0;
        std::vector<Particle> particles;
        int initialiseParticles(const std::string& casefileName);
        void setupParticleGrid();
        void calcParticleDensities();
        void calcParticleForces();
        void stepParticles(double dt);
        double getCFLTimestep(double multiplier);
        void printParameters();
        SPH(int d,
            double rd,
            double dv,
            double s,
            double ai,
            double g,
            EquationOfState *eos,
            SphKernel *kernel,
            double fps,
            double bps,
            double sl);
        ~SPH();

    private:
        int dimensions;
        double restDensity;
        double dynamicViscosity;
        double stiffness;
        double adiabaticIndex;
        double extGravity;
        EquationOfState *eos;
        SphKernel *kernel;
        double fluidParticleSize;
        double boundaryParticleSize;
        double smoothingLength;
        double fluidParticleMass;
        double boundaryParticleMass;
        int* coordOffsets;
        int gridDims[3];
        int dimFactors[3];
        double minPosition[3];
        double maxPosition[3];
        std::vector<std::vector<int>> grid;
        void initCoordOffsetArray();
        void resizeGrid();
        void projectParticlesToGrid();
        int getGridFlatIndex(int const coordArray[], int const offsetArray[]);
};

#endif
