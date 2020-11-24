// SPH particle system class
#ifndef SPH_H
#define SPH_H

#include "particle.h"
#include "eos.h"
#include "kernels.h"

#include <vector>
#include <string>
#include <memory>

namespace SphSchemes
{
    enum class BoundaryConditions
    {
        destructive,
        periodic
    };

    class SPH
    {
        public:
            double maxSpeedSquared = 0.0;
            std::vector<Particle> particles;
            int initialiseParticles(const std::string& casefileName);
            void setupParticleGrid();
            virtual double getCFLTimestep(const double multiplier);
            virtual void calcParticleDensities();
            virtual void calcParticleForces();
            virtual void stepParticles(const double dt);
            virtual void printParameters();
            SPH(int d,
                double rho0,
                double eta,
                double g,
                SphKernels::SphKernel *kernel,
                double fluidSpacing,
                double boundarySpacing,
                double h,
                BoundaryConditions bc);
            virtual ~SPH();

        protected:
            int dimensions;
            double restDensity;
            double pressureConstant;
            double dynamicViscosity;
            double extGravity;
            std::unique_ptr<SphEOS::EquationOfState> eos;
            SphKernels::SphKernel *kernel;
            double fluidParticleSize;
            double boundaryParticleSize;
            double smoothingLength;
            double fluidParticleMass;
            double boundaryParticleMass;
            std::vector<int> coordOffsets;
            int gridDims[3];
            int gridNbrhoodSize;
            int dimFactors[3];
            double minPosition[3];
            double maxPosition[3];
            double domainMin[3];
            double domainMax[3];
            int numBoundaryParticles;
            BoundaryConditions boundaryConditions;
            std::vector<std::vector<int>> grid;
            void initCoordOffsetArray();
            void resizeGrid();
            void projectParticlesToGrid();
    };

    class WCSPH: public SPH
    {
        public:
            void printParameters();
            double getCFLTimestep(const double multiplier);
            WCSPH(int d,
                double rho0,
                double eta,
                double g,
                SphKernels::SphKernel *kernel,
                double fluidSpacing,
                double boundarySpacing,
                double h,
                BoundaryConditions bc,
                double maxH,
                double rhoVar);

        private:
            double densityVariation;
            double soundSpeed;
    };

    class ThermoSPH: public SPH
    {
        public:
            void printParameters();
            ThermoSPH(int d,
                    double rho0,
                    double eta,
                    double g,
                    SphKernels::SphKernel *kernel,
                    double fluidSpacing,
                    double boundarySpacing,
                    double h,
                    BoundaryConditions bc,
                    double gamma);

        private:
            double adiabaticIndex;
    };
}

#endif
