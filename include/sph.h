// SPH particle system class
#ifndef SPH_H
#define SPH_H

#include "particle.h"
#include "eos.h"
#include "kernels.h"

#include <vector>
#include <array>
#include <string>
#include <memory>

namespace SphSchemes
{
    enum class BoundaryConditions
    {
        destructive,
        periodic
    };

    struct SphParams
    {
        int dimensions;
        double fpSize;
        double bpSize;
        double restDensity;
        double dynamicViscosity;
        double gravity;
        double smoothingLength;
        SphSchemes::BoundaryConditions boundaryConditions;
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
            virtual ~SPH() {};

        protected:
            int dimensions;
            double restDensity;
            double pressureConstant;
            double dynamicViscosity;
            double extGravity;
            std::unique_ptr<SphEOS::EquationOfState> eos;
            std::unique_ptr<SphKernels::SphKernel> kernel;
            double fluidParticleSize;
            double boundaryParticleSize;
            double fluidParticleMass;
            double boundaryParticleMass;
            double smoothingLength;
            std::vector<int> coordOffsets;
            std::array<int, 3> gridDims;
            std::array<int, 3> dimFactors;
            int gridNbrhoodSize;
            std::array<double, 3> minPosition{
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max(),
                std::numeric_limits<double>::max()
                };
            std::array<double, 3> maxPosition{
                std::numeric_limits<double>::min(),
                std::numeric_limits<double>::min(),
                std::numeric_limits<double>::min()
                };
            std::array<double, 3> domainMin;
            std::array<double, 3> domainMax;
            std::size_t numBoundaryParticles;
            BoundaryConditions boundaryConditions;
            std::vector<std::vector<int>> grid;
            SPH(SphParams params, std::unique_ptr<SphKernels::SphKernel> kernel);

        private:
            void initCoordOffsetArray();
            void resizeGrid();
            void projectParticlesToGrid();
    };

    class WCSPH: public SPH
    {
        public:
            void printParameters() override;
            double getCFLTimestep(const double multiplier) override;
            WCSPH(SphParams params,
                  std::unique_ptr<SphKernels::SphKernel> kernel,
                  double maxH,
                  double rhoVar);

        private:
            double densityVariation;
            double soundSpeed;
    };

    class ThermoSPH: public SPH
    {
        public:
            void printParameters() override;
            ThermoSPH(SphParams params,
                      std::unique_ptr<SphKernels::SphKernel> kernel,
                      double gamma);

        private:
            double adiabaticIndex;
    };
}

#endif
