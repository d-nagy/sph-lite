#include "sph.h"
#include "particle.h"

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <memory>
#include <omp.h>

#define coordOffsetsIndex(row, col, dim) ((dim*row) + col)

using SphSchemes::BoundaryConditions;
using SphSchemes::Particle;

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void SphSchemes::SPH::calcParticleDensities()
{
    std::array<int, 3> origin;
    std::array<double, 3> nbrPositionOffset;
    const double nbrRadius = kernel->getSupportRadius(smoothingLength);

    const auto nMax = particles.size();

    // #pragma omp parallel for
    for (std::size_t n=numBoundaryParticles; n<nMax; ++n)
    {
        Particle &op = particles[n];

        // Clear neighbour lists
        op.neighbours.clear();
        op.neighbourDist.clear();

        // Calculate this particle's grid coordinates
        int cellNo = op.gridCellNo;
        for (int d=dimensions-1; d>=0; d--)
        {
            int df = dimFactors[d];
            int c = cellNo / df;
            origin[d] = c;
            cellNo -= (c * df);
        }

        // Loop through all the neighbourhood grid cells
        double densitySum = 0.0;
        for (int i=0; i<gridNbrhoodSize; ++i)
        {
            // Calculate the cell index in the semi-flattened grid vector.
            cellNo = 0;
            bool validCoords = true;
            for (int d=0; d<dimensions; ++d)
            {
                int c = origin[d] + coordOffsets[coordOffsetsIndex(i, d, dimensions)];
                nbrPositionOffset[d] = 0;

                // Deal with boundary conditions
                switch (boundaryConditions)
                {
                    case BoundaryConditions::periodic:
                        if (c >= gridDims[d])
                        {
                            c %= gridDims[d];
                            nbrPositionOffset[d] += (domainMax[d] - domainMin[d]);
                        }
                        else if (c < 0)
                        {
                            c = gridDims[d] + c;
                            nbrPositionOffset[d] -= (domainMax[d] - domainMin[d]);
                        }
                        break;

                    default:
                        if (c < 0 || c >= gridDims[d])
                        {
                            validCoords = false;
                            break;
                        }
                }

                cellNo += (dimFactors[d] * c);
            }

            if (validCoords)
            {
                // Loop through all particles in the grid cell
                std::vector<int> &cell = grid[cellNo];
                for (const auto& pNbr : cell)
                {
                    Particle &nbr = particles[pNbr];
                    const double mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;

                    // Calculate distance to sample/origin particle
                    double dist = 0.0;
                    for (int d=0; d<dimensions; ++d)
                    {
                        const double dx = op.x[d] - (nbr.x[d] + nbrPositionOffset[d]);
                        dist += dx * dx;
                    }
                    dist = sqrt(dist);

                    // Check if neighbour
                    if (dist <= nbrRadius)
                    {
                        op.neighbours.push_back(pNbr);
                        op.neighbourDist.push_back(dist);
                        densitySum += (mass * kernel->W(dist, smoothingLength));
                    }
                }
            }
        }

        op.density = densitySum;
        const double eosPressure = eos->getPressure(op);
        op.pressure = (eosPressure < 0) ? 0 : eosPressure;
    }
}

// Calculate forces on particles from viscosity and pressure.
//
// Viscosity term requiring the Laplacian of the velocity
// is computed using Brookshaw's discrete Laplacian for SPH.
//
// Pressure term requiring the gradient of the pressure
// is computed using the symmetric SPH formula.
void SphSchemes::SPH::calcParticleForces()
{
    const auto nMax = particles.size();

    #pragma omp parallel for
    for (std::size_t n=numBoundaryParticles; n<nMax; ++n)
    {
        Particle &op = particles[n];

        // Reset particle accelerations
        op.a.fill(0.0);

        op.energyPerMassDerivative = 0.0;
        const double opPOverRhoSquared = op.pressure / (op.density * op.density);

        int i = -1;
        for (const auto& pNbr : op.neighbours)
        {
            Particle &nbr = particles[pNbr];
            const double dist = op.neighbourDist[++i];

            double gradWNorm = 0.0;
            const double mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;
            const double nbrPressure = (nbr.isBoundary) ? op.pressure : nbr.pressure;
            const double nbrPOverRhoSquared = nbrPressure / (nbr.density * nbr.density);

            for (int d=0; d<dimensions; ++d)
            {
                const double gradWComponent = kernel->gradWComponent(op.x[d], nbr.x[d], dist, smoothingLength);
                gradWNorm += gradWComponent * gradWComponent;

                // Add acceleration from pressure
                op.a[d] -= mass * (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;

                // Update energy per mass derivative
                op.energyPerMassDerivative += mass * (opPOverRhoSquared + nbrPOverRhoSquared) * (op.v[d] - nbr.v[d]) * gradWComponent / 2;
            }

            gradWNorm = sqrt(gradWNorm);

            if (!nbr.isBoundary && dist > 1e-12)
            {
                for (int d=0; d<dimensions; ++d)
                {
                    op.a[d] += (mass * 2 * gradWNorm * dynamicViscosity) * (nbr.v[d] - op.v[d]) / (op.density * nbr.density * dist);
                }
            }
        }
    }
}

// Move particles
// Semi-implicit Euler scheme
void SphSchemes::SPH::stepParticles(const double dt)
{
    const auto nMax = particles.size();

    #pragma omp parallel for
    for (std::size_t n=numBoundaryParticles; n<nMax; ++n)
    {
        Particle &op = particles[n];

        double speedSquared = 0.0;
        for (int d=0; d<dimensions; ++d)
        {
            // update particle accelerations due to external force (e.g. gravity)
            if (d == dimensions - 1)
                op.a[d] += extGravity;

            // Update velocities and positions
            // Semi-implicit Euler scheme
            op.v[d] += dt * op.a[d];
            op.x[d] += dt * op.v[d];

            // Update energypermass
            op.energyPerMass += dt * op.energyPerMassDerivative;
            if (op.energyPerMass < 0)
            {
                op.energyPerMass = 0;
            }

            // Handle boundary conditions
            switch (boundaryConditions)
            {
                case BoundaryConditions::periodic:
                    if (op.x[d] < domainMin[d])
                    {
                        op.x[d] += (domainMax[d] - domainMin[d]);
                    }
                    else if (op.x[d] > domainMax[d])
                    {
                        op.x[d] -= (domainMax[d] - domainMin[d]);
                    }
                    break;

                case BoundaryConditions::destructive:
                    if (op.x[d] < domainMin[d] || op.x[d] > domainMax[d])
                        op.isActive = false;
                    break;

                default:
                    break;
            }

            if (op.isActive)
            {
                if (op.x[d] < minPosition[d])
                {
                    minPosition[d] = op.x[d];
                }

                if (op.x[d] > maxPosition[d])
                {
                    maxPosition[d] = op.x[d];
                }
            }

            speedSquared += op.v[d]*op.v[d];
        }

        if (op.isActive && speedSquared > maxSpeedSquared)
        {
            #pragma omp atomic write
            maxSpeedSquared = speedSquared;
        }
    }

    if (boundaryConditions == BoundaryConditions::destructive)
    {
        particles.erase(std::remove_if(particles.begin(),
                                       particles.end(),
                                       [](const auto& p){ return !p.isActive; }),
                    particles.end());
    }
}

double SphSchemes::SPH::getCFLTimestep(const double multiplier)
{
    return multiplier * fluidParticleSize / sqrt(maxSpeedSquared);
}

double SphSchemes::WCSPH::getCFLTimestep(const double multiplier)
{
    return multiplier * fluidParticleSize / soundSpeed;
}

// Resize grid and project particles to it
void SphSchemes::SPH::setupParticleGrid()
{
    resizeGrid();
    projectParticlesToGrid();
}

// Initialise all particle positions and densities.
int SphSchemes::SPH::initialiseParticles(const std::string& casefileName)
{
    std::size_t numFluidParticles, numParticles = 0;

    std::ifstream caseFile (casefileName);
    std::string line;
    std::stringstream linestream;

    fluidParticleMass = restDensity * pow(fluidParticleSize, dimensions);
    boundaryParticleMass = restDensity * pow(boundaryParticleSize, dimensions);

    if (caseFile.is_open())
    {
        std::getline(caseFile, line);
        linestream.str(line);
        if (linestream >> numBoundaryParticles >> numFluidParticles)
            numParticles = numBoundaryParticles + numFluidParticles;

        if (!numParticles)
        {
            std::cout << "ERROR: No particles defined in case file." << '\n';
            return 0;
        }

        std::getline(caseFile, line);
        linestream.clear();
        linestream.str(line);
        linestream >> domainMin[0] >> domainMin[1] >> domainMin[2];
        linestream >> domainMax[0] >> domainMax[1] >> domainMax[2];

        // Initialise particles and set positions
        for (std::size_t p=0; p<numParticles; ++p)
        {
            Particle op = Particle();
            op.isBoundary = (p < numBoundaryParticles);
            op.density = restDensity;

            std::getline(caseFile, line);
            linestream.clear();
            linestream.str(line);
            for (int d=0; d<dimensions; ++d)
            {
                linestream >> op.x[d];

                if (op.x[d] < minPosition[d])
                    minPosition[d] = op.x[d];

                if (op.x[d] > maxPosition[d])
                    maxPosition[d] = op.x[d];
            }

            linestream >> op.isBoundary;
            if (!op.isBoundary)
            {
                linestream >> op.pressure >> op.energyPerMass;
            }

            particles.push_back(op);
        }

        caseFile.close();
    }
    else
    {
        std::cout << "ERROR: Unable to open case file." << '\n';
        return 0;
    }

    std::cout << "\tFluid particles:        " << numFluidParticles << '\n'
              << "\tBoundary particles:     " << numBoundaryParticles << '\n'
              << "\tTotal particles:        " << numParticles << '\n';

    return 1;
}

// Initialise array to hold relative coordinates of a grid cell neighbourhood
// to the central cell
void SphSchemes::SPH::initCoordOffsetArray()
{
    coordOffsets.resize(gridNbrhoodSize*dimensions);
    for (int i=0; i<gridNbrhoodSize; ++i)
    {
        for (int d=0; d<dimensions; ++d)
        {
            const int power = static_cast<int>(pow(3, dimensions-d-1));
            coordOffsets.at(coordOffsetsIndex(i, d, dimensions)) = ((i / power) % 3) - 1;
        }
    }
}

// Create grid for all particles.
void SphSchemes::SPH::resizeGrid()
{
    int gridSize = 1, dimFactor = 1;
    const double cellSize = kernel->getSupportRadius(smoothingLength);
    bool needToResize = false;

    for (int d=0; d<dimensions; ++d)
    {
        const int dim = static_cast<int>(ceil((maxPosition[d] - minPosition[d] + fluidParticleSize) / cellSize));
        if (gridDims[d] != dim) { needToResize = true; }
        gridSize *= dim;
        gridDims[d] = dim;
        dimFactors[d] = dimFactor;
        dimFactor *= dim;
    }

    if (needToResize)
    {
        grid.clear();
        grid.resize(gridSize);
    }
}

// Assign each particle to a grid cell.
void SphSchemes::SPH::projectParticlesToGrid()
{
    int p = -1;
    const double cellSize = kernel->getSupportRadius(smoothingLength);
    const double halfParticleSize = fluidParticleSize / 2;

    for (auto& cell : grid)
    {
        cell.clear();
    }

    for (auto& op : particles)
    {
        ++p;
        int cellNo = 0;
        for (int d=0; d<dimensions; ++d)
        {
            cellNo += dimFactors[d] * static_cast<int>((op.x[d] - minPosition[d] + halfParticleSize) / cellSize);
        }

        grid.at(cellNo).push_back(p);
        op.gridCellNo = cellNo;
    }
}

void SphSchemes::SPH::printParameters()
{
    std::string bcString;
    switch (boundaryConditions)
    {
        case BoundaryConditions::periodic:
            bcString = "periodic";
            break;
        case BoundaryConditions::destructive:
            bcString = "destructive";
            break;
    }

    std::cout << "STANDARD SPH PARAMETERS:" << '\n'
              << "\tDimensions:             " << dimensions << '\n'
              << "\tRest density:           " << restDensity << '\n'
              << "\tPressure constant:      " << pressureConstant << '\n'
              << "\tDynamic viscosity:      " << dynamicViscosity << '\n'
              << "\tExternal gravity:       " << extGravity << '\n'
              << "\tEquation of state:      " << *eos << '\n'
              << "\tFluid particle size:    " << fluidParticleSize << '\n'
              << "\tBoundary particle size: " << boundaryParticleSize << '\n'
              << "\tSmoothing length:       " << smoothingLength << '\n'
              << "\tSPH kernel:             " << *kernel << '\n'
              << "\tFluid particle mass:    " << fluidParticleMass << '\n'
              << "\tBoundary particle mass: " << boundaryParticleMass << '\n'
              << "\tBoundary conditions:    " << bcString << '\n';
}

void SphSchemes::WCSPH::printParameters()
{
    SPH::printParameters();

    std::cout << '\n' << "WCSPH PARAMETERS:" << '\n'
              << "\tDensity variation:      " << densityVariation << '\n'
              << "\tSound speed:            " << soundSpeed << '\n';
}

void SphSchemes::ThermoSPH::printParameters()
{
    SPH::printParameters();

    std::cout << '\n' << "ThermoSPH PARAMETERS:" << '\n'
              << "\tAdiabatic index:        " << adiabaticIndex << '\n';
}

SphSchemes::SPH::SPH(SphSchemes::SphParams params, std::unique_ptr<SphKernels::SphKernel> kernel)
{
    dimensions = params.dimensions;
    restDensity = params.restDensity;
    dynamicViscosity = params.dynamicViscosity;
    extGravity = params.gravity;
    this->kernel = std::move(kernel);
    fluidParticleSize = params.fpSize;
    boundaryParticleSize = params.bpSize;
    smoothingLength = params.smoothingLength;
    boundaryConditions = params.boundaryConditions;

    gridNbrhoodSize = static_cast<int>(pow(3, dimensions));
    fluidParticleMass = restDensity * pow(fluidParticleSize, dimensions);
    boundaryParticleMass = restDensity * pow(boundaryParticleSize, dimensions);

    initCoordOffsetArray();
}

SphSchemes::WCSPH::WCSPH(SphSchemes::SphParams params,
                         std::unique_ptr<SphKernels::SphKernel> kernel,
                         double maxH,
                         double rhoVar
                         ) : SPH(params, std::move(kernel))
{
    densityVariation = rhoVar;

    const double soundSpeedSquared = 2 * abs(extGravity) * maxH / densityVariation;
    pressureConstant = restDensity * soundSpeedSquared / 7;
    soundSpeed = sqrt(soundSpeedSquared);

    eos = std::make_unique<SphEOS::WeaklyCompressibleEOS>(pressureConstant, restDensity);
}

SphSchemes::ThermoSPH::ThermoSPH(SphSchemes::SphParams params,
                                 std::unique_ptr<SphKernels::SphKernel> kernel,
                                 double gamma
                                 ) : SPH(params, std::move(kernel))
{
    adiabaticIndex = gamma;
    eos = std::make_unique<SphEOS::IdealGasEOS>(adiabaticIndex);
}
