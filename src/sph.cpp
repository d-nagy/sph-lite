#include "sph.h"
#include "particle.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <omp.h>

#define coordOffsetsIndex(row, col, dim) ((dim*row) + col)
#define max(x, y) ((x > y) ? x : y)
#define min(x, y) ((x > y) ? y : x)

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void SPH::calcParticleDensities()
{
    int origin[dimensions];
    double densitySum, mass, dist, dx;
    int cellNo, c, df;
    double nbrRadius = kernel->getSupportRadius(smoothingLength);
    bool validCoords;
    double nbrPositionOffset[3];

    int nMax = particles.size();

    // #pragma omp parallel for private(densitySum, mass, dist, dx, cellNo, c, df, validCoords)
    for (int n=numBoundaryParticles; n<nMax; n++)
    {
        Particle &op = particles[n];

        // Clear neighbour lists
        op.neighbours.clear();
        op.neighbourDist.clear();

        // Calculate this particle's grid coordinates
        cellNo = op.gridCellNo;
        for (int d=dimensions-1; d>=0; d--)
        {
            df = dimFactors[d];
            c = cellNo / df;
            origin[d] = c;
            cellNo -= (c * df);
        }

        // Loop through all the neighbourhood grid cells
        densitySum = 0.0;
        for (int i=0; i<gridNbrhoodSize; i++)
        {
            // Calculate the cell index in the semi-flattened grid vector.
            cellNo = 0;
            validCoords = true;
            for (int d=0; d<dimensions; d++)
            {
                c = origin[d] + coordOffsets[coordOffsetsIndex(i, d, dimensions)];
                nbrPositionOffset[d] = 0;

                // Deal with boundary conditions
                switch (boundaryConditions)
                {
                    case periodic:
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
                    mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;

                    // Calculate distance to sample/origin particle
                    dist = 0;
                    for (int d=0; d<dimensions; d++)
                    {
                        dx = op.x[d] - (nbr.x[d] + nbrPositionOffset[d]);
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
        op.pressure = max(eos->getPressure(op), 0);
    }
}

// Calculate forces on particles from viscosity and pressure.
//
// Viscosity term requiring the Laplacian of the velocity
// is computed using Brookshaw's discrete Laplacian for SPH.
//
// Pressure term requiring the gradient of the pressure
// is computed using the symmetric SPH formula.
void SPH::calcParticleForces()
{
    double gradWNorm, gradWComponent, gradW[3];
    double mass, dist, opPOverRhoSquared, nbrPOverRhoSquared, nbrPressure;

    int nMax = particles.size();

    #pragma omp parallel for private(gradWNorm, gradWComponent, mass, dist, opPOverRhoSquared, nbrPOverRhoSquared, nbrPressure)
    for (int n=numBoundaryParticles; n<nMax; n++)
    {
        Particle &op = particles[n];

        // Reset particle accelerations
        for (int d=0; d<dimensions; d++)
        {
            op.a[d] = 0.0;
        }

        op.energyPerMassDerivative = 0;
        opPOverRhoSquared = op.pressure / (op.density * op.density);

        int i=0;
        for (const auto& pNbr : op.neighbours)
        {
            Particle &nbr = particles[pNbr];
            dist = op.neighbourDist[i++];

            gradWNorm = 0;
            mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;
            nbrPressure = (nbr.isBoundary) ? op.pressure : nbr.pressure;
            nbrPOverRhoSquared = nbrPressure / (nbr.density * nbr.density);

            for (int d=0; d<dimensions; d++)
            {
                gradWComponent = kernel->gradWComponent(op.x[d], nbr.x[d], dist, smoothingLength);
                gradWNorm += gradWComponent * gradWComponent;
                gradW[d] = gradWComponent;

                // Add acceleration from pressure
                op.a[d] -= mass * (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;

                // Update energy per mass derivative
                op.energyPerMassDerivative += mass * (opPOverRhoSquared + nbrPOverRhoSquared) * (op.v[d] - nbr.v[d]) * gradWComponent / 2;
            }

            gradWNorm = sqrt(gradWNorm);

            if (!nbr.isBoundary && dist > 1e-12)
            {
                for (int d=0; d<dimensions; d++)
                {
                    op.a[d] += (mass * 2 * gradWNorm * dynamicViscosity) * (nbr.v[d] - op.v[d]) / (op.density * nbr.density * dist);
                }
            }
        }
    }
}

// void SPH::addViscosityAcceleration(Particle& op, const Particle& nbr, const double (&gradW)[3], const double dist, const double mass)
// {
//     double gradWNorm = 0;
//     for (int d=0; d<dimensions; d++)
//     {
//         gradWNorm += gradW[d]*gradW[d];
//     }

//     gradWNorm = sqrt(gradWNorm);

//     if (dist > 1e-12)
//     {
//         for (int d=0; d<dimensions; d++)
//         {
//             op.a[d] += (mass * 2 * gradWNorm * dynamicViscosity) * (nbr.v[d] - op.v[d]) / (op.density * nbr.density * dist);
//         }
//     }
// }

// void WCSPH::addViscosityAcceleration(Particle& op, const Particle& nbr, const double (&gradW)[3], const double dist, const double mass)
// {
//     double dvDotdx = 0;
//     for (int d=0; d<dimensions; d++)
//     {
//         dvDotdx += (op.v[d] - nbr.v[d]) * (op.x[d] - nbr.x[d]);
//     }

//     if (dvDotdx < 0)
//     {
//         double viscousTerm = 2 * viscosityAlpha * smoothingLength * soundSpeed / (op.density + nbr.density);
//         double capitalPi = -viscousTerm * dvDotdx / (dist*dist + 0.01*smoothingLength*smoothingLength);
//         for (int d=0; d<dimensions; d++)
//         {
//             op.a[d] -= mass * capitalPi * gradW[d];
//         }
//     }
// }

// Move particles
// Semi-implicit Euler scheme
void SPH::stepParticles(double dt)
{
    double speedSquared;
    int nMax = particles.size();

    #pragma omp parallel for private(speedSquared)
    for (int n=numBoundaryParticles; n<nMax; n++)
    {
        Particle &op = particles[n];

        speedSquared = 0.0;
        for (int d=0; d<dimensions; d++)
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
            op.energyPerMass = max(op.energyPerMass, 0);

            // Handle boundary conditions
            switch (boundaryConditions)
            {
                case periodic:
                    if (op.x[d] < domainMin[d])
                    {
                        op.x[d] += (domainMax[d] - domainMin[d]);
                    }
                    else if (op.x[d] > domainMax[d])
                    {
                        op.x[d] -= (domainMax[d] - domainMin[d]);
                    }
                    break;

                case destructive:
                    if (op.x[d] < domainMin[d] || op.x[d] > domainMax[d])
                        op.isActive = false;
                    break;

                default:
                    break;
            }

            if (op.isActive)
            {
                minPosition[d] = min(minPosition[d], op.x[d]);
                maxPosition[d] = max(maxPosition[d], op.x[d]);
            }

            speedSquared += op.v[d]*op.v[d];
        }

        if (op.isActive && speedSquared > maxSpeedSquared)
        {
            #pragma omp atomic write
            maxSpeedSquared = speedSquared;
        }
    }

    if (boundaryConditions == destructive)
    {
        particles.erase(std::remove_if(particles.begin(),
                                       particles.end(),
                                       [](const Particle& p){return !p.isActive;}),
                    particles.end());
    }
}

double SPH::getCFLTimestep(double multiplier)
{
    return multiplier * fluidParticleSize / sqrt(maxSpeedSquared);
}

double WCSPH::getCFLTimestep(double multiplier)
{
    return multiplier * fluidParticleSize / soundSpeed;
}

// Resize grid and project particles to it
void SPH::setupParticleGrid()
{
    resizeGrid();
    projectParticlesToGrid();
}

// Initialise all particle positions and densities.
int SPH::initialiseParticles(const std::string& casefileName)
{
    std::ifstream caseFile (casefileName);
    int numFluidParticles, numParticles = 0;
    std::string line;
    std::stringstream linestream;
    double x;

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
            std::cout << "ERROR: No particles defined in case file." << std::endl;
            return 0;
        }

        std::getline(caseFile, line);
        linestream.clear();
        linestream.str(line);
        linestream >> domainMin[0] >> domainMin[1] >> domainMin[2];
        linestream >> domainMax[0] >> domainMax[1] >> domainMax[2];

        // Initialise particles and set positions
        for (int p=0; p<numParticles; p++)
        {
            Particle op = Particle();
            op.isBoundary = (p < numBoundaryParticles);
            op.density = restDensity;

            std::getline(caseFile, line);
            linestream.clear();
            linestream.str(line);
            for (int d=0; d<dimensions; d++)
            {
                linestream >> op.x[d];

                if (p == 0)
                {
                    minPosition[d] = op.x[d];
                    maxPosition[d] = op.x[d];
                }
                else
                {
                    if (op.x[d] < minPosition[d])
                        minPosition[d] = op.x[d];

                    if (op.x[d] > maxPosition[d])
                        maxPosition[d] = op.x[d];
                }
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
        std::cout << "ERROR: Unable to open case file." << std::endl;
        return 0;
    }

    std::cout << "\tFluid particles:        " << numFluidParticles << std::endl
              << "\tBoundary particles:     " << numBoundaryParticles << std::endl
              << "\tTotal particles:        " << numParticles << std::endl;

    return 1;
}

// Initialise array to hold relative coordinates of a grid cell neighbourhood
// to the central cell
void SPH::initCoordOffsetArray()
{
    coordOffsets = new int[gridNbrhoodSize*dimensions];
    for (int i=0; i<gridNbrhoodSize; i++)
    {
        for (int d=0; d<dimensions; d++)
        {
            int power = pow(3, dimensions-d-1);
            coordOffsets[coordOffsetsIndex(i, d, dimensions)] = ((i / power) % 3) - 1;
        }
    }
}

// Create grid for all particles.
void SPH::resizeGrid()
{
    int dim, gridSize = 1, dimFactor = 1;
    double cellSize = kernel->getSupportRadius(smoothingLength);
    bool needToResize = false;

    for (int d=0; d<dimensions; d++)
    {
        dim = ceil((maxPosition[d] - minPosition[d] + fluidParticleSize) / cellSize);
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
void SPH::projectParticlesToGrid()
{
    int cellNo, p = -1;
    double cellSize = kernel->getSupportRadius(smoothingLength);
    double halfParticleSize = fluidParticleSize / 2;

    for (auto& cell : grid)
    {
        cell.clear();
    }

    for (auto& op : particles)
    {
        ++p;
        cellNo = 0;
        for (int d=0; d<dimensions; d++)
        {
            cellNo += dimFactors[d] * (int)((op.x[d] - minPosition[d] + halfParticleSize) / cellSize);
        }

        grid.at(cellNo).push_back(p);
        op.gridCellNo = cellNo;
    }
}

void SPH::printParameters()
{
    std::string bcString;
    switch (boundaryConditions)
    {
        case periodic:
            bcString = "periodic";
            break;
        case destructive:
            bcString = "destructive";
            break;
    }

    std::cout << "STANDARD SPH PARAMETERS:" << std::endl
              << "\tDimensions:             " << dimensions << std::endl
              << "\tRest density:           " << restDensity << std::endl
              << "\tPressure constant:      " << pressureConstant << std::endl
              << "\tDynamic viscosity:      " << dynamicViscosity << std::endl
              << "\tExternal gravity:       " << extGravity << std::endl
              << "\tEquation of state:      " << *eos << std::endl
              << "\tFluid particle size:    " << fluidParticleSize << std::endl
              << "\tBoundary particle size: " << boundaryParticleSize << std::endl
              << "\tSmoothing length:       " << smoothingLength << std::endl
              << "\tSPH kernel:             " << *kernel << std::endl
              << "\tFluid particle mass:    " << fluidParticleMass << std::endl
              << "\tBoundary particle mass: " << boundaryParticleMass << std::endl
              << "\tBoundary conditions:    " << bcString << std::endl;
}

void WCSPH::printParameters()
{
    SPH::printParameters();

    std::cout << std::endl << "WCSPH PARAMETERS:" << std::endl
              << "\tDensity variation:      " << densityVariation << std::endl
              << "\tSound speed:            " << soundSpeed << std::endl;
}

void ThermoSPH::printParameters()
{
    SPH::printParameters();

    std::cout << std::endl << "ThermoSPH PARAMETERS:" << std::endl
              << "\tAdiabatic index:        " << adiabaticIndex << std::endl;
}

SPH::SPH(int d,
         double rho0,
         double eta,
         double g,
         SphKernel *kernel,
         double fluidSpacing,
         double boundarySpacing,
         double h,
         BoundaryConditions bc)
{
    dimensions = d;
    restDensity = rho0;
    dynamicViscosity = eta;
    extGravity = g;
    this->kernel = kernel;
    fluidParticleSize = fluidSpacing;
    boundaryParticleSize = boundarySpacing;
    smoothingLength = h;
    boundaryConditions = bc;

    gridNbrhoodSize = pow(3, dimensions);
    fluidParticleMass = restDensity * pow(fluidParticleSize, dimensions);
    boundaryParticleMass = restDensity * pow(boundaryParticleSize, dimensions);

    initCoordOffsetArray();
}

WCSPH::WCSPH(int d,
             double rho0,
             double eta,
             double g,
             SphKernel *kernel,
             double fluidSpacing,
             double boundarySpacing,
             double h,
             BoundaryConditions bc,
             double maxH,
             double rhoVar)
    : SPH(d, rho0, eta, g, kernel, fluidSpacing, boundarySpacing, h, bc)
{
    boundaryConditions = bc;
    densityVariation = rhoVar;

    double soundSpeedSquared = 2 * abs(g) * maxH / rhoVar;
    pressureConstant = rho0 * soundSpeedSquared / 7;
    soundSpeed = sqrt(soundSpeedSquared);

    this->eos = new WeaklyCompressibleEOS(pressureConstant, rho0);
}

ThermoSPH::ThermoSPH(int d,
                     double rho0,
                     double eta,
                     double g,
                     SphKernel *kernel,
                     double fluidSpacing,
                     double boundarySpacing,
                     double h,
                     BoundaryConditions bc,
                     double gamma)
    : SPH(d, rho0, eta, g, kernel, fluidSpacing, boundarySpacing, h, bc)
{
    adiabaticIndex = gamma;

    this->eos = new IdealGasEOS(gamma);
}

SPH::~SPH()
{
    delete eos;
    delete kernel;
    delete[] coordOffsets;
}
