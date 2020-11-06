#include "sph.h"
#include "particle.h"

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>

#define coordOffsetsIndex(row, col, dim) ((dim*row) + col)
#define max(x, y) ((x > y) ? x : y)
#define min(x, y) ((x > y) ? y : x)

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void SPH::calcParticleDensities()
{
    int origin[dimensions];
    double densitySum, mass, dist;
    int cellNo, p = -1, c, df;
    double nbrRadius = kernel->getSupportRadius(smoothingLength);

    for (auto &op : particles)
    {
        p++;
        if (op.isBoundary) continue;

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
        for (int i=0; i<pow(3, dimensions); i++)
        {
            // Calculate the cell index in the semi-flattened grid vector.
            cellNo = getGridFlatIndex(origin, &coordOffsets[coordOffsetsIndex(i, 0, dimensions)]);

            if (cellNo < 0)
                continue;

            // Loop through all particles in the grid cell
            std::vector<int> &cell = grid.at(cellNo);
            for (auto pNbr : cell)
            {
                Particle &nbr = particles.at(pNbr);
                mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;

                // Calculate distance to sample/origin particle
                dist = 0;
                for (int d=0; d<dimensions; d++)
                {
                    dist += (op.x[d] - nbr.x[d]) * (op.x[d] - nbr.x[d]);
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

        op.density = max(restDensity, densitySum);
        op.pressure = eos->getPressure(op);
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
    double gradWNorm, gradWComponent, energyPerMassDerivativeSum;
    double mass, dist, opPOverRhoSquared, nbrPOverRhoSquared, nbrPressure;

    int p = -1;
    for (auto &op : particles)
    {
        p++;
        if (op.isBoundary) continue;

        // Reset particle accelerations
        for (int d=0; d<dimensions; d++)
        {
            op.a[d] = 0.0;
        }

        /* op.energyPerMassDerivative = 0; */
        opPOverRhoSquared = op.pressure / (op.density * op.density);

        int i=0;
        for (auto pNbr : op.neighbours)
        {
            Particle &nbr = particles.at(pNbr);
            dist = op.neighbourDist[i++];

            if (&op == &nbr) continue;

            gradWNorm = 0;
            mass = (nbr.isBoundary) ? boundaryParticleMass : fluidParticleMass;
            nbrPressure = (nbr.isBoundary) ? op.pressure : nbr.pressure;
            nbrPOverRhoSquared = nbrPressure / (nbr.density * nbr.density); // Should nbr.density be op.density for boundary nbrs?

            // Add acceleration from pressure and calculate kernel gradient
            for (int d=0; d<dimensions; d++)
            {
                gradWComponent = kernel->gradWComponent(op.x[d], nbr.x[d], dist, smoothingLength);
                gradWNorm += gradWComponent * gradWComponent;
                op.a[d] -= mass * (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;
                /* op.energyPerMassDerivative += mass * (opPOverRhoSquared + nbrPOverRhoSquared) * (op.v[d] - nbr.v[d]) * gradWComponent / 2; */
            }

            gradWNorm = sqrt(gradWNorm);

            // Add acceleration from viscosity
            for (int d=0; d<dimensions; d++)
            {
                op.a[d] += (mass * 2 * gradWNorm * dynamicViscosity) * (nbr.v[d] - op.v[d]) / (op.density * nbr.density * dist);
            }
        }
    }
}

// Move particles
// Semi-implicit Euler scheme
void SPH::stepParticles(double dt)
{
    double speedSquared;
    int p = -1;
    for (auto &op : particles)
    {
        p++;
        if (op.isBoundary) continue;

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

            // Update temperature
            /* op.energyPerMass += deltaT * op.energyPerMassDerivative; */

            if (op.x[d] < minPosition[d])
                minPosition[d] = op.x[d];

            if (op.x[d] > maxPosition[d])
                maxPosition[d] = op.x[d];

            speedSquared += op.v[d]*op.v[d];
        }

        if (speedSquared > maxSpeedSquared)
            maxSpeedSquared = speedSquared;
    }
}

double SPH::getCFLTimestep(double multiplier)
{
    return multiplier * fluidParticleSize / sqrt(maxSpeedSquared);
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
    int numFluidParticles, numBoundaryParticles, numParticles = 0;
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
            std::cout << "ERROR: No particles defined in case file." << std::endl;
            return 0;
        }

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

    std::cout << "Number of fluid particles: " << numFluidParticles << std::endl
              << "Number of boundary particles: " << numBoundaryParticles << std::endl
              << "Total number of particles: " << numParticles << std::endl;

    return 1;
}

// Initialise array to hold relative coordinates of a grid cell neighbourhood
// to the central cell
void SPH::initCoordOffsetArray()
{
    int nsize = pow(3, dimensions);
    coordOffsets = new int[nsize*dimensions];
    for (int i=0; i<nsize; i++)
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

    grid.clear();

    for (int d=0; d<dimensions; d++)
    {
        dim = ceil((maxPosition[d] - minPosition[d] + fluidParticleSize) / cellSize);
        gridSize *= dim;
        gridDims[d] = dim;
        dimFactors[d] = dimFactor;
        dimFactor *= dim;
    }

    grid.resize(gridSize);
}

// Assign each particle to a grid cell.
void SPH::projectParticlesToGrid()
{
    int cellNo, p = -1;
    double cellSize = kernel->getSupportRadius(smoothingLength);

    for (auto &op : particles)
    {
        p++;
        cellNo = 0;
        for (int d=0; d<dimensions; d++)
        {
            cellNo += dimFactors[d] * (int)((op.x[d] - minPosition[d]) / cellSize);
        }

        grid.at(cellNo).push_back(p);
        op.gridCellNo = cellNo;
    }
}

// Calculate position in semi-flattened grid vector from a
// coordinate array.
int SPH::getGridFlatIndex(int const coordArray[], int const offsetArray[])
{
    int c, cellNo = 0;
    bool valid_coords = true;

    for (int d=0; d<dimensions; d++)
    {
        c = coordArray[d] + offsetArray[d];
        if (c < 0 || c >= gridDims[d]) { return -1; }
        cellNo += (dimFactors[d] * c);
    }

    return cellNo;
}

void SPH::printParameters()
{
    std::cout << "Dimensions: " << dimensions << std::endl
              << "Rest density: " << restDensity << std::endl
              << "Dynamic viscosity: " << dynamicViscosity << std::endl
              << "Stiffness coefficient: " << stiffness << std::endl
              << "Adiabatic index: " << adiabaticIndex << std::endl
              << "External gravitational acceleration: " << extGravity << std::endl
              << "Equation of state: " << *eos << std::endl
              << "Fluid particle size: " << fluidParticleSize << std::endl
              << "Boundary particle size: " << boundaryParticleSize << std::endl
              << "Smoothing length: " << smoothingLength << std::endl
              << "SPH kernel: " << *kernel << std::endl
              << "Fluid particle mass: " << fluidParticleMass << std::endl
              << "Boundary particle mass: " << boundaryParticleMass << std::endl;
}

SPH::SPH(int d,
         double rd,
         double dv,
         double s,
         double ai,
         double g,
         EquationOfState *eos,
         SphKernel *kernel,
         double fps,
         double bps,
         double sl)
{
    dimensions = d;
    restDensity = rd;
    dynamicViscosity = dv;
    stiffness = s;
    adiabaticIndex = ai;
    extGravity = g;
    this->eos = eos;
    this->kernel = kernel;
    fluidParticleSize = fps;
    boundaryParticleSize = bps;
    smoothingLength = sl;

    fluidParticleMass = restDensity * pow(fluidParticleSize, dimensions);
    boundaryParticleMass = restDensity * pow(boundaryParticleMass, dimensions);

    initCoordOffsetArray();
}

SPH::~SPH()
{
    delete eos;
    delete kernel;
    delete[] coordOffsets;
}
