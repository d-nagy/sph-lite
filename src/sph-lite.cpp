#include "kernels.h"
#include "particle.h"
#include "eos.h"
#include "vtkout.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#define coordOffsetsIndex(row, col, dim) ((dim*row) + col)
#define max(x, y) ((x > y) ? x : y)
#define min(x, y) ((x > y) ? y : x)

const std::string inputFilename = "sph.dat";
const std::string caseDir = "./cases/";
const std::string outputDir = "./output/";

int dimensions;
double minPosition[3], maxPosition[3], maxSpeedSquared = 0.0;
int gridDims[3], dimFactors[3];
int neighbourhoodSize;
unsigned int numParticles;
double t, tFinal, tPlot, deltaT, minDeltaT, deltaTPlot, cflLambda;
double fpMass, bpMass, stiffness, restDensity, dynamicViscosity, gravity = -9.81, adiabaticIndex;
double neighbourRadius, smoothingLength, fpSize, bpSize;
int *coordOffsets;
std::ofstream videoFile;
std::string caseFilename;
std::vector<Particle> particles;
std::vector<std::vector<int>> grid;
SphKernel *kernel;
EquationOfState *eos;
VTKResultsWriter *vtkout;

// Initialise all particle positions and densities.
int initialiseParticles()
{
    std::ifstream caseFile (caseDir + caseFilename);
    int numFluidParticles, numBoundaryParticles;
    std::string line;
    std::stringstream linestream;

    if (caseFile.is_open())
    {
        std::getline(caseFile, line);
        linestream.str(line);
        if (linestream >> numBoundaryParticles >> numFluidParticles)
            numParticles = numBoundaryParticles + numFluidParticles;

        if (!numParticles)
        {
            std::cout << "No particles defined in case file." << std::endl;
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
                linestream >> op.pressure;
            }

            particles.push_back(op);
        }

        caseFile.close();
    }
    else
    {
        std::cout << "Unable to open case file." << std::endl;
        return 0;
    }

    std::cout << "Number of fluid particles: " << numFluidParticles << std::endl
         << "Number of boundary particles: " << numBoundaryParticles << std::endl
         << "Total number of particles: " << numParticles << std::endl;

    return 1;
}

// Release all array memory.
void cleanupObjects()
{
    delete[] coordOffsets;
    delete kernel;
    delete eos;
    delete vtkout;
}

// Initialise array to hold relative coordinates of a grid cell neighbourhood
// to the central cell
void initCoordOffsets()
{
    coordOffsets = new int[neighbourhoodSize*dimensions];
    for (int i=0; i<neighbourhoodSize; i++)
    {
        for (int d=0; d<dimensions; d++)
        {
            int power = pow(3, dimensions-d-1);
            coordOffsets[coordOffsetsIndex(i, d, dimensions)] = ((i / power) % 3) - 1;
        }
    }
}

// Create grid for all particles.
void prepareGrid(double const minP[], double const maxP[], double cellLength)
{
    int dim, gridSize = 1, dimFactor = 1;

    grid.clear();

    for (int d=0; d<dimensions; d++)
    {
        dim = ceil((maxP[d] - minP[d] + fpSize) / cellLength);
        gridSize *= dim;
        gridDims[d] = dim;
        dimFactors[d] = dimFactor;
        dimFactor *= dim;
    }

    grid.resize(gridSize);
}

// Assign each particle to a grid cell.
void projectParticlesToGrid(double const minP[], double cellLength)
{
    int cellNo, p = -1;
    for (auto &op : particles)
    {
        p++;
        cellNo = 0;
        for (int d=0; d<dimensions; d++)
        {
            cellNo += dimFactors[d] * (int)((op.x[d] - minP[d]) / cellLength);
        }

        grid.at(cellNo).push_back(p);
        op.gridCellNo = cellNo;
    }
}

// Calculate position in semi-flattened grid vector from a
// coordinate array.
int getGridFlatIndex(int const coordArray[], int const offsetArray[])
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

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void calcParticleDensities()
{
    int origin[dimensions];
    double densitySum, mass, dist;
    int cellNo, p = -1, c, df;

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
        for (int i=0; i<neighbourhoodSize; i++)
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
                mass = (nbr.isBoundary) ? bpMass : fpMass;

                // Calculate distance to sample/origin particle
                dist = 0;
                for (int d=0; d<dimensions; d++)
                {
                    dist += (op.x[d] - nbr.x[d]) * (op.x[d] - nbr.x[d]);
                }
                dist = sqrt(dist);

                // Check if neighbour
                if (dist <= neighbourRadius)
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
void calcParticleForces()
{
    double gradWNorm, gradWComponent, temperatureDerivativeSum;
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

        /* op.temperatureDerivative = 0; */
        opPOverRhoSquared = op.pressure / (op.density * op.density);

        int i=0;
        for (auto pNbr : op.neighbours)
        {
            Particle &nbr = particles.at(pNbr);
            dist = op.neighbourDist[i++];

            if (&op == &nbr) continue;

            gradWNorm = 0;
            mass = (nbr.isBoundary) ? bpMass : fpMass;
            nbrPressure = (nbr.isBoundary) ? op.pressure : nbr.pressure;
            nbrPOverRhoSquared = nbrPressure / (nbr.density * nbr.density); // Should nbr.density be op.density for boundary nbrs?

            // Add acceleration from pressure and calculate kernel gradient
            for (int d=0; d<dimensions; d++)
            {
                gradWComponent = kernel->gradWComponent(op.x[d], nbr.x[d], dist, smoothingLength);
                gradWNorm += gradWComponent * gradWComponent;
                op.a[d] -= mass * (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;
                /* op.temperatureDerivative += mass * (opPOverRhoSquared + nbrPOverRhoSquared) * (op.v[d] - nbr.v[d]) * gradWComponent / 2; */
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
void updateParticles()
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
                op.a[d] += gravity;

            // Update velocities and positions
            // Semi-implicit Euler scheme
            op.v[d] += deltaT * op.a[d];
            op.x[d] += deltaT * op.v[d];

            // Update temperature
            /* op.temperature += deltaT * op.temperatureDerivative; */

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

// Update timestep size according to the CFL condition
void updateTimestepSize()
{
    double cflT = cflLambda * fpSize / sqrt(maxSpeedSquared);
    deltaT = min(tPlot - t, cflT);
    deltaT = max(minDeltaT, deltaT);
}

int readParameters()
{
    std::string line;
    std::ifstream inputFile (inputFilename);
    std::stringstream linestream;
    if (inputFile.is_open())
    {
        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> dimensions;
        std::cout << "Dimensions: " << dimensions << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> restDensity;
        std::cout << "Rest density: " << restDensity << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> dynamicViscosity;
        std::cout << "Dynamic viscosity: " << dynamicViscosity << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> stiffness;
        std::cout << "Stiffness coefficient: " << stiffness << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> fpSize;
        std::cout << "Fluid particle size: " << fpSize << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> bpSize;
        std::cout << "Boundary particle size: " << bpSize << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> smoothingLength;
        std::cout << "Smoothing length: " << smoothingLength << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> deltaT;
        std::cout << "Timestep: " << deltaT << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> minDeltaT;
        std::cout << "Minimum timestep: " << minDeltaT << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> cflLambda;
        std::cout << "CFL multiplier: " << cflLambda << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> deltaTPlot;
        std::cout << "Plot interval: " << deltaTPlot << std::endl;

        std::getline(inputFile, line);
        linestream.str(line);
        linestream >> tFinal;
        std::cout << "Final time: " << tFinal << std::endl;

        std::getline(inputFile, line);
        caseFilename = line;
        std::cout << "Case file: " << caseFilename << std::endl;

        inputFile.close();
        return 1;
    }
    else
    {
        std::cout << "Unable to open parameters file." << std::endl;
        return 0;
    }
}

int main(int argc, char** argv)
{
    if (readParameters())
    {
        std::cout << "All parameters read successfully." << std::endl;
    }
    else return 1;

    kernel = new CubicSplineKernel(dimensions);
    eos = new WeaklyCompressibleEOS(stiffness, restDensity);
    vtkout = new VTKResultsWriter(outputDir);

    t = 0.0;
    tPlot = t + deltaTPlot;
    neighbourhoodSize = pow(3, dimensions);
    neighbourRadius = 2 * smoothingLength;

    // Set particle masses
    fpMass = restDensity * pow(fpSize, dimensions);
    bpMass = restDensity * pow(bpSize, dimensions);

    std::cout << "Fluid particle mass: " << fpMass << std::endl
              << "Boundary particle mass: " << bpMass << std::endl;

    initialiseParticles();
    initCoordOffsets();

    prepareGrid(minPosition, maxPosition, neighbourRadius);
    projectParticlesToGrid(minPosition, neighbourRadius);
    calcParticleDensities();

    // Plot initial state
    vtkout->printSnapshot(particles);
    std::cout << "Time: " << t << "\t"
              << "deltaT: " << deltaT << "\t"
              << "max speed: " << sqrt(maxSpeedSquared) << std::endl;

    while (t <= tFinal)
    {
        prepareGrid(minPosition, maxPosition, neighbourRadius);
        projectParticlesToGrid(minPosition, neighbourRadius);
        calcParticleDensities();
        calcParticleForces();
        updateParticles();

        if (t >= tPlot)
        {
            // Plot state of the system
            vtkout->printSnapshot(particles);
            std::cout << "Time: " << t << "\t"
                      << "deltaT: " << deltaT << "\t"
                      << "max speed: " << sqrt(maxSpeedSquared) << std::endl;
            tPlot += deltaTPlot;
        }
        t += deltaT;

        updateTimestepSize();
    }

    // Plot final state of the system
    vtkout->printSnapshot(particles);
    std::cout << "Time: " << t << "\t"
              << "deltaT: " << deltaT << "\t"
              << "max speed: " << sqrt(maxSpeedSquared) << std::endl;

    // Cleanup
    cleanupObjects();

    // Testing/troubleshooting actions
    /* plotCubicKernel(smoothingLength); */
}
