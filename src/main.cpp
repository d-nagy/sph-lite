#include "kernels.h"
#include "particle.h"
#include "eos.h"
#include "vtkout.h"
#include "sph.h"

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <stdlib.h>

#define max(x, y) ((x > y) ? x : y)
#define min(x, y) ((x > y) ? y : x)

const std::string caseDir = "./cases/";
const std::string outputDir = "./output/";
std::string casefileName;

int dimensions;
double t, tFinal, tPlot, deltaT, minDeltaT, maxDeltaT, deltaTPlot, cflLambda, cflT;
double fpMass, bpMass, fpSize, bpSize;
double restDensity, dynamicViscosity, gravity, smoothingLength;
double adiabaticIndex, maxHeight, densityVariation;
BoundaryConditions boundaryConditions;
SphKernel *kernel;
SPH *sphSim;

// Read and initialise all parameters from input file, with some validation.
int readParameters(const std::string& paramfileName)
{
    std::string line, kernelString, bcString, schemeString;
    std::ifstream paramFile (paramfileName);
    std::stringstream linestream;

    if (paramFile.is_open())
    {
        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> dimensions;
        if (dimensions < 1 || dimensions > 3)
        {
            std::cout << "ERROR: Invalid number of dimensions given: " << dimensions << std::endl;
            return 0;
        }

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> restDensity;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> dynamicViscosity;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> gravity;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> fpSize;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> bpSize;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> smoothingLength;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> kernelString;
        if (kernelString.compare("cubicspline") == 0)
        {
            kernel = new CubicSplineKernel(dimensions);
        }
        else
        {
            std::cout << "ERROR: Unknown kernel specified: \"" << kernelString << "\"" << std::endl;
            return 0;
        }

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> bcString;
        if (bcString.compare("periodic") == 0)
        {
            boundaryConditions = periodic;
        }
        else if (bcString.compare("destructive") == 0)
        {
            boundaryConditions = destructive;
        }
        else
        {
            std::cout << "ERROR: Unknown boundary conditions specified: \"" << bcString << "\"" << std::endl;
            return 0;
        }

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> schemeString;
        if (schemeString.compare("wcsph") == 0)
        {
            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> maxHeight;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> densityVariation;

            sphSim = new WCSPH(dimensions,
                               restDensity,
                               dynamicViscosity,
                               gravity,
                               kernel,
                               fpSize,
                               bpSize,
                               smoothingLength,
                               boundaryConditions,
                               maxHeight,
                               densityVariation);
        }
        else if (schemeString.compare("thermosph") == 0)
        {
            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> adiabaticIndex;

            sphSim = new ThermoSPH(dimensions,
                                   restDensity,
                                   dynamicViscosity,
                                   gravity,
                                   kernel,
                                   fpSize,
                                   bpSize,
                                   smoothingLength,
                                   boundaryConditions,
                                   adiabaticIndex);
        }
        else
        {
            std::cout << "ERROR: Unknown SPH scheme specified: \"" << schemeString << "\"" << std::endl;
            return 0;
        }

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> deltaT;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> minDeltaT;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> maxDeltaT;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> cflLambda;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> deltaTPlot;

        std::getline(paramFile, line);
        linestream.str(line);
        linestream >> tFinal;

        std::getline(paramFile, line);
        casefileName = line;

        paramFile.close();
        std::cout << "All parameters read successfully." << std::endl;
        return 1;
    }
    else
    {
        std::cout << "ERROR: Unable to open parameters file." << std::endl;
        return 0;
    }
}

int main(int argc, char** argv)
{
    std::string paramfileName;

    if (argc < 2)
    {
        std::cout << "ERROR: Expected 1 argument (paramfileName), got 0." << std::endl;
        return 1;
    }

    if (readParameters(argv[1]) == 0) return 1;

    VTKResultsWriter vtkout(outputDir);

    sphSim->printParameters();
    std::cout << std::endl << "SIMULATION PARAMETERS:" << std::endl
              << "\tInitial timestep:       " << deltaT << std::endl
              << "\tMinimum timestep:       " << minDeltaT << std::endl
              << "\tMaximum timestep:       " << maxDeltaT << std::endl
              << "\tCFL multiplier:         " << cflLambda << std::endl
              << "\tPlot interval:          " << deltaTPlot << std::endl
              << "\tFinal time:             " << tFinal << std::endl
              << "\tCase file:              " << casefileName << std::endl;

    t = 0.0;
    tPlot = t + deltaTPlot;

    sphSim->initialiseParticles(caseDir + casefileName);

    sphSim->setupParticleGrid();
    sphSim->calcParticleDensities();

    // Plot initial state
    vtkout.writeSnapshot(sphSim->particles);
    std::cout << std::endl << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSim->maxSpeedSquared) << std::endl;

    while (t <= tFinal)
    {
        sphSim->setupParticleGrid();
        sphSim->calcParticleDensities();
        sphSim->calcParticleForces();
        sphSim->stepParticles(deltaT);

        if (deltaTPlot > 0 && t >= tPlot)
        {
            // Plot state of the system
            vtkout.writeSnapshot(sphSim->particles);
            std::cout << "Time: " << t << "\t"
                      << "max speed: " << sqrt(sphSim->maxSpeedSquared) << std::endl;
            tPlot += deltaTPlot;
        }
        t += deltaT;

        cflT = sphSim->getCFLTimestep(cflLambda);
        deltaT = (deltaTPlot > 0) ? min(tPlot - t, cflT) : cflT;
        deltaT = max(minDeltaT, deltaT);
        deltaT = min(maxDeltaT, deltaT);
    }

    // Plot final state of the system
    vtkout.writeSnapshot(sphSim->particles);
    std::cout << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSim->maxSpeedSquared) << std::endl;

    delete sphSim;
}
