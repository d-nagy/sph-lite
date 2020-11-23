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

const std::string caseDir = "./cases/";
const std::string outputDir = "./output/";
std::string casefileName;

double t, tFinal, tPlot, deltaT, minDeltaT, maxDeltaT, deltaTPlot, cflLambda, cflT;
SPH *sphSim;

// Read and initialise all parameters from input file, with some validation.
int readParameters(const std::string& paramfileName)
{
    int dimensions;
    double fpSize, bpSize;
    double restDensity, dynamicViscosity, gravity, smoothingLength;
    double adiabaticIndex, maxHeight, densityVariation;
    BoundaryConditions boundaryConditions;
    SphKernel *kernel;

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
            std::cout << "ERROR: Invalid number of dimensions given: " << dimensions << '\n';
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
            std::cout << "ERROR: Unknown kernel specified: \"" << kernelString << "\"" << '\n';
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
            std::cout << "ERROR: Unknown boundary conditions specified: \"" << bcString << "\"" << '\n';
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
            std::cout << "ERROR: Unknown SPH scheme specified: \"" << schemeString << "\"" << '\n';
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
        std::cout << "All parameters read successfully." << '\n';
        return 1;
    }
    else
    {
        std::cout << "ERROR: Unable to open parameters file." << '\n';
        return 0;
    }
}

int main(int argc, char** argv)
{
    std::string paramfileName;

    if (argc < 2)
    {
        std::cout << "ERROR: Expected 1 argument (paramfileName), got 0." << '\n';
        return 1;
    }

    if (readParameters(argv[1]) == 0) return 1;

    VTKResultsWriter vtkout(outputDir);

    sphSim->printParameters();
    std::cout << '\n' << "SIMULATION PARAMETERS:" << '\n'
              << "\tInitial timestep:       " << deltaT << '\n'
              << "\tMinimum timestep:       " << minDeltaT << '\n'
              << "\tMaximum timestep:       " << maxDeltaT << '\n'
              << "\tCFL multiplier:         " << cflLambda << '\n'
              << "\tPlot interval:          " << deltaTPlot << '\n'
              << "\tFinal time:             " << tFinal << '\n'
              << "\tCase file:              " << casefileName << '\n';

    t = 0.0;
    tPlot = t + deltaTPlot;

    sphSim->initialiseParticles(caseDir + casefileName);

    sphSim->setupParticleGrid();
    sphSim->calcParticleDensities();

    // Plot initial state
    vtkout.writeSnapshot(sphSim->particles);
    std::cout << '\n' << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSim->maxSpeedSquared) << '\n';

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
                      << "max speed: " << sqrt(sphSim->maxSpeedSquared) << '\n';
            tPlot += deltaTPlot;
        }
        t += deltaT;

        cflT = sphSim->getCFLTimestep(cflLambda);
        if (deltaTPlot > 0)
        {
            deltaT = (tPlot - t < cflT) ? tPlot - t : cflT;
        }
        deltaT = (deltaT < minDeltaT) ? minDeltaT : (maxDeltaT > deltaT) ? maxDeltaT : deltaT;
    }

    // Plot final state of the system
    vtkout.writeSnapshot(sphSim->particles);
    std::cout << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSim->maxSpeedSquared) << '\n';

    delete sphSim;
}
