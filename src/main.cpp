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
#include <memory>
#include <stdlib.h>

namespace Simulation
{
    std::string casefileName;

    struct TimestepData
    {
        double deltaT0;
        double tFinal;
        double deltaTPlot;
        double minDeltaT;
        double maxDeltaT;
        double cflLambda;
    } timestepData;

    SphSchemes::SphParams sphParams;
    std::unique_ptr<SphKernels::SphKernel> kernel;
    std::unique_ptr<SphSchemes::SPH> sphSystem;

    // Read and initialise all parameters from input file, with some validation.
    int configureSimulation(const std::string& paramfileName)
    {
        double adiabaticIndex, maxHeight, densityVariation;

        std::string line, kernelString, bcString, schemeString;
        std::ifstream paramFile (paramfileName);
        std::stringstream linestream;

        if (paramFile.is_open())
        {
            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.dimensions;
            if (sphParams.dimensions < 1 || sphParams.dimensions > 3)
            {
                std::cout << "ERROR: Invalid number of dimensions given: " << sphParams.dimensions << '\n';
                return 0;
            }

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.restDensity;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.dynamicViscosity;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.gravity;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.fpSize;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.bpSize;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> sphParams.smoothingLength;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> kernelString;
            if (kernelString.compare("cubicspline") == 0)
            {
                kernel = std::make_unique<SphKernels::CubicSplineKernel>(sphParams.dimensions);
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
                sphParams.boundaryConditions = SphSchemes::BoundaryConditions::periodic;
            }
            else if (bcString.compare("destructive") == 0)
            {
                sphParams.boundaryConditions = SphSchemes::BoundaryConditions::destructive;
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

                sphSystem = std::make_unique<SphSchemes::WCSPH>(sphParams,
                                                                std::move(kernel),
                                                                maxHeight,
                                                                densityVariation);
            }
            else if (schemeString.compare("thermosph") == 0)
            {
                std::getline(paramFile, line);
                linestream.str(line);
                linestream >> adiabaticIndex;

                sphSystem = std::make_unique<SphSchemes::ThermoSPH>(sphParams,
                                                                    std::move(kernel),
                                                                    adiabaticIndex);
            }
            else
            {
                std::cout << "ERROR: Unknown SPH scheme specified: \"" << schemeString << "\"" << '\n';
                return 0;
            }

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.deltaT0;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.minDeltaT;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.maxDeltaT;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.cflLambda;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.deltaTPlot;

            std::getline(paramFile, line);
            linestream.str(line);
            linestream >> timestepData.tFinal;

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
}

int main(int argc, char** argv)
{
    using Simulation::sphSystem;
    using Simulation::timestepData;

    const std::string caseDir = "./cases/";
    const std::string outputDir = "./output/";

    if (argc < 2)
    {
        std::cout << "ERROR: Expected 1 argument (paramfileName), got 0." << '\n';
        return 1;
    }

    if (Simulation::configureSimulation(argv[1]) == 0) return 1;

    SimOutput::VTKResultsWriter vtkout(outputDir);

    sphSystem->printParameters();
    std::cout << '\n' << "SIMULATION PARAMETERS:" << '\n'
              << "\tInitial timestep:       " << timestepData.deltaT0 << '\n'
              << "\tMinimum timestep:       " << timestepData.minDeltaT << '\n'
              << "\tMaximum timestep:       " << timestepData.maxDeltaT << '\n'
              << "\tCFL multiplier:         " << timestepData.cflLambda << '\n'
              << "\tPlot interval:          " << timestepData.deltaTPlot << '\n'
              << "\tFinal time:             " << timestepData.tFinal << '\n'
              << "\tCase file:              " << Simulation::casefileName << '\n';

    double t = 0.0;
    double tPlot = t + timestepData.deltaTPlot;
    double deltaT = timestepData.deltaT0;

    sphSystem->initialiseParticles(caseDir + Simulation::casefileName);

    sphSystem->setupParticleGrid();
    sphSystem->calcParticleDensities();

    // Plot initial state
    vtkout.writeSnapshot(sphSystem->particles);
    std::cout << '\n' << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSystem->maxSpeedSquared) << '\n';

    while (t <= timestepData.tFinal)
    {
        sphSystem->setupParticleGrid();
        sphSystem->calcParticleDensities();
        sphSystem->calcParticleForces();
        sphSystem->stepParticles(deltaT);

        if (timestepData.deltaTPlot > 0 && t >= tPlot)
        {
            // Plot state of the system
            vtkout.writeSnapshot(sphSystem->particles);
            std::cout << "Time: " << t << "\t"
                      << "max speed: " << sqrt(sphSystem->maxSpeedSquared) << '\n';
            tPlot += timestepData.deltaTPlot;
        }
        t += deltaT;

        double cflT = sphSystem->getCFLTimestep(timestepData.cflLambda);
        if (timestepData.deltaTPlot > 0)
        {
            deltaT = (tPlot - t < cflT) ? tPlot - t : cflT;
        }
        deltaT = (deltaT < timestepData.minDeltaT) ? timestepData.minDeltaT :
                 (timestepData.maxDeltaT > deltaT) ? timestepData.maxDeltaT : deltaT;
    }

    // Plot final state of the system
    vtkout.writeSnapshot(sphSystem->particles);
    std::cout << "Time: " << t << "\t"
              << "max speed: " << sqrt(sphSystem->maxSpeedSquared) << '\n';

}
