// VTK Output class
#ifndef VTKOUT_H
#define VTKOUT_H

#include "particle.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>

namespace SimOutput
{
    class VTKResultsWriter
    {
        public:
            VTKResultsWriter(const std::string& outputDir);
            ~VTKResultsWriter();
            void writeSnapshot(const std::vector<SphSchemes::Particle>& ps);

        private:
            std::string outputDir;
            std::ofstream videoFile;
    };
}

#endif
