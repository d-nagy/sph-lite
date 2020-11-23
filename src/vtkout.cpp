#include "vtkout.h"
#include "particle.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>

VTKResultsWriter::VTKResultsWriter(const std::string& outputDir) : outputDir(outputDir)
{
    videoFile.open(outputDir + "result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << '\n'
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << '\n'
              << "<Collection>";
}

VTKResultsWriter::~VTKResultsWriter()
{
    videoFile << "</Collection>"
              << "</VTKFile>" << '\n';
    videoFile.close();
}

void VTKResultsWriter::writeSnapshot(const std::vector<Particle>& ps)
{
    static int counter = -1;
    std::stringstream filename;
    filename << "result-" << ++counter <<  ".vtp";

    std::ofstream out(outputDir + filename.str().c_str());

    out << "<VTKFile type=\"PolyData\">" << '\n'
        << "<PolyData>" << '\n'
        << " <Piece NumberOfPoints=\"" << ps.size() << "\">" << '\n'
        << "  <Points>" << '\n'
        << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        // << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    // Output particle positions
    for (auto &p : ps)
    {
        for (int d=0; d<3; d++)
        {
            out << p.x[d] << " ";
        }
    }

    out << "   </DataArray>" << '\n'
        << "  </Points>" << '\n'
        << "  <PointData>" << '\n'
        << "   <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">" << '\n';

    // Output particle densities
    for (auto &p : ps)
    {
        out << p.density << " ";
    }

    out << "   </DataArray>" << '\n'
        << "   <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">" << '\n';

    // Output particle pressures
    for (auto &p : ps)
    {
        out << p.pressure << " ";
    }

    out << "   </DataArray>" << '\n'
        << "   <DataArray type=\"Float64\" Name=\"EnergyPerMass\" format=\"ascii\">" << '\n';

    // Output particle internal energy per mass
    for (auto &p : ps)
    {
        out << p.energyPerMass << " ";
    }

    out << "   </DataArray>" << '\n'
        << "   <DataArray type=\"Int32\" Name=\"CellNumber\" format=\"ascii\">" << '\n';

    // Output particle internal energy per mass
    for (auto &p : ps)
    {
        out << p.gridCellNo % 10 << " ";
    }

    out << "   </DataArray>" << '\n'
        << "  </PointData>" << '\n'
        << " </Piece>" << '\n'
        << "</PolyData>" << '\n'
        << "</VTKFile>"  << '\n';

    out.close();

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << '\n';
}
