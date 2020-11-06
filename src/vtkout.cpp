#include "vtkout.h"
#include "particle.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

VTKResultsWriter::VTKResultsWriter(const std::string& outputDir) : outputDir(outputDir)
{
    videoFile.open(outputDir + "result.pvd");
    videoFile << "<?xml version=\"1.0\"?>" << std::endl
              << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
              << "<Collection>";
}

VTKResultsWriter::~VTKResultsWriter()
{
    videoFile << "</Collection>"
              << "</VTKFile>" << std::endl;
    videoFile.close();
}

void VTKResultsWriter::writeSnapshot(const std::vector<Particle>& ps)
{
    static int counter = -1;
    std::stringstream filename;
    std::ofstream out(outputDir + filename.str().c_str());

    filename << "result-" << ++counter <<  ".vtp";

    out << "<VTKFile type=\"PolyData\">" << std::endl
        << "<PolyData>" << std::endl
        << " <Piece NumberOfPoints=\"" << ps.size() << "\">" << std::endl
        << "  <Points>" << std::endl
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

    out << "   </DataArray>" << std::endl
        << "  </Points>" << std::endl
        << "  <PointData>" << std::endl
        << "   <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">" << std::endl;

    // Output particle densities
    for (auto &p : ps)
    {
        out << p.density << " ";
    }

    out << "   </DataArray>" << std::endl
        << "   <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">" << std::endl;

    // Output particle pressures
    for (auto &p : ps)
    {
        out << p.pressure << " ";
    }

    out << "   </DataArray>" << std::endl
        << "  </PointData>" << std::endl
        << " </Piece>" << std::endl
        << "</PolyData>" << std::endl
        << "</VTKFile>"  << std::endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}
