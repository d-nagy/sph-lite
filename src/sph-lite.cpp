#define _USE_MATH_DEFINES

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

using namespace std;

const double cubicKernelWNormalFactors[3] = {2.0/3.0, 10.0/(7*M_PI), 1/M_PI};
const string inputFilename = "sph.dat";
const string caseDir = "./cases/";
const string outputDir = "./output/";

int dimensions;

double *minPosition, *maxPosition;
int *gridDims, *dimFactors;
int neighbourhoodSize;

unsigned int numParticles;
double t, tFinal, tPlot, deltaT, deltaTPlot;
double fpMass, bpMass, stiffness, restDensity, dynamicViscosity, gravity = -9.81;
double neighbourRadius, smoothingLength, fpSize, bpSize;
int *coordOffsets;
ofstream videoFile;
string caseFilename;

class Particle
{
    public:
        bool isBoundary = false;
        double *x, *v, *a;
        double *accelPressure, *laplacianVelocity;
        double pressure, density;
        vector<int> neighbours;
        vector<double> neighbourDist;
        unsigned int gridCellNo;
        Particle() {};
        Particle(int dims)
        {
            x = new double[dims];
            v = new double[dims];
            a = new double[dims];
            accelPressure = new double[dims];
            laplacianVelocity = new double[dims];
        }
};

vector<Particle> particles;
vector<vector<int>> grid;

// Cubic spline kernel function with compact support of size 2h
double cubicKernelW(double r, double h)
{
   double q = r/h;
   double result = cubicKernelWNormalFactors[dimensions-1] / pow(h, (double)dimensions);
   
   if (q >= 0 && q < 1)
   {
       result *= ((2 - q)*(2 - q)*(2 - q)/4.0) - ((1 - q)*(1 - q)*(1 - q));
   } 
   else if (q >= 1 && q < 2)
   {
       result *= ((2 - q)*(2 - q)*(2 - q)/4.0);
   }
   else if (q >= 2)
   {
       result *= 0;
   }

   return result;
}

// First derivative of the cubic spline kernel w.r.t. q = r/h
double cubicKernelDelW(double r, double h)
{
    double q = r/h;
    double result = cubicKernelWNormalFactors[dimensions-1] / pow(h, (double)dimensions);

    if (q >= 0 && q < 1)
    {
        result *= (3 * (((2 - q)*(2 - q)/-4.0) + ((1 - q)*(1 - q))));
    } 
    else if (q >= 1 && q < 2)
    {
        result *= (3 * ((2 - q)*(2 - q)/-4.0));
    }
    else if (q >= 2)
    {
        result *= 0;
    }

    return result;
}

// Calculate a single component of the cubic kernel gradient with respect to
// the position of a particle
inline double cubicKernelGradComponent(double x1, double x2, double r, double h)
{
    return cubicKernelDelW(r, h) * (x1 - x2) / (r * h);
}

// Plot kernel function and first derivative for testing
void plotCubicKernel(double h)
{
    cout << "h=" << h << endl;
    cout << "q," << "W," << "dW" << endl;

    double maxR = 2*h;
    double dq = 2*h / 100;
    maxR += dq;
    for (double r=0.0; r<maxR; r+=dq)
    {
        cout << r/h << "," << cubicKernelW(r, h) << "," << cubicKernelDelW(r, h) << endl;
    }
}

// Return pressure value as a function of density according to a state equation
/* inline double pressureStateEquation(double rho) */
/* { */
/*     return stiffness * (rho - restDensity); */
/* } */

/* inline double pressureStateEquation(double rho) */
/* { */
/*     return stiffness * ((rho / restDensity) - 1); */
/* } */

inline double pressureStateEquation(double rho)
{
    return stiffness * (pow((rho/restDensity), 7) - 1);
}

// Initialise all particle positions and densities.
int initialiseParticles()
{
    ifstream caseFile (caseDir + caseFilename);
    int numFluidParticles, numBoundaryParticles;
    string line;
    stringstream linestream;

    if (caseFile.is_open())
    {
        getline(caseFile, line);
        linestream.str(line);
        if (linestream >> numBoundaryParticles >> numFluidParticles)
            numParticles = numBoundaryParticles + numFluidParticles;
        
        if (!numParticles)
        {
            cout << "No particles defined in case file." << endl;
            return 0;
        }

        // Initialise particles and set positions
        for (int p=0; p<numParticles; p++)
        {
            Particle op = Particle(dimensions);
            op.isBoundary = (p < numBoundaryParticles);
            op.density = restDensity;

            getline(caseFile, line);
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
        cout << "Unable to open case file." << endl;
        return 0;
    }

    cout << "Number of fluid particles: " << numFluidParticles << endl
         << "Number of boundary particles: " << numBoundaryParticles << endl
         << "Total number of particles: " << numParticles << endl;

    return 1;
}

// Release all array memory.
void releaseParticles()
{
    delete[] coordOffsets;
    delete[] minPosition;
    delete[] maxPosition;
    delete[] gridDims;
    delete[] dimFactors;
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
void prepareGrid(double *minP, double *maxP, double cellLength)
{
    grid.clear();

    int gridSize = 1;
    int dimFactor = 1;
    for (int d=0; d<dimensions; d++)
    {
        int dim = ceil((maxP[d] - minP[d] + fpSize) / cellLength);
        gridSize *= dim;
        gridDims[d] = dim;
        dimFactors[d] = dimFactor;
        dimFactor *= dim;
    }
    
    grid.resize(gridSize);
}

// Assign each particle to a grid cell.
void projectParticlesToGrid(double *minP, double cellLength)
{
    int p = 0;
    for (auto &op : particles)
    {
        int cellNo = 0;
        for (int d=0; d<dimensions; d++)
        {
            cellNo += dimFactors[d] * (int)((op.x[d] - minP[d]) / cellLength);
        }

        if (cellNo >= grid.size())
        {
            cout << "Invalid grid cell number calculated." << endl;
            exit(EXIT_FAILURE);
        }
        grid[cellNo].push_back(p);
        op.gridCellNo = cellNo;
        p++;
    }
}

// Calculate a coordinate array from an index in the semi-flattened
// grid vector.
int* getGridCoordArray(int cellNo)
{
    static int coords[3];
    for (int d=dimensions; d>=0; d--)
    {
        int df = dimFactors[d];
        int c = cellNo / df;
        coords[d] = c;
        cellNo -= (c * df);
    }
    
    return coords;
}

// Calculate position in semi-flattened grid vector from a
// coordinate array.
int getGridFlatIndex(int* coordArray, int* offsetArray)
{
    int cellNo = 0;
    bool valid_coords = true;
    for (int d=0; d<dimensions; d++)
    {
        int c = coordArray[d] + offsetArray[d];
        if (c < 0 || c >= gridDims[d])
        {
            return -1;
        }
        cellNo += (dimFactors[d] * c); 
    }
    return cellNo;
}

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void calcParticleDensities()
{
    // Loop through all particles
    int p = -1;
    for (auto &op : particles)
    {
        p++;
        if (op.isBoundary) continue;

        // Clear neighbour lists
        op.neighbours.clear();
        op.neighbourDist.clear();

        // Calculate this particle's grid coordinates
        int* origin = getGridCoordArray(op.gridCellNo);

        // Loop through all the neighbourhood grid cells
        /* double densitySum = fpMass * cubicKernelW(0, smoothingLength); // include own contribution */
        double densitySum = 0.0;
        for (int i=0; i<neighbourhoodSize; i++)
        {
            // Calculate the cell index in the semi-flattened grid vector.
            int cellNo = getGridFlatIndex(origin, &coordOffsets[coordOffsetsIndex(i, 0, dimensions)]); 

            if (cellNo < 0)
                continue;

            // Loop through all particles in the grid cell
            vector<int> cell = grid[cellNo];
            for (auto pNbr : cell)
            {
                Particle *nbr = &particles[pNbr];
                double mass = (nbr->isBoundary) ? bpMass : fpMass; 
              
                // Calculate distance to sample/origin particle
                double dist = 0;
                for (int d=0; d<dimensions; d++)
                {
                    dist += (op.x[d] - nbr->x[d]) * (op.x[d] - nbr->x[d]);
                }
                dist = sqrt(dist);

                // Check if neighbour
                if (dist <= neighbourRadius)
                {
                    op.neighbours.push_back(pNbr);
                    op.neighbourDist.push_back(dist);
                    densitySum += (mass * cubicKernelW(dist, smoothingLength));
                }
            }
        }

        op.density = max(restDensity, densitySum);
        op.pressure = pressureStateEquation(densitySum);
    }
}

// Calculate forces on particles from viscosity, pressure
// and external forces
//
// Viscosity term requiring the Laplacian of the velocity
// is computed using Brookshaw's discrete Laplacian for SPH.
//
// Pressure term requiring the gradient of the pressure
// is computed using the symmetric SPH formula.
void calcParticleForces()
{
    for (auto &op : particles)
    {
        if (op.isBoundary) continue;

        // Reset accelPressure and laplacianVelocity arrays
        for (int d=0; d<dimensions; d++)
        {
            op.accelPressure[d] = 0.0;
            op.laplacianVelocity[d] = 0.0;
        }

        double opRho = op.density, opP = op.pressure;
        double opPOverRhoSquared = opP / (opRho * opRho);

        int i=0;
        for (auto pNbr : op.neighbours)
        {
            Particle *nbr = &particles[pNbr];
            double dist = op.neighbourDist[i++];
            
            if (&op == nbr) continue;

            // calculate pressure force per spatial axis
            // as well as kernel gradient
            double gradW[dimensions];
            double gradWNorm = 0;
            double nbrPressure = (nbr->isBoundary) ? op.pressure : nbr->pressure;
            /* double nbrDensity = (nbr->isBoundary) ? op->density : nbr->density; */
            double mass = (nbr->isBoundary) ? bpMass : fpMass;
            double nbrPOverRhoSquared = nbrPressure / (nbr->density * nbr->density); // Should nbr->density be op->density for boundary nbrs?
            for (int d=0; d<dimensions; d++)
            {
                double gradWComponent = cubicKernelGradComponent(op.x[d], nbr->x[d], dist, smoothingLength);
                op.accelPressure[d] += mass * (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;
                gradWNorm += gradWComponent * gradWComponent;
                gradW[d] = gradWComponent;
            }

            gradWNorm = sqrt(gradWNorm);

            // calculate viscosity per spatial axis
            double factor = mass * 2 * gradWNorm / (nbr->density * dist);
            for (int d=0; d<dimensions; d++)
            {
                op.laplacianVelocity[d] -= (op.v[d] - nbr->v[d]) * factor;
            }
        }
    }
}
// Move particles
// Semi-implicit Euler scheme
void updateParticles()
{
    for (auto &op : particles)
    {
        if (op.isBoundary) continue;

        for (int d=0; d<dimensions; d++)
        {
            // reset acceleration
            op.a[d] = 0.0;

            // update particle accelerations due to external force (e.g. gravity)
            if (d == dimensions - 1)
                op.a[d] += gravity;

            // update particle accelerations due to 
            // pressure force, viscosity force and external force (e.g. gravity)
            double kinematicViscosity = dynamicViscosity / op.density;
            op.a[d] += (kinematicViscosity * op.laplacianVelocity[d]); // accel. from viscosity
            op.a[d] -= op.accelPressure[d]; // accel. from pressure

            // Update velocities and positions
            // Semi-implicit Euler scheme
            op.v[d] += deltaT * op.a[d];
            op.x[d] += deltaT * op.v[d];

            if (op.x[d] < minPosition[d])
                minPosition[d] = op.x[d];

            if (op.x[d] > maxPosition[d])
                maxPosition[d] = op.x[d];
        }
    }
}

void openParaviewVideoFile()
{
  videoFile.open( outputDir + "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl
            << "<Collection>";
}

void closeParaviewVideoFile()
{
  videoFile << "</Collection>"
            << "</VTKFile>" << endl;
}

void printParaviewSnapshot()
{
    static int counter = -1;
    counter++;
    stringstream filename;
    filename << "result-" << counter <<  ".vtp";
    ofstream out( outputDir + filename.str().c_str() );
    out << "<VTKFile type=\"PolyData\" >" << endl
        << "<PolyData>" << endl
        << " <Piece NumberOfPoints=\"" << numParticles << "\">" << endl
        << "  <Points>" << endl
        << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
        // << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

    for (int p=0; p<numParticles; p++) 
    {
      Particle *op = &particles[p];
      
      int i=0;
      for (int d=0; d<dimensions; d++, i++)
      {
          out << op->x[d] << " ";
      }

      while (i < 3) 
      {
          out << "0 ";
          i++;
      }
    }

    out << "   </DataArray>" << endl
        << "  </Points>" << endl
        << "  <PointData>" << endl
        << "   <DataArray type=\"Float64\" Name=\"Density\" format=\"ascii\">" << endl;

    
    for (int p=0; p<numParticles; p++) 
    {
        out << particles[p].density << " ";
    }

    out << "   </DataArray>" << endl
        << "   <DataArray type=\"Float64\" Name=\"Pressure\" format=\"ascii\">" << endl;

    for (int p=0; p<numParticles; p++)
    {
        out << particles[p].pressure << " ";
    }

    out << "   </DataArray>" << endl
        << "  </PointData>" << endl
        << " </Piece>" << endl
        << "</PolyData>" << endl
        << "</VTKFile>"  << endl;

    videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << endl;
}

int readParameters()
{
    string line;
    ifstream inputFile (inputFilename);
    stringstream linestream;
    if (inputFile.is_open())
    {
        getline(inputFile, line);
        linestream.str(line);
        linestream >> dimensions;
        cout << "Dimensions: " << dimensions << endl;
        
        getline(inputFile, line);
        linestream.str(line);
        linestream >> restDensity;
        cout << "Rest density: " << restDensity << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> dynamicViscosity;
        cout << "Dynamic viscosity: " << dynamicViscosity << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> stiffness;
        cout << "Stiffness coefficient: " << stiffness << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> fpSize;
        cout << "Fluid particle size: " << fpSize << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> bpSize;
        cout << "Boundary particle size: " << bpSize << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> smoothingLength;
        cout << "Smoothing length: " << smoothingLength << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> deltaT;
        cout << "Timestep: " << deltaT << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> deltaTPlot;
        cout << "Plot interval: " << deltaTPlot << endl;

        getline(inputFile, line);
        linestream.str(line);
        linestream >> tFinal; 
        cout << "Final time: " << tFinal << endl;

        getline(inputFile, line);
        caseFilename = line;
        cout << "Case file: " << caseFilename << endl;

        inputFile.close();
        return 1;
    }
    else
    {
        cout << "Unable to open parameters file." << endl;
        return 0;
    }
}

int main(int argc, char** argv)
{
    if (readParameters())
    {
        cout << "All parameters read successfully." << endl;
    }
    else return 1;

    t = 0.0;
    tPlot = t + deltaTPlot;
    neighbourhoodSize = pow(3, dimensions);
    neighbourRadius = 2 * smoothingLength;
    minPosition = new double[dimensions];
    maxPosition = new double[dimensions];
    gridDims = new int[dimensions];
    dimFactors = new int[dimensions];
  
    // Set particle masses
    fpMass = restDensity * pow(fpSize, dimensions);
    bpMass = restDensity * pow(bpSize, dimensions);

    cout << "Fluid particle mass: " << fpMass << endl
         << "Boundary particle mass: " << bpMass << endl;

    initialiseParticles();
    initCoordOffsets();

    prepareGrid(minPosition, maxPosition, neighbourRadius);
    projectParticlesToGrid(minPosition, neighbourRadius);
    calcParticleDensities();

    // Plot initial state
    openParaviewVideoFile();
    printParaviewSnapshot();
    cout << "Time: " << t << endl;

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
            printParaviewSnapshot();
            cout << "Time: " << t << endl;
            tPlot += deltaTPlot;
        }
        t += deltaT;        
    }

    // Plot final state of the system
    printParaviewSnapshot();
    cout << "Time: " << t << endl;

    // Cleanup
    releaseParticles();
    closeParaviewVideoFile();

    // Testing/troubleshooting actions
    /* plotCubicKernel(smoothingLength); */    
}

