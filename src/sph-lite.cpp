#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#define coordOffsetsIndex(row, col, dim) ((dim*row) + col)

using namespace std;

const int DIM = 2;
const double cubicKernelWNormalFactors[3] = {2.0/3.0, 10.0/(7*M_PI), 1/M_PI};

double minPosition[DIM], maxPosition[DIM];
int gridDims[DIM], dimFactors[DIM];
int neighbourhoodSize = pow(3, DIM);

unsigned int numParticles;
double t, tFinal, tPlot, deltaT, deltaTPlot;
double pMass, stiffness, restDensity, dynamicViscosity, gravity = -9.81;
double neighbourRadius, smoothingLength;
int *coordOffsets;
vector<vector<int>> grid;
ofstream videoFile;

class Particle
{
    public:
        double x[DIM], v[DIM], a[DIM];
        double accelPressure[DIM], laplacianVelocity[DIM];
        double pressure, density;
        vector<int> neighbours;
        vector<double> neighbourDist;
        unsigned int gridCellNo;
};

Particle *particles;

// Cubic spline kernel function with compact support of size 2h
double cubicKernelW(double r, double h)
{
   double q = r/h;
   double result = cubicKernelWNormalFactors[DIM-1] / pow(h, (double)DIM);
   
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
    double result = cubicKernelWNormalFactors[DIM-1] / pow(h, (double)DIM);

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

// Plot kernel function and first derivative for testing
void plotCubicKernel(double h, double dq = 0.1)
{
    cout << "h=" << h << endl;
    cout << "q," << "W," << "dW" << endl;

    for (double r=0.0; r<(2*h)+1; r+=dq)
    {
        cout << r/h << "," << cubicKernelW(r, h) << "," << cubicKernelDelW(r, h) << endl;
    }
}

// Return pressure value as a function of density according to a state equation
inline double pressureStateEquation(double rho)
{
    return stiffness * ((rho/restDensity) - 1);
}

void initialiseParticles()
{
    particles = new Particle[numParticles]; // Initialise particles array
    
    // Initialise array to hold relative coordinates of a grid cell neighbourhood
    // to the central cell
    coordOffsets = new int[neighbourhoodSize*DIM];
    for (int i=0; i<neighbourhoodSize; i++)
    {
        for (int d=0; d<DIM; d++)
        {
            int power = pow(3, DIM-d-1);
            coordOffsets[coordOffsetsIndex(i, d, DIM)] = ((i / power) % 3) - 1;
        }
    }
}

void releaseParticles()
{
    delete[] particles;
    delete[] coordOffsets;
}

void prepareGrid(double *minP, double *maxP, double cellLength)
{
    int gridSize = 1;
    int dimFactor = 1;
    for (int d=0; d<DIM; d++)
    {
        int dim = ceil((maxP[d] - minP[d]) / cellLength);
        gridSize *= dim;
        gridDims[d] = dim;
        dimFactors[d] = dimFactor;
        dimFactor *= dim;
    }
    
    grid.resize(gridSize, vector<int>(numParticles / gridSize));
}

void projectParticlesToGrid(double *minP, double cellLength)
{
    for (int i=0; i<numParticles; i++)
    {
        Particle p = particles[i];
        
        int cellNo = 0;
        for (int d=0; d<DIM; d++)
        {
            cellNo += dimFactors[d] * (int)((p.x[d] - minP[d]) / cellLength);
        }

        grid[cellNo].push_back(i);
        p.gridCellNo = cellNo;
    }
}

// Calculate particle densities using SPH interpolation
// and particle pressures using state equation.
void calcParticleDensities()
{
    // Loop through all particles
    for (int p=0; p<numParticles; p++)
    {
        Particle op = particles[p];
        int originCellNo = op.gridCellNo;

        // Clear neighbour lists
        op.neighbours.clear();
        op.neighbourDist.clear();

        // Calculate this particle's grid coordinates
        int origin[DIM];
        for (int d=DIM; d>=0; d--)
        {
            int df = dimFactors[d];
            int c = originCellNo / df;
            origin[d] = c;
            originCellNo -= (c * df);
        } 

        // Loop through all the neighbourhood grid cells
        double densitySum = 0;
        for (int i=0; i<neighbourhoodSize; i++)
        {
            // Calculate the cell index in the semi-flattened grid vector.
            int cellNo = 0;
            bool valid_coords = true;
            for (int d=0; d<DIM; d++)
            {
                int c = origin[d] + coordOffsets[coordOffsetsIndex(i, d, DIM)];
                if (c < 0 || c >= gridDims[d])
                {
                    valid_coords = false;
                    break;
                }
                cellNo += (dimFactors[d] * c); 
            }

            if (!valid_coords)
                continue;

            // Loop through all particles in the grid cell
            //  - Calculate particle distance
            //  - Determine if particle is a neighbour
            //  - If neighbour:
            //      add its contribution to the density estimate
            //      store its index in the particle list
            //      store its distance from the sample particle
            vector<int> cell = grid[cellNo];
            for (auto pIndex: cell)
            {
                Particle nbr = particles[pIndex]; // potential neighbour particle
              
                // Calculate distance to sample/origin particle
                double dist = 0;
                for (int d=0; d<DIM; d++)
                {
                    dist += (op.x[d] - nbr.x[d]) * (op.x[d] - nbr.x[d]);
                }
                dist = sqrt(dist);

                // Check if neighbour
                if (dist <= neighbourRadius)
                {
                    op.neighbours.push_back(pIndex);
                    op.neighbourDist.push_back(dist);
                    densitySum += cubicKernelW(dist, smoothingLength);
                }
            }
        }

        densitySum *= pMass; // Complete density summation for equal-mass particles
        op.density = densitySum;
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
    for (int p=0; p<numParticles; p++)
    {
        Particle op = particles[p];

        // Reset accelPressure and laplacianVelocity arrays
        for (int d=0; d<DIM; d++)
        {
            op.accelPressure[d] = 0.0;
            op.laplacianVelocity[d] = 0.0;
        }

        double opRho = op.density, opP = op.pressure;
        double opPOverRhoSquared = opP / (opRho * opRho);

        int numNeighbours = op.neighbours.size();
        int *nbrPtr = (numNeighbours > 0) ? op.neighbours.data() : nullptr;
        double *nbrDistPtr = (numNeighbours > 0) ? op.neighbourDist.data() : nullptr;
        for (int i=0; i<numNeighbours; i++)
        {
            Particle nbr = particles[nbrPtr[i]];
            double dist = nbrDistPtr[i];

            // calculate pressure force per spatial axis
            // as well as kernel gradient
            double gradW[DIM];
            double gradWNorm = 0;
            double nbrPOverRhoSquared = nbr.pressure / (nbr.density * nbr.density);
            for (int d=0; d<DIM; d++)
            {
                double gradWComponent = cubicKernelDelW(dist, smoothingLength) * (op.x[d] - nbr.x[d]) / (dist * smoothingLength);
                op.accelPressure[d] += (opPOverRhoSquared + nbrPOverRhoSquared) * gradWComponent;
                gradWNorm += gradWComponent * gradWComponent;
                gradW[d] = gradWComponent;
            }

            gradWNorm = sqrt(gradWNorm);

            // calculate viscosity per spatial axis
            double factor = 2 * gradWNorm / (nbr.density * dist);
            for (int d=0; d<DIM; d++)
            {
                op.laplacianVelocity[d] -= (op.v[d] - nbr.v[d]) * factor;
            }
        }

        // update particle accelerations due to external force (e.g. gravity)
        op.a[DIM-1] += gravity;
    }
}

// Move particles
// Semi-implicit Euler scheme
void updateParticles()
{
    int p, d;
    for (int i=0; i<numParticles*DIM; i++)
    {
        p = i / 2;
        d = i % 2;
        Particle op = particles[p];

        // update particle accelerations due to 
        // pressure force, viscosity force and external force (e.g. gravity)
        double kinematicViscosity = dynamicViscosity / op.density;
        op.laplacianVelocity[d] *= pMass;
        op.a[d] += (kinematicViscosity * op.laplacianVelocity[d]); // accel. from viscosity
        op.a[d] -= op.accelPressure[d]; // accel. from pressure

        // Update velocities and positions
        // Semi-implicit Euler scheme
        op.v[d] += deltaT * op.a[d];
        op.x[d] += deltaT * op.v[d];
    }
}

void openParaviewVideoFile()
{
  videoFile.open( "result.pvd" );
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
  ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << endl
      << "<PolyData>" << endl
      << " <Piece NumberOfPoints=\"" << numParticles << "\">" << endl
      << "  <Points>" << endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"" << DIM << "\" format=\"ascii\">";
      // << "   <DataArray type=\"Float32\" NumberOfComponents=\"" << DIM << "\" format=\"ascii\">";

  for (int p=0; p<numParticles; p++) 
  {
      Particle op = particles[p];
      
      for (int d=0; d<DIM; d++)
      {
          out << op.x[d] << " ";
      }
  }

  out << "   </DataArray>" << endl
      << "  </Points>" << endl
      << " </Piece>" << endl
      << "</PolyData>" << endl
      << "</VTKFile>"  << endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << endl;
}

int main(int argc, char** argv)
{
    t = 0.0;
    tFinal = 100.0;
    deltaT = 0.01;
    deltaTPlot = 0.1;
    tPlot = t + deltaTPlot;

    initialiseParticles();
    prepareGrid(minPosition, maxPosition, neighbourRadius);

    // Plot initial state
    openParaviewVideoFile();
    printParaviewSnapshot();

    while (t <= tFinal)
    {
        projectParticlesToGrid(minPosition, neighbourRadius);
        calcParticleDensities();
        calcParticleForces();
        updateParticles();

        t += deltaT;        
        if (t >= tPlot)
        {
            // Plot state of the system
            printParaviewSnapshot();
            tPlot = t + deltaTPlot;
        }
    }

    // Plot final state of the system
    printParaviewSnapshot();

    // Cleanup
    releaseParticles();
    closeParaviewVideoFile();

    // Testing/troubleshooting actions
    // plotCubicKernel(1.0);    
}

