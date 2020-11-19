#include "eos.h"
#include "particle.h"

#include <cmath>
#include <iostream>

std::ostream& operator<<(std::ostream &out, const EquationOfState &eos)
{
    out << eos.description;
    return out;
}

WeaklyCompressibleEOS::WeaklyCompressibleEOS(double k, double rho0) : k(k), rho0(rho0)
{
    description = "Weakly Compressible";
}

double WeaklyCompressibleEOS::getPressure(Particle& p)
{
    return k * (pow((p.density/rho0), 7) - 1);
}

IdealGasEOS::IdealGasEOS(double gamma) : gamma(gamma)
{
    description = "Ideal Gas";
}

double IdealGasEOS::getPressure(Particle& p)
{
    return p.density * p.energyPerMass * (gamma - 1);
}
