#include "eos.h"
#include "particle.h"

#include <cmath>
#include <iostream>

std::ostream& SphEOS::operator<<(std::ostream &out, const SphEOS::EquationOfState &eos)
{
    out << eos.description;
    return out;
}

SphEOS::WeaklyCompressibleEOS::WeaklyCompressibleEOS(const double k, const double rho0) : k(k), rho0(rho0)
{
    description = "Weakly Compressible";
}

double SphEOS::WeaklyCompressibleEOS::getPressure(SphSchemes::Particle& p) const
{
    return k * (pow((p.density/rho0), 7) - 1);
}

SphEOS::IdealGasEOS::IdealGasEOS(const double gamma) : gamma(gamma)
{
    description = "Ideal Gas";
}

double SphEOS::IdealGasEOS::getPressure(SphSchemes::Particle& p) const
{
    return p.density * p.energyPerMass * (gamma - 1);
}
