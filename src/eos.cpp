#include "eos.h"
#include "particle.h"

#include <cmath>

double WeaklyCompressibleEOS::getPressure(Particle& p)
{
    return k * (pow((p.density/rho0), 7) - 1);
}

double IdealGasEOS::getPressure(Particle& p)
{
    return p.density * p.temperature * (gamma - 1);
}
