// SPH Equations of State
#ifndef EOS_H
#define EOS_H

#include "particle.h"

class EquationOfState
{
    public:
        virtual double getPressure(Particle& p) = 0;
        virtual ~EquationOfState() {};
};

class WeaklyCompressibleEOS: public EquationOfState
{
    public:
        WeaklyCompressibleEOS(double k, double rho0) : k(k), rho0(rho0) {};
        double getPressure(Particle& p);

    private:
        double k, rho0;
};

class IdealGasEOS: public EquationOfState
{
    public:
        IdealGasEOS(double gamma) : gamma(gamma) {};
        double getPressure(Particle& p);

    private:
        double gamma;
};

#endif
