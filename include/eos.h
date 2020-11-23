// SPH Equations of State
#ifndef EOS_H
#define EOS_H

#include "particle.h"

#include <iostream>
#include <string>

class EquationOfState
{
    public:
        virtual double getPressure(Particle& p) = 0;
        virtual ~EquationOfState() {};
        friend std::ostream& operator<<(std::ostream &out, const EquationOfState &eos);

    protected:
        std::string description;
};

class WeaklyCompressibleEOS: public EquationOfState
{
    public:
        WeaklyCompressibleEOS(const double k, const double rho0);
        double getPressure(Particle& p);

    private:
        double k, rho0;
};

class IdealGasEOS: public EquationOfState
{
    public:
        IdealGasEOS(const double gamma);
        double getPressure(Particle& p);

    private:
        double gamma;
};

#endif
