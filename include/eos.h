// SPH Equations of State
#ifndef EOS_H
#define EOS_H

#include "particle.h"

#include <iostream>
#include <string>

namespace SphEOS
{
    class EquationOfState
    {
        public:
            virtual double getPressure(SphSchemes::Particle& p) const = 0;
            virtual ~EquationOfState() {};
            friend std::ostream& operator<<(std::ostream &out, const EquationOfState &eos);

        protected:
            std::string description;
    };

    class WeaklyCompressibleEOS: public EquationOfState
    {
        public:
            WeaklyCompressibleEOS(const double k, const double rho0);
            double getPressure(SphSchemes::Particle& p) const;

        private:
            double k, rho0;
    };

    class IdealGasEOS: public EquationOfState
    {
        public:
            IdealGasEOS(const double gamma);
            double getPressure(SphSchemes::Particle& p) const;

        private:
            double gamma;
    };
}

#endif
