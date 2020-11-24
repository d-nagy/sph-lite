// SPH Kernel classes
#ifndef KERNELS_H
#define KERNELS_H

#include <cmath>
#include <string>

namespace SphKernels
{
    class SphKernel
    {
        public:
            virtual double W(const double r, const double h) = 0;
            virtual double gradWComponent(const double xi, const double xj,
                                          const double r, const double h) = 0;
            virtual inline double getSupportRadius(const double h) = 0;
            virtual ~SphKernel() {};
            friend std::ostream& operator<<(std::ostream &out, const SphKernels::SphKernel &kernel);

        protected:
            std::string description;
            virtual double delW(const double r, const double h) = 0;
            double powerH(const double h, const int d);
    };

    class CubicSplineKernel: public SphKernel
    {
        public:
            CubicSplineKernel(const int d);
            double W(const double r, const double h);
            double gradWComponent(const double xi, const double xj,
                                  const double r, const double h);
            double getSupportRadius(const double h);
            void plotKernel(const double h);

        private:
            int dims;
            const double normalFactors[3] = {2.0/3.0, 10.0/(7*M_PI), 1/M_PI};
            double delW(const double r, const double h);
    };
}

#endif
