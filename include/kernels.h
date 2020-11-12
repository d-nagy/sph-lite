// SPH Kernel classes
#ifndef KERNELS_H
#define KERNELS_H

#include <cmath>
#include <string>

class SphKernel
{
    public:
        virtual double W(double r, double h) = 0;
        virtual double gradWComponent(double xi, double xj, double r, double h) = 0;
        virtual inline double getSupportRadius(double h) = 0;
        virtual ~SphKernel() {};
        friend std::ostream& operator<<(std::ostream &out, const SphKernel &kernel);

    protected:
        std::string description;
        virtual double delW(double r, double h) = 0;
        double powerH(double h, int d);
};

class CubicSplineKernel: public SphKernel
{
    public:
        CubicSplineKernel(int d);
        double W(double r, double h);
        double gradWComponent(double xi, double xj, double r, double h);
        double getSupportRadius(double h);
        void plotKernel(double h);

    private:
        int dims;
        const double normalFactors[3] = {2.0/3.0, 10.0/(7*M_PI), 1/M_PI};
        double delW(double r, double h);
};

#endif
