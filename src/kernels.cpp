#include "kernels.h"
#include <iostream>

std::ostream& operator<<(std::ostream &out, const SphKernel &kernel)
{
    out << kernel.description;
    return out;
}

double SphKernel::powerH(double h, int d)
{
    double result = h;
    switch (d)
    {
        case 2:
            result *= h;
            break;
        case 3:
            result *= h*h;
        default:
            break;
    }
    return result;
}

CubicSplineKernel::CubicSplineKernel(int d) : dims(d)
{
    description = "Cubic Spline";
}

double CubicSplineKernel::W(double r, double h)
{
    double q = r/h;
    double result = normalFactors[dims-1] / powerH(h, dims);

    if (q >= 2)
    {
        result *= 0;
    }
    else if (q >= 1)
    {
        result *= (0.25 * (2 - q)*(2 - q)*(2 - q));
    }
    else
    {
        result *= (1 - 1.5*q*q*(1 - q/2));
    }

    return result;
}

double CubicSplineKernel::gradWComponent(double xi, double xj, double r, double h)
{
    return delW(r, h) * (xi - xj) / (r * h);
}

inline double CubicSplineKernel::getSupportRadius(double h)
{
    return 2 * h;
}

double CubicSplineKernel::delW(double r, double h)
{
    double q = r/h;
    double result = normalFactors[dims-1] / powerH(h, dims);

    if (q >= 2)
    {
        result *= 0;
    }
    else if (q >= 1)
    {
        result *= (-0.75 * ((2 - q)*(2 - q)));
    }
    else
    {
        // result *= (3 * (((2 - q)*(2 - q)/-4.0) + ((1 - q)*(1 - q))));
        result *= (0.75*q * (3*q - 4));
    }

    return result;
}

void CubicSplineKernel::plotKernel(double h)
{
    std::cout << "h=" << h << std::endl;
    std::cout << "q," << "W," << "dW" << std::endl;

    double maxR = 2*h;
    double dq = 2*h / 100;
    maxR += dq;
    for (double r=0.0; r<maxR; r+=dq)
    {
        std::cout << r/h << "," << W(r, h) << "," << delW(r, h) << std::endl;
    }
}
