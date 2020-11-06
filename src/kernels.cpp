#include "kernels.h"
#include <cmath>
#include <iostream>

CubicSplineKernel::CubicSplineKernel(int d) : dims(d) {}

double CubicSplineKernel::W(double r, double h)
{
    double q = r/h;
    double result = normalFactors[dims-1] / pow(h, (double)dims);

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

double CubicSplineKernel::gradWComponent(double xi, double xj, double r, double h)
{
    return delW(r, h) * (xi - xj) / (r * h);
}

double CubicSplineKernel::delW(double r, double h)
{
    double q = r/h;
    double result = normalFactors[dims-1] / pow(h, (double)dims);

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
