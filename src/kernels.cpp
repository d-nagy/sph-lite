#include "kernels.h"
#include <iostream>

std::ostream& SphKernels::operator<<(std::ostream &out, const SphKernels::SphKernel &kernel)
{
    out << kernel.description;
    return out;
}

double SphKernels::SphKernel::powerH(const double h, const int d)
{
    double result = h;
    switch (d)
    {
        case 2:
            return result * h;
            break;
        case 3:
            return result * h*h;
        default:
            return result;
    }
}

SphKernels::CubicSplineKernel::CubicSplineKernel(const int d) : dims(d)
{
    description = "Cubic Spline";
}

double SphKernels::CubicSplineKernel::W(const double r, const double h)
{
    const double q = r/h;
    double result = normalFactors[dims-1] / powerH(h, dims);

    if (q >= 2)
    {
        return result * 0;
    }
    else if (q >= 1)
    {
        return result * (0.25 * (2 - q)*(2 - q)*(2 - q));
    }
    else
    {
        return result * (1 - 1.5*q*q*(1 - q/2));
    }
}

double SphKernels::CubicSplineKernel::gradWComponent(const double xi, const double xj,
                                                     const double r, const double h)
{
    return (r > 1e-12) ? delW(r, h) * (xi - xj) / (r * h) : 0;
}

inline double SphKernels::CubicSplineKernel::getSupportRadius(const double h)
{
    return 2 * h;
}

double SphKernels::CubicSplineKernel::delW(const double r, const double h)
{
    if (r < 1e-12) { return 0; }

    const double q = r/h;
    double result = normalFactors[dims-1] / powerH(h, dims);

    if (q >= 2)
    {
        return result * 0;
    }
    else if (q >= 1)
    {
        return result * (-0.75 * ((2 - q)*(2 - q)));
    }
    else
    {
        return result * (0.75*q * (3*q - 4));
    }
}

void SphKernels::CubicSplineKernel::plotKernel(const double h)
{
    std::cout << "h=" << h << '\n';
    std::cout << "q," << "W," << "dW" << '\n';

    double maxR = 2*h;
    const double dq = 2*h / 100;
    maxR += dq;
    for (double r=0.0; r<maxR; r+=dq)
    {
        std::cout << r/h << "," << W(r, h) << "," << delW(r, h) << '\n';
    }
}
