#include "Colour.h"

namespace Colour
{
    void rgb2xyz( double r, double g, double b, double& x, double& y, double& z )
    {
        r = (r <= 0.04045) ? r/12.92 : std::pow((r + 0.055) / 1.055, 2.4);
        g = (g <= 0.04045) ? g/12.92 : std::pow((g + 0.055) / 1.055, 2.4);
        b = (b <= 0.04045) ? b/12.92 : std::pow((b + 0.055) / 1.055, 2.4);

        x = 0.412453*r + 0.357580*g + 0.180423*b;
        y = 0.212671*r + 0.715160*g + 0.072169*b;
        z = 0.019334*r + 0.119193*g + 0.950227*b;
    }

    void xyz2rgb( double x, double y, double z, double& r, double& g, double& b )
    {
        r = 3.240479*x - 1.537150*y - 0.498535*z;
        g = -0.969256*x + 1.875992*y + 0.041556*z;
        b = 0.055648*x - 0.204043*y  + 1.057311*z;

        r = (r <= 0.0031308) ? r*12.92 : 1.055*std::pow(r, 1.0/2.4) -0.055;
        g = (g <= 0.0031308) ? g*12.92 : 1.055*std::pow(g, 1.0/2.4) -0.055;
        b = (b <= 0.0031308) ? b*12.92 : 1.055*std::pow(b, 1.0/2.4) -0.055;
    }

    void xyz2lab( double x, double y, double z, double& l, double& a, double& b )
    {
        const static std::vector<double> D65 {0.95047, 1, 1.08883};
        const static double epsilon = std::pow(6.0/29.0, 3);
        const static double kappa = (1.0/116.0)*std::pow(29.0/3.0, 3);

        x /= D65[0];
        y /= D65[1];
        z /= D65[2];



        x = (x<=epsilon) ? kappa*x + 16.0 / 116.0 : std::pow(x, 0.3333333333333);
        y = (y<=epsilon) ? kappa*y + 16.0 / 116.0 : std::pow(y, 0.3333333333333);
        z = (z<=epsilon) ? kappa*z + 16.0 / 116.0 : std::pow(z, 0.3333333333333);

        l = 116.0 * y - 16;
        a = 500 * (x - y);
        b = 200 * (y - z);
    }

    void lab2xyz( double l, double a, double b , double& x, double& y, double& z)
    {
        const static std::vector<double> D65 {0.95047, 1, 1.08883};
        const static double epsilon = std::pow(6.0/29.0, 3);
        const static double kappa = (1.0/116.0)*std::pow(29.0/3.0, 3);

        double Lp = (l + 16.0) / 116.0;

        x = Lp + a/500;
        y = Lp;
        z = Lp - b/200;

        x = (x*x*x > epsilon) ? x*x*x : (x - 16.0/116.0)/kappa;
        y = (y*y*y > epsilon) ? y*y*y : (y - 16.0/116.0)/kappa;
        z = (z*z*z > epsilon) ? z*z*z : (z - 16.0/116.0)/kappa;

        x = D65[0] * x;
        y = D65[1] * y;
        z = D65[2] * z;
    }

    void rgb2lab(double r, double g, double b, double& l, double& a, double& bp )
    {
        rgb2xyz(r,g,b,l,a,bp);
        xyz2lab(l,a,bp,l,a,bp);
    }

    void lab2rgb(double l, double a, double bp, double& r, double& g, double& b )
    {
        lab2xyz(l,a,bp, r,g,b);
        xyz2rgb(r,g,b,r,g,b);
    }
}
