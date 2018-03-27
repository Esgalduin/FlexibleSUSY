#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "pv2.hpp"

double sqr(double x) { return x*x; }
double abssqrt(double x) { return std::sqrt(std::abs(x)); }

int main(int argc, char *argv[])
{
   if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <file>\n";
      return 0;
   }

   std::ifstream ifs(argv[1]);
   std::string line;

   std::cout.setf(std::ios_base::fixed, std::ios_base::floatfield);
   std::cout.precision(16);
   std::cout << "# A0(m0^2) A0(m1^2) B0(p^2,m0^2,m1^2) B1(p^2,m0^2,m1^2)"
      " B00(p^2,m0^2,m1^2)\n";

   while (std::getline(ifs, line)) {
      if (line.empty() || line[0] == '#')
         continue;

      std::istringstream iss(line);
      double p2 = 0., m02 = 0., m12 = 0., mu2 = 0.;
      iss >> p2 >> m02 >> m12 >> mu2;

      std::cout
         << flexiblesusy::a0(m02, mu2) << ' '
         << flexiblesusy::a0(m12, mu2) << ' '
         << flexiblesusy::b0(p2, m02, m12, mu2) << ' '
         << -flexiblesusy::b1(p2, m02, m12, mu2) << ' '
         << flexiblesusy::b22(p2, m02, m12, mu2) << ' '
         << flexiblesusy::c0(p2, m02, m12) << '\n';
   }

   return 0;
}
