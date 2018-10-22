// ====================================================================
// This file is part of FlexibleSUSY.
//
// FlexibleSUSY is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published
// by the Free Software Foundation, either version 3 of the License,
// or (at your option) any later version.
//
// FlexibleSUSY is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with FlexibleSUSY.  If not, see
// <http://www.gnu.org/licenses/>.
// ====================================================================

#include "pv2.hpp"
#include "numerics2.hpp"

#include <algorithm>
#include <complex>
#include <cmath>
#include <limits>

namespace flexiblesusy {

namespace {

constexpr double EPSTOL = 1.0e-11; ///< underflow accuracy
constexpr double TOL = 1e-4;

constexpr double sqr(double a) noexcept { return a*a; }
constexpr double pow3(double a) noexcept { return a*a*a; }
constexpr double pow6(double a) noexcept { return a*a*a*a*a*a; }
constexpr double log_abs(double a) noexcept { return std::log(std::abs(a)); }

bool is_close(double m1, double m2, double tol) noexcept
{
   using std::abs;
   const double mmax = abs(std::max(abs(m1), abs(m2)));
   const double max_tol = tol * mmax;

   if (max_tol == 0. && mmax != 0. && tol != 0.)
      return abs(m1 - m2) <= tol;

   return abs(m1 - m2) <= max_tol;
}

/// returns a/b if a/b is finite, otherwise returns numeric_limits::max()
template <typename T>
T divide_finite(T a, T b) noexcept {
   T result = a / b;
   if (!std::isfinite(result))
      result = std::numeric_limits<T>::max();
   return result;
}

double fB(const std::complex<double>& a) noexcept
{
   using flexiblesusy::fast_log;
   using std::abs;
   const double x = a.real();

   if (abs(x) < EPSTOL)
      return -1. - x + sqr(x) * 0.5;

   if (is_close(x, 1., EPSTOL))
      return -1.;

   return std::real(fast_log(1. - a) - 1. - a * fast_log(1.0 - 1.0 / a));
}

} // anonymous namespace

double a0(double m2, double q2) noexcept
{
   if (std::abs(m2) < 1e-8)
      return 0.;

   return m2 * (1.0 - log_abs(m2 / q2));
}

/**
 * B0 function with squared arguments, from hep-ph/9606211.
 * Note it returns the REAL PART ONLY.
 */
double b0(double p2, double m12, double m22, double q2) noexcept
{
   using std::abs;

   // protect against infrared divergence
   if (is_close(p2, 0., EPSTOL) && is_close(m12, 0., EPSTOL)
       && is_close(m22, 0., EPSTOL))
      return 0.;

   double ans = 0.;
   const double mMax2 = std::max(abs(m12), abs(m22));
   const double pTest = abs(divide_finite(p2, mMax2));

   if (pTest > 1e-10) {
      const double s = p2 - m22 + m12;
      const std::complex<double> iEpsilon(0., EPSTOL * mMax2);
      const std::complex<double> rt = sqrt(sqr(s) - 4.*p2*(m12 - iEpsilon));
      const std::complex<double> xPlus = 0.5 * (s + rt) / p2;
      const std::complex<double> xMinus = 0.5 * (s - rt) / p2;

      ans = -log_abs(p2 / q2) - fB(xPlus) - fB(xMinus);
   } else {
      if (is_close(m12, m22, EPSTOL)) {
         ans = -log_abs(m12 / q2);
      } else {
         const double mMin2 = std::min(abs(m12), abs(m22));

         if (mMin2 < 1.e-30) {
            ans = 1. - log_abs(mMax2 / q2);
         } else {
            ans = (m12*(1. - log_abs(m12/q2)) - m22*(1. - log_abs(m22/q2)))
               / (m12 - m22);
         }
      }
   }

   return ans;
}

/**
 * Passarino-Veltman B1 function with squared mass arguments.
 * Not symmetric w.r.t. m1 and m2
 */
double b1(double p2, double m12, double m22, double q2) noexcept
{
   using std::abs;

   // protect against infrared divergence
   if (is_close(p2, 0., EPSTOL) && is_close(m12, 0., EPSTOL)
       && is_close(m22, 0., EPSTOL))
      return 0.;

   double ans = 0.;
   const double mMax2 = std::max(abs(m12), abs(m22));
   const double pTest = abs(divide_finite(p2, mMax2));

   /// Decides level at which one switches to p=0 limit of calculations
   const double pTolerance = 1e-4;

   if (pTest > pTolerance) {
      ans = (a0(m22,q2) - a0(m12,q2) + (p2 + m12 - m22)
             * b0(p2,m12,m22,q2)) / (2.*p2);
   } else if (abs(m12) > 1.0e-30 && abs(m22) > 1.0e-30) {
      const double m14 = sqr(m12), m24 = sqr(m22);
      const double m16 = m12*m14 , m26 = m22*m24;
      const double m18 = sqr(m14), m28 = sqr(m24);
      const double p4 = sqr(p2);
      if (abs(m12 - m22) < pTolerance * mMax2) {
         ans = 0.08333333333333333*p2/m22
            + 0.008333333333333333*p4/m24
            + sqr(m12 - m22)*(0.041666666666666664/m24 +
                              0.016666666666666666*p2/m26 +
                              0.005357142857142856*p4/m28)
            + (m12 - m22)*(-0.16666666666666666/m22 -
                           0.03333333333333333*p2/m24 -
                           0.007142857142857142*p4/m26)
            - 0.5*log_abs(m22/q2);
      } else {
         const double l12 = log_abs(m12/m22);

         ans = (3*m14 - 4*m12*m22 + m24 - 2*m14*l12)/(4.*sqr(m12 - m22))
            + (p2*(4*pow3(m12 - m22)*
                   (2*m14 + 5*m12*m22 - m24) +
                   (3*m18 + 44*m16*m22 - 36*m14*m24 - 12*m12*m26 + m28)*p2
                   - 12*m14*m22*(2*sqr(m12 - m22) + (2*m12 + 3*m22)*p2)*l12))/
            (24.*pow6(m12 - m22)) - 0.5*log_abs(m22/q2);
      }
   } else {
      if (abs(m12) > abs(m22))
         ans = -0.5*log_abs(m12/q2) + 0.75;
      else
         ans = -0.5*log_abs(m22/q2) + 0.25;
   }

   return ans;
}

double b22(double p2,  double m12, double m22, double q2) noexcept
{
   using std::abs;

   // protect against infrared divergence
   if (is_close(p2, 0., EPSTOL) && is_close(m12, 0., EPSTOL)
       && is_close(m22, 0., EPSTOL))
      return 0.;

   double ans = 0.;
   const double mMax2 = std::max(abs(m12), abs(m22));
   const double pTolerance = 1.0e-10;

   if (abs(p2) < pTolerance * mMax2) {
      if (is_close(m12, m22, EPSTOL)) {
         ans = -m12 * log_abs(m12 / q2) * 0.5 + m12 * 0.5;
      } else {
         if (abs(m12) > EPSTOL && abs(m22) > EPSTOL) {
            ans = 0.375 * (m12 + m22) - 0.25 *
               (sqr(m22) * log_abs(m22 / q2) - sqr(m12) *
                log_abs(m12 / q2)) / (m22 - m12);
         } else {
            if (abs(m12) < EPSTOL) {
               ans = 0.375 * m22 - 0.25 * m22 * log_abs(m22 / q2);
            } else {
               ans = 0.375 * m12 - 0.25 * m12 * log_abs(m12 / q2);
            }
         }
      }
   } else {
      const double b0Save = b0(p2, m12, m22, q2);
      const double a01 = a0(m12, q2);
      const double a02 = a0(m22, q2);

      ans = 1.0 / 6.0 *
         (0.5 * (a01 + a02) + (m12 + m22 - 0.5 * p2)
          * b0Save + (m22 - m12) / (2.0 * p2) *
          (a02 - a01 - (m22 - m12) * b0Save) +
          m12 + m22 - p2 / 3.0);
   }

   return ans;
}

double b22bar(double p2, double m12, double m22, double q2) noexcept
{
   return b22(p2,m12,m22,q2) - 0.25*a0(m12,q2) - 0.25*a0(m22,q2);
}

double f0(double p2, double m12, double m22, double q2) noexcept
{
   return a0(m12,q2) - 2.*a0(m22,q2)
      - (2.0 * p2 + 2.0 * m12 - m22) * b0(p2,m12,m22,q2);
}

double g0(double p2, double m12, double m22, double q2) noexcept
{
   return (p2 - m12 - m22) * b0(p2,m12,m22,q2) - a0(m12,q2) - a0(m22,q2);
}

double h0(double p2, double m12, double m22, double q2) noexcept
{
   return 4. * b22(p2,m12,m22,q2) + g0(p2,m12,m22,q2);
}

double c0(double m12, double m22, double m32) noexcept
{
   double ans = 0.;

   if (is_close(m12,0.,EPSTOL) && is_close(m22,0.,EPSTOL) && is_close(m32,0.,EPSTOL)) {
      // c0 is undefined for m1 == m2 == m3 == 0
      ans = 0.;
   } else if (is_close(m22,0.,EPSTOL) && is_close(m32,0.,EPSTOL)) {
      // c0 is undefined for m2 == m3 == 0
      ans = 0.;
   } else if (is_close(m12,0.,EPSTOL) && is_close(m32,0.,EPSTOL)) {
      // c0 is undefined for m1 == m3 == 0
      ans = 0.;
   } else if (is_close(m12,0.,EPSTOL) && is_close(m22,0.,EPSTOL)) {
      // c0 is undefined for m1 == m2 == 0
      ans = 0.;
   } else if (is_close(m12,0.,EPSTOL)) {
      if (is_close(m22,m32,EPSTOL)) {
         ans = -1./m22;
      } else {
         ans = log_abs(m32/m22)/(m22 - m32);
      }
   } else if (is_close(m22,0.,EPSTOL)) {
      if (is_close(m12,m32,EPSTOL)) {
         ans = -1./m12;
      } else {
         ans = log_abs(m32/m12)/(m12 - m32);
      }
   } else if (is_close(m32,0.,EPSTOL)) {
      if (is_close(m12,m22,EPSTOL)) {
         ans = -1./m12;
      } else {
         ans = log_abs(m22/m12)/(m12 - m22);
      }
   } else if (is_close(m22, m32, EPSTOL)) {
      if (is_close(m12, m22, EPSTOL)) {
         ans = ( - 0.5 / m22 );
      } else {
         ans = ( m12 / sqr(m12-m22) * log_abs(m22/m12)
                 + 1.0 / (m12 - m22) );
      }
   } else if (is_close(m12, m22, EPSTOL)) {
      ans = ( - (1.0 + m32 / (m22-m32) * log_abs(m32/m22))
              / (m22-m32) );
   } else if (is_close(m12, m32, EPSTOL)) {
      ans = ( - (1.0 + m22 / (m32-m22) * log_abs(m22/m32))
              / (m32-m22) );
   } else {
      ans = (1.0 / (m22 - m32) *
             (m22 / (m12 - m22) *
              log_abs(m22 / m12) -
              m32 / (m12 - m32) *
              log_abs(m32 / m12)) );
   }

   return ans;
}

double d0(double m12, double m22, double m32, double m42) noexcept
{
   if (is_close(m12, m22, EPSTOL)) {
      if (is_close(m22,0.,EPSTOL)) {
         // d0 is undefined for m1 == m2 == 0
         return 0.;
      } else if (is_close(m32,0.,EPSTOL)) {
         return (-m22 + m42 - m22 * log_abs(m42/m22))/
            (m22 * sqr(m22 - m42));
      } else if (is_close(m42,0.,EPSTOL)) {
         return (-m22 + m32 - m22 * log_abs(m32/m22))/
            (m22 * sqr(m22 - m32));
      } else if (is_close(m22, m32, EPSTOL) && is_close(m22, m42, EPSTOL)) {
         return 1.0 / (6.0 * sqr(m22));
      } else if (is_close(m22, m32, EPSTOL)) {
         return (sqr(m22) - sqr(m42) + 2.0 * m42 * m22 * log_abs(m42 / m22)) /
            (2.0 * m22 * sqr(m22 - m42) * (m22 - m42));
      } else if (is_close(m22, m42, EPSTOL)) {
         return (sqr(m22) - sqr(m32) + 2.0 * m32 * m22 * log_abs(m32 / m22)) /
            (2.0 * m22 * sqr(m22 - m32) * (m22 - m32));
      } else if (is_close(m32, m42, EPSTOL)) {
         return -1.0 / sqr(m22 - m32) *
            ((m22 + m32) / (m22 - m32) * log_abs(m32 / m22) + 2.0);
      }

      return
         (m42 / sqr(m22 - m42) * log_abs(m42 / m22) +
          m42 / (m22 * (m22 - m42)) -
          m32 / sqr(m22 - m32) * log_abs(m32 / m22) -
          m32 / (m22 * (m22 - m32))) / (m32 - m42);
   }
   return (c0(m12, m32, m42) - c0(m22, m32, m42)) / (m12 - m22);
}

double d27(double m12, double m22, double m32, double m42) noexcept
{
   if (is_close(m12, m22, EPSTOL))
      m12 += TOL * 0.2;

   return (m12 * c0(m12, m32, m42) - m22 * c0(m22, m32, m42))
      / (4.0 * (m12 - m22));
}

/**
 * Derivative of B0(p^2, m1^2, m2^2, Q^2) w.r.t. p^2.
 *
 * @note Implemented only in the p^2 = 0 limit.
 *
 * @param p2 squared momentum
 * @param m2a squared mass
 * @param m2b squared mass
 *
 * @return derivative of B0 w.r.t. p^2
 */
double d1_b0(double p2, double m2a, double m2b) noexcept
{
   using std::abs;

   const double m4a = m2a * m2a;
   const double m4b = m2b * m2b;

   if ((std::abs(m2a) < 0.0001) != (std::abs(m2b) < 0.0001)) {
      return (m4a - m4b) / (2. * pow3(m2a - m2b));
   } else if (std::abs(m2a) < 0.0001 && std::abs(m2b) < 0.0001) {
      return 0.;
   } else if (std::abs(m2b - m2a) < 0.001) {
      return 1./(6. * m2a) + (m2a - m2b)/(12.* m4a);
   }

   return (m4a - m4b + 2. * m2a * m2b * std::log(m2b/m2a))
      /(2. * pow3(m2a - m2b));
}

} // namespace flexiblesusy
