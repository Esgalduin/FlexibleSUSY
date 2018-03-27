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

#include <cmath>
#include <limits>
#include <random>

#include "numerics.h"
#include "numerics2.hpp"
#include "pv2.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_pv2

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

struct Values {
   Values() = default;
   Values(double p2_, double m12_, double m22_, double q2_)
      : p2(p2_), m12(m12_), m22(m22_), q2(q2_) {}
   double p2{}, m12{}, m22{}, q2{};
};

const std::vector<Values> positive_vals = {
   Values(0.  , 1.  , 0., 1.),
   Values(1e-1, 1.  , 0., 1.),
   Values(1e-2, 1.  , 0., 1.),
   Values(1e-3, 1.  , 0., 1.),
   Values(1e-4, 1.  , 0., 1.),
   Values(1e-5, 1.  , 0., 1.),
   Values(0.  , 1.  , 1e-15, 1.),
   Values(1e-1, 1.  , 1e-15, 1.),
   Values(1e-2, 1.  , 1e-15, 1.),
   Values(1e-3, 1.  , 1e-15, 1.),
   Values(1e-4, 1.  , 1e-15, 1.),
   Values(1e-5, 1.  , 1e-15, 1.),
   Values(0.  , 1e20, 0., 1.),
   Values(1e-1, 1e20, 0., 1.),
   Values(1e-2, 1e20, 0., 1.),
   Values(1e-3, 1e20, 0., 1.),
   Values(1e-4, 1e20, 0., 1.),
   Values(1e-5, 1e20, 0., 1.),
   Values(0.  , 0.  , 1., 1.),
   Values(1e-1, 0.  , 1., 1.),
   Values(1e-2, 0.  , 1., 1.),
   Values(1e-3, 0.  , 1., 1.),
   Values(1e-4, 0.  , 1., 1.),
   Values(1e-5, 0.  , 1., 1.),
   Values(0.  , 1e20, 1., 1.),
   Values(1e-1, 1e20, 1., 1.),
   Values(1e-2, 1e20, 1., 1.),
   Values(1e-3, 1e20, 1., 1.),
   Values(1e-4, 1e20, 1., 1.),
   Values(1e-5, 1e20, 1., 1.),
   Values(0.  , 0., 1e20, 1.),
   Values(1e-1, 0., 1e20, 1.),
   Values(1e-2, 0., 1e20, 1.),
   Values(1e-3, 0., 1e20, 1.),
   Values(1e-4, 0., 1e20, 1.),
   Values(1e-5, 0., 1e20, 1.),
   Values(1.  , 1.  , 1., 1.)
};

std::vector<Values> generate_random_values(int n, double start, double stop)
{
   std::minstd_rand gen;
   std::uniform_real_distribution<double> dist(start, stop);

   const auto rand = [&dist,&gen](){
      return Values(dist(gen), dist(gen), dist(gen), std::abs(dist(gen)));
   };

   std::vector<Values> v(n);
   std::generate(begin(v), end(v), rand);

   return v;
}

template <class T>
std::vector<T> concat(const std::vector<T>& v1, const std::vector<T>& v2)
{
   std::vector<T> res = v1;
   res.reserve(v1.size() + v2.size());

   for (const auto v: v2)
      res.push_back(v);

   return res;
}

constexpr double sqr(double x) noexcept { return x*x; }
constexpr double sqrtabs(double x) noexcept { return std::sqrt(std::abs(x)); }

const double scale  = 100;
const double scale2 = sqr(scale);
const double p  = 91.0;
const double p2 = sqr(p);

BOOST_AUTO_TEST_CASE( test_ReA0 )
{
   BOOST_CHECK_EQUAL(flexiblesusy::a0(0., scale2), 0.);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(p, scale), flexiblesusy::a0(p2, scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(0, scale), flexiblesusy::a0(0., scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(1e-4, scale), flexiblesusy::a0(1e-8, scale2), 1e-15);
   BOOST_CHECK_CLOSE_FRACTION(softsusy::a0(1e-5, scale), flexiblesusy::a0(1e-10, scale2), 1e-15);
}

BOOST_AUTO_TEST_CASE( test_softsusy )
{
   const auto rand_vals = generate_random_values(10000, 0., 2000);
   const auto vals = concat(positive_vals, rand_vals);

   for (const auto v: vals) {
      const auto p2 = v.p2;
      const auto m12 = v.m12;
      const auto m22 = v.m22;
      const auto q2 = v.q2;
      const auto p = sqrtabs(p2);
      const auto m1 = sqrtabs(m12);
      const auto m2 = sqrtabs(m22);
      const auto q = sqrtabs(q2);

      BOOST_TEST_MESSAGE("testing p2 = " << p2 << ", m12 = " << m12
                         << ", m22 = " << m22 << ", q2 = " << q2);

      {
         const auto v1 = softsusy::a0(m1,q);
         const auto v2 = flexiblesusy::a0(m12,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 8e-5);
      }

      {
         const auto v1 = softsusy::b0(p,m1,m2,q);
         const auto v2 = flexiblesusy::b0(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 8e-5);
      }

      {
         const auto v1 = softsusy::b1(p,m1,m2,q);
         const auto v2 = flexiblesusy::b1(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-8);
      }

      {
         const auto v1 = softsusy::b22(p,m1,m2,q);
         const auto v2 = flexiblesusy::b22(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-7);
      }

      {
         const auto v1 = softsusy::b22bar(p,m1,m2,q);
         const auto v2 = flexiblesusy::b22bar(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-6);
      }

      {
         const auto v1 = softsusy::ffn(p,m1,m2,q);
         const auto v2 = flexiblesusy::f0(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-8);
      }

      {
         const auto v1 = softsusy::gfn(p,m1,m2,q);
         const auto v2 = flexiblesusy::g0(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-8);
      }

      {
         const auto v1 = softsusy::hfn(p,m1,m2,q);
         const auto v2 = flexiblesusy::h0(p2,m12,m22,q2);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-7);
      }

      {
         const auto v1 = softsusy::c0(p,m1,m2);
         const auto v2 = flexiblesusy::c0(p2,m12,m22);
         BOOST_CHECK_CLOSE_FRACTION(v1, v2, 1e-7);
      }

      {
         const auto v1 = softsusy::d0(p,m1,m2,q2);
         const auto v2 = flexiblesusy::d0(p2,m12,m22,q2);

         // not so precise because limits of equal masse are not implemented
         if (!(v1 < 1e-4 && v1 < 1e-4))
            BOOST_CHECK_CLOSE_FRACTION(v1, v2, 2e-3);
      }

      {
         const auto v1 = softsusy::d27(p,m1,m2,q2);
         const auto v2 = flexiblesusy::d27(p2,m12,m22,q2);

         // not so precise because limits of equal masse are not implemented
         if (!(v1 < 1e-7 && v1 < 1e-7))
            BOOST_CHECK_CLOSE_FRACTION(v1, v2, 2e-3);
      }
   }
}
