
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_SM_low_scale_constraint

#include <boost/test/unit_test.hpp>

#define private public

#include "SM_two_scale_model.hpp"
#include "SM_two_scale_low_scale_constraint.hpp"
#include "wrappers.hpp"
#include "ew_input.hpp"

using namespace flexiblesusy;
using namespace softsusy;

BOOST_AUTO_TEST_CASE( test_model_properties )
{
   BOOST_CHECK(SM_info::is_supersymmetric_model == false);
}

BOOST_AUTO_TEST_CASE( test_delta_alpha )
{
   SM<Two_scale> m;
   SM_input_parameters input;
   input.LambdaIN = 0.1;
   QedQcd qedqcd;

   const double vev = 246.;

   m.set_scale(91.);
   m.set_v(246.);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.calculate_DRbar_masses();

   SM_low_scale_constraint<Two_scale> constraint(&m, qedqcd);

   const double alpha_em = qedqcd.displayAlpha(ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(ALPHAS);
   const double scale = m.get_scale();

   const double delta_alpha_em_fs = constraint.calculate_delta_alpha_em(alpha_em);
   const double delta_alpha_s_fs  = constraint.calculate_delta_alpha_s(alpha_s);

   const double Mtop = m.get_MFu(2);

   // the extra singlet field does not couple electromagnetically and
   // does thus not contribute here (or below)
   const double delta_alpha_em =
      alpha_em / (2 * Pi) * (- 16./9. * log(Mtop/scale));

   // no MS-bar DR-bar conversion term appears here, because the SM
   // is renormalized in MS-bar
   const double delta_alpha_s  =
      alpha_s / (2 * Pi) * (- 2./3. * log(Mtop/scale));

   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_em_fs, delta_alpha_em, 1.0e-12);
   BOOST_CHECK_CLOSE_FRACTION(delta_alpha_s_fs , delta_alpha_s , 1.0e-12);
}

BOOST_AUTO_TEST_CASE( test_delta_spectrum )
{
   SM<Two_scale> m;
   SM_input_parameters input;
   input.LambdaIN = 0.1;
   QedQcd qedqcd;

   const double vev = 246.;
   const double g1 = 0.2;
   const double g2 = 0.3;
   const double g3 = 0.4;

   m.set_scale(91.);
   m.set_v(246.);
   m.set_g1(g1);
   m.set_g2(g2);
   m.set_g3(g3);
   m.set_Yu(2, 2, 165.0   * Sqrt(2.) / vev);
   m.set_Yd(2, 2, 2.9     * Sqrt(2.) / vev);
   m.set_Ye(2, 2, 1.77699 * Sqrt(2.) / vev);

   m.calculate_DRbar_masses();

   const double MVZ = 0.5 * vev * Sqrt(0.6*g1*g1 + g2*g2);

   BOOST_CHECK_CLOSE_FRACTION(m.get_MVZ(), MVZ, 1.0e-10);
}
