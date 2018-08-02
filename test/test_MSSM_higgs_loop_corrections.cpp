#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_MSSM_higgs_loop_corrections

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <cmath>

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"
#include "MSSM_mass_eigenstates.hpp"

using namespace std;
using namespace flexiblesusy;

constexpr double sqr(double x) { return x*x; }

void setup_MSSM(MSSM_mass_eigenstates &mssm)
{
   const double g1 = 3.6307180231547609E-01*std::sqrt(5./3.);
   const double g2 = 0.65056172722216188;
   const double g3 = 1.0906350505392199;
   const double scale = 1000;
   const double vd = 24.220102076678143;
   const double vu = 242.20102076678143 ;
   const double M3 = 4.61039475E+02;
   const double mhu2= -375770.57664860087;
   const double mhd2= 489612.35827468964;


   Eigen::Matrix<double,3,3>
      Yu = Eigen::Matrix<double,3,3>::Zero(),
      Yd = Eigen::Matrix<double,3,3>::Zero(),
      Ye = Eigen::Matrix<double,3,3>::Zero(),
      Tu = Eigen::Matrix<double,3,3>::Zero(),
      Td = Eigen::Matrix<double,3,3>::Zero(),
      Te = Eigen::Matrix<double,3,3>::Zero(),
      mu2 = Eigen::Matrix<double,3,3>::Zero(),
      md2 = Eigen::Matrix<double,3,3>::Zero(),
      me2 = Eigen::Matrix<double,3,3>::Zero(),
      ml2 = Eigen::Matrix<double,3,3>::Zero(),
      mq2 = Eigen::Matrix<double,3,3>::Zero();

    Yu(0,0) =   7.03427453E-06;
    Yu(0,1) =   0.00000000E+00;
    Yu(0,2) =   0.00000000E+00;
    Yu(1,0) =   0.00000000E+00;
    Yu(1,1) =   3.58407618E-03;
    Yu(1,2) =   0.00000000E+00;
    Yu(2,0) =   0.00000000E+00;
    Yu(2,1) =   0.00000000E+00;
    Yu(2,2) =   8.71236072E-01;
    Yd(0,0) =   1.20677587E-04;
    Yd(0,1) =   0.00000000E+00;
    Yd(0,2) =   0.00000000E+00;
    Yd(1,0) =   0.00000000E+00;
    Yd(1,1) =   2.28840295E-03;
    Yd(1,2) =   0.00000000E+00;
    Yd(2,0) =   0.00000000E+00;
    Yd(2,1) =   0.00000000E+00;
    Yd(2,2) =   1.22457542E-01;
    Ye(0,0) =   2.82931110E-05;
    Ye(0,1) =   0.00000000E+00;
    Ye(0,2) =   0.00000000E+00;
    Ye(1,0) =   0.00000000E+00;
    Ye(1,1) =   5.85012055E-03;
    Ye(1,2) =   0.00000000E+00;
    Ye(2,0) =   0.00000000E+00;
    Ye(2,1) =   0.00000000E+00;
    Ye(2,2) =   9.84019104E-02;
    Tu(0,0) =   -5.77790272E-03;
    Tu(0,1) =   0.00000000E+00;
    Tu(0,2) =   0.00000000E+00;
    Tu(1,0) =   0.00000000E+00;
    Tu(1,1) =   -2.44596425E+00;
    Tu(1,2) =   0.00000000E+00;
    Tu(2,0) =   0.00000000E+00;
    Tu(2,1) =   0.00000000E+00;
    Tu(2,2) =   -4.14385376E+02;
    Td(0,0) =   -1.17347423E-01;
    Td(0,1) =   0.00000000E+00;
    Td(0,2) =   0.00000000E+00;
    Td(1,0) =   0.00000000E+00;
    Td(1,1) =   -2.46428705E+00;
    Td(1,2) =   0.00000000E+00;
    Td(2,0) =   0.00000000E+00;
    Td(2,1) =   0.00000000E+00;
    Td(2,2) =   -1.10182553E+02;
    Te(0,0) =   -1.24488551E-02;
    Te(0,1) =   0.00000000E+00;
    Te(0,2) =   0.00000000E+00;
    Te(1,0) =   0.00000000E+00;
    Te(1,1) =   -2.57394779E+00;
    Te(1,2) =   0.00000000E+00;
    Te(2,0) =   0.00000000E+00;
    Te(2,1) =   0.00000000E+00;
    Te(2,2) =   -4.29498493E+01;
    mq2(0,0) =   2.15013979E+05;
    mq2(0,1) =   0.00000000E+00;
    mq2(0,2) =   0.00000000E+00;
    mq2(1,0) =   0.00000000E+00;
    mq2(1,1) =   2.15012104E+05;
    mq2(1,2) =   0.00000000E+00;
    mq2(2,0) =   0.00000000E+00;
    mq2(2,1) =   0.00000000E+00;
    mq2(2,2) =   1.65139305E+05;
    ml2(0,0) =   5.74821153E+04;
    ml2(0,1) =   0.00000000E+00;
    ml2(0,2) =   0.00000000E+00;
    ml2(1,0) =   0.00000000E+00;
    ml2(1,1) =   5.74790645E+04;
    ml2(1,2) =   0.00000000E+00;
    ml2(2,0) =   0.00000000E+00;
    ml2(2,1) =   0.00000000E+00;
    ml2(2,2) =   5.66227938E+04;
    md2(0,0) =   2.01670449E+05;
    md2(0,1) =   0.00000000E+00;
    md2(0,2) =   0.00000000E+00;
    md2(1,0) =   0.00000000E+00;
    md2(1,1) =   2.01668819E+05;
    md2(1,2) =   0.00000000E+00;
    md2(2,0) =   0.00000000E+00;
    md2(2,1) =   0.00000000E+00;
    md2(2,2) =   1.98068542E+05;
    mu2(0,0) =   2.03033322E+05;
    mu2(0,1) =   0.00000000E+00;
    mu2(0,2) =   0.00000000E+00;
    mu2(1,0) =   0.00000000E+00;
    mu2(1,1) =   2.03031160E+05;
    mu2(1,2) =   0.00000000E+00;
    mu2(2,0) =   0.00000000E+00;
    mu2(2,1) =   0.00000000E+00;
    mu2(2,2) =   1.04647252E+05;
    me2(0,0) =   4.52752497E+04;
    me2(0,1) =   0.00000000E+00;
    me2(0,2) =   0.00000000E+00;
    me2(1,0) =   0.00000000E+00;
    me2(1,1) =   4.52690899E+04;
    me2(1,2) =   0.00000000E+00;
    me2(2,0) =   0.00000000E+00;
    me2(2,1) =   0.00000000E+00;
    me2(2,2) =   4.35399223E+04;

   mssm.do_calculate_sm_pole_masses(true);
   mssm.do_calculate_bsm_pole_masses(true);
   mssm.set_pole_mass_loop_order(2);
   mssm.set_ewsb_loop_order(2);

   // set parameters
   mssm.set_g1(g1);
   mssm.set_g2(g2);
   mssm.set_g3(g3);
   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_TYu(Tu);
   mssm.set_TYd(Td);
   mssm.set_TYe(Te);
   mssm.set_mu2(mu2);
   mssm.set_md2(md2);
   mssm.set_ml2(ml2);
   mssm.set_mq2(mq2);
   mssm.set_me2(me2);
   mssm.set_scale(scale);
   mssm.set_vd(vd);
   mssm.set_vu(vu);
   mssm.set_mHu2(mhu2);
   mssm.set_mHd2(mhd2);
   mssm.set_MassG(M3);

   mssm.solve_ewsb_tree_level();
   mssm.calculate_DRbar_masses();
   mssm.calculate_vertices();
}



BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atat_atab_abab_SARAH_literature )
{
   MSSM_mass_eigenstates mssm;

   setup_MSSM(mssm);

   double mst1, mst2, thetat;
   double msb1, msb2, thetab;

   Eigen::Matrix<double,3,3>
      Yu = Eigen::Matrix<double,3,3>::Zero(),
      Yd = Eigen::Matrix<double,3,3>::Zero(),
      Ye = Eigen::Matrix<double,3,3>::Zero(),
      Tu = Eigen::Matrix<double,3,3>::Zero(),
      Td = Eigen::Matrix<double,3,3>::Zero(),
      Te = Eigen::Matrix<double,3,3>::Zero(),
      mu2 = Eigen::Matrix<double,3,3>::Zero(),
      md2 = Eigen::Matrix<double,3,3>::Zero(),
      me2 = Eigen::Matrix<double,3,3>::Zero(),
      ml2 = Eigen::Matrix<double,3,3>::Zero(),
      mq2 = Eigen::Matrix<double,3,3>::Zero();

   Yu(2,2) =    8.73058219E-01;
   Tu(2,2) =   -4.14385376E+02;
   mu2(2,2) =   1.04647252E+05;

   Yd(2,2) =    1.35439639E-01;
   Td(2,2) =   -1.10182553E+02;
   md2(2,2) =   1.98068542E+05;

   mq2(2,2) =   1.65139305E+05;

   mssm.set_g3(0);

   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_TYu(Tu);
   mssm.set_TYd(Td);
   mssm.set_TYe(Te);
   mssm.set_mu2(mu2);
   mssm.set_md2(md2);
   mssm.set_ml2(ml2);
   mssm.set_mq2(mq2);
   mssm.set_me2(me2);

   mssm.enter_gaugeless_limit();

   mssm.calculate_M2Su_3rd_generation(mst1,mst2,thetat);
   mssm.calculate_M2Sd_3rd_generation(msb1,msb2,thetab);

   Eigen::Matrix<double, 2, 2> self_energy_atat_atab_abab_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_at_at_mssm(
      sqr(mssm.get_MFu(2)), sqr(mssm.get_MFd(2)), mssm.get_M2Ah(1), mst1, mst2, msb1, msb2,
      std::sin(thetat), std::cos(thetat), std::sin(thetab), std::cos(thetab), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 1> tadpole_atat_atab_abab_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_at_at_mssm(
      sqr(mssm.get_MFu(2)), sqr(mssm.get_MFd(2)), mssm.get_M2Ah(1), mst1, mst2, msb1, msb2,
      std::sin(thetat), std::cos(thetat), std::sin(thetab), std::cos(thetab), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 2> self_energy_atat_atab_abab_sarah = (mssm.self_energy_hh_2loop(125)).real();
   Eigen::Matrix<double, 2, 1>  tadpole_atat_atab_abab_sarah;
   tadpole_atat_atab_abab_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   BOOST_CHECK_CLOSE_FRACTION(self_energy_atat_atab_abab_sarah(0,0), self_energy_atat_atab_abab_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atat_atab_abab_sarah(0,1), self_energy_atat_atab_abab_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atat_atab_abab_sarah(1,1), self_energy_atat_atab_abab_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atat_atab_abab_sarah(0), tadpole_atat_atab_abab_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atat_atab_abab_sarah(1), tadpole_atat_atab_abab_literature(1), 1e-8);

}



BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atau_atau_SARAH_literature )
{
   MSSM_mass_eigenstates mssm;

   setup_MSSM(mssm);

   double mstau1, mstau2, thetatau;
   double msv1, msv2, thetav;

   Eigen::Matrix<double,3,3>
      Yu = Eigen::Matrix<double,3,3>::Zero(),
      Yd = Eigen::Matrix<double,3,3>::Zero(),
      Ye = Eigen::Matrix<double,3,3>::Zero(),
      Tu = Eigen::Matrix<double,3,3>::Zero(),
      Td = Eigen::Matrix<double,3,3>::Zero(),
      Te = Eigen::Matrix<double,3,3>::Zero(),
      mu2 = Eigen::Matrix<double,3,3>::Zero(),
      md2 = Eigen::Matrix<double,3,3>::Zero(),
      me2 = Eigen::Matrix<double,3,3>::Zero(),
      ml2 = Eigen::Matrix<double,3,3>::Zero(),
      mq2 = Eigen::Matrix<double,3,3>::Zero();

   Ye(2,2) =    9.84019104E-02;
   Te(2,2) =   -4.29498493E+01;
   me2(2,2) =   4.35399223E+04;

   mq2(2,2) =   1.65139305E+05;
   ml2(2,2) =   5.66227938E+04;

   mssm.set_g3(0);

   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_TYu(Tu);
   mssm.set_TYd(Td);
   mssm.set_TYe(Te);
   mssm.set_mu2(mu2);
   mssm.set_md2(md2);
   mssm.set_ml2(ml2);
   mssm.set_mq2(mq2);
   mssm.set_me2(me2);

   mssm.enter_gaugeless_limit();

   mssm.calculate_M2Se_3rd_generation(mstau1,mstau2,thetatau);
   mssm.calculate_M2Sv_3rd_generation(msv1,msv2,thetav);

   Eigen::Matrix<double, 2, 2> self_energy_atau_atau_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_atau_atau_mssm(
      sqr(mssm.get_MFe(2)), mssm.get_M2Ah(1), msv2, mstau1, mstau2, std::sin(thetatau), std::cos(thetatau), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()),0);

   Eigen::Matrix<double, 2, 1> tadpole_atau_atau_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_atau_atau_mssm(
      sqr(mssm.get_MFe(2)), mssm.get_M2Ah(1), msv2, mstau1, mstau2, std::sin(thetatau), std::cos(thetatau), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 2> self_energy_atau_atau_sarah = (mssm.self_energy_hh_2loop(125)).real();
   Eigen::Matrix<double, 2, 1>  tadpole_atau_atau_sarah;
   tadpole_atau_atau_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   BOOST_CHECK_CLOSE_FRACTION(self_energy_atau_atau_sarah(0,0), self_energy_atau_atau_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atau_atau_sarah(0,1), self_energy_atau_atau_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atau_atau_sarah(1,1), self_energy_atau_atau_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atau_atau_sarah(0), tadpole_atau_atau_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atau_atau_sarah(1), tadpole_atau_atau_literature(1), 1e-8);

}

BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atas_SARAH_literature )
{
   MSSM_mass_eigenstates mssm;

   setup_MSSM(mssm);

   double mst1, mst2, thetat;

   Eigen::Matrix<double,3,3>
      Yu = Eigen::Matrix<double,3,3>::Zero(),
      Yd = Eigen::Matrix<double,3,3>::Zero(),
      Ye = Eigen::Matrix<double,3,3>::Zero(),
      Tu = Eigen::Matrix<double,3,3>::Zero(),
      Td = Eigen::Matrix<double,3,3>::Zero(),
      Te = Eigen::Matrix<double,3,3>::Zero(),
      mu2 = Eigen::Matrix<double,3,3>::Zero(),
      md2 = Eigen::Matrix<double,3,3>::Zero(),
      me2 = Eigen::Matrix<double,3,3>::Zero(),
      ml2 = Eigen::Matrix<double,3,3>::Zero(),
      mq2 = Eigen::Matrix<double,3,3>::Zero();

   Yu(2,2) =    8.73058219E-01;
   Tu(2,2) =   -4.14385376E+02;
   mu2(2,2) =   1.04647252E+05;

   mq2(2,2) =   1.65139305E+05;


   mssm.set_Yu(Yu);
   mssm.set_Yd(Yd);
   mssm.set_Ye(Ye);
   mssm.set_TYu(Tu);
   mssm.set_TYd(Td);
   mssm.set_TYe(Te);
   mssm.set_mu2(mu2);
   mssm.set_md2(md2);
   mssm.set_ml2(ml2);
   mssm.set_mq2(mq2);
   mssm.set_me2(me2);

   mssm.enter_gaugeless_limit();

   mssm.calculate_M2Su_3rd_generation(mst1,mst2,thetat);

   Eigen::Matrix<double, 2, 2> self_energy_atas_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_at_as_mssm(
      sqr(mssm.get_MFu(2)), mssm.get_MGlu(), mst1, mst2, std::sin(thetat), std::cos(thetat), sqr(mssm.get_scale()), -mssm.get_Mu(),
      mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()), mssm.get_g3(), 0);

   Eigen::Matrix<double, 2, 1> tadpole_atas_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_at_as_mssm(
      sqr(mssm.get_MFu(2)), mssm.get_MGlu(), mst1, mst2, std::sin(thetat), std::cos(thetat), sqr(mssm.get_scale()), -mssm.get_Mu(),
      mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()), mssm.get_g3());

   Eigen::Matrix<double, 2, 2> self_energy_atas_sarah = (mssm.self_energy_hh_2loop(125)).real();
   Eigen::Matrix<double, 2, 1>  tadpole_atas_sarah;
   tadpole_atas_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   mssm.set_g3(0);

   Eigen::Matrix<double, 2, 2> self_energy_atat_sarah = (mssm.self_energy_hh_2loop(125)).real();
   Eigen::Matrix<double, 2, 1>  tadpole_atat_sarah;
   tadpole_atat_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   self_energy_atas_sarah = self_energy_atas_sarah - self_energy_atat_sarah;
   tadpole_atas_sarah = tadpole_atas_sarah - tadpole_atat_sarah;

   BOOST_CHECK_CLOSE_FRACTION(self_energy_atas_sarah(0,0), self_energy_atas_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atas_sarah(0,1), self_energy_atas_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_atas_sarah(1,1), self_energy_atas_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_sarah(0), tadpole_atas_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_sarah(1), tadpole_atas_literature(1), 1e-8);

}
