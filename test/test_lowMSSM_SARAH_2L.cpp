#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowMSSM_SARAH_2L

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <cmath>

#include "mssm_twoloophiggs.hpp"
#include "mssm_twoloophiggs.h"
#include "lowMSSM_mass_eigenstates.hpp"
#include "lowMSSM_slha_io.hpp"

using namespace std;
using namespace flexiblesusy;

constexpr double sqr(double x) { return x*x; }

void setup_lowMSSM(lowMSSM_mass_eigenstates &mssm)
{
   flexiblesusy::lowMSSM_slha_io mssm_slha;

   mssm_slha.read_from_file("test/SPheno.spc.MSSM");
   mssm_slha.fill(mssm);

   /*
      At the moment, SPheno is using the solutions of the full tree-level EWSB eqs. (g1,g2 non-zero)
      to calculate the 2-loop expressions and shifts with g1=g2=0. To allow comparison we therefore
      first call solve_ewsb_tree_level() before setting g1,g2 to zero, and in particular we then call
      calculate_current_DRbar_masses() instead of calculate_DRbar_masses(), since the latter
      recalculates the tree-level EWSB parameters for the masses.
      When this gets fixed in SPheno, this should also be changed.
   */
   mssm.solve_ewsb_tree_level();

   mssm.set_g1(0);
   mssm.set_g2(0);

   mssm.calculate_current_DRbar_masses();
   mssm.calculate_vertices();
}

BOOST_AUTO_TEST_CASE( lowMSSM_SARAH_2L_SPheno_comparison )
{
   lowMSSM_mass_eigenstates mssm;
   setup_lowMSSM(mssm);

   /*
      Note: the numbers from SPheno were generated using the MSSM LesHouches file
            packaged with SARAH/SPheno. Also, the variable 'epscouplings' in
            'TwoLoopMasses/2LPole_MSSM.f90' of the generated SPheno code has to be set
            to a lower value (here 10^-12), since otherwise some expressions are discarded
   */

   std::array<std::complex<double>, 2> tadpole2L{};
   tadpole2L[0] = mssm.tadpole_hh_2loop(0);
   tadpole2L[1] = mssm.tadpole_hh_2loop(1);


   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[0]),4158292.4424673738, 3e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[1]),38787424.841057681, 3e-9);


   BOOST_CHECK_SMALL(Im(tadpole2L[0]),1e-8);
   BOOST_CHECK_SMALL(Im(tadpole2L[1]),1e-8);

   std::array<std::complex<double>, 2> tadpole2Lshift{};
   tadpole2Lshift[0] = mssm.tadpole_shift_hh_2loop(0);
   tadpole2Lshift[1] = mssm.tadpole_shift_hh_2loop(1);


   BOOST_CHECK_SMALL(Re(tadpole2Lshift[0]), 1e-10);
   BOOST_CHECK_SMALL(Re(tadpole2Lshift[1]), 1e-10);
   BOOST_CHECK_SMALL(Im(tadpole2Lshift[0]), 1e-10);
   BOOST_CHECK_SMALL(Im(tadpole2Lshift[1]), 1e-10);


   auto self_energy_hh_2L = mssm.self_energy_hh_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(0,0)),-30175.548756820193, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(0,1)), 23472.298901058326, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,0)), 23472.298901058333, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,1)), 155865.86700677092, 2e-9);



   auto self_energy_shift_hh_2L = mssm.self_energy_shift_hh_2loop(sqr(125));

   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(0,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(0,1)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(1,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(1,1)), 1e-10);



   auto self_energy_Ah_2L = mssm.self_energy_Ah_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,0)),-30175.529388820458, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,1)),-23539.703290778729, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,0)),-23539.703290778729, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,1)), 158570.19956955794, 2e-9);



   auto self_energy_shift_Ah_2L = mssm.self_energy_shift_Ah_2loop(sqr(125));

   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(0,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(0,1)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(1,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(1,1)), 1e-10);

}


BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atat_atab_abab_SARAH_literature )
{
   lowMSSM_mass_eigenstates mssm;

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

   Eigen::Matrix<double, 2, 2> self_energy_hh_atat_atab_abab_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_at_at_mssm(
      sqr(mssm.get_MFu(2)), sqr(mssm.get_MFd(2)), mssm.get_M2Ah(1), mst1, mst2, msb1, msb2,
      std::sin(thetat), std::cos(thetat), std::sin(thetab), std::cos(thetab), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 2> self_energy_Ah_atat_atab_abab_literature = flexiblesusy::mssm_twoloophiggs::self_energy_pseudoscalar_2loop_at_at_mssm(
      sqr(mssm.get_MFu(2)), sqr(mssm.get_MFd(2)), mssm.get_M2Ah(1), mst1, mst2, msb1, msb2,
      std::sin(thetat), std::cos(thetat), std::sin(thetab), std::cos(thetab), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 1> tadpole_atat_atab_abab_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_at_at_mssm(
      sqr(mssm.get_MFu(2)), sqr(mssm.get_MFd(2)), mssm.get_M2Ah(1), mst1, mst2, msb1, msb2,
      std::sin(thetat), std::cos(thetat), std::sin(thetab), std::cos(thetab), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 2> self_energy_hh_atat_atab_abab_sarah = (mssm.self_energy_hh_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 2> self_energy_Ah_atat_atab_abab_sarah = (mssm.self_energy_Ah_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 1>  tadpole_atat_atab_abab_sarah;
   tadpole_atat_atab_abab_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atat_atab_abab_sarah(0,0), self_energy_hh_atat_atab_abab_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atat_atab_abab_sarah(0,1), self_energy_hh_atat_atab_abab_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atat_atab_abab_sarah(1,1), self_energy_hh_atat_atab_abab_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atat_atab_abab_sarah(0,0), self_energy_Ah_atat_atab_abab_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atat_atab_abab_sarah(0,1), self_energy_Ah_atat_atab_abab_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atat_atab_abab_sarah(1,1), self_energy_Ah_atat_atab_abab_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atat_atab_abab_sarah(0), tadpole_atat_atab_abab_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atat_atab_abab_sarah(1), tadpole_atat_atab_abab_literature(1), 1e-8);

}



BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atau_atau_SARAH_literature )
{
   lowMSSM_mass_eigenstates mssm;

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

   Eigen::Matrix<double, 2, 2> self_energy_hh_atau_atau_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_atau_atau_mssm(
      sqr(mssm.get_MFe(2)), mssm.get_M2Ah(1), msv2, mstau1, mstau2,
      std::sin(thetatau), std::cos(thetatau), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()),0);

   Eigen::Matrix<double, 2, 2> self_energy_Ah_atau_atau_literature = flexiblesusy::mssm_twoloophiggs::self_energy_pseudoscalar_2loop_atau_atau_mssm(
      sqr(mssm.get_MFe(2)), mssm.get_M2Ah(1), msv2, mstau1, mstau2,
      std::sin(thetatau), std::cos(thetatau), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()),0);

   Eigen::Matrix<double, 2, 1> tadpole_atau_atau_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_atau_atau_mssm(
      sqr(mssm.get_MFe(2)), mssm.get_M2Ah(1), msv2, mstau1, mstau2,
      std::sin(thetatau), std::cos(thetatau), sqr(mssm.get_scale()),
      -mssm.get_Mu(), mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()));

   Eigen::Matrix<double, 2, 2> self_energy_hh_atau_atau_sarah = (mssm.self_energy_hh_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 2> self_energy_Ah_atau_atau_sarah = (mssm.self_energy_Ah_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 1>  tadpole_atau_atau_sarah;
   tadpole_atau_atau_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atau_atau_sarah(0,0), self_energy_hh_atau_atau_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atau_atau_sarah(0,1), self_energy_hh_atau_atau_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atau_atau_sarah(1,1), self_energy_hh_atau_atau_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atau_atau_sarah(0,0), self_energy_Ah_atau_atau_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atau_atau_sarah(0,1), self_energy_Ah_atau_atau_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atau_atau_sarah(1,1), self_energy_Ah_atau_atau_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atau_atau_sarah(0), tadpole_atau_atau_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atau_atau_sarah(1), tadpole_atau_atau_literature(1), 1e-8);

}

BOOST_AUTO_TEST_CASE( MSSM_higgs_2loop_atas_SARAH_literature )
{
   lowMSSM_mass_eigenstates mssm;

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

   Eigen::Matrix<double, 2, 2> self_energy_hh_atas_literature = flexiblesusy::mssm_twoloophiggs::self_energy_higgs_2loop_at_as_mssm(
      sqr(mssm.get_MFu(2)), mssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(mssm.get_scale()), -mssm.get_Mu(),
      mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()), mssm.get_g3(), 0);

   Eigen::Matrix<double, 2, 2> self_energy_Ah_atas_literature = flexiblesusy::mssm_twoloophiggs::self_energy_pseudoscalar_2loop_at_as_mssm(
      sqr(mssm.get_MFu(2)), mssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(mssm.get_scale()), -mssm.get_Mu(),
      mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()), mssm.get_g3(), 0);

   Eigen::Matrix<double, 2, 1> tadpole_atas_literature = flexiblesusy::mssm_twoloophiggs::tadpole_higgs_2loop_at_as_mssm(
      sqr(mssm.get_MFu(2)), mssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(mssm.get_scale()), -mssm.get_Mu(),
      mssm.get_vu()/mssm.get_vd(), sqr(mssm.get_vu())+sqr(mssm.get_vd()), mssm.get_g3());

   Eigen::Matrix<double, 2, 2> self_energy_hh_atas_sarah = (mssm.self_energy_hh_2loop(sqr(125))).real();
   Eigen::Matrix<double, 2, 2> self_energy_Ah_atas_sarah = (mssm.self_energy_Ah_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 1>  tadpole_atas_sarah;
   tadpole_atas_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   mssm.set_g3(0);
   mssm.calculate_vertices();

   Eigen::Matrix<double, 2, 2> self_energy_hh_atat_sarah = (mssm.self_energy_hh_2loop(sqr(125))).real();
   Eigen::Matrix<double, 2, 2> self_energy_Ah_atat_sarah = (mssm.self_energy_Ah_2loop(sqr(125))).real();

   Eigen::Matrix<double, 2, 1>  tadpole_atat_sarah;
   tadpole_atat_sarah << (mssm.tadpole_hh_2loop(0)).real()/mssm.get_vd(), (mssm.tadpole_hh_2loop(1)).real()/mssm.get_vu();

   self_energy_hh_atas_sarah -= self_energy_hh_atat_sarah;
   self_energy_Ah_atas_sarah -= self_energy_Ah_atat_sarah;

   tadpole_atas_sarah -= tadpole_atat_sarah;

   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_sarah(0,0), self_energy_hh_atas_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_sarah(0,1), self_energy_hh_atas_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_sarah(1,1), self_energy_hh_atas_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_sarah(0,0), self_energy_Ah_atas_literature(0,0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_sarah(0,1), self_energy_Ah_atas_literature(0,1), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_sarah(1,1), self_energy_Ah_atas_literature(1,1), 1e-8);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_sarah(0), tadpole_atas_literature(0), 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_sarah(1), tadpole_atas_literature(1), 1e-8);

}
