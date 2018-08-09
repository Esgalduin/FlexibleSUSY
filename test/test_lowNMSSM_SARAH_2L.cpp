#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_lowNMSSM_SARAH_2L

#include <boost/test/unit_test.hpp>
#include <Eigen/Core>
#include <cmath>

#include "nmssm_twoloophiggs.hpp"
#include "lowNMSSM_mass_eigenstates.hpp"
#include "lowNMSSM_slha_io.hpp"

using namespace std;
using namespace flexiblesusy;

constexpr double sqr(double x) { return x*x; }

void setup_lowNMSSM(lowNMSSM_mass_eigenstates &nmssm)
{
   flexiblesusy::lowNMSSM_slha_io nmssm_slha;

   nmssm_slha.read_from_file("test/SPheno.spc.lowNMSSM");
   nmssm_slha.fill(nmssm);

   /*
      At the moment, SPheno is using the solutions of the full tree-level EWSB eqs. (g1,g2 non-zero)
      to calculate the 2-loop expressions and shifts with g1=g2=0. To allow comparison we therefore
      first call solve_ewsb_tree_level() before setting g1,g2 to zero, and in particular we then call
      calculate_current_DRbar_masses() instead of calculate_DRbar_masses(), since the latter
      recalculates the tree-level EWSB parameters for the masses.
      When this gets fixed in SPheno, this should also be changed.
   */
   nmssm.solve_ewsb_tree_level();

   nmssm.set_g1(0);
   nmssm.set_g2(0);

   nmssm.calculate_current_DRbar_masses();
   nmssm.calculate_vertices();
}

BOOST_AUTO_TEST_CASE( lowNMSSM_SARAH_2L_SPheno_comparison )
{
   lowNMSSM_mass_eigenstates nmssm;
   setup_lowNMSSM(nmssm);

   /*
      Note: the numbers from SPheno were generated using the NMSSM LesHouches file
            packaged with SARAH/SPheno. Also, the variable 'epscouplings' in
            'TwoLoopMasses/2LPole_NMSSM.f90' of the generated SPheno code has to be set
            to a lower value (here 10^-12), since otherwise some expressions are discarded
   */

   std::array<std::complex<double>, 3> tadpole2L{};
   tadpole2L[0] = nmssm.tadpole_hh_2loop(0);
   tadpole2L[1] = nmssm.tadpole_hh_2loop(1);
   tadpole2L[2] = nmssm.tadpole_hh_2loop(2);

   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[0]),-7624.2312505199952, 6e-8);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[1]),-279923.66466816078, 5e-8);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[2]),-118537.18503490090, 1e-8);

   BOOST_CHECK_SMALL(Im(tadpole2L[0]),1e-8);
   BOOST_CHECK_SMALL(Im(tadpole2L[1]),1e-8);
   BOOST_CHECK_SMALL(Im(tadpole2L[2]),1e-8);


   std::array<std::complex<double>, 3> tadpole2Lshift{};
   tadpole2Lshift[0] = nmssm.tadpole_shift_hh_2loop(0);
   tadpole2Lshift[1] = nmssm.tadpole_shift_hh_2loop(1);
   tadpole2Lshift[2] = nmssm.tadpole_shift_hh_2loop(2);

   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2Lshift[0]), 2609.2074859924905, 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2Lshift[1]),-2174.9954770097429, 1e-8);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2Lshift[2]), 517023.58307811257, 1e-8);


   auto self_energy_hh_2L = nmssm.self_energy_hh_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(0,0)), 390.18053400015003, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(0,1)),-116.70340608494476, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(0,2)), 6.9507043879199699, 4e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,0)),-116.70340608494482, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,1)),-2685.7488760980996, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,2)), 345.43643705820261, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(2,0)), 6.9507043879199832, 4e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(2,1)), 345.43643705820261, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(2,2)),-2585.6478361170180, 1e-7);


   auto self_energy_shift_hh_2L = nmssm.self_energy_shift_hh_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(0,0)),-242.74124811029125, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(0,1)), 34.967175267211779, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(0,2)), 15.307001245810028, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(1,0)), 34.967175267211772, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(1,1)),-12.417810951683208, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(1,2)),-7.5460174519307479, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(2,0)), 15.307001245810048, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(2,1)),-7.5460174519307541, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_hh_2L(2,2)), 2564.0796328916099, 1e-7);


   auto self_energy_Ah_2L = nmssm.self_energy_Ah_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,0)), 334.18322284522236     , 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,1)), 65.013103856659384     , 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,2)),-0.21597970977050729     , 2e-6);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,0)), 65.013103856659370     , 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,1)),-1153.5011561948452     , 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,2)),-2.1597973538196563E-002, 2e-4);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(2,0)),-0.21597970977050807     , 2e-6);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(2,1)),-2.1597973538196480E-002, 2e-4);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(2,2)), 42.745582139784716     , 1e-7);

   /*
      The shift calculation for mass matrices of pseudoscalar Higgs in SPheno is
      missing contributions from topologies with quartic vertices. Therefore the comparison
      here is with numbers generated by FlexibleSUSY.
      This should be amended, once SPheno includes those shifts.
   */

   auto self_energy_shift_Ah_2L = nmssm.self_energy_shift_Ah_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(0,0)),-216.172452159802987, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(0,1)),-32.429789843780455, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(0,2)),-0.933747008388792, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(1,0)),-32.429789843780455, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(1,1)),-12.256150890948312, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(1,2)),-0.093374700838951, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(2,0)),-0.933747008388792, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(2,1)),-0.093374700838951, 1e-7);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_shift_Ah_2L(2,2)),-40.024044178566513, 1e-7);

}



BOOST_AUTO_TEST_CASE( lowNMSSM_SARAH_2L_literature_comparison )
{
   lowNMSSM_mass_eigenstates nmssm;

   setup_lowNMSSM(nmssm);

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

   nmssm.set_Yu(Yu);
   nmssm.set_Yd(Yd);
   nmssm.set_Ye(Ye);
   nmssm.set_TYu(Tu);
   nmssm.set_TYd(Td);
   nmssm.set_TYe(Te);
   nmssm.set_mu2(mu2);
   nmssm.set_md2(md2);
   nmssm.set_ml2(ml2);
   nmssm.set_mq2(mq2);
   nmssm.set_me2(me2);

   nmssm.enter_gaugeless_limit();

   nmssm.calculate_M2Su_3rd_generation(mst1,mst2,thetat);

   Eigen::Matrix<double, 3, 3> self_energy_hh_atas_literature = flexiblesusy::nmssm_twoloophiggs::self_energy_higgs_2loop_at_as_nmssm(
      nmssm.get_MFu(2), nmssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(nmssm.get_scale()),
      nmssm.get_vu()/nmssm.get_vd(), sqr(nmssm.get_vu())+sqr(nmssm.get_vd()),nmssm.get_Lambdax(),
      flexiblesusy::Abs(flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS())),
      sqr(nmssm.get_g3()) / (4.0 * flexiblesusy::Pi), flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS()*nmssm.get_Lambdax()));

   Eigen::Matrix<double, 3, 3> self_energy_Ah_atas_literature = flexiblesusy::nmssm_twoloophiggs::self_energy_pseudoscalar_2loop_at_as_nmssm(
      nmssm.get_MFu(2), nmssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(nmssm.get_scale()),
      nmssm.get_vu()/nmssm.get_vd(), sqr(nmssm.get_vu())+sqr(nmssm.get_vd()),nmssm.get_Lambdax(),
      flexiblesusy::Abs(flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS())),
      sqr(nmssm.get_g3()) / (4.0 * flexiblesusy::Pi), flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS()*nmssm.get_Lambdax()));

   Eigen::Matrix<double, 3, 1> tadpole_atas_literature = flexiblesusy::nmssm_twoloophiggs::tadpole_higgs_2loop_at_as_nmssm(
      sqr(nmssm.get_MFu(2)), nmssm.get_MGlu(), mst1, mst2,
      std::sin(thetat), std::cos(thetat), sqr(nmssm.get_scale()),flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS()*nmssm.get_Lambdax()),
      nmssm.get_vu()/nmssm.get_vd(),sqr(nmssm.get_vu())+sqr(nmssm.get_vd()),
      nmssm.get_g3(),flexiblesusy::Abs(flexiblesusy::Re(-0.7071067811865475*nmssm.get_vS())));

   Eigen::Matrix<double, 3, 3> self_energy_hh_atas_sarah = (nmssm.self_energy_hh_2loop(sqr(125))).real();
   Eigen::Matrix<double, 3, 3> self_energy_Ah_atas_sarah = (nmssm.self_energy_Ah_2loop(sqr(125))).real();
   Eigen::Matrix<double, 3, 1>  tadpole_atas_sarah;
   tadpole_atas_sarah << (nmssm.tadpole_hh_2loop(0)).real()/nmssm.get_vd(),
                         (nmssm.tadpole_hh_2loop(1)).real()/nmssm.get_vu(),
                         (nmssm.tadpole_hh_2loop(2)).real()/nmssm.get_vS();

   nmssm.set_g3(0);
   nmssm.calculate_vertices();

   self_energy_hh_atas_sarah -= (nmssm.self_energy_hh_2loop(sqr(125))).real();
   self_energy_Ah_atas_sarah -= (nmssm.self_energy_Ah_2loop(sqr(125))).real();

   tadpole_atas_sarah(0) -= (nmssm.tadpole_hh_2loop(0)).real()/nmssm.get_vd();
   tadpole_atas_sarah(1) -= (nmssm.tadpole_hh_2loop(1)).real()/nmssm.get_vu();
   tadpole_atas_sarah(2) -= (nmssm.tadpole_hh_2loop(2)).real()/nmssm.get_vS();


   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(0,0), self_energy_hh_atas_sarah(0,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(0,1), self_energy_hh_atas_sarah(0,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(0,2), self_energy_hh_atas_sarah(0,2), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(1,0), self_energy_hh_atas_sarah(1,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(1,1), self_energy_hh_atas_sarah(1,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(1,2), self_energy_hh_atas_sarah(1,2), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(2,0), self_energy_hh_atas_sarah(2,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(2,1), self_energy_hh_atas_sarah(2,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_hh_atas_literature(2,2), self_energy_hh_atas_sarah(2,2), 2e-10);

   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(0,0), self_energy_Ah_atas_sarah(0,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(0,1), self_energy_Ah_atas_sarah(0,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(0,2), self_energy_Ah_atas_sarah(0,2), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(1,0), self_energy_Ah_atas_sarah(1,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(1,1), self_energy_Ah_atas_sarah(1,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(1,2), self_energy_Ah_atas_sarah(1,2), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(2,0), self_energy_Ah_atas_sarah(2,0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(2,1), self_energy_Ah_atas_sarah(2,1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(self_energy_Ah_atas_literature(2,2), self_energy_Ah_atas_sarah(2,2), 1e-12);

   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_literature(0), tadpole_atas_sarah(0), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_literature(1), tadpole_atas_sarah(1), 1e-12);
   BOOST_CHECK_CLOSE_FRACTION(tadpole_atas_literature(2), tadpole_atas_sarah(2), 1e-12);


}
