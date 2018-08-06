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

   std::array<std::complex<double>, 3> tadpole2L{};
   tadpole2L[0] = mssm.tadpole_hh_2loop(0);
   tadpole2L[1] = mssm.tadpole_hh_2loop(1);


   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[0]),-4158292.4424673738, 3e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(tadpole2L[1]),-38787424.841057681, 3e-9);


   BOOST_CHECK_SMALL(Im(tadpole2L[0]),1e-8);
   BOOST_CHECK_SMALL(Im(tadpole2L[1]),1e-8);

   std::array<std::complex<double>, 3> tadpole2Lshift{};
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
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_hh_2L(1,1)), 155865.86700677092, 1e-9);



   auto self_energy_shift_hh_2L = mssm.self_energy_shift_hh_2loop(sqr(125));

   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(0,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(0,1)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(1,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_hh_2L(1,1)), 1e-10);



   auto self_energy_Ah_2L = mssm.self_energy_Ah_2loop(sqr(125));

   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,0)),-30175.529388820458, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(0,1)),-23539.703290778729, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,0)),-23539.703290778729, 1e-9);
   BOOST_CHECK_CLOSE_FRACTION(Re(self_energy_Ah_2L(1,1)), 158570.19956955794, 1e-9);



   auto self_energy_shift_Ah_2L = mssm.self_energy_shift_Ah_2loop(sqr(125));

   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(0,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(0,1)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(1,0)), 1e-10);
   BOOST_CHECK_SMALL(Re(self_energy_shift_Ah_2L(1,1)), 1e-10);

}
