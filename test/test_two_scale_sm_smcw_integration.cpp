#include "two_scale_solver.hpp"
#include "two_scale_matching.hpp"
#include "two_scale_model.hpp"
#include "sm_two_scale.hpp"
#include "sm_two_scale_experimental_constraint.hpp"
#include "smcw_two_scale.hpp"
#include "smcw_two_scale_gut_constraint.hpp"
#include "smcw_two_scale_convergence_tester.hpp"
#include "linalg.h"
#include "coupling_monitor.hpp"
#include "error.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_two_scale_sm_smcw_integration

#include <boost/test/unit_test.hpp>

using namespace flexiblesusy;
using namespace softsusy;

#define YU StandardModel<Two_scale>::YU
#define YD StandardModel<Two_scale>::YD
#define YE StandardModel<Two_scale>::YE

class Trivial_SM_SMCW_matching_condition_upwards: public Matching<Two_scale> {
public:
   Trivial_SM_SMCW_matching_condition_upwards()
      : Matching<Two_scale>()
      , sm(0)
      , smcw(0)
      {}
   virtual ~Trivial_SM_SMCW_matching_condition_upwards() {}
   virtual void match() {
      // ensure that both models are at the matching scale
      smcw->setScale(sm->get_scale());
      // copy parameters
      smcw->setYukawaMatrix(YU, sm->displayYukawaMatrix(YU));
      smcw->setYukawaMatrix(YD, sm->displayYukawaMatrix(YD));
      smcw->setYukawaMatrix(YE, sm->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         smcw->setGaugeCoupling(i, sm->displayGaugeCoupling(i));
   }
   virtual double get_scale() const {
      return 3000;
   }
   virtual void set_models(Two_scale_model* sm_, Two_scale_model* smcw_) {
      sm = cast_model<StandardModel<Two_scale>*>(sm_);
      smcw = cast_model<StandardModelCW<Two_scale>*>(smcw_);
   }
private:
   StandardModel<Two_scale>* sm;
   StandardModelCW<Two_scale>* smcw;
};

class Trivial_SM_SMCW_matching_condition_downwards: public Matching<Two_scale> {
public:
   Trivial_SM_SMCW_matching_condition_downwards()
      : Matching<Two_scale>()
      , sm(0)
      , smcw(0)
      {}
   virtual ~Trivial_SM_SMCW_matching_condition_downwards() {}
   virtual void match() {
      // ensure that both models are at the matching scale
      BOOST_REQUIRE(sm->get_scale() == smcw->get_scale());
      // copy parameters
      sm->setYukawaMatrix(YU, smcw->displayYukawaMatrix(YU));
      sm->setYukawaMatrix(YD, smcw->displayYukawaMatrix(YD));
      sm->setYukawaMatrix(YE, smcw->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         sm->setGaugeCoupling(i, smcw->displayGaugeCoupling(i));
   }
   virtual double get_scale() const {
      return 3000;
   }
   virtual void set_models(Two_scale_model* smcw_, Two_scale_model* sm_) {
      sm = cast_model<StandardModel<Two_scale>*>(sm_);
      smcw = cast_model<StandardModelCW<Two_scale>*>(smcw_);
   }
private:
   StandardModel<Two_scale>* sm;
   StandardModelCW<Two_scale>* smcw;
};

class Dynamic_SM_SMCW_matching_condition_downwards: public Matching<Two_scale> {
public:
   Dynamic_SM_SMCW_matching_condition_downwards()
      : Matching<Two_scale>()
      , sm(0)
      , smcw(0)
      , scale(3000) // initial guess
      {}
   virtual ~Dynamic_SM_SMCW_matching_condition_downwards() {}
   virtual void match() {
      // ensure that both models are at the matching scale
      BOOST_REQUIRE(sm->get_scale() == smcw->get_scale());
      // copy parameters
      sm->setYukawaMatrix(YU, smcw->displayYukawaMatrix(YU));
      sm->setYukawaMatrix(YD, smcw->displayYukawaMatrix(YD));
      sm->setYukawaMatrix(YE, smcw->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         sm->setGaugeCoupling(i, smcw->displayGaugeCoupling(i));
      update_scale();
   }
   virtual double get_scale() const {
      return scale;
   }
   virtual void set_models(Two_scale_model* smcw_, Two_scale_model* sm_) {
      sm = cast_model<StandardModel<Two_scale>*>(sm_);
      smcw = cast_model<StandardModelCW<Two_scale>*>(smcw_);
   }
   virtual void update_scale() {
      const double new_scale = smcw->calcZprimeMass();
      if (new_scale != 0.0)
         scale = new_scale;
   }

private:
   StandardModel<Two_scale>* sm;
   StandardModelCW<Two_scale>* smcw;
   double scale; ///< dynamic matching scale
};

class Dynamic_SM_SMCW_matching_condition_upwards: public Matching<Two_scale> {
public:
   Dynamic_SM_SMCW_matching_condition_upwards()
      : Matching<Two_scale>()
      , sm(0)
      , smcw(0)
      , scale(3000) // initial guess
      {}
   virtual ~Dynamic_SM_SMCW_matching_condition_upwards() {}
   virtual void match() {
      // ensure that both models are at the matching scale
      smcw->setScale(sm->get_scale());
      // copy parameters
      smcw->setYukawaMatrix(YU, sm->displayYukawaMatrix(YU));
      smcw->setYukawaMatrix(YD, sm->displayYukawaMatrix(YD));
      smcw->setYukawaMatrix(YE, sm->displayYukawaMatrix(YE));
      for (int i = 1; i <= 3; ++i)
         smcw->setGaugeCoupling(i, sm->displayGaugeCoupling(i));
   }
   virtual double get_scale() const {
      return scale;
   }
   virtual void set_models(Two_scale_model* sm_, Two_scale_model* smcw_) {
      sm = cast_model<StandardModel<Two_scale>*>(sm_);
      smcw = cast_model<StandardModelCW<Two_scale>*>(smcw_);
   }
   virtual void update_scale() {
      const double new_scale = smcw->calcZprimeMass();
      if (new_scale != 0.0)
         scale = new_scale;
   }

private:
   StandardModel<Two_scale>* sm;
   StandardModelCW<Two_scale>* smcw;
   double scale; ///< dynamic matching scale
};

BOOST_AUTO_TEST_CASE( test_trival_matching )
{
   BOOST_MESSAGE("test if trivial matching condition leaves parameters invariant");

   const double vev = 246;
   const double root2 = sqrt(2.0);
   const double mtoprun = 165;
   const double mbrun = 2.9;
   const double mtau = 1.77699;
   const double yt = mtoprun * root2 / vev;
   const double yb = mbrun * root2 / vev;
   const double ytau = mtau * root2 / vev;

   const double MZ = 91.1876;
   const double aem = 1.0 / 127.918; // at MZ
   const double sinthWsq = 0.23122;
   const double alpha1 = 5.0 * aem / (3.0 * (1.0 - sinthWsq));
   const double alpha2 = aem / sinthWsq;
   const double alpha3 = 0.1187; // at MZ

   StandardModel<Two_scale>* sm = new StandardModel<Two_scale>();
   sm->setScale(MZ);
   sm->setYukawaElement(YU, 3, 3, yt);
   sm->setYukawaElement(YD, 3, 3, yb);
   sm->setYukawaElement(YE, 3, 3, ytau);
   sm->setGaugeCoupling(1, sqrt(4 * PI * alpha1));
   sm->setGaugeCoupling(2, sqrt(4 * PI * alpha2));
   sm->setGaugeCoupling(3, sqrt(4 * PI * alpha3));

   StandardModelCW<Two_scale>* smcw = new StandardModelCW<Two_scale>();
   smcw->setScale(MZ);

   // this trivial matching condition simply forwards the parameters
   // of one model to the other
   Trivial_SM_SMCW_matching_condition_downwards mcd;
   Trivial_SM_SMCW_matching_condition_upwards mcu;

   // create convergence tester for the CW-Standard Model
   StandardModelCW_convergence_tester convergence_tester(smcw, 0.01);

   RGFlow<Two_scale> solver;
   solver.add(&mcu, sm, smcw);
   solver.add(&mcd, smcw, sm);
   solver.set_convergence_tester(&convergence_tester);

   // run two scale solver and ensure that no errors occure
   try {
      solver.solve();
   } catch (Error& e) {
      BOOST_ERROR(e.what());
   }

   // check that g4 and lambda are 0.0
   BOOST_CHECK_EQUAL(smcw->displayGaugeCoupling(4), 0.0);
   BOOST_CHECK_EQUAL(smcw->displayLambda(), 0.0);

   // check that the SM parameters are the same in both models
   for (int i = 1; i <= 3; ++i)
      BOOST_CHECK_EQUAL(smcw->displayGaugeCoupling(i),
                        sm->displayGaugeCoupling(i));

   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YU), sm->displayYukawaMatrix(YU));
   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YD), sm->displayYukawaMatrix(YD));
   BOOST_CHECK_EQUAL(smcw->displayYukawaMatrix(YE), sm->displayYukawaMatrix(YE));

   delete smcw;
   delete sm;
}

BOOST_AUTO_TEST_CASE( test_sm_smcw_constraints )
{
   BOOST_MESSAGE("test if SM and SMCW constraints are applied correctly");

   // create Standard Model and the EW constraints
   StandardModel<Two_scale> sm;
   sm.setScale(Electroweak_constants::MZ);
   StandardModel_exp_constraint sm_ew_constraint;
   sm_ew_constraint.set_model(&sm);

   // create CW-Standard Model and the GUT constraint
   StandardModelCW<Two_scale> smcw;
   smcw.setScale(Electroweak_constants::MZ);
   const double lambda_at_mgut = 1.0;
   StandardModelCWGUTConstraint smcw_gut_constraint(1.0e12, lambda_at_mgut);
   smcw_gut_constraint.set_model(&smcw);

   // create trivial matching condition
   Trivial_SM_SMCW_matching_condition_upwards mcu;
   Trivial_SM_SMCW_matching_condition_downwards mcd;

   // create convergence tester for the CW-Standard Model
   StandardModelCW_convergence_tester convergence_tester(&smcw, 1.0e-4);

   // create two scale solver
   RGFlow<Two_scale> solver;
   solver.add(&sm_ew_constraint, &sm);
   solver.add(&mcu, &sm, &smcw);
   solver.add(&smcw_gut_constraint, &smcw);
   solver.add(&mcd, &smcw, &sm);
   solver.set_convergence_tester(&convergence_tester);

   // run two scale solver and ensure that no errors occure
   try {
      solver.solve();
   } catch (Error& e) {
      BOOST_ERROR(e.what());
   }

   // to make the parameter comparison work, run sm to the same scale as smcw
   sm.run_to(smcw.get_scale());

   // check that the SM parameters are the same in both models
   for (int i = 1; i <= 3; ++i)
      BOOST_CHECK_CLOSE(smcw.displayGaugeCoupling(i),
                        sm.displayGaugeCoupling(i), 0.002);

   for (int i = 1; i <= 3; ++i) {
      for (int k = 1; k <= 3; ++k) {
         BOOST_CHECK_CLOSE(smcw.displayYukawaElement(YU, i, k),
                           sm.displayYukawaElement(YU, i, k), 0.002);
         BOOST_CHECK_CLOSE(smcw.displayYukawaElement(YD, i, k),
                           sm.displayYukawaElement(YD, i, k), 0.002);
         BOOST_CHECK_CLOSE(smcw.displayYukawaElement(YE, i, k),
                           sm.displayYukawaElement(YE, i, k), 0.002);
      }
   }

   // get the GUT scale value
   const double gut_scale = smcw_gut_constraint.get_scale();

   // run CW-Standard Model to the GUT scale and test equality of g1
   // and g2
   smcw.run_to(gut_scale);
   const double g1_at_mgut = smcw.displayGaugeCoupling(1);
   const double g2_at_mgut = smcw.displayGaugeCoupling(2);
   const double g4_at_mgut = smcw.displayGaugeCoupling(4);
   BOOST_CHECK_CLOSE(g1_at_mgut, g2_at_mgut, 1.0e-8);
   BOOST_CHECK_CLOSE(g1_at_mgut, g4_at_mgut, 1.0e-7);

   // check that the input value lambda(MGUT) was applied correctly
   const double lambda_at_mgut_output = smcw.displayLambda();
   BOOST_CHECK_CLOSE(lambda_at_mgut, lambda_at_mgut_output, 1.0e-4);
}

BOOST_AUTO_TEST_CASE( test_sm_smcw_convergence )
{
   BOOST_MESSAGE("test if two scale solver with SM and SMCW converges");

   // create Standard Model and the EW constraints
   StandardModel<Two_scale> sm;
   sm.setScale(Electroweak_constants::MZ);
   StandardModel_exp_constraint sm_ew_constraint;
   sm_ew_constraint.set_model(&sm);

   // create CW-Standard Model and the GUT constraint
   StandardModelCW<Two_scale> smcw;
   smcw.setScale(Electroweak_constants::MZ);
   const double lambda_at_mgut = 1.0;
   StandardModelCWGUTConstraint smcw_gut_constraint(1.0e12, lambda_at_mgut);
   smcw_gut_constraint.set_model(&smcw);

   // create trivial matching condition
   Trivial_SM_SMCW_matching_condition_upwards mcu;
   Trivial_SM_SMCW_matching_condition_downwards mcd;

   // create convergence tester for the CW-Standard Model
   StandardModelCW_convergence_tester convergence_tester(&smcw, 0.01);

   // create two scale solver
   RGFlow<Two_scale> solver;
   solver.set_convergence_tester(&convergence_tester);
   solver.add(&sm_ew_constraint, &sm);
   solver.add(&mcu, &sm, &smcw);
   solver.add(&smcw_gut_constraint, &smcw);
   solver.add(&mcd, &smcw, &sm);

   // run two scale solver and ensure that no errors occure
   try {
      solver.solve();
   } catch (Error& e) {
      BOOST_ERROR(e.what());
   }
   BOOST_MESSAGE("convergence after " << solver.number_of_iterations_done()
                 << " iterations");
   BOOST_CHECK_EQUAL(solver.number_of_iterations_done(), 3u);

#if 0
   // create data: all gauge couplings at different scales
   sm.run_to(Electroweak_constants::MZ);
   smcw.run_to(3000);
   const double gut_scale = smcw_gut_constraint.get_scale();

   Coupling_monitor cm;
   Gauge_coupling_getter gcg;
   cm.run(sm  , gcg, Electroweak_constants::MZ, 3000, 50, true);
   cm.run(smcw, gcg, 3000, gut_scale, 100, true);
   cm.write_to_file("running_coupling.dat");
#endif
}

BOOST_AUTO_TEST_CASE( test_sm_smcw_dynamic_convergence )
{
   BOOST_MESSAGE("test if two scale solver with SM and SMCW and dynamic MC converges");

   // create Standard Model and the EW constraints
   StandardModel<Two_scale> sm;
   sm.setScale(Electroweak_constants::MZ);
   StandardModel_exp_constraint sm_ew_constraint;
   sm_ew_constraint.set_model(&sm);

   // create CW-Standard Model and the GUT constraint
   StandardModelCW<Two_scale> smcw;
   smcw.setScale(Electroweak_constants::MZ);
   smcw.setVs(5000.0);
   const double lambda_at_mgut = 1.0;
   StandardModelCWGUTConstraint smcw_gut_constraint(1.0e12, lambda_at_mgut);
   smcw_gut_constraint.set_model(&smcw);

   // create dynamic matching condition
   Dynamic_SM_SMCW_matching_condition_upwards mcu;
   Dynamic_SM_SMCW_matching_condition_downwards mcd;

   // create convergence tester for the CW-Standard Model
   StandardModelCW_convergence_tester convergence_tester(&smcw, 0.01);

   // create two scale solver
   RGFlow<Two_scale> solver;
   solver.set_convergence_tester(&convergence_tester);
   solver.add(&sm_ew_constraint, &sm);
   solver.add(&mcu, &sm, &smcw);
   solver.add(&smcw_gut_constraint, &smcw);
   solver.add(&mcd, &smcw, &sm);

   // run two scale solver and ensure that no errors occure
   try {
      solver.solve();
   } catch (Error& e) {
      BOOST_ERROR(e.what());
   }
   BOOST_MESSAGE("convergence after " << solver.number_of_iterations_done()
                 << " iterations");
   BOOST_CHECK_EQUAL(solver.number_of_iterations_done(), 3u);

   // check that the matching scale is approx. the Z' mass
   // (assumption: smcw is currently at the matching scale)
   mcu.update_scale();
   BOOST_CHECK_CLOSE(mcu.get_scale(), smcw.calcZprimeMass(), 1.0e-8);
   BOOST_CHECK_CLOSE(mcd.get_scale(), smcw.calcZprimeMass(), 1.0e-8);

#if 0
   // create data: all gauge couplings at different scales
   const double gut_scale = smcw_gut_constraint.get_scale();
   const double matching_scale = mc.get_scale();

   sm.run_to(Electroweak_constants::MZ);
   smcw.run_to(matching_scale);

   Coupling_monitor cm;
   Gauge_coupling_getter gcg;
   cm.run(sm  , gcg, Electroweak_constants::MZ, matching_scale, 50, true);
   cm.run(smcw, gcg, matching_scale, gut_scale, 100, true);
   cm.write_to_file("running_coupling.dat");
#endif
}
