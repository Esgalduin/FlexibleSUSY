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

/**
 * @file standard_model.cpp
 * @brief implementation of the Standard_model class
 *
 */

#include "standard_model.hpp"
#include "eigen_utils.hpp"
#include "wrappers.hpp"
#include "linalg2.hpp"
#include "lowe.h"
#include "numerics2.hpp"
#include "logger.hpp"
#include "error.hpp"
#include "root_finder.hpp"
#include "fixed_point_iterator.hpp"
#include "config.h"
#include "pv.hpp"
#include "raii.hpp"
#include "thread_pool.hpp"
#include "functors.hpp"
#include "ew_input.hpp"
#include "weinberg_angle.hpp"

#include "sm_twoloophiggs.hpp"
#include "sm_threeloophiggs.hpp"
#include "sm_fourloophiggs.hpp"
#include "sm_threeloop_as.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <memory>
#include <algorithm>

namespace flexiblesusy {

namespace standard_model_info {

   const double normalization_g1 = 0.7745966692414834;
   const double normalization_g2 = 1;
   const double normalization_g3 = 1;

   const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities = {1, 1,
      3, 1, 1, 3, 3, 3, 1, 1, 1};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_names = {"VG",
      "Hp", "Fv", "Ah", "hh", "Fd", "Fu", "Fe", "VWp", "VP", "VZ"};

   const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names = {
      "g", "H^+", "\\nu", "A^0", "h", "d", "u", "e", "W^+", "\\gamma", "Z"};

   const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names = {"g1",
      "g2", "g3", "Lambdax", "Yu(0,0)", "Yu(0,1)", "Yu(0,2)", "Yu(1,0)",
      "Yu(1,1)", "Yu(1,2)", "Yu(2,0)", "Yu(2,1)", "Yu(2,2)", "Yd(0,0)", "Yd(0,1)",
      "Yd(0,2)", "Yd(1,0)", "Yd(1,1)", "Yd(1,2)", "Yd(2,0)", "Yd(2,1)", "Yd(2,2)"
      , "Ye(0,0)", "Ye(0,1)", "Ye(0,2)", "Ye(1,0)", "Ye(1,1)", "Ye(1,2)",
      "Ye(2,0)", "Ye(2,1)", "Ye(2,2)", "mu2", "v"};

   const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names = {
      "Re(Vd(0,0))", "Im(Vd(0,0))", "Re(Vd(0,1))", "Im(Vd(0,1))", "Re(Vd(0,2))",
      "Im(Vd(0,2))", "Re(Vd(1,0))", "Im(Vd(1,0))", "Re(Vd(1,1))", "Im(Vd(1,1))",
      "Re(Vd(1,2))", "Im(Vd(1,2))", "Re(Vd(2,0))", "Im(Vd(2,0))", "Re(Vd(2,1))",
      "Im(Vd(2,1))", "Re(Vd(2,2))", "Im(Vd(2,2))", "Re(Ud(0,0))", "Im(Ud(0,0))",
      "Re(Ud(0,1))", "Im(Ud(0,1))", "Re(Ud(0,2))", "Im(Ud(0,2))", "Re(Ud(1,0))",
      "Im(Ud(1,0))", "Re(Ud(1,1))", "Im(Ud(1,1))", "Re(Ud(1,2))", "Im(Ud(1,2))",
      "Re(Ud(2,0))", "Im(Ud(2,0))", "Re(Ud(2,1))", "Im(Ud(2,1))", "Re(Ud(2,2))",
      "Im(Ud(2,2))", "Re(Vu(0,0))", "Im(Vu(0,0))", "Re(Vu(0,1))", "Im(Vu(0,1))",
      "Re(Vu(0,2))", "Im(Vu(0,2))", "Re(Vu(1,0))", "Im(Vu(1,0))", "Re(Vu(1,1))",
      "Im(Vu(1,1))", "Re(Vu(1,2))", "Im(Vu(1,2))", "Re(Vu(2,0))", "Im(Vu(2,0))",
      "Re(Vu(2,1))", "Im(Vu(2,1))", "Re(Vu(2,2))", "Im(Vu(2,2))", "Re(Uu(0,0))",
      "Im(Uu(0,0))", "Re(Uu(0,1))", "Im(Uu(0,1))", "Re(Uu(0,2))", "Im(Uu(0,2))",
      "Re(Uu(1,0))", "Im(Uu(1,0))", "Re(Uu(1,1))", "Im(Uu(1,1))", "Re(Uu(1,2))",
      "Im(Uu(1,2))", "Re(Uu(2,0))", "Im(Uu(2,0))", "Re(Uu(2,1))", "Im(Uu(2,1))",
      "Re(Uu(2,2))", "Im(Uu(2,2))", "Re(Ve(0,0))", "Im(Ve(0,0))", "Re(Ve(0,1))",
      "Im(Ve(0,1))", "Re(Ve(0,2))", "Im(Ve(0,2))", "Re(Ve(1,0))", "Im(Ve(1,0))",
      "Re(Ve(1,1))", "Im(Ve(1,1))", "Re(Ve(1,2))", "Im(Ve(1,2))", "Re(Ve(2,0))",
      "Im(Ve(2,0))", "Re(Ve(2,1))", "Im(Ve(2,1))", "Re(Ve(2,2))", "Im(Ve(2,2))",
      "Re(Ue(0,0))", "Im(Ue(0,0))", "Re(Ue(0,1))", "Im(Ue(0,1))", "Re(Ue(0,2))",
      "Im(Ue(0,2))", "Re(Ue(1,0))", "Im(Ue(1,0))", "Re(Ue(1,1))", "Im(Ue(1,1))",
      "Re(Ue(1,2))", "Im(Ue(1,2))", "Re(Ue(2,0))", "Im(Ue(2,0))", "Re(Ue(2,1))",
      "Im(Ue(2,1))", "Re(Ue(2,2))", "Im(Ue(2,2))", "ZZ(0,0)", "ZZ(0,1)", "ZZ(1,0)"
      , "ZZ(1,1)"};

   const std::string model_name = "SM";

void print(std::ostream& ostr)
{
   ostr
      << "Model information\n"
      << "=================\n"
      << "Model name:                " << model_name << '\n'
      << "Is a low-energy model:     "
      << (is_low_energy_model ? "yes" : "no") << '\n'
      << "Is a supersymmetric model: "
      << (is_supersymmetric_model ? "yes" : "no") << '\n'
      << "Number of multiplets:      " << NUMBER_OF_PARTICLES << '\n'
      << "Number of parameters:      " << NUMBER_OF_PARAMETERS << '\n'
      ;

   ostr << "\n"
      "Multiplets:                ";
   for (int i = 0; i < NUMBER_OF_PARTICLES; i++) {
      ostr << particle_names[i]
           << '[' << particle_multiplicities[i] << ']';
      if (i + 1 < NUMBER_OF_PARTICLES)
         ostr << ", ";
   }

   ostr << "\n\n"
      "Parameters:                ";
   for (int i = 0; i < NUMBER_OF_PARAMETERS; i++) {
      ostr << parameter_names[i];
      if (i + 1 < NUMBER_OF_PARAMETERS)
         ostr << ", ";
   }
   ostr << '\n';
}

} // namespace standard_model_info

namespace standard_model {

const int Standard_model::numberOfParameters;

#define PHYSICAL(parameter) physical.parameter

#define HIGGS_2LOOP_CORRECTION_AT_AS       loop_corrections.higgs_at_as
#define HIGGS_2LOOP_CORRECTION_AB_AS       loop_corrections.higgs_ab_as
#define HIGGS_2LOOP_CORRECTION_AT_AT       loop_corrections.higgs_at_at
#define HIGGS_2LOOP_CORRECTION_ATAU_ATAU   loop_corrections.higgs_atau_atau
#define TOP_POLE_QCD_CORRECTION            loop_corrections.top_qcd
#define HIGGS_3LOOP_CORRECTION_AT_AS_AS    loop_corrections.higgs_at_as_as
#define HIGGS_3LOOP_CORRECTION_AB_AS_AS    loop_corrections.higgs_ab_as_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AS    loop_corrections.higgs_at_at_as
#define HIGGS_3LOOP_CORRECTION_AT_AT_AT    loop_corrections.higgs_at_at_at
#define HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS loop_corrections.higgs_at_as_as_as

Standard_model::Standard_model()
{
   set_number_of_parameters(numberOfParameters);
   set_thresholds(3);
}

Standard_model::Standard_model(double scale_, double loops_, double thresholds_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double mu2_, double v_)
   : Beta_function()
   , g1(g1_), g2(g2_), g3(g3_), Lambdax(Lambdax_), Yu(Yu_), Yd(Yd_), Ye(Ye_)
   , mu2(mu2_), v(v_)
{
   set_number_of_parameters(numberOfParameters);
   set_scale(scale_);
   set_loops(loops_);
   set_thresholds(thresholds_);
}

void Standard_model::do_force_output(bool flag)
{
   force_output = flag;
}

bool Standard_model::do_force_output() const
{
   return force_output;
}

void Standard_model::set_ewsb_loop_order(int loop_order)
{
   ewsb_loop_order = loop_order;
}

void Standard_model::set_loop_corrections(const Loop_corrections& loop_corrections_)
{
   loop_corrections = loop_corrections_;
}

const Loop_corrections& Standard_model::get_loop_corrections() const
{
   return loop_corrections;
}

void Standard_model::set_threshold_corrections(const Threshold_corrections& tc)
{
   threshold_corrections = tc;
}

const Threshold_corrections& Standard_model::get_threshold_corrections() const
{
   return threshold_corrections;
}

int Standard_model::get_number_of_ewsb_iterations() const
{
   return static_cast<int>(std::abs(-log10(ewsb_iteration_precision) * 10));
}

int Standard_model::get_number_of_mass_iterations() const
{
   return static_cast<int>(std::abs(-log10(precision) * 10));
}

void Standard_model::set_precision(double precision_)
{
   precision = precision_;
   ewsb_iteration_precision = precision_;
}

void Standard_model::set_pole_mass_loop_order(int loop_order)
{
   pole_mass_loop_order = loop_order;
}

int Standard_model::get_pole_mass_loop_order() const
{
   return pole_mass_loop_order;
}

void Standard_model::set_ewsb_iteration_precision(double precision)
{
   ewsb_iteration_precision = precision;
}

double Standard_model::get_ewsb_iteration_precision() const
{
   return ewsb_iteration_precision;
}

double Standard_model::get_precision() const
{
   return precision;
}

double Standard_model::get_ewsb_loop_order() const
{
   return ewsb_loop_order;
}

const Standard_model_physical& Standard_model::get_physical() const
{
   return physical;
}

Standard_model_physical& Standard_model::get_physical()
{
   return physical;
}

void Standard_model::set_physical(const Standard_model_physical& physical_)
{
   physical = physical_;
}

const Problems& Standard_model::get_problems() const
{
   return problems;
}

Problems& Standard_model::get_problems()
{
   return problems;
}

double Standard_model::get_MVG() const { return AbsSqrt(M2VG); }
double Standard_model::get_MHp() const { return AbsSqrt(M2Hp); }
double Standard_model::get_MAh() const { return AbsSqrt(M2Ah); }
double Standard_model::get_Mhh() const { return AbsSqrt(M2hh); }
double Standard_model::get_MVWp() const { return AbsSqrt(M2VWp); }
double Standard_model::get_MVP() const { return AbsSqrt(M2VP); }
double Standard_model::get_MVZ() const { return AbsSqrt(M2VZ); }

/**
 * Method which calculates the tadpoles at the current loop order.
 *
 * @return array of tadpoles
 */
Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::tadpole_equations() const
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole(
      Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   tadpole[0] = get_ewsb_eq_hh_1();

   if (ewsb_loop_order > 0) {
      tadpole[0] -= Re(tadpole_hh_1loop());
   }

   return tadpole;
}

/**
 * This method solves the EWSB conditions iteratively, trying several
 * root finding methods until a solution is found.
 */
int Standard_model::solve_ewsb_iteratively()
{
   auto ewsb_stepper = [this](const EWSB_vector_t& ewsb_pars) -> EWSB_vector_t {
      set_mu2(ewsb_pars(0));
      if (get_ewsb_loop_order() > 0)
         calculate_DRbar_masses();
      return ewsb_step();
   };
   auto tadpole_stepper = [this](const EWSB_vector_t& ewsb_pars) -> EWSB_vector_t {
      set_mu2(ewsb_pars(0));
      if (get_ewsb_loop_order() > 0)
         calculate_DRbar_masses();
      return tadpole_equations();
   };

   std::unique_ptr<EWSB_solver> solvers[] = {
      std::unique_ptr<EWSB_solver>(new Fixed_point_iterator<number_of_ewsb_equations, fixed_point_iterator::Convergence_tester_relative>(ewsb_stepper, get_number_of_ewsb_iterations(), fixed_point_iterator::Convergence_tester_relative(ewsb_iteration_precision))),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, get_number_of_ewsb_iterations(), ewsb_iteration_precision, Root_finder<number_of_ewsb_equations>::GSLHybridS)),
      std::unique_ptr<EWSB_solver>(new Root_finder<number_of_ewsb_equations>(tadpole_stepper, get_number_of_ewsb_iterations(), ewsb_iteration_precision, Root_finder<number_of_ewsb_equations>::GSLBroyden))
   };

   const auto x_init(ewsb_initial_guess());

   VERBOSE_MSG("\t\tSolving EWSB equations ...\n"
               "\t\tInitial guess: x_init = " << x_init.transpose());

   int status;
   for (auto& solver: solvers) {
      VERBOSE_MSG("\t\t\tStarting EWSB iteration using " << solver->name());
      status = solve_ewsb_iteratively_with(solver.get(), x_init);
      if (status == EWSB_solver::SUCCESS) {
         VERBOSE_MSG("\t\t\t" << solver->name() << " finished successfully!");
         break;
      }
#ifdef ENABLE_VERBOSE
      else {
         WARNING("\t\t\t" << solver->name() << " could not find a solution!"
                 " (requested precision: " << ewsb_iteration_precision << ")");
      }
#endif
   }

   if (status == EWSB_solver::SUCCESS) {
      problems.unflag_no_ewsb();
   } else {
      problems.flag_no_ewsb();
#ifdef ENABLE_VERBOSE
      WARNING("\t\t\tCould not find a solution to the EWSB equations!"
              " (requested precision: " << ewsb_iteration_precision << ")");
#endif
   }

   return status;
}

/**
 * Solves EWSB equations with given EWSB solver
 *
 * @param solver EWSB solver
 * @param x_init initial values
 *
 * @return status of the EWSB solver
 */
int Standard_model::solve_ewsb_iteratively_with(
   EWSB_solver* solver,
   const Eigen::Matrix<double, number_of_ewsb_equations, 1>& x_init
)
{
   const int status = solver->solve(x_init);
   const auto solution = solver->get_solution();

   mu2 = solution(0);


   return status;
}

int Standard_model::solve_ewsb_iteratively(int loop_order)
{
   // temporarily set `ewsb_loop_order' to `loop_order' and do
   // iteration
   const auto save_loop_order_raii = make_raii_save(ewsb_loop_order);
   ewsb_loop_order = loop_order;
   return solve_ewsb_iteratively();
}


int Standard_model::solve_ewsb_tree_level()
{
   int error = 0;

   const double old_mu2 = mu2;

   mu2 = Re(0.5*Lambdax*Sqr(v));

   const bool is_finite = IsFinite(mu2);

   if (!is_finite) {
      mu2 = old_mu2;
      error = 1;
   }


   return error;
}

int Standard_model::solve_ewsb_one_loop()
{
   return solve_ewsb_iteratively(1);
}

int Standard_model::solve_ewsb()
{
   VERBOSE_MSG("\t\tSolving Standard model EWSB at " << ewsb_loop_order << "-loop order");

   if (ewsb_loop_order == 0)
      return solve_ewsb_tree_level();

   return solve_ewsb_iteratively(ewsb_loop_order);
}

Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::ewsb_initial_guess()
{
   Eigen::Matrix<double, number_of_ewsb_equations, 1> x_init;

   x_init[0] = mu2;


   return x_init;
}

/**
 * Calculates EWSB output parameters including loop-corrections.
 *
 * Throws exception of type EEWSBStepFailed if new EWSB parameters are
 * inf or nan.
 *
 * @return new set of EWSB output parameters
 */
Eigen::Matrix<double, Standard_model::number_of_ewsb_equations, 1> Standard_model::ewsb_step() const
{
   std::array<double, number_of_ewsb_equations> tadpole{};
   Eigen::Matrix<double, number_of_ewsb_equations, 1> ewsb_parameters(Eigen::Matrix<double, number_of_ewsb_equations, 1>::Zero());

   if (ewsb_loop_order > 0) {
      tadpole[0] += Re(tadpole_hh_1loop());
   }

   double mu2;

   mu2 = Re((0.5*(Cube(v)*Lambdax - 2*tadpole[0]))/v);

   const bool is_finite = IsFinite(mu2);


   if (!is_finite)
      throw EEWSBStepFailed();

   ewsb_parameters[0] = mu2;


   return ewsb_parameters;
}

void Standard_model::print(std::ostream& ostr) const
{
   ostr << "========================================\n"
           "Standard model Q = " << get_scale() << "GeV \n"
           "========================================\n";
   ostr << "running parameters:\n";
   ostr << "g1 = " << g1 << '\n';
   ostr << "g2 = " << g2 << '\n';
   ostr << "g3 = " << g3 << '\n';
   ostr << "Lambdax = " << Lambdax << '\n';
   ostr << "Yu = " << Yu << '\n';
   ostr << "Yd = " << Yd << '\n';
   ostr << "Ye = " << Ye << '\n';
   ostr << "mu2 = " << mu2 << '\n';
   ostr << "v = " << v << '\n';
   ostr << "----------------------------------------\n"
           "tree-level MSbar masses:\n"
           "----------------------------------------\n";
   ostr << "M2VG = " << M2VG << '\n';
   ostr << "M2Hp = " << M2Hp << '\n';
   ostr << "MFv = " << MFv.transpose() << '\n';
   ostr << "M2Ah = " << M2Ah << '\n';
   ostr << "M2hh = " << M2hh << '\n';
   ostr << "MFd = " << MFd.transpose() << '\n';
   ostr << "MFu = " << MFu.transpose() << '\n';
   ostr << "MFe = " << MFe.transpose() << '\n';
   ostr << "M2VWp = " << M2VWp << '\n';
   ostr << "M2VP = " << M2VP << '\n';
   ostr << "M2VZ = " << M2VZ << '\n';

   ostr << "----------------------------------------\n"
           "tree-level MSbar mixing matrices:\n"
           "----------------------------------------\n";
   ostr << "Vd = " << Vd << '\n';
   ostr << "Ud = " << Ud << '\n';
   ostr << "Vu = " << Vu << '\n';
   ostr << "Uu = " << Uu << '\n';
   ostr << "Ve = " << Ve << '\n';
   ostr << "Ue = " << Ue << '\n';
   ostr << "ZZ = " << ZZ << '\n';

   physical.print(ostr);
}

/**
 * wrapper routines for passarino Veltman functions
 * @note: They take squared arguments!
 */

double Standard_model::A0(double m) const noexcept
{
   return passarino_veltman::ReA0(m, Sqr(get_scale()));
}

double Standard_model::B0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB0(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::B1(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB1(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::B00(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB00(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::B22(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReB22(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::H0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReH0(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::F0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReF0(p, m1, m2, Sqr(get_scale()));
}

double Standard_model::G0(double p, double m1, double m2) const noexcept
{
   return passarino_veltman::ReG0(p, m1, m2, Sqr(get_scale()));
}

/**
 * routine which finds the MSbar mass eigenstates and mixings.
 */
void Standard_model::calculate_DRbar_masses()
{
   const auto save_mu2_raii = make_raii_save(mu2);

   solve_ewsb_tree_level();

   calculate_M2VPVZ();
   calculate_M2VWp();
   calculate_MFe();
   calculate_MFu();
   calculate_MFd();
   calculate_M2hh();
   calculate_M2Ah();
   calculate_MFv();
   calculate_M2Hp();
   calculate_M2VG();
}

/**
 * routine which finds the pole mass eigenstates and mixings.
 */
void Standard_model::calculate_pole_masses()
{
#ifdef ENABLE_THREADS
   Thread_pool tp(std::min(std::thread::hardware_concurrency(), 9u));

   tp.run_task([this] () { calculate_M2VG_pole(); });
   tp.run_task([this] () { calculate_MFv_pole(); });
   tp.run_task([this] () { calculate_M2hh_pole(); });
   tp.run_task([this] () { calculate_M2VP_pole(); });
   tp.run_task([this] () { calculate_M2VZ_pole(); });
   tp.run_task([this] () { calculate_MFd_pole(); });
   tp.run_task([this] () { calculate_MFu_pole(); });
   tp.run_task([this] () { calculate_MFe_pole(); });
   tp.run_task([this] () { calculate_M2VWp_pole(); });

#else

   calculate_M2VG_pole();
   calculate_MFv_pole();
   calculate_M2hh_pole();
   calculate_M2VP_pole();
   calculate_M2VZ_pole();
   calculate_MFd_pole();
   calculate_MFu_pole();
   calculate_MFe_pole();
   calculate_M2VWp_pole();

#endif
}

void Standard_model::copy_DRbar_masses_to_pole_masses()
{
   PHYSICAL(M2VG) = M2VG;
   PHYSICAL(M2Hp) = M2Hp;
   PHYSICAL(MFv) = MFv;
   PHYSICAL(M2Ah) = M2Ah;
   PHYSICAL(M2hh) = M2hh;
   PHYSICAL(MFd) = MFd;
   PHYSICAL(Vd) = Vd;
   PHYSICAL(Ud) = Ud;
   PHYSICAL(MFu) = MFu;
   PHYSICAL(Vu) = Vu;
   PHYSICAL(Uu) = Uu;
   PHYSICAL(MFe) = MFe;
   PHYSICAL(Ve) = Ve;
   PHYSICAL(Ue) = Ue;
   PHYSICAL(M2VWp) = M2VWp;
   PHYSICAL(M2VP) = M2VP;
   PHYSICAL(M2VZ) = M2VZ;
   PHYSICAL(ZZ) = ZZ;
}

/**
 * Checks the pole masses for tachyons
 */
void Standard_model::check_pole_masses_for_tachyons()
{
   if (PHYSICAL(M2hh) < 0.) problems.flag_pole_tachyon(standard_model_info::hh);

}

/**
 * calculates spectrum for model once the MSbar parameters at
 * at low energies are known
 */
void Standard_model::calculate_spectrum()
{
   calculate_DRbar_masses();
   if (pole_mass_loop_order > 0)
      calculate_pole_masses();

   if (pole_mass_loop_order == 0) {
      copy_DRbar_masses_to_pole_masses();
   }

   check_pole_masses_for_tachyons();

   if (problems.have_problem() && !force_output) {
      clear_DRbar_parameters();
      physical.clear();
   }
}

void Standard_model::clear_running_parameters()
{
   g1 = 0.;
   g2 = 0.;
   g3 = 0.;
   Lambdax = 0.;
   Yu = Eigen::Matrix<double,3,3>::Zero();
   Yd = Eigen::Matrix<double,3,3>::Zero();
   Ye = Eigen::Matrix<double,3,3>::Zero();
   mu2 = 0.;
   v = 0.;
}

void Standard_model::clear_DRbar_parameters()
{
   M2VG = 0.;
   M2Hp = 0.;
   MFv = Eigen::Matrix<double,3,1>::Zero();
   M2Ah = 0.;
   M2hh = 0.;
   MFd = Eigen::Matrix<double,3,1>::Zero();
   Vd = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ud = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFu = Eigen::Matrix<double,3,1>::Zero();
   Vu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Uu = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   MFe = Eigen::Matrix<double,3,1>::Zero();
   Ve = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   Ue = Eigen::Matrix<std::complex<double>,3,3>::Zero();
   M2VWp = 0.;
   M2VP = 0.;
   M2VZ = 0.;
   ZZ = Eigen::Matrix<double,2,2>::Zero();
}

void Standard_model::clear_problems()
{
   problems.unflag_all_tachyons();
}

void Standard_model::clear()
{
   reset();
   clear_running_parameters();
   clear_DRbar_parameters();
   physical.clear();
   problems.clear();
}

std::string Standard_model::name() const
{
   return "Standard model";
}

void Standard_model::initialise_from_input(const softsusy::QedQcd& qedqcd_)
{
   const double scale = get_scale();
   auto qedqcd = qedqcd_;

   // initial guess
   qedqcd.to(qedqcd.displayPoleMZ());
   initial_guess_for_parameters(qedqcd);
   run_to(qedqcd.displayPoleMZ());

   // determine Standard model parameters iteratively
   Standard_model old(*this);

   bool converged = false;

   const double alpha_em = qedqcd.displayAlpha(softsusy::ALPHA);
   const double alpha_s  = qedqcd.displayAlpha(softsusy::ALPHAS);
   const double mz_pole  = qedqcd.displayPoleMZ();

   while (!converged) {
      calculate_DRbar_masses();

      double mz_run = mz_pole;

      if (get_thresholds() && threshold_corrections.mz > 0)
         mz_run = AbsSqrt(calculate_M2VZ_DRbar(Sqr(mz_pole)));

      double delta_alpha_em = 0.;
      double delta_alpha_s  = 0.;

      if (get_thresholds() && threshold_corrections.alpha_em > 0)
         delta_alpha_em = calculate_delta_alpha_em(alpha_em);

      if (get_thresholds() && threshold_corrections.alpha_s > 0)
         delta_alpha_s  = calculate_delta_alpha_s(alpha_s);

      const double alpha_em_drbar = alpha_em / (1.0 - delta_alpha_em);
      const double alpha_s_drbar  = alpha_s / (1.0 - delta_alpha_s);
      const double e_drbar        = Sqrt(4.0 * Pi * alpha_em_drbar);
      const double theta_w_drbar  = calculate_theta_w(qedqcd, alpha_em_drbar);

      v = Re((2 * mz_run) / Sqrt(0.6 * Sqr(g1) + Sqr(g2)));

      calculate_Yu_DRbar(qedqcd);
      calculate_Yd_DRbar(qedqcd);
      calculate_Ye_DRbar(qedqcd);
      calculate_Lambdax_DRbar();

      solve_ewsb();

      g1 = 1.2909944487358056 * e_drbar * Sec(theta_w_drbar);
      g2 = e_drbar * Csc(theta_w_drbar);
      g3 = 3.5449077018110318 * Sqrt(alpha_s_drbar);

      if (IsFinite(g1)) {
         problems.unflag_non_perturbative_parameter(standard_model_info::g1);
      } else {
         problems.flag_non_perturbative_parameter(
            standard_model_info::g1, g1, get_scale());
         g1 = Electroweak_constants::g1;
      }

      if (IsFinite(g2)) {
         problems.unflag_non_perturbative_parameter(standard_model_info::g2);
      } else {
         problems.flag_non_perturbative_parameter(
            standard_model_info::g2, g2, get_scale());
         g2 = Electroweak_constants::g2;
      }

      if (get_thresholds() && threshold_corrections.sin_theta_w > 0)
         qedqcd.setPoleMW(AbsSqrt(recalculate_m2w_pole(Sqr(qedqcd.displayPoleMW()))));

      converged = check_convergence(old);
      old = *this;
   }

   // run all SM parameters to the desired scale
   if (scale > 0.)
      run_to(scale);
}

void Standard_model::initial_guess_for_parameters(const softsusy::QedQcd& qedqcd)
{
   const double MH = input.get(Physical_input::mh_pole);
   const double mtpole = qedqcd.displayPoleMt();

   const double mu_guess = qedqcd.displayMass(softsusy::mUp);
   const double mc_guess = qedqcd.displayMass(softsusy::mCharm);
   const double mt_guess = get_thresholds() > 0 && threshold_corrections.mt > 0 ?
      qedqcd.displayMass(softsusy::mTop) - 30.0 :
      qedqcd.displayPoleMt();
   const double md_guess = qedqcd.displayMass(softsusy::mDown);
   const double ms_guess = qedqcd.displayMass(softsusy::mStrange);
   const double mb_guess = qedqcd.displayMass(softsusy::mBottom);
   const double me_guess = get_thresholds() > 0 ?
      qedqcd.displayMass(softsusy::mElectron) :
      qedqcd.displayPoleMel();
   const double mm_guess = get_thresholds() > 0 ?
      qedqcd.displayMass(softsusy::mMuon) :
      qedqcd.displayPoleMmuon();
   const double mtau_guess = qedqcd.displayMass(softsusy::mTau);

   // guess gauge couplings at mt
   const auto alpha_sm(qedqcd.guess_alpha_SM5(mtpole));

   g1 = Sqrt(4. * Pi * alpha_sm(0));
   g2 = Sqrt(4. * Pi * alpha_sm(1));
   g3 = Sqrt(4. * Pi * alpha_sm(2));

   set_scale(mtpole);

   v = Re(Electroweak_constants::vev);

   Eigen::Matrix<double,3,3> upQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   upQuarksDRbar(0,0) = mu_guess;
   upQuarksDRbar(1,1) = mc_guess;
   upQuarksDRbar(2,2) = mt_guess;
   Yu = -((1.4142135623730951*upQuarksDRbar) / v).transpose();

   Eigen::Matrix<double,3,3> downQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   downQuarksDRbar(0,0) = md_guess;
   downQuarksDRbar(1,1) = ms_guess;
   downQuarksDRbar(2,2) = mb_guess;
   Yd = ((1.4142135623730951*downQuarksDRbar)/v).transpose();

   Eigen::Matrix<double,3,3> downLeptonsDRbar(Eigen::Matrix<double,3,3>::Zero());
   downLeptonsDRbar(0,0) = me_guess;
   downLeptonsDRbar(1,1) = mm_guess;
   downLeptonsDRbar(2,2) = mtau_guess;
   Ye = ((1.4142135623730951*downLeptonsDRbar)/v).transpose();

   Lambdax = Sqr(MH) / Sqr(v);

   solve_ewsb_tree_level();
}

double Standard_model::calculate_delta_alpha_em(double alphaEm) const
{
   const double delta_alpha_em_SM = -0.28294212105225836*alphaEm*
      FiniteLog(Abs(MFu(2) / get_scale()));

   return delta_alpha_em_SM;
}

double Standard_model::calculate_delta_alpha_s(double alphaS) const
{
   const double delta_alpha_s_1loop = -0.1061032953945969*alphaS*
      FiniteLog(Abs(MFu(2) / get_scale()));

   double delta_alpha_s_2loop = 0.;
   double delta_alpha_s_3loop = 0.;

   if (get_thresholds() > 1 && threshold_corrections.alpha_s > 1) {
      sm_threeloop_as::Parameters pars;
      pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
      pars.mt   = MFu(2);
      pars.Q    = get_scale();

      const auto das_1L = sm_threeloop_as::delta_alpha_s_1loop_as(pars);
      const auto das_2L = sm_threeloop_as::delta_alpha_s_2loop_as_as(pars);

      delta_alpha_s_2loop = - das_2L + Sqr(das_1L);
   }

   if (get_thresholds() > 2 && get_threshold_corrections().alpha_s > 2) {
      sm_threeloop_as::Parameters pars;
      pars.as   = alphaS; // alpha_s(SM(5)) MS-bar
      pars.mt   = MFu(2);
      pars.Q    = get_scale();

      const auto das_1L = sm_threeloop_as::delta_alpha_s_1loop_as(pars);
      const auto das_2L = sm_threeloop_as::delta_alpha_s_2loop_as_as(pars);
      const auto das_3L = sm_threeloop_as::delta_alpha_s_3loop_as_as_as(pars);

      delta_alpha_s_3loop = - das_3L - Power3(das_1L) + 2. * das_1L * das_2L;
   }

   return delta_alpha_s_1loop + delta_alpha_s_2loop + delta_alpha_s_3loop;

}

double Standard_model::calculate_theta_w(const softsusy::QedQcd& qedqcd, double alpha_em_drbar)
{
   double theta_w = 0.;

   using namespace weinberg_angle;

   const double scale               = get_scale();
   const double mw_pole             = qedqcd.displayPoleMW();
   const double mz_pole             = qedqcd.displayPoleMZ();
   const double mt_pole             = qedqcd.displayPoleMt();
   const double mt_drbar            = MFu(2);
   const double mb_drbar            = MFd(2);
   const double mh_drbar            = AbsSqrt(M2hh);
   const double gY                  = g1 * 0.7745966692414834;
   const double pizztMZ             = Re(self_energy_VZ_1loop(Sqr(mz_pole)));
   const double piwwt0              = Re(self_energy_VWp_1loop(0.));
   const double self_energy_w_at_mw = Re(self_energy_VWp_1loop(Sqr(mw_pole)));

   Weinberg_angle::Self_energy_data se_data;
   se_data.scale    = scale;
   se_data.mt_pole  = mt_pole;
   se_data.mt_drbar = mt_drbar;
   se_data.mb_drbar = mb_drbar;
   se_data.gY       = gY;
   se_data.g2       = g2;

   double pizztMZ_corrected = pizztMZ;
   double piwwtMW_corrected = self_energy_w_at_mw;
   double piwwt0_corrected  = piwwt0;

   if (get_thresholds() > 1 && threshold_corrections.sin_theta_w > 1) {
      pizztMZ_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_z(pizztMZ, mz_pole,
            se_data);
      piwwtMW_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(
            self_energy_w_at_mw, mw_pole, se_data);
      piwwt0_corrected =
         Weinberg_angle::replace_mtop_in_self_energy_w(piwwt0, 0.,
            se_data);
   }

   Weinberg_angle::Data data;
   data.scale               = scale;
   data.alpha_em_drbar      = alpha_em_drbar;
   data.fermi_contant       = qedqcd.displayFermiConstant();
   data.self_energy_z_at_mz = pizztMZ_corrected;
   data.self_energy_w_at_mw = piwwtMW_corrected;
   data.self_energy_w_at_0  = piwwt0_corrected;
   data.mw_pole             = mw_pole;
   data.mz_pole             = mz_pole;
   data.mt_pole             = mt_pole;
   data.mh_drbar            = mh_drbar;
   data.gY                  = gY;
   data.g2                  = g2;
   data.g3                  = g3;

   Weinberg_angle weinberg;
   weinberg.disable_susy_contributions();
   weinberg.set_number_of_loops(threshold_corrections.sin_theta_w);
   weinberg.set_data(data);

   const int error = weinberg.calculate();

   theta_w = ArcSin(weinberg.get_sin_theta());

   if (error)
      problems.flag_no_sinThetaW_convergence();
   else
      problems.unflag_no_sinThetaW_convergence();

   return theta_w;
}

void Standard_model::calculate_Yu_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> upQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   upQuarksDRbar(0,0)      = qedqcd.displayMass(softsusy::mUp);
   upQuarksDRbar(1,1)      = qedqcd.displayMass(softsusy::mCharm);
   upQuarksDRbar(2,2)      = qedqcd.displayPoleMt();

   if (get_thresholds() && threshold_corrections.mt > 0)
      upQuarksDRbar(2,2) = calculate_MFu_DRbar(qedqcd.displayPoleMt(), 2);

   Yu = ((1.4142135623730951*upQuarksDRbar)/v).transpose();
}

void Standard_model::calculate_Yd_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> downQuarksDRbar(Eigen::Matrix<double,3,3>::Zero());
   downQuarksDRbar(0,0)   = qedqcd.displayMass(softsusy::mDown);
   downQuarksDRbar(1,1)   = qedqcd.displayMass(softsusy::mStrange);
   downQuarksDRbar(2,2)   = qedqcd.displayMass(softsusy::mBottom);

   if (get_thresholds() && threshold_corrections.mb > 0)
      downQuarksDRbar(2,2) = calculate_MFd_DRbar(qedqcd.displayMass(softsusy::mBottom), 2);

   Yd = ((1.4142135623730951*downQuarksDRbar)/v).transpose();
}

void Standard_model::calculate_Ye_DRbar(const softsusy::QedQcd& qedqcd)
{
   Eigen::Matrix<double,3,3> downLeptonsDRbar(Eigen::Matrix<double,3,3>::Zero());
   downLeptonsDRbar(0,0) = qedqcd.displayPoleMel();
   downLeptonsDRbar(1,1) = qedqcd.displayPoleMmuon();
   downLeptonsDRbar(2,2) = qedqcd.displayPoleMtau();

   if (get_thresholds()) {
      downLeptonsDRbar(0,0) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mElectron), 0);
      downLeptonsDRbar(1,1) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mMuon), 1);
   }

   if (get_thresholds() && threshold_corrections.mtau > 0)
      downLeptonsDRbar(2,2) = calculate_MFe_DRbar(qedqcd.displayMass(softsusy::mTau), 2);

   Ye = ((1.4142135623730951*downLeptonsDRbar)/v).transpose();
}

void Standard_model::calculate_Lambdax_DRbar()
{
   const double mh_pole_sqr = Sqr(input.get(Physical_input::mh_pole));
   double mh_msbar_sqr = mh_pole_sqr;

   if (get_thresholds() && threshold_corrections.mh > 0)
      mh_msbar_sqr = calculate_M2hh_DRbar(mh_pole_sqr);

   Lambdax = mh_msbar_sqr / Sqr(v);
}

double Standard_model::recalculate_m2w_pole(double p2)
{
   calculate_M2VWp();

   const double mw_pole_sqr = M2VWp - Re(self_energy_VWp_1loop(p2));

   if (mw_pole_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   return mw_pole_sqr;
}

double Standard_model::max_rel_diff(const Standard_model& old) const
{
   std::array<double, 12> diff{};

   diff[0] = MaxRelDiff(old.M2hh, M2hh);
   diff[1] = MaxRelDiff(old.M2VZ, M2VZ);
   for (int i = 0; i < 3; ++i)
      diff[i + 2] = MaxRelDiff(old.MFd(i), MFd(i));
   for (int i = 0; i < 3; ++i)
      diff[i + 5] = MaxRelDiff(old.MFu(i), MFu(i));
   for (int i = 0; i < 3; ++i)
      diff[i + 8] = MaxRelDiff(old.MFe(i), MFe(i));
   diff[11] = MaxRelDiff(old.M2VWp,M2VWp);

   return *std::max_element(diff.cbegin(), diff.cend());
}

bool Standard_model::check_convergence(const Standard_model& old) const
{
   return max_rel_diff(old) < precision;
}

Eigen::ArrayXd Standard_model::get() const
{
   Eigen::ArrayXd pars(numberOfParameters);

   pars(0) = g1;
   pars(1) = g2;
   pars(2) = g3;
   pars(3) = Lambdax;
   pars(4) = Yu(0,0);
   pars(5) = Yu(0,1);
   pars(6) = Yu(0,2);
   pars(7) = Yu(1,0);
   pars(8) = Yu(1,1);
   pars(9) = Yu(1,2);
   pars(10) = Yu(2,0);
   pars(11) = Yu(2,1);
   pars(12) = Yu(2,2);
   pars(13) = Yd(0,0);
   pars(14) = Yd(0,1);
   pars(15) = Yd(0,2);
   pars(16) = Yd(1,0);
   pars(17) = Yd(1,1);
   pars(18) = Yd(1,2);
   pars(19) = Yd(2,0);
   pars(20) = Yd(2,1);
   pars(21) = Yd(2,2);
   pars(22) = Ye(0,0);
   pars(23) = Ye(0,1);
   pars(24) = Ye(0,2);
   pars(25) = Ye(1,0);
   pars(26) = Ye(1,1);
   pars(27) = Ye(1,2);
   pars(28) = Ye(2,0);
   pars(29) = Ye(2,1);
   pars(30) = Ye(2,2);
   pars(31) = mu2;
   pars(32) = v;

   return pars;
}

void Standard_model::set(const Eigen::ArrayXd& pars)
{
   g1 = pars(0);
   g2 = pars(1);
   g3 = pars(2);
   Lambdax = pars(3);
   Yu(0,0) = pars(4);
   Yu(0,1) = pars(5);
   Yu(0,2) = pars(6);
   Yu(1,0) = pars(7);
   Yu(1,1) = pars(8);
   Yu(1,2) = pars(9);
   Yu(2,0) = pars(10);
   Yu(2,1) = pars(11);
   Yu(2,2) = pars(12);
   Yd(0,0) = pars(13);
   Yd(0,1) = pars(14);
   Yd(0,2) = pars(15);
   Yd(1,0) = pars(16);
   Yd(1,1) = pars(17);
   Yd(1,2) = pars(18);
   Yd(2,0) = pars(19);
   Yd(2,1) = pars(20);
   Yd(2,2) = pars(21);
   Ye(0,0) = pars(22);
   Ye(0,1) = pars(23);
   Ye(0,2) = pars(24);
   Ye(1,0) = pars(25);
   Ye(1,1) = pars(26);
   Ye(1,2) = pars(27);
   Ye(2,0) = pars(28);
   Ye(2,1) = pars(29);
   Ye(2,2) = pars(30);
   mu2 = pars(31);
   v = pars(32);
}

void Standard_model::run_to(double scale, double eps)
{
   if (eps < 0.0)
      eps = precision;
   Beta_function::run_to(scale, eps);
}

Eigen::ArrayXd Standard_model::beta() const
{
   return calc_beta().get().unaryExpr(Chop<double>(get_zero_threshold()));
}

Standard_model Standard_model::calc_beta(int loops) const
{
   double beta_g1 = 0.;
   double beta_g2 = 0.;
   double beta_g3 = 0.;
   double beta_Lambdax = 0.;
   Eigen::Matrix<double,3,3> beta_Yu = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Yd = Eigen::Matrix<double,3,3>::Zero();
   Eigen::Matrix<double,3,3> beta_Ye = Eigen::Matrix<double,3,3>::Zero();
   double beta_mu2 = 0.;
   double beta_v = 0.;

   if (loops > 0) {
      Beta_traces traces;
      calc_beta_traces(traces);

      beta_g1 += calc_beta_g1_one_loop(traces);
      beta_g2 += calc_beta_g2_one_loop(traces);
      beta_g3 += calc_beta_g3_one_loop(traces);
      beta_Lambdax += calc_beta_Lambdax_one_loop(traces);
      beta_Yu += calc_beta_Yu_one_loop(traces);
      beta_Yd += calc_beta_Yd_one_loop(traces);
      beta_Ye += calc_beta_Ye_one_loop(traces);
      beta_mu2 += calc_beta_mu2_one_loop(traces);
      beta_v += calc_beta_v_one_loop(traces);

      if (loops > 1) {
         beta_g1 += calc_beta_g1_two_loop(traces);
         beta_g2 += calc_beta_g2_two_loop(traces);
         beta_g3 += calc_beta_g3_two_loop(traces);
         beta_Lambdax += calc_beta_Lambdax_two_loop(traces);
         beta_Yu += calc_beta_Yu_two_loop(traces);
         beta_Yd += calc_beta_Yd_two_loop(traces);
         beta_Ye += calc_beta_Ye_two_loop(traces);
         beta_mu2 += calc_beta_mu2_two_loop(traces);
         beta_v += calc_beta_v_two_loop(traces);

         if (loops > 2) {
            beta_g1 += calc_beta_g1_three_loop(traces);
            beta_g2 += calc_beta_g2_three_loop(traces);
            beta_g3 += calc_beta_g3_three_loop(traces);
            beta_Lambdax += calc_beta_Lambdax_three_loop(traces);
            beta_Yu += calc_beta_Yu_three_loop(traces);
            beta_Yd += calc_beta_Yd_three_loop(traces);
            beta_Ye += calc_beta_Ye_three_loop(traces);
            beta_mu2 += calc_beta_mu2_three_loop(traces);

            if (loops > 3) {
               beta_g3 += calc_beta_g3_four_loop(traces);
               beta_Lambdax += calc_beta_Lambdax_four_loop(traces);
               beta_Yu += calc_beta_Yu_four_loop(traces);
            }
         }
      }
   }

   return Standard_model(get_scale(), loops, get_thresholds(),
                         beta_g1, beta_g2, beta_g3, beta_Lambdax, beta_Yu, beta_Yd, beta_Ye,
                         beta_mu2, beta_v);
}

Standard_model Standard_model::calc_beta() const
{
   return calc_beta(get_loops());
}

void Standard_model::calc_beta_traces(Beta_traces& traces) const
{
   const int loops = get_loops();

   if (loops > 0) {
      traces.traceYdAdjYd = Re((Yd*Yd.adjoint()).trace());
      traces.traceYeAdjYe = Re((Ye*Ye.adjoint()).trace());
      traces.traceYuAdjYu = Re((Yu*Yu.adjoint()).trace());
      traces.traceYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint())
         .trace());
      traces.traceYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint())
         .trace());
      traces.traceYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint())
         .trace());
   }

   if (loops > 1) {
      traces.traceYdAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint())
         .trace());
      traces.traceYdAdjYdYdAdjYdYdAdjYd = Re((Yd*Yd.adjoint()*Yd*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      traces.traceYdAdjYdYdAdjYuYuAdjYd = Re((Yd*Yd.adjoint()*Yd*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      traces.traceYdAdjYuYuAdjYdYdAdjYd = Re((Yd*Yu.adjoint()*Yu*Yd.adjoint(
         )*Yd*Yd.adjoint()).trace());
      traces.traceYdAdjYuYuAdjYuYuAdjYd = Re((Yd*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yd.adjoint()).trace());
      traces.traceYeAdjYeYeAdjYeYeAdjYe = Re((Ye*Ye.adjoint()*Ye*Ye.adjoint(
         )*Ye*Ye.adjoint()).trace());
      traces.traceYuAdjYuYuAdjYuYuAdjYu = Re((Yu*Yu.adjoint()*Yu*Yu.adjoint(
         )*Yu*Yu.adjoint()).trace());
   }
}


double Standard_model::calc_beta_g1_one_loop(const Beta_traces&) const
{
   double beta_g1;

   beta_g1 = Re(4.1*oneOver16PiSqr*Cube(g1));

   return beta_g1;
}

double Standard_model::calc_beta_g1_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g1;

   beta_g1 = Re(0.02*twoLoop*Cube(g1)*(199*Sqr(g1) + 5*(-5*traceYdAdjYd -
      15*traceYeAdjYe - 17*traceYuAdjYu + 27*Sqr(g2) + 88*Sqr(g3))));

   return beta_g1;
}

double Standard_model::calc_beta_g1_three_loop(const Beta_traces&) const
{
   double beta_g1;

   beta_g1 = Re(-0.000041666666666666665*threeLoop*Cube(g1)*(388613*Quad(
      g1) - 10*Sqr(g1)*(648*Lambdax + 1845*Sqr(g2) - 4384*Sqr(g3) - 8481*Sqr(Yu
      (2,2))) - 75*(3945*Quad(g2) - 6*Sqr(g2)*(-24*Lambdax + 32*Sqr(g3) + 785*
      Sqr(Yu(2,2))) + 4*(4752*Quad(g3) + 945*Quad(Yu(2,2)) - 36*Sqr(Lambdax) -
      464*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g1;
}

double Standard_model::calc_beta_g2_one_loop(const Beta_traces&) const
{
   double beta_g2;

   beta_g2 = Re(-3.1666666666666665*oneOver16PiSqr*Cube(g2));

   return beta_g2;
}

double Standard_model::calc_beta_g2_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g2;

   beta_g2 = Re(0.03333333333333333*twoLoop*Cube(g2)*(27*Sqr(g1) + 5*(-3*
      (3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu) + 35*Sqr(g2) + 72*Sqr(g3
      ))));

   return beta_g2;
}

double Standard_model::calc_beta_g2_three_loop(const Beta_traces&) const
{
   double beta_g2;

   beta_g2 = Re(0.000023148148148148147*threeLoop*Cube(g2)*(-151119*Quad(
      g1) + 270*Sqr(g1)*(24*Lambdax + 873*Sqr(g2) - 32*Sqr(g3) - 593*Sqr(Yu(2,2
      ))) + 25*(324953*Quad(g2) + 162*Sqr(g2)*(8*Lambdax + 416*Sqr(g3) - 243*
      Sqr(Yu(2,2))) + 108*(1296*Quad(g3) + 147*Quad(Yu(2,2)) - 12*Sqr(Lambdax)
      - 112*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g2;
}

double Standard_model::calc_beta_g3_one_loop(const Beta_traces&) const
{
   double beta_g3;

   beta_g3 = Re(-7*oneOver16PiSqr*Cube(g3));

   return beta_g3;
}

double Standard_model::calc_beta_g3_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_g3;

   beta_g3 = Re(-0.1*twoLoop*Cube(g3)*(-11*Sqr(g1) + 5*(-9*Sqr(g2) + 4*(
      traceYdAdjYd + traceYuAdjYu + 13*Sqr(g3)))));

   return beta_g3;
}

double Standard_model::calc_beta_g3_three_loop(const Beta_traces&) const
{
   double beta_g3;

   beta_g3 = Re(0.008333333333333333*threeLoop*Cube(g3)*(-523*Quad(g1) +
      Sqr(g1)*(-9*Sqr(g2) + 616*Sqr(g3) - 303*Sqr(Yu(2,2))) + 15*(109*Quad(g2)
      + 3*Sqr(g2)*(56*Sqr(g3) - 31*Sqr(Yu(2,2))) + 20*(13*Quad(g3) + 6*Quad(Yu(
      2,2)) - 16*Sqr(g3)*Sqr(Yu(2,2))))));

   return beta_g3;
}

double Standard_model::calc_beta_g3_four_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   double beta_g3;

   beta_g3 = Re(-2472.2837425797156*Cube(g3)*Quad(oneOver16PiSqr)*(1.*
      Power6(g3) + 0.04569149546770327*Power6(Yu(2,2)) + 0.0060672647486441755*
      Lambdax*Quad(Yu(2,2)) - 0.12603833934188147*Quad(Yu(2,2))*Sqr(g3) +
      0.16927060578749137*Quad(g3)*Sqr(Yu(2,2)) - 0.003640358849186505*Sqr(
      Lambdax)*Sqr(Yu(2,2))));

   return beta_g3;
}

double Standard_model::calc_beta_Lambdax_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_Lambdax;

   beta_Lambdax = Re(oneOver16PiSqr*(0.27*Quad(g1) + 2.25*Quad(g2) - 9*
      Lambdax*Sqr(g2) + 0.9*Sqr(g1)*(-2*Lambdax + Sqr(g2)) + 4*(-3*
      traceYdAdjYdYdAdjYd - traceYeAdjYeYeAdjYe - 3*traceYuAdjYuYuAdjYu + 3*
      traceYdAdjYd*Lambdax + traceYeAdjYe*Lambdax + 3*traceYuAdjYu*Lambdax + 3*
      Sqr(Lambdax))));

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;
   const double traceYdAdjYdYdAdjYdYdAdjYd =
      traces.traceYdAdjYdYdAdjYdYdAdjYd;
   const double traceYdAdjYdYdAdjYuYuAdjYd =
      traces.traceYdAdjYdYdAdjYuYuAdjYd;
   const double traceYdAdjYuYuAdjYdYdAdjYd =
      traces.traceYdAdjYuYuAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYuYuAdjYd =
      traces.traceYdAdjYuYuAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYeYeAdjYe =
      traces.traceYeAdjYeYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYuYuAdjYu =
      traces.traceYuAdjYuYuAdjYuYuAdjYu;

   double beta_Lambdax;

   beta_Lambdax = Re(twoLoop*(60*traceYdAdjYdYdAdjYdYdAdjYd - 24*
      traceYdAdjYdYdAdjYuYuAdjYd + 12*traceYdAdjYuYuAdjYdYdAdjYd - 12*
      traceYdAdjYuYuAdjYuYuAdjYd + 20*traceYeAdjYeYeAdjYeYeAdjYe + 60*
      traceYuAdjYuYuAdjYuYuAdjYu - 78*Cube(Lambdax) - 3*traceYdAdjYdYdAdjYd*
      Lambdax - 42*traceYdAdjYuYuAdjYd*Lambdax - traceYeAdjYeYeAdjYe*Lambdax -
      3*traceYuAdjYuYuAdjYu*Lambdax - 3.411*Power6(g1) + 38.125*Power6(g2) -
      0.125*(36*traceYdAdjYd + 12*traceYeAdjYe + 36*traceYuAdjYu + 73*Lambdax)*
      Quad(g2) + 1.5*Lambdax*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*
      traceYuAdjYu + 36*Lambdax)*Sqr(g2) - 0.015*Quad(g1)*(-60*traceYdAdjYd +
      300*traceYeAdjYe + 228*traceYuAdjYu - 629*Lambdax + 559*Sqr(g2)) - 64*
      traceYdAdjYdYdAdjYd*Sqr(g3) - 64*traceYuAdjYuYuAdjYu*Sqr(g3) + 80*
      traceYdAdjYd*Lambdax*Sqr(g3) + 80*traceYuAdjYu*Lambdax*Sqr(g3) - 72*
      traceYdAdjYd*Sqr(Lambdax) - 24*traceYeAdjYe*Sqr(Lambdax) - 72*
      traceYuAdjYu*Sqr(Lambdax) - 0.025*Sqr(g1)*(289*Quad(g2) - 6*(36*
      traceYdAdjYd + 44*traceYeAdjYe + 84*traceYuAdjYu + 39*Lambdax)*Sqr(g2) -
      4*(16*traceYdAdjYdYdAdjYd - 48*traceYeAdjYeYeAdjYe - 32*
      traceYuAdjYuYuAdjYu + 25*traceYdAdjYd*Lambdax + 75*traceYeAdjYe*Lambdax +
      85*traceYuAdjYu*Lambdax + 108*Sqr(Lambdax)))));

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_three_loop(const Beta_traces&) const
{
   double beta_Lambdax;

   beta_Lambdax = Re(0.0001*threeLoop*(-60320*Power8(g1) - 4563640*Power8
      (g2) - 40*Power6(g1)*(-14084*Lambdax + 1543*Sqr(g2) - 663*Sqr(g3) - 11117
      *Sqr(Yu(2,2))) + 20*Power6(g2)*(865483*Lambdax + 15072*Sqr(g3) + 125000*
      Sqr(Yu(2,2))) + 2*Sqr(g2)*(-968630*Cube(Lambdax) + 1482760*Power6(Yu(2,2)
      ) - 54700*Lambdax*Quad(Yu(2,2)) + 266980*Quad(Yu(2,2))*Sqr(g3) + 151443*
      Lambdax*Sqr(g3)*Sqr(Yu(2,2)) - 1797695*Sqr(Lambdax)*Sqr(Yu(2,2))) + 80*
      Quad(g2)*(7942*Quad(Yu(2,2)) - 98785*Sqr(Lambdax) - 79916*Lambdax*Sqr(Yu(
      2,2)) + Sqr(g3)*(-14286*Lambdax + 8232*Sqr(Yu(2,2)))) + 2*Quad(g1)*(
      130000*Quad(g2) + 318960*Quad(Yu(2,2)) - 927660*Sqr(Lambdax) - 748599*
      Lambdax*Sqr(Yu(2,2)) + Sqr(g3)*(-83810*Lambdax + 20320*Sqr(Yu(2,2))) + 10
      *Sqr(g2)*(61753*Lambdax + 2210*Sqr(g3) + 21254*Sqr(Yu(2,2)))) - 10*Sqr(g1
      )*(38745*Cube(Lambdax) + 151556*Power6(g2) - 135720*Power6(Yu(2,2)) +
      42030*Lambdax*Quad(Yu(2,2)) + 63869*Sqr(Lambdax)*Sqr(Yu(2,2)) - 4*Quad(g2
      )*(39819*Lambdax + 1507*Sqr(g3) + 13041*Sqr(Yu(2,2))) - 4*Sqr(g3)*(17570*
      Quad(Yu(2,2)) + 8727*Lambdax*Sqr(Yu(2,2))) + 2*Sqr(g2)*(140712*Quad(Yu(2,
      2)) + 158320*Sqr(Lambdax) - 5615*Lambdax*Sqr(Yu(2,2)) - 22772*Sqr(g3)*Sqr
      (Yu(2,2)))) - 5*(893528*Lambdax*Power6(Yu(2,2)) + 1945192*Power8(Yu(2,2))
      - 3005675*Quad(Lambdax) - 3536520*Quad(Yu(2,2))*Sqr(Lambdax) - 873000*
      Cube(Lambdax)*Sqr(Yu(2,2)) + 8*Quad(g3)*(50201*Quad(Yu(2,2)) - 178484*
      Lambdax*Sqr(Yu(2,2))) - 4*Sqr(g3)*(500988*Power6(Yu(2,2)) - 662866*
      Lambdax*Quad(Yu(2,2)) + 80385*Sqr(Lambdax)*Sqr(Yu(2,2))))));

   return beta_Lambdax;
}

double Standard_model::calc_beta_Lambdax_four_loop(const Beta_traces&) const
{
   double beta_Lambdax;

   beta_Lambdax = Re(16616.34*Power6(g3)*Quad(oneOver16PiSqr)*Quad(Yu(2,2
      )));

   return beta_Lambdax;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (oneOver16PiSqr*(Yu*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 0.85*Sqr(g1) - 2.25*Sqr(g2) - 8*Sqr(g3)) - 1.5*(Yu*
      Yd.adjoint()*Yd) + 1.5*(Yu*Yu.adjoint()*Yu))).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (twoLoop*(0.0016666666666666668*Yu*(1187*Quad(g1) + 5*Sqr(g1
      )*(75*traceYdAdjYd + 225*traceYeAdjYe + 255*traceYuAdjYu - 54*Sqr(g2) +
      152*Sqr(g3)) - 75*(46*Quad(g2) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Quad(g3) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(-43*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yu*Yd.adjoint
      ()*Yd) + 0.0125*(223*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3)))*(Yu*Yu.adjoint(
      )*Yu) + 2.75*(Yu*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - 0.25*(Yu*Yd.adjoint()
      *Yd*Yu.adjoint()*Yu) - Yu*Yu.adjoint()*Yu*Yd.adjoint()*Yd + 1.5*(Yu*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_three_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (0.00005*PROJECTOR*threeLoop*(321980*Power6(g1) + 3396580*
      Power6(g2) - 5*Quad(g2)*(21375*Lambdax + 84288*Sqr(g3) - 67960*Sqr(Yu(2,2
      ))) - 5*Quad(g1)*(5445*Lambdax + 17768*Sqr(g2) + 89276*Sqr(g3) + 97688*
      Sqr(Yu(2,2))) - 10*Sqr(g1)*(9486*Quad(g2) + 30192*Quad(g3) + 60925*Quad(
      Yu(2,2)) - 4500*Sqr(Lambdax) + Sqr(g2)*(-2925*Lambdax + 32100*Sqr(g3) -
      69658*Sqr(Yu(2,2))) + 12700*Lambdax*Sqr(Yu(2,2)) - 36148*Sqr(g3)*Sqr(Yu(2
      ,2))) + 10*Sqr(g2)*(147308*Quad(g3) + 96740*Sqr(g3)*Sqr(Yu(2,2)) + 1125*(
      -177*Quad(Yu(2,2)) + 20*Sqr(Lambdax) - 60*Lambdax*Sqr(Yu(2,2)))) + 2*(
      -45000*Cube(Lambdax) - 6193500*Power6(g3) + 586028*Power6(Yu(2,2)) +
      990000*Lambdax*Quad(Yu(2,2)) + 3637640*Quad(g3)*Sqr(Yu(2,2)) + 9375*Sqr(
      Lambdax)*Sqr(Yu(2,2)) + 10000*Sqr(g3)*(-157*Quad(Yu(2,2)) + 8*Lambdax*Sqr
      (Yu(2,2)))))*Yu(2,2)).real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yu_four_loop(const Beta_traces&) const
{
   DEFINE_PROJECTOR(3,3,3,3)

   Eigen::Matrix<double,3,3> beta_Yu;

   beta_Yu = (2308.18*PROJECTOR*Power8(g3)*Quad(oneOver16PiSqr)*Yu(2,2))
      .real();

   return beta_Yu;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (oneOver16PiSqr*(-0.25*Yd*(Sqr(g1) + 9*Sqr(g2) + 4*(-3*
      traceYdAdjYd - traceYeAdjYe - 3*traceYuAdjYu + 8*Sqr(g3))) + 1.5*(Yd*
      Yd.adjoint()*Yd) - 1.5*(Yd*Yu.adjoint()*Yu))).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (twoLoop*(-0.0016666666666666668*Yd*(127*Quad(g1) + 5*Sqr(g1
      )*(-15*(5*traceYdAdjYd + 15*traceYeAdjYe + 17*traceYuAdjYu) + 162*Sqr(g2)
      - 248*Sqr(g3)) + 75*(46*Quad(g2) - 3*Sqr(g2)*(5*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 24*Sqr(g3)) + 2*(432*Quad(g3) - 80*(
      traceYdAdjYd + traceYuAdjYu)*Sqr(g3) + 3*(9*traceYdAdjYdYdAdjYd - 2*
      traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2*
      Sqr(Lambdax))))) + 0.0125*(187*Sqr(g1) + 675*Sqr(g2) + 20*(-3*(9*
      traceYdAdjYd + 3*traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 64*Sqr(g3))
      )*(Yd*Yd.adjoint()*Yd) + 0.0125*(-79*Sqr(g1) + 5*(20*(3*traceYdAdjYd +
      traceYeAdjYe + 3*traceYuAdjYu) + 9*Sqr(g2) - 256*Sqr(g3)))*(Yd*Yu.adjoint
      ()*Yu) + 1.5*(Yd*Yd.adjoint()*Yd*Yd.adjoint()*Yd) - Yd*Yd.adjoint()*Yd*
      Yu.adjoint()*Yu - 0.25*(Yd*Yu.adjoint()*Yu*Yd.adjoint()*Yd) + 2.75*(Yd*
      Yu.adjoint()*Yu*Yu.adjoint()*Yu))).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Yd_three_loop(const Beta_traces&) const
{
   Eigen::Matrix<double,3,3> beta_Yd;

   beta_Yd = (ZEROMATRIX(3,3)).real();

   return beta_Yd;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (oneOver16PiSqr*(Ye*(3*traceYdAdjYd + traceYeAdjYe + 3*
      traceYuAdjYu - 2.25*Sqr(g1) - 2.25*Sqr(g2)) + 1.5*(Ye*Ye.adjoint()*Ye)))
      .real();

   return beta_Ye;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (twoLoop*(0.005*Ye*(1371*Quad(g1) + 5*Sqr(g1)*(25*
      traceYdAdjYd + 75*traceYeAdjYe + 85*traceYuAdjYu + 54*Sqr(g2)) - 25*(46*
      Quad(g2) - 15*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(g2) -
      2*(80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) - 3*(9*traceYdAdjYdYdAdjYd -
      2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu - 2
      *Sqr(Lambdax))))) + 0.0375*(129*Sqr(g1) + 5*(-4*(9*traceYdAdjYd + 3*
      traceYeAdjYe + 9*traceYuAdjYu + 8*Lambdax) + 45*Sqr(g2)))*(Ye*Ye.adjoint(
      )*Ye) + 1.5*(Ye*Ye.adjoint()*Ye*Ye.adjoint()*Ye))).real();

   return beta_Ye;
}

Eigen::Matrix<double,3,3> Standard_model::calc_beta_Ye_three_loop(const Beta_traces&) const
{
   Eigen::Matrix<double,3,3> beta_Ye;

   beta_Ye = (ZEROMATRIX(3,3)).real();

   return beta_Ye;
}

double Standard_model::calc_beta_mu2_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_mu2;

   beta_mu2 = Re(oneOver16PiSqr*(2*mu2*(3*traceYdAdjYd + traceYeAdjYe + 3
      *(traceYuAdjYu + Lambdax)) - 0.9*mu2*Sqr(g1) - 4.5*mu2*Sqr(g2)));

   return beta_mu2;
}

double Standard_model::calc_beta_mu2_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_mu2;

   beta_mu2 = Re(0.0025*mu2*twoLoop*(1671*Quad(g1) + 10*Sqr(g1)*(50*
      traceYdAdjYd + 150*traceYeAdjYe + 170*traceYuAdjYu + 288*Lambdax + 45*Sqr
      (g2)) - 25*(145*Quad(g2) - 12*(15*traceYdAdjYd + 5*traceYeAdjYe + 15*
      traceYuAdjYu + 48*Lambdax)*Sqr(g2) - 8*(80*(traceYdAdjYd + traceYuAdjYu)*
      Sqr(g3) - 3*(9*traceYdAdjYdYdAdjYd + 14*traceYdAdjYuYuAdjYd + 3*
      traceYeAdjYeYeAdjYe + 9*traceYuAdjYuYuAdjYu + 24*traceYdAdjYd*Lambdax + 8
      *traceYeAdjYe*Lambdax + 24*traceYuAdjYu*Lambdax + 10*Sqr(Lambdax))))));

   return beta_mu2;
}

double Standard_model::calc_beta_mu2_three_loop(const Beta_traces&) const
{
   double beta_mu2;

   beta_mu2 = Re(2*mu2*threeLoop*(128.25*Cube(Lambdax) +
      8.378314604562993*Power6(g1) + 301.7235994495886*Power6(g2) +
      154.40506064218175*Power6(Yd(2,2)) + 3.4683535473939138*Power6(Ye(2,2)) +
      154.40506064218175*Power6(Yu(2,2)) + 173.69714554123618*Lambdax*Quad(Yd(
      2,2)) + 75.89904851374538*Lambdax*Quad(Ye(2,2)) + 173.69714554123618*
      Lambdax*Quad(Yu(2,2)) - 209.24048513745396*Quad(Yd(2,2))*Sqr(g3) -
      209.24048513745396*Quad(Yu(2,2))*Sqr(g3) + 178.48396765750311*Quad(g3)*
      Sqr(Yd(2,2)) + 72.*Quad(Ye(2,2))*Sqr(Yd(2,2)) + 296.2115485137454*Quad(Yu
      (2,2))*Sqr(Yd(2,2)) + 40.19238810996319*Lambdax*Sqr(g3)*Sqr(Yd(2,2)) +
      37.125*Sqr(Lambdax)*Sqr(Yd(2,2)) + 72.*Quad(Yd(2,2))*Sqr(Ye(2,2)) + 72.*
      Quad(Yu(2,2))*Sqr(Ye(2,2)) + 12.375*Sqr(Lambdax)*Sqr(Ye(2,2)) - 54.*
      Lambdax*Sqr(Yd(2,2))*Sqr(Ye(2,2)) + Quad(g2)*(-32.25723415592714*Lambdax
      - 28.572145541236182*Sqr(g3) - 102.62651936694535*Sqr(Yd(2,2)) -
      34.208839788981784*Sqr(Ye(2,2)) - 102.62651936694535*Sqr(Yu(2,2))) + Quad
      (g1)*(-32.90278936623708*Lambdax + 9.778005607297105*Sqr(g2) -
      4.190581346047974*Sqr(g3) - 11.837888914984243*Sqr(Yd(2,2)) -
      20.542464319265452*Sqr(Ye(2,2)) - 27.721423523790055*Sqr(Yu(2,2))) +
      178.48396765750311*Quad(g3)*Sqr(Yu(2,2)) + 296.2115485137454*Quad(Yd(2,2)
      )*Sqr(Yu(2,2)) + 72.*Quad(Ye(2,2))*Sqr(Yu(2,2)) + 40.19238810996319*
      Lambdax*Sqr(g3)*Sqr(Yu(2,2)) + 37.125*Sqr(Lambdax)*Sqr(Yu(2,2)) -
      208.57214554123618*Lambdax*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 16.698731351660527
      *Sqr(g3)*Sqr(Yd(2,2))*Sqr(Yu(2,2)) - 54.*Lambdax*Sqr(Ye(2,2))*Sqr(Yu(2,2)
      ) + 10.5*Sqr(Yd(2,2))*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g2)*(
      -3.829281688145728*Quad(Yd(2,2)) + 3.223572770618091*Quad(Ye(2,2)) -
      3.829281688145728*Quad(Yu(2,2)) - 48.205536385309045*Sqr(Lambdax) -
      53.09857277061809*Lambdax*Sqr(Ye(2,2)) + Sqr(Yd(2,2))*(
      -159.29571831185424*Lambdax - 13.500000000000004*Sqr(Ye(2,2)) -
      62.83053638530905*Sqr(Yu(2,2))) - 159.29571831185424*Lambdax*Sqr(Yu(2,2))
      - 13.500000000000004*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g3)*(
      7.572145541236182*Sqr(Yd(2,2)) + 7.572145541236182*Sqr(Yu(2,2)))) + Sqr(
      g1)*(9.930778506201856*Quad(g2) + 11.916358216494471*Quad(Yd(2,2)) -
      15.051929108247235*Quad(Ye(2,2)) - 7.507690297250921*Quad(Yu(2,2)) -
      9.641107277061808*Sqr(Lambdax) - 32.86395336511993*Lambdax*Sqr(Yd(2,2)) -
      7.605285445876382*Lambdax*Sqr(Ye(2,2)) - 2.7*Sqr(Yd(2,2))*Sqr(Ye(2,2)) -
      29.849524256872694*Lambdax*Sqr(Yu(2,2)) - 23.946234141895758*Sqr(Yd(2,2)
      )*Sqr(Yu(2,2)) - 2.7*Sqr(Ye(2,2))*Sqr(Yu(2,2)) + Sqr(g3)*(
      -2.091983828751535*Sqr(Yd(2,2)) + 8.727254982244787*Sqr(Yu(2,2))) + Sqr(
      g2)*(-18.911535445876382*Lambdax - 8.109375*Sqr(Yd(2,2)) +
      10.602411385309043*Sqr(Ye(2,2)) + 11.470322300901756*Sqr(Yu(2,2))))));

   return beta_mu2;
}

double Standard_model::calc_beta_v_one_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;

   double beta_v;

   beta_v = Re(oneOver16PiSqr*(0.6*v*Sqr(g1) + v*(-3*traceYdAdjYd -
      traceYeAdjYe - 3*traceYuAdjYu + 3*Sqr(g2))));

   return beta_v;
}

double Standard_model::calc_beta_v_two_loop(const Beta_traces& traces) const
{
   const double traceYdAdjYd = traces.traceYdAdjYd;
   const double traceYeAdjYe = traces.traceYeAdjYe;
   const double traceYuAdjYu = traces.traceYuAdjYu;
   const double traceYdAdjYdYdAdjYd = traces.traceYdAdjYdYdAdjYd;
   const double traceYdAdjYuYuAdjYd = traces.traceYdAdjYuYuAdjYd;
   const double traceYeAdjYeYeAdjYe = traces.traceYeAdjYeYeAdjYe;
   const double traceYuAdjYuYuAdjYu = traces.traceYuAdjYuYuAdjYu;

   double beta_v;

   beta_v = Re(-0.00125*twoLoop*v*(1221*Quad(g1) + 10*Sqr(g1)*(122*
      traceYdAdjYd + 174*traceYeAdjYe + 242*traceYuAdjYu - 45*Sqr(g2)) - 25*(
      379*Quad(g2) - 108*(3*traceYdAdjYd + traceYeAdjYe + 3*traceYuAdjYu)*Sqr(
      g2) - 8*(80*(traceYdAdjYd + traceYuAdjYu)*Sqr(g3) - 3*(9*
      traceYdAdjYdYdAdjYd - 2*traceYdAdjYuYuAdjYd + 3*traceYeAdjYeYeAdjYe + 9*
      traceYuAdjYuYuAdjYu - 2*Sqr(Lambdax))))));

   return beta_v;
}

double Standard_model::calc_beta_v_three_loop(const Beta_traces&) const
{
   double beta_v;

   beta_v = 0;

   return beta_v;
}

double Standard_model::get_mass_matrix_VG() const
{
   const double mass_matrix_VG = Re(0);

   return mass_matrix_VG;
}

void Standard_model::calculate_M2VG()
{
   const auto mass_matrix_VG = get_mass_matrix_VG();
   M2VG = mass_matrix_VG;
}

double Standard_model::get_mass_matrix_Hp() const
{
   const double mass_matrix_Hp = Re(-mu2 + 0.5*Lambdax*Sqr(v) + 0.25*Sqr(
      g2)*Sqr(v));

   return mass_matrix_Hp;
}

void Standard_model::calculate_M2Hp()
{
   const auto mass_matrix_Hp = get_mass_matrix_Hp();
   M2Hp = mass_matrix_Hp;

   if (M2Hp < 0.) {
      problems.flag_running_tachyon(standard_model_info::Hp);
   }
}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fv() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fv;

   mass_matrix_Fv(0,0) = 0;
   mass_matrix_Fv(0,1) = 0;
   mass_matrix_Fv(0,2) = 0;
   mass_matrix_Fv(1,1) = 0;
   mass_matrix_Fv(1,2) = 0;
   mass_matrix_Fv(2,2) = 0;

   Symmetrize(mass_matrix_Fv);

   return mass_matrix_Fv;
}

void Standard_model::calculate_MFv()
{
   MFv.setConstant(0);
}

double Standard_model::get_mass_matrix_Ah() const
{
   const double mass_matrix_Ah = Re(0.25*(-4*mu2 + 2*Lambdax*Sqr(v) + Sqr
      (v)*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()))));

   return mass_matrix_Ah;
}

void Standard_model::calculate_M2Ah()
{
   const auto mass_matrix_Ah = get_mass_matrix_Ah();
   M2Ah = mass_matrix_Ah;

   if (M2Ah < 0.) {
      problems.flag_running_tachyon(standard_model_info::Ah);
   }
}

double Standard_model::get_mass_matrix_hh() const
{
   const double mass_matrix_hh = Re(-mu2 + 1.5*Lambdax*Sqr(v));

   return mass_matrix_hh;
}

void Standard_model::calculate_M2hh()
{
   const auto mass_matrix_hh = get_mass_matrix_hh();
   M2hh = mass_matrix_hh;

   if (M2hh < 0.) {
      problems.flag_running_tachyon(standard_model_info::hh);
   }
}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fd() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fd;

   mass_matrix_Fd(0,0) = 0.7071067811865475*v*Yd(0,0);
   mass_matrix_Fd(0,1) = 0.7071067811865475*v*Yd(1,0);
   mass_matrix_Fd(0,2) = 0.7071067811865475*v*Yd(2,0);
   mass_matrix_Fd(1,0) = 0.7071067811865475*v*Yd(0,1);
   mass_matrix_Fd(1,1) = 0.7071067811865475*v*Yd(1,1);
   mass_matrix_Fd(1,2) = 0.7071067811865475*v*Yd(2,1);
   mass_matrix_Fd(2,0) = 0.7071067811865475*v*Yd(0,2);
   mass_matrix_Fd(2,1) = 0.7071067811865475*v*Yd(1,2);
   mass_matrix_Fd(2,2) = 0.7071067811865475*v*Yd(2,2);

   return mass_matrix_Fd;
}

void Standard_model::calculate_MFd()
{
   const auto mass_matrix_Fd(get_mass_matrix_Fd());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fd, eigenvalue_error > precision * Abs(MFd(0)));
#else
   fs_svd(mass_matrix_Fd, MFd, Vd, Ud);
#endif

}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fu() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fu;

   mass_matrix_Fu(0,0) = 0.7071067811865475*v*Yu(0,0);
   mass_matrix_Fu(0,1) = 0.7071067811865475*v*Yu(1,0);
   mass_matrix_Fu(0,2) = 0.7071067811865475*v*Yu(2,0);
   mass_matrix_Fu(1,0) = 0.7071067811865475*v*Yu(0,1);
   mass_matrix_Fu(1,1) = 0.7071067811865475*v*Yu(1,1);
   mass_matrix_Fu(1,2) = 0.7071067811865475*v*Yu(2,1);
   mass_matrix_Fu(2,0) = 0.7071067811865475*v*Yu(0,2);
   mass_matrix_Fu(2,1) = 0.7071067811865475*v*Yu(1,2);
   mass_matrix_Fu(2,2) = 0.7071067811865475*v*Yu(2,2);

   return mass_matrix_Fu;
}

void Standard_model::calculate_MFu()
{
   const auto mass_matrix_Fu(get_mass_matrix_Fu());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fu, eigenvalue_error > precision * Abs(MFu(0)));
#else
   fs_svd(mass_matrix_Fu, MFu, Vu, Uu);
#endif

}

Eigen::Matrix<double,3,3> Standard_model::get_mass_matrix_Fe() const
{
   Eigen::Matrix<double,3,3> mass_matrix_Fe;

   mass_matrix_Fe(0,0) = 0.7071067811865475*v*Ye(0,0);
   mass_matrix_Fe(0,1) = 0.7071067811865475*v*Ye(1,0);
   mass_matrix_Fe(0,2) = 0.7071067811865475*v*Ye(2,0);
   mass_matrix_Fe(1,0) = 0.7071067811865475*v*Ye(0,1);
   mass_matrix_Fe(1,1) = 0.7071067811865475*v*Ye(1,1);
   mass_matrix_Fe(1,2) = 0.7071067811865475*v*Ye(2,1);
   mass_matrix_Fe(2,0) = 0.7071067811865475*v*Ye(0,2);
   mass_matrix_Fe(2,1) = 0.7071067811865475*v*Ye(1,2);
   mass_matrix_Fe(2,2) = 0.7071067811865475*v*Ye(2,2);

   return mass_matrix_Fe;
}

void Standard_model::calculate_MFe()
{
   const auto mass_matrix_Fe(get_mass_matrix_Fe());


#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue, eigenvalue_error);
   problems.flag_bad_mass(standard_model_info::Fe, eigenvalue_error > precision * Abs(MFe(0)));
#else
   fs_svd(mass_matrix_Fe, MFe, Ve, Ue);
#endif

}

double Standard_model::get_mass_matrix_VWp() const
{
   const double mass_matrix_VWp = Re(0.25*Sqr(g2)*Sqr(v));

   return mass_matrix_VWp;
}

void Standard_model::calculate_M2VWp()
{
   const auto mass_matrix_VWp = get_mass_matrix_VWp();
   M2VWp = mass_matrix_VWp;

   if (M2VWp < 0.) {
      problems.flag_running_tachyon(standard_model_info::VWp);
   }
}

Eigen::Matrix<double,2,2> Standard_model::get_mass_matrix_VPVZ() const
{
   Eigen::Matrix<double,2,2> mass_matrix_VPVZ;

   mass_matrix_VPVZ(0,0) = 0.15*Sqr(g1)*Sqr(v);
   mass_matrix_VPVZ(0,1) = -0.19364916731037085*g1*g2*Sqr(v);
   mass_matrix_VPVZ(1,1) = 0.25*Sqr(g2)*Sqr(v);

   Symmetrize(mass_matrix_VPVZ);

   return mass_matrix_VPVZ;
}

void Standard_model::calculate_M2VPVZ()
{
   const auto mass_matrix_VPVZ(get_mass_matrix_VPVZ());
   Eigen::Array<double,2,1> M2VPVZ;

#ifdef CHECK_EIGENVALUE_ERROR
   double eigenvalue_error;
   fs_diagonalize_hermitian(mass_matrix_VPVZ, M2VPVZ, ZZ,
      eigenvalue_error);
#else
   fs_diagonalize_hermitian(mass_matrix_VPVZ, M2VPVZ, ZZ);
#endif
   ZZ.transposeInPlace();
   normalize_to_interval(ZZ);

   M2VP = 0.;
   M2VZ = M2VPVZ(1);
}


double Standard_model::get_ewsb_eq_hh_1() const
{
   double result = Re(-(mu2*v) + 0.5*Cube(v)*Lambdax);

   return result;
}



double Standard_model::CphhHpconjHp() const
{
   const double result = -(v*Lambdax);

   return result;
}

double Standard_model::CpbargWpgZHp() const
{
   const double result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpbargZgWpconjHp() const
{
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpbargWpCgZconjHp() const
{
   const double result = 0.05*g2*v*(-5*g2*Cos(ThetaW()) + 3.872983346207417*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpbargZgWpCHp() const
{
   const double result = 0.25*g2*v*(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjHpVPVWp() const
{
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double Standard_model::CpconjHpVWpVZ() const
{
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double Standard_model::CpAhAhHpconjHp() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::CphhhhHpconjHp() const
{
   const double result = -Lambdax;

   return result;
}

double Standard_model::CpHpHpconjHpconjHp() const
{
   const double result = -2*Lambdax;

   return result;
}

std::complex<double> Standard_model::CpAhconjHpVWp() const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double Standard_model::CphhconjHpVWp() const
{
   const double result = 0.5*g2;

   return result;
}

double Standard_model::CpHpconjHpVP() const
{
   const double result = 0.5*(0.7745966692414834*g1*Cos(ThetaW()) + g2*Sin(
      ThetaW()));

   return result;
}

double Standard_model::CpHpconjHpVZ() const
{
   const double result = 0.5*(g2*Cos(ThetaW()) - 0.7745966692414834*g1*Sin(
      ThetaW()));

   return result;
}

double Standard_model::CpHpconjHpconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpHpconjHpVZVZ() const
{
   const std::complex<double> result = 0.1*(-7.745966692414834*g1*g2*Cos(ThetaW
      ())*Sin(ThetaW()) + 5*Sqr(g2)*Sqr(Cos(ThetaW())) + 3*Sqr(g1)*Sqr(Sin(ThetaW(
      ))));

   return result;
}

std::complex<double> Standard_model::CpbarFdFuconjHpPR(int gI1, int gI2) const
{
   const std::complex<double> result = SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu
      (gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFdFuconjHpPL(int gI1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,
      Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

double Standard_model::CpbarFeFvconjHpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFeFvconjHpPL(int gI1, int gI2) const
{
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,gI2))
      ;

   return result;
}

double Standard_model::CpAhAhhh() const
{
   const double result = -(v*Lambdax);

   return result;
}

std::complex<double> Standard_model::CpbargWpgWpAh() const
{
   const std::complex<double> result = std::complex<double>(0,-0.25)*v*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpbargWpCgWpCAh() const
{
   const std::complex<double> result = std::complex<double>(0,0.25)*v*Sqr(g2);

   return result;
}

double Standard_model::CpAhAhAhAh() const
{
   const double result = -3*Lambdax;

   return result;
}

double Standard_model::CpAhAhhhhh() const
{
   const double result = -Lambdax;

   return result;
}

std::complex<double> Standard_model::CpAhhhVZ() const
{
   const std::complex<double> result = std::complex<double>(0,-0.1)*(5*g2*Cos(
      ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpAhHpconjVWp() const
{
   const std::complex<double> result = std::complex<double>(0,0.5)*g2;

   return result;
}

double Standard_model::CpAhAhconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpAhAhVZVZ() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CpbarFdFdAhPR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(
      gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFdFdAhPL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*
      Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFeFeAhPR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(
      gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFeFeAhPL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*
      Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuFuAhPR(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,SUM(j1,0,2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1
      ,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFuFuAhPL(int gI1, int gI2) const
{
   const std::complex<double> result = std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))
      *Yu(j1,j2)));

   return result;
}

double Standard_model::Cphhhhhh() const
{
   const double result = -3*v*Lambdax;

   return result;
}

double Standard_model::CphhVZVZ() const
{
   const double result = 0.5*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*Sin
      (ThetaW()));

   return result;
}

double Standard_model::CphhconjVWpVWp() const
{
   const double result = 0.5*v*Sqr(g2);

   return result;
}

double Standard_model::CpbargWpgWphh() const
{
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double Standard_model::CpbargWpCgWpChh() const
{
   const double result = -0.25*v*Sqr(g2);

   return result;
}

double Standard_model::CpbargZgZhh() const
{
   const double result = -0.25*v*Sqr(g2*Cos(ThetaW()) + 0.7745966692414834*g1*
      Sin(ThetaW()));

   return result;
}

double Standard_model::Cphhhhhhhh() const
{
   const double result = -3*Lambdax;

   return result;
}

double Standard_model::CphhHpconjVWp() const
{
   const double result = -0.5*g2;

   return result;
}

double Standard_model::CphhhhconjVWpVWp() const
{
   const double result = 0.5*Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CphhhhVZVZ() const
{
   const std::complex<double> result = 0.1*(g1*Sin(ThetaW())*(7.745966692414834
      *g2*Cos(ThetaW()) + 3*g1*Sin(ThetaW())) + 5*Sqr(g2)*Sqr(Cos(ThetaW())));

   return result;
}

std::complex<double> Standard_model::CpbarFdFdhhPR(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yd(j1,j2))*Ud(gI2,j1))*Vd(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFdFdhhPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vd(
      gI2,j2))*SUM(j1,0,2,Conj(Ud(gI1,j1))*Yd(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFeFehhPR(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Ye(j1,j2))*Ue(gI2,j1))*Ve(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFeFehhPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Ve(
      gI2,j2))*SUM(j1,0,2,Conj(Ue(gI1,j1))*Ye(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuFuhhPR(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,SUM(j1,0,
      2,Conj(Yu(j1,j2))*Uu(gI2,j1))*Vu(gI1,j2));

   return result;
}

std::complex<double> Standard_model::CpbarFuFuhhPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*SUM(j2,0,2,Conj(Vu(
      gI2,j2))*SUM(j1,0,2,Conj(Uu(gI1,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpVGVGVG() const
{
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

std::complex<double> Standard_model::CpbargGgGVG() const
{
   const std::complex<double> result = std::complex<double>(0,-1)*g3;

   return result;
}

double Standard_model::CpbarFdFdVGPL(int gI1, int gI2) const
{
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double Standard_model::CpbarFdFdVGPR(int gI1, int gI2) const
{
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double Standard_model::CpbarFuFuVGPL(int gI1, int gI2) const
{
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double Standard_model::CpbarFuFuVGPR(int gI1, int gI2) const
{
   const double result = -(g3*KroneckerDelta(gI1,gI2));

   return result;
}

double Standard_model::CpVGVGVGVG1() const
{
   const double result = -16*Sqr(g3);

   return result;
}

double Standard_model::CpVGVGVGVG2() const
{
   const double result = 0;

   return result;
}

double Standard_model::CpVGVGVGVG3() const
{
   const double result = 16*Sqr(g3);

   return result;
}

double Standard_model::CpHpconjVWpVP() const
{
   const double result = 0.3872983346207417*g1*g2*v*Cos(ThetaW());

   return result;
}

double Standard_model::CpbargWpgWpVP() const
{
   const double result = g2*Sin(ThetaW());

   return result;
}

double Standard_model::CpbargWpCgWpCVP() const
{
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpHpconjHpVPVP() const
{
   const std::complex<double> result = 0.1*(g2*Sin(ThetaW())*(7.745966692414834
      *g1*Cos(ThetaW()) + 5*g2*Sin(ThetaW())) + 3*Sqr(g1)*Sqr(Cos(ThetaW())));

   return result;
}

double Standard_model::CpconjVWpVPVWp() const
{
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFdFdVPPL(int gI1, int gI2) const
{
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) - 3*g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFdFdVPPR(int gI1, int gI2) const
{
   const double result = 0.2581988897471611*g1*Cos(ThetaW())*KroneckerDelta(gI1
      ,gI2);

   return result;
}

double Standard_model::CpbarFeFeVPPL(int gI1, int gI2) const
{
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(0.7745966692414834*g1*Cos
      (ThetaW()) + g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFeFeVPPR(int gI1, int gI2) const
{
   const double result = 0.7745966692414834*g1*Cos(ThetaW())*KroneckerDelta(gI1
      ,gI2);

   return result;
}

double Standard_model::CpbarFuFuVPPL(int gI1, int gI2) const
{
   const double result = -0.16666666666666666*KroneckerDelta(gI1,gI2)*(
      0.7745966692414834*g1*Cos(ThetaW()) + 3*g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFuFuVPPR(int gI1, int gI2) const
{
   const double result = -0.5163977794943222*g1*Cos(ThetaW())*KroneckerDelta(
      gI1,gI2);

   return result;
}

double Standard_model::CpconjVWpVPVPVWp3() const
{
   const double result = -2*Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpVPVPVWp1() const
{
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpVPVPVWp2() const
{
   const double result = Sqr(g2)*Sqr(Sin(ThetaW()));

   return result;
}

double Standard_model::CpHpconjVWpVZ() const
{
   const double result = -0.3872983346207417*g1*g2*v*Sin(ThetaW());

   return result;
}

double Standard_model::CpbargWpgWpVZ() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpbargWpCgWpCVZ() const
{
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpVWpVZ() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpbarFdFdVZPL(int gI1, int gI2) const
{
   const double result = 0.16666666666666666*KroneckerDelta(gI1,gI2)*(3*g2*Cos(
      ThetaW()) + 0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFdFdVZPR(int gI1, int gI2) const
{
   const double result = -0.2581988897471611*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFeFeVZPL(int gI1, int gI2) const
{
   const double result = 0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) -
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFeFeVZPR(int gI1, int gI2) const
{
   const double result = -0.7745966692414834*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFuFuVZPL(int gI1, int gI2) const
{
   const double result = 0.03333333333333333*KroneckerDelta(gI1,gI2)*(-15*g2*
      Cos(ThetaW()) + 3.872983346207417*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFuFuVZPR(int gI1, int gI2) const
{
   const double result = 0.5163977794943222*g1*KroneckerDelta(gI1,gI2)*Sin(
      ThetaW());

   return result;
}

double Standard_model::CpbarFvFvVZPL(int gI1, int gI2) const
{
   const double result = -0.5*KroneckerDelta(gI1,gI2)*(g2*Cos(ThetaW()) +
      0.7745966692414834*g1*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbarFvFvVZPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpconjVWpVWpVZVZ1() const
{
   const double result = -2*Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpVWpVZVZ2() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpconjVWpVWpVZVZ3() const
{
   const double result = Sqr(g2)*Sqr(Cos(ThetaW()));

   return result;
}

double Standard_model::CpbargPgWpconjVWp() const
{
   const double result = -(g2*Sin(ThetaW()));

   return result;
}

double Standard_model::CpbargWpCgPconjVWp() const
{
   const double result = g2*Sin(ThetaW());

   return result;
}

double Standard_model::CpbargWpCgZconjVWp() const
{
   const double result = g2*Cos(ThetaW());

   return result;
}

double Standard_model::CpbargZgWpconjVWp() const
{
   const double result = -(g2*Cos(ThetaW()));

   return result;
}

std::complex<double> Standard_model::CpbarFdFuconjVWpPL(int gI1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      Vu(gI2,j1))*Vd(gI1,j1));

   return result;
}

double Standard_model::CpbarFdFuconjVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFeFvconjVWpPL(int gI1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7071067811865475*g2*Ve(gI1
      ,gI2),0);

   return result;
}

double Standard_model::CpbarFeFvconjVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpconjVWpconjVWpVWpVWp2() const
{
   const double result = -Sqr(g2);

   return result;
}

double Standard_model::CpconjVWpconjVWpVWpVWp1() const
{
   const double result = 2*Sqr(g2);

   return result;
}

double Standard_model::CpconjVWpconjVWpVWpVWp3() const
{
   const double result = -Sqr(g2);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Vd(gI1,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j1,0,2,Conj(Yd(j1,gO1))*Ud(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdhhPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(Vd(gI2,j2))*Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdhhPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Yd(j1,gO1))*Ud(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVGPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Ud(gI2,gO2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVGPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vd(gI2,gO1))),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVPPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.2581988897471611*g1*Cos(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVPPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      Vd(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVZPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.2581988897471611*g1*Sin(
      ThetaW())*Ud(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFdVZPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Vd(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vd(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFuconjHpPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-SUM(j2,0,2,Conj(Vu(gI2,j2))*
      Yd(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFdFuconjHpPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,SUM(j1,0,2,Conj(Yu(j1,gO1))*
      Uu(gI2,j1)),0);

   return result;
}

double Standard_model::CpbarUFdFuconjVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarUFdFuconjVWpPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      Vu(gI2,gO1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j2,0,2,Conj(Vu(gI1,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j1,0,2,Conj(Yu(j1,gO1))*Uu(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFdHpPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,SUM(j2,0,2,Conj(Vd(gI2,j2))*
      Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFdHpPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-SUM(j1,0,2,Conj(Yd(j1,gO1))*
      Ud(gI2,j1)),0);

   return result;
}

double Standard_model::CpbarUFuFdVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarUFuFdVWpPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      Vd(gI2,gO1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuhhPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(Vu(gI2,j2))*Yu(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuhhPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Yu(j1,gO1))*Uu(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVGPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Uu(gI2,gO2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVGPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-(g3*Conj(Vu(gI2,gO1))),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVPPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.5163977794943222*g1*Cos(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVPPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.12909944487358055*g1*Conj(
      Vu(gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVZPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5163977794943222*g1*Sin(
      ThetaW())*Uu(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFuFuVZPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.5*g2*Conj(Vu(gI2,gO1))*Cos
      (ThetaW()),0) + IF(gI2 < 3,0.12909944487358055*g1*Conj(Vu(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeAhPL(int gO2, int gI1) const
{
   const std::complex<double> result = IF(gO2 < 3,std::complex<double>(0.,
      0.7071067811865475)*SUM(j2,0,2,Conj(Ve(gI1,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeAhPR(int gO1, int gI1) const
{
   const std::complex<double> result = IF(gO1 < 3,std::complex<double>(0.,
      -0.7071067811865475)*SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI1,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFehhPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-0.7071067811865475*SUM(j2,0,
      2,Conj(Ve(gI2,j2))*Ye(gO2,j2)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFehhPR(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*SUM(j1,0,
      2,Conj(Ye(j1,gO1))*Ue(gI2,j1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeVPPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.7745966692414834*g1*Cos(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeVPPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.3872983346207417*g1*Conj(Ve
      (gI2,gO1))*Cos(ThetaW()),0) + IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Sin(ThetaW
      ()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeVZPR(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,-0.7745966692414834*g1*Sin(
      ThetaW())*Ue(gI2,gO2),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFeVZPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gI2 < 3,0.5*g2*Conj(Ve(gI2,gO1))*Cos(
      ThetaW()),0) + IF(gI2 < 3,-0.3872983346207417*g1*Conj(Ve(gI2,gO1))*Sin(
      ThetaW()),0);

   return result;
}

std::complex<double> Standard_model::CpbarUFeFvconjHpPL(int gO2, int gI2) const
{
   const std::complex<double> result = IF(gO2 < 3,-Ye(gO2,gI2),0);

   return result;
}

double Standard_model::CpbarUFeFvconjHpPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpbarUFeFvconjVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

double Standard_model::CpbarUFeFvconjVWpPL(int gO1, int gI2) const
{
   const double result = IF(gI2 < 3,-0.7071067811865475*g2*KroneckerDelta(gI2,
      gO1),0);

   return result;
}

double Standard_model::CpbarFvFeHpPL(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFvFeHpPR(int gO1, int gI2) const
{
   const std::complex<double> result = -SUM(j1,0,2,Conj(Ye(j1,gO1))*Ue(gI2,j1))
      ;

   return result;
}

double Standard_model::CpbarFvFeVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFvFeVWpPL(int gO1, int gI2) const
{
   const std::complex<double> result = IF(gO1 < 3,-0.7071067811865475*g2*Conj(
      Ve(gI2,gO1)),0);

   return result;
}

std::complex<double> Standard_model::CpbarFuFdHpPL(int gO2, int gI2) const
{
   const std::complex<double> result = SUM(j2,0,2,Conj(Vd(gI2,j2))*SUM(j1,0,2,
      Conj(Uu(gO2,j1))*Yu(j1,j2)));

   return result;
}

std::complex<double> Standard_model::CpbarFuFdHpPR(int gO1, int gI2) const
{
   const std::complex<double> result = -SUM(j2,0,2,SUM(j1,0,2,Conj(Yd(j1,j2))*
      Ud(gI2,j1))*Vu(gO1,j2));

   return result;
}

double Standard_model::CpbarFuFdVWpPR(int , int ) const
{
   const double result = 0;

   return result;
}

std::complex<double> Standard_model::CpbarFuFdVWpPL(int gO1, int gI2) const
{
   const std::complex<double> result = -0.7071067811865475*g2*SUM(j1,0,2,Conj(
      Vd(gI2,j1))*Vu(gO1,j1));

   return result;
}


std::complex<double> Standard_model::self_energy_Hp_1loop(double p2 ) const
{
   std::complex<double> result;

   result += -0.5*A0(M2Ah)*CpAhAhHpconjHp();
   result += -(B0(p2,M2VZ,M2VWp)*CpbargWpCgZconjHp()*CpbargZgWpCHp());
   result += -(B0(p2,M2VWp,M2VZ)*CpbargWpgZHp()*CpbargZgWpconjHp());
   result += 2*AbsSqr(CpconjHpVPVWp())*(-1 + 2*B0(p2,0,M2VWp));
   result += 2*AbsSqr(CpconjHpVWpVZ())*(-1 + 2*B0(p2,M2VWp,M2VZ));
   result += -0.5*A0(M2hh)*CphhhhHpconjHp();
   result += AbsSqr(CphhHpconjHp())*B0(p2,M2Hp,M2hh);
   result += -(A0(M2Hp)*CpHpHpconjHpconjHp());
   result += AbsSqr(CpAhconjHpVWp())*F0(p2,M2Ah,M2VWp);
   result += AbsSqr(CphhconjHpVWp())*F0(p2,M2hh,M2VWp);
   result += AbsSqr(CpHpconjHpVP())*F0(p2,M2Hp,0);
   result += AbsSqr(CpHpconjHpVZ())*F0(p2,M2Hp,M2VZ);
   result += -2*M2VWp*CpHpconjHpconjVWpVWp() + 4*A0(M2VWp)*CpHpconjHpconjVWpVWp
      ();
   result += (-M2VZ + 2*A0(M2VZ))*CpHpconjHpVZVZ();
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFuconjHpPL(gI1,gI2)) +
      AbsSqr(CpbarFdFuconjHpPR(gI1,gI2)))*G0(p2,Sqr(MFd(gI1)),Sqr(MFu(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFvconjHpPL(gI1,gI2)) +
      AbsSqr(CpbarFeFvconjHpPR(gI1,gI2)))*G0(p2,Sqr(MFe(gI1)),Sqr(MFv(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFd(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFdFuconjHpPR(gI1,gI2))*CpbarFdFuconjHpPL(gI1,gI2) + Conj(
      CpbarFdFuconjHpPL(gI1,gI2))*CpbarFdFuconjHpPR(gI1,gI2))*MFu(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFe(gI1)),Sqr(MFv(
      gI2)))*(Conj(CpbarFeFvconjHpPR(gI1,gI2))*CpbarFeFvconjHpPL(gI1,gI2) + Conj(
      CpbarFeFvconjHpPL(gI1,gI2))*CpbarFeFvconjHpPR(gI1,gI2))*MFv(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Ah_1loop(double p2 ) const
{
   std::complex<double> result;

   result += -0.5*A0(M2Ah)*CpAhAhAhAh();
   result += AbsSqr(CpAhAhhh())*B0(p2,M2hh,M2Ah);
   result += -0.5*A0(M2hh)*CpAhAhhhhh();
   result += -(A0(M2Hp)*CpAhAhHpconjHp());
   result += -(B0(p2,M2VWp,M2VWp)*Sqr(CpbargWpCgWpCAh()));
   result += -(B0(p2,M2VWp,M2VWp)*Sqr(CpbargWpgWpAh()));
   result += AbsSqr(CpAhhhVZ())*F0(p2,M2hh,M2VZ);
   result += 2*AbsSqr(CpAhHpconjVWp())*F0(p2,M2Hp,M2VWp);
   result += -2*M2VWp*CpAhAhconjVWpVWp() + 4*A0(M2VWp)*CpAhAhconjVWpVWp();
   result += (-M2VZ + 2*A0(M2VZ))*CpAhAhVZVZ();
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdAhPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdAhPR(gI1,gI2)))*G0(p2,Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeAhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeAhPR(gI1,gI2)))*G0(p2,Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuAhPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuAhPR(gI1,gI2)))*G0(p2,Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdAhPR(gI1,gI2))*CpbarFdFdAhPL(gI1,gI2) + Conj(
      CpbarFdFdAhPL(gI1,gI2))*CpbarFdFdAhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFeAhPR(gI1,gI2))*CpbarFeFeAhPL(gI1,gI2) + Conj(
      CpbarFeFeAhPL(gI1,gI2))*CpbarFeFeAhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuAhPR(gI1,gI2))*CpbarFuFuAhPL(gI1,gI2) + Conj(
      CpbarFuFuAhPL(gI1,gI2))*CpbarFuFuAhPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_hh_1loop(double p2 ) const
{
   std::complex<double> result;

   result += 0.5*AbsSqr(CpAhAhhh())*B0(p2,M2Ah,M2Ah);
   result += -0.5*A0(M2Ah)*CpAhAhhhhh();
   result += -(B0(p2,M2VWp,M2VWp)*Sqr(CpbargWpCgWpChh()));
   result += -(B0(p2,M2VWp,M2VWp)*Sqr(CpbargWpgWphh()));
   result += -(B0(p2,M2VZ,M2VZ)*Sqr(CpbargZgZhh()));
   result += 2*AbsSqr(CphhconjVWpVWp())*(-1 + 2*B0(p2,M2VWp,M2VWp));
   result += 0.5*AbsSqr(Cphhhhhh())*B0(p2,M2hh,M2hh);
   result += -0.5*A0(M2hh)*Cphhhhhhhh();
   result += -(A0(M2Hp)*CphhhhHpconjHp());
   result += AbsSqr(CphhHpconjHp())*B0(p2,M2Hp,M2Hp);
   result += AbsSqr(CphhVZVZ())*(-1 + 2*B0(p2,M2VZ,M2VZ));
   result += AbsSqr(CpAhhhVZ())*F0(p2,M2Ah,M2VZ);
   result += 2*AbsSqr(CphhHpconjVWp())*F0(p2,M2Hp,M2VWp);
   result += -2*M2VWp*CphhhhconjVWpVWp() + 4*A0(M2VWp)*CphhhhconjVWpVWp();
   result += (-M2VZ + 2*A0(M2VZ))*CphhhhVZVZ();
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdhhPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdhhPR(gI1,gI2)))*G0(p2,Sqr(MFd(gI1)),Sqr(MFd(gI2)))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFehhPL(gI1,gI2)) + AbsSqr(
      CpbarFeFehhPR(gI1,gI2)))*G0(p2,Sqr(MFe(gI1)),Sqr(MFe(gI2)))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuhhPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuhhPR(gI1,gI2)))*G0(p2,Sqr(MFu(gI1)),Sqr(MFu(gI2)))));
   result += -6*SUM(gI1,0,2,MFd(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFd(gI1)),Sqr(MFd(
      gI2)))*(Conj(CpbarFdFdhhPR(gI1,gI2))*CpbarFdFdhhPL(gI1,gI2) + Conj(
      CpbarFdFdhhPL(gI1,gI2))*CpbarFdFdhhPR(gI1,gI2))*MFd(gI2)));
   result += -2*SUM(gI1,0,2,MFe(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFe(gI1)),Sqr(MFe(
      gI2)))*(Conj(CpbarFeFehhPR(gI1,gI2))*CpbarFeFehhPL(gI1,gI2) + Conj(
      CpbarFeFehhPL(gI1,gI2))*CpbarFeFehhPR(gI1,gI2))*MFe(gI2)));
   result += -6*SUM(gI1,0,2,MFu(gI1)*SUM(gI2,0,2,B0(p2,Sqr(MFu(gI1)),Sqr(MFu(
      gI2)))*(Conj(CpbarFuFuhhPR(gI1,gI2))*CpbarFuFuhhPL(gI1,gI2) + Conj(
      CpbarFuFuhhPL(gI1,gI2))*CpbarFuFuhhPR(gI1,gI2))*MFu(gI2)));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VG_1loop(double p2 ) const
{
   std::complex<double> result;

   result += 3*AbsSqr(CpbargGgGVG())*B00(p2,M2VG,M2VG);
   result += -(AbsSqr(CpVGVGVG())*(p2 + 6*p2*B0(p2,0,0) + 15*B00(p2,0,0)));
   result += 0;
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVGPL(gI1,gI2)) +
      AbsSqr(CpbarFdFdVGPR(gI1,gI2)))*H0(p2,Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(p2
      ,Sqr(MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVGPL(gI1,
      gI2))*CpbarFdFdVGPR(gI1,gI2))));
   result += 0.5*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVGPL(gI1,gI2)) +
      AbsSqr(CpbarFuFuVGPR(gI1,gI2)))*H0(p2,Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(p2
      ,Sqr(MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVGPL(gI1,
      gI2))*CpbarFuFuVGPR(gI1,gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VP_1loop(double p2 ) const
{
   std::complex<double> result;

   result += AbsSqr(CpbargWpCgWpCVP())*B00(p2,M2VWp,M2VWp);
   result += AbsSqr(CpbargWpgWpVP())*B00(p2,M2VWp,M2VWp);
   result += -(A0(M2VWp)*(CpconjVWpVPVPVWp1() + CpconjVWpVPVPVWp2() + 4*
      CpconjVWpVPVPVWp3()));
   result += -4*AbsSqr(CpHpconjHpVP())*B00(p2,M2Hp,M2Hp);
   result += A0(M2Hp)*CpHpconjHpVPVP();
   result += 2*AbsSqr(CpHpconjVWpVP())*B0(p2,M2VWp,M2Hp);
   result += 2*M2VWp*CpconjVWpVPVPVWp3();
   result += -0.6666666666666666*AbsSqr(CpconjVWpVPVWp())*(-6*M2VWp + p2 + 3*A0
      (M2VWp) + 3*(M2VWp + 2*p2)*B0(p2,M2VWp,M2VWp) + 15*B00(p2,M2VWp,M2VWp));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVPPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdVPPR(gI1,gI2)))*H0(p2,Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(p2,Sqr(
      MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVPPL(gI1,gI2))*
      CpbarFdFdVPPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVPPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVPPR(gI1,gI2)))*H0(p2,Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(p2,Sqr(
      MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVPPL(gI1,gI2))*
      CpbarFeFeVPPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVPPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuVPPR(gI1,gI2)))*H0(p2,Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(p2,Sqr(
      MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVPPL(gI1,gI2))*
      CpbarFuFuVPPR(gI1,gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VZ_1loop(double p2 ) const
{
   std::complex<double> result;

   result += 0.5*A0(M2Ah)*CpAhAhVZVZ();
   result += -4*AbsSqr(CpAhhhVZ())*B00(p2,M2Ah,M2hh);
   result += AbsSqr(CpbargWpCgWpCVZ())*B00(p2,M2VWp,M2VWp);
   result += AbsSqr(CpbargWpgWpVZ())*B00(p2,M2VWp,M2VWp);
   result += -(A0(M2VWp)*(4*CpconjVWpVWpVZVZ1() + CpconjVWpVWpVZVZ2() +
      CpconjVWpVWpVZVZ3()));
   result += 0.5*A0(M2hh)*CphhhhVZVZ();
   result += AbsSqr(CphhVZVZ())*B0(p2,M2VZ,M2hh);
   result += -4*AbsSqr(CpHpconjHpVZ())*B00(p2,M2Hp,M2Hp);
   result += A0(M2Hp)*CpHpconjHpVZVZ();
   result += 2*AbsSqr(CpHpconjVWpVZ())*B0(p2,M2VWp,M2Hp);
   result += 2*M2VWp*CpconjVWpVWpVZVZ1();
   result += -0.6666666666666666*AbsSqr(CpconjVWpVWpVZ())*(-6*M2VWp + p2 + 3*A0
      (M2VWp) + 3*(M2VWp + 2*p2)*B0(p2,M2VWp,M2VWp) + 15*B00(p2,M2VWp,M2VWp));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFdVZPL(gI1,gI2)) + AbsSqr
      (CpbarFdFdVZPR(gI1,gI2)))*H0(p2,Sqr(MFd(gI1)),Sqr(MFd(gI2))) + 4*B0(p2,Sqr(
      MFd(gI1)),Sqr(MFd(gI2)))*MFd(gI1)*MFd(gI2)*Re(Conj(CpbarFdFdVZPL(gI1,gI2))*
      CpbarFdFdVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFeVZPL(gI1,gI2)) + AbsSqr(
      CpbarFeFeVZPR(gI1,gI2)))*H0(p2,Sqr(MFe(gI1)),Sqr(MFe(gI2))) + 4*B0(p2,Sqr(
      MFe(gI1)),Sqr(MFe(gI2)))*MFe(gI1)*MFe(gI2)*Re(Conj(CpbarFeFeVZPL(gI1,gI2))*
      CpbarFeFeVZPR(gI1,gI2))));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFuFuVZPL(gI1,gI2)) + AbsSqr
      (CpbarFuFuVZPR(gI1,gI2)))*H0(p2,Sqr(MFu(gI1)),Sqr(MFu(gI2))) + 4*B0(p2,Sqr(
      MFu(gI1)),Sqr(MFu(gI2)))*MFu(gI1)*MFu(gI2)*Re(Conj(CpbarFuFuVZPL(gI1,gI2))*
      CpbarFuFuVZPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFvFvVZPL(gI1,gI2)) + AbsSqr(
      CpbarFvFvVZPR(gI1,gI2)))*H0(p2,Sqr(MFv(gI1)),Sqr(MFv(gI2))) + 4*B0(p2,Sqr(
      MFv(gI1)),Sqr(MFv(gI2)))*MFv(gI1)*MFv(gI2)*Re(Conj(CpbarFvFvVZPL(gI1,gI2))*
      CpbarFvFvVZPR(gI1,gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_VWp_1loop(double p2 ) const
{
   std::complex<double> result;

   result += 0.5*A0(M2Ah)*CpAhAhconjVWpVWp();
   result += -4*AbsSqr(CpAhHpconjVWp())*B00(p2,M2Ah,M2Hp);
   result += AbsSqr(CpbargPgWpconjVWp())*B00(p2,M2VWp,M2VP);
   result += AbsSqr(CpbargWpCgPconjVWp())*B00(p2,M2VP,M2VWp);
   result += AbsSqr(CpbargWpCgZconjVWp())*B00(p2,M2VZ,M2VWp);
   result += AbsSqr(CpbargZgWpconjVWp())*B00(p2,M2VWp,M2VZ);
   result += -(A0(M2VWp)*(CpconjVWpconjVWpVWpVWp1() + 4*CpconjVWpconjVWpVWpVWp2
      () + CpconjVWpconjVWpVWpVWp3()));
   result += 0;
   result += AbsSqr(CphhconjVWpVWp())*B0(p2,M2VWp,M2hh);
   result += 0.5*A0(M2hh)*CphhhhconjVWpVWp();
   result += -4*AbsSqr(CphhHpconjVWp())*B00(p2,M2hh,M2Hp);
   result += A0(M2Hp)*CpHpconjHpconjVWpVWp();
   result += AbsSqr(CpHpconjVWpVP())*B0(p2,0,M2Hp);
   result += AbsSqr(CpHpconjVWpVZ())*B0(p2,M2VZ,M2Hp);
   result += 2*M2VWp*CpconjVWpconjVWpVWpVWp2();
   result += -(AbsSqr(CpconjVWpVPVWp())*(-2*M2VWp + 0.6666666666666666*p2 + A0(
      M2VWp) + (M2VWp + 4*p2)*B0(p2,M2VWp,0) + 10*B00(p2,M2VWp,0)));
   result += 0.5*(2*M2VZ*CpconjVWpVWpVZVZ1() - A0(M2VZ)*(4*CpconjVWpVWpVZVZ1()
      + CpconjVWpVWpVZVZ2() + CpconjVWpVWpVZVZ3()));
   result += -(AbsSqr(CpconjVWpVWpVZ())*(-2*(M2VWp + M2VZ - 0.3333333333333333*
      p2) + A0(M2VWp) + A0(M2VZ) + (M2VWp + M2VZ + 4*p2)*B0(p2,M2VZ,M2VWp) + 10*
      B00(p2,M2VZ,M2VWp)));
   result += 3*SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFdFuconjVWpPL(gI1,gI2)) +
      AbsSqr(CpbarFdFuconjVWpPR(gI1,gI2)))*H0(p2,Sqr(MFd(gI1)),Sqr(MFu(gI2))) + 4*
      B0(p2,Sqr(MFd(gI1)),Sqr(MFu(gI2)))*MFd(gI1)*MFu(gI2)*Re(Conj(
      CpbarFdFuconjVWpPL(gI1,gI2))*CpbarFdFuconjVWpPR(gI1,gI2))));
   result += SUM(gI1,0,2,SUM(gI2,0,2,(AbsSqr(CpbarFeFvconjVWpPL(gI1,gI2)) +
      AbsSqr(CpbarFeFvconjVWpPR(gI1,gI2)))*H0(p2,Sqr(MFe(gI1)),Sqr(MFv(gI2))) + 4*
      B0(p2,Sqr(MFe(gI1)),Sqr(MFv(gI2)))*MFe(gI1)*MFv(gI2)*Re(Conj(
      CpbarFeFvconjVWpPL(gI1,gI2))*CpbarFeFvconjVWpPR(gI1,gI2))));

   return result * oneOver16PiSqr;

}

std::complex<double> Standard_model::self_energy_Fd_1loop_1(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarUFdFdAhPL(gO2,gI1)
      )*CpbarUFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarUFdFdhhPL(gO2,gI2)
      )*CpbarUFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPR
      (gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),M2VZ))*Conj(
      CpbarUFdFdVZPR(gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarUFdFuconjHpPL(gO2,
      gI2))*CpbarUFdFuconjHpPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_1(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PR(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarUFdFdAhPR(gO2
      ,gI1))*CpbarUFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarUFdFdhhPR(gO2
      ,gI2))*CpbarUFdFdhhPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPL(gO2,gI2))*CpbarUFdFdVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPL(
      gO2,gI2))*CpbarUFdFdVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VZ))*Conj(CpbarUFdFdVZPL
      (gO2,gI2))*CpbarUFdFdVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarUFdFuconjHpPR
      (gO2,gI2))*CpbarUFdFuconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarUFdFuconjVWpPL(gO2,gI2))*CpbarUFdFuconjVWpPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PR(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PL(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarUFdFdAhPL(gO2
      ,gI1))*CpbarUFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarUFdFdhhPL(gO2
      ,gI2))*CpbarUFdFdhhPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),0))*
      Conj(CpbarUFdFdVGPR(gO2,gI2))*CpbarUFdFdVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),0))*Conj(CpbarUFdFdVPPR(
      gO2,gI2))*CpbarUFdFdVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VZ))*Conj(CpbarUFdFdVZPR
      (gO2,gI2))*CpbarUFdFdVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarUFdFuconjHpPL
      (gO2,gI2))*CpbarUFdFuconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarUFdFuconjVWpPR(gO2,gI2))*CpbarUFdFuconjVWpPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PL(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_1(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPL(gO2,gI1)
      )*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPL(gO2,gI2)
      )*CpbarUFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPL(gO2,gI2)
      )*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -5.333333333333333*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR
      (gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),M2VZ))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_1(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPR(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPR(gO2
      ,gI2))*CpbarUFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPR(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPL(gO2,gI2))*CpbarUFuFuVGPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarUFuFuVZPL
      (gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PR(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPL(gO2
      ,gI2))*CpbarUFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPL(gO1,gI2));
   result += -1.3333333333333333*SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*
      Conj(CpbarUFuFuVGPR(gO2,gI2))*CpbarUFuFuVGPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarUFuFuVZPR
      (gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PL(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_1(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarUFeFeAhPL(gO2,gI1)
      )*CpbarUFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarUFeFehhPL(gO2,gI2)
      )*CpbarUFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPR
      (gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFe(gI2)),M2VZ))*Conj(
      CpbarUFeFeVZPR(gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarUFeFvconjHpPL(gO2,
      gI2))*CpbarUFeFvconjHpPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_1(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PR(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarUFeFeAhPR(gO2
      ,gI1))*CpbarUFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarUFeFehhPR(gO2
      ,gI2))*CpbarUFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPL(
      gO2,gI2))*CpbarUFeFeVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VZ))*Conj(CpbarUFeFeVZPL
      (gO2,gI2))*CpbarUFeFeVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarUFeFvconjHpPR
      (gO2,gI2))*CpbarUFeFvconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarUFeFvconjVWpPL(gO2,gI2))*CpbarUFeFvconjVWpPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PR(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PL(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarUFeFeAhPL(gO2
      ,gI1))*CpbarUFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarUFeFehhPL(gO2
      ,gI2))*CpbarUFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),0))*Conj(CpbarUFeFeVPPR(
      gO2,gI2))*CpbarUFeFeVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VZ))*Conj(CpbarUFeFeVZPR
      (gO2,gI2))*CpbarUFeFeVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarUFeFvconjHpPL
      (gO2,gI2))*CpbarUFeFvconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarUFeFvconjVWpPR(gO2,gI2))*CpbarUFeFvconjVWpPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PL(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fv_1loop_1(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI2,0,2,B0(p2,Sqr(MFe(gI2)),M2Hp)*Conj(CpbarFvFeHpPL(gO2,gI2))
      *CpbarFvFeHpPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFe(gI2)),M2VWp))*Conj(
      CpbarFvFeVWpPR(gO2,gI2))*CpbarFvFeVWpPL(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFv(gI2)),M2VZ))*Conj(
      CpbarFvFvVZPR(gO2,gI2))*CpbarFvFvVZPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fv_1loop_1(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_1(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fv_1loop_PR(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2Hp)*Conj(CpbarFvFeHpPR(gO2,
      gI2))*CpbarFvFeHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VWp))*Conj(
      CpbarFvFeVWpPL(gO2,gI2))*CpbarFvFeVWpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VZ))*Conj(CpbarFvFvVZPL(
      gO2,gI2))*CpbarFvFvVZPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fv_1loop_PR(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PR(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fv_1loop_PL(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2Hp)*Conj(CpbarFvFeHpPL(gO2,
      gI2))*CpbarFvFeHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VWp))*Conj(
      CpbarFvFeVWpPR(gO2,gI2))*CpbarFvFeVWpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VZ))*Conj(CpbarFvFvVZPR(
      gO2,gI2))*CpbarFvFvVZPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fv_1loop_PL(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fv_1loop_PL(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarFdFdAhPL(gO2,gI1))
      *CpbarFdFdAhPR(gO1,gI1)*MFd(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarFdFdhhPL(gO2,gI2))
      *CpbarFdFdhhPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),M2VZ))*Conj(
      CpbarFdFdVZPR(gO2,gI2))*CpbarFdFdVZPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarFdFuconjHpPL(gO2,
      gI2))*CpbarFdFuconjHpPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_1_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_1_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarFdFdAhPR(gO2,
      gI1))*CpbarFdFdAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarFdFdhhPR(gO2,
      gI2))*CpbarFdFdhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VZ))*Conj(CpbarFdFdVZPL(
      gO2,gI2))*CpbarFdFdVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarFdFuconjHpPR(
      gO2,gI2))*CpbarFdFuconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarFdFuconjVWpPL(gO2,gI2))*CpbarFdFuconjVWpPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PR_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PR_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fd_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFd(gI1)),M2Ah)*Conj(CpbarFdFdAhPL(gO2,
      gI1))*CpbarFdFdAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2hh)*Conj(CpbarFdFdhhPL(gO2,
      gI2))*CpbarFdFdhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VZ))*Conj(CpbarFdFdVZPR(
      gO2,gI2))*CpbarFdFdVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2Hp)*Conj(CpbarFdFuconjHpPL(
      gO2,gI2))*CpbarFdFuconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VWp))*Conj(
      CpbarFdFuconjVWpPR(gO2,gI2))*CpbarFdFuconjVWpPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fd_1loop_PL_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fd_1loop_PL_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarFeFeAhPL(gO2,gI1))
      *CpbarFeFeAhPR(gO1,gI1)*MFe(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarFeFehhPL(gO2,gI2))
      *CpbarFeFehhPR(gO1,gI2)*MFe(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFe(gI2)),M2VZ))*Conj(
      CpbarFeFeVZPR(gO2,gI2))*CpbarFeFeVZPL(gO1,gI2)*MFe(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarFeFvconjHpPL(gO2,
      gI2))*CpbarFeFvconjHpPR(gO1,gI2)*MFv(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2)*MFv(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_1_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_1_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarFeFeAhPR(gO2,
      gI1))*CpbarFeFeAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarFeFehhPR(gO2,
      gI2))*CpbarFeFehhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VZ))*Conj(CpbarFeFeVZPL(
      gO2,gI2))*CpbarFeFeVZPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarFeFvconjHpPR(
      gO2,gI2))*CpbarFeFvconjHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarFeFvconjVWpPL(gO2,gI2))*CpbarFeFvconjVWpPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PR_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PR_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fe_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFe(gI1)),M2Ah)*Conj(CpbarFeFeAhPL(gO2,
      gI1))*CpbarFeFeAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFe(gI2)),M2hh)*Conj(CpbarFeFehhPL(gO2,
      gI2))*CpbarFeFehhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFe(gI2)),M2VZ))*Conj(CpbarFeFeVZPR(
      gO2,gI2))*CpbarFeFeVZPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFv(gI2)),M2Hp)*Conj(CpbarFeFvconjHpPL(
      gO2,gI2))*CpbarFeFvconjHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFv(gI2)),M2VWp))*Conj(
      CpbarFeFvconjVWpPR(gO2,gI2))*CpbarFeFvconjVWpPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fe_1loop_PL_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fe_1loop_PL_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarFuFuAhPL(gO2,gI1))
      *CpbarFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarFuFdHpPL(gO2,gI2))
      *CpbarFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarFuFuhhPL(gO2,gI2))
      *CpbarFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR(
      gO2,gI2))*CpbarFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),M2VZ))*Conj(
      CpbarFuFuVZPR(gO2,gI2))*CpbarFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_1_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarFuFuAhPR(gO2,
      gI1))*CpbarFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarFuFdHpPR(gO2,
      gI2))*CpbarFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarFuFdVWpPL(gO2,gI2))*CpbarFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarFuFuhhPR(gO2,
      gI2))*CpbarFuFuhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPL(gO2
      ,gI2))*CpbarFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarFuFuVZPL(
      gO2,gI2))*CpbarFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PR_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarFuFuAhPL(gO2,
      gI1))*CpbarFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarFuFdHpPL(gO2,
      gI2))*CpbarFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarFuFdVWpPR(gO2,gI2))*CpbarFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarFuFuhhPL(gO2,
      gI2))*CpbarFuFuhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarFuFuVPPR(gO2
      ,gI2))*CpbarFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarFuFuVZPR(
      gO2,gI2))*CpbarFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PL_heavy_rotated(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy_rotated(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_1_heavy(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += SUM(gI1,0,2,B0(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPL(gO2,gI1)
      )*CpbarUFuFuAhPR(gO1,gI1)*MFu(gI1));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPL(gO2,gI2)
      )*CpbarUFuFdHpPR(gO1,gI2)*MFd(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2)*MFd(gI2));
   result += SUM(gI2,0,2,B0(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPL(gO2,gI2)
      )*CpbarUFuFuhhPR(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR
      (gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2)*MFu(gI2));
   result += -4*SUM(gI2,0,2,(-0.5 + B0(p2,Sqr(MFu(gI2)),M2VZ))*Conj(
      CpbarUFuFuVZPR(gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2)*MFu(gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_1_heavy(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_1_heavy(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PR_heavy(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPR(gO2
      ,gI1))*CpbarUFuFuAhPR(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPR(gO2
      ,gI2))*CpbarUFuFdHpPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPL(gO2,gI2))*CpbarUFuFdVWpPL(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPR(gO2
      ,gI2))*CpbarUFuFuhhPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPL(
      gO2,gI2))*CpbarUFuFuVPPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarUFuFuVZPL
      (gO2,gI2))*CpbarUFuFuVZPL(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PR_heavy(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PR_heavy(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::self_energy_Fu_1loop_PL_heavy(double p2 , int gO1, int gO2) const
{
   std::complex<double> result;

   result += -0.5*SUM(gI1,0,2,B1(p2,Sqr(MFu(gI1)),M2Ah)*Conj(CpbarUFuFuAhPL(gO2
      ,gI1))*CpbarUFuFuAhPL(gO1,gI1));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFd(gI2)),M2Hp)*Conj(CpbarUFuFdHpPL(gO2
      ,gI2))*CpbarUFuFdHpPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFd(gI2)),M2VWp))*Conj(
      CpbarUFuFdVWpPR(gO2,gI2))*CpbarUFuFdVWpPR(gO1,gI2));
   result += -0.5*SUM(gI2,0,2,B1(p2,Sqr(MFu(gI2)),M2hh)*Conj(CpbarUFuFuhhPL(gO2
      ,gI2))*CpbarUFuFuhhPL(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),0))*Conj(CpbarUFuFuVPPR(
      gO2,gI2))*CpbarUFuFuVPPR(gO1,gI2));
   result += -SUM(gI2,0,2,(0.5 + B1(p2,Sqr(MFu(gI2)),M2VZ))*Conj(CpbarUFuFuVZPR
      (gO2,gI2))*CpbarUFuFuVZPR(gO1,gI2));

   return result * oneOver16PiSqr;

}

Eigen::Matrix<std::complex<double>,3,3> Standard_model::self_energy_Fu_1loop_PL_heavy(double p2) const
{
   Eigen::Matrix<std::complex<double>,3,3> self_energy;

   for (int i = 0; i < 3; i++)
      for (int k = 0; k < 3; k++)
         self_energy(i, k) = self_energy_Fu_1loop_PL_heavy(p2, i, k);

   return self_energy;
}

std::complex<double> Standard_model::tadpole_hh_1loop() const
{
   std::complex<double> result;

   result += -0.5*A0(M2Ah)*CpAhAhhh();
   result += A0(M2VWp)*CpbargWpCgWpChh();
   result += A0(M2VWp)*CpbargWpgWphh();
   result += A0(M2VZ)*CpbargZgZhh();
   result += -0.5*A0(M2hh)*Cphhhhhh();
   result += -(A0(M2Hp)*CphhHpconjHp());
   result += -2*M2VWp*CphhconjVWpVWp() + 4*A0(M2VWp)*CphhconjVWpVWp();
   result += (-M2VZ + 2*A0(M2VZ))*CphhVZVZ();
   result += 6*SUM(gI1,0,2,A0(Sqr(MFd(gI1)))*(CpbarFdFdhhPL(gI1,gI1) +
      CpbarFdFdhhPR(gI1,gI1))*MFd(gI1));
   result += 2*SUM(gI1,0,2,A0(Sqr(MFe(gI1)))*(CpbarFeFehhPL(gI1,gI1) +
      CpbarFeFehhPR(gI1,gI1))*MFe(gI1));
   result += 6*SUM(gI1,0,2,A0(Sqr(MFu(gI1)))*(CpbarFuFuhhPL(gI1,gI1) +
      CpbarFuFuhhPR(gI1,gI1))*MFu(gI1));

   return result * oneOver16PiSqr;

}




double Standard_model::self_energy_hh_2loop(double p2) const
{
   using namespace flexiblesusy::sm_twoloophiggs;

   const double mb = MFd(2);
   const double mt = MFu(2);
   const double mtau = MFe(2);
   const double yb = Yd(2,2);
   const double yt = Yu(2,2);
   const double ytau = Ye(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_2LOOP_CORRECTION_AT_AT) {
      self_energy -= delta_mh_2loop_at_at_sm(p2, scale, mt, yt, mb);
   }

   if (HIGGS_2LOOP_CORRECTION_AT_AS) {
      self_energy -= delta_mh_2loop_at_as_sm(p2, scale, mt, yt, gs);
   }

   if (HIGGS_2LOOP_CORRECTION_AB_AS) {
      self_energy -= delta_mh_2loop_ab_as_sm(p2, scale, mb, yb, gs);
   }

   if (HIGGS_2LOOP_CORRECTION_ATAU_ATAU) {
      self_energy -= delta_mh_2loop_atau_atau_sm(p2, scale, mtau, ytau);
   }

   return self_energy;
}

double Standard_model::self_energy_hh_3loop() const
{
   using namespace flexiblesusy::sm_threeloophiggs;

   const double mt = MFu(2);
   const double yt = Yu(2,2);
   const double gs = g3;
   const double mh2 = M2hh;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_3LOOP_CORRECTION_AT_AT_AT) {
      self_energy -= delta_mh_3loop_at_at_at_sm(scale, mt, yt, mh2);
   }

   if (HIGGS_3LOOP_CORRECTION_AT_AT_AS) {
      self_energy -= delta_mh_3loop_at_at_as_sm(scale, mt, yt, gs);
   }

   if (HIGGS_3LOOP_CORRECTION_AT_AS_AS) {
      self_energy -= delta_mh_3loop_at_as_as_sm(scale, mt, yt, gs);
   }

   return self_energy;
}

double Standard_model::self_energy_hh_4loop() const
{
   using namespace flexiblesusy::sm_fourloophiggs;

   const double mt = MFu(2);
   const double yt = Yu(2,2);
   const double gs = g3;
   const double scale = get_scale();
   double self_energy = 0.;

   if (HIGGS_4LOOP_CORRECTION_AT_AS_AS_AS) {
      self_energy -= delta_mh_4loop_at_as_as_as_sm(scale, mt, yt, gs);
   }

   return self_energy;
}




void Standard_model::calculate_M2VG_pole()
{
   // diagonalization with medium precision
   PHYSICAL(M2VG) = 0.;
}

void Standard_model::calculate_MFv_pole()
{
   // diagonalization with medium precision
   PHYSICAL(MFv).setConstant(0.);
}

void Standard_model::calculate_M2hh_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::hh))
      return;

   // diagonalization with high precision
   const auto number_of_mass_iterations = get_number_of_mass_iterations()
      ;
   int iteration = 0;
   double diff = 0.0;
   decltype(M2hh) old_M2hh(M2hh), new_M2hh(M2hh);

   do {
      const double M_tree(get_mass_matrix_hh());
      const double p2 = old_M2hh;
      double self_energy = Re(self_energy_hh_1loop(p2));
      if (pole_mass_loop_order > 1)
         self_energy += self_energy_hh_2loop(p2);
      if (pole_mass_loop_order > 2)
         self_energy += self_energy_hh_3loop();
      if (pole_mass_loop_order > 3)
         self_energy += self_energy_hh_4loop();
      const double mass_sqr = M_tree - self_energy;

      PHYSICAL(M2hh) = mass_sqr;

      new_M2hh = PHYSICAL(M2hh);
      diff = MaxRelDiff(new_M2hh, old_M2hh);
      old_M2hh = new_M2hh;
      iteration++;
   } while (diff > precision
            && iteration < number_of_mass_iterations);

   if (diff > precision)
      problems.flag_no_pole_mass_convergence(standard_model_info::hh);
   else
      problems.unflag_no_pole_mass_convergence(standard_model_info::hh);
}

void Standard_model::calculate_M2VP_pole()
{
   // diagonalization with medium precision
   PHYSICAL(M2VP) = 0.;
}

void Standard_model::calculate_M2VZ_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VZ))
      return;

   // diagonalization with medium precision
   const double M_tree(M2VZ);
   const double p2 = M2VZ;
   const double self_energy = Re(self_energy_VZ_1loop(p2));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VZ);

   PHYSICAL(M2VZ) = mass_sqr;
}

void Standard_model::calculate_MFd_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fd());
   for (int es = 0; es < 3; ++es) {
      const double p2 = AbsSqr(MFd(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fd_1loop_1(p2));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fd_1loop_PL(p2));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fd_1loop_PR(p2));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vd) mix_Vd;
      decltype(Ud) mix_Ud;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fd, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vd, mix_Ud);
   #endif
      if (es == 0) {
         PHYSICAL(Vd) = mix_Vd;
         PHYSICAL(Ud) = mix_Ud;
      }
      PHYSICAL(MFd(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_MFu_pole()
{
   // diagonalization with medium precision
   double qcd_1l = 0.;

   {
      const double currentScale = get_scale();
      qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(2))/Sqr(
         currentScale)))*Sqr(g3);
   }

   double qcd_2l = 0.;

   if (pole_mass_loop_order > 1 && TOP_POLE_QCD_CORRECTION > 0) {
      const double currentScale = get_scale();
      qcd_2l = -0.005284774766427138*Quad(g3) - 0.0032348537833770956*
         Log(Sqr(currentScale)/Sqr(MFu(2)))*Quad(g3) - 0.0008822328500119351*
         Quad(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(2))));
   }

   double qcd_3l = 0.;

   if (pole_mass_loop_order > 2 && TOP_POLE_QCD_CORRECTION > 1) {
      const double currentScale = get_scale();
      qcd_3l = -0.00003352082872926087*Power6(g3)*(35.702577217116016
         + 1.*Cube(Log(Sqr(currentScale)/Sqr(MFu(2)))) + 15.387410814884797*Log
         (Sqr(currentScale)/Sqr(MFu(2))) + 5.378787878787879*Sqr(Log(Sqr(
         currentScale)/Sqr(MFu(2)))));
   }

   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fu());
   for (int es = 0; es < 3; ++es) {
      const double p2 = AbsSqr(MFu(es));
      Eigen::Matrix<double,3,3> self_energy_1;
      Eigen::Matrix<double,3,3> self_energy_PL;
      Eigen::Matrix<double,3,3> self_energy_PR;
      for (int i1 = 0; i1 < 3; ++i1) {
         for (int i2 = 0; i2 < 3; ++i2) {
            if (i1 == 2 && i2 == 2) {
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1loop_1_heavy(p2,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_1loop_PL_heavy(p2,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_1loop_PR_heavy(p2,i1,i2));
            } else {
               self_energy_1(i1,i2)  = Re(
                  self_energy_Fu_1loop_1(p2,i1,i2));
               self_energy_PL(i1,i2) = Re(
                  self_energy_Fu_1loop_PL(p2,i1,i2));
               self_energy_PR(i1,i2) = Re(
                  self_energy_Fu_1loop_PR(p2,i1,i2));
            }
         }
      }
      Eigen::Matrix<double,3,3> delta_M(- self_energy_PR * M_tree -
         M_tree * self_energy_PL - self_energy_1);
      delta_M(2,2) -= M_tree(2,2) * (qcd_1l + qcd_2l + qcd_3l);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Vu) mix_Vu;
      decltype(Uu) mix_Uu;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fu, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Vu, mix_Uu);
   #endif
      if (es == 0) {
         PHYSICAL(Vu) = mix_Vu;
         PHYSICAL(Uu) = mix_Uu;
      }
      PHYSICAL(MFu(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_MFe_pole()
{
   // diagonalization with medium precision
   const Eigen::Matrix<double,3,3> M_tree(get_mass_matrix_Fe());
   for (int es = 0; es < 3; ++es) {
      const double p2 = AbsSqr(MFe(es));
      const Eigen::Matrix<double,3,3> self_energy_1  = Re(
         self_energy_Fe_1loop_1(p2));
      const Eigen::Matrix<double,3,3> self_energy_PL = Re(
         self_energy_Fe_1loop_PL(p2));
      const Eigen::Matrix<double,3,3> self_energy_PR = Re(
         self_energy_Fe_1loop_PR(p2));
      const Eigen::Matrix<double,3,3> delta_M(- self_energy_PR *
         M_tree - M_tree * self_energy_PL - self_energy_1);
      const Eigen::Matrix<double,3,3> M_loop(M_tree + delta_M);
      Eigen::Array<double,3,1> eigen_values;
      decltype(Ve) mix_Ve;
      decltype(Ue) mix_Ue;
   #ifdef CHECK_EIGENVALUE_ERROR
      double eigenvalue_error;
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue, eigenvalue_error);
      problems.flag_bad_mass(standard_model_info::Fe, eigenvalue_error > precision
         * Abs(eigen_values(0)));
   #else
      fs_svd(M_loop, eigen_values, mix_Ve, mix_Ue);
   #endif
      if (es == 0) {
         PHYSICAL(Ve) = mix_Ve;
         PHYSICAL(Ue) = mix_Ue;
      }
      PHYSICAL(MFe(es)) = Abs(eigen_values(es));
   }
}

void Standard_model::calculate_M2VWp_pole()
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VWp))
      return;

   // diagonalization with medium precision
   const double M_tree(M2VWp);
   const double p2 = M2VWp;
   const double self_energy = Re(self_energy_VWp_1loop(p2));
   const double mass_sqr = M_tree - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   PHYSICAL(M2VWp) = mass_sqr;
}

double Standard_model::calculate_M2VWp_pole(double p2)
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VWp))
      return 0.;

   const double self_energy = Re(self_energy_VWp_1loop(p2));
   const double mass_sqr = M2VWp - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VWp);

   return mass_sqr;
}

double Standard_model::calculate_M2VZ_pole(double p2)
{
   if (!force_output && problems.is_running_tachyon(standard_model_info::VZ))
      return 0.;

   const double self_energy = Re(self_energy_VZ_1loop(p2));
   const double mass_sqr = M2VZ - self_energy;

   if (mass_sqr < 0.)
      problems.flag_pole_tachyon(standard_model_info::VZ);

   return mass_sqr;
}


double Standard_model::calculate_MFv_DRbar(double, int) const
{
   return 0.0;
}

double Standard_model::calculate_MFe_DRbar(double m_sm_msbar, int idx) const
{
   const double p2 = AbsSqr(m_sm_msbar);
   const double self_energy_1  = Re(self_energy_Fe_1loop_1_heavy_rotated(
      p2, idx, idx));
   const double self_energy_PL = Re(self_energy_Fe_1loop_PL_heavy_rotated
      (p2, idx, idx));
   const double self_energy_PR = Re(self_energy_Fe_1loop_PR_heavy_rotated
      (p2, idx, idx));
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mf_1loop = - self_energy_1/m_sm_drbar -
      self_energy_PL - self_energy_PR;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mf_1loop);

   return m_susy_drbar;
}

double Standard_model::calculate_MFu_DRbar(double m_pole, int idx) const
{
   const double p2 = AbsSqr(m_pole);
   const double self_energy_1  = Re(self_energy_Fu_1loop_1_heavy_rotated(
      p2, idx, idx));
   const double self_energy_PL = Re(self_energy_Fu_1loop_PL_heavy_rotated
      (p2, idx, idx));
   const double self_energy_PR = Re(self_energy_Fu_1loop_PR_heavy_rotated
      (p2, idx, idx));

   const double currentScale = get_scale();
   double qcd_1l = 0., qcd_2l = 0., qcd_3l = 0.;

   qcd_1l = -0.008443431970194815*(4. - 3.*Log(Sqr(MFu(idx))/Sqr(
      currentScale)))*Sqr(g3);

   if (get_thresholds() > 1 && threshold_corrections.mt > 1) {
      const double q_2l = 0.005284774766427138*Quad(g3) +
         0.0032348537833770956*Log(Sqr(currentScale)/Sqr(MFu(idx)))*Quad(g3) +
         0.0008822328500119351*Quad(g3)*Sqr(Log(Sqr(currentScale)/Sqr(MFu(idx))
         ));

      qcd_2l = -q_2l + qcd_1l * qcd_1l;
   }

   if (get_thresholds() > 2 && threshold_corrections.mt > 2) {
      qcd_3l = -0.0008783313853540776*Power6(g3) -
         5.078913443827405e-6*Cube(Log(Sqr(currentScale)/Sqr(MFu(idx))))*Power6
         (g3) - 0.00011624292515063117*Log(Sqr(currentScale)/Sqr(MFu(idx)))*
         Power6(g3) - 0.00002183932780845784*Power6(g3)*Sqr(Log(Sqr(
         currentScale)/Sqr(MFu(idx))));
   }

   const double m_susy_drbar = m_pole + self_energy_1 + m_pole * (
      self_energy_PL + self_energy_PR + qcd_1l + qcd_2l + qcd_3l);

   return m_susy_drbar;
}

double Standard_model::calculate_MFd_DRbar(double m_sm_msbar, int idx) const
{
   const double p2 = SignedAbsSqr(m_sm_msbar);
   const double self_energy_1  = Re(self_energy_Fd_1loop_1_heavy_rotated(
      p2, idx, idx));
   const double self_energy_PL = Re(self_energy_Fd_1loop_PL_heavy_rotated
      (p2, idx, idx));
   const double self_energy_PR = Re(self_energy_Fd_1loop_PR_heavy_rotated
      (p2, idx, idx));
   const double m_tree = MFd(idx);
   const double drbar_conversion = 1;
   const double m_sm_drbar = m_sm_msbar * drbar_conversion;
   const double delta_mb_1loop = - self_energy_1/m_tree - self_energy_PL
      - self_energy_PR;
   double qcd_2l = 0.;

   const double m_susy_drbar = m_sm_drbar / (1.0 + delta_mb_1loop +
      qcd_2l);

   return m_susy_drbar;
}

double Standard_model::calculate_M2VP_DRbar(double)
{
   return 0.0;
}

double Standard_model::calculate_M2VZ_DRbar(double m_pole_sqr)
{
   const double p2 = m_pole_sqr;
   const double self_energy = Re(self_energy_VZ_1loop(p2));
   const double mass_sqr = m_pole_sqr + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(standard_model_info::VZ);
      return m_pole_sqr;
   }

   return mass_sqr;
}

double Standard_model::calculate_M2VWp_DRbar(double m_pole_sqr)
{
   const double p2 = m_pole_sqr;
   const double self_energy = Re(self_energy_VWp_1loop(p2));
   const double mass_sqr = m_pole_sqr + self_energy;

   if (mass_sqr < 0.) {
      problems.flag_pole_tachyon(standard_model_info::VWp);
      return m_pole_sqr;
   }

   return mass_sqr;
}

double Standard_model::calculate_M2hh_DRbar(double m_pole_sqr)
{
   const double p2 = m_pole_sqr;
   const double self_energy = Re(self_energy_hh_1loop(p2));
   const double tadpole = Re(tadpole_hh_1loop());
   const double mass_sqr = m_pole_sqr + self_energy - tadpole/v;

   if (mass_sqr < 0.) {
      problems.flag_running_tachyon(standard_model_info::hh);
      return m_pole_sqr;
   }

   return mass_sqr;
}

double Standard_model::ThetaW() const
{
   return ArcCos(Abs(ZZ(0,0)));
}

std::ostream& operator<<(std::ostream& ostr, const Standard_model& model)
{
   model.print(ostr);
   return ostr;
}

} // namespace standard_model

} // namespace flexiblesusy
