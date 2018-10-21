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
 * @file standard_model.hpp
 *
 * @brief contains class for Standard model running and self-energies
 *
 */

#ifndef STANDARD_MODEL_H
#define STANDARD_MODEL_H

#include "betafunction.hpp"
#include "standard_model_physical.hpp"
#include "loop_corrections.hpp"
#include "threshold_corrections.hpp"
#include "error.hpp"
#include "problems.hpp"
#include "config.h"
#include "physical_input.hpp"

#include <array>
#include <iosfwd>
#include <string>

#include <Eigen/Core>

namespace softsusy {
   class QedQcd;
} // namespace softsusy

namespace flexiblesusy {

class EWSB_solver;

namespace standard_model_info {
   enum Particles : int { VG, Hp, Fv, Ah, hh, Fd, Fu, Fe, VWp, VP, VZ,
      NUMBER_OF_PARTICLES };

   enum Masses : int { M2VG, M2Hp, MFv_1, MFv_2, MFv_3, M2Ah, M2hh, MFd_1,
      MFd_2, MFd_3, MFu_1, MFu_2, MFu_3, MFe_1, MFe_2, MFe_3, M2VWp, M2VP, M2VZ,
      NUMBER_OF_MASSES };

   enum Parameters : int { g1, g2, g3, Lambdax, Yu0_0, Yu0_1, Yu0_2, Yu1_0,
      Yu1_1, Yu1_2, Yu2_0, Yu2_1, Yu2_2, Yd0_0, Yd0_1, Yd0_2, Yd1_0, Yd1_1, Yd1_2,
      Yd2_0, Yd2_1, Yd2_2, Ye0_0, Ye0_1, Ye0_2, Ye1_0, Ye1_1, Ye1_2, Ye2_0, Ye2_1
      , Ye2_2, mu2, v, NUMBER_OF_PARAMETERS };

   enum Mixings : int { ReVd0_0, ImVd0_0, ReVd0_1, ImVd0_1, ReVd0_2, ImVd0_2,
      ReVd1_0, ImVd1_0, ReVd1_1, ImVd1_1, ReVd1_2, ImVd1_2, ReVd2_0, ImVd2_0,
      ReVd2_1, ImVd2_1, ReVd2_2, ImVd2_2, ReUd0_0, ImUd0_0, ReUd0_1, ImUd0_1,
      ReUd0_2, ImUd0_2, ReUd1_0, ImUd1_0, ReUd1_1, ImUd1_1, ReUd1_2, ImUd1_2,
      ReUd2_0, ImUd2_0, ReUd2_1, ImUd2_1, ReUd2_2, ImUd2_2, ReVu0_0, ImVu0_0,
      ReVu0_1, ImVu0_1, ReVu0_2, ImVu0_2, ReVu1_0, ImVu1_0, ReVu1_1, ImVu1_1,
      ReVu1_2, ImVu1_2, ReVu2_0, ImVu2_0, ReVu2_1, ImVu2_1, ReVu2_2, ImVu2_2,
      ReUu0_0, ImUu0_0, ReUu0_1, ImUu0_1, ReUu0_2, ImUu0_2, ReUu1_0, ImUu1_0,
      ReUu1_1, ImUu1_1, ReUu1_2, ImUu1_2, ReUu2_0, ImUu2_0, ReUu2_1, ImUu2_1,
      ReUu2_2, ImUu2_2, ReVe0_0, ImVe0_0, ReVe0_1, ImVe0_1, ReVe0_2, ImVe0_2,
      ReVe1_0, ImVe1_0, ReVe1_1, ImVe1_1, ReVe1_2, ImVe1_2, ReVe2_0, ImVe2_0,
      ReVe2_1, ImVe2_1, ReVe2_2, ImVe2_2, ReUe0_0, ImUe0_0, ReUe0_1, ImUe0_1,
      ReUe0_2, ImUe0_2, ReUe1_0, ImUe1_0, ReUe1_1, ImUe1_1, ReUe1_2, ImUe1_2,
      ReUe2_0, ImUe2_0, ReUe2_1, ImUe2_1, ReUe2_2, ImUe2_2, ZZ0_0, ZZ0_1, ZZ1_0,
      ZZ1_1, NUMBER_OF_MIXINGS };

   extern const double normalization_g1;
   extern const double normalization_g2;
   extern const double normalization_g3;

   extern const std::array<int, NUMBER_OF_PARTICLES> particle_multiplicities;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_names;
   extern const std::array<std::string, NUMBER_OF_PARTICLES> particle_latex_names;
   extern const std::array<std::string, NUMBER_OF_PARAMETERS> parameter_names;
   extern const std::array<std::string, NUMBER_OF_MIXINGS> particle_mixing_names;
   extern const std::string model_name;
   constexpr bool is_low_energy_model = false;
   constexpr bool is_supersymmetric_model = false;

   class Standard_model_particle_names : public Names {
   public:
      virtual ~Standard_model_particle_names() = default;
      virtual const std::string& get(int index) const override {
         return particle_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARAMETERS;
      }
   };

   class Standard_model_parameter_names : public Names {
   public:
      virtual ~Standard_model_parameter_names() = default;
      virtual const std::string& get(int index) const override {
         return parameter_names[index];
      }
      virtual int size() const override {
         return NUMBER_OF_PARTICLES;
      }
   };

   const Standard_model_particle_names  particle_names_getter{};
   const Standard_model_parameter_names parameter_names_getter{};

} // namespace standard_model_info

namespace standard_model {

template <class T>
class StandardModel;

/**
 * @class Standard_model
 * @brief model class with routines for SM running and self-energies
 */
class Standard_model : public Beta_function {
public:

   Standard_model();
   Standard_model(double scale_, double loops_, double thresholds_
   , double g1_, double g2_, double g3_, double Lambdax_, const Eigen::Matrix<
   double,3,3>& Yu_, const Eigen::Matrix<double,3,3>& Yd_, const Eigen::Matrix<
   double,3,3>& Ye_, double mu2_, double v_);
   Standard_model(const Standard_model&) = default;
   Standard_model(Standard_model&&) = default;

   virtual ~Standard_model() = default;

   Standard_model& operator=(const Standard_model&) = default;
   Standard_model& operator=(Standard_model&&) = default;

   /// number of EWSB equations
   static const int number_of_ewsb_equations = 1;

   void calculate_DRbar_masses();
   void calculate_pole_masses();
   void check_pole_masses_for_tachyons();
   void do_force_output(bool);
   bool do_force_output() const;
   void set_ewsb_iteration_precision(double);
   void set_ewsb_loop_order(int);
   void set_loop_corrections(const Loop_corrections&);
   const Loop_corrections& get_loop_corrections() const;
   void set_threshold_corrections(const Threshold_corrections&);
   const Threshold_corrections& get_threshold_corrections() const;
   void set_pole_mass_loop_order(int);
   int get_pole_mass_loop_order() const;
   void set_physical(const Standard_model_physical&);
   double get_ewsb_iteration_precision() const;
   double get_ewsb_loop_order() const;
   const Standard_model_physical& get_physical() const;
   Standard_model_physical& get_physical();
   const Problems& get_problems() const;
   Problems& get_problems();
   int solve_ewsb_tree_level();
   int solve_ewsb_one_loop();
   int solve_ewsb();            ///< solve EWSB at ewsb_loop_order level

   virtual Eigen::ArrayXd beta() const override;
   virtual Eigen::ArrayXd get() const override;
   void print(std::ostream& out = std::cerr) const;
   virtual void set(const Eigen::ArrayXd&) override;

   Standard_model calc_beta() const;
   Standard_model calc_beta(int) const;
   void clear();
   void clear_running_parameters();
   void clear_DRbar_parameters();
   void clear_problems();

   void calculate_spectrum();
   std::string name() const;
   virtual void run_to(double scale, double eps = -1.0) override;
   void set_precision(double);
   double get_precision() const;

   void set_physical_input(const Physical_input& input_) { input = input_; }
   const Physical_input& get_physical_input() const { return input; }
   Physical_input& get_physical_input() { return input; }

   void initialise_from_input(const softsusy::QedQcd&);

   void set_g1(double g1_) { g1 = g1_; }
   void set_g2(double g2_) { g2 = g2_; }
   void set_g3(double g3_) { g3 = g3_; }
   void set_Lambdax(double Lambdax_) { Lambdax = Lambdax_; }
   void set_Yu(const Eigen::Matrix<double,3,3>& Yu_) { Yu = Yu_; }
   void set_Yu(int i, int k, double value) { Yu(i,k) = value; }
   void set_Yd(const Eigen::Matrix<double,3,3>& Yd_) { Yd = Yd_; }
   void set_Yd(int i, int k, double value) { Yd(i,k) = value; }
   void set_Ye(const Eigen::Matrix<double,3,3>& Ye_) { Ye = Ye_; }
   void set_Ye(int i, int k, double value) { Ye(i,k) = value; }
   void set_mu2(double mu2_) { mu2 = mu2_; }
   void set_v(double v_) { v = v_; }

   double get_g1() const { return g1; }
   double get_g2() const { return g2; }
   double get_g3() const { return g3; }
   double get_Lambdax() const { return Lambdax; }
   const Eigen::Matrix<double,3,3>& get_Yu() const { return Yu; }
   double get_Yu(int i, int k) const { return Yu(i,k); }
   const Eigen::Matrix<double,3,3>& get_Yd() const { return Yd; }
   double get_Yd(int i, int k) const { return Yd(i,k); }
   const Eigen::Matrix<double,3,3>& get_Ye() const { return Ye; }
   double get_Ye(int i, int k) const { return Ye(i,k); }
   double get_mu2() const { return mu2; }
   double get_v() const { return v; }

   double get_M2VG() const { return M2VG; }
   double get_MVG() const;
   double get_M2Hp() const { return M2Hp; }
   double get_MHp() const;
   const Eigen::Array<double,3,1>& get_MFv() const { return MFv; }
   double get_MFv(int i) const { return MFv(i); }
   double get_M2Ah() const { return M2Ah; }
   double get_MAh() const;
   double get_M2hh() const { return M2hh; }
   double get_Mhh() const;
   const Eigen::Array<double,3,1>& get_MFd() const { return MFd; }
   double get_MFd(int i) const { return MFd(i); }
   const Eigen::Array<double,3,1>& get_MFu() const { return MFu; }
   double get_MFu(int i) const { return MFu(i); }
   const Eigen::Array<double,3,1>& get_MFe() const { return MFe; }
   double get_MFe(int i) const { return MFe(i); }
   double get_M2VWp() const { return M2VWp; }
   double get_MVWp() const;
   double get_M2VP() const { return M2VP; }
   double get_MVP() const;
   double get_M2VZ() const { return M2VZ; }
   double get_MVZ() const;

   const Eigen::Matrix<std::complex<double>,3,3>& get_Vd() const { return Vd; }
   const std::complex<double>& get_Vd(int i, int k) const { return Vd(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ud() const { return Ud; }
   const std::complex<double>& get_Ud(int i, int k) const { return Ud(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Vu() const { return Vu; }
   const std::complex<double>& get_Vu(int i, int k) const { return Vu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Uu() const { return Uu; }
   const std::complex<double>& get_Uu(int i, int k) const { return Uu(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ve() const { return Ve; }
   const std::complex<double>& get_Ve(int i, int k) const { return Ve(i,k); }
   const Eigen::Matrix<std::complex<double>,3,3>& get_Ue() const { return Ue; }
   const std::complex<double>& get_Ue(int i, int k) const { return Ue(i,k); }
   const Eigen::Matrix<double,2,2>& get_ZZ() const { return ZZ; }
   double get_ZZ(int i, int k) const { return ZZ(i,k); }

   double get_mass_matrix_VG() const;
   void calculate_M2VG();
   double get_mass_matrix_Hp() const;
   void calculate_M2Hp();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fv() const;
   void calculate_MFv();
   double get_mass_matrix_Ah() const;
   void calculate_M2Ah();
   double get_mass_matrix_hh() const;
   void calculate_M2hh();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fd() const;
   void calculate_MFd();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fu() const;
   void calculate_MFu();
   Eigen::Matrix<double,3,3> get_mass_matrix_Fe() const;
   void calculate_MFe();
   double get_mass_matrix_VWp() const;
   void calculate_M2VWp();
   Eigen::Matrix<double,2,2> get_mass_matrix_VPVZ() const;
   void calculate_M2VPVZ();

   double get_ewsb_eq_hh_1() const;
   double CphhHpconjHp() const;
   double CpbargWpgZHp() const;
   double CpbargZgWpconjHp() const;
   double CpbargWpCgZconjHp() const;
   double CpbargZgWpCHp() const;
   double CpconjHpVPVWp() const;
   double CpconjHpVWpVZ() const;
   double CpAhAhHpconjHp() const;
   double CphhhhHpconjHp() const;
   double CpHpHpconjHpconjHp() const;
   std::complex<double> CpAhconjHpVWp() const;
   double CphhconjHpVWp() const;
   double CpHpconjHpVP() const;
   double CpHpconjHpVZ() const;
   double CpHpconjHpconjVWpVWp() const;
   std::complex<double> CpHpconjHpVZVZ() const;
   std::complex<double> CpbarFdFuconjHpPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFuconjHpPL(int gI1, int gI2) const;
   double CpbarFeFvconjHpPR(int, int) const;
   std::complex<double> CpbarFeFvconjHpPL(int gI1, int gI2) const;
   double CpAhAhhh() const;
   std::complex<double> CpbargWpgWpAh() const;
   std::complex<double> CpbargWpCgWpCAh() const;
   double CpAhAhAhAh() const;
   double CpAhAhhhhh() const;
   std::complex<double> CpAhhhVZ() const;
   std::complex<double> CpAhHpconjVWp() const;
   double CpAhAhconjVWpVWp() const;
   std::complex<double> CpAhAhVZVZ() const;
   std::complex<double> CpbarFdFdAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFeAhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuAhPL(int gI1, int gI2) const;
   double Cphhhhhh() const;
   double CphhVZVZ() const;
   double CphhconjVWpVWp() const;
   double CpbargWpgWphh() const;
   double CpbargWpCgWpChh() const;
   double CpbargZgZhh() const;
   double Cphhhhhhhh() const;
   double CphhHpconjVWp() const;
   double CphhhhconjVWpVWp() const;
   std::complex<double> CphhhhVZVZ() const;
   std::complex<double> CpbarFdFdhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFdFdhhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFeFehhPL(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPR(int gI1, int gI2) const;
   std::complex<double> CpbarFuFuhhPL(int gI1, int gI2) const;
   std::complex<double> CpVGVGVG() const;
   std::complex<double> CpbargGgGVG() const;
   double CpbarFdFdVGPL(int gI1, int gI2) const;
   double CpbarFdFdVGPR(int gI1, int gI2) const;
   double CpbarFuFuVGPL(int gI1, int gI2) const;
   double CpbarFuFuVGPR(int gI1, int gI2) const;
   double CpVGVGVGVG1() const;
   double CpVGVGVGVG2() const;
   double CpVGVGVGVG3() const;
   double CpHpconjVWpVP() const;
   double CpbargWpgWpVP() const;
   double CpbargWpCgWpCVP() const;
   std::complex<double> CpHpconjHpVPVP() const;
   double CpconjVWpVPVWp() const;
   double CpbarFdFdVPPL(int gI1, int gI2) const;
   double CpbarFdFdVPPR(int gI1, int gI2) const;
   double CpbarFeFeVPPL(int gI1, int gI2) const;
   double CpbarFeFeVPPR(int gI1, int gI2) const;
   double CpbarFuFuVPPL(int gI1, int gI2) const;
   double CpbarFuFuVPPR(int gI1, int gI2) const;
   double CpconjVWpVPVPVWp3() const;
   double CpconjVWpVPVPVWp1() const;
   double CpconjVWpVPVPVWp2() const;
   double CpHpconjVWpVZ() const;
   double CpbargWpgWpVZ() const;
   double CpbargWpCgWpCVZ() const;
   double CpconjVWpVWpVZ() const;
   double CpbarFdFdVZPL(int gI1, int gI2) const;
   double CpbarFdFdVZPR(int gI1, int gI2) const;
   double CpbarFeFeVZPL(int gI1, int gI2) const;
   double CpbarFeFeVZPR(int gI1, int gI2) const;
   double CpbarFuFuVZPL(int gI1, int gI2) const;
   double CpbarFuFuVZPR(int gI1, int gI2) const;
   double CpbarFvFvVZPL(int gI1, int gI2) const;
   double CpbarFvFvVZPR(int, int) const;
   double CpconjVWpVWpVZVZ1() const;
   double CpconjVWpVWpVZVZ2() const;
   double CpconjVWpVWpVZVZ3() const;
   double CpbargPgWpconjVWp() const;
   double CpbargWpCgPconjVWp() const;
   double CpbargWpCgZconjVWp() const;
   double CpbargZgWpconjVWp() const;
   std::complex<double> CpbarFdFuconjVWpPL(int gI1, int gI2) const;
   double CpbarFdFuconjVWpPR(int, int) const;
   std::complex<double> CpbarFeFvconjVWpPL(int gI1, int gI2) const;
   double CpbarFeFvconjVWpPR(int, int) const;
   double CpconjVWpconjVWpVWpVWp2() const;
   double CpconjVWpconjVWpVWpVWp1() const;
   double CpconjVWpconjVWpVWpVWp3() const;
   std::complex<double> CpbarUFdFdAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFdFdAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFdFdhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFdVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFdVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFdFuconjHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFdFuconjHpPR(int gO1, int gI2) const;
   double CpbarUFdFuconjVWpPR(int, int) const;
   std::complex<double> CpbarUFdFuconjVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFuFuAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFuFdHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFdHpPR(int gO1, int gI2) const;
   double CpbarUFuFdVWpPR(int, int) const;
   std::complex<double> CpbarUFuFdVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuhhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuhhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVGPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVGPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFuFuVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFuFuVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeAhPL(int gO2, int gI1) const;
   std::complex<double> CpbarUFeFeAhPR(int gO1, int gI1) const;
   std::complex<double> CpbarUFeFehhPL(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFehhPR(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVPPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVPPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFeVZPR(int gO2, int gI2) const;
   std::complex<double> CpbarUFeFeVZPL(int gO1, int gI2) const;
   std::complex<double> CpbarUFeFvconjHpPL(int gO2, int gI2) const;
   double CpbarUFeFvconjHpPR(int, int) const;
   double CpbarUFeFvconjVWpPR(int, int) const;
   double CpbarUFeFvconjVWpPL(int gO1, int gI2) const;
   double CpbarFvFeHpPL(int, int) const;
   std::complex<double> CpbarFvFeHpPR(int gO1, int gI2) const;
   double CpbarFvFeVWpPR(int, int) const;
   std::complex<double> CpbarFvFeVWpPL(int gO1, int gI2) const;
   std::complex<double> CpbarFuFdHpPL(int gO2, int gI2) const;
   std::complex<double> CpbarFuFdHpPR(int gO1, int gI2) const;
   double CpbarFuFdVWpPR(int, int) const;
   std::complex<double> CpbarFuFdVWpPL(int gO1, int gI2) const;

   std::complex<double> self_energy_Hp_1loop(double p2 ) const;
   std::complex<double> self_energy_Ah_1loop(double p2 ) const;
   std::complex<double> self_energy_hh_1loop(double p2 ) const;
   std::complex<double> self_energy_VG_1loop(double p2 ) const;
   std::complex<double> self_energy_VP_1loop(double p2 ) const;
   std::complex<double> self_energy_VZ_1loop(double p2 ) const;
   std::complex<double> self_energy_VWp_1loop(double p2 ) const;
   std::complex<double> self_energy_Fd_1loop_1(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1(double p2) const;
   std::complex<double> self_energy_Fd_1loop_PR(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR(double p2) const;
   std::complex<double> self_energy_Fd_1loop_PL(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL(double p2) const;
   std::complex<double> self_energy_Fu_1loop_1(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PR(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PL(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL(double p2) const;
   std::complex<double> self_energy_Fe_1loop_1(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1(double p2) const;
   std::complex<double> self_energy_Fe_1loop_PR(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR(double p2) const;
   std::complex<double> self_energy_Fe_1loop_PL(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL(double p2) const;
   std::complex<double> self_energy_Fv_1loop_1(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_1(double p2) const;
   std::complex<double> self_energy_Fv_1loop_PR(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PR(double p2) const;
   std::complex<double> self_energy_Fv_1loop_PL(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fv_1loop_PL(double p2) const;
   std::complex<double> self_energy_Fd_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_1_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fd_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PR_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fd_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fd_1loop_PL_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fe_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_1_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fe_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PR_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fe_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fe_1loop_PL_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy_rotated(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy_rotated(double p2) const;
   std::complex<double> self_energy_Fu_1loop_1_heavy(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_1_heavy(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PR_heavy(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PR_heavy(double p2) const;
   std::complex<double> self_energy_Fu_1loop_PL_heavy(double p2 , int gO1, int gO2) const;
   Eigen::Matrix<std::complex<double>,3,3> self_energy_Fu_1loop_PL_heavy(double p2) const;

   std::complex<double> tadpole_hh_1loop() const;

   /// calculates the tadpoles at current loop order
   Eigen::Matrix<double, number_of_ewsb_equations, 1> tadpole_equations() const;

   /// calculates Higgs 2-loop self-energy
   double self_energy_hh_2loop(double p) const;
   /// calculates Higgs 3-loop self-energy
   double self_energy_hh_3loop() const;
   /// calculates Higgs 4-loop self-energy
   double self_energy_hh_4loop() const;

   void calculate_M2VG_pole();
   void calculate_MFv_pole();
   void calculate_M2hh_pole();
   void calculate_M2VP_pole();
   void calculate_M2VZ_pole();
   void calculate_MFd_pole();
   void calculate_MFu_pole();
   void calculate_MFe_pole();
   void calculate_M2VWp_pole();
   double calculate_M2VWp_pole(double);
   double calculate_M2VZ_pole(double);

   double calculate_MFv_DRbar(double, int) const;
   double calculate_MFe_DRbar(double, int) const;
   double calculate_MFu_DRbar(double, int) const;
   double calculate_MFd_DRbar(double, int) const;
   double calculate_M2VP_DRbar(double);
   double calculate_M2VZ_DRbar(double);
   double calculate_M2VWp_DRbar(double);
   double calculate_M2hh_DRbar(double);

   double ThetaW() const;

   double calculate_delta_alpha_em(double alphaEm) const;
   double calculate_delta_alpha_s(double alphaS) const;
   void calculate_Lambdax_DRbar();
   double calculate_theta_w(const softsusy::QedQcd&, double alpha_em_drbar);
   void calculate_Yu_DRbar(const softsusy::QedQcd&);
   void calculate_Yd_DRbar(const softsusy::QedQcd&);
   void calculate_Ye_DRbar(const softsusy::QedQcd&);
   double recalculate_m2w_pole(double);
   double max_rel_diff(const Standard_model& old) const;

protected:

   // Running parameters
   double g1{};
   double g2{};
   double g3{};
   double Lambdax{};
   Eigen::Matrix<double,3,3> Yu{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Yd{Eigen::Matrix<double,3,3>::Zero()};
   Eigen::Matrix<double,3,3> Ye{Eigen::Matrix<double,3,3>::Zero()};
   double mu2{};
   double v{};

private:

   static const int numberOfParameters = 33;

   struct Beta_traces {
      double traceYdAdjYd{};
      double traceYeAdjYe{};
      double traceYuAdjYu{};
      double traceYdAdjYdYdAdjYd{};
      double traceYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYu{};
      double traceYdAdjYuYuAdjYd{};
      double traceYdAdjYdYdAdjYdYdAdjYd{};
      double traceYdAdjYdYdAdjYuYuAdjYd{};
      double traceYdAdjYuYuAdjYdYdAdjYd{};
      double traceYdAdjYuYuAdjYuYuAdjYd{};
      double traceYeAdjYeYeAdjYeYeAdjYe{};
      double traceYuAdjYuYuAdjYuYuAdjYu{};
   };
   void calc_beta_traces(Beta_traces&) const;

   double calc_beta_g1_one_loop(const Beta_traces&) const;
   double calc_beta_g1_two_loop(const Beta_traces&) const;
   double calc_beta_g1_three_loop(const Beta_traces&) const;
   double calc_beta_g2_one_loop(const Beta_traces&) const;
   double calc_beta_g2_two_loop(const Beta_traces&) const;
   double calc_beta_g2_three_loop(const Beta_traces&) const;
   double calc_beta_g3_one_loop(const Beta_traces&) const;
   double calc_beta_g3_two_loop(const Beta_traces&) const;
   double calc_beta_g3_three_loop(const Beta_traces&) const;
   double calc_beta_g3_four_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_one_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_two_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_three_loop(const Beta_traces&) const;
   double calc_beta_Lambdax_four_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_three_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yu_four_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Yd_three_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_one_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_two_loop(const Beta_traces&) const;
   Eigen::Matrix<double,3,3> calc_beta_Ye_three_loop(const Beta_traces&) const;
   double calc_beta_mu2_one_loop(const Beta_traces&) const;
   double calc_beta_mu2_two_loop(const Beta_traces&) const;
   double calc_beta_mu2_three_loop(const Beta_traces&) const;
   double calc_beta_v_one_loop(const Beta_traces&) const;
   double calc_beta_v_two_loop(const Beta_traces&) const;
   double calc_beta_v_three_loop(const Beta_traces&) const;

   using EWSB_vector_t = Eigen::Matrix<double,number_of_ewsb_equations,1>;

   class EEWSBStepFailed : public Error {
   public:
      virtual ~EEWSBStepFailed() = default;
      virtual std::string what() const override { return "Could not perform EWSB step."; }
   };

   int ewsb_loop_order{4};
   int pole_mass_loop_order{4};
   bool force_output{false};      ///< switch to force output of pole masses
   double precision{1e-4};        ///< RG running precision
   double ewsb_iteration_precision{1e-5};
   Standard_model_physical physical{}; ///< contains the pole masses and mixings
   Problems problems{standard_model_info::model_name,
                     &standard_model_info::particle_names_getter,
                     &standard_model_info::parameter_names_getter};
   Loop_corrections loop_corrections{}; ///< used loop pole mass corrections
   Threshold_corrections threshold_corrections{}; ///< used low-energy threshold corrections
   Physical_input input{};

   int get_number_of_ewsb_iterations() const;
   int get_number_of_mass_iterations() const;
   int solve_ewsb_iteratively();
   int solve_ewsb_iteratively(int);
   int solve_ewsb_iteratively_with(EWSB_solver*, const Eigen::Matrix<double, number_of_ewsb_equations, 1>&);
   int solve_ewsb_tree_level_custom();
   EWSB_vector_t ewsb_initial_guess();
   EWSB_vector_t ewsb_step() const;
   void copy_DRbar_masses_to_pole_masses();

   void initial_guess_for_parameters(const softsusy::QedQcd&);
   bool check_convergence(const Standard_model& old) const;

   // Passarino-Veltman loop functions
   double A0(double) const noexcept;
   double B0(double, double, double) const noexcept;
   double B1(double, double, double) const noexcept;
   double B00(double, double, double) const noexcept;
   double B22(double, double, double) const noexcept;
   double H0(double, double, double) const noexcept;
   double F0(double, double, double) const noexcept;
   double G0(double, double, double) const noexcept;

   // DR-bar masses
   double M2VG{};
   double M2Hp{};
   Eigen::Array<double,3,1> MFv{Eigen::Array<double,3,1>::Zero()};
   double M2Ah{};
   double M2hh{};
   Eigen::Array<double,3,1> MFd{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFu{Eigen::Array<double,3,1>::Zero()};
   Eigen::Array<double,3,1> MFe{Eigen::Array<double,3,1>::Zero()};
   double M2VWp{};
   double M2VP{};
   double M2VZ{};

   // DR-bar mixing matrices
   Eigen::Matrix<std::complex<double>,3,3> Vd{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ud{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Vu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Uu{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ve{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<std::complex<double>,3,3> Ue{Eigen::Matrix<std::complex<double>,3,3>::Zero()};
   Eigen::Matrix<double,2,2> ZZ{Eigen::Matrix<double,2,2>::Zero()};


};

std::ostream& operator<<(std::ostream&, const Standard_model&);

} // namespace standard_model

} // namespace flexiblesusy

#endif
