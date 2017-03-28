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

#ifndef EIGEN_UTILS_H
#define EIGEN_UTILS_H

#include <Eigen/Core>
#include <cassert>
#include <complex>
#include <iomanip>
#include <sstream>
#include <string>
#include <limits>

namespace flexiblesusy {

/**
 * @class LRS_tensor
 * @brief collection of 3 NxN matrices
 */
template <int N, typename Scalar = std::complex<double> >
struct LRS_tensor {
public:
   using Matrix_t = Eigen::Matrix<Scalar,N,N>;
   using Tuple_t = std::tuple<Matrix_t,Matrix_t,Matrix_t>;

   LRS_tensor()
      : lrs(std::make_tuple(Matrix_t::Zero(),Matrix_t::Zero(),Matrix_t::Zero())) {}
   explicit LRS_tensor(const Tuple_t& t)
      : lrs(t) {}

   Matrix_t L() const { return std::get<0>(lrs); }
   Matrix_t R() const { return std::get<1>(lrs); }
   Matrix_t S() const { return std::get<2>(lrs); }
   Tuple_t tuple() const { return lrs; }
   static LRS_tensor PL() { return LRS_tensor(std::make_tuple(Matrix_t::Ones(),Matrix_t::Zero(),Matrix_t::Zero())); }
   static LRS_tensor PR() { return LRS_tensor(std::make_tuple(Matrix_t::Zero(),Matrix_t::Ones(),Matrix_t::Zero())); }
   static LRS_tensor PS() { return LRS_tensor(std::make_tuple(Matrix_t::Zero(),Matrix_t::Zero(),Matrix_t::Ones())); }

   LRS_tensor operator+(const LRS_tensor& rhs)
   {
      return LRS_tensor(std::make_tuple((L()+rhs.L()).eval(), (R()+rhs.R()).eval(), (S()+rhs.S()).eval()));
   }

   template <typename T>
   LRS_tensor operator*(T rhs)
   {
      return LRS_tensor(std::make_tuple((rhs*L()).eval(), (rhs*R()).eval(), (rhs*S()).eval()));
   }

   template <typename T>
   friend LRS_tensor operator*(T lhs, const LRS_tensor& rhs)
   {
      return LRS_tensor(std::make_tuple((lhs*rhs.L()).eval(), (lhs*rhs.R()).eval(), (lhs*rhs.S()).eval()));
   }

private:
   Tuple_t lrs; ///< L, R, S components
};

template <typename Derived>
int closest_index(double mass, const Eigen::ArrayBase<Derived>& v)
{
   int pos;
   typename Derived::PlainObject tmp;
   tmp.setConstant(mass);

   (v - tmp).abs().minCoeff(&pos);

   return pos;
}

template <class BinaryOp, class Derived>
Derived binary_map(
   const Eigen::DenseBase<Derived>& a, const Eigen::DenseBase<Derived>& b, BinaryOp op)
{
   typename Derived::PlainObject result(a.rows(), b.cols());

   assert(a.rows() == b.rows());
   assert(a.cols() == b.cols());

   for (int k = 0; k < a.cols(); k++)
      for (int i = 0; i < a.rows(); i++)
         result(i,k) = op(a(i,k), b(i,k));

   return result;
}

/**
 * Divides a by b element wise.  If the quotient is not finite, it is
 * set to zero.
 *
 * @param a numerator
 * @param b denominator
 */
template <class Derived>
Derived div_safe(
   const Eigen::DenseBase<Derived>& a, const Eigen::DenseBase<Derived>& b)
{
   using Scalar = typename Derived::Scalar;

   return binary_map(a, b, [](Scalar x, Scalar y){
         const Scalar q = x / y;
         return std::isfinite(q) ? q : Scalar{};
      });
}

/**
 * Calls eval() on any Eigen expression.  If the argument is not an
 * Eigen expression, the argument is returned.
 *
 * @param expr expression
 *
 * @return evaluated expression or argument
 */
template <class Derived>
auto Eval(const Eigen::DenseBase<Derived>& expr) -> decltype(expr.eval())
{
   return expr.eval();
}

inline char                 Eval(char                        expr) { return expr; }
inline short                Eval(short                       expr) { return expr; }
inline int                  Eval(int                         expr) { return expr; }
inline long                 Eval(long                        expr) { return expr; }
inline unsigned short       Eval(unsigned short              expr) { return expr; }
inline unsigned int         Eval(unsigned int                expr) { return expr; }
inline unsigned long        Eval(unsigned long               expr) { return expr; }
inline double               Eval(double                      expr) { return expr; }
inline std::complex<double> Eval(const std::complex<double>& expr) { return expr; }

/**
 * The element of v, which is closest to mass, is moved to the
 * position idx.
 *
 * @param idx new index of the mass eigenvalue
 * @param mass mass to compare against
 * @param v vector of masses
 * @param z corresponding mixing matrix
 */
template <typename DerivedArray, typename DerivedMatrix>
void move_goldstone_to(int idx, double mass, Eigen::ArrayBase<DerivedArray>& v,
                       Eigen::MatrixBase<DerivedMatrix>& z)
{
   int pos = closest_index(mass, v);
   if (pos == idx)
      return;

   const int sign = (idx - pos) < 0 ? -1 : 1;
   int steps = std::abs(idx - pos);

   // now we shuffle the states
   while (steps--) {
      const int new_pos = pos + sign;
      v.row(new_pos).swap(v.row(pos));
      z.row(new_pos).swap(z.row(pos));
      pos = new_pos;
   }
}

namespace {
template <class T>
struct Is_not_finite {
   bool operator()(T x) const noexcept { return !std::isfinite(x); }
};
} // anonymous namespace

/**
 * Returns all elements from src, which are not close to the elements
 * in cmp.  The returned vector will have the length (src.size() -
 * cmp.size()).
 *
 * @param src source vector
 * @param cmp vector with elements to compare against
 * @return vector with elements of src not close to cmp
 */
template<class Real, int Nsrc, int Ncmp>
Eigen::Array<Real,Nsrc - Ncmp,1> remove_if_equal(
   const Eigen::Array<Real,Nsrc,1>& src,
   const Eigen::Array<Real,Ncmp,1>& cmp)
{
   Eigen::Array<Real,Nsrc,1> non_equal(src);
   Eigen::Array<Real,Nsrc - Ncmp,1> dst;

   for (int i = 0; i < Ncmp; i++) {
      const int idx = closest_index(cmp(i), non_equal);
      non_equal(idx) = std::numeric_limits<double>::infinity();
   }

   std::remove_copy_if(non_equal.data(), non_equal.data() + Nsrc,
                       dst.data(), Is_not_finite<Real>());

   return dst;
}

/**
 * @brief reorders vector v according to ordering in vector v2
 * @param v vector with elementes to be reordered
 * @param v2 vector with reference ordering
 */
template<class Real, int N>
void reorder_vector(
   Eigen::Array<Real,N,1>& v,
   const Eigen::Array<Real,N,1>& v2)
{
   Eigen::PermutationMatrix<N> p;
   p.setIdentity();
   std::sort(p.indices().data(), p.indices().data() + p.indices().size(),
             [&v2] (int i, int j) { return v2[i] < v2[j]; });

#if EIGEN_VERSION_AT_LEAST(3,1,4)
   v.matrix().transpose() *= p.inverse();
#else
   Eigen::Map<Eigen::Matrix<Real,N,1> >(v.data()).transpose() *= p.inverse();
#endif
}

/**
 * @brief reorders vector v according to ordering of diagonal elements in mass_matrix
 * @param v vector with elementes to be reordered
 * @param matrix matrix with diagonal elements with reference ordering
 */
template<class Derived>
void reorder_vector(
   Eigen::Array<double,Eigen::MatrixBase<Derived>::RowsAtCompileTime,1>& v,
   const Eigen::MatrixBase<Derived>& matrix)
{
   reorder_vector(v, matrix.diagonal().array().eval());
}

template<class Derived>
std::string print_scientific(const Eigen::DenseBase<Derived>& v,
                             int number_of_digits = std::numeric_limits<typename Derived::Scalar>::digits10 + 1)
{
   std::ostringstream sstr;

   for (int k = 0; k < v.cols(); k++) {
      for (int i = 0; i < v.rows(); i++) {
         sstr << std::setprecision(number_of_digits)
              << std::scientific << v(i,k) << ' ';
      }
      sstr << '\n';
   }

   return sstr.str();
}

} // namespace flexiblesusy

#endif
