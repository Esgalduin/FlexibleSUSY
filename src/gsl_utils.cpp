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

#include "gsl_utils.hpp"
#include "gsl_vector.hpp"
#include <cstddef>
#include <cmath>

namespace flexiblesusy {

/**
 * Returns true if GSL vector contains only finite elements (neither
 * nan nor inf), false otherwise.
 *
 * @param x GSL vector
 * @return true if vector contains only finite elements, false otherwise.
 */
bool is_finite(const gsl_vector* x)
{
   const std::size_t length = x->size;
   bool is_finite = true;

   for (std::size_t i = 0; i < length; i++)
      is_finite = is_finite && std::isfinite(gsl_vector_get(x, i));

   return is_finite;
}

/**
 * Returns true if GSL_vector contains only finite elements (neither
 * nan nor inf), false otherwise.
 *
 * @param x GSL vector
 * @return true if vector contains only finite elements, false otherwise.
 */
bool is_finite(const GSL_vector& v)
{
   const std::size_t length = v.size();
   bool finite = true;

   for (std::size_t i = 0; i < length; i++)
      finite = finite && std::isfinite(v[i]);

   return finite;
}

/**
 * Returns an Eigen array which contains the elements of the given GSL
 * vector.
 *
 * @param v GSL vector
 * @return Eigen array
 */
Eigen::ArrayXd to_eigen_array(const gsl_vector* v)
{
   return to_eigen_vector(v);
}

/**
 * Returns an Eigen array which contains the elements of the given GSL
 * vector.
 *
 * @param v GSL vector
 * @return Eigen array
 */
Eigen::ArrayXd to_eigen_array(const GSL_vector& v)
{
   return to_eigen_vector(v);
}

/**
 * Returns an Eigen array which contains the elements of the given GSL
 * vector.
 *
 * @param v GSL vector
 * @return Eigen vector
 */
Eigen::VectorXd to_eigen_vector(const gsl_vector* v)
{
   return to_eigen_vector(GSL_vector(v));
}

/**
 * Returns an Eigen array which contains the elements of the given GSL
 * vector.
 *
 * @param v GSL vector
 * @return Eigen vector
 */
Eigen::VectorXd to_eigen_vector(const GSL_vector& v)
{
   Eigen::VectorXd v2(v.size());

   for (std::size_t i = 0; i < v.size(); i++)
      v2(i) = v[i];

   return v2;
}

GSL_vector to_GSL_vector(const Eigen::VectorXd& v)
{
   GSL_vector v2(v.rows());

   for (std::size_t i = 0; i < v.rows(); i++)
      v2[i] = v(i);

   return v2;
}

GSL_vector to_GSL_vector(const gsl_vector* v)
{
   return GSL_vector(v);
}

/**
 * Copies values from an Eigen array to a GSL vector.
 *
 * @param src Eigen array
 * @param dst GSL vector
 */
void copy(const Eigen::ArrayXd& src, gsl_vector* dst)
{
   const std::size_t dim = src.rows();

   assert(dim == dst->size);

   for (std::size_t i = 0; i < dim; i++)
      gsl_vector_set(dst, i, src(i));
}

}
