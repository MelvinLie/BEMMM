// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACESINGLELAYERPOTENTIAL_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACESINGLELAYERPOTENTIAL_H_

namespace Bembel {
// forward declaration of class LaplaceSingleLayerPotential in order to define
// traits
template <typename LinOp>
class LaplaceSingleLayerPotential;

template <typename LinOp>
struct PotentialTraits<LaplaceSingleLayerPotential<LinOp>> {
  typedef Eigen::VectorXd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup Laplace
 */
template <typename LinOp>
class LaplaceSingleLayerPotential
    : public PotentialBase<LaplaceSingleLayerPotential<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceSingleLayerPotential() {}
  Eigen::Matrix<
      typename PotentialReturnScalar<
          typename LinearOperatorTraits<LinOp>::Scalar, double>::Scalar,
      1, 1>
  evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();

    // evaluate kernel
    auto kernel = evaluateKernel(point, x_f);

    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);

    // integrand without basis functions
    auto integrand = kernel * cauchy_value * x_kappa * ws;

    return integrand;
  }

  //////////////////////////////////////////////////////////////////////////////
  //    Add for post-processing of magnetic measurement data
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<double,1, 3>
  evaluateIntegrandDerivative_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();

    // evaluate kernel
    auto kernel = evaluateKernelGradient(point, x_f);

    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);


    // integrand without basis functions
    auto integrand = kernel * cauchy_value(0) * x_kappa * ws;

    return integrand;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    Add for post-processing of magnetic measurement data
  //////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Fundamental solution of Laplace problem
   */
  double evaluateKernel(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y) const {
    return 1. / 4. / BEMBEL_PI / (x - y).norm();
  }

  /**
   * \brief Gradient of fundamental solution of Laplace problem
   */
  Eigen::Matrix<double,1, 3> evaluateKernelGradient(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y) const {

    Eigen::Matrix<double,1, 3> v_ret;
    v_ret(0,0) = y[0] - x[0];
    v_ret(0,1) = y[1] - x[1];
    v_ret(0,2) = y[2] - x[2];

    return v_ret / 4. / BEMBEL_PI / pow((x - y).norm(), 3);
  }

};

}  // namespace Bembel
#endif
