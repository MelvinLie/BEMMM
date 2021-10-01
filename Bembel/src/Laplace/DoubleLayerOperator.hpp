// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYEROPERATOR_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYEROPERATOR_H_

namespace Bembel {
// forward declaration of class LaplaceDoubleLayerOperator in order to define
// traits
class LaplaceDoubleLayerOperator;

template <>
struct LinearOperatorTraits<LaplaceDoubleLayerOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = -2,
    //Form = DifferentialForm::Continuous,
    Form = DifferentialForm::Discontinuous, //tested with X_1
    NumberOfFMMComponents = 1
  };
};

/**
 * \ingroup Laplace
 */
class LaplaceDoubleLayerOperator
    : public LinearOperatorBase<LaplaceDoubleLayerOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceDoubleLayerOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<LaplaceDoubleLayerOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get quadrature weights
    auto ws = p1(2);
    auto wt = p2(2);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute normal vectors direction
    //auto n_x = x_f_dx.cross(x_f_dy);
    auto n_y = y_f_dx.cross(y_f_dy);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();
    auto y_kappa = n_y.norm();

    // scale normal vector
    n_y /= y_kappa;

    // integrand without basis functions
    auto integrand = evaluateKernel(x_f, y_f, n_y) * x_kappa * y_kappa * ws * wt;

    //std::cout << "Integrand = " << integrand << std::endl;

    // multiply basis functions with integrand and add to intval, this is an
    // efficient implementation of
    //(*intval) += super_space.basisInteraction(s, t) * evaluateKernel(x_f, y_f)
    //* x_kappa * y_kappa * ws * wt;
    super_space.addScaledBasisInteraction(intval, integrand, s, t);

    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {
    // get evaluation points on unit square
    auto s = p1.segment<2>(0);
    auto t = p2.segment<2>(0);

    // get points on geometry and tangential derivatives
    auto x_f = p1.segment<3>(3);
    auto x_f_dx = p1.segment<3>(6);
    auto x_f_dy = p1.segment<3>(9);
    auto y_f = p2.segment<3>(3);
    auto y_f_dx = p2.segment<3>(6);
    auto y_f_dy = p2.segment<3>(9);

    // compute normal vectors direction
    //auto n_x = x_f_dx.cross(x_f_dy);
    auto n_y = y_f_dx.cross(y_f_dy);

    // compute surface measures from tangential derivatives
    auto x_kappa = x_f_dx.cross(x_f_dy).norm();
    auto y_kappa = n_y.norm();

    // scale normal vector
    n_y /= y_kappa;

    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = evaluateKernel(x_f, y_f, n_y) * x_kappa * y_kappa;

    return intval;
  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
  double evaluateKernel(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y,
                        const Eigen::Vector3d &n_y) const {
    Eigen::Vector3d c;

    c(0) = x(0) - y(0);
    c(1) = x(1) - y(1);
    c(2) = x(2) - y(2);

    return ((c(0) * n_y(0) + c(1) * n_y(1) + c(2) * n_y(2)) / (4 * BEMBEL_PI *
                  pow(c(0) * c(0) + c(1) * c(1) + c(2) * c(2), 1.5)));

  }
};

}  // namespace Bembel
#endif
