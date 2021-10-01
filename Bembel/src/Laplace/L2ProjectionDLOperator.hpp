// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_L2PROJECTIONDLLELAYEROPERATOR_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_L2PROJECTIONDLLAYEROPERATOR_H_

namespace Bembel {
// forward declaration of class LaplaceDoubleLayerOperator in order to define
// traits
class L2ProjectionDLOperator;

template <>
struct LinearOperatorTraits<L2ProjectionDLOperator> {
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
class L2ProjectionDLOperator
    : public LinearOperatorBase<L2ProjectionDLOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  L2ProjectionDLOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<L2ProjectionDLOperator>::Scalar,
          Eigen::Dynamic, Eigen::Dynamic> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);



    return;
  }

  Eigen::Matrix<double, 1, 1> evaluateFMMInterpolation_impl(
      const SurfacePoint &p1, const SurfacePoint &p2) const {

    // interpolation
    Eigen::Matrix<double, 1, 1> intval;
    intval(0) = 0;

    return intval;
  }

};

}  // namespace Bembel
#endif
