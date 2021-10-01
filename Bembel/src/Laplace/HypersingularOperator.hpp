// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEHYPERSINGULARLAYEROPERATOR_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEHYPERSINGULARLAYEROPERATOR_H_


namespace Bembel {
// forward declaration of class LaplaceDoubleLayerOperator in order to define
// traits
class LaplaceHypersingularOperator;

template <>
struct LinearOperatorTraits<LaplaceHypersingularOperator> {
  typedef Eigen::VectorXd EigenType;
  typedef Eigen::VectorXd::Scalar Scalar;
  enum {
    OperatorOrder = 1,   //H^{1/2} -> H^{-1/2}
    Form = DifferentialForm::Continuous,
    NumberOfFMMComponents = 2
  };
};

/**
 * \ingroup Laplace
 */
class LaplaceHypersingularOperator
    : public LinearOperatorBase<LaplaceHypersingularOperator> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceHypersingularOperator() {}
  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p1, const SurfacePoint &p2,
      Eigen::Matrix<
          typename LinearOperatorTraits<LaplaceHypersingularOperator>::Scalar,
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


    // integrand without basis functions
    auto integrand = evaluateKernel(x_f, y_f) * ws * wt;

    super_space.addScaledCurlBasisInteraction(intval, integrand, s, t,
                                              x_f_dx,x_f_dy,
                                              y_f_dx,y_f_dy);

    return;
  }

  Eigen::Matrix<double, 2, 2> evaluateFMMInterpolation_impl(
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

    // compute surface measures from tangential derivatives
    //auto x_kappa = x_f_dx.cross(x_f_dy).norm();
    //auto y_kappa = y_f_dx.cross(y_f_dy).norm();

    double kernel = evaluateKernel(x_f, y_f);// * x_kappa * y_kappa;

    // interpolation
    Eigen::Matrix<double, 2, 2> intval;
    intval.setZero();

    //An interpolation-based fast multipole method...
    //page 18
    intval(0, 0) =  kernel * x_f_dy.dot(y_f_dy);
    intval(0, 1) =  -kernel * x_f_dy.dot(y_f_dx);
    intval(1, 0) =  -kernel * x_f_dx.dot(y_f_dy);
    intval(1, 1) =  kernel * x_f_dx.dot(y_f_dx);


    //std::cout << "evaluateFMMInterpolation_impl"<< std::endl;

    return intval;
  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
   double evaluateKernel(const Eigen::Vector3d &x,
                         const Eigen::Vector3d &y) const {
     return 1. / 4. / BEMBEL_PI / (x - y).norm();
   }



};
/**
 * \brief The Laplace hypersingular operator requires a special treatment of the
 * moment matrices of the FMM due to the involved derivatives on the ansatz
 * functions.
 */
namespace H2Multipole {

template <typename InterpolationPoints>
struct Moment2D<InterpolationPoints, LaplaceHypersingularOperator> {
  static std::vector<Eigen::MatrixXd> compute2DMoment(
      const SuperSpace<LaplaceHypersingularOperator> &super_space,
      const int cluster_level, const int cluster_refinements,
      const int number_of_points) {


    Eigen::MatrixXd moment_dx = moment2DComputer<
        Moment1DDerivative<InterpolationPoints, LaplaceHypersingularOperator>,
        Moment1D<InterpolationPoints, LaplaceHypersingularOperator>>(
        super_space, cluster_level, cluster_refinements, number_of_points);

   Eigen::MatrixXd moment_dy = moment2DComputer<
            Moment1D<InterpolationPoints, LaplaceHypersingularOperator>,
            Moment1DDerivative<InterpolationPoints, LaplaceHypersingularOperator>>(
            super_space, cluster_level, cluster_refinements, number_of_points);


    std::vector<Eigen::MatrixXd> vector_of_moments;

    Eigen::MatrixXd moment1(moment_dx.rows() + moment_dy.rows(), moment_dx.cols());
    moment1 << moment_dx, moment_dy;


    int refine_lvl = super_space.get_refinement_level();
    int polynomial_degree = super_space.get_polynomial_degree();

    //This factor must come from the derivative we are taking in the moment matrices
    double fac = std::pow(2.,refine_lvl);

    moment1 /= fac;
    vector_of_moments.push_back(moment1);


    return vector_of_moments;

  }


};

}

}  // namespace Bembel

#endif
