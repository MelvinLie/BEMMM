// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_VECTORPOTENTIALSL_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_VECTORPOTENTIALSL_H_

namespace Bembel {
// forward declaration of class LaplaceVectorPotentialSL in order to define
// traits
template <typename LinOp>
class LaplaceVectorPotentialSL;

template <typename LinOp>
struct PotentialTraits<LaplaceVectorPotentialSL<LinOp>> {
  typedef Eigen::VectorXd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 3;
};

/**
 * \ingroup Laplace
 */
template <typename LinOp>
class LaplaceVectorPotentialSL
    : public PotentialBase<LaplaceVectorPotentialSL<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceVectorPotentialSL() {}

  Eigen::Matrix<double,3, 1>
  evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);

    // evaluate kernel
    auto kernel = evaluateKernel(point, y_f);

    // assemble surface curl
    Eigen::Vector2d curl_v = fun_ev.evaluateCurl(element, p);

    // Here we need to compute J * ( d2 , -d1 )^T * N_j * v_j
    Eigen::Matrix<double,3, 1> g;
    g(0,0) = y_f_dx(0) * curl_v(0) + y_f_dy(0) * curl_v(1);
    g(1,0) = y_f_dx(1) * curl_v(0) + y_f_dy(1) * curl_v(1);
    g(2,0) = y_f_dx(2) * curl_v(0) + y_f_dy(2) * curl_v(1);

    // integrand without basis functions
    Eigen::Matrix<double, 3, 1> integrand = kernel * g * ws ;

    return integrand;
  }

  Eigen::Matrix<double,3, 3>
  evaluateIntegrandDerivative_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);


    // evaluate kernel gradient
    Eigen::Matrix<double,1, 3> kernel_der = evaluateKernelGradient(point, y_f);

    // assemble surface curl
    Eigen::Vector2d curl_v = fun_ev.evaluateCurl(element, p);

    // Here we need to compute J * ( d2 , -d1 )^T * N_j * v_j
    Eigen::Matrix<double,1, 3> g;
    g(0) = y_f_dx(0) * curl_v(0) + y_f_dy(0) * curl_v(1);
    g(1) = y_f_dx(1) * curl_v(0) + y_f_dy(1) * curl_v(1);
    g(2) = y_f_dx(2) * curl_v(0) + y_f_dy(2) * curl_v(1);


    //We then find:
    Eigen::Matrix<double,3, 3> integrand;
    integrand.row(0) = kernel_der * g(0)  * ws ;
    integrand.row(1) = kernel_der * g(1)  * ws ;
    integrand.row(2) = kernel_der * g(2)  * ws ;

    return integrand;
  }


  Eigen::Matrix<std::complex<double>,3, Eigen::Dynamic>
  evaluateLocalHarmonicExpansion(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p,
                         const int L) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);


    //std::cout << "point = "<< point << std::endl;
    // approximate the inverse distance
    Eigen::VectorXcd approx_kernel= evaluateHarmExpKernel(point,y_f,L);


    // assemble surface curl
    Eigen::Vector2d curl_v = fun_ev.evaluateCurl(element, p);

    // Here we need to compute J * ( d2 , -d1 )^T * N_j * v_j
    Eigen::Matrix<double,3, 1> g;
    g(0,0) = y_f_dx(0) * curl_v(0) + y_f_dy(0) * curl_v(1);
    g(1,0) = y_f_dx(1) * curl_v(0) + y_f_dy(1) * curl_v(1);
    g(2,0) = y_f_dx(2) * curl_v(0) + y_f_dy(2) * curl_v(1);

    // integrand without basis functions
    Eigen::MatrixXcd integrand(3,(L+1)*(L+1));

    integrand.row(0) = approx_kernel * g(0) * ws;
    integrand.row(1) = approx_kernel * g(1) * ws;
    integrand.row(2) = approx_kernel * g(2) * ws;

    return integrand;
  }
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

   /**
    * \brief Kernel for the expansion of the inverse distance
    */
   Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> evaluateHarmExpKernel(const Eigen::Vector3d &x,
                         const Eigen::Vector3d &y, const int L) const {

     //Coordinates are shifted to x
     Eigen::Vector3d d = y - x;
     double r = d.norm();
     double theta = std::acos(d(2)/r);
     double phi = std::atan2(d(1),d(0));

     int num_coeffs = (L+1)*(L+1);

     //std::cout << "Num coeffs = "<<num_coeffs << std::endl;
     Eigen::Matrix<std::complex<double>,Eigen::Dynamic,1> retval;
     retval.resize(num_coeffs,1);

     //counter for coefficients
     int k = 0;
     for(int l = 0; l <= L; ++l){

       for (int m = -l; m <= l; ++m){

         retval(k,0) = std::conj(Ylm(l, m, theta, phi))/std::pow(r,l+1)/(2.*l + 1.);
         //complex conjugate
         //retval(k,0).imag(-1.*retval(k).imag());

         //increment
         ++k;
       }

     }

     return retval;

   }



};

}  // namespace Bembel
#endif
