// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYERPOTENTIAL_H_
#define BEMBEL_LINEAROPERATOR_LAPLACE_LAPLACEDOUBLELAYERPOTENTIAL_H_

namespace Bembel {
// forward declaration of class LaplaceDoubleLayerPotential in order to define
// traits
template <typename LinOp>
class LaplaceDoubleLayerPotential;

template <typename LinOp>
struct PotentialTraits<LaplaceDoubleLayerPotential<LinOp>> {
  typedef Eigen::VectorXd::Scalar Scalar;
  static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup Laplace
 */
template <typename LinOp>
class LaplaceDoubleLayerPotential
    : public PotentialBase<LaplaceDoubleLayerPotential<LinOp>, LinOp> {
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  LaplaceDoubleLayerPotential() {}
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
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);

    // compute normal vector direction
    auto n_y = y_f_dx.cross(y_f_dy);

    // compute surface measures from tangential derivatives
    auto y_kappa = n_y.norm();

    // scale normal vector
    n_y /= y_kappa;

    // evaluate kernel
    auto kernel = evaluateKernel(point, y_f, n_y);

    // assemble Galerkin solution
    auto cauchy_value = fun_ev.evaluate(element, p);

    // integrand without basis functions
    auto integrand = kernel * cauchy_value * y_kappa * ws * element.get_h() * element.get_h();

    return integrand;
  }

  void evaluateIntegrand_impl(const SuperSpace<LinOp> &super_space,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, 1> *intval) const {
    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);

    // compute normal vector direction
    auto n_y = y_f_dx.cross(y_f_dy);

    // compute surface measures from tangential derivatives
    auto y_kappa = n_y.norm();

    // scale normal vector
    n_y /= y_kappa;

    // evaluate kernel
    auto kernel = evaluateKernel(point, y_f, n_y);

    //evaluate basis functions
    super_space.addScaledBasis(intval , kernel * y_kappa * ws * element.get_h(), s);//

    //super_space.addScaledBasis(&basis_eval , kappa * ws * element.get_h() , s);
  }

  /*
  Evaluates the matrix entries to compute the derivative of the magnetic double
  layer potential, at position point, for all basis functions on the element
  */
  void evaluateIntegrandDerivative_mat(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // evaluate kernel gradient
    Eigen::Vector3d kernel_der = evaluateKernelGradient(point, x_f);


    // compute the surface curl for all the basis function on this element
    super_space.addScaledVecTimesCurlBasisInteraction(intval, kernel_der, ws,s, x_f_dx,x_f_dy);

  }

  /*
  Evaluates the matrix entries to compute the derivative of the magnetic double
  layer potential, at position point, for all basis functions on the element
  */
  void evaluateIntegrandDerivative_far_field(const SuperSpace<LinOp> &super_space,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, 3> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto y_f = p.segment<3>(3);
    auto y_f_dx = p.segment<3>(6);
    auto y_f_dy = p.segment<3>(9);

    // compute normal vector direction
    auto n_y = y_f_dx.cross(y_f_dy);

    // compute surface measures from tangential derivatives
    auto y_kappa = n_y.norm();

    // scale normal vector
    n_y /= y_kappa;

    // evaluate kernel
    Eigen::Vector3d kernel = evaluateGradKernel(point, y_f, n_y);


    //evaluate basis functions
    Eigen::Matrix<double,-1,1> tmp_val(I2,1);
    tmp_val.setZero();
    super_space.addScaledBasis(&tmp_val , y_kappa * ws * element.get_h(), s);//

    //increment
    intval->block(0,0,I2,1) +=  kernel(0) * tmp_val;
    intval->block(0,1,I2,1) +=  kernel(1) * tmp_val;
    intval->block(0,2,I2,1) +=  kernel(2) * tmp_val;
    //super_space.addScaledBasis(&basis_eval , kappa * ws * element.get_h() , s);

  }


  /*
  Evaluates the matrix entries to compute the second derivative of the magnetic double
  layer potential, at position point, for all basis functions on the element
  */
  void evaluateIntegrandSecondDerivative_mat(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // evaluate kernel gradient
    Eigen::Matrix3d kernel_der = evaluateKernelDerivative2(point, x_f);

    Eigen::MatrixXd tmp_val(3,I2);


    // compute the surface curl for all the basis function on this element
    tmp_val.setZero();
    super_space.addScaledVecTimesCurlBasisInteraction(&tmp_val, kernel_der.row(0), ws,s, x_f_dx,x_f_dy);
    intval->block(0,0,3,I2) += tmp_val;

    tmp_val.setZero();
    super_space.addScaledVecTimesCurlBasisInteraction(&tmp_val, kernel_der.row(1), ws,s, x_f_dx,x_f_dy);
    intval->block(3,0,3,I2) += tmp_val;

    tmp_val.setZero();
    super_space.addScaledVecTimesCurlBasisInteraction(&tmp_val, kernel_der.row(2), ws,s, x_f_dx,x_f_dy);
    intval->block(6,0,3,I2) += tmp_val;

  }

  /*

  */
  void evaluateFMMMoments(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> *intval,
                         const int L = 5  ) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    //int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);


    // compute normal vector
    Eigen::Vector3d n = x_f_dx.cross(x_f_dy);

    // compute surface measure from tangential derivatives
    double kappa = n.norm();

    n /= kappa;

    //compute derivatives of solid Harmonics
    Eigen::MatrixXcd R_p = Rlm_p_alt(L, x_f, point);
    Eigen::VectorXcd nR_p = n(0)*R_p.col(0) + n(1)*R_p.col(1) + n(2)*R_p.col(2);

    //evaluate basis functions
    Eigen::Matrix<double,-1,1> basis_eval;
    basis_eval.resize((P+1)*(P+1),1);
    basis_eval.setZero();
    super_space.addScaledBasis(&basis_eval , kappa * ws * element.get_h() , s);

    (*intval) += nR_p.conjugate() * basis_eval.transpose();

    return;
  }

  void evaluateFMMMoments_surface_currents(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> *intval,
                         const int L = 5  ) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    //int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);


    // compute normal vector
    Eigen::Vector3d n = x_f_dx.cross(x_f_dy);

    // compute surface measure from tangential derivatives
    double kappa = n.norm();

    n /= kappa;

    //compute derivatives of solid Harmonics
    Eigen::MatrixXcd R_p = Rlm_p_alt(L, x_f, point);
    Eigen::VectorXcd nR_p = n(0)*R_p.col(0) + n(1)*R_p.col(1) + n(2)*R_p.col(2);

    //evaluate basis functions
    Eigen::Matrix<double,-1,1> basis_eval;
    basis_eval.resize((P+1)*(P+1),1);
    basis_eval.setZero();
    super_space.addScaledBasis(&basis_eval , kappa * ws * element.get_h() , s);

    //return outer product

    //std::cout<< "Intval = " << (*intval).rows() << " x " << (*intval).cols() << std::endl;
    //std::cout<< "nR_p = " << nR_p.rows() << " x " << nR_p.cols() << std::endl;
    //std::cout<< "basis_eval = " << basis_eval.rows() << " x " << basis_eval.cols() << std::endl;
    //Eigen::MatrixXcd prod =    nR_p * basis_eval.transpose();
    //std::cout<< "prod = " << prod.rows() << " x " << prod.cols() << std::endl;

    (*intval) += nR_p.conjugate() * basis_eval.transpose();

    return;
  }

  /*
  Evaluates the matrix entries to compute the derivative of the magnetic double
  layer potential, at position point, for one basis function on the element
  */
  void evaluateIntegrandDerivative_single_basis(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval, int k) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // evaluate kernel gradient
    Eigen::Vector3d kernel_der = evaluateKernelGradient(point, x_f);


    // compute the surface curl for one basis function on this element
    super_space.addScaledVecTimesCurlBasisInteraction_single_basis(intval, kernel_der, ws,s, x_f_dx,x_f_dy,k);

  }

  /*
  Evaluates the matrix entries to compute the derivative of the magnetic double
  layer potential, at position point, for one a linear combination of Bernstein
  Polynomials
  */
  void evaluateIntegrandDerivative_BSpline(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval, double *weights, int *brnst_fcns, int num_fcns) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    // get quadrature weights
    auto ws = p(2);

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // evaluate kernel gradient
    Eigen::Vector3d kernel_der = evaluateKernelGradient(point, x_f);


    // compute the surface curl for one basis function on this element
    super_space.addScaledVecTimesCurlBasisInteraction_BSpline(intval, kernel_der, ws,s, x_f_dx,x_f_dy,weights,brnst_fcns, num_fcns);

  }

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


    // evaluate kernel gradient
    Eigen::Matrix<double,1, 3> kernel_der = evaluateKernelGradient(point, x_f);

    // assemble surface curl
    Eigen::Vector2d curl_v = fun_ev.evaluateCurl(element, p);

    // Here we need to compute J * ( d2 , -d1 )^T * N_j * v_j
    Eigen::Matrix<double,1, 3> g;
    g(0) = x_f_dx(0) * curl_v(0) + x_f_dy(0) * curl_v(1);
    g(1) = x_f_dx(1) * curl_v(0) + x_f_dy(1) * curl_v(1);
    g(2) = x_f_dx(2) * curl_v(0) + x_f_dy(2) * curl_v(1);


    //We then find:
    Eigen::Matrix<double,1, 3> integrand;
    integrand(0) = (kernel_der(1) * g(2) - kernel_der(2) * g(1)) * ws;
    integrand(1) = (kernel_der(2) * g(0) - kernel_der(0) * g(2)) * ws;
    integrand(2) = (kernel_der(0) * g(1) - kernel_der(1) * g(0)) * ws;

    return integrand;
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
    /*
    c.x = x.x - y.x;
    c.y = x.y - y.y;
    c.z = x.z - y.z;

    return ((c.x * n_y.x + c.y * n_y.y + c.z * n_y.z) / (4 * BEMBEL_PI *
                                  pow(c.x * c.x + c.y * c.y + c.z * c.z, 1.5)));

    */
  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
  Eigen::Vector3d evaluateGradKernel(const Eigen::Vector3d &x,
                        const Eigen::Vector3d &y,
                        const Eigen::Vector3d &n_y) const {

    Eigen::Vector3d c;

    c(0) = x(0) - y(0);
    c(1) = x(1) - y(1);
    c(2) = x(2) - y(2);

    double dist = c.norm();
    double proj = n_y.transpose()*c;

    Eigen::Vector3d ret_vec = n_y/dist/dist/dist;

    ret_vec -= 3.*proj*c/dist/dist/dist/dist/dist;


    return ret_vec / 4. / BEMBEL_PI;

  }

  /**
   * \brief Fundamental solution of Laplace problem
   */
   double evaluateSingleLayerKernel(const Eigen::Vector3d &x,
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

   Eigen::Matrix<double,3, 3> evaluateKernelDerivative2(const Eigen::Vector3d &x,
                         const Eigen::Vector3d &y) const {

     Eigen::Matrix<double,3, 3> v_ret;
     Eigen::VectorXd d = x - y;
     double dist = d.norm();
     double dist_sq = dist*dist;

     v_ret(0,0) = 1. - 3.*d(0)*d(0)/dist_sq;
     v_ret(0,1) = -3.*d(0)*d(1)/dist_sq;
     v_ret(0,2) = -3.*d(0)*d(2)/dist_sq;

     v_ret(1,0) = v_ret(0,1);
     v_ret(1,1) = 1. - 3.*d(1)*d(1)/dist_sq;
     v_ret(1,2) = -3.*d(1)*d(2)/dist_sq;

     v_ret(2,0) = v_ret(0,2);
     v_ret(2,1) = v_ret(1,2);
     v_ret(2,2) = 1. - 3.*d(2)*d(2)/dist_sq;

     return -1.*v_ret / 4. / BEMBEL_PI / dist_sq/dist;
   }


};

}  // namespace Bembel
#endif
