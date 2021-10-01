// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_POTENTIAL_DISCRETEMVP_H_
#define BEMBEL_POTENTIAL_DISCRETEMVP_H_

namespace Bembel {
/**
 *  \ingroup Potential
 *  \brief DiscretePotential
 *  \todo  Add a documentation
 */
class DiscreteMVP {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  DiscreteMVP() {}
  DiscreteMVP(const AnsatzSpace<LaplaceHypersingularOperator> &ansatz_space) {
    init_DiscretePotential(ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscretePotential
  //////////////////////////////////////////////////////////////////////////////

  void init_DiscretePotential(const AnsatzSpace<LaplaceHypersingularOperator> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    fun_ev_ = FunctionEvaluator<LaplaceHypersingularOperator>(ansatz_space_);
    /**
     * \todo obtain this from ansatz space
     */
    deg_ = ansatz_space_.get_polynomial_degree() + 1;
    return;
  }


  //////////////////////////////////////////////////////////////////////////////
  //    compute
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<double,Eigen::Dynamic, 3>
  evaluate(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    auto FunctionSpaceVectorDimension = 1;
    auto OutputDimension = 3;

    GaussSquare<Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];

    auto super_space = ansatz_space_.get_superspace();

    auto element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();

    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    Eigen::Matrix<double,Eigen::Dynamic,3> potential;
    potential.resize(points.rows(),3);
    potential.setZero();

#pragma omp parallel
    {
      Eigen::Matrix<double,Eigen::Dynamic,3>  my_potential;
      my_potential.resize(points.rows(),3);
      my_potential.setZero();
      for (auto element = element_tree.cpbegin();
           element != element_tree.cpend(); ++element) {
#pragma omp single nowait
        {
          SurfacePoint qp;
          for (auto j = 0; j < Q.w_.size(); ++j) {
            super_space.map2surface(
                *element, Q.xi_.col(j),
                element->get_h() * element->get_h() * Q.w_(j), &qp);
            for (auto i = 0; i < points.rows(); ++i) {
              my_potential.row(i) += pot_.evaluateIntegrand_impl(
                  fun_ev_, *element, points.row(i), qp);
            }
          }
        }
      }
#pragma omp critical
      potential += my_potential;
    }
    return potential;
  }


  //////////////////////////////////////////////////////////////////////////////
  //    Evaluate Gradient of Potential
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<double, Eigen::Dynamic, 3>
  evaluate_der(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
    auto FunctionSpaceVectorDimension = 1;

    GaussSquare<Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];

    auto super_space = ansatz_space_.get_superspace();

    auto element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();

    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    int pot_dim = 3;

    Eigen::MatrixXd  potential_der;

    potential_der.resize(points.rows()*pot_dim,3);

    potential_der.setZero();

  #pragma omp parallel
    {
      Eigen::Matrix<double,Eigen::Dynamic,3> my_derivative;

      my_derivative.resize(points.rows()*pot_dim,3);
      my_derivative.setZero();
      for (auto element = element_tree.cpbegin();
           element != element_tree.cpend(); ++element) {
  #pragma omp single nowait
        {
          SurfacePoint qp;
          for (auto j = 0; j < Q.w_.size(); ++j) {
            super_space.map2surface(
                *element, Q.xi_.col(j),
                element->get_h() * element->get_h() * Q.w_(j), &qp);
            for (auto i = 0; i < points.rows(); ++i) {
              my_derivative.block(i*pot_dim,0,pot_dim,3) += pot_.evaluateIntegrandDerivative_impl(fun_ev_, *element, points.row(i), qp);
            }
          }
        }
      }
  #pragma omp critical
      potential_der += my_derivative;
    }
    return potential_der;
  }

  //////////////////////////////////////////////////////////////////////////////
  //    Compute local harmonic expansion
  //////////////////////////////////////////////////////////////////////////////
  Eigen::Matrix<std::complex<double>, Eigen::Dynamic,  Eigen::Dynamic>
  compute_harmonic_expansion(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points, const int L) {

    GaussSquare<Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];

    auto super_space = ansatz_space_.get_superspace();

    auto element_tree = super_space.get_mesh().get_element_tree();
    auto number_of_elements = element_tree.get_number_of_elements();

    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);

    //number of expansion coefficients
    int num_coeffs = (L+1)*(L+1);

    //return matrix:
    // we have 3 rows for each point. In each of these rows, we store the coefficients for one component:
    // M = [c_lm_x^T,
    //      c_lm_y^T,
    //      c_lm_z^T]
    Eigen::MatrixXcd harmonic_expansion;

    harmonic_expansion.resize(points.rows()*3,num_coeffs);


    harmonic_expansion.setZero();

  #pragma omp parallel
    {
      Eigen::MatrixXcd my_expansion;

      my_expansion.resize(points.rows()*3,num_coeffs);
      my_expansion.setZero();
      for (auto element = element_tree.cpbegin();
           element != element_tree.cpend(); ++element) {
  #pragma omp single nowait
        {
          SurfacePoint qp;
          for (auto j = 0; j < Q.w_.size(); ++j) {
            super_space.map2surface(
                *element, Q.xi_.col(j),
                element->get_h() * element->get_h() * Q.w_(j), &qp);
            for (auto i = 0; i < points.rows(); ++i) {
              my_expansion.block(i*3,0,3,num_coeffs) += pot_.evaluateLocalHarmonicExpansion(fun_ev_, *element, points.row(i), qp, L);
            }
          }
        }
      }
  #pragma omp critical
      harmonic_expansion += my_expansion;
    }
    return harmonic_expansion;
  }


  Eigen::Matrix<double, Eigen::Dynamic,  3> evaluate_harmonic_expansion(const Eigen::Matrix<double,1,3> &ref_pos, const Eigen::MatrixXd &points, const Eigen::MatrixXcd &coeffs) {

    int num_coeffs = coeffs.cols();
    int L = std::sqrt(num_coeffs) - 1;

    int num_particles = points.rows();

    //Eigen::MatrixXd ret_mat;
    Eigen::MatrixXcd tmp_mat;

    tmp_mat.resize(num_particles,3);
    tmp_mat.setZero();

    Eigen::Vector3d d;
    double theta,phi,r,r_pow;

    std::complex<double> sp_harm;

    for(int i = 0; i < num_particles; ++i){

      int k = 0;

      d = points.row(i)-ref_pos.row(0);

      r = d.norm();


      if (r > 0.){
        theta = std::acos(d(2)/r);
        phi = std::atan2(d(1),d(0));
      }
      else{
        theta = 0.;
        phi = 0.;
      }

      r_pow = 1.;

      for(int l = 0; l <= L; ++l){
        for(int m = -l ; m <= l ; ++m){

          sp_harm = Ylm(l, m, theta, phi);

          tmp_mat(i,0) += r_pow*coeffs(0,k)*sp_harm;
          tmp_mat(i,1) += r_pow*coeffs(1,k)*sp_harm;
          tmp_mat(i,2) += r_pow*coeffs(2,k)*sp_harm;

          ++k;

        }
        r_pow *= r;
      }

    }
    return tmp_mat.real();

  }

  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_cauchy_data(
      const Eigen::Matrix<double,Eigen::Dynamic, 1> &cauchy_data) {
    fun_ev_.set_function(cauchy_data);
  }
  void set_degree(const int &deg) { deg_ = deg; }



  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////

 private:
  int deg_;
  LaplaceVectorPotentialSL<LaplaceHypersingularOperator> pot_;
  AnsatzSpace<LaplaceHypersingularOperator> ansatz_space_;
  FunctionEvaluator<LaplaceHypersingularOperator> fun_ev_;



};  // namespace Bembel

}  // namespace Bembel
#endif
