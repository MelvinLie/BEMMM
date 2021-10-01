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
  DiscretePotential() {}
  DiscretePotential(const AnsatzSpace<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> &ansatz_space) {
    init_DiscretePotential(ansatz_space);
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscretePotential
  //////////////////////////////////////////////////////////////////////////////
  void init_DiscretePotential(const AnsatzSpace<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    fun_ev_ = FunctionEvaluator<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>>(ansatz_space_);
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
  //    Compute Gradient of Potential Matrix
  // here we iterate over the elements and call the function
  // evaluateIntegrandDerivative_mat
  //////////////////////////////////////////////////////////////////////////////
  void compute_grad_evaluation_matrix(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dx,
                                    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dy,
                                  Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dz) {

      //number of measurements
      int M = points.rows();

      //projector matrix
      auto projector = ansatz_space_.get_transformation_matrix();


      //number of basis functions
      int N_large = projector.rows();
      int N = ansatz_space_.get_number_of_dofs();

      Dx = Eigen::MatrixXd(M,N_large);
      Dy = Eigen::MatrixXd(M,N_large);
      Dz = Eigen::MatrixXd(M,N_large);

      GaussSquare<Constants::maximum_quadrature_degree> GS;



      auto super_space = ansatz_space_.get_superspace();
      auto element_tree = super_space.get_mesh().get_element_tree();

      auto number_of_elements = element_tree.get_number_of_elements();

      auto polynomial_degree = super_space.get_polynomial_degree();
      auto polynomial_degree_plus_one_squared =
          (polynomial_degree + 1) * (polynomial_degree + 1);

      auto Q = GS[deg_];

#pragma omp parallel for
for (int m = 0; m < M; ++m){

      //iterate over all elements
      for (auto element = element_tree.cpbegin();
               element != element_tree.cpend(); ++element){

              //point on surface for integratio
              SurfacePoint qp;

              //temporal storage for matrix blocks
              Eigen::MatrixXd tmp_D;
              tmp_D.resize(3,polynomial_degree_plus_one_squared);

              tmp_D.setZero();

             //Gaussian integtation
             for (auto j = 0; j < Q.w_.size(); ++j) {
                super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                pot_.evaluateIntegrandDerivative_mat(super_space, *element, points.row(m),qp, &tmp_D);

              }
                Dx.block(m,polynomial_degree_plus_one_squared * element->id_,
                               1,polynomial_degree_plus_one_squared) = tmp_D.row(0);

                Dy.block(m,polynomial_degree_plus_one_squared * element->id_,
                               1,polynomial_degree_plus_one_squared) = tmp_D.row(1);

                Dz.block(m,polynomial_degree_plus_one_squared * element->id_,
                               1,polynomial_degree_plus_one_squared) = tmp_D.row(2);

            }
        }



    Dx = Dx * projector;
    Dy = Dy * projector;
    Dz = Dz * projector;

    return;
  }

  //this function imposes a vanishing mean gauge to the matrices Dx,Dy and Dz
  Eigen::VectorXd impose_zero_mean_gauge(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dx,
                              Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dy,
                              Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &Dz){

      //number of unknowns
      int N = Dx.cols();

      //number of evaluations
      int M = Dx.rows();

      //function for stability vector
      std::function<double(Eigen::Vector3d)> unit_function = [](Eigen::Vector3d in){
        return 1.;
      };

      //compute stability
      DiscreteLinearForm<DirichletTrace<double>, LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> disc_stability(ansatz_space_);
      disc_stability.get_linear_form().set_function(unit_function);
      disc_stability.compute();
      Eigen::VectorXd s = disc_stability.get_discrete_linear_form();

      //compute gauge matrices
      Eigen::MatrixXd Sx = Dx.col(0) * s.segment(1,N-1).transpose();
      Eigen::MatrixXd Sy = Dy.col(0) * s.segment(1,N-1).transpose();
      Eigen::MatrixXd Sz = Dz.col(0) * s.segment(1,N-1).transpose();

      //arithmetic gauging
      Dx.block(0,1,M,N-1) -= Sx/s(0);
      Dy.block(0,1,M,N-1) -= Sy/s(0);
      Dz.block(0,1,M,N-1) -= Sz/s(0);


      //delete first column
      Dx.block(0,0,M,N-1) = Dx.block(0,1,M,N-1);
      Dy.block(0,0,M,N-1) = Dy.block(0,1,M,N-1);
      Dz.block(0,0,M,N-1) = Dz.block(0,1,M,N-1);

      Dx.conservativeResize(M,N-1);
      Dy.conservativeResize(M,N-1);
      Dz.conservativeResize(M,N-1);

      return s;

  }

  Eigen::VectorXd compute_stability_vector(){

      //function for stability vector
      std::function<double(Eigen::Vector3d)> unit_function = [](Eigen::Vector3d in){
        return 1.;
      };

      //compute stability
      DiscreteLinearForm<DirichletTrace<double>, LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> disc_stability(ansatz_space_);
      disc_stability.get_linear_form().set_function(unit_function);
      disc_stability.compute();
      Eigen::VectorXd s = disc_stability.get_discrete_linear_form();

      return s;

  }

  //this function imposes a vanishing mean gauge to the matrices Dx,Dy and Dz
  Eigen::VectorXd recover_full_solution(Eigen::VectorXd v,Eigen::VectorXd s){

      //number of unknowns
      int N = v.size()+1;

      Eigen::VectorXd v_new;

      v_new.resize(N);

      v_new.segment(1,N-1) = v;

      v_new(0) = -1*(s.segment(1,N-1).transpose() * v)(0)/s(0);

      return v_new;

  }


  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_cauchy_data(
      const Eigen::Matrix<typename LinearOperatorTraits<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>>::Scalar,
                          Eigen::Dynamic, 1> &cauchy_data) {
    fun_ev_.set_function(cauchy_data);
  }
  void set_degree(const int &deg) { deg_ = deg; }
  //////////////////////////////////////////////////////////////////////////////
  //    getter
  //////////////////////////////////////////////////////////////////////////////
  Derived &get_potential() { return pot_; }
  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  int deg_;
  Derived pot_;
  AnsatzSpace<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> ansatz_space_;
  FunctionEvaluator<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>> fun_ev_;


};  // namespace Bembel

}  // namespace Bembel
#endif
