// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_DISCRETE_SENSOR_H_
#define BEMBEL_DISCRETE_SENSOR_H_

#include <Eigen/Dense>

namespace Bembel {
  template <typename Derived, typename LinOp>
  class DiscreteSensor{
   public:
    //////////////////////////////////////////////////////////////////////////////
    //    constructors
    //////////////////////////////////////////////////////////////////////////////
    DiscreteSensor() {}
    DiscreteSensor(const AnsatzSpace<LinOp> &ansatz_space) {
      init_DiscreteSensor(ansatz_space);
    }
    //////////////////////////////////////////////////////////////////////////////
    //    init_Sensor
    //////////////////////////////////////////////////////////////////////////////
    void init_DiscreteSensor(const AnsatzSpace<LinOp> &ansatz_space) {

      ansatz_space_ = ansatz_space;

      fun_ev_ = FunctionEvaluator<LinOp>(ansatz_space_);

      /**
       * \todo obtain this from ansatz space
       */
      deg_ = ansatz_space_.get_polynomial_degree() + 1;
      return;
    }

    //////////////////////////////////////////////////////////////////////////////
    //    compute
    //////////////////////////////////////////////////////////////////////////////
    Eigen::Matrix<double, Eigen::Dynamic, SensorTraits<Derived>::OutputSpaceDimension>
    evaluate(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points) {
      auto FunctionSpaceVectorDimension =
          getFunctionSpaceVectorDimension<LinearOperatorTraits<LinOp>::Form>();
      auto OutputDimension = SensorTraits<Derived>::OutputSpaceDimension;

      GaussSquare<Constants::maximum_quadrature_degree> GS;
      auto Q = GS[deg_];

      auto super_space = ansatz_space_.get_superspace();

      auto element_tree = super_space.get_mesh().get_element_tree();
      auto number_of_elements = element_tree.get_number_of_elements();

      auto polynomial_degree = super_space.get_polynomial_degree();
      auto polynomial_degree_plus_one_squared =
          (polynomial_degree + 1) * (polynomial_degree + 1);

      //Synthetic Measurement Data
      Eigen::Matrix<double,Eigen::Dynamic,
                    SensorTraits<Derived>::OutputSpaceDimension>
          data;
      data.resize(points.rows(),SensorTraits<Derived>::OutputSpaceDimension);
      data.setZero();

  #pragma omp parallel
      {
        Eigen::Matrix<double,Eigen::Dynamic,
                      SensorTraits<Derived>::OutputSpaceDimension>
            my_data;
        my_data.resize(points.rows(),SensorTraits<Derived>::OutputSpaceDimension);
        my_data.setZero();
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
                my_data.row(i) += sensor_.evaluateIntegrand_impl(
                    fun_ev_, *element, points.row(i), qp);
              }
            }
          }
        }
  #pragma omp critical
        data += my_data;
      }
      return data;
    }


    //////////////////////////////////////////////////////////////////////////////
    //    Compute the Measurement Matrix
    //////////////////////////////////////////////////////////////////////////////
    void compute_measurement_matrix(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H,
                                          bool arithmetic_gauge = false) {

        //number of measurements
        int M = points.rows();

        //Output space dimension
        int dim_out = SensorTraits<Derived>::OutputSpaceDimension;

        //projector matrix
        auto projector = ansatz_space_.get_transformation_matrix();


        //number of basis functions
        int N_large = projector.rows();
        int N = ansatz_space_.get_number_of_dofs();

        //H matrix on local basis
        Eigen::MatrixXd H_l = Eigen::MatrixXd(M*dim_out,N_large);
        //H matrix on global basis
        H = Eigen::MatrixXd(M*dim_out,N);

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;



  //#pragma omp parallel
            {
  //#pragma omp single
              {


  #pragma omp parallel for
  for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_H;

        tmp_H.resize(dim_out,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){



              tmp_H.setZero();


               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrand(super_space, *element, points.row(m),qp, &tmp_H);
                }

                //#pragma omp critical
                {
                  for(int i = 0 ; i < dim_out ; ++i){
                    H_l.block(m+i*M,polynomial_degree_plus_one_squared * element->id_,
                                   1,polynomial_degree_plus_one_squared) = tmp_H.row(i);
                  }
                }
            }
            }
          }
        }



      //go back to continuous basis
      for(int i = 0 ; i < dim_out ; ++i){
        H.block(i*M,0,M,N) = H_l.block(i*M,0,M,N_large) * projector;
      }

      if(arithmetic_gauge){
        stability_vec_ = impose_zero_mean_gauge(H);
      }

      return;
    }

    //////////////////////////////////////////////////////////////////////////////
    //    Compute the Measurement Matrix
    //////////////////////////////////////////////////////////////////////////////
    void compute_measurement_matrix_low_ram(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H,
                                          bool arithmetic_gauge = false) {

        //number of measurements
        int M = points.rows();

        //Output space dimension
        int dim_out = SensorTraits<Derived>::OutputSpaceDimension;

        //projector matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
        projector = ansatz_space_.get_transformation_matrix();


        //number of basis functions

        int N = ansatz_space_.get_number_of_dofs();


        //H matrix on global basis
        H = Eigen::MatrixXd(M*dim_out,N);
        H.setZero();

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;



    //#pragma omp parallel
            {
    //#pragma omp single
              {


    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_H;

        tmp_H.resize(dim_out,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){



              tmp_H.setZero();


               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrand(super_space, *element, points.row(m),qp, &tmp_H);
                }

                //#pragma omp critical
                {

                  for(int i = 0 ; i < dim_out ; ++i){

                      //the projector provides the connectivity to assemble the
                      //small matrix directly.

                      //iterate over the Bernstein basis functions
                      for(int k = 0; k < polynomial_degree_plus_one_squared; ++k){

                        //iterate over the projectors column, corresponding to the
                        //Bernstein basis function k
                        for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(projector,polynomial_degree_plus_one_squared * element->id_ + k); it; ++it)
                        {

                          H(m+i*M,it.col()) += it.value()*tmp_H(i,k);

                        }
                      }
                  }
                }
            }
            }
          }
        }

      if(arithmetic_gauge){
        stability_vec_ = impose_zero_mean_gauge(H);
      }

      return;
    }

    /*--------------------------------------------------------------------------
    Computation of H matrices for UQ
    --------------------------------------------------------------------------*/
    void compute_H_matrices(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dtH_1,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dtH_2,
                                        const int offset) {

        //number of measurements
        int M = points.rows();


        //projector matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
        projector = ansatz_space_.get_transformation_matrix();

        //number of basis functions
        int N = ansatz_space_.get_number_of_dofs();

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;


    //#pragma omp parallel
            {
    //#pragma omp single
              {


    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_H;

        tmp_H.resize(3,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){



              tmp_H.setZero();


               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrands_UQ(super_space, *element, points.row(m),qp, &tmp_H);
                }

                //#pragma omp critical
                {

                  //the projector provides the connectivity to assemble the
                  //small matrix directly.

                  //iterate over the Bernstein basis functions
                  for(int k = 0; k < polynomial_degree_plus_one_squared; ++k){

                    //iterate over the projectors column, corresponding to the
                    //Bernstein basis function k
                    for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(projector,polynomial_degree_plus_one_squared * element->id_ + k); it; ++it)
                    {

                      H(     m+offset , it.col() )  += it.value()*tmp_H(0,k);
                      dtH_1( m+offset , it.col() )  += it.value()*tmp_H(1,k);
                      dtH_2( m+offset , it.col() )  += it.value()*tmp_H(2,k);

                    }
                  }

                }
            }
            }
          }
        }

      return;
    }


    void compute_H(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H,
                              const int offset) {

        //number of measurements
        int M = points.rows();

        //projector matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
        projector = ansatz_space_.get_transformation_matrix();

        //number of basis functions
        int N = ansatz_space_.get_number_of_dofs();

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;


    //#pragma omp parallel
            {
    //#pragma omp single
              {


    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_H;

        tmp_H.resize(1,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){


              tmp_H.setZero();
               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrand(super_space, *element, points.row(m),qp, &tmp_H);
                }


                {

                  //the projector provides the connectivity to assemble the
                  //small matrix directly.
                  //iterate over the Bernstein basis functions
                  for(int k = 0; k < polynomial_degree_plus_one_squared; ++k){

                    //iterate over the projectors column, corresponding to the
                    //Bernstein basis function k
                    for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(projector,polynomial_degree_plus_one_squared * element->id_ + k); it; ++it)
                    {
                      H(m+offset,it.col()) += it.value()*tmp_H(0,k);
                    }
                  }

                }
            }
            }
          }
        }

      return;
    }

    //////////////////////////////////////////////////////////////////////////////
    //    Compute the Measurement Derivative Matrix
    //////////////////////////////////////////////////////////////////////////////
    void compute_measurement_derivative_matrix(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH,
                                          bool arithmetic_gauge = false) {

        //number of measurements
        int M = points.rows();

        //Output space dimension
        int dim_out = 3*SensorTraits<Derived>::OutputSpaceDimension;

        //projector matrix
        auto projector = ansatz_space_.get_transformation_matrix();


        //number of basis functions
        int N_large = projector.rows();
        int N = ansatz_space_.get_number_of_dofs();

        //dH matrix on local basis
        Eigen::MatrixXd dH_l = Eigen::MatrixXd(M*dim_out,N_large);

        //H matrix on global basis
        dH.resize(M*dim_out,N);

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;


    //#pragma omp parallel
            {
    //#pragma omp single
              {


    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_dH;

        tmp_dH.resize(dim_out,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){



              tmp_dH.setZero();

               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrandDerivative(super_space, *element, points.row(m),qp, &tmp_dH);
                }

                //#pragma omp critical
                {
                  for(int i = 0 ; i < dim_out ; ++i){
                    dH_l.block(m+i*M,polynomial_degree_plus_one_squared * element->id_,
                                   1,polynomial_degree_plus_one_squared) = tmp_dH.row(i);
                  }
                }
            }
            }
          }
        }



      //go back to continuous basis
      for(int i = 0 ; i < dim_out ; ++i){
        dH.block(i*M,0,M,N) = dH_l.block(i*M,0,M,N_large) * projector;
      }

      if(arithmetic_gauge){
        stability_vec_ = impose_zero_mean_gauge(dH);
      }

      return;
    }

    //////////////////////////////////////////////////////////////////////////////
    //    Compute all matrices needed for UQ at once. We avoid evaluating the
    //    surface multiple times
    //////////////////////////////////////////////////////////////////////////////
    void compute_UQ_matrices(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dtH_1,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dtH_2,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_x,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_y,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_z,
                                        const int offset = 0) {

        //number of measurements
        int M = points.rows();

        //projector matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
        projector = ansatz_space_.get_transformation_matrix();

        //number of basis functions
        int N = ansatz_space_.get_number_of_dofs();

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;



    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_H;
        Eigen::MatrixXd tmp_dH;

        tmp_H.resize(3,polynomial_degree_plus_one_squared);
        tmp_dH.resize(3,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){

              tmp_H.setZero();
              tmp_dH.setZero();

               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrands_UQ(super_space, *element, points.row(m),qp, &tmp_H);
                  sensor_.evaluateIntegrandDerivative(super_space, *element, points.row(m),qp, &tmp_dH);


                }

                //iterate over the Bernstein basis functions
                for(int k = 0; k < polynomial_degree_plus_one_squared; ++k){

                  //iterate over the projectors column, corresponding to the
                  //Bernstein basis function k
                  for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(projector,polynomial_degree_plus_one_squared * element->id_ + k); it; ++it)
                  {

                    H(     m+offset , it.col() )  += it.value()*tmp_H(0,k);
                    dtH_1( m+offset , it.col() )  += it.value()*tmp_H(1,k);
                    dtH_2( m+offset , it.col() )  += it.value()*tmp_H(2,k);

                    dH_x(m + offset , it.col()) += it.value()*tmp_dH(0,k);
                    dH_y(m + offset , it.col()) += it.value()*tmp_dH(1,k);
                    dH_z(m + offset , it.col()) += it.value()*tmp_dH(2,k);

                  }
              }
          }
        }

      return;
    }

  //////////////////////////////////////////////////////////////////////////////
  //    Compute the derivatives of the measurement operation directly
  //////////////////////////////////////////////////////////////////////////////
  Eigen::MatrixXd compute_voltage_and_derivatives(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        const Eigen::VectorXd &v) {

    //number of measurements
    int M = points.rows();

    //projector matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
    projector = ansatz_space_.get_transformation_matrix();

    //go discontinuous in v
    Eigen::VectorXd v_disc = (projector * v).eval();

    //number of basis functions
    int N = ansatz_space_.get_number_of_dofs();

    GaussSquare<Constants::maximum_quadrature_degree> GS;

    auto Q = GS[deg_];

    auto super_space = ansatz_space_.get_superspace();
    auto element_tree = super_space.get_mesh().get_element_tree();

    auto number_of_elements = element_tree.get_number_of_elements();

    auto polynomial_degree = super_space.get_polynomial_degree();
    auto I2 = (polynomial_degree + 1) * (polynomial_degree + 1);

    //return matrix
    Eigen::MatrixXd ret_mat(M,6);
    ret_mat.setZero();

    std::shared_ptr<Bembel::ElementTreeMemory> el_memory_ = ansatz_space_.get_superspace().get_mesh().get_element_tree().root().memory_;

    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

      //temporal storage for matrix blocks
      Eigen::MatrixXd tmp_H;
      Eigen::MatrixXd tmp_dH;

      tmp_H.resize(3,I2);
      tmp_dH.resize(3,I2);

      //point on surface for integratio
      SurfacePoint qp;

      //iterate over all elements
      for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){

          tmp_H.setZero();
          tmp_dH.setZero();

          //Gaussian integtation
          for (auto j = 0; j < Q.w_.size(); ++j) {

              super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

              sensor_.evaluateIntegrands_UQ(super_space, *element, points.row(m),qp, &tmp_H);
              sensor_.evaluateIntegrandDerivative(super_space, *element, points.row(m),qp, &tmp_dH);



          }


          //position of this element in memory
          int k0 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_begin(*element));

          ret_mat.block(m,0,1,3) += (tmp_dH*v_disc.segment(k0,I2)).transpose();
          ret_mat.block(m,3,1,3) += (tmp_H*v_disc.segment(k0,I2)).transpose();


        }
      }

    return ret_mat;
  }


    //////////////////////////////////////////////////////////////////////////////
    //    Compute the Measurement Derivative Matrix
    //////////////////////////////////////////////////////////////////////////////
    void compute_measurement_derivative_matrix_low_ram(const Eigen::Matrix<double, Eigen::Dynamic, 3> &points,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_x,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_y,
                                        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &dH_z,
                                        const int row_offset = 0) {

        //number of measurements
        int M = points.rows();

        //projector matrix
        Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
        projector = ansatz_space_.get_transformation_matrix();

        //number of basis functions
        int N = ansatz_space_.get_number_of_dofs();

        GaussSquare<Constants::maximum_quadrature_degree> GS;

        auto Q = GS[deg_];

        auto super_space = ansatz_space_.get_superspace();
        auto element_tree = super_space.get_mesh().get_element_tree();

        auto number_of_elements = element_tree.get_number_of_elements();

        auto polynomial_degree = super_space.get_polynomial_degree();
        auto polynomial_degree_plus_one_squared =
            (polynomial_degree + 1) * (polynomial_degree + 1);

        int k = 0;



    #pragma omp parallel for
    for (int m = 0; m < M; ++m){

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_dH;

        tmp_dH.resize(3,polynomial_degree_plus_one_squared);

        //point on surface for integratio
        SurfacePoint qp;

        //iterate over all elements
        for (auto element = element_tree.cpbegin();
                 element != element_tree.cpend(); ++element){

              tmp_dH.setZero();

               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.evaluateIntegrandDerivative(super_space, *element, points.row(m),qp, &tmp_dH);
                }

                //iterate over the Bernstein basis functions
                for(int k = 0; k < polynomial_degree_plus_one_squared; ++k){

                  //iterate over the projectors column, corresponding to the
                  //Bernstein basis function k
                  for (Eigen::SparseMatrix<double,Eigen::RowMajor>::InnerIterator it(projector,polynomial_degree_plus_one_squared * element->id_ + k); it; ++it)
                  {

                    dH_x(m+row_offset,it.col()) = it.value()*tmp_dH(0,k);
                    dH_y(m+row_offset,it.col()) = it.value()*tmp_dH(1,k);
                    dH_z(m+row_offset,it.col()) = it.value()*tmp_dH(2,k);

                  }
              }
          }
        }

      return;
    }

    //this function imposes a vanishing mean gauge to the matrix H
    Eigen::VectorXd impose_zero_mean_gauge(Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &H){

        //number of unknowns
        int N = H.cols();

        //number of evaluations
        int M = H.rows();

        //function for stability vector
        std::function<double(Eigen::Vector3d)> unit_function = [](Eigen::Vector3d in){
          return 1.;
        };

        //compute stability
        DiscreteLinearForm<DirichletTrace<double>, LinOp> disc_stability(ansatz_space_);
        disc_stability.get_linear_form().set_function(unit_function);
        disc_stability.compute();
        Eigen::VectorXd s = disc_stability.get_discrete_linear_form();

        //compute gauge matrices
        Eigen::MatrixXd S = H.col(0) * s.segment(1,N-1).transpose();

        //arithmetic gauging
        H.block(0,1,M,N-1) -= S/s(0);

        //delete first column
        H.block(0,0,M,N-1) = H.block(0,1,M,N-1);

        H.conservativeResize(M,N-1);

        return s;

    }

    Eigen::VectorXd recover_full_solution(Eigen::VectorXd v,Eigen::VectorXd s){

        //number of unknowns
        int N = v.size()+1;

        Eigen::VectorXd v_new;

        v_new.resize(N);

        v_new.segment(1,N-1) = v;

        v_new(0) = -1*(s.segment(1,N-1).transpose() * v)(0)/s(0);

        return v_new;

    }

    double get_sensitivity(){

      return sensor_.get_sensitivity();

    }

    Eigen::VectorXcd compute_local_evaluation_vector(const int L,const Eigen::Vector3d &eval_pos, const Eigen::Vector3d &cell_center){

      // number of coefficients
      int num_coeffs = (L+1)*(L+1);

      Eigen::VectorXcd ret_vec(num_coeffs);
      ret_vec.setZero();

      // evaluate the gradient of the solid harmonics
      Eigen::MatrixXcd grad_R_lm = Rlm_p_alt(L,eval_pos,cell_center);//grad_R_lm(num_multipoles_, eval_pos - cell_center);

      //tranfer function of this sensor
      std::vector<double> tf = sensor_.get_transfer_function();

      //orientation vector
      Eigen::Vector3d n = sensor_.get_orientation_vector();

      ret_vec = n(0)*grad_R_lm.block(0,0,num_coeffs,1)
                + n(1)*grad_R_lm.block(0,1,num_coeffs,1)
                + n(2)*grad_R_lm.block(0,2,num_coeffs,1);

      return -1.*tf[1]*ret_vec;

    }

    //////////////////////////////////////////////////////////////////////////////
    //    setter
    //////////////////////////////////////////////////////////////////////////////
    void set_cauchy_data(
        const Eigen::Matrix<typename LinearOperatorTraits<LinOp>::Scalar,
                            Eigen::Dynamic, 1> &cauchy_data) {
      fun_ev_.set_function(cauchy_data);
    }
    void set_degree(const int &deg) { deg_ = deg; }
    //////////////////////////////////////////////////////////////////////////////
    //    getter
    //////////////////////////////////////////////////////////////////////////////
    Derived &get_sensor() { return sensor_; }
    //////////////////////////////////////////////////////////////////////////////
    //    private member variables
    //////////////////////////////////////////////////////////////////////////////
   private:
    int deg_;
    Derived sensor_;
    AnsatzSpace<LinOp> ansatz_space_;
    FunctionEvaluator<LinOp> fun_ev_;
    Eigen::VectorXd stability_vec_;


  };  // namespace Bembel

  }  // namespace Bembel
  #endif
