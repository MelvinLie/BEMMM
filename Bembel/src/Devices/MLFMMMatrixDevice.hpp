// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_DEVICES_MLFMMMATRIXDEVICE_H__
#define __BEMBEL_DEVICES_MLFMMMATRIXDEVICE_H__


/**
 *  \class LRMatrix
 *  \brief Low Rank Matrix class, which extends the EigenBase class.
 *
 *  The idea is to provide an easy to use interface to the H2-matrix
 *  from the fast boundary element method. At the moment, we inherit the
 *  traits of an Eigen::SparseMatrix, since this seems to be the minimum
 *  properties for a derived object to ensure that the matrix-vector
 *  multiplication can be specialised for H2Matrix.
 *  In particular, this allows for the use of the Eigen iterative solvers
 *  with a Hierarchical matrix.
 *
 *  \todo Maybe, we find something better then the SparsMatrix traits
 *        in the future
 **/

namespace Eigen {
/// forward definition of the LRMatrix Class in order to define traits
template <typename Derived, typename LinOp>
class MLFMMMatrixDevice;
/// inherit the traits from the Eigen::SparseMatrix class
namespace internal {
template <typename Derived, typename LinOp>
struct traits<MLFMMMatrixDevice<Derived,LinOp>> : public internal::traits<SparseMatrix<double>> {};
}  // namespace internal



// actual definition of the class
template <typename Derived, typename LinOp>
class MLFMMMatrixDevice : public EigenBase<MLFMMMatrixDevice<Derived,LinOp>> {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// Eigen related things
  //////////////////////////////////////////////////////////////////////////////
  // Required typedefs, constants and so on.
  typedef double Scalar;
  typedef typename NumTraits<double>::Real RealScalar;
  typedef Index StorageIndex;
  enum {
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    IsRowMajor = false,
    Flags = NestByRefBit
  };
  // Minimum specialisation of EigenBase methods
  Index rows() const { return rows_; }
  Index cols() const { return cols_; }
  // Definition of the matrix multiplication
  template <typename Rhs>
  Product<MLFMMMatrixDevice<Derived,LinOp>, Rhs, AliasFreeProduct> operator*(
      const MatrixBase<Rhs>& x) const {
    return Product<MLFMMMatrixDevice<Derived,LinOp>, Rhs, AliasFreeProduct>(*this, x.derived());
  }
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  MLFMMMatrixDevice() {}
  /**
   * \brief Assemble H2-Matrix for linear operator linOp and AnsatzSpace
   * ansatz_space with number_of_points interpolation points in one direction of
   * the unit square (standard is number_of_points=16)
   */
  //Make this a template class to be able to handle different devices. So far we
  //are using the Hall probe array
  //template <typename Derived>
  void init_MLFMMMatrix(Bembel::HallProbeArray<Derived,LinOp> *device, Bembel::MeasurementData *meas_data, int max_tree_lvl, const Eigen::Matrix<double,2,3> &bbox) {

    // device
    device_ = device;
    // ansatz space
    ansatz_space_ = device->get_ansatz_space_ptr();
    // measurement data
    meas_data_ = meas_data;
    //discrete Potential
    potential_ = device_->get_sensor(0).get_sensor().get_potential_ptr();

    // some reoccuring variables
    I_ = ansatz_space_->get_polynomial_degree();
    I2_ = (I_+1)*(I_+1);
    num_coeffs_ = (num_multipoles_+1)*(num_multipoles_+1);
    int num_meas = meas_data->get_number_of_measurements();
    num_sensors_ = device_->get_number_of_sensors();

    // dimensions of the matrix
    //cols_ = ansatz_space_->get_transformation_matrix().rows();
    cols_ = ansatz_space_->get_number_of_dofs();
    cols_disc_ = ansatz_space_->get_transformation_matrix().rows();
    rows_ = num_meas * num_sensors_;


    //sensor positions
    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> positions = device_->get_sensor_positions();

    // get superspace for integration
    Bembel::SuperSpace<LinOp> super_space = ansatz_space_->get_superspace();


    // get element tree from SuperSpace
    //auto element_tree = super_space.get_mesh().get_element_tree();

    // maximum recursion level for tree generation
    max_tree_lvl_ = max_tree_lvl;

    // setup element tree
    el_tree_.set_min_numel(min_numel_);
    el_tree_.init_ElementOctTree(ansatz_space_,max_tree_lvl,bbox);
    //el_tree_.print_tree();


    // setup the measurement tree
    meas_data->set_min_numel(min_numel_);
    meas_data->set_bounding_box(bbox);
    meas_data->init_cluster_tree(max_tree_lvl);
    //meas_data->print_cluster_tree();

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //measurement positions
    Eigen::MatrixXd meas_pos = meas_data->get_positions();

    //pointer to element cluster memory. With this pointer we access the elements for integration
    std::vector<Bembel::ElementTreeNode>::const_iterator mesh_it = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();



    //--------------------------------------------------------------------------
    //Compute moment matrices
    //--------------------------------------------------------------------------

    //quadrature degree
    deg_ = ansatz_space_->get_polynomial_degree() + 1;

    //number of basis functions on one element
    int polynomial_degree_plus_one_squared = deg_*deg_;



    //number of multipole coefficients
    int num_coeffs = (num_multipoles_+1)*(num_multipoles_+1);


    //Integration rule
    Bembel::GaussSquare<Bembel::Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];


    //compute moment matrices
    m2m_matrices_.fill(num_multipoles_, max_tree_lvl, bbox, Bembel::MultipoleMomentContainer::SOURCE);

    //compute local matrices
    l2l_matrices_.fill(num_multipoles_, max_tree_lvl, bbox, Bembel::MultipoleMomentContainer::LOCAL);

    //diameter of box
    compute_box_diam(bbox);

    //--------------------------------------------------------------------------
    //setup interaction region and near field
    //--------------------------------------------------------------------------
    //std::cout << "Make interaction region and near field" << std::endl;
    make_ir_and_nf();


    //number of DoFs
    int num_DoFs = ansatz_space_->get_number_of_dofs();

    //projector matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
    projector = ansatz_space_->get_transformation_matrix();

    Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> projector_cmplx;
    projector_cmplx =  projector.cast<std::complex<double>>();

    //number of DoFs in Bernstein basis
    int num_DoFs_disc = projector.rows();

    //--------------------------------------------------------------------------
    //Make leaf moment matrices
    //--------------------------------------------------------------------------
    //std::cout << "Make leaf moment matrices" << std::endl;
    for(int l = max_tree_lvl_+1; l >= 0 ; --l){

      //std::vector<Bembel::ElementOctTreeNode>::iterator el_it;

#pragma omp parallel for
      for(std::vector<Bembel::ElementOctTreeNode>::iterator el_it = el_tree_mem->lbegin(l); el_it < el_tree_mem->lend(l) ; ++el_it){

        //mesh tree
        std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = ansatz_space_->get_superspace().get_mesh().get_element_tree().root().memory_;

        //surface point
        SurfacePoint qp;

        //moments are computed for the leafs only (#sons = 0)
        if((*el_it).sons_.size() == 0){

          //integration result
          Eigen::MatrixXcd int_mat;
          int_mat.resize(num_coeffs,I2_);

          //make space for multipole moment matrix
          (*el_it).moment_mat_.resize(num_coeffs,num_DoFs);

          //running index
          //int k = 0;

          //make space for temporal sparse matrix on discontinuous basis
          Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> tmp_l(num_coeffs,num_DoFs_disc);

          //make tripletList to fill sparse moment matrix
          typedef Eigen::Triplet<std::complex<double>> T;
          std::vector<T> tripletList;
          tripletList.reserve(num_coeffs*(*el_it).indices_.size()*I2_);

          //iterate over the elements of this cluster
          for(int el_indx: (*el_it).indices_){


            //zero integration results
            int_mat.setZero();

            for (auto j = 0; j < Q.w_.size(); ++j) {
              super_space.map2surface(mesh_it[el_indx], Q.xi_.col(j),mesh_it[el_indx].get_h() * mesh_it[el_indx].get_h() * Q.w_(j), &qp);

              potential_->evaluateFMMMoments(super_space,mesh_it[el_indx], (*el_it).center_, qp, &int_mat,num_multipoles_);

            }

            //position of this element in memory
            int k0 = I2_*std::distance(mesh_memory->cpbegin(),mesh_memory->cluster_begin(mesh_it[el_indx]));

            //append to tripletList
            for( int k = 0; k < I2_; ++k){
              for( int m = 0; m < num_coeffs; ++m){
                  tripletList.push_back(T(m , k0 + k , int_mat(m,k)/4./BEMBEL_PI ));
              }
            }

          }

          //make sparse matrix
          tmp_l.setFromTriplets(tripletList.begin(), tripletList.end());

          //setup moment matrix
          (*el_it).moment_mat_ = (tmp_l * projector_cmplx).eval();

          //remeber the position of this leaf in the memory
          #pragma omp critical
          {
            source_leafs_.push_back(el_it->pos_);
          }

        }

      }

    }

    //--------------------------------------------------------------------------
    //Make leafs local matrices
    //--------------------------------------------------------------------------
    //std::cout << "Make leafs local matrices" << std::endl;
    for(int l = max_tree_lvl_+1; l >= 2 ; --l){

#pragma omp parallel for
      for(std::vector<Bembel::MeasurementTreeNode>::iterator meas_it = meas_tree_mem->lbegin(l); meas_it < meas_tree_mem->lend(l) ; meas_it++){


        //difference vector
        Eigen::Vector3d d;

        //locals are computed for the leafs only (#sons = 0)
        if(meas_it->sons_.size() == 0){

          //make space for local matrix
          meas_it->local_mat_.resize(rows_,num_coeffs);


          //make tripletList to fill sparse moment matrix
          typedef Eigen::Triplet<std::complex<double>> T;
          std::vector<T> tripletList;
          tripletList.reserve(meas_it->indices_.size()*num_coeffs);

          //iterate over the measurements of this cell
          for(int meas_indx: meas_it->indices_){

            //d = meas_pos.row(meas_indx).transpose() - meas_it->center_;

            //Eigen::VectorXcd Rlm = Bembel::Rlm_alt
            Eigen::MatrixXcd eval_mat = device_->compute_local_evaluation_matrix(num_multipoles_,meas_pos.row(meas_indx),meas_it->center_);

            for(int k = 0; k < eval_mat.rows(); ++ k){
              for(int s = 0; s < num_sensors_ ; ++s){
                tripletList.push_back(T(meas_indx + s*num_meas , k , eval_mat(k,s) ));
              }
            }

          }
          meas_it->local_mat_.setFromTriplets(tripletList.begin(), tripletList.end());

        //remeber the position of this leaf in the memory
        #pragma omp critical
        {
          local_leafs_.push_back(meas_it->pos_);
        }

        }

      }
    }

    //--------------------------------------------------------------------------
    //Make sparse near field matrix
    //--------------------------------------------------------------------------
    //std::cout << "Make sparse near field matrix" << std::endl;

    //projector matrix
    //Eigen::SparseMatrix<double,Eigen::RowMajor> projector; //we need it row major for what we will do lateron
    //projector = ansatz_space_->get_transformation_matrix();

    //number of BSpline basis functions
    int N = ansatz_space_->get_number_of_dofs();

    H_near_.resize(rows_,N);
    H_near_.setZero();

    // we are integrating over the hyptersingular kernel, we need to increase
    // the quadrature degree drastically!
    Q = GS[deg_ + 2];

    #pragma omp parallel
    {

      Eigen::SparseMatrix<double,Eigen::RowMajor> my_H(rows_,N);

      //local and source indices
      int i_l,i_s;


      int k0,k1;

      //sensor transfer function
      std::vector<double> tf;

      //surface point
      SurfacePoint qp;

      //mesh tree
      std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = ansatz_space_->get_superspace().get_mesh().get_element_tree().root().memory_;

      #pragma omp for
      for(int i = 0; i < near_field_interaction_list_.size() ; ++i){

        //indices of interaction
        i_l = near_field_interaction_list_[i][1];
        i_s = near_field_interaction_list_[i][0];

        //interacting cells
        Bembel::MeasurementTreeNode *meas_node = meas_tree_mem->get_element_ptr(i_l);
        Bembel::ElementOctTreeNode *source_node = el_tree_mem->get_element_ptr(i_s);


        //tripletList to construct sparse matrix
        typedef Eigen::Triplet<double> T;
        std::vector<T> tripletList;
        tripletList.reserve(num_sensors_*num_near_field_interactions_);

        Eigen::SparseMatrix<double,Eigen::RowMajor> this_H(rows_,cols_disc_);

        for(int source_indx : source_node->indices_){

          Eigen::MatrixXd int_mat, int_vec ; //int_vec_p;
          int_mat.resize(I2_,num_sensors_);
          int_vec.resize(I2_,1);



          //std::cout << "source index = " << source_indx << std::endl;

          for(int meas_indx : meas_node->indices_){


            int_mat.setZero();

            //compute dense matrix block
            for (auto j = 0; j < Q.w_.size(); ++j) {


              super_space.map2surface(mesh_it[source_indx], Q.xi_.col(j),mesh_it[source_indx].get_h()  * mesh_it[source_indx].get_h()  * Q.w_(j), &qp);//

              for(int s = 0 ; s < num_sensors_ ; ++s){

                int_vec.setZero();

                //potential_->evaluateIntegrand_impl(super_space,mesh_it[source_indx], meas_pos.row(meas_indx).transpose(), qp, &int_vec_0);
                //potential_->evaluateIntegrand_impl(super_space,mesh_it[source_indx], meas_pos.row(meas_indx).transpose() + x_p, qp, &int_vec_p);

                //potential_->evaluateIntegrandDerivative_mat(super_space,mesh_it[source_indx],meas_pos.row(meas_indx).transpose(),qp,&int_mat_der);

                //potential_->evaluateIntegrandDerivative_far_field(super_space,mesh_it[source_indx],meas_pos.row(meas_indx).transpose(),qp,&int_mat_der);
                /**/
                device_->get_sensor(s).get_sensor().evaluateIntegrand_far_field(super_space,mesh_it[source_indx],
                                                                meas_pos.row(meas_indx).transpose() + positions[s],
                                                                qp,
                                                                &int_vec);

                //std::cout << "mesh parameter = " << mesh_it[source_indx].get_h() << std::endl;
                //std::cout << "analytical = " << int_vec.transpose() << std::endl;
                //std::cout << "approx = " << -1.*(int_vec_p - int_vec_0).transpose()/1e-4 << std::endl;
                //std::cout << "der = " << int_mat_der.col(0).transpose() << std::endl;
                //int_mat.block(0,s,I2_,1) -= (int_vec_p - int_vec_0)/1e-4;//int_vec.block(0,0,1,I2_).transpose();
                //int_mat.block(0,s,I2_,1) += int_vec.block(0,0,1,I2_).transpose();

                int_mat.block(0,s,I2_,1) += int_vec;

              }
            }


            k0 = I2_*std::distance(mesh_memory->cpbegin(),mesh_memory->cluster_begin(mesh_it[source_indx]));
            //k1 = I2_*std::distance(mesh_memory->cpbegin(),mesh_memory->cluster_end(mesh_it[source_indx]));

            for(int s = 0 ; s < num_sensors_ ; ++s){
              for( int k = 0; k < I2_; ++k){


                double sens = device_->get_sensor(s).get_sensor().get_transfer_function()[1];

                tripletList.push_back(T(meas_indx + s*num_meas, k0 + k , sens*int_mat(k,s) ));

              }
            }
          }

        }

        this_H.setFromTriplets(tripletList.begin(), tripletList.end());

        //go continuous
        this_H = (this_H * projector).eval();

        //add to my_H
        my_H += this_H;

      }

      #pragma omp critical
      {
        H_near_ += my_H;
      }

    }
    /*
    */

  }

  /*
  double count_RAM(){

    //Moment matrices
    double ram_moments_matrices = count_moment_matrix_RAM();
    //std::cout << "RAM of moment matrices = " << ram_moments_matrices << " GB"<< std::endl;

    double m2m_ram = m2m_matrices_.count_RAM();
    //std::cout << "RAM of m2m matrices = " << m2m_ram << " GB"<< std::endl;

    double l2l_ram = l2l_matrices_.count_RAM();
    //std::cout << "RAM of l2l matrices = " << l2l_ram << " GB"<< std::endl;

    double m2l_ram = count_m2l_ram();
    //std::cout << "RAM of m2l matrices = " << m2l_ram << " GB"<< std::endl;

    double locals_ram = count_locals_RAM();
    //std::cout << "RAM of local matrices = " << locals_ram << " GB"<< std::endl;

    double near_field_ram = H_near_.nonZeros()*8e-9;
    //std::cout << "RAM of near field matrix = " << near_field_ram << " GB"<< std::endl;

    return ram_moments_matrices + m2m_ram + l2l_ram + m2l_ram + locals_ram + near_field_ram;

  }
  */
  void set_number_of_multipoles(const int num_multipoles){

    num_multipoles_ = num_multipoles;
    num_coeffs_ = (num_multipoles_+1)*(num_multipoles_+1);

  }
  /*
  void setup_normal_equations(const Eigen::SparseMatrix<double> R){

    //
  }
  */
  Eigen::VectorXd matvec(const Eigen::VectorXd &rhs){

    //for time measurement
    std::chrono::time_point<std::chrono::steady_clock> t_start, t_end;
    std::chrono::duration<double> t_el;



    //projector matrix
    //auto projector = ansatz_space_->get_transformation_matrix();

    //transform to discontinuous space
    //Eigen::VectorXd rhs_l = projector * rhs;

    //std::cout << "compute leaf moments" << std::endl;

    //t_start =  std::chrono::steady_clock::now();

    compute_leaf_moments(rhs);

    //t_end =  std::chrono::steady_clock::now();
    //t_el = t_end - t_start;

    //std::cout << "Elapsed time = " << t_el.count() << " sec" << std::endl;

    //std::cout << "upward pass" << std::endl;

    //t_start =  std::chrono::steady_clock::now();

    upward_pass();

    //t_end =  std::chrono::steady_clock::now();
    //t_el = t_end - t_start;

    //std::cout << "Elapsed time = " << t_el.count() << " sec" << std::endl;

    //std::cout << "downward pass" << std::endl;

    //t_start =  std::chrono::steady_clock::now();

    downward_pass();

    //t_end =  std::chrono::steady_clock::now();
    //t_el = t_end - t_start;

    //std::cout << "Elapsed time = " << t_el.count() << " sec" << std::endl;

    //std::cout << "evaluate near field" << std::endl;

    //t_start =  std::chrono::steady_clock::now();

    Eigen::VectorXd near_field = evaluate_near_field(rhs);

    //t_end =  std::chrono::steady_clock::now();
    //t_el = t_end - t_start;

    //std::cout << "near field max coeff = " << near_field.maxCoeff() << std::endl;
    //std::cout << "Elapsed time = " << t_el.count() << " sec" << std::endl;

    //std::cout << "evaluate far field" << std::endl;

    //t_start =  std::chrono::steady_clock::now();


    Eigen::VectorXd far_field = evaluate_far_field();

    //t_end =  std::chrono::steady_clock::now();
    //t_el = t_end - t_start;

    //std::cout << "Elapsed time = " << t_el.count() << " sec" << std::endl;

    //return evaluate_multipole_expansion(meas_data_->get_positions(),rhs);
    return  near_field + far_field;


  }

  Eigen::VectorXd matvec_transpose(const Eigen::VectorXd &rhs){

    std::cout << "compute_leaf_moments_transpose" << std::endl;
    compute_leaf_moments_transpose(rhs);

    std::cout << "upward_pass_transpose" << std::endl;
    upward_pass_transpose();

    std::cout << "downward_pass_transpose" << std::endl;
    downward_pass_transpose();

    std::cout << "evaluate_near_field_transpose" << std::endl;
    Eigen::VectorXd near_field = evaluate_near_field_transpose(rhs);

    std::cout << "evaluate_far_field_transpose" << std::endl;
    Eigen::VectorXd far_field = evaluate_far_field_transpose();

    //std::cout << "near_field = " << near_field.rows() << " x " << near_field.cols() << std::endl;
    //std::cout << "far_field = " << far_field.rows() << " x " << far_field.cols() << std::endl;
    return near_field + far_field;

  }

  void matvec(const Eigen::VectorXd &rhs, Eigen::VectorXd *near_field_out, Eigen::VectorXd *far_field_out){

    //projector matrix
    auto projector = ansatz_space_->get_transformation_matrix();

    //transform to discontinuous space
    Eigen::VectorXd rhs_l = projector * rhs;

    //std::cout << "compute leaf moments" << std::endl;
    compute_leaf_moments(rhs);

    //std::cout << "upward pass" << std::endl;
    upward_pass();

    //std::cout << "downward pass" << std::endl;
    downward_pass();

    //std::cout << "evaluate near field" << std::endl;
    (*near_field_out) = evaluate_near_field(rhs);

    //std::cout << "evaluate far field" << std::endl;
    (*far_field_out) = evaluate_far_field();

    return;

  }

  /*
  Eigen::VectorXd evaluate_multipole_expansion(const Eigen::MatrixXd &pos, const Eigen::VectorXd &density){

        //number of evaluations
        int num_eval = pos.rows();

        //make space for potential
        Eigen::VectorXd ret_vec;
        ret_vec.resize(num_eval);
        ret_vec.setZero();

        //get access to the ansatz space
        //Bembel::AnsatzSpace<LinOp> *ansatz_space = potential_->get_ansatz_space();

        //polynomial_degree_plus_one_squared
        int I = ansatz_space_->get_polynomial_degree();
        int I2_ = (I+1)*(I+1);

        //projector matrix
        auto projector = ansatz_space_->get_transformation_matrix();

        //transform to discontinuous space
        Eigen::VectorXd density_l = projector * density;

        //element tree memory
        std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

        //mesh tree
        std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = ansatz_space_->get_superspace().get_mesh().get_element_tree().root().memory_;

        //pointer to element cluster memory. With this pointer we access the elements for integration
        std::vector<Bembel::ElementTreeNode>::const_iterator mesh_it = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();

        //number of multipole coefficients
        int num_coeffs = (num_multipoles_+1)*(num_multipoles_+1);

        //product of Moments and rhs
        Eigen::VectorXcd Mv(num_coeffs);

        //evaluation matrix
        Eigen::VectorXcd S_lm(num_coeffs);

        //to get start and end intervals
        int k0,k1;

        //distance vector
        Eigen::Vector3d d;

        //spherical coordinates
        double r,t,p;

        for(int l = max_tree_lvl_+1; l >= 0 ; --l){


          //std::cout << "level = " << l << std::endl;

          std::vector<Bembel::ElementOctTreeNode>::iterator el_it;

          for(el_it = el_tree_mem->lbegin(l); el_it != el_tree_mem->lend(l) ; el_it++)
          {

            //moments were computed for the leafs only (#sons = 0)
            if((*el_it).sons_.size() == 0){

              //running index
              int k = 0;


              Mv =  (*el_it).moment_mat_ * density;


            for(int i = 0 ; i < num_eval; ++i){

              S_lm = Bembel::Slm_alt(num_multipoles_,  pos.row(i), (*el_it).center_);//Bembel::Slm(num_multipoles_,  pos.row(i), (*el_it).center_);


              Eigen::VectorXcd tmp = S_lm.transpose()*Mv;

              ret_vec(i) += tmp(0).real();

            }
            }
          }
        }

        return ret_vec;///4./BEMBEL_PI;
      }

  */
  Eigen::VectorXd evaluate_near_field(const Eigen::VectorXd &rhs){

    //std::cout << "H_near = " << H_near_.rows() << " x " << H_near_.cols() << std::endl;
    //std::cout << H_near_ << std::endl;
    //Eigen::VectorXd my_result

    Eigen::VectorXd result = H_near_*rhs;

    return result;


  }

  Eigen::VectorXd evaluate_near_field_transpose(const Eigen::VectorXd &rhs){

    //std::cout << "H_near = " << H_near_.rows() << " x " << H_near_.cols() << std::endl;
    //std::cout << H_near_ << std::endl;
    //Eigen::VectorXd my_result
    Eigen::VectorXd result = H_near_.transpose()*rhs;

    return result;


  }
  /*
    void test_yourself(){

      //this test case is used to validate the multipole translations

      int max_tree_lvl = 2;


      //we define some bounding box
      Eigen::Matrix<double,2,3> bbox;
      bbox <<  -1.1, -1.1, -1.1,
                1.1, 1.1,  1.1;


      //center of source multipole expansion
      Eigen::Vector3d y_s, y_f, x_s, x_f;
      y_s  << -0.9625,
              -0.9625,
              -0.9625;

      //center of source father cell
      y_f  << -0.825,
              -0.825,
              -0.825;

      //center of local father cell
      x_f  << 0.825,
              0.825,
              0.825;

      //center of target cell
      x_s  << 0.9625,
              0.9625,
              0.9625;

      //we define some points around x_s for verification
      Eigen::MatrixXd eval_pos(10,3);
      eval_pos.col(0) = Eigen::VectorXd::LinSpaced(10,x_s(0)+0.125,x_s(0)-0.125);
      eval_pos.col(1) = Eigen::VectorXd::LinSpaced(10,x_s(1)-0.125,x_s(1)+0.125);
      eval_pos.col(2) = Eigen::VectorXd::LinSpaced(10,x_s(2)-0.125,x_s(2)+0.125);

      //multipoles at source
      Eigen::VectorXcd O((num_multipoles_+1)*(num_multipoles_+1));
      O.setZero();

      O(9)  = -0.19764235;
      O(10) =  0.48412292;
      O(11) = -0.45927933;
      O(12) = -0.1767767;
      O(13) = -0.45927933;
      O(14) =  0.48412292;
      O(15) =  -0.19764235;

      std::cout << "multipole moment container for m2m translations" << std::endl;
      Bembel::MultipoleMomentContainer m2m_container;
      m2m_container.fill(num_multipoles_, max_tree_lvl, bbox);

      std::cout << "multipoles at source father cell" << std::endl;
      Eigen::VectorXcd M = m2m_container.m2m(2, y_s-y_f, O);

      std::cout << "rotation matrix for ml2 transformation" << std::endl;
      Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> rot_mat = make_m2l_rot(y_f,x_f);

      std::cout << "rotate M" << std::endl;
      Eigen::VectorXcd M_r = rot_mat * M;

      std::cout << "transfer matrix for ml2 transformation" << std::endl;
      Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> transl_mat = make_m2l_tansl((y_f-x_f).norm());

      std::cout << "translate multipoles" << std::endl;
      Eigen::VectorXcd L_r = transl_mat * M_r;

      std::cout << "rotate back" << std::endl;
      Eigen::VectorXcd L_f = rot_mat.adjoint()  * L_r;

      std::cout << "multipole moment container for l2l translations" << std::endl;
      Bembel::MultipoleMomentContainer l2l_container;
      l2l_container.fill(num_multipoles_, max_tree_lvl, bbox, Bembel::MultipoleMomentContainer::LOCAL);

      std::cout << "locals at son cell" << std::endl;
      Eigen::VectorXcd L = l2l_container.l2l(2, x_f-x_s, L_f);

      //compute the potential due to these sources
      Eigen::VectorXd pot_gt = Bembel::evaluate_multipole_expansion(eval_pos, y_s, O);
      //compute the potential due to these sources
      Eigen::VectorXd pot_m2m = Bembel::evaluate_multipole_expansion(eval_pos, y_f, M);
      //compute the potential due to these sources
      Eigen::VectorXd pot_m2l = Bembel::evaluate_local_expansion(eval_pos, x_f, L_f,true);
      //compute the potential due to these sources
      Eigen::VectorXd pot_l2l = Bembel::evaluate_local_expansion(eval_pos, x_s, L);

      std::cout << "Pot GT = " << pot_gt.transpose() << std::endl;
      std::cout << "Pot m2m = "<< pot_m2m.transpose() << std::endl;
      std::cout << "Pot m2l = "<< pot_m2l.transpose() << std::endl;
      std::cout << "Pot l2l = "<< pot_l2l.transpose() << std::endl;

      Eigen::MatrixXd out_mat(eval_pos.rows(),7);
      out_mat.block(0,0,eval_pos.rows(),3) = eval_pos;
      out_mat.block(0,3,eval_pos.rows(),1) = pot_gt;
      out_mat.block(0,4,eval_pos.rows(),1) = pot_m2m;
      out_mat.block(0,5,eval_pos.rows(),1) = pot_m2l;
      out_mat.block(0,6,eval_pos.rows(),1) = pot_l2l;

      // define the format you want, you only need one instance of this...
      const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

      std::ofstream out_file_eval("../../output/thesis/fmm/accelerated/test_yourself.csv");
      out_file_eval << out_mat.format(CSVFormat);
      out_file_eval.close();
      }

  */

  void set_min_numel(const int min_numel){

    min_numel_ = min_numel;

  }




  Eigen::MatrixXcd get_moment_matrix(const int indx){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    Bembel::ElementOctTreeNode* el_node = el_tree_mem->get_element_ptr(source_leafs_[indx]);

    return el_node->moment_mat_;

  }

  Eigen::SparseMatrix<double,Eigen::RowMajor> get_H_near(){

    return H_near_;

  }

  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  // we declare functionality which has not been implemented (yet)
  // to be private
  MLFMMMatrixDevice(const MLFMMMatrixDevice<Derived,LinOp>& FMM);
  MLFMMMatrixDevice(MLFMMMatrixDevice<Derived,LinOp>&& FMM);
  MLFMMMatrixDevice& operator=(const MLFMMMatrixDevice<Derived,LinOp>& FMM);
  MLFMMMatrixDevice& operator=(MLFMMMatrixDevice<Derived,LinOp>&& FMM);

  //potential operator
  Bembel::HallProbeArray<Derived,LinOp> *device_;

  //Sparse near field matrix
  Eigen::SparseMatrix<double,Eigen::RowMajor> H_near_;


  Index rows_;
  Index cols_;
  Index cols_disc_;

  int num_multipoles_ = 5;

  //minumum tree level, at this level, all interactions are computed as ml2
  int min_tree_level_ = 1;

  //number of near field interactions
  int num_near_field_interactions_ = 0;

  Bembel::ElementOctTree<LinOp> el_tree_;
  Bembel::MeasurementData *meas_data_;

  int max_tree_lvl_;
  int min_numel_ = 10;
  int num_coeffs_;
  int I_;
  int I2_;
  int num_sensors_;

  //container for m2m matrices
  Bembel::MultipoleMomentContainer m2m_matrices_;
  //container for l2l matrices
  Bembel::MultipoleMomentContainer l2l_matrices_;

  //container for m2l matrices. This is a unique list. The measurement and element trees
  //store the indices for the transformation matrix in this list
  //std::vector<Eigen::MatrixXcd,Eigen::aligned_allocator<Eigen::MatrixXcd>> m2l_matrices_;
  std::vector<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>,
              Eigen::aligned_allocator<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>>> m2l_matrices_;

  std::vector<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>,
              Eigen::aligned_allocator<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>>> m2l_tansl_;

  std::vector<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>,
              Eigen::aligned_allocator<Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor>>> m2l_rot_;

  //difference vectors corresponding to the m2l transfer in the unique list
  std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> d_vecs_;
  //difference vectors corresponding to the m2l rotations in the unique list
  std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> d_rot_vecs_;

  //quad degree
  int deg_ ;

  //distance from center to corner of box for level 0
  double diam_;

  //cell size
  double hx_,hy_,hz_;

  //diagonal of box
  Eigen::Vector3d diag_;

  //near field interaction list (el_indx , meas_indx)
  std::vector<std::vector<int>> near_field_interaction_list_;

  Bembel::AnsatzSpace<LinOp> *ansatz_space_;

  //discrete Potential
  Derived *potential_;

  std::vector<int> local_leafs_;
  std::vector<int> source_leafs_;

  //Eigen::MatrixXd *meas_pos_;

  /*
  double count_m2l_ram(){

    int count = 0;

    //rotation matrices
    for(int i = 0 ; i < m2l_rot_.size() ; ++i){

      count += m2l_rot_[i].nonZeros();
    }
    //transfer matrices
    for(int i = 0 ; i < m2l_tansl_.size() ; ++i){

      count += m2l_tansl_[i].nonZeros();
    }

    return count * 8e-9;


  }

  void m2m(Eigen::VectorXcd *M_p, const Eigen::VectorXcd &M, const Eigen::SparseMatrix<std::complex<double>> &R){

    (*M_p) += R*M;

  }
  */
  void compute_box_diam(const Eigen::Matrix<double,2,3> &bbox){

    diag_ = (bbox.row(1)-bbox.row(0)).transpose();
    diam_ = diag_.norm();

    hx_ = bbox(1,0) - bbox(0,0);
    hy_ = bbox(1,1) - bbox(0,1);
    hz_ = bbox(1,2) - bbox(0,2);

  }

  /*



  //moving from x to x_p
  Eigen::SparseMatrix<std::complex<double>> make_m2m_trafo(const Eigen::Vector3d &x_p, const Eigen::Vector3d &x){

    Eigen::VectorXcd R = Bembel::Rlm(num_multipoles_, x-x_p);

    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tripletList;
    tripletList.reserve(num_multipoles_*num_multipoles_/2);


    //running indices
    int row = 0;
    int k_p;

    for(int l = 0; l <= num_multipoles_; ++l){
      for(int m = -l; m <= l ; ++m){

        int k = 0;
        for(int l_p = 0; l_p <= l; ++l_p){
          for(int m_p = -l_p; m_p <= l_p ; ++m_p){


            k_p = make_k_p(l,m,l_p,m_p);
            if(k_p > -1){

              //std::cout << "( " << row << " , "  << k_p << " , " << R(k) << " )" <<std::endl;
              tripletList.push_back(T(row , k_p , R(k) ));
            }
            ++k;

          }
        }

        ++row;
      }
    }

    int num_coeffs = (num_multipoles_+1)*(num_multipoles_+1);

    Eigen::SparseMatrix<std::complex<double>> mat(num_coeffs,num_coeffs);

    mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return mat;


  }

  int make_k_p(const int l,const int m,const int l_p,const int m_p){

    //bottom index
    int b = l-l_p;
    if(b < 0){
      return -1;
    }

    //top index
    int t = m-m_p;

    if(std::abs(t) > b){
      return -1;
    }

    //k = l**2 + l + m
    return b*b + b + t;




  }


  Eigen::VectorXcd m2l(const Eigen::Vector3d &x_l , const Eigen::Vector3d &x_c, const Eigen::VectorXcd &moments){

    //return vector
    Eigen::VectorXcd ret_vec(num_coeffs_);
    ret_vec.setZero();


    //irregular solid harmonics
    Eigen::VectorXcd S_lm = Bembel::Slm(num_multipoles_,  x_l, x_c);
    //Eigen::VectorXcd S_lm(moments.size());
    //running indices
    int k;
    int k_s = 0;
    int k_p;

    double fac;

    for(int l = 0; l <= num_multipoles_; ++l){

      fac = std::pow(-1.,l);

      for(int m = -l; m <= l; ++m){

        int k_p = 0;

        for(int lp = 0; lp <= num_multipoles_; ++lp){

          for(int mp = -lp; mp <= lp; ++mp){


            k = make_k_m2l( l, m, lp, mp);


            if(k > -1){

              ret_vec(k_s) += fac*std::conj(S_lm(k))*moments(k_p);

            }
            ++k_p;
          }
        }
        ++k_s;
      }
    }
    return ret_vec;

  }
  */
  void make_ir_and_nf(){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //distance vector
    Eigen::Vector3d d_vec;
    //distance vector of fathers
    Eigen::Vector3d d_vec_f;

    //reset near field interaction count
    num_near_field_interactions_ = 0;

    int I2 = (ansatz_space_->get_polynomial_degree() + 1)*(ansatz_space_->get_polynomial_degree() + 1);

    //distance
    double d;

    //father cells
    Bembel::ElementOctTreeNode el_father;
    Bembel::MeasurementTreeNode meas_father;

    for(int l = 2; l <= max_tree_lvl_+1 ; ++l){

      //std::cout << "level = " << l << std::endl;

      std::vector<Bembel::ElementOctTreeNode>::iterator el_it;
      std::vector<Bembel::MeasurementTreeNode>::iterator meas_it;

      for(el_it = el_tree_mem->lbegin(l); el_it != el_tree_mem->lend(l) ; el_it++)
      {


        el_father = el_tree_mem->get_element((*el_it).father_);

        for(meas_it = meas_tree_mem->lbegin(l); meas_it != meas_tree_mem->lend(l) ; meas_it++){


              meas_father = meas_tree_mem->get_element((*meas_it).father_);

              int c = compare_cells(l, (*el_it).center_ ,(*meas_it).center_,el_father.center_,meas_father.center_);


              if(c == 0){
                  //dense operation

                  //only leafs can interact as near field.
                  if((el_it->sons_.size() == 0) || (meas_it->sons_.size() == 0)){

                    (*meas_it).near_field_.push_back((*el_it).pos_);
                    (*el_it).near_field_.push_back((*meas_it).pos_);

                    near_field_interaction_list_.push_back({(*el_it).pos_,(*meas_it).pos_});

                    num_near_field_interactions_ += (*meas_it).indices_.size() * (*el_it).indices_.size() * I2;
                  }
              }
              else if( c == 1 ){
                //compressed operation

                (*meas_it).interaction_region_.push_back((*el_it).pos_);
                (*el_it).interaction_region_.push_back((*meas_it).pos_);

                //Eigen::MatrixXcd tmp = make_m2l_matrix(meas_it->center_ - el_it->center_);

                //source - local
                //we are searching for the m2l matrix based on the difference vector
                Eigen::Vector3d d_tmp = el_it->center_ - meas_it->center_;

                //find rotation matrix
                int rot_indx = find_rot_matrix(d_tmp);

                //if not found, we append this rotation to the unique list
                if(rot_indx == -1){
                    m2l_rot_.push_back(make_m2l_rot(el_it->center_,meas_it->center_));
                    d_rot_vecs_.push_back(d_tmp);

                    rot_indx = m2l_rot_.size() - 1;
                }

                //find transfer matrix
                int m2l_indx = find_m2l_matrix(d_tmp);

                //if not found, we append this rotation to the unique list
                if(m2l_indx == -1){
                    m2l_tansl_.push_back(make_m2l_tansl(d_tmp.norm()));
                    d_vecs_.push_back(d_tmp);

                    m2l_indx = m2l_tansl_.size() - 1;
                }

                (*meas_it).interaction_m2l_index_.push_back(m2l_indx);
                (*meas_it).interaction_m2l_rot_index_.push_back(rot_indx);

                (*el_it).interaction_m2l_index_.push_back(m2l_indx);
                (*el_it).interaction_m2l_rot_index_.push_back(rot_indx);

              }
        }
      }
    }
  }
  /*
  int make_k_m2l(const int l,const int m,const int l_p,const int m_p){

    //bottom index
    int b = l+l_p;
    if(b > num_multipoles_){
      return -1;
    }

    //top index
    int t = m+m_p;

    if(std::abs(t) > b){
      return -1;
    }

    //k = l**2 + l + m
    return b*b + b + t;


  }

  */

  void upward_pass(){

      //element tree memory
      std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

      //std::cout << "upward pass" << std::endl;

      //we start from the second deepest level
      //we iterate up until level 1, this level is twice bisected and contains: 64
      //cells.
      for(int l = max_tree_lvl_; l > 1 ; --l){

        //std::cout << "level = " << l << std::endl;

  #pragma omp parallel for
        for(std::vector<Bembel::ElementOctTreeNode>::iterator  el_it = el_tree_mem->lbegin(l); el_it < el_tree_mem->lend(l) ; ++el_it){

          //std::cout << "Cell center = " << (*el_it).center_ << std::endl;

          //transformation is only performed to father cells
          if(el_it->sons_.size() > 0){

            //std::cout << "\tis father" << std::endl;

            //pointer to son
            Bembel::ElementOctTreeNode *son;

            //make space for moments in father
            el_it->moments_.resize(num_coeffs_);
            el_it->moments_.setZero();

            for(int i = 0; i < el_it->sons_.size(); ++i){


              son = el_tree_mem->get_element_ptr(el_it->sons_[i]);


              m2m_matrices_.m2m(l, son->center_ - el_it->center_, son->moments_, &el_it->moments_);

            }

          }
        }
      }
  }


  void compute_leaf_moments(const Eigen::VectorXd &rhs){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //mesh tree
    //std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = potential_->get_ansatz_space()->get_superspace().get_mesh().get_element_tree().root().memory_;

    //here we run upwards, and look for the leafs in the element tree
    //the first level in the tree is labeled with -1. max_tree_lvl = 2 means that
    //there are the levels -1 , 0 , 1 and 2 so four levels. The last level is
    //at position max_tree_lvl_ +1 in the memory.
    for(int l = max_tree_lvl_+1; l > 1 ; --l){

    //std::cout << "level = " << l << std::endl;
    #pragma omp parallel
    {

    #pragma omp for
      for(std::vector<Bembel::ElementOctTreeNode>::iterator  el_it = el_tree_mem->lbegin(l); el_it < el_tree_mem->lend(l) ; el_it++)
      {


        //moments are computed for the leafs only (#sons = 0)
        if(el_it->sons_.size() == 0){

          //running index
          int k = 0;

          //for start and end in mesh memory
          int k0,k1;

          //we could do this also in the matrix init
          el_it->moments_.resize(num_coeffs_);

          el_it->moments_ =  el_it->moment_mat_ * rhs;


        }
        }
      }
    }
  }


  void compute_leaf_moments_transpose(const Eigen::VectorXd &rhs){

    #pragma omp parallel
    {

      //meas tree memory
      std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();


      //measurement tree node
      Bembel::MeasurementTreeNode *meas_node;

      #pragma omp for
      for(int i = 0; i < local_leafs_.size(); ++i){

        meas_node = meas_tree_mem->get_element_ptr(local_leafs_[i]);

        //std::cout << "meas node center = " << meas_node->center_ << std::endl;
        //std::cout << "meas_node->local_mat_ = " << meas_node->local_mat_.rows() << " x " << meas_node->local_mat_.cols() << std::endl;
        //std::cout << "rhs = " << rhs.rows() << " x " << rhs.cols() << std::endl;
        meas_node->locals_ = meas_node->local_mat_.adjoint()*rhs;

      }

    }
  }

  void upward_pass_transpose(){

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //we start from the second deepest level
    //we iterate up until level 1, this level is twice bisected and contains: 64
    //cells.
    for(int l = max_tree_lvl_; l > 1 ; --l){

      #pragma omp parallel for
        for(std::vector<Bembel::MeasurementTreeNode>::iterator  meas_it = meas_tree_mem->lbegin(l); meas_it < meas_tree_mem->lend(l) ; ++meas_it){

          //std::cout << "Cell center = " << (*meas_it).center_ << std::endl;

          //transformation is only performed to father cells
          if(meas_it->sons_.size() > 0){

            //pointer to son
            Bembel::MeasurementTreeNode *son;

            //make space for locals in father
            meas_it->locals_.resize(num_coeffs_);
            meas_it->locals_.setZero();

            for(int i = 0; i < meas_it->sons_.size(); ++i){

              son = meas_tree_mem->get_element_ptr(meas_it->sons_[i]);

              //adjoint should be equal to swapping the positions
              l2l_matrices_.l2l(l, meas_it->center_ - son->center_ , son->locals_, &meas_it->locals_,true);

            }

        }
      }
    }
  }

  void downward_pass_transpose(){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //std::cout << "max_tree_lvl_ = " << max_tree_lvl_ << std::endl;
    //we start from the second level. This level is at memory position 2
    for(int l = 2; l <= max_tree_lvl_+1 ; ++l){

    //std::cout << "level = " << l << std::endl;
    #pragma omp parallel for
    for(std::vector<Bembel::ElementOctTreeNode>::iterator  el_it = el_tree_mem->lbegin(l); el_it < el_tree_mem->lend(l) ; ++el_it){


        //std::cout << "Meas center = " << (*meas_it).center_ << std::endl;

        //pointer to element cell
        Bembel::MeasurementTreeNode *meas_cell;

        //make space for locals
        el_it->moments_.resize(num_coeffs_);
        //set locals to zero
        el_it->moments_.setZero();

        //std::cout << "downward_pass" << std::endl;


        //m2l transformations in the interaction region
        for(int i = 0; i < el_it->interaction_region_.size(); ++i){

          meas_cell = meas_tree_mem->get_element_ptr(el_it->interaction_region_[i]);

          //std::cout << "Interacting source center = " << (*el_cell).center_<< std::endl;

          //compute harmonics
          //meas_it->locals_ += m2l(meas_it->center_,el_cell->center_,el_cell->moments_);

          //based on stored matrices
          //meas_it->locals_ += m2l_matrices_[meas_it->interaction_m2l_index_[i]]*el_cell->moments_;

          //std::cout << "m2l_tansl_.size() = " << m2l_tansl_.size() << std::endl;
          //std::cout << "m2l_rot_.size() = " << m2l_rot_.size() << std::endl;
          //std::cout << "rot index = " << meas_it->interaction_m2l_rot_index_[i] << std::endl;
          //std::cout << "transl index = " << meas_it->interaction_m2l_index_[i] << std::endl;

          //rotate
          Eigen::VectorXcd tmp = m2l_rot_[el_it->interaction_m2l_rot_index_[i]] * meas_cell->locals_;

          tmp = (m2l_tansl_[el_it->interaction_m2l_index_[i]].adjoint() * tmp).eval();

          //rotate back
          el_it->moments_ += (m2l_rot_[el_it->interaction_m2l_rot_index_[i]].adjoint() * tmp).eval();


        }
        if(l > 2){

        //std::cout << "l2l" << std::endl;
        //l2l transformations for the far field

        //std::cout << "\tl2l" << std::endl;


        //pointer to father cell
        Bembel::ElementOctTreeNode *father = el_tree_mem->get_element_ptr(el_it->father_);

        //std::cout << "\t\tfather center = " << father->center_.transpose() << std::endl;
        //std::cout << "\t\tson center = " << meas_it->center_.transpose() << std::endl;

        //l2l transformation
        //adjoint should be equal to swapping the positions
        m2m_matrices_.m2m(l-1, el_it->center_ - father->center_, father->moments_, &el_it->moments_,true);

        }
      }
    }
  }

  double count_moment_matrix_RAM(){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //mesh tree
    //std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = potential_->get_ansatz_space()->get_superspace().get_mesh().get_element_tree().root().memory_;

    //memory in GB
    double ram;

    //here we run upwards, and look for the leafs in the element tree
    //the first level in the tree is labeled with -1. max_tree_lvl = 2 means that
    //there are the levels -1 , 0 , 1 and 2 so four levels. The last level is
    //at position max_tree_lvl_ +1 in the memory.
    for(int l = max_tree_lvl_+1; l > 1 ; --l){

    //std::cout << "level = " << l << std::endl;
    //#pragma omp parallel
    {

      //count
      int my_count = 0;

      //mesh iterator
      std::vector<Bembel::ElementTreeNode>::const_iterator mesh_it = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();

    //#pragma omp for
      for(std::vector<Bembel::ElementOctTreeNode>::iterator  el_it = el_tree_mem->lbegin(l); el_it < el_tree_mem->lend(l) ; el_it++)
      {

        //moments are computed for the leafs only (#sons = 0)
        if(el_it->sons_.size() == 0){

          my_count += el_it->moment_mat_.nonZeros();

          }
        }

        //#pragma omp critical
        {
          ram += my_count * 8e-9;
        }
      }
    }
    return ram;
  }

  void downward_pass(){

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //std::cout << "max_tree_lvl_ = " << max_tree_lvl_ << std::endl;
    //we start from the second level. This level is at memory position 2
    for(int l = 2; l <= max_tree_lvl_+1 ; ++l){

    //std::cout << "level = " << l << std::endl;
    #pragma omp parallel for
      for(std::vector<Bembel::MeasurementTreeNode>::iterator meas_it = meas_tree_mem->lbegin(l); meas_it < meas_tree_mem->lend(l) ; meas_it++){


        //std::cout << "Meas center = " << (*meas_it).center_ << std::endl;

        //pointer to element cell
        Bembel::ElementOctTreeNode *el_cell;

        //make space for locals
        meas_it->locals_.resize(num_coeffs_);
        //set locals to zero
        meas_it->locals_.setZero();


        //m2l transformations in the interaction region
        for(int i = 0; i < meas_it->interaction_region_.size(); ++i){




          el_cell = el_tree_mem->get_element_ptr(meas_it->interaction_region_[i]);

          //std::cout << "Interacting source center = " << (*el_cell).center_<< std::endl;

          //compute harmonics
          //meas_it->locals_ += m2l(meas_it->center_,el_cell->center_,el_cell->moments_);

          //based on stored matrices
          //meas_it->locals_ += m2l_matrices_[meas_it->interaction_m2l_index_[i]]*el_cell->moments_;

          //std::cout << "m2l_tansl_.size() = " << m2l_tansl_.size() << std::endl;
          //std::cout << "m2l_rot_.size() = " << m2l_rot_.size() << std::endl;
          //std::cout << "rot index = " << meas_it->interaction_m2l_rot_index_[i] << std::endl;
          //std::cout << "transl index = " << meas_it->interaction_m2l_index_[i] << std::endl;

          //rotate
          Eigen::VectorXcd tmp = m2l_rot_[meas_it->interaction_m2l_rot_index_[i]] * el_cell->moments_;

          tmp = (m2l_tansl_[meas_it->interaction_m2l_index_[i]] * tmp).eval();

          //rotate back
          meas_it->locals_ += (m2l_rot_[meas_it->interaction_m2l_rot_index_[i]].adjoint() * tmp).eval();


        }
        if(l > 2){


          //l2l transformations for the far field

          //std::cout << "\tl2l" << std::endl;


          //pointer to father cell
          Bembel::MeasurementTreeNode *father = meas_tree_mem->get_element_ptr(meas_it->father_);


          if(father->locals_.norm() > 1e-12){
            //l2l transformation
            l2l_matrices_.l2l(l-1, father->center_ - meas_it->center_, father->locals_, &meas_it->locals_);
          }
        }
      }
    }
  }
  /*
  //Eigen::MatrixXcd make_m2l_matrix(const Eigen::Vector3d &d){
  Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> make_m2l_matrix(const Eigen::Vector3d &d){
    //Eigen::MatrixXcd ret_mat(num_coeffs_,num_coeffs_);

    Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> ret_mat;
    ret_mat.resize(num_coeffs_,num_coeffs_);
    ret_mat.setZero();

    //irregular solid harmonics
    Eigen::VectorXcd S_lm = Bembel::Slm(num_multipoles_, d);

    //running indices
    int k;
    int k_s = 0;
    int k_p;

    double fac;


    //tripletList to construct sparse matrix
    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tripletList;
    tripletList.reserve(num_coeffs_*num_coeffs_);

    for(int l = 0; l <= num_multipoles_; ++l){

      fac = std::pow(-1.,l);

      for(int m = -l; m <= l; ++m){

        int k_p = 0;

        for(int lp = 0; lp <= num_multipoles_; ++lp){

          for(int mp = -lp; mp <= lp; ++mp){


            k = make_k_m2l( l, m, lp, mp);


            if(k > -1){

              tripletList.push_back(T(k_s ,k_p, fac*std::conj(S_lm(k))));
              //ret_mat(k_s,k_p) += fac*std::conj(S_lm(k));

            }
            ++k_p;
          }
        }
        ++k_s;
      }
    }

    //std::cout << "make sparse matrix" << std::endl;
    ret_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << ret_mat << std::endl;
    return ret_mat;
  }
  */
  //Make a moment to local translation matrix for translate along z
  Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> make_m2l_tansl(const double dist){

    Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> ret_mat;
    ret_mat.resize(num_coeffs_,num_coeffs_);
    ret_mat.setZero();

    int row;
    int col;

    std::complex<double> num;
    std::complex<double> den;

    //imaginary number as helper
    std::complex<double> I(0.0,1.0);

    //tripletList to construct sparse matrix
    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tripletList;
    tripletList.reserve(num_coeffs_*num_coeffs_);

    for(int j = 0; j <= num_multipoles_; ++j){
      for(int k = -j; k <= j; ++k){

        for(int n = 0; n <= num_multipoles_; ++n){
          for(int m = -n; m <= n ; ++m){

            if(m == k){

              row = j*j + j + k;
              col = n*n + n + m;

              num = std::pow(-1.,-std::abs(m))*Bembel::A_nm(n,m)*Bembel::A_nm(j,k)*Bembel::Ynm_alt(0.,0.,j+n,0);
              den = std::pow(-1.,n)* Bembel::A_nm(j+n,0)*std::pow(dist,j+n+1);

              tripletList.push_back(T(row ,col, num/den ));

            }
          }
        }
      }
    }


    //std::cout << "make sparse matrix" << std::endl;
    ret_mat.setFromTriplets(tripletList.begin(), tripletList.end());
    //std::cout << ret_mat << std::endl;
    return ret_mat;
  }

  //Make a moment to local rotation matrix
  Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> make_m2l_rot(const Eigen::Vector3d &x_s,const Eigen::Vector3d &x_l){

    Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> ret_mat;
    ret_mat.resize(num_coeffs_,num_coeffs_);
    ret_mat.setZero();

    //source - local
    Eigen::Vector3d d = x_s - x_l;
    double r = d.norm();

    //algles
    double alpha = std::atan2(d(1),d(0));
    double beta = std::acos(d(2)/r);

    ret_mat =  Bembel::make_rotation_matrix(num_multipoles_,beta,alpha);

    return ret_mat;
  }

  int compare_cells(const int level, const Eigen::Vector3d &xs_1,const Eigen::Vector3d &xs_2,
                                      const Eigen::Vector3d &xf_1,const Eigen::Vector3d &xf_2){

    Eigen::Vector3d diff_s = xs_1 - xs_2;
    Eigen::Vector3d diff_f = xf_1 - xf_2;




    if((std::fabs(diff_s(0)) <= hx_ / std::pow(2,level)  * 1.001 )
        && (std::fabs(diff_s(1)) <= hy_ / std::pow(2,level) * 1.001 )
        && (std::fabs(diff_s(2)) <= hz_ / std::pow(2,level) * 1.001) ){

          //std::cout << "found adjacent cells" << std::endl;
          //adjacent
          return 0;

        }
    else if(level == 2){

          //std::cout << "all other cells on level 2 are interacting" << std::endl;
          //interaction region, all are interacting in level 2
          return 1;
    }
    else if((std::fabs(diff_f(0)) <= hx_ / std::pow(2,level-1) * 1.001)
        && (std::fabs(diff_f(1)) <= hy_ / std::pow(2,level-1) * 1.001 )
        && (std::fabs(diff_f(2)) <= hz_ / std::pow(2,level-1) * 1.001 )){

      //std::cout << "found interacting cells" << std::endl;
      //interaction region
      return 1;

    }
    else{

      //std::cout << "this is far field" << std::endl;
      //far field
      return 2;

    }

  }


  Eigen::VectorXd evaluate_far_field(){

      Eigen::VectorXd result(meas_data_->get_number_of_measurements() * num_sensors_);
      result.setZero();

      #pragma omp parallel
      {

      //meas tree memory
      std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();


      //measurement tree node
      Bembel::MeasurementTreeNode *meas_node;

      Eigen::VectorXd my_result(meas_data_->get_number_of_measurements() * num_sensors_ );
      my_result.setZero();


      #pragma omp for
      for(int i = 0; i < local_leafs_.size(); ++i){

        meas_node = meas_tree_mem->get_element_ptr(local_leafs_[i]);

        //std::cout << "meas node center = " << meas_node->center_ << std::endl;

        //Eigen::VectorXcd matvec = meas_node->local_mat_*meas_node->locals_;

        my_result += (meas_node->local_mat_*meas_node->locals_).real();



        //for(int j = 0; j < meas_node->indices_.size(); ++j){

        //  my_result(meas_node->indices_[j]) += matvec(j).real();

        //}

      }

      #pragma omp critical
      {
        result += my_result;
      }

    }

    //std::cout << "Far field norm = " << result.norm() << std::endl;
    return result;///4 ./BEMBEL_PI;


  }

  Eigen::VectorXd evaluate_far_field_transpose(){

      Eigen::VectorXd result(cols_);
      result.setZero();

      //std::cout << "cols = " << cols_ << std::endl;
      //element tree memory
      std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

      //element tree node
      Bembel::ElementOctTreeNode *el_node;

      Eigen::VectorXd my_result(cols_);
      my_result.setZero();


      #pragma omp for
      for(int i = 0; i < source_leafs_.size(); ++i){

        el_node = el_tree_mem->get_element_ptr(source_leafs_[i]);

        //std::cout << "meas node center = " << meas_node->center_ << std::endl;

        //Eigen::VectorXcd matvec = meas_node->local_mat_*meas_node->locals_;
        //std::cout << el_node->moment_mat_.rows() << " x " << el_node->moment_mat_.cols() << std::endl;
        //std::cout << el_node->moments_.rows() << " x " << el_node->moments_.cols() << std::endl;

        my_result += (el_node->moment_mat_.adjoint()*el_node->moments_).real();

        //std::cout << my_result.rows() << " x " << my_result.cols() << std::endl;

        //for(int j = 0; j < meas_node->indices_.size(); ++j){

        //  my_result(meas_node->indices_[j]) += matvec(j).real();

        //}

      }

      #pragma omp critical
      {
        result += my_result;
      }

    return result;///4 ./BEMBEL_PI;


  }
  /*
  double count_locals_RAM(){

    double ram = 0.;

    //#pragma omp parallel
    {

      int my_count = 0;

      //meas tree memory
      std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

      //measurement tree node
      Bembel::MeasurementTreeNode *meas_node;

      //#pragma omp for
      for(int i = 0; i < local_leafs_.size(); ++i){

        meas_node = meas_tree_mem->get_element_ptr(local_leafs_[i]);
        my_count += meas_node->local_mat_.nonZeros(); //meas_node->local_mat_.rows()*meas_node->local_mat_.cols();

      }
      //#pragma omp critical
      {
        ram += my_count*8e-9;
      }

    }
    return ram;
  }
  */
  int find_m2l_matrix(const Eigen::Vector3d &d){

    int ret_val = -1;

    double dist = d.norm();


    for(int i = 0; i < d_vecs_.size(); ++i){


      if( std::abs(d_vecs_[i].norm() - dist) < 1e-10 ) {

        ret_val = i;
        break;
      }
    }
    return ret_val;
  }

  //source - local
  int find_rot_matrix(const Eigen::Vector3d &d){

    int ret_val = -1;

    //std::cout << "find_rot_matrix" << std::endl;

    for(int i = 0; i < d_rot_vecs_.size(); ++i){

      //std::cout << "d_rot = " << d_rot_vecs_[i].transpose() << std::endl;
      //std::cout << "d = " << d.transpose() << std::endl;

      double proj = d_rot_vecs_[i].transpose() * d;

      //std::cout << "proj = " << proj << std::endl;

      if( std::abs(proj/d.norm()/d_rot_vecs_[i].norm()  - 1.) < 1e-10){

        ret_val = i;
        break;
      }
    }
    return ret_val;
  }
  /*

  int count_nonzeros(){

    int count = 0;

    //element tree memory
    std::shared_ptr<Bembel::ElementOctTreeMemory> el_tree_mem = el_tree_.get_tree_memory();

    //meas tree memory
    std::shared_ptr<Bembel::MeasurementTreeMemory> meas_tree_mem = meas_data_->get_tree_memory();

    //pointer to element cluster memory. With this pointer we access the elements for integration
    std::vector<Bembel::ElementTreeNode>::const_iterator mesh_it = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();

    //local and source indices
    int i_l,i_s;

    int k0,k1;

    //mesh tree
    std::shared_ptr<Bembel::ElementTreeMemory> mesh_memory = potential_->get_ansatz_space()->get_superspace().get_mesh().get_element_tree().root().memory_;

    //#pragma omp for
    for(int i = 0; i < near_field_interaction_list_.size() ; ++i){

      //indices of interaction
      i_l = near_field_interaction_list_[i][1];
      i_s = near_field_interaction_list_[i][0];

      //interacting cells
      Bembel::MeasurementTreeNode *meas_node = meas_tree_mem->get_element_ptr(i_l);
      Bembel::ElementOctTreeNode *source_node = el_tree_mem->get_element_ptr(i_s);

      for(int source_indx : source_node->indices_){
        for(int meas_indx : meas_node->indices_){

          k0 = I2_*std::distance(mesh_memory->cpbegin(),mesh_memory->cluster_begin(mesh_it[source_indx]));
          k1 = I2_*std::distance(mesh_memory->cpbegin(),mesh_memory->cluster_end(mesh_it[source_indx]));

          count += k1-k0;

        }
      }
    }
    return count;
  }

};


}  // namespace Eigen

/**
 * \brief Implementation of H2Matrix * Eigen::DenseVector through a
 * specialization of internal::generic_product_impl
 */
 /*
namespace internal {
template <typename Rhs, typename Device >
struct generic_product_impl<LRMatrix<Device>, Rhs, SparseShape, DenseShape,
                            GemvProduct>  // GEMV stands for matrix-vector
    : generic_product_impl_base<LRMatrix<Device>, Rhs,
                                generic_product_impl<LRMatrix<Device>, Rhs>> {
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const LRMatrix<Device>& lhs,
                            const Rhs& rhs, const double& alpha) {

      assert(alpha == double(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);

      if(lhs.is_squared()){

        dst += lhs.squared_product(rhs);

      }
      else if(lhs.is_transposed()){

        dst += lhs.get_device()->mat_vec_transpose(rhs);
      }
      else{

        //std::cout << "!!!!!!!!!!!!!!" << std::endl;
        dst += lhs.get_device()->mat_vec(rhs);

      }




  }

  */
};


//}  // namespace internal
}  // namespace Eigen

#endif
