// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_ACAMATRIX_H_
#define BEMBEL_ACAMATRIX_H_

#include <Eigen/Dense>


namespace Bembel {

 //template <typename Scalar>
 template<typename Device, typename LinOp>
 class ACAMatrix {

 public:
   //constructor
   ACAMatrix(){


     eta_ = 1.6;
     max_tree_lvl_ = 3;
     min_cluster_level_ = 1;

   }

   void set_parameters(double eta, int max_tree_level, int min_cluster_level, double eps = 1e-4){
     eta_ = eta;
     max_tree_lvl_ = max_tree_level;
     min_cluster_level_ = min_cluster_level;
     eps_ = eps;


   }

   double get_eta(){

     return eta_;
   }

   //template<typename Derived, typename LinOp>
   void init_ACAMatrix(const DiscreteSensor<Device,LinOp> &sensor, MeasurementData *meas_data, AnsatzSpace<LinOp> *ansatz_space){


     //sensor
     sensor_ = sensor;


     //get ansatz_space
     ansatz_space_ = ansatz_space;

     //get super space
     auto super_space = ansatz_space->get_superspace();


     //get polynomail degree
     int polynomial_degree = ansatz_space->get_polynomial_degree();
     int polynomial_degree_plus_one_squared = (polynomial_degree + 1) * (polynomial_degree + 1);


     //quadrature degree
     deg_ = polynomial_degree + 1;

     //Gauss quadrature rule
     Bembel::GaussSquare<Bembel::Constants::maximum_quadrature_degree> GS;
     auto Q = GS[deg_];

     //Get measurement data pointer
     pos_ = meas_data->get_position_ptr();

     //total number of measurements
     num_meas_ = meas_data->get_number_of_measurements();

     //init measurement cluster
     meas_data->init_cluster_tree(max_tree_lvl_);

     //setup block cluster tree
     meas_bct_.set_measurement_data(meas_data);
     meas_bct_.set_parameters(eta_, max_tree_lvl_, min_cluster_level_);
     meas_bct_.setup_element_tree(*ansatz_space);

     //initialize block cluster tree
     meas_bct_.init_BlockClusterTree();

     //how many leaves are there
     int num_leaves = meas_bct_.leaf_pointers_->size();

     //pointer to element cluster memory
     el_memory_ = meas_bct_.get_element_cluster()->memory_;

     #pragma omp parallel for
     for (int l = 0; l < num_leaves; ++l){

       //tmp integers for matrix size
       int tmp_rows, tmp_cols;

       //pointers to the element and measurement clusters
       const ElementTreeNode *this_el;
       const MeasurementTreeNode *this_meas;

       this_el = (*meas_bct_.leaf_pointers_)[l]->get_element_cluster();
       this_meas = (*meas_bct_.leaf_pointers_)[l]->get_measurement_cluster();


       tmp_rows = this_meas->indices_.size();

       tmp_cols = std::distance(el_memory_->cluster_begin(*this_el),
                                      el_memory_->cluster_end(*this_el)) *
                        polynomial_degree_plus_one_squared;

       if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){
         //LOW RANK
         //std::cout << "Low rank" << std::endl;

         //temporal storage for low rank matrix
         Eigen::MatrixXd L, R;

         assemble_low_rank_matrix_block(super_space,*this_el,*this_meas,  Q, &L,&R);

         (*meas_bct_.leaf_pointers_)[l]->leaf_->set_L(L);
         (*meas_bct_.leaf_pointers_)[l]->leaf_->set_R(R);

         (*meas_bct_.leaf_pointers_)[l]->leaf_->set_low_rank_flag(true);



       }
       else{
         //DENSE
         //std::cout << "Dense" << std::endl;

         //temporal storage for dense matrices
         Eigen::MatrixXd F;

         //assemble dense matrix block
         assemble_dense_matrix_block(super_space,*this_el,*this_meas,Q,&F); //, this_el, this_meas, GS, &F
         (*meas_bct_.leaf_pointers_)[l]->leaf_->set_F(F);

         (*meas_bct_.leaf_pointers_)[l]->leaf_->set_low_rank_flag(false);

       }

     }

  }

  int get_compressed_number_of_elements(){

    Eigen::MatrixXd tmp;
    int retval = 0;

    //how many leaves are there
    int num_leaves = meas_bct_.leaf_pointers_->size();

    for (int l = 0; l < num_leaves; ++l){

      if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){
        //LOW RANK
        tmp = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_L();
        retval += tmp.rows() * tmp.cols();
        tmp = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_R();
        retval += tmp.rows() * tmp.cols();

      }
      else{
        //DENSE
        tmp = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_F();
        retval += tmp.rows() * tmp.cols();
      }

    }

    return retval;
  }

  Eigen::MatrixXd get_sparsity_pattern(){

    //number of patches
    int num_patches = ansatz_space_->get_number_of_patches();
    //refinement level
    int refine_lvl = ansatz_space_->get_refinement_level();

    //get polynomail degree
    int I = ansatz_space_->get_polynomial_degree();
    int I2 = (I + 1) * (I + 1);

    //number of degrees of freedom of large matrix
    int num_DoFs_large = I2*num_patches*std::pow(4,refine_lvl);

    //pointer to element cluster
    const ElementTreeNode *this_el_cluster;

    //for segment in rhs
    int k0, k1;

    //for measurement indices
    std::vector<int> meas_indx;

    //return matrix
    Eigen::MatrixXd ret_val;
    ret_val.resize(num_meas_,num_DoFs_large);

    //how many leaves are there
    int num_leaves = meas_bct_.leaf_pointers_->size();

    for (int l = 0; l < num_leaves; ++l){

      this_el_cluster = (*meas_bct_.leaf_pointers_)[l]->get_element_cluster();


      k0 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_begin(*this_el_cluster));
      k1 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_end(*this_el_cluster));

      meas_indx = (*meas_bct_.leaf_pointers_)[l]->meas_cluster_->indices_;

      if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){

        for(int m = 0; m < meas_indx.size();++m){
          for(int k = 0; k < k1-k0 ; ++k){
            ret_val(meas_indx[m],k0+k) = (*meas_bct_.leaf_pointers_)[l]->lvl_;
          }
        }


      }
      else{
        for(int m = 0; m < meas_indx.size();++m){
          for(int k = 0; k < k1-k0 ; ++k){
            ret_val(meas_indx[m],k0+k) = -1;
          }
        }
      }


    }

    return ret_val;

  }



   void mat_vec_prod(Eigen::VectorXd *res,const Eigen::VectorXd &rhs){

    res->resize(num_meas_);
    res->setZero();

    //number of leaves
    int num_leaves = meas_bct_.leaf_pointers_->size();

    //polynomial degree plus one squared
    int I2 = meas_bct_.polynomial_degree_plus_one_squared_;


    // go discontinuous in rhs
    Eigen::VectorXd long_rhs = (ansatz_space_->get_transformation_matrix() * rhs).eval();


    Eigen::VectorXd sum;
    sum.resize(num_meas_);
    sum.setZero();

    //res->resize(num_meas_);
    //res->setZero();


    //#pragma omp parallel
    {


    //#pragma omp for nowait
    #pragma omp parallel for
    for (int l = 0; l < num_leaves; ++l){

      //temporal storage
      Eigen::VectorXd tmp_res(num_meas_);

      //for segment in rhs
      int k0, k1, tmp;

      k0 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_begin(*(*meas_bct_.leaf_pointers_)[l]->get_element_cluster()));
      k1 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_end(*(*meas_bct_.leaf_pointers_)[l]->get_element_cluster()));

      std::vector<int> meas_indx = (*meas_bct_.leaf_pointers_)[l]->meas_cluster_->indices_;

      if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){

        tmp = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_R().rows();

        //LOW RANK
        tmp_res.block(0,0,tmp,1) = ((*meas_bct_.leaf_pointers_)[l]->leaf_->get_R() *  long_rhs.segment(k0,k1-k0)).eval();
        tmp_res.segment(0,meas_indx.size()) = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_L() * tmp_res.segment(0,tmp);

      }
      else{

        tmp_res.segment(0,meas_indx.size()) = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_F() *  long_rhs.segment(k0,k1-k0);

      }

      #pragma omp critical
      for(int m = 0; m < meas_indx.size();++m){
        sum(meas_indx[m]) += tmp_res(m);
      }

    }

    }
    (*res) = sum;
   }
   //!!!!!!!!!!! CHECK THIS CODE. NOT WORKING
   void mat_vec_transpose(Eigen::VectorXd *res,const Eigen::VectorXd &rhs){


    int num_DoFs = ansatz_space_->get_number_of_dofs();

    const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

    //res->resize(num_DoFs);
    //res->setZero();


    //number of leaves
    int num_leaves = meas_bct_.leaf_pointers_->size();

    //polynomial degree plus one squared
    int I2 = meas_bct_.polynomial_degree_plus_one_squared_;


    // go discontinuous in rhs
    Eigen::VectorXd long_lhs;
    long_lhs.resize(ansatz_space_->get_transformation_matrix().rows());
    long_lhs.setZero();

    //#pragma omp parallel
    {

    //#pragma omp for nowait
    #pragma omp parallel for
    for (int l = 0; l < num_leaves; ++l){


      //for segment in lhs
      int k0, k1, tmp;

      k0 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_begin(*(*meas_bct_.leaf_pointers_)[l]->get_element_cluster()));
      k1 = I2*std::distance(el_memory_->cpbegin(),el_memory_->cluster_end(*(*meas_bct_.leaf_pointers_)[l]->get_element_cluster()));

      std::vector<int> meas_indx = (*meas_bct_.leaf_pointers_)[l]->meas_cluster_->indices_;


      //temporal storage
      Eigen::VectorXd tmp_res(k1-k0);

      if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){

        //LOW RANK

        //tmp = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_L().cols();

        //tmp_res.segment(0,tmp) = prod_L_transpose((*meas_bct_.leaf_pointers_)[l]->leaf_->get_L(), rhs, meas_indx);//.segment(0,tmp)

        tmp_res = (*meas_bct_.leaf_pointers_)[l]->leaf_->get_R().transpose()
                                              * mat_prod_transpose((*meas_bct_.leaf_pointers_)[l]->leaf_->get_L(), rhs, meas_indx);



      }
      else{

        tmp_res = mat_prod_transpose((*meas_bct_.leaf_pointers_)[l]->leaf_->get_F(), rhs, meas_indx);


      }

      #pragma omp critical
      long_lhs.segment(k0,k1-k0) += tmp_res;


    }

    }

    (*res) += ansatz_space_->get_transformation_matrix().transpose() * long_lhs;
   }

  int get_number_of_dense_leafs(){

       int retval = 0;
       //how many leaves are there
       int num_leaves = meas_bct_.leaf_pointers_->size();
       for (int l = 0; l < num_leaves; ++l){

         if((*meas_bct_.leaf_pointers_)[l]->cc_ != 1){
           ++retval;
         }
       }
       return retval;
  }

  int get_number_of_low_rank_leafs(){

       int retval = 0;
       //how many leaves are there
       int num_leaves = meas_bct_.leaf_pointers_->size();
       for (int l = 0; l < num_leaves; ++l){

         if((*meas_bct_.leaf_pointers_)[l]->cc_ == 1){
           ++retval;
         }
       }
       return retval;
  }

   template <class T>
   void assemble_low_rank_matrix_block(const T& super_space,
                                    const ElementTreeNode &el_cluster,
                                    const MeasurementTreeNode &meas_cluster,
                                    Bembel::Quadrature<2>& Q,
                                    Eigen::MatrixXd* L, Eigen::MatrixXd* R ){

      //polynomial degree plus one squared
      int I2 = super_space.get_polynomial_degree_plus_one_squared();

      //number of measurements in this cluster
      int M = meas_cluster.indices_.size();

      //number of basis functions in this cluster
      int K = std::distance(el_memory_->cluster_begin(el_cluster),el_memory_->cluster_end(el_cluster)) * I2;



      //minimum dimension
      int min_dim = M;
      if (K > M) min_dim = K;

      //Low rank aproximations
      Eigen::MatrixXd U,V;
      Eigen::MatrixXd A_col,A_row;

      U.resize(M,K);
      V.resize(M,K);
      U.setZero();
      V.setZero();

      Eigen::VectorXd r_row, r_col;
      r_row.resize(K);
      r_col.resize(M);

      Eigen::Vector3d this_pos;

      //std::cout << "Data Management" << std::endl;

      std::vector<int> I, J;
      int k = 0;
      //$i^\asterix$
      int is = 0;
      int js;
      double delta_k = 0;
      double Af_sq = 0.;


      //------------------------------------------------------
      //ACA assembly
      //------------------------------------------------------

      while(true){

        this_pos = pos_->row(meas_cluster.indices_[is]);


        generate_row_component(super_space,el_cluster,this_pos,Q, &A_row);

        residual_row(&r_row,A_row,U,V,is,k);


        //$ j^\asterix = \text{argmax}_j |R_{k-1}(i^\asterix,j)| $
         //js = find_max(r_row.array().abs());
         js = find_max_constrained(r_row.array().abs(),J);
         if (js < 0) break;
         J.push_back(js);


         //$ \delta_k $
         delta_k = r_row(js);

         r_row /= delta_k;

         if (fabs(delta_k) < 1e-12){
           if (I.size() == min_dim - 1){
             break;
           }
         }
         else{
           generate_col_component(super_space,js,el_cluster,meas_cluster,Q, &A_col);

           residual_col(&r_col,A_col,U,V,js,k);



           U.col(k) = r_col;
           V.row(k) = r_row;

           ++k;

           I.push_back(is);


           is = find_max_constrained(r_col.array().abs(),I);
           if (is < 0) break;


         }

         //std::cout << "k = " << k << std::endl;

         if(check_stop(U, V, &Af_sq,k) || (k >= min_dim) ){
           //std::cout << "break" << std::endl;
           break;
         }

       }

       (*L) = U.block(0,0,M,k);
       (*R) = V.block(0,0,k,K);



   }

   MeasurementBockClusterTree get_block_cluster_tree(){

     return meas_bct_;
   }


 private:
   double eta_;                             //admissibility criterion
   double eps_;                             //stopping criterion for ACA
   double xi_;                              //threshold to adapt quadrature degree
   int max_tree_lvl_;
   int min_cluster_level_;
   int deg_;
   int num_meas_;
   MeasurementBockClusterTree meas_bct_;
   std::shared_ptr<Bembel::ElementTreeMemory> el_memory_;
   Eigen::Vector3d tmp_derivative_;
   SurfacePoint qp_;
   Eigen::MatrixXd *pos_;
   DiscreteSensor<Device,LinOp> sensor_;
   AnsatzSpace<LinOp> *ansatz_space_;

   enum cross_type{ ROW, COLUMN };

   Eigen::VectorXd mat_prod_transpose(const Eigen::MatrixXd &L, const Eigen::VectorXd &rhs, const std::vector<int> indices){

     Eigen::VectorXd ret_val(L.cols());
     ret_val.setZero();

     for(int k = 0 ; k < L.cols(); ++k){
       for(int m = 0 ; m < L.rows(); ++m){

         ret_val(k) += L(m,k) * rhs(indices[m]);

       }
     }

     return ret_val;
   }


   //template <typename Derived, typename LinOp, class T>
   template <class T>
   void assemble_dense_matrix_block(const T& super_space,
                                    const ElementTreeNode &el_cluster,
                                    const MeasurementTreeNode &meas_cluster,
                                    Bembel::Quadrature<2>& Q,
                                    Eigen::MatrixXd* F){

        //polynomial degree plus one squared
        int I2 = super_space.get_polynomial_degree_plus_one_squared();

        //number of measurements in this cluster
        int M = meas_cluster.indices_.size();

        //number of basis functions in this cluster
        int K = std::distance(el_memory_->cluster_begin(el_cluster),el_memory_->cluster_end(el_cluster)) * I2;

        //make space for matrices
        F->resize(M,K);

        //init matrices
        F->setZero();

        //temporal storage for matrix blocks
        Eigen::MatrixXd tmp_D;
        tmp_D.resize(1,I2);


        //surface point for quadratures
        SurfacePoint qp;

        //running index for basis functions
        int k;

        //iterate over measurements
       for (int m = 0; m < meas_cluster.indices_.size(); ++m) {

        k = 0;

        //iterate over elements
        for (auto element = el_memory_->cluster_begin(el_cluster);
             element != el_memory_->cluster_end(el_cluster); ++element) {

                 tmp_D.setZero();

                   //Gaussian integtation
                   for (auto j = 0; j < Q.w_.size(); ++j) {
                      super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                      sensor_.get_sensor().evaluateIntegrand(super_space, *element,pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D);
                    }

                    F->block(m,I2 * k,1,I2) = tmp_D.row(0);

                  ++k;
                  }

                }
          }

   void residual_col(Eigen::VectorXd *res, const Eigen::MatrixXd &A_col,
                                const Eigen::MatrixXd &U,
                                const Eigen::MatrixXd &V, int j, int k){


      (*res) = A_col.col(0);
      if(k>0){

        (*res) -= U.block(0,0,U.rows(),k) * V.block(0,j,k,1);
      }


   }

   void residual_row(Eigen::VectorXd *res, const Eigen::MatrixXd &A_row,
                                const Eigen::MatrixXd &U,
                                const Eigen::MatrixXd &V, int i, int k){


      (*res) = A_row.row(0);
      if(k>0){
        (*res) -= (U.block(i,0,1,k) * V.block(0,0,k,V.cols())).transpose();
    }
   }



   bool check_stop(const Eigen::MatrixXd &U, const Eigen::MatrixXd &V, double *Af_sq_k, int k){

      int M = U.rows();
      int K = V.cols();
      double Af_sq = 0.;
      double uk_sq = U.col(k-1).squaredNorm();
      double vk_sq = V.row(k-1).squaredNorm();

      if (k>1){
        for (int j = 0; j < k-1 ; ++j){
          Af_sq += fabs((U.col(k-1).transpose() * U.col(j)).eval() * (V.row(j) *  V.row(k-1).transpose()).eval());
        }
      }

      Af_sq *= 2.;

      Af_sq +=  fabs(uk_sq * vk_sq) + (*Af_sq_k);

      (*Af_sq_k) = Af_sq;

      if (fabs(uk_sq*vk_sq) < eps_*eps_* fabs(Af_sq)){
        return true;
      }
      else{

        return false;
      }

   }

   template <class T>
   void generate_row(const T& super_space, const ElementTreeNode &el_cluster,
                      const Eigen::Vector3d &position, Bembel::Quadrature<2>& Q,
                      Eigen::MatrixXd* F){

     //polynomial degree plus one squared
     int I2 = super_space.get_polynomial_degree_plus_one_squared();


     //number of basis functions in this cluster
     int K = std::distance(el_memory_->cluster_begin(el_cluster),el_memory_->cluster_end(el_cluster)) * I2;

     //make space for matrices
     F->resize(1,K);


     //init matrices
     F->setZero();

     //temporal storage for matrix blocks
     Eigen::VectorXd tmp_D;


     tmp_D.resize(I2);


     //surface point for quadratures
     SurfacePoint qp;

     //running index for basis functions
     int k;

     k = 0;

     //iterate over elements
     for (auto element = el_memory_->cluster_begin(el_cluster);
          element != el_memory_->cluster_end(el_cluster); ++element) {

              tmp_D.setZero();

                //Gaussian integtation
                for (auto j = 0; j < Q.w_.size(); ++j) {
                   super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                   sensor_.get_sensor().evaluateIntegrand(super_space, *element,position ,qp, &tmp_D);
                 }

                 F->block(0,I2 * k,1,I2) = tmp_D;

               ++k;
     }


  }

  template <class T>
  void generate_row_component(const T& super_space, const ElementTreeNode &el_cluster,
                     const Eigen::Vector3d &position, Bembel::Quadrature<2>& Q,
                     Eigen::MatrixXd* F){ //, const int component

    //polynomial degree plus one squared
    int I2 = super_space.get_polynomial_degree_plus_one_squared();



    //number of basis functions in this cluster
    int K = std::distance(el_memory_->cluster_begin(el_cluster),el_memory_->cluster_end(el_cluster)) * I2;

    //make space for matrices
    F->resize(1,K);

    //init matrices
    F->setZero();

    //temporal storage for matrix blocks
    Eigen::MatrixXd tmp_D;
    tmp_D.resize(1,I2);

    //surface point for quadratures
    SurfacePoint qp;

    //running index for basis functions
    int k;

    k = 0;

    //iterate over elements
    for (auto element = el_memory_->cluster_begin(el_cluster);
         element != el_memory_->cluster_end(el_cluster); ++element) {

             tmp_D.setZero();

               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.get_sensor().evaluateIntegrand(super_space, *element, position ,qp, &tmp_D);

                }

                F->block(0,I2 * k,1,I2) = tmp_D;

               ++k;
    }
  }

  template <class T>
  void generate_col_component(const T& super_space, const int basis_index,const ElementTreeNode &el_cluster,
                     const MeasurementTreeNode &meas_cluster, Bembel::Quadrature<2>& Q,
                     Eigen::MatrixXd* F){ //this_tf[j]

    //polynomial degree plus one squared
    int I2 = super_space.get_polynomial_degree_plus_one_squared();


    //number of measurements in this cluster
    int M = meas_cluster.indices_.size();

    //make space for matrices
    F->resize(M,1);

    //init matrices
    F->setZero();


    int el_indx    = basis_index / I2;
    int basis_indx = basis_index % I2;

    //temporal storage for matrix blocks
    double tmp_D;

    //surface point for quadratures
    SurfacePoint qp;

    auto element = get_cluster_element(el_cluster,el_indx);


    //iterate over measurements
   for (int m = 0; m < meas_cluster.indices_.size(); ++m) {

             tmp_D = 0.;


               //Gaussian integtation
               for (auto j = 0; j < Q.w_.size(); ++j) {
                  super_space.map2surface(*element, Q.xi_.col(j), element->get_h() * element->get_h() * Q.w_(j), &qp);

                  sensor_.get_sensor().evaluateIntegrand_single_basis(super_space, *element,pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D,basis_indx);

                };

                (*F)(m,0) = tmp_D;

    }
  }



    int find_zero(Eigen::VectorXd vec){

      int i = 0;
      for( ; i < vec.size(); ++i) if (vec[i] == 0) break;

      return i;

    }

    int find_max(Eigen::VectorXd vec){

      double max_val = 0;
      int max_coeff = 0;
      int i = 0;
      for( ; i < vec.size(); ++i) {

        if (vec[i] > max_val){
          max_val = vec[i];
          max_coeff = i;
        }

      }

      return max_coeff;

    }

    int find_max_constrained(Eigen::VectorXd vec,std::vector<int> constr){

      double max_val = 0;
      int max_coeff = -1;
      bool in_constr = false;
      int i = 0;

      for( ; i < vec.size(); ++i) {

        in_constr = false;
        if (vec[i] > max_val){
          for (int ic = 0; ic < constr.size() ; ++ic){

            if (constr[ic] == i){

              in_constr = true;
            }
          }

          if ( in_constr == false){
            max_val = vec[i];
            max_coeff = i;
          }

        }

      }

      return max_coeff;
    }

    bool linear_dependent(const Eigen::MatrixXd &ab, const Eigen::MatrixXd &ruv){

      double sum_abs_ab = 0;
      double sum_abs_ruv = 0;

      bool is_lin_dep = false;

      for(int i = 0; i < 3 ; ++i){



        sum_abs_ab = ab.row(i).array().abs().sum();
        sum_abs_ruv = ruv.row(i).array().abs().sum();

        //std::cout << i << " ab = " << sum_abs_ab << std::endl;
        //std::cout << i << " ruv = " << sum_abs_ruv << std::endl;


        if ((sum_abs_ab < 1e-12) || (sum_abs_ruv < 1e-12)){

          is_lin_dep = true;
          break;
        }

      }
      return is_lin_dep;
    }

    std::vector<ElementTreeNode>::const_iterator get_cluster_element(const ElementTreeNode &el_cluster,const int el_indx){

      auto element = el_memory_->cluster_begin(el_cluster);

      int e = 0;

      //find the element
      for (;element != el_memory_->cluster_end(el_cluster); ++element) {
             if(e == el_indx) break;
              ++e;
      }
      return element;
    }



 };


}  // namespace Bembel
#endif
