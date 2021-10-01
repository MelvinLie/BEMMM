// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_ACAMATRIX_BF_H_
#define BEMBEL_ACAMATRIX_BF_H_

#include <Eigen/Dense>
#include <chrono>

namespace Bembel {



 //template <typename Scalar>
 template<typename Device, typename LinOp>
 class ACAMatrixBF {

 public:
   //constructor
   ACAMatrixBF(){
     eta_ = 1.6;
     max_tree_lvl_ = 3;
     min_cluster_level_ = 1;
     eps_ = 1e-5;

   }
   ACAMatrixBF(const DiscreteSensor<Device,LinOp> &sensor, MeasurementData *meas_data, AnsatzSpace<LinOp> *ansatz_space,double eta = 1.6, int max_tree_lvl = 3, int min_cluster_level = 1, double eps = 1e-5){
     eta_ = eta;
     max_tree_lvl_ = max_tree_lvl;
     min_cluster_level_ = min_cluster_level;
     eps_  = eps;

     init_ACAMatrix(sensor,meas_data, ansatz_space);

   }


   void set_parameters(double eta, int max_tree_level, int min_cluster_level, double eps = 1e-5){
     eta_ = eta;
     max_tree_lvl_ = max_tree_level;
     min_cluster_level_ = min_cluster_level;
     eps_ = eps;


   }

   double get_eta(){

     return eta_;
   }

   //template<typename Derived, typename LinOp>
   //TO DO MAKE DEVICE ARGUMENT POINTER
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


     //Get measurement data pointer
     pos_ = meas_data->get_position_ptr();

     //std::cout << meas_data->get_positions() << std::endl;

     //total number of measurements
     num_meas_ = meas_data->get_number_of_measurements();

     //total number of dofs
     num_dofs_ = ansatz_space->get_number_of_dofs();


     //setup block cluster tree
     bct_ = BaseFcnMeasBockClusterTree<LinOp>(meas_data,*ansatz_space_,max_tree_lvl_);

     //initialize block cluster tree
     bct_.init_BlockClusterTree();

     //bct_.print_block_cluster_tree();

     //how many leaves are there
     int num_leaves = bct_.leaf_pointers_->size();

     //lookup table for basis splines
     BSplines_ = bct_.get_BaseFcnTree_ptr()->make_BSpline_table();


     //for(std::vector<BSplineBasisCombination>::iterator it = BSplines_[0].begin();it != BSplines_[0].end()

     //pointer to element cluster memory
     element_it_ = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();


//#pragma omp parallel
     {
     //Gauss quadrature rule
     Bembel::GaussSquare<Bembel::Constants::maximum_quadrature_degree> GS;
     auto Q = GS[deg_];

     #pragma omp parallel for
     for (int l = 0; l < num_leaves; ++l){

       //temporal storage for dense matrices
       Eigen::MatrixXd F;

       //temporal storage for low rank matrix
       Eigen::MatrixXd L, R;

       //pointers to the element and measurement clusters
       const BaseFcnTreeNode *this_dofs;
       const MeasurementTreeNode *this_meas;

       this_dofs = (*bct_.leaf_pointers_)[l]->get_dof_cluster();
       this_meas = (*bct_.leaf_pointers_)[l]->get_measurement_cluster();

       //std::cout << "Block = ("<< this_meas->indices_.size() << " x " << this_dofs->indices_.size() << ")" <<std::endl;

       if((*bct_.leaf_pointers_)[l]->cc_ == 1){
         //LOW RANK
         //std::cout << "Low rank" << std::endl;
         assemble_low_rank_matrix_block(super_space,*this_dofs,*this_meas,Q, BSplines_, &L,&R);

         (*bct_.leaf_pointers_)[l]->leaf_->set_L(L);
         (*bct_.leaf_pointers_)[l]->leaf_->set_R(R);

         (*bct_.leaf_pointers_)[l]->leaf_->set_low_rank_flag(true);


       }
       else{
         //DENSE
         //std::cout << "Dense" << std::endl;
         //assemble dense matrix block

         assemble_dense_matrix_block(super_space,*this_dofs,*this_meas,Q,BSplines_,&F); //, this_el, this_meas, GS, &F

         (*bct_.leaf_pointers_)[l]->leaf_->set_F(F);

         (*bct_.leaf_pointers_)[l]->leaf_->set_low_rank_flag(false);


         }
       }
     }
  }


  void test_yourself(){

    std::cout << "Test Yourself" << std::endl;
    //------------------------------------
    //Mesh
    //------------------------------------
    //Load a unit sphere
    Geometry geometry("../sphere.dat");

    int refinement_level = 2;
    int polynomial_degree = 2;

    // Build ansatz space
    AnsatzSpace<LaplaceHypersingularOperator> ansatz_space(
        geometry, refinement_level, polynomial_degree);

    //get ansatz_space
    ansatz_space_ = &ansatz_space;

    //get super space
    auto super_space = ansatz_space.get_superspace();


    BaseFcnTreeNode dof_cluster;

    std::ifstream file;
    std::string current_line;
    std::string current_element;

    //lookup table for basis splines
    BSplines_.resize(150);

    file.open("../../data/tests/test_aca_bf.csv");

    int dof_ctr;
    int element_indx;
    int loc_indx;
    double weight;

    while(std::getline(file, current_line)){


      if(current_line.substr(0, 3).compare("DOF") == 0){

        dof_ctr = atoi(current_line.substr(3, 6).c_str());
        dof_cluster.indices_.push_back(dof_ctr);
      }
      else{
        std::stringstream current_data(current_line);

        std::getline(current_data,current_element,',');
        element_indx = atoi(current_element.c_str());
        std::getline(current_data,current_element,',');
        loc_indx= atoi(current_element.c_str());
        std::getline(current_data,current_element,',');
        weight = atof(current_element.c_str());

        BSplines_[dof_ctr].append_bernstein_basis(element_indx,loc_indx,weight);
      }


    }

    file.close();

    int num_smpls = 5;
    Eigen::MatrixXd test_pos(num_smpls*num_smpls*num_smpls,3);
    Eigen::VectorXd test_line = Eigen::VectorXd::LinSpaced(num_smpls,0.6,0.65);

    for(int i = 0; i < num_smpls; ++i){
      for(int j = 0; j < num_smpls; ++j){
        for(int k = 0; k < num_smpls; ++k){
          test_pos.row(i*num_smpls*num_smpls+j*num_smpls+k) = Eigen::Vector3d(test_line(i), test_line(j), test_line(k));
        }
      }
    }

    pos_ = &test_pos;

    MeasurementTreeNode meas_cluster;
    for (int i = 0; i < test_pos.rows(); ++i) meas_cluster.append_index(i);

    //quadrature degree
    deg_ =  3;

    //pointer to element cluster memory
    element_it_ = ansatz_space_->get_superspace().get_mesh().get_element_tree().cpbegin();

    //Gauss quadrature rule
    Bembel::GaussSquare<Bembel::Constants::maximum_quadrature_degree> GS;
    auto Q = GS[deg_];

    Eigen::MatrixXd F,L,R;

    std::chrono::time_point<std::chrono::steady_clock> t_start, t_end;
    std::chrono::duration<double> t_dense,t_aca;


    t_start =  std::chrono::steady_clock::now();
    assemble_dense_matrix_block(super_space,dof_cluster,meas_cluster,Q,BSplines_,&F);
    t_end =  std::chrono::steady_clock::now();
    t_dense = t_end - t_start;

    t_start =  std::chrono::steady_clock::now();
    assemble_low_rank_matrix_block(super_space,dof_cluster,meas_cluster,Q,BSplines_, &L,&R);
    t_end =  std::chrono::steady_clock::now();
    t_aca = t_end - t_start;

    int len_dense = F.rows()*F.cols();
    int len_aca = R.rows()*R.cols()+L.rows()*L.cols();


    std::cout << "                    |     dense      |      aca       |"  << std::endl;
    std::cout << "--------------------|----------------|----------------|"  << std::endl;
    std::cout << "Elapsed time [sec.] | " << std::setw(14) << t_dense.count() << " | "<< std::setw(14) << t_aca.count() << " |" << std::endl;
    std::cout << "Number of elements  | " << std::setw(14) << len_dense << " | "<< std::setw(14) << len_aca << " |" << std::endl;

    Eigen::VectorXd test_rhs = Eigen::VectorXd::Random(F.cols());

    Eigen::VectorXd matvec_dense = F*test_rhs;
    Eigen::VectorXd matvec_aca = L*R*test_rhs;

    std::cout << "Compression Ratio = " << ((double) len_aca)/len_dense << std::endl;
    std::cout << "Max Error = " << (matvec_dense-matvec_aca).array().abs().maxCoeff() << std::endl;

    file.open("../../data/tests/dense_mat_aca_test.csv");

    int row_indx = 0;
    int col_indx = 0;

    Eigen::MatrixXd D_in(num_smpls*num_smpls*num_smpls,152);
    Eigen::MatrixXd D_bm(num_smpls*num_smpls*num_smpls,dof_cluster.indices_.size());

    while(std::getline(file, current_line)){
      std::stringstream current_data(current_line);
      col_indx = 0;
      while(std::getline(current_data, current_element,',')){
        //std::cout << row_indx << " x " << col_indx << std::endl;
        D_in(row_indx,col_indx) = atof(current_element.c_str());
        ++col_indx;
      }
      ++row_indx;
    }
    for(int i = 0; i < dof_cluster.indices_.size() ;++i){
      D_bm.col(i) = D_in.col(dof_cluster.indices_[i]);
    }
    file.close();

    std::cout << "Benchmark" << std::endl;
    std::cout << D_bm.block(0,0,10,10) << std::endl;
    std::cout << "Dense" << std::endl;
    std::cout << F.block(0,0,10,10) << std::endl;
    std::cout << "ACA" << std::endl;
    std::cout << (L*R).block(0,0,10,10) << std::endl;

    /*
    SHOULD BE
    BM
    -0.000242935 -0.000131203 -8.04244e-05 -0.000307821 -0.000164054 -0.000131203 -0.000589154 -0.000307821 -0.000242935 -0.000263118
    -0.000238529 -0.000128994 -7.91495e-05 -0.000302698 -0.000161526 -0.000129287 -0.000580286 -0.000303548  -0.00023974 -0.000259731
    -0.000234185  -0.00012681 -7.78869e-05  -0.00029762 -0.000159014 -0.000127383 -0.000571433 -0.000299276 -0.000236542 -0.000256342
    -0.000229902 -0.000124653 -7.66368e-05 -0.000292587 -0.000156519 -0.000125489 -0.000562599 -0.000295006 -0.000233343 -0.000252954
    -0.000225681 -0.000122523 -7.53996e-05 -0.000287602 -0.000154043 -0.000123606 -0.000553787  -0.00029074 -0.000230146 -0.000249568
     -0.00023974 -0.000129287 -7.91495e-05 -0.000303548 -0.000161526 -0.000128994 -0.000580286 -0.000302698 -0.000238529  -0.00025802
    -0.000235415 -0.000127125 -7.79042e-05 -0.000298521 -0.000159051 -0.000127125 -0.000571589 -0.000298521 -0.000235415 -0.000254728
     -0.00023115 -0.000124987 -7.66707e-05 -0.000293537 -0.000156593 -0.000125266 -0.000562904 -0.000294344 -0.000232298 -0.000251434
    -0.000226945 -0.000122875 -7.54493e-05 -0.000288597 -0.000154152 -0.000123418 -0.000554236 -0.000290169  -0.00022918  -0.00024814
      -0.0002228 -0.000120788 -7.42401e-05 -0.000283703 -0.000151729  -0.00012158 -0.000545591 -0.000285997 -0.000226062 -0.000244848
    */
  }

  int get_nonzeros_in_projector(){

    int nonzeros = 0;

    for(int i = 0; i < (*element_index_table_).size(); ++i){
      std::cout << "DoF " << i << ": ";
      for(int j = 0; j < (*element_index_table_)[i].size(); ++j){

        std::cout << (*element_index_table_)[i][j] << " , ";

      }
      std::cout << "\n" << std::endl;

      nonzeros += (*element_index_table_)[i].size();
    }

    return nonzeros;

  }

  Eigen::MatrixXd get_dense(){

    Eigen::MatrixXd ret_val(num_meas_,num_dofs_);
    ret_val.setZero();
    Eigen::MatrixXd tmp;
    std::vector<int> meas_indx,dof_indx;

    //pointer to element cluster
    const BaseFcnTreeNode *this_dof_cluster;

    int num_leaves = bct_.leaf_pointers_->size();



    for (int l = 0; l < num_leaves; ++l){


      this_dof_cluster = (*bct_.leaf_pointers_)[l]->get_dof_cluster();

      dof_indx = (*bct_.leaf_pointers_)[l]->dof_cluster_->indices_;
      meas_indx = (*bct_.leaf_pointers_)[l]->meas_cluster_->indices_;


      tmp.resize(meas_indx.size(),dof_indx.size());

      if((*bct_.leaf_pointers_)[l]->cc_ == 1){

        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_L()*(*bct_.leaf_pointers_)[l]->leaf_->get_R();
      }
      else{

        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_F();

      }

      for(int i = 0; i < meas_indx.size(); ++i){
        for(int j = 0; j < dof_indx.size(); ++j){

          ret_val(meas_indx[i],dof_indx[j]) = tmp(i,j);

        }
      }

    }
    return ret_val;

  }


  int get_compressed_number_of_elements(){

    Eigen::MatrixXd tmp;
    int retval = 0;

    //how many leaves are there
    int num_leaves = bct_.leaf_pointers_->size();

    for (int l = 0; l < num_leaves; ++l){

      if((*bct_.leaf_pointers_)[l]->cc_ == 1){
        //LOW RANK
        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_L();
        retval += tmp.rows() * tmp.cols();
        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_R();
        retval += tmp.rows() * tmp.cols();

      }
      else{
        //DENSE
        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_F();
        retval += tmp.rows() * tmp.cols();
      }

    }

    return retval;
  }

  int get_number_of_dense_leafs(){

    Eigen::MatrixXd tmp;
    int retval = 0;
    //how many leaves are there
    int num_leaves = bct_.leaf_pointers_->size();
    for (int l = 0; l < num_leaves; ++l){

      if((*bct_.leaf_pointers_)[l]->cc_ != 1){
        ++retval;
      }
    }
    return retval;
  }

  int get_number_of_low_rank_leafs(){

    Eigen::MatrixXd tmp;
    int retval = 0;
    //how many leaves are there
    int num_leaves = bct_.leaf_pointers_->size();
    for (int l = 0; l < num_leaves; ++l){

      if((*bct_.leaf_pointers_)[l]->cc_ == 1){
        ++retval;
      }
    }
    return retval;
  }


   void mat_vec_prod(Eigen::VectorXd *res,const Eigen::VectorXd &rhs){

    res->resize(num_meas_);
    res->setZero();

    //number of leaves
    int num_leaves = bct_.leaf_pointers_->size();

    //polynomial degree plus one squared
    int I2 = bct_.polynomial_degree_plus_one_squared_;

    Eigen::VectorXd sum;
    sum.resize(num_meas_);
    sum.setZero();


    #pragma omp parallel
    {

    //temporal storage
    Eigen::VectorXd tmp_rhs(num_dofs_);
    tmp_rhs.setZero();

    Eigen::VectorXd tmp_res(num_meas_);
    tmp_res.setZero();

    //for measurement indices
    std::vector<int> meas_indx;

    //dof indices
    std::vector<int> dof_indx;

    //pointer to element cluster
    const BaseFcnTreeNode *this_dof_cluster;

    //for segment in rhs
    int k0, k1, tmp;

    #pragma omp for nowait
    for (int l = 0; l < num_leaves; ++l){

      this_dof_cluster = (*bct_.leaf_pointers_)[l]->get_dof_cluster();

      dof_indx = (*bct_.leaf_pointers_)[l]->dof_cluster_->indices_;
      meas_indx = (*bct_.leaf_pointers_)[l]->meas_cluster_->indices_;

      for(int i = 0; i < dof_indx.size(); ++i){
         tmp_rhs(i) = rhs(dof_indx[i]);
       }

      if((*bct_.leaf_pointers_)[l]->cc_ == 1){

        tmp = (*bct_.leaf_pointers_)[l]->leaf_->get_R().rows();

        //LOW RANK

        tmp_res.block(0,0,tmp,1) = ((*bct_.leaf_pointers_)[l]->leaf_->get_R() *  tmp_rhs.segment(0,dof_indx.size())).eval();

        tmp_res.segment(0,meas_indx.size()) = (*bct_.leaf_pointers_)[l]->leaf_->get_L() * tmp_res.segment(0,tmp);

        //This takes very long!
        //tmp_res.segment(0,meas_indx.size()) = (*meas_bct_.leaf_pointers_)[l]->leafX_->get_L() * ((*meas_bct_.leaf_pointers_)[l]->leafX_->get_R() *  long_rhs.segment(k0,k1-k0)).eval();


      }
      else{

        tmp_res.segment(0,meas_indx.size()) = (*bct_.leaf_pointers_)[l]->leaf_->get_F() *  tmp_rhs.segment(0,dof_indx.size());

      }

      #pragma omp critical
      for(int m = 0; m < meas_indx.size();++m){
        sum(meas_indx[m]) += tmp_res(m);
      }

    }

    }
    (*res) = sum;
   }



   template <class T>
   void assemble_low_rank_matrix_block(const T& super_space,
                                    const BaseFcnTreeNode &dof_cluster,
                                    const MeasurementTreeNode &meas_cluster,
                                    Bembel::Quadrature<2>& Q,
                                    std::vector<BSplineBasisCombination> BSplines_,
                                    Eigen::MatrixXd* L, Eigen::MatrixXd* R ){

      //polynomial degree plus one squared
      int I2 = super_space.get_polynomial_degree_plus_one_squared();

      //number of measurements in this cluster
      int M = meas_cluster.indices_.size();

      //number of basis functions in this cluster
      int K = dof_cluster.indices_.size();

      //std::cout << "Block Size = ( " << M <<" x " << K << " )" << std::endl;

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


        generate_row_component(super_space,dof_cluster,this_pos,Q,BSplines_, &A_row);

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

           generate_col_component(super_space,js,dof_cluster,meas_cluster,Q,BSplines_, &A_col);

           residual_col(&r_col,A_col,U,V,js,k);



           U.col(k) = r_col;
           V.row(k) = r_row;

           ++k;

           I.push_back(is);


           is = find_max_constrained(r_col.array().abs(),I);
           if (is < 0) break;


         }



         if(check_stop(U, V, &Af_sq,k) || (k >= min_dim) ){
           //std::cout << "break" << std::endl;
           break;
         }

       }

       (*L) = U.block(0,0,M,k);
       (*R) = V.block(0,0,k,K);



   }


 private:
   double eta_;                             //admissibility criterion
   double eps_ = 1e-5;                      //stopping criterion for ACA
   double xi_;                              //threshold to adapt quadrature degree
   int max_tree_lvl_;
   int min_cluster_level_;
   int deg_;
   int num_meas_;
   int num_dofs_;
   BaseFcnMeasBockClusterTree<LinOp> bct_;
   std::shared_ptr<Bembel::ElementTreeMemory> el_memory_;
   Eigen::Vector3d tmp_derivative_;
   SurfacePoint qp_;
   Eigen::MatrixXd *pos_;
   DiscreteSensor<Device,LinOp> sensor_;
   AnsatzSpace<LinOp> *ansatz_space_;
   std::vector<ElementTreeNode>::const_iterator element_it_;

   std::vector<BSplineBasisCombination> BSplines_;

   enum cross_type{ ROW, COLUMN };

   std::vector<std::vector<int>> *element_index_table_;
   std::vector<std::vector<int>> *local_index_table_;
   std::vector<std::vector<double>> *weight_table_;




   //template <typename Derived, typename LinOp, class T>
   template <class T>
   void assemble_dense_matrix_block(const T& super_space,
                                    const BaseFcnTreeNode &dof_cluster,
                                    const MeasurementTreeNode &meas_cluster,
                                    Bembel::Quadrature<2>& Q,
                                    std::vector<BSplineBasisCombination> BSplines_,
                                    Eigen::MatrixXd* F){


        //polynomial degree plus one squared
        int I2 = super_space.get_polynomial_degree_plus_one_squared();

        //number of measurements in this cluster
        int M = meas_cluster.indices_.size();

        //number of basis functions in this cluster
        int K = dof_cluster.indices_.size();

        //make space for matrices
        F->resize(M,K);

        //init matrices
        F->setZero();

        //temporal storage for matrix blocks
        //Eigen::MatrixXd tmp_D(1,I2);
        double tmp_D;

        int dof_index;
        int element_index;

        //surface point for quadratures
        SurfacePoint qp;

        int num_fcns;

        //iterate over measurements
        for (int m = 0; m < meas_cluster.indices_.size(); ++m) {

          //iterate over dofs
          for(int i = 0; i < K; ++i){


              dof_index = dof_cluster.indices_[i];

              //iterate over discontinuous basis functions
              //for(std::vector<local_bernstein_combination>::iterator spline_it = BSplines_[dof_index].begin(); spline_it != BSplines_[dof_index].end(); ++spline_it){
              for(std::vector<local_bernstein_combination>::iterator spline_it = BSplines_[dof_index].begin(); spline_it != BSplines_[dof_index].end(); ++spline_it){

                num_fcns = (*spline_it).bernstein_basis_ids.size();

                element_index = (*spline_it).element_id;

                //tmp_D.setZero();
                tmp_D = 0.;
                //Gaussian integtation
                for (int k = 0; k < Q.w_.size(); ++k) {

                    //It took me a day to get this, but here is the problem of this class!
                    // The most expansive part of the assembly is map2surface.... Therefore we would like to call it the minimum amount of times!
                    // For this reason, it is superior in every case to iterate over the elements.
                    // Put this whole class in the trash...
                    super_space.map2surface(element_it_[element_index], Q.xi_.col(k), element_it_[element_index].get_h() * element_it_[element_index].get_h() * Q.w_(k), &qp);

                    sensor_.get_sensor().evaluateIntegrand_BSpline(super_space, element_it_[element_index],
                            pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D,&(*spline_it).bernstein_weights[0],&(*spline_it).bernstein_basis_ids[0],num_fcns);

                  }


                  (*F)(m,i) += tmp_D;


              }
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
  void generate_row_component(const T& super_space, const BaseFcnTreeNode &dof_cluster,
                     const Eigen::Vector3d &position, Bembel::Quadrature<2>& Q,
                     std::vector<BSplineBasisCombination> BSplines_,
                     Eigen::MatrixXd* F){ //, const int component


    //polynomial degree plus one squared
    int I2 = super_space.get_polynomial_degree_plus_one_squared();

    //number of basis functions in this cluster
    int K = dof_cluster.indices_.size();

    //make space for matrices
    F->resize(1,K);

    //init matrices
    F->setZero();

    //temporal storage for matrix blocks
    //Eigen::MatrixXd tmp_D(1,I2);
    double tmp_D = 0.;

    //surface point for quadratures
    SurfacePoint qp;

    int dof_index;
    int element_index;

    int num_fcns;

    //iterate over dofs
    for(int i = 0; i < K; ++i){

        dof_index = dof_cluster.indices_[i];

        //iterate over discontinuous basis functions
        for(std::vector<local_bernstein_combination>::iterator spline_it = BSplines_[dof_index].begin(); spline_it != BSplines_[dof_index].end(); ++spline_it){


          num_fcns = (*spline_it).bernstein_basis_ids.size();

          element_index = (*spline_it).element_id;

          //tmp_D.setZero();
          tmp_D =  0.;
          //Gaussian integtation
          for (int k = 0; k < Q.w_.size(); ++k) {

            super_space.map2surface(element_it_[element_index], Q.xi_.col(k), element_it_[element_index].get_h() * element_it_[element_index].get_h() * Q.w_(k), &qp);



            sensor_.get_sensor().evaluateIntegrand_BSpline(super_space, element_it_[element_index],
                    position ,qp, &tmp_D,&(*spline_it).bernstein_weights[0],&(*spline_it).bernstein_basis_ids[0],num_fcns);

          }
          (*F)(0,i) += tmp_D;
      }

    }
    /*

    //polynomial degree plus one squared
    int I2 = super_space.get_polynomial_degree_plus_one_squared();

    //number of basis functions in this cluster
    int K = dof_cluster.indices_.size();

    //make space for matrices
    F->resize(1,K);

    //init matrices
    F->setZero();

    //temporal storage for matrix blocks
    double tmp_D;

    //surface point for quadratures
    SurfacePoint qp;

    int dof_index;
    int element_index;
    int local_index;
    //int global_index;

    //iterate over dofs
    for(int i = 0; i < K; ++i){

      dof_index = dof_cluster.indices_[i];

      //iterate over discintinuous basis functions
      for(int j = 0; j < (*element_index_table_)[dof_index].size(); ++j){


        element_index = (*element_index_table_)[dof_index][j];
        local_index = (*local_index_table_)[dof_index][j];


        tmp_D = 0.;

        //Gaussian integtation
        for (int k = 0; k < Q.w_.size(); ++k) {


           super_space.map2surface(element_it_[element_index], Q.xi_.col(k), element_it_[element_index].get_h() * element_it_[element_index].get_h() * Q.w_(k), &qp);
           sensor_.get_sensor().evaluateIntegrand_single_basis(super_space, element_it_[element_index], position ,qp, &tmp_D,local_index);

         }

         (*F)(0,i) += (*weight_table_)[dof_index][j]*tmp_D;

      }

    }
    */

  }


  template <class T>
  void generate_col_component(const T& super_space, const int basis_index,const BaseFcnTreeNode &BaseFcnTreeNode,
                     const MeasurementTreeNode &meas_cluster, Bembel::Quadrature<2>& Q,
                     std::vector<BSplineBasisCombination> BSplines_,
                     Eigen::MatrixXd* F){ //this_tf[j]

    //polynomial degree plus one squared
    int I2 = super_space.get_polynomial_degree_plus_one_squared();


    //number of measurements in this cluster
    int M = meas_cluster.indices_.size();

    //make space for matrices
    F->resize(M,1);

    //init matrices
    F->setZero();

    int dof_index = BaseFcnTreeNode.indices_[basis_index];

    int element_index;
    int local_index;
    //int global_index;

    //temporal storage for matrix blocks
    //Eigen::MatrixXd tmp_D(1,I2);
    double tmp_D;

    //surface point for quadratures
    SurfacePoint qp;

    int num_fcns;

    dof_index = BaseFcnTreeNode.indices_[basis_index];

    //iterate over discontinuous basis functions
    for(std::vector<local_bernstein_combination>::iterator spline_it = BSplines_[dof_index].begin(); spline_it != BSplines_[dof_index].end(); ++spline_it){

      element_index = (*spline_it).element_id;

      num_fcns = (*spline_it).bernstein_basis_ids.size();

      //iterate over measurements
      for (int m = 0; m < meas_cluster.indices_.size(); ++m) {

          //tmp_D.setZero();
          tmp_D = 0.;

          //Gaussian integtation
          for (int k = 0; k < Q.w_.size(); ++k) {

            super_space.map2surface(element_it_[element_index], Q.xi_.col(k), element_it_[element_index].get_h() * element_it_[element_index].get_h() * Q.w_(k), &qp);



            sensor_.get_sensor().evaluateIntegrand_BSpline(super_space, element_it_[element_index],
                    pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D,&(*spline_it).bernstein_weights[0],&(*spline_it).bernstein_basis_ids[0],num_fcns);

            //It would be more efficient to only sum up the basis functions we need here!
            //sensor_.get_sensor().evaluateIntegrand(super_space, element_it_[element_index], pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D);

          }
          (*F)(m,0) += tmp_D;
          //for(int j = 0; j < (*spline_it).bernstein_basis_ids.size(); ++j){
          //  (*F)(m,0) += (*spline_it).bernstein_weights[j]*tmp_D(0,(*spline_it).bernstein_basis_ids[j]);
          //}
      }
    }
    /*
    //iterate over measurements
   for (int m = 0; m < meas_cluster.indices_.size(); ++m) {



             //iterate over discintinuous basis functions
             for(int j = 0; j < (*element_index_table_)[dof_index].size(); ++j){

               element_index = (*element_index_table_)[dof_index][j];
               local_index = (*local_index_table_)[dof_index][j];

               //global_index = element_index+local_index;
               tmp_D = 0.;

               //Gaussian integtation
               for (auto k = 0; k < Q.w_.size(); ++k) {

                 super_space.map2surface(element_it_[element_index], Q.xi_.col(k), element_it_[element_index].get_h() * element_it_[element_index].get_h() * Q.w_(k), &qp);

                 sensor_.get_sensor().evaluateIntegrand_single_basis(super_space, element_it_[element_index], pos_->row(meas_cluster.indices_[m]) ,qp, &tmp_D,local_index);

                };

                (*F)(m,0) += (*weight_table_)[dof_index][j]*tmp_D;
            }

    }
    */
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


 };


}  // namespace Bembel
#endif
