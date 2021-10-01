// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_THREEAXISHALLPROBE_H_
#define BEMBEL_THREEAXISHALLPROBE_H_

#include <Eigen/Dense>


namespace Bembel {

/**
 *  \ingroup ThreeAxisHallProbe
 *  \brief Helper class for the treatment of Hall probe measurement data
 */
template <typename Pot, typename LinOp>
class ThreeAxisHallProbe {

 public:
   //////////////////////////////////////////////////////////////////////////////
   //    constructors
   //////////////////////////////////////////////////////////////////////////////
   ThreeAxisHallProbe() {}

   ThreeAxisHallProbe(const AnsatzSpace<LinOp> &ansatz_space) {
     init_ThreeAxisHallProbe(ansatz_space);
   }

   //////////////////////////////////////////////////////////////////////////////
   //    init_ThreeAxisHallProbe
   //////////////////////////////////////////////////////////////////////////////
   void init_ThreeAxisHallProbe(const AnsatzSpace<LinOp> &ansatz_space) {

     ansatz_space_ = ansatz_space;

     deg_ = ansatz_space_.get_polynomial_degree() + 1;

     disc_pot_.init_DiscretePotential(ansatz_space);

     return;
   }



    void load_parameters(const std::string &orientation_file,
                          const std::string &tranfer_fcn_file,const char del = ','){

      std::ifstream file;
      file.open(orientation_file);
      std::string current_line;
      std::string current_element;
      std::stringstream current_data;
      int i = 0;

      std::vector<double> tmp_vector;

      if (!file) {
        std::cerr << "File " << orientation_file << " doesn't exist!" << std::endl;
        exit(1);
      }
      // first two rows are irrelevant
      getline(file, current_line);
      getline(file, current_line);

      //covert this line
      getline(file, current_line);
      current_data << current_line;

      while(std::getline(current_data,current_element,del)){
          n_u_(i) = atof(current_element.c_str());
          ++i;
      }
      current_data.str("");

      i = 0;
      //covert this line
      getline(file, current_line);
      std::stringstream n_v_line(current_line);
      while(std::getline(n_v_line,current_element,del)){
          n_v_(i) = atof(current_element.c_str());
          ++i;
      }

      i = 0;
      //covert this line
      getline(file, current_line);
      std::stringstream n_w_line(current_line);
      while(std::getline(n_w_line,current_element,del)){
          n_w_(i) = atof(current_element.c_str());
          ++i;
      }

      //this line is irrelevant
      getline(file, current_line);

      i = 0;
      //covert this line
      getline(file, current_line);
      std::stringstream d_u_line(current_line);
      while(std::getline(d_u_line,current_element,del)){
          d_u_(i) = atof(current_element.c_str())*1e-3;
          ++i;
      }

      i = 0;
      //covert this line
      getline(file, current_line);
      std::stringstream d_v_line(current_line);
      while(std::getline(d_v_line,current_element,del)){
          d_v_(i) = atof(current_element.c_str())*1e-3;
          ++i;
      }

      i = 0;
      //covert this line
      getline(file, current_line);
      std::stringstream d_w_line(current_line);
      while(std::getline(d_w_line,current_element,del)){
          d_w_(i) = atof(current_element.c_str())*1e-3;
          ++i;
      }

      file.close();

      //sensitivities
      q_u_.clear(); q_v_.clear(); q_w_.clear();

      file.open(tranfer_fcn_file);

      // first two rows are irrelevant
      getline(file, current_line);
      getline(file, current_line);
      getline(file, current_line);
      std::stringstream qu_line(current_line);
      while(std::getline(qu_line,current_element,del)){
          q_u_.push_back(atof(current_element.c_str()));
      }

      getline(file, current_line);
      getline(file, current_line);
      getline(file, current_line);
      std::stringstream qv_line(current_line);
      while(std::getline(qv_line,current_element,del)){
          q_v_.push_back(atof(current_element.c_str()));
      }

      getline(file, current_line);
      getline(file, current_line);
      getline(file, current_line);
      std::stringstream qw_line(current_line);
      while(std::getline(qw_line,current_element,del)){
          q_w_.push_back(atof(current_element.c_str()));
      }

      file.close();


    }

    Eigen::Matrix3d get_orientation_vecs(){

      Eigen::Matrix3d ret_mat;
      ret_mat.row(0) = n_u_;
      ret_mat.row(1) = n_v_;
      ret_mat.row(2) = n_w_;

      return ret_mat;
    }

    Eigen::Matrix3d get_probe_positions(){

      Eigen::Matrix3d ret_mat;
      ret_mat.row(0) = d_u_;
      ret_mat.row(1) = d_v_;
      ret_mat.row(2) = d_w_;

      return ret_mat;
    }

    void get_sensitivities(std::vector<double> *qu,
                            std::vector<double> *qv,
                            std::vector<double> *qw){

      (*qu).resize(q_u_.size());
      (*qv).resize(q_v_.size());
      (*qw).resize(q_w_.size());

      for(int i = 0; i < q_u_.size() ; ++i)   (*qu)[i] = q_u_[i];
      for(int i = 0; i < q_v_.size() ; ++i)   (*qv)[i] = q_v_[i];
      for(int i = 0; i < q_w_.size() ; ++i)   (*qw)[i] = q_w_[i];

    }

    //setup parameters
    void set_parameters(const Eigen::VectorXd &absolute_orientation,
                        const Eigen::MatrixXd &probe_positions,
                         const std::vector<double> &q_u,
                         const std::vector<double> &q_v,
                         const std::vector<double> &q_w) {

        //angles
        alpha_v_ = absolute_orientation(0);
        alpha_w_ = absolute_orientation(1);
        beta_u_  = absolute_orientation(2);
        beta_w_  = absolute_orientation(3);
        gamma_u_ = absolute_orientation(4);
        gamma_v_ = absolute_orientation(5);

        //positions
        d_u_ = probe_positions.col(0);
        d_v_ = probe_positions.col(1);
        d_w_ = probe_positions.col(2);

        //sensitivities
        q_u_.clear(); q_v_.clear(); q_w_.clear();

        for (int i = 0 ; i < q_u.size() ; ++i){
          q_u_.push_back(q_u[i]);
          q_v_.push_back(q_v[i]);
          q_w_.push_back(q_w[i]);
        }


        //orientation vectors
        n_u_(0) =     std::cos(gamma_u_)*std::cos(beta_u_);
        n_u_(1) =     std::sin(gamma_u_)*std::cos(beta_u_);
        n_u_(2) = -1.*std::sin(beta_u_);

        n_v_(0) = -1.*std::sin(gamma_v_)*std::cos(alpha_v_);
        n_v_(1) =     std::cos(gamma_v_)*std::cos(alpha_v_);
        n_v_(2) =     std::sin(alpha_v_);

        n_w_(0) =     std::sin(beta_w_)*std::cos(alpha_w_);
        n_w_(1) = -1.*std::sin(alpha_w_);
        n_w_(2) =     std::cos(beta_w_)*std::cos(alpha_w_);

        }

    //Computes the dense design matrix for measurements at positions pos.
    //Here only the linear part of the transfer function is considered
    //For small problems only
    Eigen::MatrixXd compute_design_matrix(const Eigen::MatrixXd &pos){

      //number of measuremements
      int M = pos.rows();

      //number of DoFs
      int N = ansatz_space_.get_number_of_dofs();
      //note that we will gauge for zero mean arithmetically. So the Matrices
      //will be of size (M x N-1)

      //B field matrices
      Eigen::MatrixXd Bx,By,Bz;

      //make space for B matrices
      Bx.resize(M,N); By.resize(M,N); Bz.resize(M,N);

      //make space for design matrix
      Eigen::MatrixXd D;
      D.resize(3*M,N-1);

      //Helper unit matrix
      Eigen::MatrixXd unit_mat = Eigen::MatrixXd::Ones(M,1);

      //u-axis:
      Eigen::MatrixXd tmp_pos;
      tmp_pos.resize(M,3);
      tmp_pos.col(0) = pos.col(0) + d_u_(0)*unit_mat.col(0);
      tmp_pos.col(1) = pos.col(1) + d_u_(1)*unit_mat.col(0);
      tmp_pos.col(2) = pos.col(2) + d_u_(2)*unit_mat.col(0);

      //Assemble B matrices
      disc_pot_.compute_grad_evaluation_matrix(tmp_pos,Bx,By,Bz);

      //Zero mean gauge
      stability_vec_ = disc_pot_.impose_zero_mean_gauge(Bx,By,Bz);

      D.block(0,0,M,N-1) = q_u_[1]*(n_u_[0]*Bx
                                  + n_u_[1]*By
                                  + n_u_[2]*Bz);

      //v-axis:
      tmp_pos.col(0) = pos.col(0) + d_v_(0)*unit_mat.col(0);
      tmp_pos.col(1) = pos.col(1) + d_v_(1)*unit_mat.col(0);
      tmp_pos.col(2) = pos.col(2) + d_v_(2)*unit_mat.col(0);

      //Assemble B matrices
      disc_pot_.compute_grad_evaluation_matrix(tmp_pos,Bx,By,Bz);

      //Zero mean gauge
      disc_pot_.impose_zero_mean_gauge(Bx,By,Bz);

      D.block(M,0,M,N-1) = q_v_[1]*(n_v_[0]*Bx
                                  + n_v_[1]*By
                                  + n_v_[2]*Bz);

      //w-axis:
      tmp_pos.col(0) = pos.col(0) + d_w_(0)*unit_mat.col(0);
      tmp_pos.col(1) = pos.col(1) + d_w_(1)*unit_mat.col(0);
      tmp_pos.col(2) = pos.col(2) + d_w_(2)*unit_mat.col(0);

      //Assemble B matrices
      disc_pot_.compute_grad_evaluation_matrix(tmp_pos,Bx,By,Bz);

      //Zero mean gauge
      disc_pot_.impose_zero_mean_gauge(Bx,By,Bz);

      D.block(2*M,0,M,N-1) = q_w_[1]*(n_w_[0]*Bx
                                  + n_w_[1]*By
                                  + n_w_[2]*Bz);


      return D;
    }


    //Computes the compresses design matrix for measurements at positions pos.
    //Here only we initialize the matrix to compure n.B. The nonlinearity is
    //applied in the function ...
    void init_aca_matrix(const Eigen::MatrixXd &pos){


    }

    void rotate_probe(double Alpha,double Beta,double Gamma){

         Eigen::Matrix3d Rot_X;

          Rot_X  <<   1,         0      ,       0        ,
                      0,    std::cos(Alpha), -1.*std::sin(Alpha),
                      0,    std::sin(Alpha),     std::cos(Alpha);

          Eigen::Matrix3d Rot_Y;

          Rot_Y <<   std::cos(Beta),    0    ,     std::sin(Beta),
                            0      ,    1    ,        0          ,
                  -1*std::sin(Beta),    0    ,     std::cos(Beta);

          Eigen::Matrix3d Rot_Z;

          Rot_Z <<   std::cos(Gamma),-1.*std::sin(Gamma) ,        0 ,
                     std::sin(Gamma),    std::cos(Gamma) ,        0 ,
                            0       ,           0        ,        1 ;

          n_u_ = Rot_X*Rot_Y*Rot_Z*n_u_;
          n_v_ = Rot_X*Rot_Y*Rot_Z*n_v_;
          n_w_ = Rot_X*Rot_Y*Rot_Z*n_w_;

          d_u_ = Rot_X*Rot_Y*Rot_Z*d_u_;
          d_v_ = Rot_X*Rot_Y*Rot_Z*d_v_;
          d_w_ = Rot_X*Rot_Y*Rot_Z*d_w_;

        }


        Eigen::Vector3d get_offsets(){

          Eigen::Vector3d ret_val;

          ret_val(0) = q_u_[0];
          ret_val(1) = q_v_[0];
          ret_val(2) = q_w_[0];

          return  ret_val;
        }

        Eigen::VectorXd recover_full_solution(Eigen::VectorXd v){

          return disc_pot_.recover_full_solution(v,stability_vec_);
        }

        Eigen::VectorXd get_stabilization_vector(){

          return stability_vec_;

        }

  void print_parameters(){

    std::cout << "n_u = " << n_u_.transpose() << std::endl;
    std::cout << "n_v = " << n_v_.transpose() << std::endl;
    std::cout << "n_w = " << n_w_.transpose() << std::endl;

    std::cout << "d_u_ = " << d_u_.transpose() << std::endl;
    std::cout << "d_v_ = " << d_v_.transpose() << std::endl;
    std::cout << "d_w_ = " << d_w_.transpose() << std::endl;

    std::cout << "q_u_ = " << q_u_[0] ;
    for(int j = 1 ; j < q_u_.size() ; ++j)  std::cout << " + " <<  q_u_[j] << " * B**" << j;
    std::cout << std::endl;

    std::cout << "q_v_ = " << q_v_[0] ;
    for(int j = 1 ; j < q_v_.size() ; ++j)  std::cout << " + " <<  q_v_[j] << " * B**" << j;
    std::cout << std::endl;

    std::cout << "q_w_ = " << q_w_[0] ;
    for(int j = 1 ; j < q_w_.size() ; ++j)  std::cout << " + " <<  q_w_[j] << " * B**" << j;
    std::cout << std::endl;

  }


private:

  //Ansatz space
  AnsatzSpace<LinOp> ansatz_space_;

  //disctrete operator
  DiscretePotential<Pot,LinOp> disc_pot_ ;

  //stability vector
  Eigen::VectorXd stability_vec_;

  //polynomial degree
  int deg_;

  //absolute probe orientation
  double alpha_v_ = 0;
  double alpha_w_ = 0;
  double beta_u_  = 0;
  double beta_w_  = 0;
  double gamma_u_ = 0;
  double gamma_v_ = 0;

  //relative probe orientation
  double phi_uv_ = 0;
  double phi_uw_ = 0;
  double phi_vw_ = 0;

  //absolute orientation of measuring head
  double Alpha_ = 0.;
  double Beta_   = 0.;
  double Gamma_   = 0.;

  //absolute orientation vectors
  Eigen::Vector3d n_u_;
  Eigen::Vector3d n_v_;
  Eigen::Vector3d n_w_;

  //probe positions
  Eigen::Vector3d d_u_;
  Eigen::Vector3d d_v_;
  Eigen::Vector3d d_w_;

  //absolute rotation matrices
  Eigen::Matrix3d Rx_, Ry_, Rz_;

  //transfer functions polynomial coefficients
  std::vector<double> q_u_,q_v_,q_w_;

  //stages positioning accuracy
  double acc_x_ = 3e-6; //[m]
  double acc_y_ = 5e-6; //[m]
  double acc_z_ = 3e-6; //[m]

  void update_absolute_orientation(){

    //these angles are defined in the local uvw coordinate system
    //by definition, this system has:
    //$\vec{n}_u = \left(\cos(\beta_u),0,-\sin(\beta_u)\right)^T $
    //$\vec{n}_v = \left(-\sin(\gamma_v)\cos(\alpha_v),\cos(\alpha_v)\cos(\gamma_v),\sin(\alpha_v)\right)^T $
    //$\vec{n}_w = \left(0,0,1\right)^T $
    //a_w = b_w = g_u = 0
    double a_v, b_u, g_v;

    b_u = phi_uw_ - BEMBEL_PI/2;
    a_v = BEMBEL_PI/2 - phi_vw_;

    g_v = std::asin(-1*(std::cos(phi_uv_)+std::sin(b_u)*std::sin(a_v))/std::cos(b_u)/std::cos(a_v));

    Rx_ <<     1     ,         0       ,      0          ,
               0     , std::cos(Alpha_),-std::sin(Alpha_),
               0     , std::sin(Alpha_), std::cos(Alpha_);

    Ry_ << std::cos(Beta_) ,         0       , std::sin(Beta_) ,
                 0         ,         1       ,      0          ,
          -std::sin(Beta_) ,         0       , std::cos(Beta_) ;

    Rz_ << std::cos(Gamma_) , -std::sin(Gamma_) ,   0 ,
           std::sin(Gamma_) ,  std::cos(Gamma_) ,   0 ,
                   0        ,         0         ,   1 ;

    // orientation vectors in (u,v,w) coordinates
    Eigen::Vector3d nu_p,nv_p,nw_p;
    nu_p.setZero();
    nv_p.setZero();
    nw_p.setZero();

    nw_p(2) = 1.;

    nu_p(0) =  std::cos(b_u);
    nu_p(1) = -std::sin(b_u);

    nv_p(0) =  -std::sin(g_v)*std::cos(a_v);
    nv_p(1) =   std::cos(g_v)*std::cos(a_v);
    nv_p(2) =   std::sin(a_v);

    // orientation vectors in (x,y,z) coordinates
    n_u_ = Rz_ * Ry_ * Rx_ * nu_p;
    n_v_ = Rz_ * Ry_ * Rx_ * nv_p;
    n_w_ = Rz_ * Ry_ * Rx_ * nw_p;

    //absolute probe orientation angles

  }

 };


}  // namespace Bembel
#endif



/*
//constructor based on a calibration file
ThreeAxisHallProbe(const std::string &rel_orientation_file,
                   const std::string &nonlinearity_file, const char del = ",") {

   std::ifstream file;
   file.open(rel_orientation_file);
   std::string current_line;
   std::string current_element;

   std::vector<double> tmp_vector;

   if (!file) {
     std::cerr << "File " << rel_orientation_file << " doesn't exist!";
   }
   else{
     // first row is irrelevant
     getline(file, current_line);

     //second row stored relative orientation
     getline(file, current_line);

     //covert this line
     while(std::getline(current_line,current_element,del)){
         tmp_vector.push_back(atof(current_element.c_str()));

     }

     //setup the relative probe orientation
     phi_uv_ = tmp_vector[0]*BEMBEL_PI/180;
     phi_uw_ = tmp_vector[1]*BEMBEL_PI/180;
     phi_vw_ = tmp_vector[2]*BEMBEL_PI/180;

     //update the absolute probe orientation
     update_absolute_orientation();


   }
 }
 */
