// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_PERTURBATIONCOV_H_
#define BEMBEL_PERTURBATIONCOV_H_

#include <Eigen/Dense>


namespace Bembel {

/**
 *  \ingroup MeasurementData
 *  \brief Helper class that is used in order to input and store measurement data
 */

 // forward declaration of class MeasurementData
 class PerturbationCovariance;


 class PerturbationCovariance {

 public:
  //constructors
  PerturbationCovariance() {}

  //constructor reading the covariance information from file
  PerturbationCovariance(const std::string &file_name, const int axis_id) {

    read_from_file(file_name,axis_id);

  }

  //read the covariance information from file
  void read_from_file(const std::string &file_name, const int axis_id) {

    Eigen::MatrixXd data = read_csv(file_name,true);

    //number of samples available
    num_samples_ = data.rows();

    //max lag
    max_lag_ = data(num_samples_-1,0);

    //to identify the arm axis orienation
    axis_id_ = axis_id;

    //covariance vectors
    cov_x_ = data.col(1);
    cov_y_ = data.col(2);
    cov_z_ = data.col(3);
    cov_t1_ = data.col(4);
    cov_t2_ = data.col(5);
    cov_wt_1_ = data.col(6);
    cov_wt_2_ = data.col(7);


  }

  Eigen::MatrixXd get_prec(const int num_meas){

    //make space for matrix
    Eigen::MatrixXd L(5*num_meas,5*num_meas);
    L.setZero();

    //std::cout << "cov_x_ = " << cov_x_.size() << std::endl;
    //std::cout << "num_meas = " << num_meas << std::endl;
    //diagonal blocks
    L.block(0,0,num_meas,num_meas) = make_toeplitz(cov_x_,num_meas).inverse();
    L.block(num_meas,num_meas,num_meas,num_meas) = make_toeplitz(cov_y_,num_meas).inverse();
    L.block(2*num_meas,2*num_meas,num_meas,num_meas) = make_toeplitz(cov_z_,num_meas).inverse();
    L.block(3*num_meas,3*num_meas,num_meas,num_meas) = make_toeplitz(cov_t1_,num_meas).inverse();
    L.block(4*num_meas,4*num_meas,num_meas,num_meas) = make_toeplitz(cov_t2_,num_meas).inverse();


    return L;

  }

  Eigen::MatrixXd get_cov(const int num_meas){

    //make space for matrix
    Eigen::MatrixXd D(5*num_meas,5*num_meas);
    D.setZero();

    //diagonal blocks
    D.block(0,0,num_meas,num_meas) = make_toeplitz(cov_x_,num_meas);
    D.block(num_meas,num_meas,num_meas,num_meas) = make_toeplitz(cov_y_,num_meas);
    D.block(2*num_meas,2*num_meas,num_meas,num_meas) = make_toeplitz(cov_z_,num_meas);
    D.block(3*num_meas,3*num_meas,num_meas,num_meas) = make_toeplitz(cov_t1_,num_meas);
    D.block(4*num_meas,4*num_meas,num_meas,num_meas) = make_toeplitz(cov_t2_,num_meas);

    //found better performance without this
    if(false){
      //off diagonal blocks
      if(axis_id_ == 0){
        //arm is orientated along x
        //vibration id in xy and xz plane


        D.block(num_meas,3*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);
        D.block(3*num_meas,num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);

        D.block(2*num_meas,4*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);
        D.block(4*num_meas,2*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);


      }
      else if(axis_id_ == 1){
        //arm is orientated along y
        //vibration id in yx and yz plane

        D.block(0,3*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);
        D.block(3*num_meas,0,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);

        D.block(2*num_meas,4*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);
        D.block(4*num_meas,2*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);


      }
      else if(axis_id_ == 2){
        //arm is orientated along z
        //vibration id in zx and zy plane


        D.block(num_meas,3*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);
        D.block(3*num_meas,num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_1_,num_meas);

        D.block(0,3*num_meas,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);
        D.block(3*num_meas,0,num_meas,num_meas) =  make_toeplitz(cov_wt_2_,num_meas);

      }

    }
    return D;
  }

  void save_to_file(const std::string &filename, const int num_samples, const int prec = 8){


    Eigen::MatrixXd cov = get_cov(num_samples);

    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(prec, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream out_file(filename);
    out_file << cov.format(CSVFormat);
    out_file.close();

    return;

  }


  private:

    //vectors to store all relevant information about the perturbation covariance
    Eigen::VectorXd cov_x_;
    Eigen::VectorXd cov_y_;
    Eigen::VectorXd cov_z_;
    Eigen::VectorXd cov_t1_;
    Eigen::VectorXd cov_t2_;
    Eigen::VectorXd cov_wt_1_;
    Eigen::VectorXd cov_wt_2_;

    //number of samples
    int num_samples_;

    //maximum lag
    double max_lag_;

    //mapper axis id
    int axis_id_;


    Eigen::MatrixXd make_toeplitz(const Eigen::VectorXd &vec, const int num_meas){

      Eigen::MatrixXd ret_mat(num_meas,num_meas);
      ret_mat.setZero();

      Eigen::VectorXd vec_rev = vec.segment(0,num_meas).reverse().eval();

      for(int i = 0 ; i < num_meas; ++i){


        ret_mat.block(i,i,1,num_meas-i) = vec.segment(0,num_meas-i).transpose();
        ret_mat.block(i,0,1,i) = vec_rev.segment(num_meas-i-1,i).transpose();

      }
      return ret_mat;

    }


 };


}  // namespace Bembel
#endif
