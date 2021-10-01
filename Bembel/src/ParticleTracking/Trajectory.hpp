// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_TRACKING_TRAJECTORY_H_
#define BEMBEL_TRACKING_TRAJECTORY_H_

namespace Bembel {
/**
 *  \ingroup Potential
 *  \brief DiscretePotential
 *  \todo  Add a documentation
 */
class Trajectory {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  Trajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen){
     setup_trajectory(start_pos,direction,arclen);
  }

  Trajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen,double curvature){
     setup_trajectory(start_pos,direction,arclen,curvature);
  }


  void setup_trajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen){

    q0_ = start_pos;
    dir_ = direction;
    dir_ /= dir_.norm();
    arclen_  = arclen;
    is_curved_ = false;

  }
  void setup_trajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen,double curvature){

    q0_ = start_pos;
    dir_ = direction;
    dir_ /= dir_.norm();
    arclen_  = arclen;
    curvature_ = curvature;
    is_curved_ = true;


    phi0_ =  std::atan2(dir_(2),dir_(0)) + BEMBEL_PI/2.;


    ctr_ = q0_;

    ctr_(0) += dir_(2)/curvature;
    ctr_(2) -= dir_(0)/curvature;


  }

  Eigen::Vector3d get_pos(double s){

    Eigen::Vector3d ret_val;


    if(is_curved_){

      ret_val = ctr_;

      ret_val(0) += std::cos(s*curvature_ - phi0_)/curvature_;
      ret_val(2) -= std::sin(s*curvature_ - phi0_)/curvature_;

    }
    else{

      ret_val = s*dir_ + q0_;
    }

    return ret_val;

  }

  Eigen::Vector3d get_dir(double s){

    Eigen::Vector3d ret_val;
    ret_val.setZero();

    if(is_curved_){

      ret_val(0) = -std::sin(s*curvature_ - phi0_);
      ret_val(2) = -std::cos(s*curvature_ - phi0_);

    }
    else{

      ret_val = dir_;
    }

    return ret_val;

  }

  Eigen::Vector3d get_end(){

    return get_pos(arclen_);

  }

  Eigen::Vector3d get_end_direction(){

    return get_dir(arclen_);

  }

  Eigen::Matrix3d get_unit_vectors(double s){

    Eigen::Matrix3d ret_val;
    ret_val.setZero();

    ret_val.col(2) = get_dir(s);
    ret_val(1,1) = 1.;
    ret_val.col(0) = ret_val.col(1).cross(ret_val.col(2));

    return ret_val;

  }

  double get_curvature(){
    return curvature_;
  }

 private:
    Eigen::Vector3d q0_;
    Eigen::Vector3d ctr_;
    double arclen_;
    Eigen::Vector3d dir_;
    double curvature_;
    bool is_curved_;
    double phi0_;

};  // namespace Bembel

}  // namespace Bembel
#endif
