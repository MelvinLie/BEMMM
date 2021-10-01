// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_TRACKING_REFERENCETRAJECTORY_H_
#define BEMBEL_TRACKING_REFERENCETRAJECTORY_H_

namespace Bembel {
/**
 *  \ingroup Potential
 *  \brief DiscretePotential
 *  \todo  Add a documentation
 */
class ReferenceTrajectory {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  ReferenceTrajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen){

    arclen_.push_back(0.);
    path_.push_back(Trajectory(start_pos,direction,arclen));
    arclen_.push_back(arclen);

  }

  ReferenceTrajectory(const Eigen::Vector3d &start_pos,const Eigen::Vector3d &direction,
                        double arclen,double curvature){

    arclen_.push_back(0.);
    path_.push_back(Trajectory(start_pos,direction,arclen,curvature));
    arclen_.push_back(arclen);

  }

  void add_straight_section(double arclen){

    Eigen::Vector3d start_pos = path_.back().get_end();
    Eigen::Vector3d start_direction = path_.back().get_end_direction();


    path_.push_back(Trajectory(start_pos,start_direction,arclen));
    arclen_.push_back(arclen+arclen_.back());

  }

  void add_curve(double arclen,double curvature){

    Eigen::Vector3d start_pos = path_.back().get_end();
    Eigen::Vector3d start_direction = path_.back().get_end_direction();

    path_.push_back(Trajectory(start_pos,start_direction,arclen,curvature));
    arclen_.push_back(arclen+arclen_.back());

  }

  Eigen::Vector3d get_pos(double s){

    Eigen::Vector3d ret_val;

    for (int i = 0; i < arclen_.size()-1; ++i){
      if((s >= arclen_[i]) && (s < arclen_[i+1])){

        ret_val = path_[i].get_pos(s-arclen_[i]);

      }

    }

    return ret_val;

  }

  Eigen::Vector3d get_end(){

    return get_pos(get_pathlength());

  }

  Eigen::Vector3d get_dir(double s){

    Eigen::Vector3d ret_val;

    for (int i = 0; i < arclen_.size()-1; ++i){
      if((s >= arclen_[i]) && (s < arclen_[i+1])){

        ret_val = path_[i].get_dir(s-arclen_[i]);

      }

    }

    return ret_val;

  }

  Eigen::Matrix3d get_unit_vectors(double s){

    Eigen::Matrix3d ret_val;

    for (int i = 0; i < arclen_.size()-1; ++i){
      if((s >= arclen_[i]) && (s < arclen_[i+1])){

        ret_val = path_[i].get_unit_vectors(s-arclen_[i]);

      }

    }

    return ret_val;

  }

  double get_pathlength(){

    return arclen_.back();

  }

  double get_curvature(double s){

    double ret_val = 0.;

    for (int i = 0; i < arclen_.size()-1; ++i){
      if((s >= arclen_[i]) && (s < arclen_[i+1])){

        ret_val = path_[i].get_curvature();

      }

    }
    return ret_val;
  }

 private:
    std::vector<Trajectory> path_;
    std::vector<double> arclen_;

};  // namespace Bembel

}  // namespace Bembel
#endif
