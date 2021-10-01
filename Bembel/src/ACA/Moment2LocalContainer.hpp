// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_ACA_MOMENT2LOCALCONTAINER_H_
#define BEMBEL_ACA_MOMENT2LOCALCONTAINER_H_


namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief This class organizes an element structure on a Geometry object and
 * handles refinement.
 *
 * \todo Describe the ElementOctTree
 */
class Moment2LocalContainer {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  //ElementOctTree(AnsatzSpace<Derived> *ansatz_space) {}
  Moment2LocalContainer(){ };

  void fill(const int num_multipoles, const int max_level, const Eigen::Matrix<double,2,3> &bbox){

    num_multipoles_ = num_multipoles;

    Eigen::Vector3d center = 0.5*(bbox.row(0)+bbox.row(1)).transpose();

    dx_ = (bbox(1,0) - bbox(0,0))/4.;
    dy_ = (bbox(1,1) - bbox(0,1))/4.;
    dz_ = (bbox(1,2) - bbox(0,2))/4.;

    D_.setZero();

    Eigen::Vector3d this_pos;

    double hx = dx_;
    double hy = dy_;
    double hz = dz_;

    for(int depth = 0; depth <= max_level - 2 ; ++depth ){

      //running index
      int row = 0;
      std::vector<Eigen::MatrixXcd,Eigen::aligned_allocator<Eigen::MatrixXcd>> tmp;

      for(int i = -3; i <= 3; ++i){
        if(std::abs(i) < 2) continue;

        for(int j = -3; j < 3; ++j){
          if(std::abs(j) < 2) continue;

          for(int k = -3; k < 3; ++k){
            if(std::abs(k) < 2) continue;

            this_pos(0) = i*dx_;
            this_pos(1) = j*dy_;
            this_pos(2) = k*dz_;

            //We store the positions in the second level in the D matrix
            if(depth == 0){
              D_.row(row) = this_pos.transpose();
            }

            tmp.push_back(make_transformation_matrix(this_pos));

            ++row;

          }
        }
      }
      transformation_matrices_.push_back(tmp);

      //refine
      hx *= 0.5;
      hy *= 0.5;
      hz *= 0.5;

    }

  }

  Eigen::MatrixXcd get_transformation_matrix(const int level,const int index){

    return transformation_matrices_[level][index];

  }



  Eigen::VectorXcd m2l(const int level,const int index, const Eigen::VectorXcd &moments_in){

    return get_transformation_matrix(level, index)*moments_in;

  }

  int get_index(const int level, const Eigen::Vector3d &x){

    int ret_val = -1;

    double scale_fac = 1./(level-1.);

    Eigen::Vector3d d;

    for(int i = 0 ; i < 316 ; ++i){

      d = x - scale_fac*D_.row(i).transpose();

      if(d.norm() < 1e-10){
        return i;
      }
    }
    return ret_val;
  }


  private:

    std::vector<std::vector<Eigen::MatrixXcd,
                Eigen::aligned_allocator<Eigen::MatrixXcd>>> transformation_matrices_;

    int levels_;
    int num_multipoles_;
    double dx_,dy_,dz_;

    Eigen::Matrix<double,316,3> D_;

    Eigen::MatrixXcd make_transformation_matrix(const Eigen::Vector3d &x){


      Eigen::MatrixXcd ret_mat;
      ret_mat.setZero();

      //irregular solid harmonics
      Eigen::VectorXcd S_lm = Slm(num_multipoles_,  x);

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

                ret_mat(k_s,k_p) += fac*std::conj(S_lm(k));

              }
              ++k_p;
            }
          }
          ++k_s;
        }
      }

      return ret_mat;

    }

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
  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif
