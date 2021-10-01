// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_DEVICES_NearFieldMatrixBlock_H__
#define __BEMBEL_DEVICES_NearFieldMatrixBlock_H__


/**
 *  \class UniqueList
 *  \brief
 *  \todo
 **/

namespace Bembel {

class NearFieldMatrixBlock {
 public:
   NearFieldMatrixBlock(const Eigen::MatrixXd &in_mat, const std::vector<int> indices) {

     mat_ = in_mat;
     dof_table_ = indices;
   }


   Eigen::VectorXd matvec(const Eigen::VectorXd &rhs){

     return mat_*dof_table_.unaryExpr(rhs);

   }




 private:
   Eigen::MatrixXd mat_;
   Eigen::VectorXi dof_table_;


};


}  // namespace Bembel

#endif
