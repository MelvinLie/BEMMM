// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_BINARYTREENODE_H_
#define BEMBEL_BINARYTREENODE_H_


namespace Bembel {

/**
 *  \ingroup MeasurementClusterNode
 *  \brief Cluster of measurements
 */

 // forward declaration of memory is necessary here
 struct BinaryTreeMemory;


 class BinaryTreeNode {

 public:
  //constructor
  BinaryTreeNode() { }

  void append_index(const int index){
    indices_.push_back(index);
  }


  void set_memory(std::shared_ptr<BinaryTreeMemory> memory) {
    memory_ = memory;
    return;
  }

  void setup_cluster_box(Eigen::MatrixXd bbox) {
    bbox_ = bbox;
    center_ = bbox.colwise().mean();
    diam_ = (bbox.row(0) - bbox.row(1)).norm();
    return;
  }

  void set_center(Eigen::Vector3d in){
    center_ = in;
  }




   std::vector<int> indices_;
   Eigen::MatrixXd bbox_;
   std::shared_ptr<BinaryTreeMemory> memory_;
   int father_;
   double diam_;
   std::vector<int> sons_;
   int pos_,level_;
   Eigen::Vector3d center_;


};

}  // namespace Bembel
#endif
