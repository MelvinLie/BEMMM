// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_MEASUREMENTCLUSTERNODE_H_
#define BEMBEL_MEASUREMENTCLUSTERNODE_H_


namespace Bembel {

/**
 *  \ingroup MeasurementClusterNode
 *  \brief Cluster of measurements
 */

 // forward declaration of memory is necessary here
 struct MeasurementTreeMemory;


 class MeasurementTreeNode {

 public:
  //constructor
  MeasurementTreeNode() { }

  void append_index(const int index){
    indices_.push_back(index);
  }


  void set_memory(std::shared_ptr<MeasurementTreeMemory> memory) {
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

  void set_meas_table(std::vector<int> in){

      meas_table_ = Eigen::Map<Eigen::VectorXi>(&(in.data()[0]),indices_.size());
  }

  void make_meas_table(){
      meas_table_ = Eigen::Map<Eigen::VectorXi>(&(indices_.data()[0]),indices_.size());
  }


   std::vector<int> indices_;
   std::vector<int> interaction_region_;
   std::vector<int> interaction_m2l_index_;
   std::vector<int> interaction_m2l_rot_index_;
   std::vector<int> near_field_;
   Eigen::VectorXi meas_table_; //(we could in future remove indices_ since it stores the same info)
   Eigen::MatrixXd bbox_;
   std::shared_ptr<MeasurementTreeMemory> memory_;
   int father_;
   double diam_;
   std::vector<int> sons_;
   int pos_,level_;
   Eigen::Vector3d center_;

   //local FM expansion coefficients
   Eigen::VectorXcd locals_;
   //Eigen::MatrixXcd local_mat_;
   Eigen::SparseMatrix<std::complex<double>,Eigen::RowMajor> local_mat_;

};

}  // namespace Bembel
#endif
