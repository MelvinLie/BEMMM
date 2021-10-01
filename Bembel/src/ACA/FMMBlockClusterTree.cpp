// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_FMMBLOCKCLUSTERTREE_H_
#define BEMBEL_FMMBLOCKCLUSTERTREE_H_

#include <Eigen/Dense>


namespace Bembel {

 template <typename LinOp>
 class FMMBockClusterTree {

 public:
   //constructor
   FMMBockClusterTree(){};


   void init_BlockClusterTree(AnsatzSpace<LinOp> *ansatz_space,MeasurementData *meas_data, int max_tree_lvl, const Eigen::Matrix<double,2,3> &bbox){

     max_tree_lvl_ = max_tree_lvl;

     el_tree_.set_min_numel(min_numel_);
     el_tree_.init_ElementOctTree(ansatz_space,max_tree_lvl,bbox);
     el_tree_.print_tree();

     //meas_data_ = meas_data;

     meas_data_->set_min_numel(min_numel_);
     meas_data_->set_bounding_box(bbox);
     meas_data_->init_cluster_tree(max_tree_lvl_);

     //meas_data_->print_cluster_tree();

     /*
     std::vector<int> levels = meas_data_->get_levels();


     std::cout << "levels :" << std::endl;

     for (int i = 0; i < levels.size(); ++i){

       std::cout << levels[i] << " , ";

     }
     std::cout << std::endl;
     */
   }

 private:

   ElementOctTree<LinOp> el_tree_;
   MeasurementData *meas_data_;

   int max_tree_lvl_;
   int min_numel_ = 10;

}; //FMMBockClusterTree

}  // namespace Bembel
#endif
