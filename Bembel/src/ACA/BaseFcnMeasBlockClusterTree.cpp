// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_BASEFCNMEASBLOCKCLUSTERTREE_H_
#define BEMBEL_BASEFCNMEASBLOCKCLUSTERTREE_H_

#include <Eigen/Dense>


namespace Bembel {

 template <typename Derived>
 class BaseFcnMeasBockClusterTree {

 public:
   //constructor
   BaseFcnMeasBockClusterTree() {};

   BaseFcnMeasBockClusterTree(MeasurementData *meas_data, AnsatzSpace<Derived> ansatz_space,int max_level){

     max_level_ = max_level;

     //measurement tree
     meas_data->init_cluster_tree(max_level_);
     meas_cluster_ = std::addressof(meas_data->get_tree_memory()->get_root());


     //dof tree
     auto element_tree = ansatz_space.get_superspace().get_mesh().get_element_tree();
     el_cluster_ = std::addressof(element_tree.root());

     // set parameters for matrix assembly
     //max_level_ = el_cluster_->memory_->max_level_;
     dof_tree_ = BaseFcnTree<Derived>(&ansatz_space,max_level_);
     dof_cluster_ = std::addressof(dof_tree_.get_tree_memory()->get_root());


     polynomial_degree_ =  ansatz_space.get_superspace().get_polynomial_degree();
     polynomial_degree_plus_one_squared_ = (polynomial_degree_ + 1)
                                                  *(polynomial_degree_ + 1);

  };


   void init_BlockClusterTree(){

     // set up leaf_pointers
     leaf_pointers_ = std::make_shared<std::vector<BaseFcnMeasBockClusterTree*>>();
     leaf_pointers_->clear();

     is_leaf_ = false;

     lvl_ = 0;
     int counter = 0;


     appendSubtree(dof_cluster_,meas_cluster_, &counter);
     update_leaf_pointers();


   }


   void appendSubtree(const BaseFcnTreeNode *dof_cluster,
                      const MeasurementTreeNode *meas_cluster,
                      int *counter){


        std::shared_ptr<BaseFcnTreeMemory> dof_mem = dof_cluster->memory_;
        std::shared_ptr<MeasurementTreeMemory> meas_mem = meas_cluster->memory_;

        if (lvl_ == 0) cc_ = 0;
        else{
            cc_ = compareCluster(*dof_cluster, *meas_cluster);
        }

        if (cc_ == 0) {
          //Refine

          //std::cout << "refine"<< std::endl;
          // reserve memory for sons_
          sons_.resize(dof_cluster->sons_.size(), meas_cluster->sons_.size());

          for (auto j = 0; j < sons_.cols(); ++j){
            //std::cout << "j = " << j << std::endl;
            for (auto i = 0; i < sons_.rows(); ++i) {
              //std::cout << "i = " << i << std::endl;

              //this is not a leaf
              is_leaf_ = false;

              sons_(i, j).set_level(lvl_+1);

              sons_(i, j).set_parameters( eta_, max_level_, min_cluster_level_);

              const BaseFcnTreeNode &dof_son = dof_mem->son(*dof_cluster, i);
              sons_(i, j).dof_cluster_ = std::addressof(dof_son);

              const MeasurementTreeNode &meas_son = meas_mem->son(*meas_cluster, j);
              //std::cout << "address of son " << std::addressof(meas_son) << std::endl;
              sons_(i, j).meas_cluster_ = std::addressof(meas_son);

              //set leaf pointers
              sons_(i, j).set_leaf_pointers(leaf_pointers_);

              sons_(i, j).appendSubtree(std::addressof(dof_son), std::addressof(meas_son),counter);

            }
          }
        }
          else {
           is_leaf_ = true;
           leaf_ctr_ = *counter;
           ++(*counter);
           leaf_ = std::make_shared<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>();
         }

         return;


    }

    void print_block_cluster_tree(int max_level = -1){

      if(max_level == -1) max_level = lvl_;
      std::stringstream indent;

      for (int l = 0; l < max_level; ++l) indent << "\t";

      int num_dofs = dof_cluster_->indices_.size();
      int num_meas = meas_cluster_->indices_.size();

      int rows_son = sons_.rows();
      int cols_son = sons_.cols();

      std::cout << indent.str() << "-----------" << std::endl;
      std::cout << indent.str() << "Leaf number (" << leaf_ctr_<< ")" << std::endl;
      std::cout << indent.str() << "block size = ( " << num_meas << " , "<< num_dofs << " )" << std::endl;
      std::cout << indent.str() << "Number of sons = " << rows_son*cols_son << std::endl;
      std::cout << indent.str() << "dofs = ";
      for(int i = 0; i < num_dofs ; ++i) std::cout << dof_cluster_->indices_[i] << " , ";
      std::cout << std::endl;
      std::cout << indent.str() << "meas = ";
      for(int i = 0; i < num_meas ; ++i) std::cout << meas_cluster_->indices_[i] << " , ";
      std::cout << std::endl;

      if (is_leaf_){
        std::cout << indent.str() << "This is a ";
        if (cc_ == 1){
          std::cout << "low rank leaf" << std::endl;
        }
        else{
          std::cout  << "dense leaf" << std::endl;
        }
      }
      else{
        std::cout << indent.str() << "This is no leaf" << std::endl;
      }

      std::cout << indent.str() << "-----------" << std::endl;

      ++leaf_ctr_;

      for (int j = 0; j < sons_.cols(); ++j){
        //std::cout << "j = " << j << std::endl;
        for (int i = 0; i < sons_.rows(); ++i) {
          sons_(i,j).print_block_cluster_tree();

      }
    }


    }

    void set_eta(double eta){

      eta_ = eta;
    }

    void set_parameters(double eta,int max_level, int min_cluster_level){

      eta_ = eta;
      max_level_ = max_level;
      min_cluster_level_ = min_cluster_level;

    }


    const BaseFcnTreeNode *get_dof_cluster() const {
      return dof_cluster_;
    }

    const MeasurementTreeNode *get_measurement_cluster() const {
      return meas_cluster_;
    }


    void set_level(int level){
      lvl_ = level;
    }


    void set_leaf_pointers(std::shared_ptr<std::vector<BaseFcnMeasBockClusterTree *>> lp){

      leaf_pointers_ = lp;
    }

   std::vector<int> get_measurement_indices(){

     return meas_cluster_->indices_;

   }

   BaseFcnTree<Derived> *get_BaseFcnTree_ptr(){

     return &dof_tree_;
   }



   //members
   double eta_ = 1.6;
   MeasurementData *meas_data_;
   BaseFcnTree<Derived> dof_tree_;
   //ElementTree *element_tree_;
   const BaseFcnTreeNode *dof_cluster_;
   const ElementTreeNode *el_cluster_;
   const MeasurementTreeNode *meas_cluster_;
   //std::shared_ptr<ElementTreeMemory> el_mem_;
   //std::shared_ptr<MeasurementTreeMemory> meas_mem_;
   int max_level_;
   int min_cluster_level_;
   int polynomial_degree_;
   int polynomial_degree_plus_one_squared_;
   //int ffield_deg_ = 2;
   int lvl_;
   int cc_;
   int leaf_ctr_ = -1;
   int rows_,cols_;
   bool is_leaf_;
   std::shared_ptr<std::vector<BaseFcnMeasBockClusterTree *>> leaf_pointers_;
   GenericMatrix<BaseFcnMeasBockClusterTree> sons_;
   std::shared_ptr<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> leaf_;

 private:
  /**
   *  \brief determines admissible clusters
   **/


  int compareCluster(const BaseFcnTreeNode &dof_cluster,
                     const MeasurementTreeNode &meas_cluster) {

    int retval = 0;

    double dist = (dof_cluster.center_ - meas_cluster.center_).norm();

    double max_diam = dof_cluster.diam_ ;
    if (meas_cluster.diam_ > max_diam)  max_diam = meas_cluster.diam_;

    //std::cout << "dist = " << dist << std::endl;
    //std::cout << "max_diam = " << max_diam << std::endl;

    if (eta_*dist < max_diam){
        //refine
        retval = 0;
    }
    else{
        //low rank
        //std::cout << "Low Rank" << std::endl;
        retval = 1;

    }

    //std::cout << "lhs " << max_level_ - dof_cluster.level_  << std::endl;
    //std::cout << "rhs " << min_cluster_level_ << std::endl;
    if(max_level_ == lvl_){//(retval == 0) &&
        //dense
        //std::cout << "Dense" << std::endl;
        retval = 2;
    }

    return retval;
  }


  void update_leaf_pointers() {
    // make sure we have the root of the block cluster tree calling this method
    if (dof_cluster_->pos_ != 0 || meas_cluster_->pos_ != 0) {
      return;
    } else {
      leaf_pointers_->clear();
      update_leaf_pointers_recursion(*this);
    }
  }
  void update_leaf_pointers_recursion(BaseFcnMeasBockClusterTree &child) {

    if (child.cc_ == 0) {

      for (auto j = 0; j < child.sons_.cols(); ++j){
        for (auto i = 0; i < child.sons_.rows(); ++i){
          update_leaf_pointers_recursion(child.sons_(i, j));
          }
      }
    } else

      leaf_pointers_->push_back(std::addressof(child));
    return;
  }




 };


}  // namespace Bembel
#endif
