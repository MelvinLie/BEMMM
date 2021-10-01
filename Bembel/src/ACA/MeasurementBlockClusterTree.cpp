// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_MEASUREMENTBLOCKCLUSTERTREE2_H_
#define BEMBEL_MEASUREMENTBLOCKCLUSTERTREE2_H_

#include <Eigen/Dense>


namespace Bembel {

 //template <typename Scalar>
 class MeasurementBockClusterTree {

 public:
   //constructor
   MeasurementBockClusterTree(){   };


   void set_measurement_data(MeasurementData *meas_data){
     //meas_data_ = meas_data_;
     meas_cluster_ = std::addressof(meas_data->get_tree_memory()->get_root());
     //meas_mem_ = meas_cluster_->memory_;
   }

   template <typename Derived>
   void setup_element_tree(const AnsatzSpace<Derived> &ansatz_space){
     auto element_tree =
         ansatz_space.get_superspace().get_mesh().get_element_tree();
     el_cluster_ = std::addressof(element_tree.root());
     // get memory structure from element tree
     //el_mem_ = el_cluster_->memory_;
     // set parameters for matrix assembly
     max_level_ = el_cluster_->memory_->max_level_;

     polynomial_degree_ =
         ansatz_space.get_superspace().get_polynomial_degree();
     polynomial_degree_plus_one_squared_ = (polynomial_degree_ + 1)
                                             *(polynomial_degree_ + 1);
   }


   void init_BlockClusterTree(){

     // set up leaf_pointers
     leaf_pointers_ = std::make_shared<std::vector<MeasurementBockClusterTree *>>();
     leaf_pointers_->clear();

     is_leaf_ = false;

     //leafX_ = std::make_shared<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>();
     //leafY_ = std::make_shared<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>();
     //leafZ_ = std::make_shared<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>>();

     lvl_ = 0;
     int counter = 0;

     appendSubtree(el_cluster_,meas_cluster_, &counter);
     update_leaf_pointers();


   }

   void appendSubtree(const ElementTreeNode *el_cluster,
                      const MeasurementTreeNode *meas_cluster,
                      int *counter){


        std::shared_ptr<ElementTreeMemory> el_mem = el_cluster->memory_;
        std::shared_ptr<MeasurementTreeMemory> meas_mem = meas_cluster->memory_;



        cc_ = compareCluster(*el_cluster, *meas_cluster);


        if (cc_ == 0) {
          //Refine

          //std::cout << "refine"<< std::endl;
          // reserve memory for sons_
          sons_.resize(el_cluster->sons_.size(), meas_cluster->sons_.size());

          for (auto j = 0; j < sons_.cols(); ++j){
            //std::cout << "j = " << j << std::endl;
            for (auto i = 0; i < sons_.rows(); ++i) {
              //std::cout << "i = " << i << std::endl;

              //this is not a leaf
              is_leaf_ = false;

              sons_(i, j).set_level(lvl_+1);
              sons_(i, j).set_parameters( eta_, max_level_, min_cluster_level_);

              const ElementTreeNode &el_son = el_mem->son(*el_cluster, i);
              sons_(i, j).el_cluster_ = std::addressof(el_son);

              const MeasurementTreeNode &meas_son = meas_mem->son(*meas_cluster, j);
              //std::cout << "address of son " << std::addressof(meas_son) << std::endl;
              sons_(i, j).meas_cluster_ = std::addressof(meas_son);

              //set leaf pointers
              sons_(i, j).set_leaf_pointers(leaf_pointers_);

              sons_(i, j).appendSubtree(std::addressof(el_son), std::addressof(meas_son),counter);

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


    void print_block_cluster_tree(){


      std::stringstream indent;
      for (int l = 0; l < lvl_; ++l) indent << "\t";

      std::shared_ptr<ElementTreeMemory> el_mem = el_cluster_->memory_;

      int num_elements = std::distance(el_mem->cluster_begin(*el_cluster_),el_mem->cluster_end(*el_cluster_));
      int num_meas = meas_cluster_->indices_.size();

      int rows_son = sons_.rows();
      int cols_son = sons_.cols();

      std::cout << indent.str() << "-----------" << std::endl;
      std::cout << indent.str() << "Leaf number (" << leaf_ctr_<< ")" << std::endl;
      std::cout << indent.str() << "block size = ( " << num_meas << " , "<< num_elements << " )" << std::endl;
      std::cout << indent.str() << "Number of sons = " << rows_son*cols_son << std::endl;

      if (is_leaf_){
        std::cout << indent.str() << "This is a ";
        if (leaf_->is_low_rank()){
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

    const ElementTreeNode *get_element_cluster() const {
      return el_cluster_;
    }

    const MeasurementTreeNode *get_measurement_cluster() const {
      return meas_cluster_;
    }

    void set_level(int level){
      lvl_ = level;
    }

    void set_leaf_pointers(std::shared_ptr<std::vector<MeasurementBockClusterTree *>> lp){

      leaf_pointers_ = lp;
    }

    int get_row_start_index() {
      return polynomial_degree_plus_one_squared_*
             std::distance(el_cluster_->get_memory()->cpbegin(),
                           el_cluster_->get_memory()->cluster_begin(*el_cluster_));
    }
    int get_row_end_index() {

      return polynomial_degree_plus_one_squared_ *
             std::distance(el_cluster_->get_memory()->cpbegin(),
                           el_cluster_->get_memory()->cluster_end(*el_cluster_));
    }

   std::vector<int> get_measurement_indices(){

     return meas_cluster_->indices_;

   }

   void set_dof_table(std::vector<int> list){

     Eigen::Map<Eigen::VectorXi> tmp(&list[0], list.size());
     dof_table_ = tmp;

   }

   Eigen::VectorXi get_dof_table(){

     return dof_table_;
   }


   //members
   double eta_ = 1.6;
   MeasurementData *meas_data_;
   ElementTree *element_tree_;
   const ElementTreeNode *el_cluster_;
   const MeasurementTreeNode *meas_cluster_;
   //std::shared_ptr<ElementTreeMemory> el_mem_;
   //std::shared_ptr<MeasurementTreeMemory> meas_mem_;
   int max_level_;
   int min_cluster_level_;                        // a element tree leaf has 4^min_cluster_level_ elements
   int polynomial_degree_;
   int polynomial_degree_plus_one_squared_;
   //int ffield_deg_ = 2;
   int lvl_;
   int cc_;
   int leaf_ctr_ = -1;
   int rows_,cols_;
   bool is_leaf_;
   std::shared_ptr<std::vector<MeasurementBockClusterTree *>> leaf_pointers_;
   GenericMatrix<MeasurementBockClusterTree> sons_;
   std::shared_ptr<TreeLeaf<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>> leaf_;
   Eigen::VectorXi dof_table_;

 private:
  /**
   *  \brief determines admissible clusters
   **/

  int compareCluster(const ElementTreeNode &el_cluster,
                     const MeasurementTreeNode &meas_cluster) {

    int retval = 0;

    //the root node in el_cluster has inf radius
    if(std::isinf(el_cluster.radius_ )) return 0;

    double dist = (el_cluster.midpoint_ - meas_cluster.center_).norm();

    double max_diam = el_cluster.radius_ ;
    if (meas_cluster.diam_ > max_diam)  max_diam = meas_cluster.diam_;

    if (eta_*dist < max_diam){
        //refine
        retval = 0;
    }
    else{
        //low rank
        retval = 1;
    }

    //std::cout << "lhs " << max_level_ - el_cluster.level_  << std::endl;
    //std::cout << "rhs " << min_cluster_level_ << std::endl;
    //if((max_level_ - el_cluster.level_ <=  min_cluster_level_)){//(retval == 0) &&
    if(lvl_ ==  max_level_-1){
        //dense
        retval = 2;
    }

    return retval;
  }

  void update_leaf_pointers() {
    // make sure we have the root of the block cluster tree calling this method
    if (el_cluster_->id_ != -1 || meas_cluster_->pos_ != 0) {
      return;
    } else {
      leaf_pointers_->clear();
      update_leaf_pointers_recursion(*this);
    }
  }
  void update_leaf_pointers_recursion(MeasurementBockClusterTree &child) {

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
