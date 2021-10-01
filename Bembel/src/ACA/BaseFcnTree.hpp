// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_CLUSTERTREE_BASEFCNTREE_H_
#define BEMBEL_CLUSTERTREE_BASEFCNTREE_H_


namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief This class organizes an element structure on a Geometry object and
 * handles refinement.
 *
 * \todo Describe the BaseFcnTree
 */
 template <typename Derived>
class BaseFcnTree {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  BaseFcnTree() {}
  BaseFcnTree(AnsatzSpace<Derived> *ansatz_space, unsigned int max_level) {
    init_BaseFcnTree(ansatz_space, max_level);
  }
  //////////////////////////////////////////////////////////////////////////////
  /// init
  //////////////////////////////////////////////////////////////////////////////
  void init_BaseFcnTree(AnsatzSpace<Derived> *ansatz_space, unsigned int max_level) {

    // initialise the data fields
    ansatz_space_ = ansatz_space;

    //num dofs
    num_dofs_ = ansatz_space_->get_number_of_dofs();

    T_matrix_ = ansatz_space_->get_transformation_matrix();

    init_dof_table();
    //std::cout << "Cog 0 = " << get_cog(0) << std::endl;

    //allocate memory pointer
    mem_ = std::make_shared<BaseFcnTreeMemory>();
    //set max level
    mem_->max_level_ = max_level;
    //allocate a pointer to the tree nodes
    mem_->memory_ = std::make_shared<std::vector<BaseFcnTreeNode>>();
    //allocate memory for the hole tree
    mem_->memory_->resize(mem_->cumNumElements(max_level));


    //container for domain bisection
    std::vector<Eigen::MatrixXd> domain_bisection;

    //compute the root bounding box
    init_bounding_box();

    //get reference to the root node
    BaseFcnTreeNode &root = mem_->get_root();
    //initialize the root node measurement indices
    for (int i = 0; i < num_dofs_ ; ++i) root.append_index(i);
    //set memory pointer of root node
    root.set_memory(mem_);
    //set root nodes center and bounding box
    root.setup_cluster_box(bounding_box_);
    //set level and position in memory of root node
    root.pos_ = 0;
    root.level_ = -1; //root is defined as level -1

    //setup level counter to -1
    lvl_ctr_ = -1;

    //setup the total node counter
    node_ctr_ = 1;

    generate_cluster_tree(0,1);

    mem_->son(root,0);

    return;
  }

  std::shared_ptr<BaseFcnTreeMemory> get_tree_memory(){

    return mem_;
  }

  std::vector<std::vector<int>> *get_element_index_table_ptr(){
    return &element_index_table_;
  }

  std::vector<std::vector<int>> *get_local_index_table_ptr(){
    return &local_index_table_;
  }

  std::vector<std::vector<double>> *get_weight_table_ptr(){
    return &weight_table_;
  }

  /*
  const std::vector<local_bernstein_combination> *get_local_bernstein_combination_ptr(){


  }
  */


  void print_tree(){

    BaseFcnTreeNode tmp_node;

    if((*mem_->memory_).size() == 0){

      std::cout << "Initialize a cluster tree first!" << std::endl;
    }
    else{
        for (int n = 0; n < node_ctr_; ++n){
          tmp_node = (*mem_->memory_)[n];
          std::cout << "------------------------------" << std::endl;
          std::cout << "Node (" << n << "):" << std::endl;
          std::cout << "\tpos =" << tmp_node.pos_ << std::endl;
          std::cout << "\tCenter = ( " << tmp_node.center_(0) << " , ";
          std::cout  << tmp_node.center_(1) << " , " << tmp_node.center_(2) << " )" << std::endl;
          std::cout << "\tDiam = " << tmp_node.diam_ << std::endl;
          std::cout << "\tBbox = ( " << tmp_node.bbox_(0,0) << " , ";
          std::cout  << tmp_node.bbox_(0,1) << " , " << tmp_node.bbox_(0,2) << " )" << std::endl;
          std::cout << "\t       ( " << tmp_node.bbox_(1,0) << " , ";
          std::cout  << tmp_node.bbox_(1,0) << " , " << tmp_node.bbox_(1,2) << " )" << std::endl;
          std::cout << "\tNumber of dofs = " << tmp_node.indices_.size() << std::endl;
          std::cout << "\tLevel = " << tmp_node.level_ << std::endl;
          std::cout << "\tNumber of sons = " << tmp_node.sons_.size() << std::endl;
          std::cout << "\tsons = ";
          for (int i = 0 ; i <  tmp_node.sons_.size() ; ++i) std::cout << tmp_node.sons_[i] << " , ";
          std::cout << std::endl;
          std::cout << "\tDoFs = ";
                    for (int i = 0 ; i <  tmp_node.indices_.size() ; ++i) std::cout << tmp_node.indices_[i] << " , ";
          std::cout << std::endl;
        }

    }
  }

  std::vector<BSplineBasisCombination> make_BSpline_table(){

    std::vector<BSplineBasisCombination> bspline_table;
    bspline_table.resize(num_dofs_);


    int I2 = (ansatz_space_->get_polynomial_degree()+1)*(ansatz_space_->get_polynomial_degree()+1);


    std::vector<std::vector<double>> weight_table;
    std::vector<std::vector<int>> element_index_table;
    std::vector<std::vector<int>> local_index_table;

    //resize tables
    element_index_table.resize(num_dofs_);
    local_index_table.resize(num_dofs_);
    weight_table.resize(num_dofs_);

    int el_indx;

    init_element_id_table();

    for (int k=0; k < T_matrix_.outerSize(); ++k){

      for (Eigen::SparseMatrix<double>::InnerIterator it(T_matrix_,k); it; ++it)
      {

        el_indx = it.row() / I2;

        bspline_table[k].append_bernstein_basis(element_id_table_[el_indx],it.row() % I2,it.value());

      }
    }

    /*
    for (int i = 0; i < bspline_table.size() ; ++i ){

      std::cout << "*******************" << std::endl;
      std::cout << "DoF "  << i << std::endl;
      bspline_table[i].print();

    }
    */

    return bspline_table;

  }

  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  AnsatzSpace<Derived> *ansatz_space_;
  std::shared_ptr<BaseFcnTreeMemory> mem_;
  int number_of_basis_fcns_;
  int max_level_;
  int num_dofs_;
  int number_of_points_;
  Eigen::Matrix<double,2,3> bounding_box_;
  int lvl_ctr_;
  int node_ctr_;
  Eigen::SparseMatrix<double> T_matrix_;
  std::vector<std::vector<double>> weight_table_;
  std::vector<std::vector<int>> element_index_table_;
  std::vector<std::vector<int>> local_index_table_;
  std::vector<std::vector<local_bernstein_combination>> local_bernstein_table_;
  std::vector<int> element_id_table_;

  void init_bounding_box(){

    PatchVector geometry = ansatz_space_->get_superspace().get_geometry();

    double x_min = std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity();
    double z_min = std::numeric_limits<double>::infinity();
    double x_max = -1*std::numeric_limits<double>::infinity();
    double y_max = -1*std::numeric_limits<double>::infinity();
    double z_max = -1*std::numeric_limits<double>::infinity();

    bounding_box_.setZero();

    //corner points in 3d
    Eigen::Vector3d c1,c2,c3,c4;
    Eigen::Vector3d cntr;
    Eigen::Vector3d eval_ptr;

    int steps_eval = 10;
    double u,v;


    for(int j = 0; j < geometry.size(); ++j){

      for(int i = 0; i < steps_eval;++i){

        u =  ((double) i)/(steps_eval-1);

        for(int k = 0; k < steps_eval;++k){

          v =  ((double) k)/(steps_eval-1);

          cntr = geometry[j].eval(u,v);


          if (cntr(0) < x_min) x_min = cntr(0);
          if (cntr(1) < y_min) y_min = cntr(1);
          if (cntr(2) < z_min) z_min = cntr(2);

          if (cntr(0) > x_max) x_max = cntr(0);
          if (cntr(1) > y_max) y_max = cntr(1);
          if (cntr(2) > z_max) z_max = cntr(2);
        }

      }


    }
    Eigen::Matrix<double,2,3> bounding_box;

    bounding_box_(0,0) = x_min;
    bounding_box_(1,0) = x_max;
    bounding_box_(0,1) = y_min;
    bounding_box_(1,1) = y_max;
    bounding_box_(0,2) = z_min;
    bounding_box_(1,2) = z_max;


    //margin for outer bounding box
    double marg_ratio = 0.1;
    double marg = marg_ratio*(bounding_box_.row(1)-bounding_box_.row(0)).maxCoeff();

    bounding_box_(0,0) -= marg;
    bounding_box_(0,1) -= marg;
    bounding_box_(0,2) -= marg;
    bounding_box_(1,0) += marg;
    bounding_box_(1,1) += marg;
    bounding_box_(1,2) += marg;

    return;

  }

  /*
    */


  void generate_cluster_tree(int mem_from,int mem_to){

    //increment level counter
    lvl_ctr_ += 1;
    //number of current root nodes
    int num_root = mem_to - mem_from;
    //position of first son in memory = mem_to;
    //initialize son counter
    int son_ctr = 0;
    //bisected indices list
    std::vector<std::vector<int>> bs_idx_list;
    //container for domain bisection
    std::vector<Eigen::MatrixXd> domain_bisection;

    //iterate over all given nodes
    for (int n = 0; n < num_root; ++n){

      //bisect the domain and separate the index list of this node
      bs_idx_list = bisect_index_list((*mem_->memory_)[mem_from + n].bbox_,
                                        (*mem_->memory_)[mem_from + n].indices_,
                                        &domain_bisection);

      for (int s = 0; s < bs_idx_list.size(); ++s){
        //it can be that there are no dofs in this subsection!
        if (bs_idx_list[s].size() == 0) continue;
        else{
          //there are some measurements here, this is a node
          //increment the total node counter
          node_ctr_++;
          //initialize the new node
          (*mem_->memory_)[mem_to + son_ctr].indices_ = bs_idx_list[s];
          //give it access to memory
          (*mem_->memory_)[mem_to + son_ctr].set_memory(mem_);
          //setup his father
          (*mem_->memory_)[mem_to + son_ctr].father_ = mem_from + n;
          //this the geometric information
          (*mem_->memory_)[mem_to + son_ctr].setup_cluster_box(domain_bisection[s]);
          //mark level in level indicator
          (*mem_->memory_)[mem_to + son_ctr].level_ = lvl_ctr_;
          //position of this node in memory_
          (*mem_->memory_)[mem_to + son_ctr].pos_ = node_ctr_-1;
          //add this position for the son to father
          (*mem_->memory_)[mem_from + n].sons_.push_back(mem_to + son_ctr);
          //increment son counter
          son_ctr++;
        }

        }


      }
    if (lvl_ctr_ < mem_->max_level_){
      //refine further
      generate_cluster_tree(mem_to, mem_to + son_ctr);

    }

    }



    std::vector<std::vector<int>> bisect_index_list(Eigen::MatrixXd box, std::vector<int> index_list, std::vector<Eigen::MatrixXd> *bisection){

      double bound_marg = 1e-8;

      double x_sect[3] = {box(0,0)-bound_marg,0.5*(box(1,0)+box(0,0)),box(1,0)+bound_marg};
      double y_sect[3] = {box(0,1)-bound_marg,0.5*(box(1,1)+box(0,1)),box(1,1)+bound_marg};
      double z_sect[3] = {box(0,2)-bound_marg,0.5*(box(1,2)+box(0,2)),box(1,2)+bound_marg};

      std::vector<std::vector<int>> ret_list;

      //we separate into 8 subdomains
      ret_list.resize(8);
      bisection->resize(8);
      for (int i = 0; i < 8; ++i) (*bisection)[i].resize(2,3);

      //temporal index list with remaining indices
      std::vector<int> rem_indx_list;// = index_list;

      //center of gravity of basis function
      Eigen::Vector3d cog;

      //index counter
      int current_index;

      //domain counter
      int l = 0;
      for (int i = 0; i < 2; ++i){

        for (int j = 0; j < 2; ++j){
          for(int k = 0; k < 2; ++k){

              (*bisection)[l](0,0) = x_sect[i];
              (*bisection)[l](1,0) = x_sect[i+1];
              (*bisection)[l](0,1) = y_sect[j];
              (*bisection)[l](1,1) = y_sect[j+1];
              (*bisection)[l](0,2) = z_sect[k];
              (*bisection)[l](1,2) = z_sect[k+1];



              for(int m = 0; m < index_list.size(); m++){
                current_index = index_list[m];


                //std::cout << "dof = " << current_index<< std::endl;
                //std::cout << "\tcog = " << cog.transpose() << std::endl;

                cog = get_cog(current_index);



                if((x_sect[i] <= cog(0)) & (cog(0) < x_sect[i+1]) &
                   (y_sect[j] <= cog(1)) & (cog(1) < y_sect[j+1]) &
                   (z_sect[k] <= cog(2)) & (cog(2) < z_sect[k+1])){
                     //std::cout << "\ttake" << std::endl;

                     //append index to index list
                     ret_list[l].push_back(current_index);

                   }
                   else{
                    //std::cout << "\treject" << std::endl;
                    //add index to the remaining ones
                    rem_indx_list.push_back(current_index);

                   }

              }
              //hand over remaining index list to avoid double indices
              index_list = rem_indx_list;
              rem_indx_list.clear();
              l++;


              //for(int ii = 0; ii < index_list.size(); ++ii) std::cout << index_list[ii] << " , ";
              //std::cout << std::endl;
          }
        }
      }
      if(rem_indx_list.size() != 0){
        std::cout << "!!!Something went wrong! Not all DoFs could be clustered." << std::endl;
        std::cout << "Remaining DoFs = ";
        for(int i = 0; i < rem_indx_list.size(); ++i){
          std::cout << rem_indx_list[i] << " , ";
        }
        std::cout << std::endl;
      }
      return ret_list;
    }



    void init_dof_table(){



      int I2 = (ansatz_space_->get_polynomial_degree()+1)*(ansatz_space_->get_polynomial_degree()+1);



      //resize tables
      element_index_table_.resize(num_dofs_);
      local_index_table_.resize(num_dofs_);
      weight_table_.resize(num_dofs_);

      int el_indx;

      init_element_id_table();

      for (int k=0; k < T_matrix_.outerSize(); ++k){


        for (Eigen::SparseMatrix<double>::InnerIterator it(T_matrix_,k); it; ++it)
        {

          el_indx = it.row() / I2;

          element_index_table_[k].push_back(element_id_table_[el_indx]);
          local_index_table_[k].push_back(it.row() % I2);
          weight_table_[k].push_back(it.value());
        }
      }
    }

    void init_element_id_table(){


      //get the element tree
      ElementTree el_tree = ansatz_space_->get_superspace().get_mesh().get_element_tree();

      element_id_table_.clear();

      //std::cout << "Element -> ID table = ";
      for (auto element = el_tree.cpbegin(); element != el_tree.cpend(); ++element){

        element_id_table_.push_back(element->id_);

        //std::cout << element_id_table_.back() << " , ";
      }
      //std::cout << std::endl;

      std::sort(element_id_table_.begin(),element_id_table_.end());
    }

    Eigen::Vector3d get_cog(int dof_indx){

      Eigen::Vector3d cog;
      cog.setZero();

      //get the element tree
      ElementTree el_tree = ansatz_space_->get_superspace().get_mesh().get_element_tree();

      //pointer to first element in lowest level of element tree
      std::vector<ElementTreeNode>::const_iterator element = el_tree.cpbegin();

      int num_elements = element_index_table_[dof_indx].size();


      for(int i = 0; i < num_elements; ++i){

        //std::cout << "element index = " << element_index_table_[dof_indx][i] << ", midpoint = " << element[element_index_table_[dof_indx][i]].midpoint_ << std::endl;

        cog += element[element_index_table_[dof_indx][i]].midpoint_;

      }

      cog /= num_elements;

      return cog;

    }




  //////////////////////////////////////////////////////////////////////////////
};
}  // namespace Bembel
#endif
