// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_MEASUREMENTDATA_H_
#define BEMBEL_MEASUREMENTDATA_H_

#include <Eigen/Dense>


namespace Bembel {

/**
 *  \ingroup MeasurementData
 *  \brief Helper class that is used in order to input and store measurement data
 */

 // forward declaration of class MeasurementData
 class MeasurementData;


 class MeasurementData {

 public:
  //constructors
  MeasurementData() {}

  MeasurementData(const std::string &file_name, const char del = '\t') {

    std::ifstream file;
    int numLines = 0;
    std::string current_line;
    std::string current_element;

    std::vector<std::vector<double>> tmp_data;

    file.open(file_name);

    if (!file) {
      std::cerr << "File " << file_name << " doesn't exist!";
      exit(1);
    }

    // first row is irrelevant
    getline(file, current_line);


    while ( std::getline(file, current_line) ){

      std::stringstream current_data(current_line);
      std::vector<double> tmp_vector;

      while(std::getline(current_data,current_element,del)){
        tmp_vector.push_back(atof(current_element.c_str()));

      }
      tmp_data.push_back(tmp_vector);



    }

    file.close();

    int num_moves = tmp_data[0].size()/6;
    int num_samples = tmp_data.size();
    num_meas_ = num_moves*num_samples;


    pos_.resize(num_meas_,3);
    Y_.resize(num_meas_,3);

    for (int m = 0; m < num_moves; ++m){
      for (int s = 0; s < num_samples; ++s){

        pos_(m*num_samples+s,0) = tmp_data[s][m*6];
        pos_(m*num_samples+s,1) = tmp_data[s][m*6+1];
        pos_(m*num_samples+s,2) = tmp_data[s][m*6+2];

        Y_(m*num_samples+s,0) = tmp_data[s][m*6+3];
        Y_(m*num_samples+s,1) = tmp_data[s][m*6+4];
        Y_(m*num_samples+s,2) = tmp_data[s][m*6+5];

      }
    }


  }

void set_data(const Eigen::MatrixXd &pos,const Eigen::MatrixXd &measurements){

  pos_.resize(pos_.rows(),pos_.cols());
  Y_.resize(Y_.rows(),Y_.cols());

  pos_ = pos;
  Y_ = measurements;

  num_meas_ = pos_.rows();

}

void set_data(const Eigen::MatrixXd &pos,const Eigen::MatrixXd &measurements, const Eigen::VectorXd &move_nums ){

  pos_.resize(pos_.rows(),pos_.cols());
  Y_.resize(Y_.rows(),Y_.cols());

  pos_ = pos;
  Y_ = measurements;

  num_meas_ = pos_.rows();

  //the separator stores the indices of begin and end of a move
  move_separator_.clear();

  move_separator_.push_back(0);

  for(int i = 1; i < move_nums.size() ; ++i){
    std::cout << "move num = " << move_nums(i) << std::endl;
    if(move_nums(i) >  move_separator_.back()){
      move_separator_.push_back(i);
    }
  }
  move_separator_.push_back(move_nums.size());
}

void set_data(const Eigen::MatrixXd &pos, const Eigen::MatrixXd &measurements, const Eigen::VectorXd &move_nums , const Eigen::VectorXd &move_ids ){

  pos_.resize(pos_.rows(),pos_.cols());
  Y_.resize(Y_.rows(),Y_.cols());

  pos_ = pos;
  Y_ = measurements;

  num_meas_ = pos_.rows();

  //the separator stores the indices of begin and end of a move
  move_separator_.clear();

  move_ids_.clear();
  move_ids_.push_back(move_ids(0));

  move_separator_.push_back(0);

  int move_last = move_nums(0);

  for(int i = 1; i < move_nums.size() ; ++i){

    if(move_nums(i) >  move_last){

      move_separator_.push_back(i);
      move_ids_.push_back(move_ids(i));
      move_last = move_nums(i);

    }
  }
  move_separator_.push_back(move_nums.size());
}

//read measurements from a formatted data file
void set_data(const std::string &filename){

  //read measurement file
  Eigen::MatrixXd in_data = read_csv(filename,true);
  //number of measurements
  num_meas_ = in_data.rows();

  //make space for positions and data
  pos_.resize(num_meas_,3);
  Y_.resize(num_meas_,3);

  //positions are stored in the first 3 columns
  pos_ = in_data.block(0,0,num_meas_,3);
  //data is stored in column 3 to 6
  Y_ = in_data.block(0,3,num_meas_,3);
  //the separator stores the indices of begin and end of a move
  move_separator_.clear();
  //the move id stores an idetifier for each move, in this way, velocities,
  //move directions can be distriguished. It is stored in column 8
  move_ids_.clear();
  move_ids_.push_back(in_data(0,7));

  move_separator_.push_back(0);

  int move_cnt = in_data(0,6);

  for(int i = 1; i < num_meas_ ; ++i){

    if(in_data(i,6) > move_cnt){
      move_separator_.push_back(i);
      move_ids_.push_back(in_data(i,7));

      move_cnt = in_data(i,6);
    }
  }
  move_separator_.push_back(num_meas_);
}

void read_data_in_rows(const std::string &file_name, const char del = ','){

  std::ifstream file;
  int numLines = 0;
  std::string current_line;
  std::string current_element;

  std::vector<std::vector<double>> tmp_data;

  file.open(file_name);

  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!";
    exit(1);
  }

  // first row is irrelevant
  getline(file, current_line);


  while ( std::getline(file, current_line) ){

    std::stringstream current_data(current_line);
    std::vector<double> tmp_vector;

    while(std::getline(current_data,current_element,del)){
      tmp_vector.push_back(atof(current_element.c_str()));

    }
    tmp_data.push_back(tmp_vector);

  }
  file.close();

  num_meas_ = tmp_data.size();

  pos_.resize(num_meas_,3);
  Y_.resize(num_meas_,3);

  for (int m = 0; m < num_meas_; ++m){

      pos_(m,0) = tmp_data[m][0];
      pos_(m,1) = tmp_data[m][1];
      pos_(m,2) = tmp_data[m][2];

      Y_(m,0) = tmp_data[m][3];
      Y_(m,1) = tmp_data[m][4];
      Y_(m,2) = tmp_data[m][5];

  }

}

void read_run(const std::string &file_name,const int run, const char del = '\t'){

  std::ifstream file;
  int numLines = 0;
  std::string current_line;
  std::string current_element;

  std::vector<std::vector<double>> tmp_data;

  file.open(file_name);

  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!" << std::endl;
    exit(1);
  }

  //jump over all these runs
  for(int i = 0; i < 6*(run-1) ; ++i){
    getline(file, current_line);
  }
  //jump read this run
  for(int i = 0; i < 6 ; ++i){
    getline(file, current_line);

    std::stringstream current_data(current_line);
    std::vector<double> tmp_vector;

    //fist element is irrelevant
    std::getline(current_data,current_element,del);

    while(std::getline(current_data,current_element,del)){
      tmp_vector.push_back(atof(current_element.c_str()));

    }
    tmp_data.push_back(tmp_vector);

  }

  file.close();

  int num_moves = 1;
  int num_samples = tmp_data[0].size()-1;
  num_meas_ = num_samples;


  pos_.resize(num_samples,3);
  Y_.resize(num_samples,3);

  for (int s = 0; s < num_samples; ++s){

      pos_(s,0) = tmp_data[0][s];
      pos_(s,1) = tmp_data[1][s];
      pos_(s,2) = tmp_data[2][s];

      Y_(s,0) = tmp_data[3][s];
      Y_(s,1) = tmp_data[4][s];
      Y_(s,2) = tmp_data[5][s];

  }



}

int number_of_runs_raw_data(const std::string &file_name, const char del = '\t'){

  std::ifstream file;
  int numLines = 0;
  std::string current_line;
  std::string current_element;


  file.open(file_name);

  if (!file) {
    std::cerr << "File " << file_name << " doesn't exist!" << std::endl;
    exit(1);
  }


  while ( std::getline(file, current_line) ){


    ++numLines;


  }

  file.close();

  return numLines/6;
}

void rotate_measurements(double Alpha,double Beta,double Gamma){

  Eigen::Matrix3d Rot_X;

  Rot_X  <<   1,         0      ,       0        ,
              0,    std::cos(Alpha), -1.*std::sin(Alpha),
              0,    std::sin(Alpha),     std::cos(Alpha);

  Eigen::Matrix3d Rot_Y;

  Rot_Y <<   std::cos(Beta),    0    ,     std::sin(Beta),
                    0      ,    1    ,        0          ,
          -1*std::sin(Beta),    0    ,     std::cos(Beta);

  Eigen::Matrix3d Rot_Z;

  Rot_Z <<   std::cos(Gamma),-1.*std::sin(Gamma) ,        0 ,
             std::sin(Gamma),    std::cos(Gamma) ,        0 ,
                    0       ,           0        ,        1 ;

  Eigen::MatrixXd pos_rot = (Rot_X*Rot_Y*Rot_Z)*pos_.transpose().eval();

  pos_ = pos_rot.transpose();

}

void subtract_offsets(Eigen::Vector3d offsets){

  for(int i = 0; i < Y_.rows() ; ++i){
      Y_.row(i) -= offsets;
  }


}

Eigen::MatrixXd get_positions(){
  return pos_;
}
Eigen::MatrixXd get_measurements(){
  return Y_;
}
int get_number_of_measurements(){
  return num_meas_;
}


void filter_nth(int n){

  std::vector<std::vector<double>> pos_vec;
  std::vector<std::vector<double>> y_vec;

  for (int m = 0; m < num_meas_ ; ++m){
    if(m % n == 0){
      pos_vec.push_back({pos_(m,0),pos_(m,1),pos_(m,2)});
      y_vec.push_back({Y_(m,0),Y_(m,1),Y_(m,2)});
    }
  }

  num_meas_ = pos_vec.size();

  pos_.resize(num_meas_,3);
  Y_.resize(num_meas_,3);

  for (int m = 0; m < num_meas_ ; ++m){

    pos_(m,0) = pos_vec[m][0];
    pos_(m,1) = pos_vec[m][1];
    pos_(m,2) = pos_vec[m][2];

    Y_(m,0) = y_vec[m][0];
    Y_(m,1) = y_vec[m][1];
    Y_(m,2) = y_vec[m][2];

  }

}

void filter_sphere(double radius,Eigen::Vector3d center){

  std::vector<std::vector<double>> pos_vec;
  std::vector<std::vector<double>> y_vec;

  double r_tmp;


  for (int m = 0; m < num_meas_ ; ++m){

    r_tmp = std::sqrt((pos_(m,0)-center(0))*(pos_(m,0)-center(0))
                      + (pos_(m,1)-center(1))*(pos_(m,1)-center(1))
                      + (pos_(m,2)-center(2))*(pos_(m,2)-center(2)));


    if (r_tmp < radius){


      pos_vec.push_back({pos_(m,0),pos_(m,1),pos_(m,2)});
      y_vec.push_back({Y_(m,0),Y_(m,1),Y_(m,2)});

    }

  }
  num_meas_ = pos_vec.size();

  pos_.resize(num_meas_,3);
  Y_.resize(num_meas_,3);

  for (int m = 0; m < num_meas_ ; ++m){

    pos_(m,0) = pos_vec[m][0];
    pos_(m,1) = pos_vec[m][1];
    pos_(m,2) = pos_vec[m][2];

    Y_(m,0) = y_vec[m][0];
    Y_(m,1) = y_vec[m][1];
    Y_(m,2) = y_vec[m][2];

  }

}

void custom_filter(const std::function<bool(Eigen::Vector3d)> &function){

  std::vector<std::vector<double>> pos_vec;
  std::vector<std::vector<double>> y_vec;

  for (int m = 0; m < num_meas_ ; ++m){

    if (function(pos_.row(m))==true){

      pos_vec.push_back({pos_(m,0),pos_(m,1),pos_(m,2)});
      y_vec.push_back({Y_(m,0),Y_(m,1),Y_(m,2)});

    }

  }
  num_meas_ = pos_vec.size();

  pos_.resize(num_meas_,3);
  Y_.resize(num_meas_,3);

  for (int m = 0; m < num_meas_ ; ++m){

    pos_(m,0) = pos_vec[m][0];
    pos_(m,1) = pos_vec[m][1];
    pos_(m,2) = pos_vec[m][2];

    Y_(m,0) = y_vec[m][0];
    Y_(m,1) = y_vec[m][1];
    Y_(m,2) = y_vec[m][2];

  }

}


Eigen::Vector3d get_positions_cog(){

  return pos_.colwise().mean();

}

Eigen::MatrixXd *get_position_ptr(){

  return &pos_;

}

void translate_positions(Eigen::Vector3d t){

  pos_.col(0) = pos_.col(0).array() + t(0);
  pos_.col(1) = pos_.col(1).array() + t(1);
  pos_.col(2) = pos_.col(2).array() + t(2);

}

void scale_positions(double ratio){

  pos_ *= ratio;

}

Eigen::VectorXd flatten_measurements(){

  Eigen::VectorXd ret_val;

  ret_val.resize(3*num_meas_);
  ret_val.segment(0,num_meas_)          = Y_.col(0);
  ret_val.segment(num_meas_  ,num_meas_) = Y_.col(1);
  ret_val.segment(2*num_meas_,num_meas_) = Y_.col(2);

  return ret_val;

}

Eigen::MatrixXd extract_validation_set(const int n){

    //filtered data
    std::vector<std::vector<double>> pos_vec;
    std::vector<std::vector<double>> y_vec;

    //validation data
    std::vector<std::vector<double>> pos_val;
    std::vector<std::vector<double>> y_val;

    for (int m = 0; m < num_meas_ ; ++m){
      if(m % n != 0){
        pos_vec.push_back({pos_(m,0),pos_(m,1),pos_(m,2)});
        y_vec.push_back({Y_(m,0),Y_(m,1),Y_(m,2)});
      }
      else{
        pos_val.push_back({pos_(m,0),pos_(m,1),pos_(m,2)});
        y_val.push_back({Y_(m,0),Y_(m,1),Y_(m,2)});
      }
    }

    //number of remaining measurements
    num_meas_ = pos_vec.size();

    //size of validation set
    int num_val = pos_val.size();

    //remaining measurement data
    pos_.resize(num_meas_,3);
    Y_.resize(num_meas_,3);

    for (int m = 0; m < num_meas_ ; ++m){

      pos_(m,0) = pos_vec[m][0];
      pos_(m,1) = pos_vec[m][1];
      pos_(m,2) = pos_vec[m][2];

      Y_(m,0) = y_vec[m][0];
      Y_(m,1) = y_vec[m][1];
      Y_(m,2) = y_vec[m][2];

    }

    //validation data
    Eigen::MatrixXd val_data(num_val,6);
    val_data.setZero();

    for (int m = 0; m < num_val ; ++m){

      val_data(m,0) = pos_val[m][0];
      val_data(m,1) = pos_val[m][1];
      val_data(m,2) = pos_val[m][2];

      val_data(m,3) = y_val[m][0];
      val_data(m,4) = y_val[m][1];
      val_data(m,5) = y_val[m][2];

    }

    return val_data;

}

Eigen::MatrixXd get_bounding_box(){

  //bounding box is a 2 x 3 matrix with:
  //[ x0 , y0 , z0 ]
  //[ x1 , y1 , z1 ]
  if (bounding_box_.size() == 0) compute_bounding_box();
  return bounding_box_;

}

std::vector<std::vector<int>> bisect_index_list(Eigen::MatrixXd box, std::vector<int> index_list,std::vector<Eigen::MatrixXd> *bisection){

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

            if((x_sect[i] <= pos_(current_index,0)) & (pos_(current_index,0) < x_sect[i+1]) &
               (y_sect[j] <= pos_(current_index,1)) & (pos_(current_index,1) < y_sect[j+1]) &
               (z_sect[k] <= pos_(current_index,2)) & (pos_(current_index,2) < z_sect[k+1])){
                 //std::cout << "\t\t\t\ttake" << std::endl;

                 //append index to index list
                 ret_list[l].push_back(current_index);

               }
               else{
                //std::cout << "\t\t\t\treject" << std::endl;
                //add index to the remaining ones
                rem_indx_list.push_back(current_index);

               }

          }
          //hand over remaining index list to avoid double indices
          index_list = rem_indx_list;
          rem_indx_list.clear();
          l++;
      }
    }
  }
  return ret_list;
}

std::vector<std::vector<int>> binary_bisect_index_list(Eigen::MatrixXd box, std::vector<int> index_list,std::vector<Eigen::MatrixXd> *bisection){

  double bound_marg = 1e-8;

  double x_sect[3] = {box(0,0)-bound_marg,0.5*(box(1,0)+box(0,0)),box(1,0)+bound_marg};
  double y_sect[3] = {box(0,1)-bound_marg,0.5*(box(1,1)+box(0,1)),box(1,1)+bound_marg};
  double z_sect[3] = {box(0,2)-bound_marg,0.5*(box(1,2)+box(0,2)),box(1,2)+bound_marg};

  std::vector<std::vector<int>> ret_list;

  //we separate into 2 subdomains
  ret_list.resize(2);
  bisection->resize(2);
  for (int i = 0; i < 2; ++i) (*bisection)[i].resize(2,3);

  //initialize bisected bounding boxes
  (*bisection)[0](0,0) = x_sect[0];
  (*bisection)[0](1,0) = x_sect[2];
  (*bisection)[0](0,1) = y_sect[0];
  (*bisection)[0](1,1) = y_sect[2];
  (*bisection)[0](0,2) = z_sect[0];
  (*bisection)[0](1,2) = z_sect[2];

  (*bisection)[1](0,0) = x_sect[0];
  (*bisection)[1](1,0) = x_sect[2];
  (*bisection)[1](0,1) = y_sect[0];
  (*bisection)[1](1,1) = y_sect[2];
  (*bisection)[1](0,2) = z_sect[0];
  (*bisection)[1](1,2) = z_sect[2];

  //find longest side
  double dim = z_sect[2] - z_sect[0];
  int bisect_index = 2;

  if (y_sect[2] - y_sect[0] > dim) {
    bisect_index = 1;
    dim = y_sect[2] - y_sect[0];
    (*bisection)[0](1,1) = y_sect[1];
    (*bisection)[1](0,1) = y_sect[1];

  }
  else if (x_sect[2] - x_sect[0] > dim) {
    bisect_index = 0;
    dim = x_sect[2] - x_sect[0] ;
    (*bisection)[0](1,0) = x_sect[1];
    (*bisection)[1](0,0) = x_sect[1];
  }
  else{
    (*bisection)[0](1,2) = z_sect[1];
    (*bisection)[1](0,2) = z_sect[1];
  }

  //temporal index list with remaining indices
  std::vector<int> rem_indx_list;// = index_list;

  //index counter
  int current_index;


  //domain counter
  for (int l = 0; l < 2; ++l){

          for(int m = 0; m < index_list.size(); m++){
            current_index = index_list[m];

            if(((*bisection)[l](0,0) <= pos_(current_index,0)) & (pos_(current_index,0) < (*bisection)[l](1,0)) &
               ((*bisection)[l](0,1) <= pos_(current_index,1)) & (pos_(current_index,1) < (*bisection)[l](1,1)) &
               ((*bisection)[l](0,2) <= pos_(current_index,2)) & (pos_(current_index,2) < (*bisection)[l](1,2))){
                 //std::cout << "\t\t\t\ttake" << std::endl;

                 //append index to index list
                 ret_list[l].push_back(current_index);

               }
               else{
                //std::cout << "\t\t\t\treject" << std::endl;
                //add index to the remaining ones
                rem_indx_list.push_back(current_index);

               }

          }
          //hand over remaining index list to avoid double indices
          index_list = rem_indx_list;
          rem_indx_list.clear();
      }

  return ret_list;
}


void init_cluster_tree(int max_level){


  //container for domain bisection
  std::vector<Eigen::MatrixXd> domain_bisection;


  //compute the root bounding box
  if(bounding_box_set_ == false){
    compute_bounding_box();
  }


  //allocate MeasurementTreeMemory pointer
  cluster_tree_memory_ = std::make_shared<MeasurementTreeMemory>();
  //set max level
  cluster_tree_memory_->max_level_ = max_level;
  //allocate a pointer to the tree nodes
  cluster_tree_memory_->memory_ = std::make_shared<std::vector<MeasurementTreeNode>>();
  //allocate memory for the hole tree
  cluster_tree_memory_->memory_->resize(cluster_tree_memory_->cumNumElements(max_level));


  //get reference to the root node
  MeasurementTreeNode &root = cluster_tree_memory_->get_root();
  //initialize the root node measurement indices
  for (int i = 0; i < num_meas_ ; ++i) root.append_index(i);
  root.make_meas_table();
  //set memory pointer of root node
  root.set_memory(cluster_tree_memory_);
  //set root nodes center and bounding box
  root.setup_cluster_box(bounding_box_);
  //set level and position in memory of root node
  root.pos_ = 0;
  root.level_ = -1; //root is defined as level -1

  //setup level counter to -1
  lvl_ctr_ = -1;

  //setup the total node counter
  node_ctr_ = 1;

  //first level starts with root node
  cluster_tree_memory_->levels_.push_back(0);
  cluster_tree_memory_->levels_.push_back(1);


  //std::cout << "Generate Cluster Tree" << std::endl;
  generate_cluster_tree(0,1);

  cluster_tree_memory_->son(root,0);

  //sort level vector
  std::sort(cluster_tree_memory_->levels_.begin(),cluster_tree_memory_->levels_.end());
  //sort_positions();

  }

  void init_binary_tree(int max_level){

    //container for domain bisection
    std::vector<Eigen::MatrixXd> domain_bisection;

    //compute the root bounding box
    compute_bounding_box();

    //allocate MeasurementTreeMemory pointer
    binary_tree_memory_ = std::make_shared<BinaryTreeMemory>();
    //set max level
    binary_tree_memory_->max_level_ = max_level;
    //allocate a pointer to the tree nodes
    binary_tree_memory_->memory_ = std::make_shared<std::vector<BinaryTreeNode>>();
    //allocate memory for the hole tree
    binary_tree_memory_->memory_->resize(binary_tree_memory_->cumNumElements(max_level));


    //get reference to the root node
    BinaryTreeNode &root = binary_tree_memory_->get_root();
    //initialize the root node measurement index set
    for (int i = 0; i < num_meas_ ; ++i) root.append_index(i);
    //root.make_meas_table();
    //set memory pointer of root node
    root.set_memory(binary_tree_memory_);
    //set root nodes center and bounding box
    root.setup_cluster_box(bounding_box_);
    //set level and position in memory of root node
    root.pos_ = 0;
    root.level_ = -1; //root is defined as level -1

    //setup level counter to -1
    lvl_ctr_ = -1;

    //setup the total node counter
    node_ctr_ = 1;

    generate_binary_tree(0,1);

    binary_tree_memory_->son(root,0);

    //sort_positions();

    }


void print_cluster_tree(){

  MeasurementTreeNode tmp_node;

  if((*cluster_tree_memory_->memory_).size() == 0){

    std::cout << "Initialize a cluster tree first!" << std::endl;
  }
  else{
      for (int n = 0; n < node_ctr_; ++n){
        tmp_node = (*cluster_tree_memory_->memory_)[n];
        std::cout << "------------------------------" << std::endl;
        std::cout << "Node (" << n << "):" << std::endl;
        std::cout << "\tpos =" << tmp_node.pos_ << std::endl;
        std::cout << "\tCenter = ( " << tmp_node.center_(0) << " , ";
        std::cout  << tmp_node.center_(1) << " , " << tmp_node.center_(2) << " )" << std::endl;
        std::cout << "\tDiam = " << tmp_node.diam_ << std::endl;
        std::cout << "\tBbox = ( " << tmp_node.bbox_(0,0) << " , ";
        std::cout  << tmp_node.bbox_(0,1) << " , " << tmp_node.bbox_(0,2) << " )" << std::endl;
        std::cout << "\t       ( " << tmp_node.bbox_(1,0) << " , ";
        std::cout  << tmp_node.bbox_(1,2) << " , " << tmp_node.bbox_(1,2) << " )" << std::endl;
        std::cout << "\tNumber of measurements = " << tmp_node.indices_.size() << std::endl;
        std::cout << "\tMeasurement table = " << tmp_node.meas_table_.transpose() << std::endl;
        std::cout << "\tMeasurement indices = ";
        for (int i = 0; i < tmp_node.indices_.size(); ++i) std::cout << tmp_node.indices_[i] << " ";
        std::cout << std::endl;
        std::cout << "\tLevel = " << tmp_node.level_ << std::endl;
        std::cout << "\tNumber of sons = " << tmp_node.sons_.size() << std::endl;
        std::cout << "\tsons = ";
        for (int i = 0 ; i <  tmp_node.sons_.size() ; ++i) std::cout << tmp_node.sons_[i] << " , ";
        std::cout << std::endl;
      }

  }
}

void sort_positions(){

  //sorted positions
  Eigen::MatrixXd sorted_pos;
  sorted_pos.resize(pos_.rows(),3);

  //index table
  std::vector<int> index_table;
  index_table.resize(pos_.rows());

  //temporal node storage
  MeasurementTreeNode tmp_node;

  //measurement counter
  int meas_cnt = 0;

  if((*cluster_tree_memory_->memory_).size() == 0){

    std::cout << "Initialize a cluster tree first!" << std::endl;
  }
  else{
      for (int n = 0; n < node_ctr_; ++n){
        tmp_node = (*cluster_tree_memory_->memory_)[n];

        //we sort corresponding to leafs of the tree
        if(tmp_node.sons_.size() == 0){

          for(int m = 0; m < tmp_node.indices_.size(); ++m){

            //copy this row in sorted matrix
            sorted_pos.row(meas_cnt) = pos_.row(tmp_node.indices_[m]);
            //save the change in index table
            index_table[tmp_node.indices_[m]] = meas_cnt;

            meas_cnt++;
            //tmp_node.indices_[m] = m+meas_cnt;
          }
        }
      }
      pos_ = sorted_pos;

      //apply this change to all nodes
      for (int n = 0; n < node_ctr_; ++n){
        tmp_node = (*cluster_tree_memory_->memory_)[n];

        for(int m = 0; m < tmp_node.indices_.size(); ++m){
          (*cluster_tree_memory_->memory_)[n].indices_[m] = index_table[tmp_node.indices_[m]];
        }
      }

  }
}

std::shared_ptr<MeasurementTreeMemory> get_tree_memory(){

  return cluster_tree_memory_;
}

void set_bounding_box(const Eigen::Matrix<double,2,3> &bbox){
  bounding_box_set_ = true;
  bounding_box_ = bbox;
}

void set_min_numel(int min_numel){

  min_numel_ = min_numel;

}

std::vector<int> get_move_separator(){

  return move_separator_;

}

std::vector<int> get_move_ids(){

  return move_ids_;

}

Eigen::MatrixXd get_move(const int move){

  int num_positions = move_separator_[move+1] - move_separator_[move];

  Eigen::MatrixXd ret_mat(num_positions,6);

  ret_mat.block(0,0,num_positions,3) = pos_.block(move_separator_[move],0,num_positions,3);
  ret_mat.block(0,3,num_positions,3) = Y_.block(move_separator_[move],0,num_positions,3);

  return ret_mat;

}

private:
  Eigen::MatrixXd pos_;
  Eigen::MatrixXd Y_;
  Eigen::MatrixXd bounding_box_;
  std::shared_ptr<MeasurementTreeMemory> cluster_tree_memory_;
  std::shared_ptr<BinaryTreeMemory> binary_tree_memory_;

  //move separator
  std::vector<int> move_separator_;
  std::vector<int> move_ids_;

  //number of measurements
  int num_meas_;

  //to count the level in generate_cluster_tree
  int lvl_ctr_;

  //to count the current number of nodes in generate_cluster_tree
  int node_ctr_;

  //minimum number of measurements in a cell with sons
  int min_numel_ = 1;

  int min_tree_level_ = 2;

  bool bounding_box_set_ = false;


  void compute_bounding_box(){

    //bounding box is a 2 x 3 matrix with:
    //[ x0 , y0 , z0 ]
    //[ x1 , y1 , z1 ]
    bounding_box_.resize(2,3);

    bounding_box_(0,0) = pos_.col(0).minCoeff();
    bounding_box_(1,0) = pos_.col(0).maxCoeff();
    bounding_box_(0,1) = pos_.col(1).minCoeff();
    bounding_box_(1,1) = pos_.col(1).maxCoeff();
    bounding_box_(0,2) = pos_.col(2).minCoeff();
    bounding_box_(1,2) = pos_.col(2).maxCoeff();

  }


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

      if((lvl_ctr_ < min_tree_level_) || (*cluster_tree_memory_->memory_)[mem_from + n].indices_.size() > min_numel_){

        //number of measurements in this node
        //int num_node_meas = (*cluster_tree_memory_->memory_)[mem_from + n].indices_.size();

        //bisect the domain and separate the index list of this node
        bs_idx_list = bisect_index_list((*cluster_tree_memory_->memory_)[mem_from + n].bbox_,
                                          (*cluster_tree_memory_->memory_)[mem_from + n].indices_,
                                          &domain_bisection);

        for (int s = 0; s < bs_idx_list.size(); ++s){
          //it can be that there are no measurements in this subsection!
          if (bs_idx_list[s].size() == 0) continue;
          else{
            //there are some measurements here, this is a node
            //increment the total node counter
            node_ctr_++;
            //initialize the new node
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].indices_ = bs_idx_list[s];
            //initialize the measurement table (we could in future remove indices_ since it stores the same info)
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].set_meas_table(bs_idx_list[s]);
            //give it access to memory
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].set_memory(cluster_tree_memory_);
            //setup his father
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].father_ = mem_from + n;
            //this the geometric information
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].setup_cluster_box(domain_bisection[s]);
            //mark level in level indicator
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].level_ = lvl_ctr_;
            //position of this node in memory_
            (*cluster_tree_memory_->memory_)[mem_to + son_ctr].pos_ = node_ctr_-1;
            //add this position for the son to father
            (*cluster_tree_memory_->memory_)[mem_from + n].sons_.push_back(mem_to + son_ctr);
            //increment son counter
            son_ctr++;
          }

        }

      }


    }
    if (lvl_ctr_ < cluster_tree_memory_->max_level_){
      //refine further
      generate_cluster_tree(mem_to, mem_to + son_ctr);

    }
    //save a reference to the first element of next level in levels_
    cluster_tree_memory_->levels_.push_back(mem_to + son_ctr);

    }

    void generate_binary_tree(int mem_from,int mem_to){

      //increment level counter
      lvl_ctr_ += 1;
      //number of current root nodes
      int num_root = mem_to - mem_from;
      //position of first son in memory = mem_to;
      //initialize son counter
      int son_ctr = 0;
      //bisected index list
      std::vector<std::vector<int>> bs_idx_list;
      //container for domain bisection
      std::vector<Eigen::MatrixXd> domain_bisection;

      //iterate over all given nodes
      for (int n = 0; n < num_root; ++n){

        //number of measurements in this node
        //int num_node_meas = (*cluster_tree_memory_->memory_)[mem_from + n].indices_.size();

        //bisect the domain and separate the index list of this node
        bs_idx_list = binary_bisect_index_list((*binary_tree_memory_->memory_)[mem_from + n].bbox_,
                                          (*binary_tree_memory_->memory_)[mem_from + n].indices_,
                                          &domain_bisection);

        for (int s = 0; s < bs_idx_list.size(); ++s){
          //it can be that there are no measurements in this subsection!
          if (bs_idx_list[s].size() == 0) continue;
          else{
            //there are some measurements here, this is a node
            //increment the total node counter
            node_ctr_++;
            //initialize the new node
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].indices_ = bs_idx_list[s];
            //give it access to memory
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].set_memory(binary_tree_memory_);
            //setup his father
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].father_ = mem_from + n;
            //this the geometric information
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].setup_cluster_box(domain_bisection[s]);
            //mark level in level indicator
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].level_ = lvl_ctr_;
            //position of this node in memory_
            (*binary_tree_memory_->memory_)[mem_to + son_ctr].pos_ = node_ctr_-1;
            //add this position for the son to father
            (*binary_tree_memory_->memory_)[mem_from + n].sons_.push_back(mem_to + son_ctr);
            //increment son counter
            son_ctr++;
          }

          }


        }
      if (lvl_ctr_ < binary_tree_memory_->max_level_){
        //refine further

        generate_binary_tree(mem_to, mem_to + son_ctr);

      }

      }

 };


}  // namespace Bembel
#endif
