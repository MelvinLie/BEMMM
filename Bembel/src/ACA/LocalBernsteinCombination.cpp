// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LOCALBERNSTEINCOMBINATION_H_
#define BEMBEL_LOCALBERNSTEINCOMBINATION_H_

#include <Eigen/Dense>
#include <chrono>

namespace Bembel {

struct local_bernstein_combination{

  int element_id;

  std::vector<int> bernstein_basis_ids;
  std::vector<double> bernstein_weights;

};

class BSplineBasisCombination {

  public:
    BSplineBasisCombination(){};

    void append_bernstein_basis(int element_id,int bernstein_basis_id,double bernstein_weight){

      if(loc_bernstein_table_.size() == 0){

        local_bernstein_combination loc_brnst;
        loc_brnst.element_id = element_id;
        loc_brnst.bernstein_basis_ids.push_back(bernstein_basis_id);
        loc_brnst.bernstein_weights.push_back(bernstein_weight);

        loc_bernstein_table_.push_back(loc_brnst);

      }
      else{

        bool found = false;

        for(int i = 0; i < loc_bernstein_table_.size(); ++i){


          if(loc_bernstein_table_[i].element_id == element_id){
            found = true;

            loc_bernstein_table_[i].bernstein_basis_ids.push_back(bernstein_basis_id);
            loc_bernstein_table_[i].bernstein_weights.push_back(bernstein_weight);
            }
          }
          if(found == false){

            local_bernstein_combination loc_brnst;
            loc_brnst.element_id = element_id;
            loc_brnst.bernstein_basis_ids.push_back(bernstein_basis_id);
            loc_brnst.bernstein_weights.push_back(bernstein_weight);

            loc_bernstein_table_.push_back(loc_brnst);

        }
      }
    }

    void print(){
      for(int i = 0; i < loc_bernstein_table_.size(); ++i){
        std::cout << "------------------"<< std::endl;
        std::cout << "Element " <<  loc_bernstein_table_[i].element_id << std::endl;
        for(int j = 0; j <  loc_bernstein_table_[i].bernstein_basis_ids.size(); ++j){
          std::cout << "( " << loc_bernstein_table_[i].bernstein_basis_ids[j] << " , "<< loc_bernstein_table_[i].bernstein_weights[j] << " )" << std::endl;
        }
      }
    }

    const std::vector<local_bernstein_combination>::iterator begin(){
      return  loc_bernstein_table_.begin();
    }

    const std::vector<local_bernstein_combination>::iterator end(){
      return  loc_bernstein_table_.end();
    }

    std::vector<local_bernstein_combination> *get_table_ptr(){
      return &loc_bernstein_table_;
    }

  private:
    std::vector<local_bernstein_combination> loc_bernstein_table_;

};



}  // namespace Bembel
#endif
