// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_MEASUREMENTTREEMEMORY_H_
#define BEMBEL_MEASUREMENTTREEMEMORY_H_


namespace Bembel {

/**
 *  \ingroup MeasurementTreeMemory
 *  \brief This struct keeps track of the measurement tree memory
 */


struct MeasurementTreeMemory {


  MeasurementTreeNode &get_root() {
    return (*memory_)[0];
  }

  std::vector<int>::size_type nsons(MeasurementTreeNode *etn) {
    return etn->sons_.size();
  }

  void test(MeasurementTreeNode etn,int i){
    MeasurementTreeNode test = (*memory_)[etn.sons_[i]];
    //std::cout << "Test " << etn.sons_[i] <<std::endl;
  }

  MeasurementTreeNode &son(MeasurementTreeNode &etn,  std::vector<int>::size_type id) {//
    //std::cout << etn.sons_[id] << std::endl;
    return (*memory_)[etn.sons_[id]];
  }

  const MeasurementTreeNode &son(const MeasurementTreeNode &etn,
                             std::vector<int>::size_type id) const {
    //std::cout << "\tPos of son = " << etn.sons_[id] << std::endl;
    return (*memory_)[etn.sons_[id]];
  }

  //gives the cumulative sum of nodes of the tree up to level l
  int cumNumElements(int l) const {
    //root is lvl -1
    if (l == -1) return 1;
    //Maximum number of nodes is given by: \sum_{l=-1}^{L} 8^l = (8^{L+2}-1)/7
    else return (std::pow(8,l+2)-1)/7;
  }

  std::vector<MeasurementTreeNode>::iterator lbegin(unsigned int level) {
    return (*memory_).begin() + levels_[level];
  }

  std::vector<MeasurementTreeNode>::iterator lend(unsigned int level) {
    return (*memory_).begin() + levels_[level + 1];
  }

  MeasurementTreeNode &get_element(std::vector<MeasurementTreeNode>::size_type id) {
    return (*memory_)[id];
  }

  MeasurementTreeNode *get_element_ptr(std::vector<MeasurementTreeNode>::size_type id) {
    return &(*memory_)[id];
  }



  std::shared_ptr<std::vector<MeasurementTreeNode>> memory_;
  int max_level_;

  std::vector<int> levels_;


};

}  // namespace Bembel
#endif
