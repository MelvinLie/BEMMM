// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_BASEFCNTREEMEMORY_H_
#define BEMBEL_BASEFCNTREEMEMORY_H_


namespace Bembel {

/**
 *  \ingroup BaseFcnTreeMemory
 *  \brief This struct keeps track of the base function tree memory
 */


struct BaseFcnTreeMemory {


  BaseFcnTreeNode &get_root() {
    return (*memory_)[0];
  }

  std::vector<int>::size_type nsons(BaseFcnTreeNode *etn) {
    return etn->sons_.size();
  }

  void test(BaseFcnTreeNode etn,int i){
    BaseFcnTreeNode test = (*memory_)[etn.sons_[i]];
    std::cout << "Test " << etn.sons_[i] <<std::endl;
  }

  BaseFcnTreeNode &son(BaseFcnTreeNode &etn,  std::vector<int>::size_type id) {//
    //std::cout << etn.sons_[id] << std::endl;
    return (*memory_)[etn.sons_[id]];
  }

  const BaseFcnTreeNode &son(const BaseFcnTreeNode &etn,
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

  std::shared_ptr<std::vector<BaseFcnTreeNode>> memory_;
  int max_level_;


};

}  // namespace Bembel
#endif
