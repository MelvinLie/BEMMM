// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.

#ifndef BEMBEL_ACA_ELEMENTOCTTREEMEMORY_H_
#define BEMBEL_ACA_ELEMENTOCTTREEMEMORY_H_

namespace Bembel {

/**
 *  \ingroup ClusterTree
 *  \brief this nice struct keeps track of the entire memory management of the
 *         element tree. In fact, due to this struct, the tree is now easily
 *         copyable.
 */
struct ElementOctTreeMemory {
  //////////////////////////////////////////////////////////////////////////////
  std::vector<int>::size_type nsons(ElementOctTreeNode *etn) {
    return etn->sons_.size();
  }

  ElementOctTreeNode &son(ElementOctTreeNode &etn, std::vector<int>::size_type id) {
    return (*memory_)[etn.sons_[id]];
  }

  const ElementOctTreeNode &son(const ElementOctTreeNode &etn,
                             std::vector<int>::size_type id) const {
    return (*memory_)[etn.sons_[id]];
  }

  /*
  ElementOctTreeNode &adjcent(ElementOctTreeNode &etn,
                           std::vector<int>::size_type id) {
    return (*memory_)[etn.adjcents_[id]];
  }


  const ElementOctTreeNode &adjcent(const ElementOctTreeNode &etn,
                                 std::vector<int>::size_type id) const {
    return (*memory_)[etn.adjcents_[id]];
  }
  */
  int cumNumElements(int l) const {
    //root is lvl -1
    if (l == -1) return 1;
    //Maximum number of nodes is given by: \sum_{l=-1}^{L} 8^l = (8^{L+2}-1)/7
    else return (std::pow(8,l+2)-1)/7;
  }

  ElementOctTreeNode &get_root() { return (*memory_)[0]; }

  const ElementOctTreeNode &get_root() const { return (*memory_)[0]; }

  ElementOctTreeNode &get_element(std::vector<ElementOctTreeNode>::size_type id) {
    return (*memory_)[id];
  }

  ElementOctTreeNode *get_element_ptr(std::vector<ElementOctTreeNode>::size_type id) {
    return &(*memory_)[id];
  }

  const ElementOctTreeNode &get_element(
      std::vector<ElementOctTreeNode>::size_type id) const {
    return (*memory_)[id];
  }

  //////////////////////////////////////////////////////////////////////////////
  /// iterators
  //////////////////////////////////////////////////////////////////////////////
  /*
  std::vector<ElementTreeNode>::const_iterator cluster_begin(
      const ElementTreeNode &etn) const {
    const ElementTreeNode *left = std::addressof(etn);
    while (left->sons_.size()) left = std::addressof(son(*left, 0));
    assert(left->level_ == max_level_ && "panels on different levels");
    size_t inc = cumNumElements(max_level_ - 1) + left->id_;
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cluster_end(
      const ElementTreeNode &etn) const {
    const ElementTreeNode *right = std::addressof(etn);
    while (right->sons_.size())
      right = std::addressof(son(*right, right->sons_.size() - 1));
    assert(right->level_ == max_level_ && "panels on different levels");
    size_t inc = cumNumElements(max_level_ - 1) + right->id_ + 1;
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cpbegin() const {
    size_t inc = cumNumElements(max_level_ - 1);
    return (*memory_).cbegin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator cpend() const {
    return (*memory_).cend();
  }

  std::vector<ElementTreeNode>::iterator pbegin() {
    size_t inc = cumNumElements(max_level_ - 1);
    return (*memory_).begin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator pend() {
    return (*memory_).end();
  }
  */
  std::vector<ElementOctTreeNode>::iterator lbegin(unsigned int level) {
    return (*memory_).begin() + levels_[level];
  }

  std::vector<ElementOctTreeNode>::iterator lend(unsigned int level) {
    return (*memory_).begin() + levels_[level + 1];
  }
  /*
  std::vector<ElementTreeNode>::const_iterator clbegin(
      unsigned int level) const {
    size_t inc = cumNumElements(int(level) - 1);
    return (*memory_).cbegin() + inc;
  }

  std::vector<ElementTreeNode>::const_iterator clend(unsigned int level) const {
    size_t inc = cumNumElements(int(level));
    return (*memory_).cbegin() + inc;
  }
  */
  //////////////////////////////////////////////////////////////////////////////
  /// member variables
  //////////////////////////////////////////////////////////////////////////////
  std::shared_ptr<std::vector<ElementOctTreeNode>> memory_;
  int number_of_patches_;
  int max_level_;

  //table which stores the indeces of the first element in every level of the tree
  std::vector<int> levels_;

};
}  // namespace Bembel
#endif
