// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_DEVICES_UNIQUELIST_H__
#define __BEMBEL_DEVICES_UNIQUELIST_H__


/**
 *  \class UniqueList
 *  \brief
 *  \todo
 **/

namespace Bembel {

class UniqueList {
 public:
   UniqueList() {}

   void append(const int row,const int col, const double val){

     bool found = false;
     for(int i = 0 ; i < size_; ++i){
       if((rows_[i] == row) && (cols_[i] == col)){
         vals_[i] += val;
         found = true;
         break;
       }
     }
     if(found == false){
       rows_.push_back(row);
       cols_.push_back(col);
       vals_.push_back(val);
       size_ += 1;
     }
   }

   int size(){
     return size_;
   }

   int get_row(const int i){
     return rows_[i];
   }

   int get_col(const int i){
     return cols_[i];
   }

   double get_val(const int i){
     return vals_[i];
   }

 private:
   int size_ = 0;
   std::vector<int> rows_;
   std::vector<int> cols_;
   std::vector<double> vals_;


};


}  // namespace Bembel

#endif
