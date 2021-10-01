// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_GEOMETRY_GEOMETRY_H_
#define BEMBEL_GEOMETRY_GEOMETRY_H_

namespace Bembel {
/**
 *  \ingroup Geometry
 *  \brief this class wraps a GeometryVector and provides some basic
 *         functionality, like reading Geometry files
 */
class Geometry {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    Constructors
  //////////////////////////////////////////////////////////////////////////////
  Geometry() {}
  Geometry(const std::string &filename) { init_Geometry(filename); }
  Geometry(Geometry &&other) { geometry_ = std::move(other.geometry_); }
  // though we are using a shared pointer, we are creating an actual
  // copy here. might be useful if we want to modify the geometry object
  Geometry(const Geometry &other) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = *(other.geometry_);
  }
  Geometry(const PatchVector &in) {
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = in;
  }
  Geometry &operator=(Geometry other) {
    std::swap(geometry_, other.geometry_);
    return *this;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_Geometry
  //////////////////////////////////////////////////////////////////////////////
  inline void init_Geometry(const std::string &filename) {
    // Note that the Shredder is required. The order of ansatz functions allows
    // to be chosen higher than the smoothness of the NÙRBS mappings. Thus, we
    // need to shredder the geometry mappings to have Bezier patches. You can
    // achieve the higher regularity by changing coefficients in the projector.
    auto tmp = Bembel::PatchShredder(Bembel::LoadGeometryFile(filename));
    geometry_ = std::make_shared<PatchVector>();
    *geometry_ = tmp;
  }
  //////////////////////////////////////////////////////////////////////////////
  //    getters
  //////////////////////////////////////////////////////////////////////////////
  const PatchVector &get_geometry() const { return *geometry_; }
  PatchVector &get_geometry() { return *geometry_; }
  const std::shared_ptr<PatchVector> get_geometry_ptr() const {
    return geometry_;
  }
  std::shared_ptr<PatchVector> get_geometry_ptr() { return geometry_; }
  int get_number_of_patches() const { return geometry_->size(); };

  void scale_geometry(double ratio){

    for(int j = 0; j < (*geometry_).size(); ++j){

      //std::cout << "Patch " << j << std::endl;
      for(int k = 0; k < (*geometry_)[j].data_.size(); ++k){

          if ((k+1) % 4 != 0){

            (*geometry_)[j].data_[k] *= ratio;
          }
        }
    }
  }

  void shift_geometry(const Eigen::Vector3d shift){


    for(int j = 0; j < (*geometry_).size(); ++j){


      for(int k = 0; k < (*geometry_)[j].data_.size(); ++k){

        //std::cout << (*geometry_)[j].data_[k] << std::endl;
          if ((k+1) % 4 != 0){
            int rem = k % 4;

            (*geometry_)[j].data_[k] += shift[rem]*(*geometry_)[j].data_[k + 3 - rem];
          }
        }
    }
  }

  /*
  void scale_geometry_along_axis(const double ratio,const int axis){

    for(int j = 0; j < (*geometry_).size(); ++j){


      //std::cout << "Patch " << j << std::endl;
      //std::cout << "(*geometry_)[j].data_ :" << std::endl;
      //int i = 0;
      for(int k = 0; k < (*geometry_)[j].data_.size(); ++k){


          //std::cout << (*geometry_)[j].data_[k] << " , ";
          //std::cout << (*geometry_)[j].data_[k]  << " , ";

          //if (i == axis){//((k+1) % 4 != 0){

            if ((k+1) % 4 != 0){
              if( k % 4 == axis){
                  (*geometry_)[j].data_[k] *= ratio;
              }

            }
          //}

          //if ((k+1) % 4 == 0){
          //    i += 1;
          //    std::cout << std::endl;
          //}
        }
    }

  }
  */

  Eigen::MatrixXd get_bounding_box() const {

    double x_min = std::numeric_limits<double>::infinity();
    double y_min = std::numeric_limits<double>::infinity();
    double z_min = std::numeric_limits<double>::infinity();
    double x_max = -1*std::numeric_limits<double>::infinity();
    double y_max = -1*std::numeric_limits<double>::infinity();
    double z_max = -1*std::numeric_limits<double>::infinity();

    //corner points in 3d
    Eigen::Vector3d c1,c2,c3,c4;
    Eigen::Vector3d cntr;
    Eigen::Vector3d eval_ptr;

    int steps_eval = 10;
    double u,v;


    for(int j = 0; j < (*geometry_).size(); ++j){

      for(int i = 0; i < steps_eval;++i){

        u =  i/(steps_eval-1);

        for(int k = 0; k < steps_eval;++k){

          v =  k/(steps_eval-1);

          cntr = (*geometry_)[j].eval(u,v);

          if (cntr(0) < x_min) x_min = cntr(0);
          if (cntr(1) < y_min) y_min = cntr(1);
          if (cntr(2) < z_min) z_min = cntr(2);

          if (cntr(0) > x_max) x_max = cntr(0);
          if (cntr(1) > y_max) y_max = cntr(1);
          if (cntr(2) > z_max) z_max = cntr(2);
        }

      }


    }
    Eigen::MatrixXd bounding_box(2,3);

    bounding_box(0,0) = x_min;
    bounding_box(1,0) = x_max;
    bounding_box(0,1) = y_min;
    bounding_box(1,1) = y_max;
    bounding_box(0,2) = z_min;
    bounding_box(1,2) = z_max;

    return bounding_box;

  }

  //////////////////////////////////////////////////////////////////////////////
  //    private member variables
  //////////////////////////////////////////////////////////////////////////////
 private:
  std::shared_ptr<PatchVector> geometry_;
};
}  // namespace Bembel
#endif
