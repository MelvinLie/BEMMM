// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef __BEMBEL_LRMATRIX_LRMATRIX_H__
#define __BEMBEL_LRMATRIX_LRMATRIX_H__


/**
 *  \class LRMatrix
 *  \brief Low Rank Matrix class, which extends the EigenBase class.
 *
 *  The idea is to provide an easy to use interface to the H2-matrix
 *  from the fast boundary element method. At the moment, we inherit the
 *  traits of an Eigen::SparseMatrix, since this seems to be the minimum
 *  properties for a derived object to ensure that the matrix-vector
 *  multiplication can be specialised for H2Matrix.
 *  In particular, this allows for the use of the Eigen iterative solvers
 *  with a Hierarchical matrix.
 *
 *  \todo Maybe, we find something better then the SparsMatrix traits
 *        in the future
 **/
namespace Eigen {
/// forward definition of the LRMatrix Class in order to define traits
template <typename Device>
class LRMatrix;
/// inherit the traits from the Eigen::SparseMatrix class
namespace internal {
template <typename Device>
struct traits<LRMatrix<Device>> : public internal::traits<SparseMatrix<double>> {};
}  // namespace internal

// actual definition of the class
template <typename Device>
class LRMatrix : public EigenBase<LRMatrix<Device>> {
 public:
  //////////////////////////////////////////////////////////////////////////////
  /// Eigen related things
  //////////////////////////////////////////////////////////////////////////////
  // Required typedefs, constants and so on.
  typedef double Scalar;
  typedef typename NumTraits<double>::Real RealScalar;
  typedef Index StorageIndex;
  enum {
    ColsAtCompileTime = Dynamic,
    MaxColsAtCompileTime = Dynamic,
    IsRowMajor = false,
    Flags = NestByRefBit
  };
  // Minimum specialisation of EigenBase methods
  Index rows() const { return rows_; }
  Index cols() const { return cols_; }
  // Definition of the matrix multiplication
  template <typename Rhs>
  Product<LRMatrix, Rhs, AliasFreeProduct> operator*(
      const MatrixBase<Rhs>& x) const {
    return Product<LRMatrix, Rhs, AliasFreeProduct>(*this, x.derived());
  }
  //////////////////////////////////////////////////////////////////////////////
  /// constructors
  //////////////////////////////////////////////////////////////////////////////
  LRMatrix() {}
  /**
   * \brief Assemble H2-Matrix for linear operator linOp and AnsatzSpace
   * ansatz_space with number_of_points interpolation points in one direction of
   * the unit square (standard is number_of_points=16)
   */
  //Make this a template class to be able to handle different devices. So far we
  //are using the Hall probe array
  //template <typename Derived>
  //<typename Pot,typename LinOp>
  void init_LRMatrix(Device *device) {//HallProbeArray *device

    cols_ = device->get_number_of_dofs();
    rows_ = device->get_number_of_measurements()*device->get_number_of_sensors();

    device_ = device;

    stability_vec = device->get_stability_vector();

    if(device->gauged_formulation_enabled()){
      set_gauging(true);
    }

  }

  void setup_inverse_problem(const Eigen::SparseMatrix<double> &R_inv ){

    meas_prec_ = R_inv;

    squared_ = true;

    if(gauging_){

      cols_ = device_->get_number_of_dofs() - 1;
      rows_ = device_->get_number_of_dofs() - 1;
    }
    else{

      cols_ = device_->get_number_of_dofs();
      rows_ = device_->get_number_of_dofs();
    }


  }



  void setup_inverse_problem(const Eigen::SparseMatrix<double> &R_inv,
                              const Eigen::SparseMatrix<double> &P ){
    meas_prec_ = R_inv;
    prior_prec_ = P;
    squared_ = true;
    prior_ = true;

  }

  Eigen::VectorXd make_rhs(const Eigen::VectorXd &rhs){

    return device_->mat_vec_transpose(rhs);

  }

  //////////////////////////////////////////////////////////////////////////////
  /// getter
  //////////////////////////////////////////////////////////////////////////////
  /*
  Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> get_dense() const {
    Eigen::Matrix<ScalarT, Eigen::Dynamic, Eigen::Dynamic> dense(rows(),
                                                                 cols());
    for (int i = 0; i < cols(); ++i) {
      Eigen::VectorXd unit(cols());
      unit.setZero();
      unit(i) = 1.;
      dense.col(i) = (*this) * unit;
    }
    return dense;
  }
  */
  Eigen::MatrixXd get_dense() const{

    return test_matrix_;

  }

  Device *get_device() const{

    return device_;

  }

  Eigen::VectorXd squared_product(const Eigen::VectorXd in) const{

    //std::cout << "in = " << in.rows() << " x " << in.cols() << std::endl;

    Eigen::VectorXd tmp = meas_prec_*device_->mat_vec(in);

    //std::cout << "tmp = " << tmp.rows() << " x " << tmp.cols() << std::endl;

    Eigen::VectorXd ret_val = device_->mat_vec_transpose(tmp);

    //if prior_
    // just add the contribution from the prior
    // not implemeted yet

    return ret_val;

  }

  bool is_squared() const {

    return squared_;
  }

  bool is_transposed() const {

    return transposed_;
  }

  void set_test(){

    test_matrix_ << 1.,2.,3.,
                    2.,2.,2.,
                    3.,2.,1.;

    test_vec_ << 1,
                  1,
                  1;

    test_ = true;
    squared_ = true;

    rows_ = 3;
    cols_ = 3;

  }

  void set_squared(bool enable){

    squared_ = enable;
    if(enable == false){

      transposed_ = false;

      cols_ = device_->get_number_of_dofs();
      rows_ = device_->get_number_of_measurements()*device_->get_number_of_sensors();


    }

  }


  void set_transposed(){

    transposed_ = true;

    rows_ = device_->get_number_of_dofs();
    cols_ = device_->get_number_of_measurements()*device_->get_number_of_sensors();

  }



  //////////////////////////////////////////////////////////////////////////////
  /// private members
  //////////////////////////////////////////////////////////////////////////////
 private:
  // we declare functionality which has not been implemented (yet)
  // to be private
  LRMatrix(const LRMatrix& LR);
  LRMatrix(LRMatrix&& LR);
  LRMatrix& operator=(const LRMatrix& LR);
  LRMatrix& operator=(LRMatrix&& LR);


  Device *device_;

  Index rows_;
  Index cols_;

  bool squared_ = false;
  bool transposed_ = false;
  bool prior_ = false;
  bool gauging_ = true;

  bool test_ = false;

  Eigen::Matrix3d test_matrix_;
  Eigen::Vector3d test_vec_;

  Eigen::SparseMatrix<double> prior_prec_;
  Eigen::SparseMatrix<double> meas_prec_;

  Eigen::VectorXd stability_vec;


  void set_gauging(bool enable){

    gauging_ = enable;

    if(enable){
      cols_ = device_->get_number_of_dofs() - 1;

      if(squared_){
        rows_ = device_->get_number_of_dofs() - 1;
      }
    }
  }

};  // namespace Eigen

/**
 * \brief Implementation of H2Matrix * Eigen::DenseVector through a
 * specialization of internal::generic_product_impl
 */
namespace internal {
template <typename Rhs, typename Device >
struct generic_product_impl<LRMatrix<Device>, Rhs, SparseShape, DenseShape,
                            GemvProduct>  // GEMV stands for matrix-vector
    : generic_product_impl_base<LRMatrix<Device>, Rhs,
                                generic_product_impl<LRMatrix<Device>, Rhs>> {
  template <typename Dest>
  static void scaleAndAddTo(Dest& dst, const LRMatrix<Device>& lhs,
                            const Rhs& rhs, const double& alpha) {

      assert(alpha == double(1) && "scaling is not implemented");
      EIGEN_ONLY_USED_FOR_DEBUG(alpha);

      if(lhs.is_squared()){

        dst += lhs.squared_product(rhs);

      }
      else if(lhs.is_transposed()){

        dst += lhs.get_device()->mat_vec_transpose(rhs);
      }
      else{

        //std::cout << "!!!!!!!!!!!!!!" << std::endl;
        dst += lhs.get_device()->mat_vec(rhs);

      }




  }
};

}  // namespace internal
}  // namespace Eigen

#endif
