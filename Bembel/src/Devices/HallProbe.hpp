// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_DEVICES_HALLPROBE_H_
#define BEMBEL_DEVICES_HALLPROBE_H_

namespace Bembel {
// forward declaration of class HallProbe in order to define
// traits
template <typename Derived,typename LinOp>
class HallProbe;

template <typename Derived,typename LinOp>
struct SensorTraits<HallProbe<Derived,LinOp>> {
  static constexpr int OutputSpaceDimension = 1;
};

/**
 * \ingroup Devices
 */
template <typename Derived,typename LinOp>
class HallProbe {


 public:
  HallProbe() {
    n_(0) = 1;
    n_(1) = 0;
    n_(2) = 0;

    tf_.resize(2);
    tf_.push_back(0.);
    tf_.push_back(1.);

    //by default we assume that the mapper arm is orientated along z
    setup_transversal_vectors(2);

  }
  HallProbe(Eigen::Vector3d n, std::vector<double> tf, const int arm_axis_id = 2) {

    n_ = n;
    tf_ = tf;

    setup_transversal_vectors(arm_axis_id);


  }

  void set_parameters(Eigen::Vector3d n, std::vector<double> tf, const int arm_axis_id = 2){

    n_ = n;
    tf_ = tf;

    setup_transversal_vectors(arm_axis_id);

    return;

  }

  void get_parameters(Eigen::Vector3d *n, std::vector<double> *tf){

    n = &n_;
    tf = &tf_;

    return;

  }

  Eigen::Vector3d get_orientation(){

    return n_;
  }

  std::vector<double> get_transfer_function(){

    return tf_;
  }

  Eigen::VectorXd apply_transfer_function(const Eigen::VectorXd &in){

    Eigen::VectorXd ret_val, power;
    ret_val.resize(in.size());
    ret_val.setOnes();
    ret_val *= tf_[0];

    power.resize(in.size());
    power.setOnes();

    for(int i = 1 ; i < tf_.size() ; ++i){
      power = power.cwiseProduct(in).eval();

      ret_val += tf_[i]*power;

    }

    return ret_val;
  }


  double evaluateIntegrand_impl(const FunctionEvaluator<LinOp> &fun_ev,
                         const ElementTreeNode &element,
                         const Eigen::Vector3d &point,
                         const SurfacePoint &p) const {

    // evaluate gradient
    Eigen::Matrix<double,1, 3> gradient = disc_pot_.evaluateIntegrandDerivative_impl(fun_ev, element, point, p);

    double integrand = tf_[0];

    double nB = gradient(0,0)*n_(0) + gradient(0,1)*n_(1) +   gradient(0,2)*n_(2);

    for(int i = 1 ; i < tf_.size() ; ++i) integrand += tf_[i]*std::pow(nB,i);

    return integrand;
  }

  /*
  Evaluates the matrix entries for the evaluation of n.B. The nonlinearity of
  must be computed on n.(B_mat v)
  */
  void evaluateIntegrand(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::MatrixXd intvec;
    intvec.resize(3,I2);
    intvec.setZero();

    //compute gradient entries
    //disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);
    disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);

    //Eigen::MatrixXd nB;
    //nB.resize(1,I2);
    //nB.setZero();

    (*intval) += intvec.row(0)*n_(0) + intvec.row(1)*n_(1) + intvec.row(2)*n_(2);


    //for(int i = 1 ; i < tf_.size() ; ++i) (*intval) += std::pow(tf_[i],i)*nB;

  }


  /*
  Evaluates the matrix entries for the evaluation of n.B. The nonlinearity of
  must be computed on n.(B_mat v)
  Here, we integrate over the double layer potential. Take care to increase the
  quadrature degree when using this function.
  */
  void evaluateIntegrand_far_field(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::Matrix<double,-1,3> intvec;
    intvec.resize(I2,3);
    intvec.setZero();

    //compute gradient entries
    //disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);
    disc_pot_.evaluateIntegrandDerivative_far_field(super_space,element,point,p,&intvec);

    //Eigen::MatrixXd nB;
    //nB.resize(1,I2);
    //nB.setZero();

    (*intval) -= intvec.col(0)*n_(0) + intvec.col(1)*n_(1) + intvec.col(2)*n_(2);


    //for(int i = 1 ; i < tf_.size() ; ++i) (*intval) += std::pow(tf_[i],i)*nB;

  }


  /*
  Evaluates the matrix entries for the evaluation of n.B as well as t1.B  and
  t2.B. t1 and t2 are tangential vectors. These expressions are needed for UQ.
  The nonlinearity of must be computed on n.(B_mat v)
  */
  void evaluateIntegrands_UQ(const SuperSpace<LinOp> &super_space,
                          const ElementTreeNode &element,
                          const Eigen::Vector3d &point,
                          const SurfacePoint &p,
                          Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::MatrixXd intvec;
    intvec.resize(3,I2);
    intvec.setZero();

    //compute gradient entries
    disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);

    //Eigen::MatrixXd nB;
    //nB.resize(1,I2);
    //nB.setZero();


    (*intval).row(0) += intvec.row(0)*n_(0)  + intvec.row(1)*n_(1)  + intvec.row(2)*n_(2);
    (*intval).row(1) += intvec.row(0)*t1_(0) + intvec.row(1)*t1_(1) + intvec.row(2)*t1_(2);
    (*intval).row(2) += intvec.row(0)*t2_(0) + intvec.row(1)*t2_(1) + intvec.row(2)*t2_(2);


    //for(int i = 1 ; i < tf_.size() ; ++i) (*intval) += std::pow(tf_[i],i)*nB;

  }



  /*
  Evaluates the matrix entries for the evaluation of n.B. The nonlinearity of
  must be computed on n.(B_mat v), based on a single basis function
  */
  void evaluateIntegrand_single_basis(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         double *intval, int k) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::MatrixXd intvec;
    intvec.resize(3,I2);
    intvec.setZero();

    //compute gradient entries
    //disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);
    disc_pot_.evaluateIntegrandDerivative_single_basis(super_space,element,point,p,&intvec,k);

    //Eigen::MatrixXd nB;
    //nB.resize(1,I2);
    //nB.setZero();

    (*intval) += (intvec.row(0)*n_(0) + intvec.row(1)*n_(1) + intvec.row(2)*n_(2))(0);


    //for(int i = 1 ; i < tf_.size() ; ++i) (*intval) += std::pow(tf_[i],i)*nB;

  }

  /*
  Evaluates the matrix entries for the evaluation of n.B. The nonlinearity of
  must be computed on n.(B_mat v). This function is used to integrate the
  BSpline as a linear combination of Bernstein Polynomials directly
  */
  void evaluateIntegrand_BSpline(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         double *intval, double *weights, int *brnst_fcns, int num_fcns) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::MatrixXd intvec;
    intvec.resize(3,1);
    intvec.setZero();

    //compute gradient entries
    //disc_pot_.evaluateIntegrandDerivative_mat(super_space,element,point,p,&intvec);
    disc_pot_.evaluateIntegrandDerivative_BSpline(super_space,element,point,p,&intvec,weights,brnst_fcns,num_fcns);

    //std::cout << "Intvec " << intvec << std::endl;

    (*intval) += (intvec.row(0)*n_(0) + intvec.row(1)*n_(1) + intvec.row(2)*n_(2))(0);



  }

  /*
  Evaluates the matrix entries for the evaluation of n.B. The nonlinearity of
  must be computed on n.(B_mat v)
  */
  void evaluateIntegrandDerivative(const SuperSpace<LinOp> &super_space,const ElementTreeNode &element,
                         const Eigen::Vector3d &point,const SurfacePoint &p,
                         Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> *intval) const {

    //polynomial degree
    int P = super_space.get_polynomial_degree();
    int I2 = (P+1)*(P+1);

    Eigen::MatrixXd intmat;
    intmat.resize(9,I2);
    intmat.setZero();

    //compute gradient entries
    disc_pot_.evaluateIntegrandSecondDerivative_mat(super_space,element,point,p,&intmat);

    //Eigen::MatrixXd nB;
    //nB.resize(1,I2);
    //nB.setZero();

    (*intval).row(0) += intmat.row(0)*n_(0) + intmat.row(1)*n_(1) + intmat.row(2)*n_(2);
    (*intval).row(1) += intmat.row(3)*n_(0) + intmat.row(4)*n_(1) + intmat.row(5)*n_(2);
    (*intval).row(2) += intmat.row(6)*n_(0) + intmat.row(7)*n_(1) + intmat.row(8)*n_(2);

    //for(int i = 1 ; i < tf_.size() ; ++i) (*intval) += std::pow(tf_[i],i)*nB;

  }

  Derived &get_potential(){

    return disc_pot_;
  }

  Derived *get_potential_ptr(){

    return &disc_pot_;
  }

  void rotate_probe(double Alpha,double Beta,double Gamma){

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

        n_ = Rot_X*Rot_Y*Rot_Z*n_;
        n_ = Rot_X*Rot_Y*Rot_Z*n_;
        n_ = Rot_X*Rot_Y*Rot_Z*n_;

  }
  void reset_offset(){

       tf_[0] = 0.;

  }

  double get_sensitivity(){

    return tf_[1];

  }

  void setup_transversal_vectors(const int arm_axis_id){

    /*
    Eigen::Matrix3d Rot_1, Rot_2;

    //helper matrices for unit rotations
    Eigen::Matrix3d Rot_x, Rot_y, Rot_z;

    //rot 90 deg around x
    Rot_x << 1., 0.,  0.,
              0., 0., -1.,
              0., 1.,  0. ;

    //rot 90 deg around y
    Rot_y << 0., 0.,  1.,
            0., 1.,  0.,
            -1., 0.,  0. ;

    //rot 90 deg around z
    Rot_z << 0., -1.,  0.,
            1.,  0.,  0.,
            0.,  0.,  1. ;
    */
    // arm is orientated along x
    if (arm_axis_id == 0){
      //std::cout << "X" << std::endl;
      //theta 1 is in xy plane
      //90 deg around z
      //Rot_1 = Rot_z;
      t1_ << -n_(1),
              n_(0),
              0.;

      //theta 2 is in xz plane
      //90 deg around y
      //Rot_2 = Rot_y;
      t2_ <<  n_(2),
              0.,
              -n_(0);

    }
    // arm is orientated along y
    if (arm_axis_id == 1){
      //std::cout << "Y" << std::endl;
      //theta 1 is in yz plane
      //90 deg around x
      //Rot_1 = Rot_x;
      t1_ <<  0.,
              -n_(2),
               n_(1);


      //theta 2 is in xy plane
      //90 deg around z
      //Rot_2 = Rot_z;
      t2_ <<  -n_(1),
              n_(0),
              0.;

    }
    // arm is orientated along z
    if (arm_axis_id == 2){
      //std::cout << "Z" << std::endl;
      //theta 1 is in yz plane
      //90 deg around x
      //Rot_1 = Rot_x;
      t1_ <<  0.,
              -n_(2),
               n_(1);

      //theta 2 is in xz plane
      //90 deg around y
      //Rot_2 = Rot_y;
      t2_ <<  n_(2),
              0.,
              -n_(0);

    }

    //t1_ = Rot_1 * n_;
    //t2_ = Rot_2 * n_;

    //std::cout << "t_1 = " << t1_.transpose() << std::endl;
    //std::cout << "t_2 = " << t2_.transpose() << std::endl;

  }

  Eigen::Vector3d get_orientation_vector(){

    return n_;

  }


 private:
   //Orientation
   Eigen::Vector3d n_;
   //Rotations of n_ around the two transversal axes of a mapper arm
   Eigen::Vector3d t1_,t2_;
   //Transfer function
   std::vector<double> tf_;
   //discrete Potential
   Derived disc_pot_;



};

}  // namespace Bembel
#endif
