// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_TRACKING_DIRECTRUNGEKUTTA_H_
#define BEMBEL_TRACKING_DIRECTRUNGEKUTTA_H_

namespace Bembel {
/**
 *  \ingroup Potential
 *  \brief DiscretePotential
 *  \todo  Add a documentation
 */
template <typename LinOp>
class DirectRungeKutta {
 public:
  //////////////////////////////////////////////////////////////////////////////
  //    constructors
  //////////////////////////////////////////////////////////////////////////////
  DirectRungeKutta() {}
  DirectRungeKutta(const AnsatzSpace<LinOp> &ansatz_space) {

    init_Tracker(ansatz_space);

  }
  //////////////////////////////////////////////////////////////////////////////
  //    init_DiscretePotential
  //////////////////////////////////////////////////////////////////////////////
  void init_Tracker(const AnsatzSpace<LinOp> &ansatz_space) {
    ansatz_space_ = ansatz_space;
    fun_ev_ = FunctionEvaluator<LinOp>(ansatz_space_);
    /**
     * \todo obtain this from ansatz space
     */
    deg_ = ansatz_space_.get_polynomial_degree() + 1;

    pot_ = DiscretePotential<LaplaceVectorPotentialSL<LaplaceHypersingularOperator>,
                      LaplaceHypersingularOperator>(ansatz_space);

    //RK4 is implemented
    ButcherTable.resize(5,5);


    //FUTUR WORK: Make different integration schemes available! So far RK4 is
    //hard coded.
    //initialize the butcher table
    /*
    ButcherTable(1,0) = 0.5;
    ButcherTable(2,0) = 0.5;
    ButcherTable(3,0) = 1.;

    ButcherTable(1,1) = 0.5;
    ButcherTable(4,1) = 1./6.;

    ButcherTable(2,2) = 0.5;
    ButcherTable(4,2) = 1./3.;

    ButcherTable(3,3) = 1.;
    ButcherTable(4,3) = 1./3.;

    ButcherTable(4,4) = 1./6.;
    */

    return;
  }


  //////////////////////////////////////////////////////////////////////////////
  //    setter
  //////////////////////////////////////////////////////////////////////////////
  void set_cauchy_data(const Eigen::Matrix<double,Eigen::Dynamic, 1> &cauchy_data) {
    fun_ev_.set_function(cauchy_data);
    pot_.set_cauchy_data(cauchy_data);
  }

  void set_initial_condition(const Eigen::Matrix<double,Eigen::Dynamic,6> &W_0){

    num_particles_ = W_0.rows();
    W_0_ = W_0;

  }

  void set_initial_condition(const Eigen::Matrix<double,Eigen::Dynamic,6> &W_0,const Eigen::VectorXd &p_s){

    num_particles_ = W_0.rows();
    W_0_ = W_0;

    p_vec_ = p_s;

  }

  void set_reference_trajectory(ReferenceTrajectory *traj){

    traj_ = traj;

  }

  void set_step_size(double step_size){

    step_ = step_size;
  }

  Eigen::MatrixXd integrate(double step_size,double s_max){

    //w = [x,y,t,Px,Py,Pt]

    int num_steps = s_max/step_size-1;

    std::cout << "*************************************" << std::endl;
    std::cout << "Launching Runge Kutta Integration" << std::endl;
    std::cout << "*************************************" << std::endl;
    std::cout << "\tNumber of particles = " << num_particles_ << std::endl;
    std::cout << "\tNumber of RK steps = " << num_steps << std::endl;

    Eigen::MatrixXd ret_val;
    ret_val.resize(num_particles_,6);

    Eigen::MatrixXd W_n, W_np1;
    W_n.resize(num_particles_,6); W_np1.resize(num_particles_,6);


    //make space for results
    X_.resize(num_particles_,num_steps);
    Y_.resize(num_particles_,num_steps);
    Px_.resize(num_particles_,num_steps);
    Py_.resize(num_particles_,num_steps);
    T_.resize(num_particles_,num_steps);
    Pt_.resize(num_particles_,num_steps);
    s_.resize(num_steps);

    //for RK4
    Eigen::VectorXd k1,k2,k3,k4;

    double s_n = 0.;

    W_n = W_0_;


    for (int i = 0; i < num_steps; ++i){


      for(int j = 0; j < num_particles_; ++j){

        if(test_enable_){
          k1 = f_dot(W_n.row(j) , s_n);
          k2 = f_dot(W_n.row(j) + step_size*k1.transpose()/2, s_n + step_size/2.);
          k3 = f_dot(W_n.row(j) + step_size*k2.transpose()/2, s_n + step_size/2.);
          k4 = f_dot(W_n.row(j) + step_size*k3.transpose(),   s_n + step_size);
        }
        else{
          p_ = p_vec_(j);

          k1 = w_dot(W_n.row(j) , s_n);
          k2 = w_dot(W_n.row(j) + step_size*k1.transpose()/2, s_n + step_size/2.);
          k3 = w_dot(W_n.row(j) + step_size*k2.transpose()/2, s_n + step_size/2.);
          k4 = w_dot(W_n.row(j) + step_size*k3.transpose(),   s_n + step_size);
        }


        W_np1.row(j) = W_n.row(j) + step_size/6.*(k1.transpose()
                                                  + 2*k2.transpose()
                                                   + 2*k3.transpose()
                                                    + k4.transpose());

      }

      //increment
      s_n += step_size;
      W_n = W_np1;

      //store variables
      X_.col(i) = W_np1.col(0);
      Y_.col(i) = W_np1.col(1);
      T_.col(i) = W_np1.col(2);
      Px_.col(i) = W_np1.col(3);
      Py_.col(i) = W_np1.col(4);
      Pt_.col(i) = W_np1.col(5);
      s_(i) = s_n;



    }

    return W_np1;

  }

  void save_trajectories(const std::string filename){

    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", "\n");

    Eigen::MatrixXd tmp_mat;
    tmp_mat.resize(X_.cols(),1+6*num_particles_);

    tmp_mat.col(0) = s_;

    for(int j = 0; j < num_particles_; ++j){

      tmp_mat.col(1+6*j) = X_.row(j);
      tmp_mat.col(2+6*j) = Y_.row(j);
      tmp_mat.col(3+6*j) = T_.row(j);
      tmp_mat.col(4+6*j) = Px_.row(j);
      tmp_mat.col(5+6*j) = Py_.row(j);
      tmp_mat.col(6+6*j) = Pt_.row(j);

    }
    std::ofstream file_path(filename);
    file_path << tmp_mat.format(CSVFormat);

  }

  void set_test(bool enable){

    test_enable_ = enable;

  }

  void extract_integrated_field(int multipole_order,double R,int num_samples_s = 500){

    Eigen::Vector3d start_pos = traj_->get_pos(0.);
    Eigen::Vector3d end_pos = traj_->get_end();

    int num_angles = 3*multipole_order;
    double s_0 = start_pos(2);
    double s_1 = end_pos(2);
    double angle_step = 2.*BEMBEL_PI/num_angles;
    double s_step = (s_1-s_0)/num_samples_s;

    Eigen::MatrixXd pos;
    pos.resize(num_samples_s,3);

    Eigen::MatrixXd Axys;
    Eigen::VectorXd As_i(num_angles);
    As_i.setZero();


    for(int i = 0 ; i < num_angles ; ++i){

      //std::cout << "phi = " << angle_step*i << " rad" << std::endl;

      //Integrate here over the curved boundary in future!!!
      pos.block(0,0,num_samples_s,1) = R*std::cos(angle_step*i)*Eigen::VectorXd::Ones(num_samples_s);
      pos.block(0,1,num_samples_s,1) = R*std::sin(angle_step*i)*Eigen::VectorXd::Ones(num_samples_s);
      pos.block(0,2,num_samples_s,1) = Eigen::VectorXd::LinSpaced(num_samples_s,s_0,s_1);

      Axys = pot_.evaluate_der(pos);

      for(int j = 0; j < num_samples_s-1 ; ++j){
        As_i(i) += 0.5*(Axys(j+1,2)+Axys(j,2))*s_step;
      }


    }
    Eigen::FFT<double> fft;
    Eigen::VectorXcd Cn;
    Eigen::VectorXd in;
    fft.fwd(Cn, As_i);

    Cn *= 2./num_angles;

    std::cout << "Size Cn = " << Cn.size() << std::endl;


    //Scale Cn's by radius
    for (int n = 0; n < num_angles ; ++n){
      std::cout << Cn(n) << std::endl;

    }

    //double Cref = std::sqrt(Cn(1).real()*Cn(1).real() + Cn(1).imag()*Cn(1).imag());
    //double angle = 0.;

  }


 private:
  int deg_;
  //tracking is (so far) always based on the mvp
  DiscretePotential<LaplaceVectorPotentialSL<LinOp>,LinOp> pot_;
  AnsatzSpace<LinOp> ansatz_space_;
  FunctionEvaluator<LinOp> fun_ev_;
  ReferenceTrajectory* traj_;
  int num_particles_;
  double step_;
  Eigen::MatrixXd X_,Y_,Px_,Py_,T_,Pt_;
  Eigen::VectorXd s_;
  bool test_enable_ = false;

  Eigen::Matrix<double,Eigen::Dynamic,6> W_0_;
  //particle momentum
  Eigen::VectorXd p_vec_;
  double p_;


  //speed of light
  double c_ = 299792458;     //[m/s]

  //elementary charge
  double e_ = 1.602176634e-19;  //[As]

  //mass
  double m_;

  Eigen::MatrixXd ButcherTable;

  //This is a test funcion, to test the numerical integration scheme
  Eigen::VectorXd f_dot(const Eigen::VectorXd &w, double s){

    Eigen::VectorXd ret_val;
    ret_val.resize(6);
    ret_val.setZero();


    ret_val(0) = w(0)*std::cos(s);

    return ret_val;

  }


  Eigen::VectorXd w_dot(const Eigen::VectorXd &w, double s){

    //Frenet–Serret-Frame
    double h = traj_->get_curvature(s);
    Eigen::Matrix<double,1,3> r = traj_->get_pos(s).transpose();
    Eigen::Matrix3d E = traj_->get_unit_vectors(s);

    //Compute MVP
    Eigen::MatrixXd A = pot_.evaluate(r).eval();
    Eigen::MatrixXd dA = pot_.evaluate_der(r).eval();



    //transform to Frenet–Serret-Frame
    A *= E;
    dA = dA.transpose()*E;

    //MVP components in Frenet–Serret-Frame
    double Ax = A(0);
    double Ay = A(1);
    double As = A(2);

    //MVP derivatives in Frenet–Serret-Frame
    double Axx = dA(0,0);
    double Ayx = dA(0,1);
    double Asx = dA(0,2);
    double Axy = dA(1,0);
    double Ayy = dA(1,1);
    double Asy = dA(1,2);


    //return value
    //w = [x,y,t,Px,Py,Pt]
    Eigen::VectorXd w_dot;
    w_dot.resize(6);

    //sqrt in denominator
    double den = std::sqrt(p_*p_ - (w(3)-e_*Ax)*(w(3)-e_*Ax)- (w(4)-e_*Ay)*(w(4)-e_*Ay));


    //dx/ds
    w_dot(0) = (1+h*w(0))*(w(3)-e_*Ax)/den;

    //dy/ds
    w_dot(1) = (1+h*w(0))*(w(4)-e_*Ay)/den;

    //Px/ds
    w_dot(3) = (1+h*w(0))*e_*Asx + h*e_*As + h*den + (1+h*w(0))*((w(3)-e_*Ax)*e_*Axx + (w(4)-e_*Ay)*e_*Ayx)/den;

    //Pz/ds
    w_dot(4) = (1+h*w(0))*e_*Asy + (1+h*w(0))*((w(3)-e_*Ax)*e_*Axy + (w(4)-e_*Ay)*e_*Ayy)/den;

    //t/ds
    w_dot(2) = m_*(1+h*w(0))/den;

    //Pt/ds
    w_dot(5) = 0.;

    return w_dot;

  }


};  // namespace Bembel

}  // namespace Bembel
#endif
