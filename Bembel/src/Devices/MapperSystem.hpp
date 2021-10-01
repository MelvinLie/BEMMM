// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_DEVICES_MAPPERSYSTEM_H_
#define BEMBEL_DEVICES_MAPPERSYSTEM_H_


#include <Eigen/StdVector>

namespace Bembel {

/**
 * \ingroup Devices
 */
template <typename Pot,typename LinOp>
class MapperSystem {

 public:
  MapperSystem(const AnsatzSpace<LinOp> &ansatz_space) {

    ansatz_space_ = ansatz_space;

    z_.push_back(0.);
    num_elements_ = 0;

  }

  /*****************************************************************************
  Read mapper arm geometry from file.

  *****************************************************************************/
  void read_geometry_from_file(const std::string &filename){

    char del = ',';
    std::ifstream file;
    file.open(filename);
    std::string current_line;
    std::string current_element;
    std::stringstream current_data;

    std::vector<double> tmp_vec;
    double I;
    int segs;

    if (!file) {
      std::cerr << "File " << filename << " doesn't exist!" << std::endl;
      exit(1);
    }

    // first  row is irrelevant
    getline(file, current_line);

    while(getline(file, current_line)){

      tmp_vec.clear();
      std::stringstream tmp_line(current_line);

      while(std::getline(tmp_line,current_element,del)){
          tmp_vec.push_back(atof(current_element.c_str()));
      }

      segs = (int) tmp_vec[5];
      //so far only hollow cylinders are implemented
      I = 0.25*BEMBEL_PI*(std::pow(tmp_vec[3]/1e3,4) - std::pow(tmp_vec[2]/1e3,4));

      for(int i = 0; i < segs;++i){

        add_segment(tmp_vec[1]*1e9, I, tmp_vec[4], tmp_vec[0]/segs*1e-3);

      }
    }

    file.close();



  }
  void print_geometry(){

    std::cout << "Euler-Bernoulli Beam with of length " << get_length()*1e3 <<" mm" << std::endl;
    std::cout << "|    E [GPa]   |    I   [mm^4]   |    M    [kg/m]    |    l   [mm]     |" << std::endl;
    std::cout << "|--------------|-----------------|-------------------|-----------------|" << std::endl;
    for(int i = 0; i < num_elements_; ++i){
      std::cout << "| " << std::setw(12) << E_[i]*1e-9  << " | " << std::setw(15) << I_[i]*1e12 << " | " << std::setw(17) << m_[i] << " | " << std::setw(15) << h_[i] << " |" << std::endl;
    }

  }

  double get_length(){

    return z_.back();

  }


  void add_segment(double e_mod, double inertia, double mass, double length){

    E_.push_back(e_mod);
    I_.push_back(inertia);
    m_.push_back(mass);
    h_.push_back(length);
    z_.push_back(z_.back() + length);

    num_elements_ += 1;

  }

  void add_segment(double e_mod, double r_i, double r_o, double mass, double length){

    E_.push_back(e_mod);
    I_.push_back(BEMBEL_PI/4*(r_o*r_o*r_o*r_o + r_i*r_i*r_i*r_i));
    m_.push_back(mass);
    h_.push_back(length);
    z_.push_back(z_.back() + length);


    num_elements_ += 1;

  }

  void set_damping_coeffs(double l, double u){

    l_ = l;
    u_ = u;

  }


  void initialize_galerkin_matrices(){


    num_dofs_ = num_elements_*2;

    M_ = Eigen::SparseMatrix<double>(num_dofs_,num_dofs_);
    C_ = Eigen::SparseMatrix<double>(num_dofs_,num_dofs_);
    K_ = Eigen::SparseMatrix<double>(num_dofs_,num_dofs_);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tM,tK;
    tM.reserve(6*num_dofs_+4);
    tK.reserve(6*num_dofs_+4);

    //index of upper left corner
    int ulc;

    //std::cout << "Number of DoFs FEM = " << num_dofs << std::endl;

    double m_fac,k_fac;

    for(int i = 0; i < num_elements_; ++i){

      m_fac = m_[i]*h_[i]/420;


      k_fac = E_[i]*I_[i]/h_[i]/h_[i]/h_[i];


      if(i == 0){

        tM.push_back(T(0,0, 156*m_fac));
        tM.push_back(T(0,1, -22*m_fac*h_[i]));
        tM.push_back(T(1,0, -22*m_fac*h_[i]));
        tM.push_back(T(1,1,  4*m_fac*h_[i]*h_[i]));

        //std::cout <<"m = " <<m_[i]<< std::endl;
        //std::cout <<"h = " <<h_[i]<< std::endl;
        //std::cout <<"first element = " <<156*m_fac<< std::endl;

        tK.push_back(T(0,0,  12*k_fac));
        tK.push_back(T(0,1, -6*k_fac*h_[i]));
        tK.push_back(T(1,0, -6*k_fac*h_[i]));
        tK.push_back(T(1,1,  4*k_fac*h_[i]*h_[i]));

      }
      else{

        ulc = 2*(i-1);

        tM.push_back(T(ulc  ,ulc  , 156*m_fac));
        tM.push_back(T(ulc  ,ulc+1, 22*m_fac*h_[i]));
        tM.push_back(T(ulc  ,ulc+2, 54*m_fac));
        tM.push_back(T(ulc  ,ulc+3,-13*m_fac*h_[i]));

        tM.push_back(T(ulc+1,ulc  , 22*m_fac*h_[i]));
        tM.push_back(T(ulc+1,ulc+1, 4*m_fac*h_[i]*h_[i]));
        tM.push_back(T(ulc+1,ulc+2, 13*m_fac*h_[i]));
        tM.push_back(T(ulc+1,ulc+3,-3*m_fac*h_[i]*h_[i]));

        tM.push_back(T(ulc+2,ulc  , 54*m_fac));
        tM.push_back(T(ulc+2,ulc+1, 13*m_fac*h_[i]));
        tM.push_back(T(ulc+2,ulc+2, 156*m_fac));
        tM.push_back(T(ulc+2,ulc+3,-22*m_fac*h_[i]));

        tM.push_back(T(ulc+3,ulc  ,-13*m_fac*h_[i]));
        tM.push_back(T(ulc+3,ulc+1,-3*m_fac*h_[i]*h_[i]));
        tM.push_back(T(ulc+3,ulc+2,-22*m_fac*h_[i]));
        tM.push_back(T(ulc+3,ulc+3, 4*m_fac*h_[i]*h_[i]));


        tK.push_back(T(ulc  ,ulc  , 12*k_fac));
        tK.push_back(T(ulc  ,ulc+1, 6*k_fac*h_[i]));
        tK.push_back(T(ulc  ,ulc+2,-12*k_fac));
        tK.push_back(T(ulc  ,ulc+3, 6*k_fac*h_[i]));

        tK.push_back(T(ulc+1,ulc  , 6*k_fac*h_[i]));
        tK.push_back(T(ulc+1,ulc+1, 4*k_fac*h_[i]*h_[i]));
        tK.push_back(T(ulc+1,ulc+2,-6*k_fac*h_[i]));
        tK.push_back(T(ulc+1,ulc+3, 2*k_fac*h_[i]*h_[i]));

        tK.push_back(T(ulc+2,ulc  ,-12*k_fac));
        tK.push_back(T(ulc+2,ulc+1, -6*k_fac*h_[i]));
        tK.push_back(T(ulc+2,ulc+2, 12*k_fac));
        tK.push_back(T(ulc+2,ulc+3, -6*k_fac*h_[i]));

        tK.push_back(T(ulc+3,ulc  , 6*k_fac*h_[i]));
        tK.push_back(T(ulc+3,ulc+1, 2*k_fac*h_[i]*h_[i]));
        tK.push_back(T(ulc+3,ulc+2, -6*k_fac*h_[i]));
        tK.push_back(T(ulc+3,ulc+3, 4*k_fac*h_[i]*h_[i]));
      }

    }

    M_.setFromTriplets(tM.begin(), tM.end());


    K_.setFromTriplets(tK.begin(), tK.end());

    C_ = 2.*u_*M_ + l_*K_;

    //vector of additional forces
    f_.resize(num_dofs_/2);
    f_.setZero();

  }


  void setup_load_vector(){

    p_.resize(num_dofs_);
    p_.setZero();

    //setup load vector
    for(int i = 0; i < num_elements_; ++i){
      if(i == 0){
        p_(0) = -g_*m_[i]*h_[i]*0.5;
        p_(1) =  g_*m_[i]*h_[i]*h_[i]/12;
      }

      else{

        p_(2*(i-1))   += -g_*m_[i]*h_[i]*0.5;
        p_(2*(i-1)+1) += -g_*m_[i]*h_[i]*h_[i]/12;
        p_(2*(i-1)+2) += -g_*m_[i]*h_[i]*0.5;
        p_(2*(i-1)+3) +=  g_*m_[i]*h_[i]*h_[i]/12;
      }
    }

    //apply additional forces
    for(int i = 0 ; i < num_dofs_; i = i + 2){
      p_(i) += f_(i/2);
    }

  }

  void apply_force_to_node(const int node,const double f){

    //std::cout << "num dofs / 2 = " << num_dofs_/2 << std::endl;
    f_.resize(num_dofs_/2);
    f_(node) = f;

  }


  void reset_external_forces(){

    f_.setZero();

  }

  Eigen::VectorXd compute_static_displacement(){

    setup_load_vector();


    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    //solver.analyzePattern(K_);
    solver.compute(K_);

    Eigen::VectorXd w = solver.solve(p_);


    Eigen::VectorXd ret_val;
    ret_val.resize(w.size()+2);
    ret_val.setZero();
    ret_val.segment(2,w.size()) = w;

    return ret_val;

  }

  Eigen::VectorXd evaluate_displacement(const Eigen::VectorXd &x, const Eigen::VectorXd &z){

    Eigen::VectorXd w;

    w.resize(z.size());
    w.setZero();

    int seg_indx = -1;

    double z_hat;



    for (int i = 0; i < z.size(); ++i){



      //find segment
      seg_indx = -1;

      for( int j = 0; j < num_elements_; ++j){

        if((z[i] >= z_[j] ) && (z[i] < z_[j+1])){
          seg_indx = j;
        }


      }

      if ((seg_indx < 0) && (std::fabs(z[i] - z_.back())) < 1e-14){
        seg_indx = num_elements_-1;
      }


      if(seg_indx >= 0){

      z_hat = z[i] - z_[seg_indx];


      w[i] = x[2*seg_indx]*(2*std::pow(z_hat,3) - 3*h_[seg_indx]*std::pow(z_hat,2) + std::pow(h_[seg_indx],3))/std::pow(h_[seg_indx],3);
      w[i] += x[2*seg_indx+1]*(h_[seg_indx]*std::pow(z_hat,3) - 2*std::pow(h_[seg_indx],2)*std::pow(z_hat,2) + std::pow(h_[seg_indx],3)*z_hat)/std::pow(h_[seg_indx],3);
      w[i] += x[2*seg_indx+2]*(-2.*std::pow(z_hat,3) + 3*h_[seg_indx]*std::pow(z_hat,2))/std::pow(h_[seg_indx],3);
      w[i] += x[2*seg_indx+3]*(h_[seg_indx]*std::pow(z_hat,3) - std::pow(h_[seg_indx],2)*std::pow(z_hat,2))/std::pow(h_[seg_indx],3);


      }

    }

    return w;

  }

  void setup_support_motion(const std::function<double(double)> &a_fcn,
                            const std::function<double(double)> &a_dot_fcn,
                            const std::function<double(double)> &a_ddot_fcn) {

    a_      = a_fcn;
    a_dot_  = a_dot_fcn;
    a_ddot_ = a_ddot_fcn;

  }

  Eigen::SparseMatrix<std::complex<double>>  get_complex_system_matrix(double omega){

    Eigen::SparseMatrix<std::complex<double>> A(num_dofs_,num_dofs_);

    std::complex<double> j(0.0,1.0);//std::sqrt(-1);

    //std::cout << j << std::endl;
    std::complex<double> tmp;

    typedef Eigen::Triplet<std::complex<double>> T;
    std::vector<T> tA;
    tA.reserve(3*M_.nonZeros());

    for (int k=0; k<M_.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(M_,k); it; ++it)
      {

        tmp = -omega*omega*it.value();
        tA.push_back(T(it.row()  ,it.col()  , tmp));
      }
    }

    for (int k=0; k<C_.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(C_,k); it; ++it)
      {

        tmp = j*omega*it.value();

        tA.push_back(T(it.row()  ,it.col()  , tmp));
      }
    }

    for (int k=0; k<K_.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(K_,k); it; ++it)
      {

        tmp = it.value();
        tA.push_back(T(it.row()  ,it.col()  , tmp));
      }
    }

    A.setFromTriplets(tA.begin(), tA.end());

    return A;

  }

  void modal_analysis(Eigen::VectorXd *freq, Eigen::MatrixXd *modes){

    //matrices
    initialize_galerkin_matrices();

    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd> solver;


    //Unfortulately we need to get dense here...
    solver.compute(Eigen::MatrixXd(K_),Eigen::MatrixXd(M_));

    //std::cout << solver.alphas() << std::endl;
    //std::cout << "betas" << std::endl;
    //std::cout << solver.betas() << std::endl;

    (*freq).resize(solver.eigenvalues().size());
    (*modes).resize(solver.eigenvalues().size(),solver.eigenvalues().size());

    Eigen::VectorXd omegas = solver.eigenvalues().real();
    //std::cout << omegas << std::endl;
    (*freq) = omegas.cwiseSqrt();
    (*freq) *= 1./2./BEMBEL_PI;
    (*modes) = solver.eigenvectors().real();

  }

  Eigen::MatrixXcd run_bode_analysis(const Eigen::MatrixXd &f){

    int num_freq = f.size();

    Eigen::MatrixXcd ret_val(num_dofs_+2,num_freq);
    ret_val.setZero();

    //matrices
    initialize_galerkin_matrices();
    Eigen::SparseMatrix<std::complex<double>> A;
    A.resize(num_dofs_,num_dofs_);

    double omega;
    std::complex<double> j(0.0,1.0);

    Eigen::VectorXcd g(num_dofs_);
    g.setZero();

    Eigen::VectorXcd x_n(num_dofs_);

    //solver
    //Eigen::SimplicialLLT<Eigen::SparseMatrix<std::complex<double>>> solver;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXcd> solver;


    for(int n = 0; n < num_freq; ++n){

      //angular frequency
      omega = 2*BEMBEL_PI*f(n);

      //system matrix
      A.setZero();
      A = get_complex_system_matrix(omega);

      //rhs
      g.setZero();
      g(0) = -m_[0]*h_[0]/420*(54.*(-omega*omega+2.*u_*j*omega));
      g(0) += E_[0]*I_[0]/std::pow(h_[0],3)*(12.*(1.+l_*j*omega));

      g(1) =  m_[0]*h_[0]/420*(13.*h_[0]*(-omega*omega+2*u_*j*omega));
      g(1) -= E_[0]*I_[0]/std::pow(h_[0],2)*(6.*(1.+l_*j*omega));


      //solve
      solver.compute(Eigen::MatrixXcd(A));
      x_n = solver.solve(g);

      ret_val.block(2,n,num_dofs_,1) = x_n;
      ret_val(0,n) = 1.;

    }

    return ret_val;

  }

  void setup_initial_condition(const Eigen::VectorXd &w_init){

    init_set_ = true;

    w_init_ = w_init;
  }

  void set_gravity(double g){

    g_ = g;
  }


  Eigen::MatrixXd run_newmark(double dt,double T_max,double beta = 1./6., double gamma = 0.5){

    //matrices
    initialize_galerkin_matrices();

    //load vector
    setup_load_vector();

    //initial deformation
    Eigen::VectorXd w_s;

    if(init_set_){
      //compute static deformation
      w_s = w_init_;
    }
    else{
      //compute static deformation
      w_s = compute_static_displacement();
    }


    //initial condition
    Eigen::VectorXd w_dot_s(w_s.size());
    Eigen::VectorXd w_ddot_s(w_s.size());
    w_dot_s.setZero();
    w_ddot_s.setZero();

    for(int i = 0; i < w_s.size()/2 ; ++i){
      w_s[2*i] += a_(0.);
      w_dot_s[2*i] += a_dot_(0.);
      w_ddot_s[2*i] += a_ddot_(0.);
    }

    //system matrix
    Eigen::SparseMatrix<double> A = M_/beta/dt/dt + gamma*C_/beta/dt + K_;

    //Cholesky factorization
    //Undamped case -> Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver;
    solver.analyzePattern(A);
    solver.factorize(A);

    //solver.compute(A);

    //make space for data
    int num_steps = T_max/dt;

    Eigen::MatrixXd ret_val(num_steps,num_dofs_+2);
    ret_val.row(0) = w_s;


    Eigen::VectorXd b(num_dofs_);

    Eigen::VectorXd g(num_dofs_);
    g.setZero();

    Eigen::VectorXd t1(num_dofs_);
    Eigen::VectorXd t2(num_dofs_);

    Eigen::VectorXd w_n = w_s.segment(2,num_dofs_);
    Eigen::VectorXd w_dot_n = w_dot_s.segment(2,num_dofs_);
    Eigen::VectorXd w_ddot_n = w_ddot_s.segment(2,num_dofs_);

    Eigen::VectorXd w_np1(num_dofs_);
    Eigen::VectorXd w_dot_np1(num_dofs_);
    Eigen::VectorXd w_ddot_np1(num_dofs_);


    //iteration
    double t;

    for (int n = 1; n < num_steps; ++n){

      t = n*dt;


      //support condition
      g(0) = -m_[0]*h_[0]/420*(54*(a_ddot_(t)+2*u_*a_dot_(t)));
      g(0) += E_[0]*I_[0]/std::pow(h_[0],3)*(12*(a_(t)+l_*a_dot_(t)));

      g(1) =  m_[0]*h_[0]/420*(13*h_[0]*(a_ddot_(t)+2*u_*a_dot_(t)));
      g(1) -= E_[0]*I_[0]/std::pow(h_[0],2)*(6*(a_(t)+l_*a_dot_(t)));

      t1 = w_n/beta/dt/dt + w_dot_n/beta/dt + (0.5/beta - 1.)*w_ddot_n;
      t2 = gamma*w_n/beta/dt - w_dot_n*(1-gamma/beta) - dt*w_ddot_n*(1-0.5*gamma/beta);

      b = p_ + g + M_*t1 + C_*t2;


      //update
      w_np1 = solver.solve(b);

      Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>>   solver_tmp;
      solver_tmp.analyzePattern(A);
      solver_tmp.factorize(A);

      Eigen::VectorXd w_test = solver_tmp.solve(b);


      w_dot_np1 = gamma*(w_np1 - w_n)/beta/dt + w_dot_n*(1-gamma/beta) + dt*(1-0.5*gamma/beta)*w_ddot_n;
      w_ddot_np1 = (w_np1 - w_n)/beta/dt/dt - w_dot_n/beta/dt - (0.5/beta-1)*w_ddot_n;

      ret_val.block(n,2,1,num_dofs_) = w_np1.transpose();
      ret_val(n,0) = a_(t);

      //increment
      w_n = w_np1;
      w_dot_n = w_dot_np1;
      w_ddot_n = w_ddot_np1;



    }

    return ret_val;

  }

  Eigen::SparseMatrix<double> get_M(){

    return M_;
  }

  Eigen::SparseMatrix<double> get_K(){

    return K_;
  }


  Eigen::SparseMatrix<double> get_C(){

    return C_;
  }



 private:
   //Hall Probes
   std::vector<DiscreteSensor<HallProbe<Pot,LinOp>,LinOp>> sensors_;
   //positions
   std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> positions_;
   //ansatz space
   AnsatzSpace<LinOp> ansatz_space_;
   //aca matrices
   std::vector<ACAMatrix<HallProbe<Pot,LinOp>,LinOp>> aca_matrices_;
   //number of evaluations
   int num_eval_;
   //flag for gauged formulation
   bool enable_gauge_;
   //stability vector
   Eigen::VectorXd stability_vec_;

   //beam segments
   //elastic modulus
   std::vector<double> E_;
   //Second moment of area
   std::vector<double> I_;
   //segment length
   std::vector<double> h_;
   //knots
   std::vector<double> z_;
   //Mass per unit length
   std::vector<double> m_;
   //number of elements
   int num_elements_ = 0;

   //Rayleigh damping
   double l_ = 0.;
   double u_ = 0.;

   //earths gravity
   double g_ = 9.80665; //[m/s^2]

   int num_dofs_;

   //Galerkin matrices
   Eigen::SparseMatrix<double> M_, K_, C_;

   //Load vector
   Eigen::VectorXd p_;

   //additional force vector
   Eigen::VectorXd f_;

   //initial condition
   Eigen::VectorXd w_init_;

   //flag initial condition
   bool init_set_ = false;

   std::function<double(double)> a_;
   std::function<double(double)> a_dot_;
   std::function<double(double)> a_ddot_;

};

}  // namespace Bembel
#endif
