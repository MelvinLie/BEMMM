// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_DEVICES_HALLPROBEARRAY_H_
#define BEMBEL_DEVICES_HALLPROBEARRAY_H_


#include <Eigen/StdVector>

namespace Bembel {

/**
 * \ingroup Devices
 */
template <typename Pot,typename LinOp>
class HallProbeArray {

 public:
  HallProbeArray(const AnsatzSpace<LinOp> &ansatz_space) {

    ansatz_space_ = ansatz_space;

    move_ids_.clear();

    Pert_Covs_.clear();

  }

  void load_parameter_from_file(const std::string &filename){

    char del = ',';
    std::ifstream file;
    file.open(filename);
    std::string current_line;
    std::string current_element;
    std::stringstream current_data;
    int i = 0;
    int k;
    int num_probes;


    std::vector<double> tmp_vector;

    std::vector<double> q;
    Eigen::Vector3d n;
    Eigen::Vector3d p;

    noise_stdev_.clear();

    if (!file) {
      std::cerr << "File " << filename << " doesn't exist!" << std::endl;
      exit(1);
    }

    // first  two rows is irrelevant
    getline(file, current_line);
    getline(file, current_line);

    //number of probes
    //covert this line
    getline(file, current_line);
    num_probes = atoi(current_line.c_str());

    //irrelevant line
    getline(file, current_line);

    for(; i < num_probes ; ++i){

      q.clear();

      //irrelevant line
      getline(file, current_line);
      //irrelevant line
      getline(file, current_line);

      getline(file, current_line);
      std::stringstream q_line(current_line);
      while(std::getline(q_line,current_element,del)){
          q.push_back(atof(current_element.c_str()));
      }

      //irrelevant line
      getline(file, current_line);
      //irrelevant line
      getline(file, current_line);

      getline(file, current_line);
      std::stringstream p_line(current_line);

      k = 0;
      while(std::getline(p_line,current_element,del)){
          p(k) = atof(current_element.c_str())*1e-3;
          ++k;
      }

      //irrelevant line
      getline(file, current_line);
      //irrelevant line
      getline(file, current_line);

      getline(file, current_line);
      std::stringstream n_line(current_line);

      k = 0;
      while(std::getline(n_line,current_element,del)){
          n(k) = atof(current_element.c_str());
          ++k;
      }

      //irrelevant line
      getline(file, current_line);
      //irrelevant line
      getline(file, current_line);



      //noise standard deviation
      double noise_std;
      getline(file, current_line);
      std::stringstream std_line(current_line);

      while(std::getline(std_line,current_element,del)){
          noise_std = atof(current_element.c_str());
      }

      //irrelevant line
      getline(file, current_line);
      //irrelevant line
      getline(file, current_line);

      append_probe(n,p,q,noise_std);

    }

    file.close();


  }

  void append_probe(const Eigen::Vector3d &n,  const Eigen::Vector3d &pos, std::vector<double> tf, const double std = 0.) {

    positions_.push_back(pos);

    DiscreteSensor<HallProbe<Pot,LinOp>,LinOp> disc_sensor(ansatz_space_);

    disc_sensor.get_sensor().set_parameters(n,tf);

    sensors_.push_back(disc_sensor);

    noise_stdev_.push_back(std);


    //aca_matrices_.resize(sensors_.size());

  }

  void print_parameters(){

    std::cout << "Array of " << sensors_.size() << " Hall probes" << std::endl;

    Eigen::Vector3d this_n;
    std::vector<double> this_tf;

    for (int i = 0 ; i < sensors_.size() ; ++i){

      this_n = sensors_[i].get_sensor().get_orientation();
      this_tf = sensors_[i].get_sensor().get_transfer_function();

      std::cout << "Probe " << i+1 << ":" << std::endl;
      std::cout << "\tn = ( " << this_n(0) << " , "<< this_n(1)<< " , "<< this_n(2)<< " )" << std::endl;
      std::cout << "\tpos = ( " << positions_[i](0) << " , "<< positions_[i](1)<< " , "<< positions_[i](2)<< " )" << std::endl;
      std::cout << "\tnoise stdev = " << noise_stdev_[i] << std::endl;
      std::cout << "\ttf = " << this_tf[0] ;
      for(int j = 1 ; j < this_tf.size() ; ++j)  std::cout << " + " <<  this_tf[j] << " * B**" << j;
      std::cout << std::endl;

    }

  }

  Eigen::VectorXd evaluate_H(const Eigen::VectorXd &v){

    int M = H_mat.rows();
    int N = H_mat.cols();

    return H_mat.block(0,1,M,N-1) * v.segment(1,N-1);

  }

  Eigen::VectorXd evaluate_H_tilde(const Eigen::VectorXd &v){

    int M = H_mat.rows();
    int N = H_mat.cols();

    return H_tilde_.block(0,1,M,N-1) * v.segment(1,N-1);

  }


  Eigen::MatrixXd evaluate_grad_H(const Eigen::VectorXd &v){

    int M = H_mat.rows();
    int N = H_mat.cols();

    Eigen::MatrixXd ret_mat(M,3);

    ret_mat.col(0) = dHx_mat.block(0,1,M,N-1)*v.segment(1,N-1);
    ret_mat.col(1) = dHy_mat.block(0,1,M,N-1)*v.segment(1,N-1);
    ret_mat.col(2) = dHz_mat.block(0,1,M,N-1)*v.segment(1,N-1);

    return ret_mat;

  }

  Eigen::VectorXd evaluate_grad_H(const Eigen::VectorXd &v,
                                    const Eigen::VectorXd &dw_x,
                                    const Eigen::VectorXd &dw_y,
                                    const Eigen::VectorXd &dw_z){

    int M = H_mat.rows()/3;
    int N = H_mat.cols();

    Eigen::VectorXd dH_dx = dHx_mat.block(0,1,3*M,N-1)*v.segment(1,N-1);
    Eigen::VectorXd dH_dy = dHy_mat.block(0,1,3*M,N-1)*v.segment(1,N-1);
    Eigen::VectorXd dH_dz = dHz_mat.block(0,1,3*M,N-1)*v.segment(1,N-1);

    Eigen::VectorXd ret_vec(3*M);
    ret_vec.setZero();

    ret_vec.segment(0,M)     +=  (dH_dx.segment(0,M).array()   * dw_x.array()).matrix();



    ret_vec.segment(M,M)     +=  (dH_dx.segment(M,M).array()   * dw_x.array()).matrix();
    ret_vec.segment(2*M,M)   +=  (dH_dx.segment(2*M,M).array() * dw_x.array()).matrix();

    ret_vec.segment(0,M)     +=  (dH_dy.segment(0,M).array()   * dw_y.array()).matrix();
    ret_vec.segment(M,M)     +=  (dH_dy.segment(M,M).array()   * dw_y.array()).matrix();
    ret_vec.segment(2*M,M)   +=  (dH_dy.segment(2*M,M).array() * dw_y.array()).matrix();

    ret_vec.segment(0,M)     +=  (dH_dz.segment(0,M).array()   * dw_z.array()).matrix();
    ret_vec.segment(M,M)     +=  (dH_dz.segment(M,M).array()   * dw_z.array()).matrix();
    ret_vec.segment(2*M,M)   +=  (dH_dz.segment(2*M,M).array() * dw_z.array()).matrix();

    return ret_vec;

  }


  Eigen::MatrixXd evaluate_dH_dt(const Eigen::VectorXd &v){

    int M = H_mat.rows();
    int N = H_mat.cols();

    Eigen::MatrixXd ret_mat(M,2);

    ret_mat.col(0) = dtH_1.block(0,1,M,N-1)*v.segment(1,N-1);
    ret_mat.col(1) = dtH_2.block(0,1,M,N-1)*v.segment(1,N-1);

    return ret_mat;

  }

  Eigen::VectorXd evaluate_dH_dt(const Eigen::VectorXd &v,
                                    const Eigen::VectorXd &dt_1,
                                    const Eigen::VectorXd &dt_2){

    int M = H_mat.rows()/3;
    int N = H_mat.cols();

    Eigen::VectorXd dH_dt1 = dtH_1.block(0,1,3*M,N-1)*v.segment(1,N-1);
    Eigen::VectorXd dH_dt2 = dtH_2.block(0,1,3*M,N-1)*v.segment(1,N-1);

    Eigen::VectorXd ret_vec(3*M);
    ret_vec.setZero();

    ret_vec.segment(0,M)     +=  (dH_dt1.segment(0,M).array()   * dt_1.array()).matrix();
    ret_vec.segment(M,M)     +=  (dH_dt1.segment(M,M).array()   * dt_1.array()).matrix();
    ret_vec.segment(2*M,M)   +=  (dH_dt1.segment(2*M,M).array() * dt_1.array()).matrix();

    ret_vec.segment(0,M)     +=  (dH_dt2.segment(0,M).array()   * dt_2.array()).matrix();
    ret_vec.segment(M,M)     +=  (dH_dt2.segment(M,M).array()   * dt_2.array()).matrix();
    ret_vec.segment(2*M,M)   +=  (dH_dt2.segment(2*M,M).array() * dt_2.array()).matrix();

    return ret_vec;

  }


  //here we assemble all matrices needed for UQ.
  void prepare_UQ_matrices(const Eigen::MatrixXd &pos,
                           const int arm_axis_id){

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //set mapper axis in Hall Probes
    for(int i = 0 ; i < num_probes ; ++i){
      sensors_[i].get_sensor().setup_transversal_vectors(arm_axis_id);
    }

    compute_UQ_matrices(pos);

    //make a discrete potential
    DiscretePotential<Pot,LinOp> disc_pot;
    disc_pot.init_DiscretePotential(ansatz_space_);

    //impose a zero mean gauge gently
    //use block(0,1,M,N-1) lateron
    disc_pot.impose_zero_mean_gauge_low_ram(H_mat,dtH_1,dtH_2);
    disc_pot.impose_zero_mean_gauge_low_ram(dHx_mat,dHy_mat,dHz_mat);

    //std::cout << "H_mat = " << H_mat.block(0,0,5,5) << std::endl;


    arm_axis_id_ = arm_axis_id;

  }


  void compute_UQ_matrices(const Eigen::MatrixXd &pos){

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();
    //note that we will gauge for zero mean arithmetically.

    //make space for all matrices
    H_mat.resize(num_probes*M , N);
    H_mat.setZero();

    dtH_1.resize(num_probes*M , N);
    dtH_1.setZero();

    dtH_2.resize(num_probes*M , N);
    dtH_2.setZero();

    dHx_mat.resize(num_probes*M , N);
    dHx_mat.setZero();

    dHy_mat.resize(num_probes*M , N);
    dHy_mat.setZero();

    dHz_mat.resize(num_probes*M , N);
    dHz_mat.setZero();

    //evaluation positions
    Eigen::MatrixXd sensor_pos;
    sensor_pos.resize(M,3);

    //sensor transfer function
    std::vector<double> tf;

    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){

      //shift the positions to Hall element
      sensor_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      sensor_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      sensor_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);

      //Assemble matrices
      sensors_[i].compute_UQ_matrices(sensor_pos, H_mat,
                                                    dtH_1,
                                                    dtH_2,
                                                    dHx_mat,
                                                    dHy_mat,
                                                    dHz_mat,
                                                    M*i);

      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();

      //apply tranfer function of this sensor
      H_mat.block(i*M,0,M,N) *= tf[1];
      dtH_1.block(i*M,0,M,N) *= tf[1];
      dtH_2.block(i*M,0,M,N) *= tf[1];

      dHx_mat.block(i*M,0,M,N) *= tf[1];
      dHy_mat.block(i*M,0,M,N) *= tf[1];
      dHz_mat.block(i*M,0,M,N) *= tf[1];
    }
  }

  //here we avoid the assembly of dense UQ matrices and directly compute the derivatives.
  Eigen::MatrixXd compute_derivatives(const Eigen::MatrixXd &pos,
                                      const Eigen::VectorXd &v,
                                      const int arm_axis_id){

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    //set mapper axis in Hall Probes
    for(int i = 0 ; i < num_probes ; ++i){
      sensors_[i].get_sensor().setup_transversal_vectors(arm_axis_id);
    }

    //this matrix stores ( dHx, dHy, dHz, H,  dHt1, dHt2 )
    Eigen::MatrixXd derivative_mat(num_probes*M,6);
    derivative_mat.setZero();


    //temporal storage for derivatives
    Eigen::MatrixXd der_tmp;
    der_tmp.resize(M , 6);

    //evaluation positions
    Eigen::MatrixXd sensor_pos;
    sensor_pos.resize(M,3);

    //Helper unit vector
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);

    //sensor transfer function
    std::vector<double> tf;

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){

      //shift the positions to Hall element
      sensor_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      sensor_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      sensor_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);

      //zero temporal storage
      der_tmp.setZero();

      //compute derivatives
      der_tmp = sensors_[i].compute_voltage_and_derivatives(sensor_pos, v);

      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();

      //write into derivative matrix
      derivative_mat.block(i*M,0,M,6) = tf[1]*der_tmp;


    }

    return derivative_mat;

  }


  //give measurement data as a pointer in future
  Eigen::SparseMatrix<double,Eigen::RowMajor> compute_R_decorrelated(MeasurementData meas,
                                                                      const Eigen::VectorXd &v,
                                                                      const int arm_axis_id){

    //compute the derivatives of the measurement operation
    std::cout << " computing derivatives...";
    Eigen::MatrixXd H_der = compute_derivatives(meas.get_positions(),v,arm_axis_id);
    std::cout << " ... done!" << std::endl;

    //storage scheme H_der: ( dHx, dHy, dHz, H,  dHt1, dHt2 )
    Eigen::DiagonalMatrix<double, -1> dHx,dHy,dHz,dHt1,dHt2;

    //move separator
    std::vector<int> move_sep = meas.get_move_separator();

    //move identifier
    std::vector<int> move_ids = meas.get_move_ids();

    //setup the move precision matrices
    //setup_precision_matrices(move_sep, move_ids);

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //number of measurements
    int num_meas = meas.get_number_of_measurements();

    //count nonzeros in global matrix
    int nonzeros = cout_nonzero_cov(move_sep,num_probes);

    std::cout << " nonzeros in cov = " << nonzeros << std::endl;
    std::cout << " sparsity of cov = " << (double) nonzeros / num_probes/num_meas/num_probes/num_meas << std::endl;


    //make tripletList to fill sparse moment matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(nonzeros);

    //iterate over moves
    #pragma omp parallel for
    for(int m = 0; m < move_sep.size()-1 ; ++m){

      //number of positions in this move
      int num_pos = move_sep[m+1] - move_sep[m];

      //position of this move in memory
      int move_pos = get_move_index(move_ids[m]);

      //get the covariance matrix of this move
      Eigen::MatrixXd cov = get_cov(move_pos,num_pos);

      Eigen::MatrixXd this_cov(3*num_pos,3*num_pos);
      this_cov.setZero();

      //get local observation operator
      Eigen::SparseMatrix<double,Eigen::RowMajor> dH_loc = get_dH_loc(H_der,move_sep[m],num_pos,num_meas);


      //compute local covatiance
      this_cov = dH_loc * cov * dH_loc.transpose();

      //sort this covariance into the global matrix
      #pragma omp critical
      {
        for(int i = 0 ; i < num_probes; ++i){
          for(int j = 0 ; j < num_probes ; ++j){

            for(int k = 0; k < num_pos; ++k){
              for(int l = 0; l < num_pos; ++l){

                tripletList.push_back(T(move_sep[m] + k + i*num_meas , move_sep[m] + l + j*num_meas  , this_cov(k + i*num_pos, l + j*num_pos )));
              }
            }
          }
        }
      }
    }

    //make sparse matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> ret_mat(num_probes*num_meas,num_probes*num_meas);
    ret_mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return ret_mat;
  }

  int  cout_nonzero_cov(const std::vector<int> move_sep, const int num_probes){

    int count = 0;

    //iterate over moves
    for(int m = 0; m < move_sep.size()-1 ; ++m){

      //number of positions in this move
      int num_pos = move_sep[m+1] - move_sep[m];

      count += num_probes*num_probes*num_pos*num_pos;

    }

    return count;


  }

  Eigen::SparseMatrix<double,Eigen::RowMajor> get_dH_loc(const Eigen::MatrixXd &H_der,
                                                          const int m_from,
                                                          const int num_pos,
                                                          const int num_meas){

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //make tripletList to fill sparse moment matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(num_pos*num_pos*5);

    for(int i = 0 ; i < num_probes; ++i){
      for(int m = 0; m < num_pos; ++m){


        tripletList.push_back(T(m + num_pos*i , m             , H_der(m_from + m + i*num_meas, 0 )));
        tripletList.push_back(T(m + num_pos*i , m + num_pos   , H_der(m_from + m + i*num_meas, 1 )));
        tripletList.push_back(T(m + num_pos*i , m + 2*num_pos , H_der(m_from + m + i*num_meas, 2 )));

        tripletList.push_back(T(m + num_pos*i , m + 3*num_pos , H_der(m_from + m + i*num_meas, 4 )));
        tripletList.push_back(T(m + num_pos*i , m + 4*num_pos , H_der(m_from + m + i*num_meas, 5 )));
      }

    }

    //make sparse matrix
    Eigen::SparseMatrix<double,Eigen::RowMajor> ret_mat(num_probes*num_pos,5*num_pos);
    ret_mat.setFromTriplets(tripletList.begin(), tripletList.end());

    return ret_mat;
  }

  void compute_H_matrices(const Eigen::MatrixXd &pos){

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();
    //note that we will gauge for zero mean arithmetically.

    //make space for H matrices
    H_mat.resize(3*M,N);
    H_mat.setZero();

    dtH_1.resize(3*M,N);
    dtH_1.setZero();

    dtH_2.resize(3*M,N);
    dtH_2.setZero();

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //evaluation positions
    Eigen::MatrixXd sensor_pos;
    sensor_pos.resize(M,3);

    //sensor transfer function
    std::vector<double> tf;

    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){

      //shift the positions to Hall element
      sensor_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      sensor_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      sensor_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);

      //Assemble Bn part
      sensors_[i].compute_H_matrices(sensor_pos, H_mat, dtH_1, dtH_2,M*i);

      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();

      //apply tranfer function of this sensor
      H_mat.block(i*M,0,M,N) *= tf[1];
      dtH_1.block(i*M,0,M,N) *= tf[1];
      dtH_2.block(i*M,0,M,N) *= tf[1];


    }
  }

  void compute_H(const Eigen::MatrixXd &pos){

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();
    //note that we will gauge for zero mean arithmetically.

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //make space for H matrices
    H_mat.resize(num_probes*M,N);
    H_mat.setZero();

    //evaluation positions
    Eigen::MatrixXd sensor_pos;
    sensor_pos.resize(M,3);

    //sensor transfer function
    std::vector<double> tf;

    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){


      //shift the positions to Hall element
      sensor_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      sensor_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      sensor_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);


      //Assemble Bn part
      sensors_[i].compute_H(sensor_pos, H_mat,M*i);


      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();


      //apply tranfer function of this sensor
      H_mat.block(i*M,0,M,N) *= tf[1];

    }

    //make a discrete potential
    DiscretePotential<Pot,LinOp> disc_pot;
    disc_pot.init_DiscretePotential(ansatz_space_);

    //impose a zero mean gauge gently
    //!!!!!!!! use block(0,1,M,N-1) lateron  !!!!!!!!
    disc_pot.impose_zero_mean_gauge_low_ram(H_mat);
  }

  void compute_Q(const Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> &llt,
                            Eigen::MatrixXd &Q){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    Q.resize(M,N-1);
    Q = llt.solve(H_mat.block(0,1,M,N-1));

    return;

  }


  void compute_A(const Eigen::MatrixXd &Q,
                  Eigen::MatrixXd &A){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    A.resize(N-1,N-1);

    A = H_mat.block(0,1,M,N-1).transpose()*Q;

    return;

  }

  Eigen::MatrixXd compute_A(const Eigen::SparseMatrix<double> &R_inv){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_mat.block(0,1,M,N-1).transpose()*R_inv*H_mat.block(0,1,M,N-1);

  }

  Eigen::MatrixXd compute_A_perturbed(const Eigen::SparseMatrix<double> &R_inv){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_tilde_.block(0,1,M,N-1).transpose()*R_inv*H_tilde_.block(0,1,M,N-1);

  }


  Eigen::MatrixXd compute_A(const Eigen::SparseMatrix<double> &R_inv,
                            const Eigen::VectorXd &pert){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    //Perturbed observation matrix
    Eigen::MatrixXd H_tilde(M,N);

    compute_H_tilde(pert,H_tilde);


    return H_tilde.block(0,1,M,N-1).transpose()*R_inv*H_tilde.block(0,1,M,N-1);

  }

  Eigen::MatrixXd compute_rhs(const Eigen::SparseMatrix<double> &R_inv,
                              const Eigen::VectorXd &y){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_mat.block(0,1,M,N-1).transpose() * R_inv * y;

  }

  Eigen::MatrixXd compute_rhs_perturbed(const Eigen::SparseMatrix<double> &R_inv,
                              const Eigen::VectorXd &y){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_tilde_.block(0,1,M,N-1).transpose() * R_inv * y;

  }

  Eigen::MatrixXd compute_rhs_perturbed(const Eigen::SparseMatrix<double> &R_inv,
                              const Eigen::VectorXd &y,
                              const Eigen::VectorXd &L_eps){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_tilde_.block(0,1,M,N-1).transpose() * (R_inv * y + L_eps);

  }

  Eigen::MatrixXd compute_rhs(const Eigen::SparseMatrix<double> &R_inv,
                              const Eigen::VectorXd &y,
                              const Eigen::VectorXd &L_eps){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_mat.block(0,1,M,N-1).transpose() * (R_inv * y + L_eps);

  }

  Eigen::VectorXd compute_rhs(const Eigen::VectorXd &u,
                              const Eigen::VectorXd &L_eps){

    //number of measuremements
    int M = H_mat.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    return H_mat.block(0,1,M,N-1).transpose() * (u + L_eps);

  }

  void read_arm_transfer_function_from_csv(const std::string filename,
                                            Eigen::VectorXcd &Tw,
                                            Eigen::VectorXcd &Tt,
                                            Eigen::VectorXd &omegas){

    //read the file
    Eigen::MatrixXd in_data = read_csv(filename,true);

    //how many frequencies are there?
    int num_freqs = in_data.rows();

    //make space for vectors
    Tw.resize(num_freqs);
    Tt.resize(num_freqs);
    omegas.resize(num_freqs);

    //immaginary unit
    std::complex<double> I(0,1);

    //assign data
    Tw = in_data.col(1) + I*in_data.col(2);
    Tt = in_data.col(3) + I*in_data.col(4);
    omegas = 2.*BEMBEL_PI*in_data.col(0);

    return;

  }

  void append_arm_transfer_function(const int move_axis_id,
                                  const Eigen::VectorXcd Tw_x,
                                  const Eigen::VectorXcd Tw_y,
                                  const Eigen::VectorXcd Tw_z,
                                  const Eigen::VectorXcd Tt_x,
                                  const Eigen::VectorXcd Tt_y,
                                  const Eigen::VectorXcd Tt_z,
                                  const Eigen::VectorXd omega_x,
                                  const Eigen::VectorXd omega_y,
                                  const Eigen::VectorXd omega_z,
                                  const double f_s,
                                  const int max_samples_per_move){

    //number of frequencies
    int num_freq_x = omega_x.rows();
    int num_freq_y = omega_y.rows();
    int num_freq_z = omega_z.rows();


    //make space for new transfer matrices
    Eigen::MatrixXd tmp_Fw_x,tmp_Fw_y,tmp_Fw_z,tmp_Ft_x,tmp_Ft_y,tmp_Ft_z;

    //the position has one more for static offsets
    //every frequency has two DoFs, except for the DC offset, to encode phase and
    //amplitude
    tmp_Fw_x.resize(max_samples_per_move,2*num_freq_x + 1);
    tmp_Fw_y.resize(max_samples_per_move,2*num_freq_y + 1);
    tmp_Fw_z.resize(max_samples_per_move,2*num_freq_z + 1);


    tmp_Ft_x.resize(max_samples_per_move,2*num_freq_x + 1);
    tmp_Ft_y.resize(max_samples_per_move,2*num_freq_y + 1);
    tmp_Ft_z.resize(max_samples_per_move,2*num_freq_z + 1);


    //time
    double t = 0.;

    //time between two triggers
    double t_s = 1/f_s;

    //fill transfer matrices
    for(int i = 0; i < max_samples_per_move; ++i, t += t_s){

      //static displacements
      tmp_Fw_x(i,0) = 1.;
      tmp_Fw_y(i,0) = 1.;
      tmp_Fw_z(i,0) = 1.;

      tmp_Ft_x(i,0) = 1.;
      tmp_Ft_y(i,0) = 1.;
      tmp_Ft_z(i,0) = 1.;

      //x axis
      for(int j = 0; j < num_freq_x; ++j){


        tmp_Fw_x(i,2*j+1) = Tw_x(j).real()*std::cos(omega_x(j)*t)
                      - Tw_x(j).imag()*std::sin(omega_x(j)*t);

        tmp_Fw_x(i,2*j+2) = Tw_x(j).real()*std::cos(omega_x(j)*t)
                      + Tw_x(j).imag()*std::sin(omega_x(j)*t);

        tmp_Ft_x(i,2*j+1) = Tt_x(j).real()*std::cos(omega_x(j)*t)
                      - Tt_x(j).imag()*std::sin(omega_x(j)*t);

        tmp_Ft_x(i,2*j+2) = Tt_x(j).real()*std::cos(omega_x(j)*t)
                      + Tt_x(j).imag()*std::sin(omega_x(j)*t);

        }


      //y axis
      for(int j = 0; j < num_freq_y; ++j){

        tmp_Fw_y(i,2*j+1) = Tw_y(j).real()*std::cos(omega_y(j)*t)
                      - Tw_y(j).imag()*std::sin(omega_y(j)*t);

        tmp_Fw_y(i,2*j+2) = Tw_y(j).real()*std::cos(omega_y(j)*t)
                        + Tw_y(j).imag()*std::sin(omega_y(j)*t);

        tmp_Ft_y(i,2*j+1) = Tt_y(j).real()*std::cos(omega_y(j)*t)
                      - Tt_y(j).imag()*std::sin(omega_y(j)*t);

        tmp_Ft_y(i,2*j+2) = Tt_y(j).real()*std::cos(omega_y(j)*t)
                      + Tt_y(j).imag()*std::sin(omega_z(j)*t);

        }
      //z axis
      for(int j = 0; j < num_freq_z; ++j){


        tmp_Fw_z(i,2*j+1) = Tw_z(j).real()*std::cos(omega_z(j)*t)
                      - Tw_z(j).imag()*std::sin(omega_z(j)*t);

        tmp_Fw_z(i,2*j+2) = Tw_z(j).real()*std::cos(omega_z(j)*t)
                      + Tw_z(j).imag()*std::sin(omega_z(j)*t);

        tmp_Ft_z(i,2*j+1) = Tt_z(j).real()*std::cos(omega_z(j)*t)
                    - Tt_z(j).imag()*std::sin(omega_z(j)*t);

        tmp_Ft_z(i,2*j+2) = Tt_z(j).real()*std::cos(omega_z(j)*t)
                      + Tt_z(j).imag()*std::sin(omega_z(j)*t);

      }

    }

    //save move direction
    move_directions.push_back(move_axis_id);

    //save transfer matrices
    Fw_x.push_back(tmp_Fw_x);
    Fw_y.push_back(tmp_Fw_y);
    Fw_z.push_back(tmp_Fw_z);
    Ft_x.push_back(tmp_Ft_x);
    Ft_y.push_back(tmp_Ft_y);
    Ft_z.push_back(tmp_Ft_z);

  }


  //compute the measurements spatial derivatives
  void compute_dH_matrices(const Eigen::MatrixXd &pos){

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();

    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //make space for matrix
    dHx_mat.resize(num_probes*M,N), dHy_mat.resize(num_probes*M,N), dHz_mat.resize(num_probes*M,N);

    //evaluation positions
    Eigen::MatrixXd eval_pos;
    eval_pos.resize(M,3);

    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);
    Eigen::MatrixXd I_N = Eigen::MatrixXd::Ones(1,N);

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){

      //shift the positions to Hall element
      eval_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      eval_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      eval_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);

      sensors_[i].compute_measurement_derivative_matrix_low_ram(eval_pos, dHx_mat, dHy_mat, dHz_mat,i*M);

    }

    return;
  }

  Eigen::MatrixXd assemble_jacobian(const Eigen::MatrixXd &pos,
                                    const Eigen::VectorXd v_0,
                                    const bool use_small_matrix = true){

    //The Jacobian is given by:
    //   J(v_0) = ( sum_i i q_i (H v_0)^(i-1) x I_(1 x N) ) o H
    //where i iterates over the sensors in the array.
    // x is the Kronecker and o is the Hadamard product

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();
    if (enable_gauge_){
      N -= 1;
    }

    //note that we will gauge for zero mean arithmetically. So the Matrices
    //will be of size (M x N-1)

    //Bn matrix
    Eigen::MatrixXd H;
    //make space for Bn matrix contributions
    H.resize(M,N);

    //H.v0 matrix
    Eigen::VectorXd Hv0;
    Hv0.resize(M);

    //left side vector
    Eigen::VectorXd sHv0;
    sHv0.resize(M);


    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //make space for jacobian
    Eigen::MatrixXd J;
    J.resize(num_probes*M,N);
    J.setZero();


    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);
    Eigen::MatrixXd I_N = Eigen::MatrixXd::Ones(1,N);

    //evaluation positions
    Eigen::MatrixXd eval_pos;
    eval_pos.resize(M,3);

    //sensor transfer function
    std::vector<double> tf;


    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){


      //shift the positions to Hall element
      eval_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      eval_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      eval_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);



      //Assemble Bn matrix
      //#pragma omp parallel
      if(use_small_matrix)
      {
        sensors_[i].compute_measurement_matrix_low_ram(eval_pos,H,enable_gauge_);
      }
      else
      {
        sensors_[i].compute_measurement_matrix(eval_pos,H,enable_gauge_);
      }



      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();



      //Bn component (we could spare this multiplication for linear sensors...)
      if (enable_gauge_)  Hv0 = H * v_0.segment(1,N);
      else                Hv0 = H * v_0;

      sHv0.setZero();

      //we now iterate over the nonlinear terms
      for(int j = 1; j < tf.size() ; ++j){
        sHv0 += j*tf[j]*Hv0.array().pow(j-1).matrix();
      }


      J.block(i*M,0,M,N) = ((sHv0 * I_N.row(0)).array() * H.array()).matrix();

    }


    return J;

  }

  void assemble_spatial_derivatives(const Eigen::MatrixXd &pos,
                                    const Eigen::VectorXd v_0,
                                    Eigen::MatrixXd *Jx,
                                    Eigen::MatrixXd *Jy,
                                    Eigen::MatrixXd *Jz ) {

    //Same as the function above, but we compute the spatial derivatives of the
    //measurement operation. This is achieved by deriving the greens kernel.

    //number of measuremements
    int M = pos.rows();

    //number of DoFs
    int N = ansatz_space_.get_number_of_dofs();
    if (enable_gauge_){
      N -= 1;
    }

    //all the matrices below are 3 x M since we need to compute the three spatial
    //derivatives
    //Bn matrix
    Eigen::MatrixXd dH;
    //make space for Bn matrix contributions
    dH.resize(3*M,N);

    //dH.v0 vector
    Eigen::VectorXd dHv0;
    dHv0.resize(3*M);

    //left side vector
    Eigen::VectorXd sdHv0;
    sdHv0.resize(3*M);


    //number of Hall probes in array
    int num_probes =  sensors_.size();

    //make space for jacobians
    Jx->resize(num_probes*M,N);
    Jx->setZero();

    Jy->resize(num_probes*M,N);
    Jy->setZero();

    Jz->resize(num_probes*M,N);
    Jz->setZero();

    //Helper unit vectors
    Eigen::MatrixXd I_M = Eigen::MatrixXd::Ones(M,1);
    Eigen::MatrixXd I_N = Eigen::MatrixXd::Ones(1,N);

    //evaluation positions
    Eigen::MatrixXd eval_pos;
    eval_pos.resize(M,3);

    //sensor transfer function
    std::vector<double> tf;

    //iterate over the probes in the array
//#pragma omp parallel for
    for(int i = 0 ; i < num_probes ; ++i){

      //shift the positions to Hall element

      eval_pos.col(0) = pos.col(0) + positions_[i](0)*I_M.col(0);
      eval_pos.col(1) = pos.col(1) + positions_[i](1)*I_M.col(0);
      eval_pos.col(2) = pos.col(2) + positions_[i](2)*I_M.col(0);

      //Assemble dBn matrix
      //#pragma omp parallel
      {
        sensors_[i].compute_measurement_derivative_matrix(eval_pos,dH,enable_gauge_);
      }

      //tranfer function of this sensor
      tf = sensors_[i].get_sensor().get_transfer_function();

      //Bn component (we could spare this multiplication for linear sensors...)
      if (enable_gauge_)  dHv0 = dH * v_0.segment(1,N);
      else                dHv0 = dH * v_0;

      sdHv0.setZero();

      //THIS IS WRONG FOR NONLINEAR SENSORS CHECK THE EQUATIONS
      //we now iterate over the nonlinear terms
      for(int j = 1; j < tf.size() ; ++j){

        sdHv0 += j*tf[j]*dHv0.array().pow(j-1).matrix();
      }


      Jx->block(i*M,0,M,N) = ((sdHv0.segment(0,M)   * I_N.row(0)).array() * dH.block(0,0,M,N).array()).matrix();
      Jy->block(i*M,0,M,N) = ((sdHv0.segment(M,M)   * I_N.row(0)).array() * dH.block(M,0,M,N).array()).matrix();
      Jz->block(i*M,0,M,N) = ((sdHv0.segment(2*M,M) * I_N.row(0)).array() * dH.block(2*M,0,M,N).array()).matrix();


    }
    return;
  }

  Eigen::VectorXd get_offset_vector(int M){


    //numner of Hall probes in array
    int num_probes =  sensors_.size();

    //return vector
    Eigen::VectorXd ret_vec;
    ret_vec.resize(num_probes*M);

    //Helper unit vectors
    Eigen::VectorXd i_M(M);
    i_M.setOnes();

    //iterate over the probes in the array
    for(int i = 0 ; i < num_probes ; ++i){

      //tranfer function of this sensor
      std::vector<double> tf = sensors_[i].get_sensor().get_transfer_function();

      ret_vec.segment(i*M,M) = tf[0]*i_M;

    }

    return ret_vec;

  }

  void init_aca_matrix(MeasurementData *meas,int refine_lvl,double acc = 1e-4,double dist = 1.6){


    aca_matrices_.resize(sensors_.size());


    num_eval_ = meas->get_number_of_measurements();

    meas->init_cluster_tree(refine_lvl);



    for (int i = 0 ; i < sensors_.size() ; ++i){

      aca_matrices_[i].set_parameters(dist,refine_lvl, 1, acc);

      meas->translate_positions(positions_[i]);

      aca_matrices_[i].init_ACAMatrix(sensors_[i],meas,&ansatz_space_);

      meas->translate_positions(-1*positions_[i]);


    }

    if(h_.size() == 0){

      h_ = make_h_vector();
    }

  }



  Eigen::VectorXd mat_vec(const Eigen::VectorXd &rhs){


    Eigen::VectorXd ret_val,tmp_res,tmp_rhs;
    ret_val.resize(num_eval_*sensors_.size());

    int num_DoFs = rhs.size();


    if(enable_gauge_){

      tmp_rhs.resize(num_DoFs+1);
      tmp_rhs.setZero();
      tmp_rhs.segment(1,num_DoFs) = rhs;


      double sv = (stability_vec_.segment(1,num_DoFs).transpose() * rhs)(0);
      tmp_rhs(0) = -1*sv/stability_vec_(0);

    }
    else{

      tmp_rhs.resize(num_DoFs);
      tmp_rhs = rhs;
    }


    for (int i = 0 ; i < sensors_.size() ; ++i){

      aca_matrices_[i].mat_vec_prod(&tmp_res,tmp_rhs);


      ret_val.segment(i*num_eval_,num_eval_) = tmp_res;

    }

    return ret_val;

  }

  Eigen::VectorXd mat_vec_transpose(const Eigen::VectorXd &rhs){


    Eigen::VectorXd ret_val,tmp_res;

    int num_DoFs = ansatz_space_.get_number_of_dofs();

    //tmp_res.resize(num_DoFs);
    //tmp_res.setZero();

    ret_val.resize(num_DoFs);
    ret_val.setZero();


    for (int i = 0 ; i < sensors_.size() ; ++i){
      aca_matrices_[i].mat_vec_transpose(&ret_val,rhs,i*num_eval_);

    }


    if(enable_gauge_){


      //double sv = stability_vec_.transpose() * ret_val;
      //double delta = sv / ret_val.size(); //stability_vec_.array().sum();

      //std::cout << "h = " << h_.rows() << " x " << h_.cols() << std::endl;
      //std::cout << "rhs = " << rhs.rows() << " x " << rhs.cols() << std::endl;
      //std::cout << "stability_vec_.segment(1,num_DoFs-1) = " << stability_vec_.segment(1,num_DoFs-1).rows() << " x " << stability_vec_.segment(1,num_DoFs-1).cols() << std::endl;

      double c = h_.transpose()*rhs.eval();

      Eigen::VectorXd correction = c * stability_vec_.segment(1,num_DoFs-1);

      //std::cout << "correction = " << correction.rows() << " x " << correction.cols() << std::endl;
      //std::cout << "ret_val.segment(1,num_DoFs-1) = " << ret_val.segment(1,num_DoFs-1).rows() << " x " << ret_val.segment(1,num_DoFs-1).cols() << std::endl;
      //ret_val -= sv*stability_vec_;

      ret_val = (ret_val.segment(1,num_DoFs-1) - correction).eval();//- delta*Eigen::VectorXd::Ones(num_DoFs-1)).eval() ;

    }

    return ret_val;




  }

  Eigen::VectorXd compute_voltages(const Eigen::VectorXd &rhs){

    Eigen::VectorXd ret_val,tmp_res,tmp_rhs;
    ret_val.resize(num_eval_*sensors_.size());

    int num_DoFs = rhs.size();

    if(enable_gauge_){

      tmp_rhs.resize(num_DoFs+1);
      tmp_rhs.setZero();
      tmp_rhs.segment(1,num_DoFs) = rhs;

      double sv = (stability_vec_.segment(1,num_DoFs).transpose() * rhs)(0);
      tmp_rhs(0) = -1*sv/stability_vec_(0);

    }
    else{

      tmp_rhs.resize(num_DoFs);
      tmp_rhs = rhs;
    }

    for (int i = 0 ; i < sensors_.size() ; ++i){

      //n.B
      aca_matrices_[i].mat_vec_prod(&tmp_res,tmp_rhs);

      //nonlinearity
      tmp_res = sensors_[i].get_sensor().apply_transfer_function(tmp_res);

      ret_val.segment(i*num_eval_,num_eval_) = tmp_res;

    }
    return ret_val;

  }


  DiscreteSensor<HallProbe<Pot,LinOp>,LinOp> &get_sensor(int i){

    return sensors_[i];
  }

  int get_number_of_sensors(){

    return sensors_.size();
  }

  Eigen::SparseMatrix<double> get_R_inv(const int num_meas){

    int number_of_sensors  = get_number_of_sensors();

    Eigen::SparseMatrix<double> R_inv(number_of_sensors*num_meas,number_of_sensors*num_meas);
    R_inv.setIdentity();
    for(int i = 0 ; i < number_of_sensors; ++i){
      R_inv.block(i*num_meas,i*num_meas,num_meas,num_meas) /= noise_stdev_[i]*noise_stdev_[i];

    }
    return R_inv;
  }

  Eigen::SparseMatrix<double> get_R(const int num_meas){

    int number_of_sensors  = get_number_of_sensors();

    Eigen::SparseMatrix<double> R(number_of_sensors*num_meas,number_of_sensors*num_meas);
    R.setIdentity();
    for(int i = 0 ; i < number_of_sensors; ++i){
      R.block(i*num_meas,i*num_meas,num_meas,num_meas) *= noise_stdev_[i]*noise_stdev_[i];

    }
    return R;
  }

  ACAMatrixRed<HallProbe<Pot,LinOp>,LinOp> *get_aca_matrix(int sensor_id){

    return &aca_matrices_[sensor_id];
  }


  void rotate(double Alpha,double Beta,double Gamma){

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


      for(int i = 0; i < sensors_.size() ; ++i){
        sensors_[i].get_sensor().rotate_probe(Alpha,Beta,Gamma);
        positions_[i] = Rot_X*Rot_Y*Rot_Z*positions_[i];
      }
  }

  void reset_offset_voltages(){

    for(int i = 0; i < sensors_.size() ; ++i){
      sensors_[i].get_sensor().reset_offset();
    }
  }

  void set_gauged_formulation(bool enable){


    if( enable ){



      DiscretePotential<Pot,LinOp> disc_pot;
      disc_pot.init_DiscretePotential(ansatz_space_);

      stability_vec_ = disc_pot.compute_stability_vector();


      if((h_.size() == 0) && aca_matrices_.size() > 0){

        //compute H_{:,1}
        h_ = make_h_vector();

      }

      enable_gauge_ = enable;

    }
    /*
    else{

      std::cout << "Tried to setup gauge condition for tranpose mapping, but no matrices were found. Pleas initialize the matrices before enabling gauging!" << std::endl;

    }
    */

  }

  Eigen::VectorXd recover_full_solution(const Eigen::VectorXd v_in){

    //number of unknowns
    int N = v_in.size()+1;

    //std::cout << "N = " << N << std::endl;

    Eigen::VectorXd v_new;

    v_new.resize(N);

    v_new.segment(1,N-1) = v_in;

    //std::cout << "stability_vec_ = " << stability_vec_.size() << std::endl;

    v_new(0) = -1*(stability_vec_.segment(1,N-1).transpose() * v_in)(0)/stability_vec_(0);

    return v_new;
  }

  Eigen::MatrixXd recover_full_solution(const Eigen::MatrixXd V_in){

    //number of unknowns
    int N = V_in.rows()+1;

    //number of ensembles
    int E = V_in.cols();


    Eigen::MatrixXd V_new;

    V_new.resize(N,E);

    V_new.block(1,0,N-1,E) = V_in;


    for(int i = 0 ; i < E ; ++i){

      V_new(0,i) = -1.*(stability_vec_.segment(1,N-1).transpose() * V_in.col(i))(0)/stability_vec_(0);
    }


    return V_new;
  }

  int get_compressed_number_of_elements(){

    int ret_val = 0;

    for(int i = 0; i < sensors_.size() ; ++i){
      ret_val += aca_matrices_[i].get_compressed_number_of_elements();
    }
    return ret_val;
  }

  int get_number_of_dofs(){

    return  ansatz_space_.get_number_of_dofs();
  }

  int get_number_of_dofs_disc(){

    return  ansatz_space_.get_transformation_matrix().rows();
  }

  int get_number_of_measurements(){

    return  num_eval_;
  }

  Eigen::VectorXd get_stability_vector(){

    if(enable_gauge_ == false){

      DiscretePotential<Pot,LinOp> disc_pot;

      disc_pot.init_DiscretePotential(ansatz_space_);

      stability_vec_ = disc_pot.compute_stability_vector();
    }

      return stability_vec_;
  }

  bool gauged_formulation_enabled(){

    return enable_gauge_;
  }

  Eigen::VectorXd make_h_vector(){

    //compute H_{:,1}
    int num_DoFs = ansatz_space_.get_number_of_dofs();

    Eigen::VectorXd test_vec(num_DoFs);
    test_vec.setZero();
    test_vec(0) = 1.;

    bool switch_gauge_flag = false;

    if(enable_gauge_){
      enable_gauge_ = false;
    }
    else{
      get_stability_vector();
    }


    Eigen::VectorXd ret_vec = mat_vec(test_vec)/stability_vec_(0);

    if(switch_gauge_flag) enable_gauge_ = true;

    return ret_vec;

  }

  AnsatzSpace<LinOp> *get_ansatz_space_ptr(){

    return &ansatz_space_;
  }

  Eigen::VectorXd compute_dx(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();

    return Fw_x[move_id]*a.segment(0,num_x);

  }
  Eigen::VectorXd compute_dy(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();
    int num_y = Fw_y[move_id].cols();

    return Fw_y[move_id]*a.segment(num_x,num_y);

  }
  Eigen::VectorXd compute_dz(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();
    int num_y = Fw_y[move_id].cols();
    int num_z = Fw_z[move_id].cols();

    return Fw_z[move_id]*a.segment(num_x+num_y,num_z);

  }

  Eigen::VectorXd compute_tx(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();
    int num_y = Fw_y[move_id].cols();
    int num_z = Fw_z[move_id].cols();

    return Ft_x[move_id]*a.segment(num_x + num_y + num_z,num_x);

  }
  Eigen::VectorXd compute_ty(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();
    int num_y = Fw_y[move_id].cols();
    int num_z = Fw_z[move_id].cols();

    return Ft_y[move_id]*a.segment(2*num_x + num_y + num_z , num_y);

  }
  Eigen::VectorXd compute_tz(const int move_id, const Eigen::VectorXd &a){

    int num_x = Fw_x[move_id].cols();
    int num_y = Fw_y[move_id].cols();
    int num_z = Fw_z[move_id].cols();


    return Ft_z[move_id]*a.segment(2*num_x + 2*num_y + num_z, num_z);

  }

  Eigen::MatrixXd compute_P_matrix(const int move_id, const Eigen::MatrixXd &dH){

    //number of measurements
    int num_meas = Fw_x[move_id].rows();
    //int num_meas = dH.rows();

    //number of DoFs for motion
    int num_ax = Fw_x[move_id].cols();
    int num_ay = Fw_y[move_id].cols();
    int num_az = Fw_z[move_id].cols();


    //helper diagonal matrix
    Eigen::DiagonalMatrix<double, -1> dHx_x(dH.block(0,0,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHx_y(dH.block(num_meas,0,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHx_z(dH.block(2*num_meas,0,num_meas,1));

    Eigen::DiagonalMatrix<double, -1> dHy_x(dH.block(0,1,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHy_y(dH.block(num_meas,1,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHy_z(dH.block(2*num_meas,1,num_meas,1));

    Eigen::DiagonalMatrix<double, -1> dHz_x(dH.block(0,2,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHz_y(dH.block(num_meas,2,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHz_z(dH.block(2*num_meas,2,num_meas,1));

    Eigen::DiagonalMatrix<double, -1> dHt1_x(dH.block(0,3,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHt1_y(dH.block(num_meas,3,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHt1_z(dH.block(2*num_meas,3,num_meas,1));

    Eigen::DiagonalMatrix<double, -1> dHt2_x(dH.block(0,4,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHt2_y(dH.block(num_meas,4,num_meas,1));
    Eigen::DiagonalMatrix<double, -1> dHt2_z(dH.block(2*num_meas,4,num_meas,1));

    /*
    DISCRETE PARAMETERS
    Eigen::MatrixXd P(3*num_meas , 5*num_meas);

    P.block(0,0,num_meas,num_meas) = dHx_x;
    P.block(num_meas,0,num_meas,num_meas) = dHx_y;
    P.block(2*num_meas,0,num_meas,num_meas) = dHx_z;

    P.block(0,num_meas,num_meas,num_meas) = dHy_x;
    P.block(num_meas,num_meas,num_meas,num_meas) = dHy_y;
    P.block(2*num_meas,num_meas,num_meas,num_meas) = dHy_z;

    P.block(0,2*num_meas,num_meas,num_meas) = dHz_x;
    P.block(num_meas,2*num_meas,num_meas,num_meas) = dHz_y;
    P.block(2*num_meas,2*num_meas,num_meas,num_meas) = dHz_z;

    P.block(0,3*num_meas,num_meas,num_meas) = dHt1_x;
    P.block(num_meas,3*num_meas,num_meas,num_meas) = dHt1_y;
    P.block(2*num_meas,3*num_meas,num_meas,num_meas) = dHt1_z;

    P.block(0,4*num_meas,num_meas,num_meas) = dHt2_x;
    P.block(num_meas,4*num_meas,num_meas,num_meas) = dHt2_y;
    P.block(2*num_meas,4*num_meas,num_meas,num_meas) = dHt2_z;
    */

    /*
    W AND T INDEPENDENT
    */
    Eigen::MatrixXd P(3*num_meas , 2*num_ax + 2*num_ay + 2*num_az );
    P.setZero();

    P.block(0,0,num_meas,num_ax)          = dHx_x*Fw_x[move_id];
    P.block(num_meas,0,num_meas,num_ax)   = dHx_y*Fw_x[move_id];
    P.block(2*num_meas,0,num_meas,num_ax) = dHx_z*Fw_x[move_id];

    P.block(0,num_ax,num_meas,num_ay)          = dHy_x*Fw_y[move_id];
    P.block(num_meas,num_ax,num_meas,num_ay)   = dHy_y*Fw_y[move_id];
    P.block(2*num_meas,num_ax,num_meas,num_ay) = dHy_z*Fw_y[move_id];

    P.block(0,num_ax+num_ay,num_meas,num_az)          = dHz_x*Fw_z[move_id];
    P.block(num_meas,num_ax+num_ay,num_meas,num_az)   = dHz_y*Fw_z[move_id];
    P.block(2*num_meas,num_ax+num_ay,num_meas,num_az) = dHz_z*Fw_z[move_id];

    int col_begin = num_ax+num_ay+num_az;

    if(arm_axis_id_ == 0){
      //std::cout << "x axis" << std::endl;

      //theta 1 is in xy plane
      P.block(0,col_begin,num_meas,num_ay)          = dHt1_x*Ft_y[move_id];
      P.block(num_meas,col_begin,num_meas,num_ay)   = dHt1_y*Ft_y[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_ay) = dHt1_z*Ft_y[move_id];


      col_begin +=  num_ay;


      //theta 2 is in xz plane
      P.block(0,col_begin,num_meas,num_az)          = dHt2_x*Ft_z[move_id];
      P.block(num_meas,col_begin,num_meas,num_az)   = dHt2_y*Ft_z[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_az) = dHt2_z*Ft_z[move_id];

    }
    else if(arm_axis_id_ == 1){
      //std::cout << "y axis" << std::endl;
      //theta 1 is in yz plane
      P.block(0,col_begin,num_meas,num_az)          = dHt1_x*Ft_z[move_id];
      P.block(num_meas,col_begin,num_meas,num_az)   = dHt1_y*Ft_z[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_az) = dHt1_z*Ft_z[move_id];

      col_begin +=  num_az;

      //theta 2 is in xy plane
      P.block(0,col_begin,num_meas,num_ax)          = dHt2_x*Ft_x[move_id];
      P.block(num_meas,col_begin,num_meas,num_ax)   = dHt2_y*Ft_x[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_ax) = dHt2_z*Ft_x[move_id];
    }
    else{
      //std::cout << "z axis" << std::endl;
      //theta 1 is in yz plane
      P.block(0,col_begin,num_meas,num_ay)          = dHt1_x*Ft_y[move_id];
      P.block(num_meas,col_begin,num_meas,num_ay)   = dHt1_y*Ft_y[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_ay) = dHt1_z*Ft_y[move_id];

      col_begin +=  num_ay;

      //theta 2 is in xz plane
      P.block(0,col_begin,num_meas,num_ax)          = dHt2_x*Ft_x[move_id];
      P.block(num_meas,col_begin,num_meas,num_ax)   = dHt2_y*Ft_x[move_id];
      P.block(2*num_meas,col_begin,num_meas,num_ax) = dHt2_z*Ft_x[move_id];
    }



    /*
    W AND T LINKED WITH TF

    Eigen::MatrixXd P(3*num_meas , num_ax + num_ay + num_az);
    P.setZero();

    P.block(0,0,num_meas,num_ax)          = dHx_x*Fw_x[move_id];
    P.block(num_meas,0,num_meas,num_ax)   = dHx_y*Fw_x[move_id];
    P.block(2*num_meas,0,num_meas,num_ax) = dHx_z*Fw_x[move_id];

    P.block(0,num_ax,num_meas,num_ay)          = dHy_x*Fw_y[move_id];
    P.block(num_meas,num_ax,num_meas,num_ay)   = dHy_y*Fw_y[move_id];
    P.block(2*num_meas,num_ax,num_meas,num_ay) = dHy_z*Fw_y[move_id];

    P.block(0,num_ax+num_ay,num_meas,num_az)          = dHz_x*Fw_z[move_id];
    P.block(num_meas,num_ax+num_ay,num_meas,num_az)   = dHz_y*Fw_z[move_id];
    P.block(2*num_meas,num_ax+num_ay,num_meas,num_az) = dHz_z*Fw_z[move_id];


        if(arm_axis_id_ == 0){
          std::cout << "x axis" << std::endl;
          //theta 1 is in xy plane
          P.block(0,num_ax+1,num_meas,num_ay-1)          += dHt1_x*Ft_y[move_id];
          P.block(num_meas,num_ax+1,num_meas,num_ay-1)   = dHt1_y*Ft_y[move_id];
          P.block(2*num_meas,num_ax+1,num_meas,num_ay-1) = dHt1_z*Ft_y[move_id];

          //theta 2 is in xz plane
          P.block(0,num_ax+num_ay+1,num_meas,num_az-1)          += dHt2_x*Ft_z[move_id];
          P.block(num_meas,num_ax+num_ay+1,num_meas,num_az-1)   += dHt2_y*Ft_z[move_id];
          P.block(2*num_meas,num_ax+num_ay+1,num_meas,num_az-1) += dHt2_z*Ft_z[move_id];

        }
        else if(arm_axis_id_ == 1){
          std::cout << "y axis" << std::endl;
          //theta 1 is in yz plane
          P.block(0,num_ax+num_ay+1,num_meas,num_az-1)          += dHt1_x*Ft_z[move_id];
          P.block(num_meas,num_ax+num_ay+1,num_meas,num_az-1)   += dHt1_y*Ft_z[move_id];
          P.block(2*num_meas,num_ax+num_ay+1,num_meas,num_az-1) += dHt1_z*Ft_z[move_id];

          //theta 2 is in xy plane
          P.block(0,1,num_meas,num_ax-1)          += dHt2_x*Ft_x[move_id];
          P.block(num_meas,1,num_meas,num_ax-1)   += dHt2_y*Ft_x[move_id];
          P.block(2*num_meas,1,num_meas,num_ax-1) += dHt2_z*Ft_x[move_id];
        }
        else{
          std::cout << "z axis" << std::endl;
          //theta 1 is in yz plane
          P.block(0,num_ax+1,num_meas,num_ay-1)          += dHt1_x*Ft_y[move_id];
          P.block(num_meas,num_ax+1,num_meas,num_ay-1)   += dHt1_y*Ft_y[move_id];
          P.block(2*num_meas,num_ax+1,num_meas,num_ay-1) += dHt1_z*Ft_y[move_id];

          //theta 2 is in xz plane
          P.block(0,1,num_meas,num_ax-1)          += dHt2_x*Ft_x[move_id];
          P.block(num_meas,1,num_meas,num_ax-1)   += dHt2_y*Ft_x[move_id];
          P.block(2*num_meas,1,num_meas,num_ax-1) += dHt2_z*Ft_x[move_id];
        }
    */

    return P;

  }

  void save_UQ_matrices(const std::string out_dir){

    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream out_file_1(out_dir + "H_mat.csv");
    out_file_1 << H_mat.format(CSVFormat);
    out_file_1.close();

    std::ofstream out_file_2(out_dir + "dtH_1.csv");
    out_file_2 << dtH_1.format(CSVFormat);
    out_file_2.close();

    std::ofstream out_file_3(out_dir + "dtH_2.csv");
    out_file_3 << dtH_2.format(CSVFormat);
    out_file_3.close();

    std::ofstream out_file_4(out_dir + "dHx_mat.csv");
    out_file_4 << dHx_mat.format(CSVFormat);
    out_file_4.close();

    std::ofstream out_file_5(out_dir + "dHy_mat.csv");
    out_file_5 << dHy_mat.format(CSVFormat);
    out_file_5.close();

    std::ofstream out_file_6(out_dir + "dHz_mat.csv");
    out_file_6 << dHz_mat.format(CSVFormat);
    out_file_6.close();
  }

  void read_UQ_matrices(const std::string out_dir){

    H_mat = read_csv(out_dir + "H_mat.csv",false);
    dtH_1 = read_csv(out_dir + "dtH_1.csv",false);
    dtH_2 = read_csv(out_dir + "dtH_2.csv",false);
    dHx_mat = read_csv(out_dir + "dHx_mat.csv",false);
    dHy_mat = read_csv(out_dir + "dHy_mat.csv",false);
    dHz_mat = read_csv(out_dir + "dHz_mat.csv",false);


  }

  void set_arm_axis_orientation(const int axis_id){

    arm_axis_id_ = axis_id;

    int num_probes =  sensors_.size();

    for(int i = 0 ; i < num_probes ; ++i){
      sensors_[i].get_sensor().setup_transversal_vectors(axis_id);
    }

  }

  void append_move_property(const int move_id, std::string move_cov_file){

    move_ids_.push_back(move_id);

    if(arm_axis_id_ == -1){

      std::cout << " WARNING arm axis is not defined!" << std::endl;
    }

    PerturbationCovariance this_pert_cov(move_cov_file,arm_axis_id_);

    //std::cout << "move_id = " << move_id << std::endl;
    //std::cout << "move_cov_file = " << move_cov_file << std::endl;

    Pert_Covs_.push_back(this_pert_cov);



    return;
  }

  /*
  void setup_move_properties(const std::vector<int> move_ids,
                            const std::vector<std::string> move_covs){

    //number of moves
    int num_moves = move_covs.size();
    //make space for perturbation covariance matrices
    Pert_Covs_.resize(num_moves);


    for(int i = 0; i < num_moves; ++i){

      Pert_Covs_[i].read_from_file(move_covs[i], arm_axis_id_);

      //Pert_Covs_[i].save_to_file("../../../../../Measurements/vibrations_3D_mapper/august_2021/check_" + std::to_string(i) + ".csv");
    }

    move_ids_ = move_ids;
    return;
  }
  */

  int get_move_index(const int move_id){

    for(int i = 0 ; i < move_ids_.size(); ++i){

      if(move_ids_[i] == move_id) return i;

    }
    return -1;
  }

  /*
  void compute_sparse_perturbation_matrix(const Eigen::VectorXd &v){

    //spatial derivatives
    Eigen::MatrixXd dH_dr = evaluate_grad_H(v);
    //angular derivatives
    Eigen::MatrixXd dH_dt = evaluate_dH_dt(v);

    int M_large = dH_dr.rows();
    int M_small = M_large/3;

    // filling the matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(5*M_large);


    for(int m = 0; m < M_small ; ++m)
    {

      for(int i = 0; i < 3 ; ++i){

        tripletList.push_back(T(m+i*M_small,m,dH_dr(m+i*M_small,0)));
        tripletList.push_back(T(m+i*M_small,m+M_small,dH_dr(m+i*M_small,1)));
        tripletList.push_back(T(m+i*M_small,m+2*M_small,dH_dr(m+i*M_small,2)));

        tripletList.push_back(T(m+i*M_small,m+3*M_small,dH_dt(m+i*M_small,0)));
        tripletList.push_back(T(m+i*M_small,m+4*M_small,dH_dt(m+i*M_small,1)));

      }

    }
    dH_.resize(M_large,5*M_small);
    dH_.setFromTriplets(tripletList.begin(), tripletList.end());

  }
  */
  Eigen::SparseMatrix<double> compute_sparse_perturbation_matrix(const Eigen::VectorXd &v){

    //spatial derivatives
    Eigen::MatrixXd dH_dr = evaluate_grad_H(v);
    //angular derivatives
    Eigen::MatrixXd dH_dt = evaluate_dH_dt(v);

    int M_large = dH_dr.rows();
    int M_small = M_large/3;

    // filling the matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(5*M_large);


    for(int m = 0; m < M_small ; ++m)
    {

      for(int i = 0; i < 3 ; ++i){

        tripletList.push_back(T(m+i*M_small,m,dH_dr(m+i*M_small,0)));
        tripletList.push_back(T(m+i*M_small,m+M_small,dH_dr(m+i*M_small,1)));
        tripletList.push_back(T(m+i*M_small,m+2*M_small,dH_dr(m+i*M_small,2)));

        tripletList.push_back(T(m+i*M_small,m+3*M_small,dH_dt(m+i*M_small,0)));
        tripletList.push_back(T(m+i*M_small,m+4*M_small,dH_dt(m+i*M_small,1)));

      }

    }
    //dH_.resize(M_large,5*M_small);
    //dH_.setFromTriplets(tripletList.begin(), tripletList.end());

    //return dH_;

    Eigen::SparseMatrix<double> dH;

    dH.resize(M_large,5*M_small);
    dH.setFromTriplets(tripletList.begin(), tripletList.end());

    return dH;
  }

  Eigen::SparseMatrix<double> compute_sparse_perturbation_matrix(const int index_from,
                                        const int num_samples,
                                        const Eigen::MatrixXd &dH_dr,
                                        const Eigen::MatrixXd &dH_dt){

    // res[k] = dH_dr[k,0] * dw_x + dH_dr[k,1] * dw_y + dH_dr[k,2] * dw_z + dH_dt[k,0] * d_t1 + dH_dt[k,1] * d_t2
    // so the matrix involved is of size (num_probes*num_samples) x 5*num_samples with num_probes*5*num_samples entries

    //number of probes
    int num_probes = sensors_.size();

    //number of measurements in total
    int num_meas = dH_dr.rows()/num_probes;

    // filling the matrix
    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(5*num_probes*num_samples);


    for(int m = 0; m < num_samples ; ++m)
    {

      for(int i = 0; i < num_probes ; ++i){//

        tripletList.push_back(T(m+i*num_samples , m              ,dH_dr(m+index_from+i*num_meas , 0)));
        tripletList.push_back(T(m+i*num_samples , m+num_samples  ,dH_dr(m+index_from+i*num_meas , 1)));
        tripletList.push_back(T(m+i*num_samples , m+2*num_samples,dH_dr(m+index_from+i*num_meas , 2)));

        tripletList.push_back(T(m+i*num_samples , m+3*num_samples,dH_dt(m+index_from+i*num_meas  , 0)));
        tripletList.push_back(T(m+i*num_samples , m+4*num_samples,dH_dt(m+index_from+i*num_meas  , 1)));

      }

    }
    Eigen::SparseMatrix<double> dH;

    dH.resize(num_probes*num_samples,5*num_samples);
    dH.setFromTriplets(tripletList.begin(), tripletList.end());

    return dH;
  }

  Eigen::VectorXd compute_MAP_perturbations(const int index_from,
                                        const int num_samples,
                                        const int move_id,
                                        const Eigen::VectorXd residuals,
                                        const Eigen::MatrixXd &dH_dr,
                                        const Eigen::MatrixXd &dH_dt,
                                        const double delta){

    //number of probes
    int num_probes = sensors_.size();

    //number of total measurements
    int num_meas = residuals.size()/num_probes;

    //position of this move in memory
    int move_pos = get_move_index(move_id);


    //std::cout << "compute_sparse_perturbation_matrix" << std::endl;
    //compute derivative matrix
    Eigen::SparseMatrix<double> dH = compute_sparse_perturbation_matrix(index_from,
                                                                        num_samples,
                                                                        dH_dr,
                                                                        dH_dt);



    Eigen::SparseMatrix<double> R_inv = get_R_inv(num_samples);
    //with this we can compute the A matrix
    //we need to apply the sensor precision on the diagonal
    Eigen::MatrixXd A = dH.transpose()*R_inv*dH;


    //Eigen::MatrixXd Cov = Pert_Covs_[move_pos].get_cov(num_samples);
    //std::ofstream out_file_Cov("~/cernbox/Measurements/vibrations_3D_mapper/august_2021/verify_recovery_algorithm/output/Cov.csv");
    //out_file_Cov << Cov.format(CSVFormat);
    //out_file_Cov.close();

    //std::cout << "compute the prior for regularisation"  << std::endl;
    //we now compute the prior for regularisation
    Eigen::MatrixXd L = delta*Pert_Covs_[move_pos].get_prec(num_samples);

    //and now we construct the right hand side
    Eigen::VectorXd rhs(num_probes*num_samples);
    for(int i = 0; i < num_probes; ++i){
      rhs.segment(i*num_samples,num_samples) = residuals.segment(index_from+i*num_meas,num_samples);
    }
    rhs = dH.transpose()*R_inv*rhs;



    /*
    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream out_file_A("../../../../../Measurements/synthetic/cylinder/A.csv");
    out_file_A << A.format(CSVFormat);
    out_file_A.close();

    std::ofstream out_file_rhs("../../../../../Measurements/synthetic/cylinder/rhs.csv");
    out_file_rhs << A.format(CSVFormat);
    out_file_rhs.close();

    std::ofstream out_file_L("../../../../../Measurements/synthetic/cylinder/L.csv");
    out_file_L << L.format(CSVFormat);
    out_file_L.close();
    */

    //finally we solve for the perturbations
    Eigen::VectorXd d = (A + L).lu().solve(rhs);

    return d;
  }


  Eigen::VectorXd sample_from_prior_perturbation(const int num_samples,const int move_id){

    //position of this move in memory
    int move_pos = get_move_index(move_id);

    Eigen::MatrixXd cov = Pert_Covs_[move_pos].get_cov(num_samples).eval();
    Eigen::VectorXd p_mean(cov.rows());
    p_mean.setZero();

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigenSolver(cov);
    Eigen::MatrixXd normTransform = eigenSolver.eigenvectors()
                   * eigenSolver.eigenvalues().cwiseSqrt().asDiagonal();


    Eigen::VectorXd sample = generate_gaussian_samples(1, p_mean, cov, std::rand()).col(0);//


    return sample;
  }

  void setup_precision_matrices(const std::vector<int> move_separator,
                                      const std::vector<int> move_ids){

    //clear all precision matrices in memory
    prec_list_.clear();

    //clear move prec table
    move_prec_table_.clear();


    for(int m = 0; m < move_separator.size() - 1; ++m){

      //number of positions in this move
      int num_positions = move_separator[m+1] - move_separator[m];

      //find out if this matrix was computed already
      int pos_in_memory = get_position_in_move_prec_table(num_positions,move_ids[m]);

      if(pos_in_memory == -1){

        //we must make it

        //position of this move in memory
        int move_pos = get_move_index(move_ids[m]);

        //position of this move in memory
        prec_list_.push_back(Pert_Covs_[move_pos].get_prec(num_positions));

        move_prec_table_.push_back({ num_positions , move_ids[m] });

        Pert_Covs_[move_pos].get_cov(num_positions);

      }

    }

    return;


  }

  Eigen::MatrixXd get_cov(const int move_pos,const int num_positions){

      std::cout << "move pos = " << move_pos << std::endl;

      return Pert_Covs_[move_pos].get_cov(num_positions);

  }

  int get_position_in_move_prec_table(const int num_moves, const int move_id){


    for(int i = 0; i < move_prec_table_.size(); ++i){

      if((move_prec_table_[i][0] == num_moves) && (move_prec_table_[i][1] == move_id)){

        return i;
      }
    }
    return -1;

  }


  Eigen::VectorXd sample_perturbations(const int index_from,
                                        const int num_samples,
                                        const int move_id,
                                        const Eigen::VectorXd residuals,
                                        const Eigen::MatrixXd &dH_dr,
                                        const Eigen::MatrixXd &dH_dt,
                                        const double delta){



    //number of probes
    int num_probes = sensors_.size();

    //number of total measurements
    int num_meas = residuals.size()/num_probes;

    //position of this move in memory
    //std::cout << "move id = " << move_id << std::endl;
    int move_pos = get_move_index(move_id);

    //std::cout << "compute_sparse_perturbation_matrix" << std::endl;
    //compute derivative matrix
    Eigen::SparseMatrix<double> dH = compute_sparse_perturbation_matrix(index_from,
                                                                        num_samples,
                                                                        dH_dr,
                                                                        dH_dt);

    //std::cout << "compute inverse measurement covariance" << std::endl;
    //helper identity matrix
    Eigen::SparseMatrix<double> I_n(5*num_samples,5*num_samples);
    I_n.setIdentity();

    Eigen::SparseMatrix<double> R_inv = get_R_inv(num_samples);

    //Cholesky factorization of R_inv
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt_of_R_inv(R_inv);

    //std::cout << "compute A" << std::endl;
    //with this we can compute the A matrix
    //we need to apply the sensor precision on the diagonal
    Eigen::MatrixXd A = dH.transpose()*R_inv*dH;


    //we now compute the prior for regularisation
    int pos_in_memory = get_position_in_move_prec_table(num_samples,move_id);

    Eigen::MatrixXd L = delta*prec_list_[pos_in_memory];//Pert_Covs_[move_pos].get_cov(num_samples);

    //std::cout << "compute rhs" << std::endl;
    //and now we construct the right hand side
    Eigen::VectorXd rhs(num_probes*num_samples);
    for(int i = 0; i < num_probes; ++i){
      rhs.segment(i*num_samples,num_samples) = residuals.segment(index_from+i*num_meas,num_samples);
    }


    //helper identity matrix
    Eigen::SparseMatrix<double> I_m(3*num_samples,3*num_samples);
    I_m.setIdentity();

    Eigen::VectorXd eps = generate_gaussian_samples_sparse(1, Eigen::VectorXd::Zero(3*num_samples), I_m, std::rand()).col(0);

    Eigen::VectorXd b = dH.transpose()*(R_inv*rhs + llt_of_R_inv.matrixL()*eps);

    // sample from the posterior by solving a randomized linear equation system
    Eigen::VectorXd d = (A + L).lu().solve(b);

    return d;

  }


  std::vector<double> get_noise_stdev(){
    return noise_stdev_;
  }



  void save_uq_matrices(const std::string output_dir){

    write_binary((output_dir + "H_mat.bin").c_str(), H_mat);
    write_binary((output_dir + "dHdx_mat.bin").c_str(), dHx_mat);
    write_binary((output_dir + "dHdy_mat.bin").c_str(), dHy_mat);
    write_binary((output_dir + "dHdz_mat.bin").c_str(), dHz_mat);
    write_binary((output_dir + "dHdt1_mat.bin").c_str(), dtH_1);
    write_binary((output_dir + "dHdt2_mat.bin").c_str(), dtH_2);

    /*
    // define the format you want, you only need one instance of this...
    const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

    std::ofstream out_file_H(output_dir + "H_mat.csv");
    out_file_H << H_mat.format(CSVFormat);
    out_file_H.close();

    std::ofstream out_file_dHx(output_dir + "dHdx_mat.csv");
    out_file_dHx << dHx_mat.format(CSVFormat);
    out_file_dHx.close();

    std::ofstream out_file_dHy(output_dir + "dHdy_mat.csv");
    out_file_dHy << dHy_mat.format(CSVFormat);
    out_file_dHy.close();

    std::ofstream out_file_dHz(output_dir + "dHdz_mat.csv");
    out_file_dHz << dHz_mat.format(CSVFormat);
    out_file_dHz.close();

    std::ofstream out_file_dtH_1(output_dir + "dHdt1_mat.csv");
    out_file_dtH_1 << dtH_1.format(CSVFormat);
    out_file_dtH_1.close();

    std::ofstream out_file_dtH_2(output_dir + "dHdt2_mat.csv");
    out_file_dtH_2 << dtH_2.format(CSVFormat);
    out_file_dtH_2.close();
    */
  }

  void load_uq_matrices(const std::string input_dir,const int arm_axis_id){

    read_binary((input_dir + "H_mat.bin").c_str(), H_mat);
    read_binary((input_dir + "dHdx_mat.bin").c_str(), dHx_mat);
    read_binary((input_dir + "dHdy_mat.bin").c_str(), dHy_mat);
    read_binary((input_dir + "dHdz_mat.bin").c_str(), dHz_mat);
    read_binary((input_dir + "dHdt1_mat.bin").c_str(), dtH_1);
    read_binary((input_dir + "dHdt2_mat.bin").c_str(), dtH_2);

    /*
    H_mat = read_csv(input_dir +  "H_mat.csv");
    dHx_mat = read_csv(input_dir +  "dHdx_mat.csv");
    dHy_mat = read_csv(input_dir +  "dHdy_mat.csv");
    dHz_mat = read_csv(input_dir +  "dHdz_mat.csv");
    dtH_1 = read_csv(input_dir +  "dHdt1_mat.csv");
    dtH_2 = read_csv(input_dir +  "dHdt2_mat.csv");
    */
    arm_axis_id_ = arm_axis_id;

  }

  Eigen::MatrixXd get_H(){

    return H_mat;
  }

  Eigen::MatrixXd get_dHx(){

    return dHx_mat;
  }

  Eigen::MatrixXd get_dHy(){

    return dHy_mat;
  }

  Eigen::MatrixXd get_dHz(){

    return dHz_mat;
  }

  Eigen::MatrixXd get_dHt1(){

    return dtH_1;
  }

  Eigen::MatrixXd get_dHt2(){

    return dtH_2;
  }


  Eigen::MatrixXcd compute_local_evaluation_matrix(const int L,const Eigen::Vector3d &eval_pos, const Eigen::Vector3d &cell_center){

    // number of coefficients
    int num_coeffs = (L+1)*(L+1);

    Eigen::MatrixXcd ret_mat(num_coeffs,sensors_.size());
    ret_mat.setZero();



    for (int i = 0; i < sensors_.size() ; ++i){

      ret_mat.block(0,i,num_coeffs,1) = sensors_[i].compute_local_evaluation_vector(L, eval_pos + positions_[i], cell_center );
    }

    return ret_mat;

  }

  std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> get_sensor_positions(){

    return positions_;
  }

  void compute_H_tilde(const Eigen::VectorXd &pert){

    int num_meas = pert.size()/5;
    int num_DoFs = H_mat.cols();
    int num_probes = sensors_.size();

    H_tilde_.resize(3*num_meas,num_DoFs);
    H_tilde_.setZero();

    //update
    /*
      */
    for(int m = 0; m < num_meas; ++m){
      for(int n = 0; n < num_DoFs-1 ; ++n){

        for(int p = 0; p < num_probes; ++p){

          H_tilde_(m+p*num_meas,n+1) = H_mat(m+p*num_meas,n+1)
                          + dHx_mat(m+p*num_meas,n+1)*pert(m)
                          + dHy_mat(m+p*num_meas,n+1)*pert(m+num_meas)
                          + dHz_mat(m+p*num_meas,n+1)*pert(m+2*num_meas);
                          + dtH_1(m+p*num_meas,n+1)*pert(m+3*num_meas)
                          + dtH_2(m+p*num_meas,n+1)*pert(m+4*num_meas);
        }

      }
    }


  }


  void compute_H_tilde(const Eigen::VectorXd &pert,
                        Eigen::MatrixXd &H_tilde){

    int num_meas = pert.size()/5;
    int num_DoFs = H_mat.cols();
    int num_probes = sensors_.size();




    //update
    /*
      */
    for(int m = 0; m < num_meas; ++m){
      for(int n = 0; n < num_DoFs-1 ; ++n){

        for(int p = 0; p < num_probes; ++p){

          H_tilde(m+p*num_meas,n+1) = H_mat(m+p*num_meas,n+1)
                          + dHx_mat(m+p*num_meas,n+1)*pert(m)
                          + dHy_mat(m+p*num_meas,n+1)*pert(m+num_meas)
                          + dHz_mat(m+p*num_meas,n+1)*pert(m+2*num_meas);
                          + dtH_1(m+p*num_meas,n+1)*pert(m+3*num_meas)
                          + dtH_2(m+p*num_meas,n+1)*pert(m+4*num_meas);
        }

      }
    }


  }

 private:
   //Hall Probes
   std::vector<DiscreteSensor<HallProbe<Pot,LinOp>,LinOp> , Eigen::aligned_allocator<DiscreteSensor<HallProbe<Pot,LinOp>,LinOp>> > sensors_;
   //positions
   std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d>> positions_;
   //ansatz space
   AnsatzSpace<LinOp> ansatz_space_;
   //aca matrices
   std::vector<ACAMatrixRed<HallProbe<Pot,LinOp>,LinOp>, Eigen::aligned_allocator<ACAMatrixRed<HallProbe<Pot,LinOp>,LinOp>>> aca_matrices_;
   //number of evaluations
   int num_eval_;
   //flag for gauged formulation
   bool enable_gauge_ = false;
   //stability vector
   Eigen::VectorXd stability_vec_;
   //H_[:,1]/s_1
   Eigen::VectorXd h_;

   //orientation of arm axis
   int arm_axis_id_ = -1;

   //Transfer matrices for displacements for different moves
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Fw_x;
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Fw_y;
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Fw_z;

   //Transfer matrices for slopes and different moves
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Ft_x;
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Ft_y;
   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> Ft_z;

   //Sparse matrix for derivatives of measurement operation
   Eigen::SparseMatrix<double> dH_;

   std::vector<PerturbationCovariance,Eigen::aligned_allocator<PerturbationCovariance>> Pert_Covs_;

   std::vector<Eigen::MatrixXd,Eigen::aligned_allocator<Eigen::MatrixXd>> prec_list_;

   std::vector<std::vector<int>> move_prec_table_;

   //move directions vector
   std::vector<int> move_directions;

   //to identity different move types
   std::vector<int> move_ids_;

   //noise std
   std::vector<double> noise_stdev_;


   //Matrices to H and the two dtH
   Eigen::MatrixXd H_mat,dtH_1,dtH_2;

   Eigen::MatrixXd H_tilde_;

   //Measurement derivative matrix
   Eigen::MatrixXd dHx_mat,dHy_mat,dHz_mat;




};

}  // namespace Bembel
#endif
