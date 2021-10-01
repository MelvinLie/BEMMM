#include <iostream>
#include <fstream>
#include <random>


#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>

#include <Bembel/AnsatzSpace>
#include <Bembel/Geometry>
#include <Bembel/H2Matrix>
#include <Bembel/Laplace>
#include <Bembel/LinearForm>
#include <Bembel/IO>
#include <Bembel/Devices>
#include <Bembel/Measurements>
#include <Bembel/ACA>


#include <chrono>


int main() {
  using namespace Bembel;
  using namespace Eigen;

  // some general variables


  // for time measurement
  std::chrono::time_point<std::chrono::steady_clock> t_start, t_end;
  std::chrono::duration<double> t_el;

  // for csv output
  const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");

  //mapper arm is orientated along the x axis
  int mapper_arm_axis_id = 0;


  std::cout << "*****************************************************" << std::endl;
  std::cout << "estimate RMS positioning error" << std::endl;
  std::cout << "*****************************************************" << std::endl;


  //****************************************************************************
  // Input filenames
  //****************************************************************************
  std::string geo_filename, meas_filename, sensor_props_filename, input_dir, output_dir, sol_filename;

  std::string input_parameter_file = "input_example_curved_domain.xml";
  IO::InputParameter params_in(input_parameter_file);

  params_in.print_parameters();

  geo_filename = params_in.mesh_filename_;
  meas_filename = params_in.meas_filename_;
  sensor_props_filename = params_in.sensor_props_filename_;
  output_dir = params_in.output_directory_;
  sol_filename = params_in.solution_filename_;
  input_dir = params_in.input_directory_;

  int p = params_in.p_[0];
  int h = params_in.h_[0];

  //****************************************************************************
  // Boundary mesh
  //****************************************************************************
  //Load geometry
  Geometry geometry(geo_filename);

  //****************************************************************************
  // Ansatz space
  //****************************************************************************
  AnsatzSpace<LaplaceHypersingularOperator> ansatz_space(geometry, h, p);

  int num_DoFs = ansatz_space.get_number_of_dofs();
  std::cout << " number of DoFs = " << num_DoFs << std::endl;

  //****************************************************************************
  // Measurement data
  //****************************************************************************
  //read csv file
  Eigen::MatrixXd in_data = read_csv(meas_filename,true);

  //Load measurement data
  MeasurementData meas;

  meas.set_data(in_data.block(0,0,in_data.rows(),3),in_data.block(0,3,in_data.rows(),3), in_data.block(0,6,in_data.rows(),1), in_data.block(0,7,in_data.rows(),1));


  //number of measurements
  int num_meas = meas.get_number_of_measurements();


  std::cout << " number of xyz measurements = " << num_meas << std::endl;
  std::cout << " total number of measurements = " << 3*num_meas << std::endl;

  //****************************************************************************
  // Make Hall probe array
  //****************************************************************************
  std::cout << " setup Hall probe array" << std::endl;
  HallProbeArray<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                    LaplaceHypersingularOperator> hall_array(ansatz_space);

  hall_array.set_gauged_formulation(true);
  hall_array.load_parameter_from_file(sensor_props_filename);

  hall_array.print_parameters();

  //intercept corrected measurement vector
  Eigen::VectorXd y = meas.flatten_measurements() - hall_array.get_offset_vector(num_meas);


  std::cout << " setup move properties" << std::endl;
  //setup move properties
  hall_array.append_move_property(0,input_dir + "D_move_x.csv");
  hall_array.append_move_property(1,input_dir + "D_move_y.csv");
  hall_array.append_move_property(2,input_dir + "D_move_z.csv");

  //number of sensors
  int num_sensors = hall_array.get_number_of_sensors();

  //****************************************************************************
  // Load least squares solution
  //****************************************************************************
  std::cout << " read least squares solution" << std::endl;
  Eigen::VectorXd v_ls = read_csv(sol_filename,true).col(0);

  std::cout << " v_ls = " << v_ls.rows() << " x " << v_ls.cols() << std::endl;
  //sparse perturbation matrix
  //Eigen::SparseMatrix<double> dH = hall_array.compute_sparse_perturbation_matrix(v_ls);

  //std::cout << " perturbation matrix = " << dH.rows() << " x " << dH.cols() << std::endl;


  //****************************************************************************
  // Load evaluation positions and assemble field evaluation
  //****************************************************************************
  std::cout << " read evaluation positions" << std::endl;
  Eigen::MatrixXd eval_pos = read_csv(input_dir + "curved_domain_evaluation_pos.csv",true);

  //make discrete potential
  DiscretePotential<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                    LaplaceHypersingularOperator> disc_pot(ansatz_space);

  int num_eval = eval_pos.rows();

  //****************************************************************************
  // Compute observation operator
  //****************************************************************************
  std::cout << " compute H" << std::endl;

  t_start =  std::chrono::steady_clock::now();

  hall_array.compute_H(meas.get_positions());

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;

  //****************************************************************************
  // Compute covariance matrix
  //****************************************************************************

  std::cout << " compute covariance matrix" << std::endl;

  t_start =  std::chrono::steady_clock::now();

  //Eigen::MatrixXd H_der = hall_array.compute_derivatives(meas.get_positions(),v_ls,mapper_arm_axis_id);
  //hall_array.prepare_UQ_matrices(meas.get_positions(),mapper_arm_axis_id);
  //hall_array.compute_H(meas.get_positions());

  Eigen::SparseMatrix<double,Eigen::RowMajor> R_d = hall_array.compute_R_decorrelated(meas,v_ls,mapper_arm_axis_id);

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;

  Eigen::SparseMatrix<double,Eigen::RowMajor> R_y = hall_array.get_R(num_meas);

  //****************************************************************************
  // Compute sparse LLT of covariance matrix
  //****************************************************************************
  std::cout << " sparse Cholesky facotrization" << std::endl;

  t_start =  std::chrono::steady_clock::now();

  Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> llt(R_d + R_y);//
  Eigen::SparseMatrix<double> L = llt.matrixL();//

  SparseLU<SparseMatrix<double>, COLAMDOrdering<int> >  L_solver;
  L_solver.analyzePattern(L);
  L_solver.factorize(L);


  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;

  //****************************************************************************
  // Compute Q
  //****************************************************************************
  std::cout << " compute Q" << std::endl;

  Eigen::MatrixXd Q;

  t_start =  std::chrono::steady_clock::now();

  hall_array.compute_Q(llt,Q);

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;

  //****************************************************************************
  // Compute A
  //****************************************************************************
  std::cout << " compute A" << std::endl;

  Eigen::MatrixXd A;

  t_start =  std::chrono::steady_clock::now();

  hall_array.compute_A(Q,A);

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;


  //****************************************************************************
  // Compute u
  //****************************************************************************
  std::cout << " compute u" << std::endl;

  Eigen::VectorXd u = llt.solve(y);



  //****************************************************************************
  // Compute eps
  //****************************************************************************
  std::cout << " compute eps" << std::endl;

  //helper identity matrix
  Eigen::SparseMatrix<double> I(3*num_meas,3*num_meas);
  I.setIdentity();

  int num_samples = 10;

  Eigen::MatrixXd eps = generate_gaussian_samples_sparse_covariance(num_samples, Eigen::VectorXd::Zero(3*num_meas),I);

  //****************************************************************************
  // Compute v
  //****************************************************************************
  std::cout << " compute v" << std::endl;

  Eigen::MatrixXd V(num_DoFs,num_samples);
  Eigen::VectorXd L_eps;
  Eigen::VectorXd v(num_DoFs-1);
  Eigen::VectorXd rhs;
  Eigen::VectorXd v_mean(num_DoFs);

  Eigen::MatrixXd Bx_mat(num_eval,num_samples);
  Eigen::MatrixXd By_mat(num_eval,num_samples);
  Eigen::MatrixXd Bz_mat(num_eval,num_samples);
  v_mean.setZero();

  t_start =  std::chrono::steady_clock::now();

  for(int i = 0; i < num_samples; ++i){

    //normal transformation
    L_eps = L_solver.solve(eps.col(i));
    //right hand side
    rhs = hall_array.compute_rhs(u,L_eps);
    //solve
    v = A.lu().solve(rhs);
    V.col(i) = hall_array.recover_full_solution(v);

    v_mean += V.col(i);

    //evaluate
    //set cauchy data
    disc_pot.set_cauchy_data(V.col(i));

    Eigen::MatrixXd B = disc_pot.evaluate_der(eval_pos).eval();
    Bx_mat.col(i) = B.col(0);
    By_mat.col(i) = B.col(1);
    Bz_mat.col(i) = B.col(2);

    //std::cout << "******************************" << std::endl;
    //std::cout << "sample " << i << std::endl;
    //std::cout << "******************************" << std::endl;
    //std::cout << V.col(i) << std::endl;
    //std::cout << "******************************" << std::endl;

  }

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << " elapsed time = " << t_el.count() << " sec" << std::endl;

  v_mean /= num_samples;

  //compute residuals
  Eigen::VectorXd residuals = hall_array.evaluate_H(v_mean) - y;

  double rms_error = std::sqrt((residuals.segment(0,num_meas).array().pow(2)
                    + residuals.segment(num_meas,num_meas).array().pow(2)
                    + residuals.segment(2*num_meas,num_meas).array().pow(2)).sum()/3./num_meas);

  std::cout << " RMS error " << rms_error << " V" << std::endl;



  std::ofstream out_file_V(output_dir + "samples.csv");
  out_file_V << V.format(CSVFormat);
  out_file_V.close();

  std::ofstream out_file_v_mean(output_dir + "v_mean.csv");
  out_file_v_mean << v_mean.format(CSVFormat);
  out_file_v_mean.close();

  std::ofstream out_file_Bx(output_dir + "Bx_samples.csv");
  out_file_Bx << Bx_mat.format(CSVFormat);
  out_file_Bx.close();

  std::ofstream out_file_By(output_dir + "By_samples.csv");
  out_file_By << By_mat.format(CSVFormat);
  out_file_By.close();

  std::ofstream out_file_Bz(output_dir + "Bz_samples.csv");
  out_file_Bz << Bz_mat.format(CSVFormat);
  out_file_Bz.close();

  return 0;
}
