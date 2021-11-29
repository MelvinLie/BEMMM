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

  std::cout << "*****************************************************" << std::endl;
  std::cout << "Run Gibbs sampler for UQ" << std::endl;
  std::cout << "*****************************************************" << std::endl;

  //for time measurement
  std::chrono::time_point<std::chrono::steady_clock> t_start, t_end;
  std::chrono::duration<double> t_el;

  // define the format you want, you only need one instance of this...
  const static Eigen::IOFormat CSVFormat(16, Eigen::DontAlignCols, ", ", "\n");


  //****************************************************************************
  // Input parameters
  //****************************************************************************
  std::string geo_filename, meas_filename, sensor_props_filename, input_dir, output_dir, sol_filename;

  std::string input_parameter_file = "input_example_gibbs_sampling.xml";
  IO::InputParameter params_in(input_parameter_file);

  params_in.print_parameters();

  geo_filename = params_in.mesh_filename_;
  meas_filename = params_in.meas_filename_;
  sensor_props_filename = params_in.sensor_props_filename_;
  output_dir = params_in.output_directory_;
  sol_filename = params_in.solution_filename_;
  input_dir = params_in.input_directory_;

  //polynomial degree
  int p = params_in.p_[0];
  //refinement level
  int h = params_in.h_[0];

  //number of samples to generate
  int num_samples = 100;

  //arm axis id. the arm is aligned along the x axis
  int arm_axis_id = 0;

  //regularisation parameter for perturbation estimation
  double delta = 1.;

  //****************************************************************************
  // Domain of Interest
  //****************************************************************************
  //Load geometry
  Geometry geometry(geo_filename);

  //****************************************************************************
  // Ansatz space
  //****************************************************************************
  AnsatzSpace<LaplaceHypersingularOperator> ansatz_space(geometry, h, p);

  int num_DoFs = ansatz_space.get_number_of_dofs();
  std::cout << "Number of DoFs = " << num_DoFs << std::endl;

  //****************************************************************************
  // Measurement data
  //****************************************************************************
  //the measurement file is formatted as follows:
  // x , y , z , meas number , meas is , Ux , Uy , Uz
  Eigen::MatrixXd meas_data = read_csv(meas_filename,true);

  //Load measurement data
  MeasurementData meas;

  meas.set_data(meas_data.block(0,0,meas_data.rows(),3),meas_data.block(0,3,meas_data.rows(),3), meas_data.block(0,6,meas_data.rows(),1), meas_data.block(0,7,meas_data.rows(),1) );


  //number of measurements
  int num_meas = meas.get_number_of_measurements();

  std::cout << "Number of xyz measurements = " << num_meas << std::endl;
  std::cout << "Total number of measurements = " << 3*num_meas << std::endl;

  //****************************************************************************
  // Make Hall probe array
  //****************************************************************************
  HallProbeArray<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                    LaplaceHypersingularOperator> hall_array(ansatz_space);

  hall_array.set_gauged_formulation(true);
  hall_array.load_parameter_from_file(sensor_props_filename);
  hall_array.set_arm_axis_orientation(arm_axis_id);

  hall_array.print_parameters();

  //intercept corrected measurement vector
  Eigen::VectorXd y = meas.flatten_measurements() - hall_array.get_offset_vector(num_meas);

  //hall array error matrices
  Eigen::SparseMatrix<double> R_inv = hall_array.get_R_inv(num_meas);
  Eigen::SparseMatrix<double> R = hall_array.get_R(num_meas);

  //Cholesky factorization of R_inv
  Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > llt_of_R_inv(R_inv);
  Eigen::SparseMatrix<double> L = llt_of_R_inv.matrixL();

  //setup move properties
  hall_array.append_move_property(0,input_dir + "cov_x_v_30_dist_2_mm.csv");
  hall_array.append_move_property(1,input_dir + "cov_y_v_30_dist_2_mm.csv");
  hall_array.append_move_property(2,input_dir + "cov_z_v_30_dist_2_mm.csv");
  //the third move was verified by Leica and taken with a higher sampling frequency
  hall_array.append_move_property(3,input_dir + "cov_x_v_30_dist_0.3_mm.csv");

  //number of sensors
  int num_sensors = hall_array.get_number_of_sensors();

  //****************************************************************************
  // Prepare UQ matrices
  //****************************************************************************
  bool load_uq_matrices = true;
  bool save_uq_matrices = false;

  if(load_uq_matrices){

    hall_array.load_uq_matrices(output_dir,arm_axis_id);
  }
  else{
    std::cout << "Prepare UQ-matrices" << std::endl;

    t_start =  std::chrono::steady_clock::now();

    hall_array.prepare_UQ_matrices(meas.get_positions(),arm_axis_id);

    t_end =  std::chrono::steady_clock::now();
    t_el = t_end - t_start;

    std::cout << "\telapsed time  = " << t_el.count()/60. << " min"<< std::endl;
  }
  if(save_uq_matrices){

    hall_array.save_uq_matrices(output_dir);
  }

  std::cout << "get move sepatators" << std::endl;
  //move separator
  std::vector<int> move_sep = meas.get_move_separator();
  //move identifier
  std::vector<int> move_ids = meas.get_move_ids();

  std::cout << "setup the precision matrices" << std::endl;
  hall_array.setup_precision_matrices(move_sep, move_ids);

  //****************************************************************************
  // Make space for result data
  //****************************************************************************
  std::cout << "Make space for result data" << std::endl;


  Eigen::MatrixXd v_container(num_DoFs,num_samples); v_container.setZero();

  //if perturbations are of interest

  Eigen::MatrixXd wx_container(num_meas,num_samples); wx_container.setZero();
  Eigen::MatrixXd wy_container(num_meas,num_samples); wy_container.setZero();
  Eigen::MatrixXd wz_container(num_meas,num_samples); wz_container.setZero();
  Eigen::MatrixXd t1_container(num_meas,num_samples); t1_container.setZero();
  Eigen::MatrixXd t2_container(num_meas,num_samples); t2_container.setZero();
  /*
  */


  //current field solution
  Eigen::VectorXd this_v(num_DoFs); this_v.setZero();
  //current perturbations
  Eigen::VectorXd this_p(5*num_meas); this_p.setZero();


  //****************************************************************************
  // Gibbs sampling
  //****************************************************************************
  //helper identity matrix
  Eigen::SparseMatrix<double> I(3*num_meas,3*num_meas);
  I.setIdentity();

  //helper identity matrix
  Eigen::SparseMatrix<double> I_N(num_DoFs-1,num_DoFs-1);
  I_N.setIdentity();


  //container for derivatives
  Eigen::MatrixXd dH_dr,dH_dt;

  //container for residuals
  Eigen::VectorXd r;

  Eigen::VectorXd v_mean(num_DoFs);
  v_mean.setZero();

  //****************************************************************************
  // position of validation move is saved for comparison to Leica
  // we make a container to store the results
  //****************************************************************************
  int validation_move = 60;

  //containter for validation data
  Eigen::MatrixXd val_meas_data =  meas.get_move(validation_move);

  int num_meas_val = move_sep[validation_move+1] - move_sep[validation_move];

  // we save out the vertical displacement, for validation versus Leica
  Eigen::MatrixXd val_data(num_meas_val,num_samples + 2);


  val_data.block(0,0,num_meas_val,2) = read_csv(input_dir + "optical_measurement.csv",true);


  std::cout << "Start Gibbs sampling" << std::endl;

  //plot solution
  VTKSurfaceExport writer(geometry, h);

  FunctionEvaluator<LaplaceHypersingularOperator> evaluator_v(ansatz_space);

  std::function<double(int, const Eigen::Vector2d &)> density_v =
      [&](int patch_number,
          const Eigen::Vector2d &reference_domain_point) {
        return evaluator_v.evaluateOnPatch(patch_number,
                                           reference_domain_point)(0);
                                          };


  //make discrete potential for debugging
  DiscretePotential<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                  LaplaceHypersingularOperator> disc_pot(ansatz_space);



  Eigen::MatrixXd Bx_container(num_meas_val,num_samples); Bx_container.setZero();
  Eigen::MatrixXd By_container(num_meas_val,num_samples); By_container.setZero();
  Eigen::MatrixXd Bz_container(num_meas_val,num_samples); Bz_container.setZero();

  Eigen::MatrixXd res_x_container(num_meas,num_samples); res_x_container.setZero();
  Eigen::MatrixXd res_y_container(num_meas,num_samples); res_y_container.setZero();
  Eigen::MatrixXd res_z_container(num_meas,num_samples); res_z_container.setZero();

  // construct prior
  //Eigen::VectorXd v_0 = read_csv(input_dir + "prior.csv",true).col(0).segment(1,num_DoFs-1);
  // prior covariance matrix
  //Eigen::SparseMatrix<double> Q(num_DoFs-1,num_DoFs-1);
  //Q.setIdentity();

  //double delta = 0.;

  t_start =  std::chrono::steady_clock::now();



  for(int i = 0; i < num_samples; ++i){

    std::cout << " generating sample v_" << i << std::endl;
    //--------------------------------------------------------------------------
    //Sample v
    //--------------------------------------------------------------------------

    hall_array.compute_H_tilde(this_p);

    //compute A
    Eigen::MatrixXd A = hall_array.compute_A_perturbed(R_inv);


    //compute rhs
    Eigen::VectorXd eps = generate_gaussian_samples_sparse(1, 0.*y, I, std::rand()).col(0);

    std::cout << " compute_rhs_perturbed"  << std::endl;
    Eigen::VectorXd rhs = hall_array.compute_rhs_perturbed(R_inv,y, L*eps);


    //rhs += delta*Q*v_0;

    //eps = generate_gaussian_samples_sparse(1, 0.*v_0, I_N, std::rand()).col(0);
    //rhs += std::sqrt(delta)*Q*eps;

    //draw sample
    Eigen::VectorXd v_s = (A).lu().solve(rhs);

    //go back to full vector
    this_v = hall_array.recover_full_solution(v_s);

    v_mean += this_v;

    std::cout << " residual"  << std::endl;
    Eigen::VectorXd r = meas.flatten_measurements() - hall_array.evaluate_H_tilde(this_v) - hall_array.get_offset_vector(num_meas);


    evaluator_v.set_function(this_v);
    writer.addDataSet("v_" + std::to_string(i) , density_v);

    disc_pot.set_cauchy_data(this_v);

    Eigen::MatrixXd B = disc_pot.evaluate_der(val_meas_data.block(0,0,num_meas_val,3)).eval();


    Bx_container.col(i) = B.block(0,0,num_meas_val,1);
    By_container.col(i) = B.block(0,1,num_meas_val,1);
    Bz_container.col(i) = B.block(0,2,num_meas_val,1);


    res_x_container.col(i) = r.segment(0,num_meas);
    res_y_container.col(i) = r.segment(num_meas,num_meas);
    res_z_container.col(i) = r.segment(2*num_meas,num_meas);


    //--------------------------------------------------------------------------
    //Sample perturbations
    //--------------------------------------------------------------------------

    Eigen::VectorXd this_y = hall_array.evaluate_H_tilde(this_v) + hall_array.get_offset_vector(num_meas);


    //compute derivatives
    dH_dr = hall_array.evaluate_grad_H(this_v);
    dH_dt = hall_array.evaluate_dH_dt(this_v);

    std::cout << " generating sample d_" << i << std::endl;

    //iterate over moves
    //#pragma omp parallel for
    for(int m = 0; m < move_sep.size()-1 ; ++m){

      //number of positions in this move
      int num_pos = move_sep[m+1] - move_sep[m];


      Eigen::VectorXd my_pert = hall_array.sample_perturbations(move_sep[m],num_pos,move_ids[m],r,dH_dr,dH_dt,delta);

      //critical
      //#pragma omp critical
      {
        this_p.segment(move_sep[m],num_pos)             = my_pert.segment(0,num_pos);
        this_p.segment(move_sep[m]+num_meas,num_pos)    = my_pert.segment(num_pos,num_pos);
        this_p.segment(move_sep[m]+2*num_meas,num_pos)  = my_pert.segment(2*num_pos,num_pos);
        this_p.segment(move_sep[m]+3*num_meas,num_pos)  = my_pert.segment(3*num_pos,num_pos);
        this_p.segment(move_sep[m]+4*num_meas,num_pos)  = my_pert.segment(4*num_pos,num_pos);
      }

      if(m == validation_move){
        val_data.block(0,i + 2,num_meas_val,1) = my_pert.segment(num_pos,num_pos);
      }

      /*
      else{
        //test, we take only the valdation move
        this_p.segment(move_sep[m],5*num_pos)             *= 0.;
        this_p.segment(move_sep[m]+num_meas,5*num_pos)    *= 0.;
        this_p.segment(move_sep[m]+2*num_meas,5*num_pos)  *= 0.;
        this_p.segment(move_sep[m]+3*num_meas,5*num_pos)  *= 0.;
        this_p.segment(move_sep[m]+4*num_meas,5*num_pos)  *= 0.;
      }
      */

    }


    //store samples if needed

    v_container.col(i) = this_v;
    wx_container.col(i) = this_p.segment(0,num_meas);
    wy_container.col(i) = this_p.segment(num_meas,num_meas);
    wz_container.col(i) = this_p.segment(2*num_meas,num_meas);
    t1_container.col(i) = this_p.segment(3*num_meas,num_meas);
    t2_container.col(i) = this_p.segment(4*num_meas,num_meas);
    /*
    */
  }

  t_end =  std::chrono::steady_clock::now();
  t_el = t_end - t_start;

  std::cout << "elapsed time for " << num_samples << " samples: " << t_el.count() << " sec." <<  std::endl;




  std::cout << "save out everything" << std::endl;


  //save out everything

  std::ofstream out_file_Wx(output_dir + "Wx.csv");
  out_file_Wx << wx_container.format(CSVFormat);
  out_file_Wx.close();

  std::ofstream out_file_Wy(output_dir + "Wy.csv");
  out_file_Wy << wy_container.format(CSVFormat);
  out_file_Wy.close();

  std::ofstream out_file_Wz(output_dir + "Wz.csv");
  out_file_Wz << wz_container.format(CSVFormat);
  out_file_Wz.close();

  std::ofstream out_file_Ty(output_dir + "Ty.csv");
  out_file_Ty << t1_container.format(CSVFormat);
  out_file_Ty.close();

  std::ofstream out_file_Tz(output_dir + "Tz.csv");
  out_file_Tz << t2_container.format(CSVFormat);
  out_file_Tz.close();
  /*
  */

  std::ofstream out_file_V(output_dir + "V.csv");
  out_file_V << v_container.format(CSVFormat);
  out_file_V.close();

  std::ofstream out_file_val(output_dir + "val.csv");
  out_file_val << val_data.format(CSVFormat);
  out_file_val.close();


  std::ofstream out_file_Bx(output_dir + "Bx.csv");
  out_file_Bx << Bx_container.format(CSVFormat);
  out_file_Bx.close();

  std::ofstream out_file_By(output_dir + "By.csv");
  out_file_By << By_container.format(CSVFormat);
  out_file_By.close();

  std::ofstream out_file_Bz(output_dir + "Bz.csv");
  out_file_Bz << Bz_container.format(CSVFormat);
  out_file_Bz.close();

  std::ofstream out_file_res_x(output_dir + "res_x.csv");
  out_file_res_x << res_x_container.format(CSVFormat);
  out_file_res_x.close();

  std::ofstream out_file_res_y(output_dir + "res_y.csv");
  out_file_res_y << res_y_container.format(CSVFormat);
  out_file_res_y.close();

  std::ofstream out_file_res_z(output_dir + "res_z.csv");
  out_file_res_z << res_z_container.format(CSVFormat);
  out_file_res_z.close();

  evaluator_v.set_function(v_mean/num_samples);
  writer.addDataSet("v_mean", density_v);

  writer.writeToFile(output_dir +  "v.vtp");



  return 0;
}
