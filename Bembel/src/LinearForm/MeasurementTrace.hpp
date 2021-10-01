// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_LINEARFORM_MEASUREMENTTRACE_H_
#define BEMBEL_LINEARFORM_MEASUREMENTTRACE_H_

namespace Bembel {

template <typename Scalar>
class MeasurementTrace;

template <typename ScalarT>
struct LinearFormTraits<MeasurementTrace<ScalarT>> {
  typedef ScalarT Scalar;
};

/**
 *  \ingroup LinearForm
 *  \brief This class provides an implementation of the Neumann trace operator
 * and a corresponding method to evaluate the linear form corresponding to the
 * right hand side of the system via quadrature.
 */
template <typename Scalar>
class MeasurementTrace : public LinearFormBase<MeasurementTrace<Scalar>, Scalar> {
 public:
  MeasurementTrace() {}

  double read_data(const std::string &filename, int N_in){

    N = N_in;
    z.clear();
    An.clear();
    Bn.clear();

    An.resize(N);
    Bn.resize(N);

    std::ifstream file_stream;
    std::string current_line;
    std::string current_element;

    file_stream.open(filename);
    if (!file_stream) {
      std::cerr << "File " << filename << " doesn't exist!" << std::endl;
      exit(1);
    }

    for( std::string current_line; getline( file_stream, current_line ); )
    {
      std::stringstream current_data(current_line);
      std::vector<double> tmp_vector;

      while(std::getline(current_data,current_element,',')){
        tmp_vector.push_back(atof(current_element.c_str()));

      }
      z.push_back(tmp_vector[1]);
      for(int i = 0; i < N ; ++i){
        An[i].push_back(tmp_vector[2+2*i]);
        Bn[i].push_back(tmp_vector[3+2*i]);
      }

    }
    file_stream.close();

  }


  template <class T>
  void evaluateIntegrand_impl(
      const T &super_space, const SurfacePoint &p,
      Eigen::Matrix<Scalar, Eigen::Dynamic, 1> *intval) const {
    auto polynomial_degree = super_space.get_polynomial_degree();
    auto polynomial_degree_plus_one_squared =
        (polynomial_degree + 1) * (polynomial_degree + 1);


    // get evaluation points on unit square
    auto s = p.segment<2>(0);

    auto h = 1. / (1 << super_space.get_refinement_level());

    //std::cout << "h = "<<  h << std::endl;
    // get quadrature weights
    auto ws = p(2);//*h;//*h

    // get points on geometry and tangential derivatives
    auto x_f = p.segment<3>(3);
    auto x_f_dx = p.segment<3>(6);
    auto x_f_dy = p.segment<3>(9);

    // compute normal vector direction
    auto n = x_f_dx.cross(x_f_dy);

    // compute surface measures from tangential derivatives
    auto x_kappa = n.norm();

    //scale normal vector WE WOULD NOT NEED THIS IF WE WOULD REMOVE x_kappa BELOW
    n /= x_kappa;

    //evaluare gradient function

    double func_eval;
    if(std::abs(n(2)) > 0.97){
      func_eval = 0;
    }
    else{
      func_eval = evaluate(x_f);
    }


    // integrand without basis functions
    auto integrand = func_eval * x_kappa * ws; //x_kappa



    // multiply basis functions with integrand
    super_space.addScaledBasis(intval, integrand, s);


    return;
  };


 private:

   std::vector<double> z;
   std::vector<std::vector<double>> An,Bn;
   int N;


   double evaluate(const Eigen::Vector3d &in) const {

     double phi = std::atan2(in(1),in(0));
     double B_ret = 0.;
     double An_int,Bn_int;

     //bool interval_hit = false;

       for (int i = 0; i < z.size()-1; ++i){


         if ((z[i] < in(2) ) && (z[i+1] > in(2) )){

           //interval_hit = true;
           for(int n = 1; n < N+1 ; ++n){
             //linear interpolation
             An_int = (An[n-1][i+1]-An[n-1][i])*(in(2)-z[i])/(z[i+1]-z[i]) + An[n-1][i];
             Bn_int = (Bn[n-1][i+1]-Bn[n-1][i])*(in(2)-z[i])/(z[i+1]-z[i]) + Bn[n-1][i];

             B_ret += An_int*std::cos(n*phi) - Bn_int*std::sin(n*phi);
           }
         }
       }

       //if(interval_hit == false) std::cout << "not hit" << std::endl;
       return B_ret;
   };


};
}  // namespace Bembel

#endif
