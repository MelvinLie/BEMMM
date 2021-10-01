// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_SPHERICALHARM_H_
#define BEMBEL_SPHERICALHARM_H_

#include <iostream>
#include <fstream>

#include <sys/sysinfo.h>
#include <unistd.h>

//Maybe for the future to speed up the computations
//https://www.boost.org/doc/libs/1_67_0/libs/math/doc/html/math_toolkit/sf_poly/legendre.html
#include <boost/math/special_functions/legendre.hpp>

namespace Bembel {

  double factorial(int I){

    double ret_val = 1.;

    for(int i = 2 ; i < I + 1 ; ++i){
      ret_val *= i;
    }
    return ret_val;

  }

  /*
  double associated_legendre(double x, int l,int m){

    //https://www.osti.gov/servlets/purl/4370075
    double p_lm = 0.;
    int m_abs = m; if(m_abs < 0) m_abs *= -1;
    double x_abs = x; if(x_abs < 0) x_abs *= -1;
    double arg = (1.-x_abs)/2.;


    for(int n = 0; n < l+1; ++n){

      //std::cout << std::endl;
      //std::cout << "n = " << n << std::endl;
      //std::cout << "fac_lpn = " << factorial(l+n) << std::endl;
      //std::cout << "fac_lmn = " << factorial(l-n) << std::endl;
      //std::cout << "fac_mpn = " << factorial(m_abs+n) << std::endl;
      //std::cout << "fac_n = " << factorial(n) << std::endl;
      //std::cout << "arg = " << std::pow(arg,n) << std::endl;
      //std::cout << "pow_m1 = " << std::pow(-1.,n) << std::endl;

      //think about recursively computing the factorials. Would help recuding fops
      p_lm += factorial(l+n)/factorial(l-n)*std::pow(-1.,n)/factorial(m_abs+n)/factorial(n)*std::pow(arg,n);

    }
    p_lm *= std::pow((1.-x_abs)/(1.+x_abs),0.5*m_abs);

    if(x < 0){
      p_lm *= std::pow(-1.,l+m_abs);
    }

    if(m > 0){
      //std::cout << "l = " << l <<std::endl;
      //std::cout << "m = " << m <<std::endl;
      //std::cout << "(l+m)! = " << factorial(l+m) <<std::endl;
      //std::cout << "(l-m)! = " << factorial(l-m) <<std::endl;

      p_lm *= std::pow(-1.,m_abs)*factorial(l+m)/factorial(l-m);

    }


    return p_lm;

    }
    */

    double associated_legendre(const int l,const int m,const double x){

      //https://www.osti.gov/servlets/purl/4370075
      double p_lm = 0.;
      int m_abs = m; if(m_abs < 0) m_abs *= -1;
      double x_abs = x; if(x_abs < 0) x_abs *= -1;

      double arg = 1.;//(1.-x_abs)/2.;
      double arg_inc = (1.-x_abs)/2.;

      double fac_lpn = factorial(l);
      double fac_lmn = fac_lpn;
      double fac_mpn = factorial(m_abs);
      double fac_n = 1.;
      double pow_m1 = 1.;

      p_lm = fac_lpn/fac_lmn*pow_m1/fac_mpn/fac_n*arg;

      for(int n = 1; n < l+1; ++n){

        fac_lpn *= (l+n);
        fac_lmn /= (l-n+1);
        fac_mpn *= (m_abs+n);
        fac_n *= n;
        arg *= arg_inc;
        pow_m1 *= -1;

        p_lm += fac_lpn/fac_lmn*pow_m1/fac_mpn/fac_n*arg;


      }
      p_lm *= std::pow((1.-x_abs)/(1.+x_abs),0.5*m_abs);

      if(x < 0){
        p_lm *= std::pow(-1.,l+m_abs);
      }

      if(m > 0){


        p_lm *= std::pow(-1.,m_abs)*factorial(l+m)/factorial(l-m);

      }


      return p_lm;

    }
    //Here we compute the derivative of the associated legendre polynomial linked
    //with the polar angle theta by: P_l^m \circ \cos(\theta)
    double associated_legendre_derivative(const int l,const int m,const double theta){
    //https://math.stackexchange.com/questions/3369949/derivative-of-normalized-associated-legendre-function-at-the-limits-of-x-1
    //https://math.stackexchange.com/questions/391672/derivative-of-associated-legendre-polynomials-at-x-pm-1

     double ret_val = 0.;
     double x = std::cos(theta);
     double m_abs = fabs(m);


     //Deal with singularity at \theta = 0
     if(fabs(theta) < 1e-8){
       if(m_abs == 1){
         ret_val = -0.5*l*(l+1);

         if(m < 0){


           ret_val *= std::pow(-1.,m_abs)*factorial(l-m_abs)/factorial(l+m_abs);

         }
       }
     }
     //Deal with singularity at \theta = \pm \pi
     else if(fabs(theta)  > BEMBEL_PI - 1e-8){


       if(m_abs == 1){
         ret_val = 0.5*std::pow(-1,l+1)*l*(l+1);


         if(m < 0){


           ret_val *= std::pow(-1.,m_abs)*factorial(l-m_abs)/factorial(l+m_abs);

        }
       }
     }
     else{

       ret_val = ((l+1-m)*associated_legendre(l+1,m,x) - (l+1)*x*associated_legendre(l,m,x))/(x*x - 1);
       ret_val *= -1.*std::sin(theta);
     }


     return ret_val;

    }

    std::complex<double> Ylm(const int l, const int m,const double theta,const double phi){

      double P_lm = associated_legendre(l , m, std::cos(theta));
      double N_lm = std::sqrt((2.*l+1.)/4./BEMBEL_PI*factorial(l-m)/factorial(l+m));

      //std::cout << "P = " << P_lm << std::endl;
      //std::cout << "N = "<< N_lm << std::endl;
      //std::cout << "factorial(l-m) = "<< factorial(l-m) << std::endl;
      //std::cout << "factorial(l+m) = "<< factorial(l+m) << std::endl;
      //std::cout << "m = "<< m << std::endl;
      //std::cout << "l = "<< l << std::endl;

      std::complex<double> Y_lm;


      Y_lm.real(P_lm*N_lm*std::cos(m*phi));
      Y_lm.imag(P_lm*N_lm*std::sin(m*phi));

      //std::cout << Y_lm << std::endl;

      return Y_lm;

    }


    std::complex<double> dYlm_dt(const int l, const int m,const double theta,const double phi){

      double dP_lm_dt = associated_legendre_derivative(l , m, theta);
      double N_lm = std::sqrt((2.*l+1.)/4./BEMBEL_PI*factorial(l-m)/factorial(l+m));

      //std::cout << "P = " << P_lm << std::endl;
      //std::cout << "N = "<< N_lm << std::endl;
      //std::cout << "factorial(l-m) = "<< factorial(l-m) << std::endl;
      //std::cout << "factorial(l+m) = "<< factorial(l+m) << std::endl;
      //std::cout << "m = "<< m << std::endl;
      //std::cout << "l = "<< l << std::endl;

      std::complex<double> dY_lm_dt;


      dY_lm_dt.real(dP_lm_dt*N_lm*std::cos(m*phi));
      dY_lm_dt.imag(dP_lm_dt*N_lm*std::sin(m*phi));

      //std::cout << Y_lm << std::endl;

      return dY_lm_dt;

    }

    std::complex<double> dYlm_dp(const int l, const int m,const double theta,const double phi){

      double P_lm = associated_legendre(l , m, std::cos(theta));
      double N_lm = std::sqrt((2.*l+1.)/4./BEMBEL_PI*factorial(l-m)/factorial(l+m));

      //std::cout << "P = " << P_lm << std::endl;
      //std::cout << "N = "<< N_lm << std::endl;
      //std::cout << "factorial(l-m) = "<< factorial(l-m) << std::endl;
      //std::cout << "factorial(l+m) = "<< factorial(l+m) << std::endl;
      //std::cout << "m = "<< m << std::endl;
      //std::cout << "l = "<< l << std::endl;

      std::complex<double> dY_lm_dp;


      dY_lm_dp.real(-m*P_lm*N_lm*std::sin(m*phi));
      dY_lm_dp.imag(m*P_lm*N_lm*std::cos(m*phi));

      //std::cout << Y_lm << std::endl;

      return dY_lm_dp;

    }

    Eigen::MatrixXcd evaluate_harmonic_expansion(const int L, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      Eigen::MatrixXcd ret_mat;
      ret_mat.resize(num_coeffs,2);

      //running index
      int k = 0;



      //slow implementation
      for(int l = 0; l <= L; ++l){

        for(int m = -l; m <= l ; ++m){


          ret_mat(k,0) = Ylm(l,m,theta,phi);
          ret_mat(k,1) = dYlm_dt(l,m,theta,phi);

          ++k;

        }
      }

      return ret_mat;
    }

    //Problematic at this point is only the component R_phi for theta -> 0 and pi
    // From Hosbitals rule we get:
    // lim theta -> 0 R_phi = \partial_\theta Y_lm / \cos(\theta) \|_\theta = 0
    // lim theta -> pi R_phi = \partial_\theta Y_lm / \cos(\theta) \|_\theta = pi
    Eigen::MatrixXcd evaluate_derivatives_of_solid_harmonics(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      Eigen::MatrixXcd sph_harms =  evaluate_harmonic_expansion(L,theta,phi);

      Eigen::VectorXcd Rr(num_coeffs);
      Eigen::VectorXcd Rt(num_coeffs);
      Eigen::VectorXcd Rp(num_coeffs);
      Rr.setZero();
      Rt.setZero();
      Rp.setZero();

      //denumerator for S_phi. Again we need to take care of singulatities.
      // From Hospitals rule, we can compute the limit for theta -> 0
      // lim = d_Y_lm_dt / cos(theta) |_theta = 0
      double den_Rp = std::sin(theta);
      int index_Rp = 0;

      if(fabs(theta) < 1e-8){
        den_Rp = 1.;
        index_Rp = 1;

      }
      else if(fabs(theta)  > BEMBEL_PI - 1e-8){
        den_Rp = -1.;
        index_Rp = 1;
      }

      int k = 1;

      for(int l = 1; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          Rr(k) = l*std::pow(r,l-1)*sph_harms(k,0);
          Rt(k) = std::pow(r,l-1)*sph_harms(k,1);
          Rp(k) = 1.*m*I*std::pow(r,l-1)/den_Rp*sph_harms(k,index_Rp);

          ++k;

        }
      }


      Eigen::MatrixXcd ret_mat(num_coeffs,3);

      ret_mat.col(0) = std::sin(theta)*std::cos(phi)*Rr
                            + std::cos(theta)*std::cos(phi)*Rt
                            - std::sin(phi)*Rp;

      ret_mat.col(1) = std::sin(theta)*std::sin(phi)*Rr
                            + std::cos(theta)*std::sin(phi)*Rt
                            + std::cos(phi)*Rp;

      ret_mat.col(2) = std::cos(theta)*Rr - std::sin(theta)*Rt;


      return ret_mat;


    }




    Eigen::VectorXd recursive_legendre_evaluation(const int L, const int m, const double x){

      //number of coefficients along m
      int num_coeffs = L + 1 - m;

      //return values
      Eigen::VectorXd ret_vec(num_coeffs);

      ret_vec(0) = boost::math::legendre_p(m, m, x);

      if(num_coeffs > 1){
        ret_vec(1) = boost::math::legendre_p(m + 1, m, x);

        for (int i = 2; i < num_coeffs; ++i){
          ret_vec(i) = boost::math::legendre_next(m + i - 1, m,x, ret_vec(i-1), ret_vec(i-2));
        }
      }


      return ret_vec;

    }


    Eigen::VectorXd recurse_all_legendre_polynomials(const int L, const double x){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //return values
      Eigen::VectorXd ret_vec(num_coeffs);

      //temporal container for polynomials
      Eigen::VectorXd p_container;



      for(int m = 0; m <= L; ++m){

        p_container = recursive_legendre_evaluation(L, m, x);

        //running variable
        int k = 0;

        for (int l = m; l <= L ; ++l){

          ret_vec(l*l + l + m) = p_container(k);
          if(m > 0){
            ret_vec(l*l + l - m) = std::pow(-1.,m)*factorial(l-m)/factorial(l+m)*ret_vec(l*l + l + m) ;
          }
          ++k;

         }

      }

      return ret_vec;

    }

    Eigen::MatrixXd recurse_all_legendre_polynomials_and_derivative(const int L, const double theta){

      //argument for legendre polynomial
      double x = std::cos(theta);

      //number of coefficients. We calculate one more to get the derivative
      int num_coeffs = (L+1)*(L+1);

      //return values
      Eigen::MatrixXd ret_mat(num_coeffs,2);
      ret_mat.setZero();

      //temporal container for polynomials
      Eigen::VectorXd p_container;

      //temporal container for derivatives
      Eigen::VectorXd dp_container;

      //Helper vector
      Eigen::VectorXd l_helper;

      //Helper unit array
      Eigen::ArrayXd unit_helper;

      //value of derivative for singularities
      Eigen::VectorXd sing_value;



      for(int m = 0; m <= L; ++m){

        //legendre polynomials
        p_container = recursive_legendre_evaluation(L+1, m, x);

        //array with l counting up
        l_helper = Eigen::VectorXd::LinSpaced(p_container.size(),m,L+1);


        //Deal with the singularity at \theta = 0
        if(fabs(theta) < 1e-8){
          if(fabs(m) == 1){
            sing_value = -0.5*l_helper.array()*(l_helper.array() + 1.);
          }
          else{
            sing_value = 0.*l_helper;
          }
          dp_container = sing_value.segment(0,l_helper.size()-1);
        }
        //Deal with the singularity at \theta = \pm \pi
        else if(fabs(theta)  > BEMBEL_PI - 1e-8){
          if(fabs(m) == 1){
            unit_helper = -1*Eigen::ArrayXd::Ones(l_helper.size());

            sing_value = 0.5*unit_helper.pow(l_helper.array() + 1.)*l_helper.array()*(l_helper.array() + 1.);

          }
          else{
            sing_value = 0.*l_helper;
          }
          dp_container = sing_value.segment(0,l_helper.size()-1);
        }
        else{

          dp_container = (l_helper.segment(0,p_container.size()-1).array() + 1. - m )*p_container.segment(1,p_container.size()-1).array()
                        - (l_helper.segment(0,p_container.size()-1).array() + 1.) * x * p_container.segment(0,p_container.size()-1).array();
          dp_container *= -1.*std::sin(theta)/(x*x - 1.);

        }

        //running variable
        int k = 0;


        for (int l = m; l <= L ; ++l){


          ret_mat(l*l + l + m,0) = p_container(k);
          ret_mat(l*l + l + m,1) = dp_container(k);

          if(m > 0){
            ret_mat(l*l + l - m,0) = std::pow(-1.,m)*factorial(l-m)/factorial(l+m)*ret_mat(l*l + l + m,0) ;
            ret_mat(l*l + l - m,1) = std::pow(-1.,m)*factorial(l-m)/factorial(l+m)*ret_mat(l*l + l + m,1) ;
          }

          ++k;

         }

      }


      return ret_mat;

    }


    //Here we make use of the recursive computation of Legendre polynomials. Use
    //ALWAYS this function. The function above was implemented to benchmark the
    //recursive stuff.
    Eigen::MatrixXcd evaluate_harmonic_expansion_fast(const int L, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      Eigen::MatrixXcd ret_mat;
      ret_mat.resize(num_coeffs,2);

      //running index
      int k = 0;

      //scaling factor
      double N_lm;

      //recusive computation of all polynomials and their first derivatives
      Eigen::MatrixXd p_lm = recurse_all_legendre_polynomials_and_derivative(L, theta);

      //slow implementation
      for(int l = 0; l <= L; ++l){

        for(int m = -l; m <= l ; ++m){


          N_lm = std::sqrt((2.*l+1.)/4./BEMBEL_PI*factorial(l-m)/factorial(l+m));

          //Y_lm
          ret_mat(k,0) = p_lm(k,0)*N_lm*(std::cos(m*phi) + I*std::sin(m*phi));

          //d_Y_lm_dt
          ret_mat(k,1) = p_lm(k,1)*N_lm*(std::cos(m*phi) + I*std::sin(m*phi));

          ++k;

        }
      }

      return ret_mat;
    }

    //Regular solid harmonics R_lm in spherical coordinates
    //We follow the definitions in
    //Yijun Liu. Fast Multipole Boundary Element Method: Theory and Applications in Engineering. Cambridge University Press, 2009. doi: 10.1017/CBO9780511605345.
    //rendering easy formulations for the shift relations
    Eigen::VectorXcd Rlm(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd R_lm(num_coeffs);
      R_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          fac = 1./factorial(l+m);

          R_lm(k) = fac*p_lm(k)*std::pow(r,l)*(std::cos(m*phi) + I*std::sin(m*phi));

          ++k;

        }
      }

      return R_lm;
    }

    //Regular solid harmonics R_lm in spherical coordinates
    //We follow the definitions in
    //Yijun Liu. Fast Multipole Boundary Element Method: Theory and Applications in Engineering. Cambridge University Press, 2009. doi: 10.1017/CBO9780511605345.
    //rendering easy formulations for the shift relations
    Eigen::VectorXcd Rlm(const int L, const Eigen::Vector3d &x){

      double r = x.norm();
      double theta = 0.;
      if(r > 1e-12){
        theta = std::acos(x(2)/r);
      }
      double phi = std::atan2(x(1),x(0));

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd R_lm(num_coeffs);
      R_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          fac = 1./factorial(l+m);

          R_lm(k) = fac*p_lm(k)*std::pow(r,l)*(std::cos(m*phi) + I*std::sin(m*phi));

          ++k;

        }
      }

      return R_lm;
    }

    //Regular solid harmonics R_lm in spherical coordinates
    //We follow the definitions in
    //Yijun Liu. Fast Multipole Boundary Element Method: Theory and Applications in Engineering. Cambridge University Press, 2009. doi: 10.1017/CBO9780511605345.
    //rendering easy formulations for the shift relations
    /*
    Eigen::VectorXcd R_lm_rpt(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      Eigen::MatrixXcd sph_harms =  evaluate_harmonic_expansion_fast(L,theta,phi);

      Eigen::VectorXcd R_rtp(num_coeffs);
      R_rtp.setZero();

      int k = 0;

      double sqrt_4pi = std::sqrt(4.*BEMBEL_PI);
      double fac;

      for(int l = 0; l <= L; ++l){

        fac = sqrt_4pi/std::sqrt(2.*l+1);

        for(int m = -l; m <= l ; ++m){

          R_rtp(k) = fac*sph_harms(k,0)*std::pow(r,l);

          ++k;

        }
      }

      return R_rtp;
    }
    */

    //Spatial derivatives of the solid harmonics R_lm in spherical coordinates
    Eigen::MatrixXcd Rlm_p_rpt(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      //Eigen::MatrixXcd sph_harms =  evaluate_harmonic_expansion_fast(L,theta,phi);
      //recusive computation of all polynomials and their first derivatives
      Eigen::MatrixXd p_lm = recurse_all_legendre_polynomials_and_derivative(L, theta);

      Eigen::MatrixXcd R_rtp(num_coeffs,3);
      R_rtp.setZero();


      //denominator for R_phi. Again we need to take care of singulatities.
      // From Hospitals rule, we can compute the limit for theta -> 0
      // lim = d_Y_lm_dt / cos(theta) |_theta = 0
      double den_Rp = std::sin(theta);
      int index_Rp = 0;

      if(fabs(theta) < 1e-8){
        den_Rp = 1.;
        index_Rp = 1;

      }
      else if(fabs(theta)  > BEMBEL_PI - 1e-8){
        den_Rp = -1.;
        index_Rp = 1;
      }

      int k = 1;

      //double sqrt_4pi = std::sqrt(4.*BEMBEL_PI);
      double fac;

      for(int l = 1; l <= L; ++l){

        //fac = sqrt_4pi/std::sqrt(2.*l+1);

        for(int m = -l; m <= l ; ++m){

          fac = 1./factorial(l+m);

          R_rtp(k,0) = fac*l*std::pow(r,l-1)*p_lm(k,0)*(std::cos(m*phi) + I*std::sin(m*phi));
          R_rtp(k,1) = fac*std::pow(r,l-1)*p_lm(k,1)*(std::cos(m*phi) + I*std::sin(m*phi));
          R_rtp(k,2) = fac*m*std::pow(r,l-1)/den_Rp*p_lm(k,index_Rp)*(-1.*std::sin(m*phi) + I*std::cos(m*phi));

          ++k;

        }
      }

      return R_rtp;


    }


    //Spatial derivatives of the solid harmonics R_lm
    Eigen::MatrixXcd Rlm_p(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);


      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      //Derivatives in spherical coordinates
      Eigen::MatrixXcd R_rtp = Rlm_p_rpt(L,r,theta,phi);

      Eigen::MatrixXcd ret_mat(num_coeffs,3);

      ret_mat.col(0) = std::sin(theta)*std::cos(phi)*R_rtp.col(0)
                            + std::cos(theta)*std::cos(phi)*R_rtp.col(1)
                            - std::sin(phi)*R_rtp.col(2);

      ret_mat.col(1) = std::sin(theta)*std::sin(phi)*R_rtp.col(0)
                            + std::cos(theta)*std::sin(phi)*R_rtp.col(1)
                            + std::cos(phi)*R_rtp.col(2);

      ret_mat.col(2) = std::cos(theta)*R_rtp.col(0) - std::sin(theta)*R_rtp.col(1);


      return ret_mat;


    }

    //Irregular solid harmonics S_lm in spherical coordinates
    Eigen::VectorXcd Slm(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd S_lm(num_coeffs);
      S_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          fac = factorial(l-m);

          S_lm(k) = fac*p_lm(k)*(std::cos(m*phi) + I*std::sin(m*phi))/std::pow(r,l+1);

          ++k;

        }
      }

      return S_lm;


    }

    //Irregular solid harmonics S_lm in spherical coordinates evaluated around y at x
    Eigen::VectorXcd Slm(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y){

      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      return Slm(L, r, theta, phi);


    }

    //Irregular solid harmonics S_lm in spherical coordinates at d
    Eigen::VectorXcd Slm(const int L, const Eigen::VectorXd &d){

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      return Slm(L, r, theta, phi);


    }

    //**************************************************************************
    //  Convention of A new version of the Fast Multipole
    //                Method for the Laplace equation in three
    //                dimensions
    //
    // Be careful, the definition of the spherical harmoics is different here!
    // We denote this "alternative" definition by the suffic _alt
    //**************************************************************************

    // A spherica harmonic in the old definition relates to this defintion by
    /*
      \tilde{Y}_n^m(\theta,\varphi) = \begin{cases}
        (-1)^m \sqrt{\dfrac{(n-m)!}{(n+m)!}}P_n^m(\cos(\theta))\exp(jm\varphi), \qquad m < 0
        \sqrt{\dfrac{(n-m)!}{(n+m)!}}P_n^m(\cos(\theta))\exp(jm\varphi), \qquad m \leq 0
      \end{cases}
    */

    std::complex<double> Ynm_alt(const double theta,const double phi,const int n, const int m){

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      double pnm = associated_legendre(n, std::abs(m) , std::cos(theta));

      return std::sqrt(factorial(n - std::abs(m)) / factorial(n + std::abs(m))) * pnm * (std::cos(m*phi) + I*std::sin(m*phi));

    }


    //Regular solid harmonics R_lm in spherical coordinates
    //We follow the definitions in
    //A new version of the Fast Multipole  Method for the Laplace equation in three dimensions
    //rendering easy formulations for the shift relations
    Eigen::VectorXcd Rlm_alt(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd R_lm(num_coeffs);
      R_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          if(m < 0){
            fac = std::pow(-1.,m)*std::sqrt(factorial(l-m)/factorial(l+m));
          }
          else{
            fac = std::sqrt(factorial(l-m)/factorial(l+m));
          }


          R_lm(k) = fac*p_lm(k)*std::pow(r,l)*(std::cos(m*phi) + I*std::sin(m*phi));

          ++k;

        }
      }

      return R_lm;
    }

    //Spatial derivatives of the solid harmonics R_lm in spherical coordinates
    // this function returns the harmonics in reverse order! Such that the matrix
    // can directly be used for multipole expansions
    Eigen::MatrixXcd Rlm_p_rpt_alt(const int L, const double r, const double theta, const double phi, const bool reverse_m = false){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      //Eigen::MatrixXcd sph_harms =  evaluate_harmonic_expansion_fast(L,theta,phi);
      //recusive computation of all polynomials and their first derivatives
      Eigen::MatrixXd p_lm = recurse_all_legendre_polynomials_and_derivative(L, theta);

      Eigen::MatrixXcd R_rtp(num_coeffs,3);
      R_rtp.setZero();


      //denominator for R_phi. Again we need to take care of singulatities.
      // From Hospitals rule, we can compute the limit for theta -> 0
      // lim = d_Y_lm_dt / cos(theta) |_theta = 0
      double den_Rp = std::sin(theta);
      int index_Rp = 0;

      if(fabs(theta) < 1e-8){
        den_Rp = 1.;
        index_Rp = 1;

      }
      else if(fabs(theta)  > BEMBEL_PI - 1e-8){
        den_Rp = -1.;
        index_Rp = 1;
      }

      int k = 1;

      //double sqrt_4pi = std::sqrt(4.*BEMBEL_PI);
      double fac;

      if(reverse_m){
        for(int l = 1; l <= L; ++l){

          //fac = sqrt_4pi/std::sqrt(2.*l+1);

          for(int m = l; m >= -l ; --m){

            if(m < 0){
              fac = std::pow(-1.,m)*std::sqrt(factorial(l-m)/factorial(l+m));
            }
            else{
              fac = std::sqrt(factorial(l-m)/factorial(l+m));
            }


            R_rtp(k,0) = fac*l*std::pow(r,l-1)*p_lm(k,0)*(std::cos(m*phi) + I*std::sin(m*phi));
            R_rtp(k,1) = fac*std::pow(r,l-1)*p_lm(k,1)*(std::cos(m*phi) + I*std::sin(m*phi));
            R_rtp(k,2) = fac*m*std::pow(r,l-1)/den_Rp*p_lm(k,index_Rp)*(-1.*std::sin(m*phi) + I*std::cos(m*phi));

            ++k;

          }
        }
      }
      else{
        for(int l = 1; l <= L; ++l){

          //fac = sqrt_4pi/std::sqrt(2.*l+1);

          for(int m = -l; m <= l ; ++m){

            if(m < 0){
              fac = std::pow(-1.,m)*std::sqrt(factorial(l-m)/factorial(l+m));
            }
            else{
              fac = std::sqrt(factorial(l-m)/factorial(l+m));
            }


            R_rtp(k,0) = fac*l*std::pow(r,l-1)*p_lm(k,0)*(std::cos(m*phi) + I*std::sin(m*phi));
            R_rtp(k,1) = fac*std::pow(r,l-1)*p_lm(k,1)*(std::cos(m*phi) + I*std::sin(m*phi));
            R_rtp(k,2) = fac*m*std::pow(r,l-1)/den_Rp*p_lm(k,index_Rp)*(-1.*std::sin(m*phi) + I*std::cos(m*phi));

            ++k;

          }
        }
      }

      return R_rtp;


    }


    //Spatial derivatives of the solid harmonics R_lm
    Eigen::MatrixXcd Rlm_p_alt(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y, const bool reverse_m = false){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);


      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = 0.;
      if (r > 1e-10){
        theta = std::acos(d(2)/r);
      }
      double phi = std::atan2(d(1),d(0));

      //Derivatives in spherical coordinates
      Eigen::MatrixXcd R_rtp = Rlm_p_rpt_alt(L,r,theta,phi,reverse_m);

      Eigen::MatrixXcd ret_mat(num_coeffs,3);

      ret_mat.col(0) = std::sin(theta)*std::cos(phi)*R_rtp.col(0)
                            + std::cos(theta)*std::cos(phi)*R_rtp.col(1)
                            - std::sin(phi)*R_rtp.col(2);

      ret_mat.col(1) = std::sin(theta)*std::sin(phi)*R_rtp.col(0)
                            + std::cos(theta)*std::sin(phi)*R_rtp.col(1)
                            + std::cos(phi)*R_rtp.col(2);

      ret_mat.col(2) = std::cos(theta)*R_rtp.col(0) - std::sin(theta)*R_rtp.col(1);


      return ret_mat;


    }


    //Regular solid harmonics R_lm in spherical coordinates
    Eigen::VectorXcd Rlm_alt(const int L, const Eigen::Vector3d &x, bool print = false){

      double r = x.norm();
      double theta = 0.;
      if(r > 1e-12){
        theta = std::acos(x(2)/r);
      }
      double phi = std::atan2(x(1),x(0));

      if(print){

        std::cout << "x = " << x.transpose() << std::endl;
        std::cout << "theta = " << theta/BEMBEL_PI << " pi" << std::endl;
        std::cout << "phi = " << phi/BEMBEL_PI << " pi" << std::endl;
      }

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd R_lm(num_coeffs);
      R_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          if(m < 0){
            fac = std::pow(-1.,m)*std::sqrt(factorial(l-m)/factorial(l+m));
          }
          else{
            fac = std::sqrt(factorial(l-m)/factorial(l+m));
          }

          R_lm(k) = fac*p_lm(k)*std::pow(r,l)*(std::cos(m*phi) + I*std::sin(m*phi));

          ++k;

        }
      }

      return R_lm;
    }

    //Irregular solid harmonics S_lm in spherical coordinates
    Eigen::VectorXcd Slm_alt(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      Eigen::VectorXd p_lm = recurse_all_legendre_polynomials( L, std::cos(theta) );

      Eigen::VectorXcd S_lm(num_coeffs);
      S_lm.setZero();

      int k = 0;
      double fac;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          if(m < 0){
            fac = std::pow(-1.,m)*std::sqrt(factorial(l-m)/factorial(l+m));
          }
          else{
            fac = std::sqrt(factorial(l-m)/factorial(l+m));
          }


          S_lm(k) = fac*p_lm(k)*(std::cos(m*phi) + I*std::sin(m*phi))/std::pow(r,l+1);

          ++k;

        }
      }

      return S_lm;


    }

    //Irregular solid harmonics S_lm in spherical coordinates evaluated around y at x
    Eigen::VectorXcd Slm_alt(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y){

      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      return Slm_alt(L, r, theta, phi);


    }

    double A_nm(const int n, const int m){

      return std::pow(-1.,n)/std::sqrt(factorial(n-m)*factorial(n+m));
    }

    double wigner_d_small(const double theta,const int j,const int m, const int mp){

      double fac1 = std::sqrt( factorial(j+m)*factorial(j-m)*factorial(j+mp)*factorial(j-mp) );

      int s_min = 0;
      if((m-mp) > s_min) s_min = m-mp;

      int s_max = j+m;
      if((j-mp) < s_max) s_max = j-mp;

      //std::cout << "\tfac1 = " << fac1 << std::endl;
      //std::cout << "\ts_min = " << s_min << std::endl;
      //std::cout << "\ts_max = " << s_max << std::endl;
      //std::cout << "\tj = " << j << std::endl;
      //std::cout << "\tm = " << m << std::endl;
      //std::cout << "\tmp = " << mp << std::endl;


      double tmp_sum = 0.;

      for(int s = s_min; s <= s_max; ++s){

        double num = std::pow(-1.,mp-m+s) * std::pow(std::cos(0.5*theta),2*j+m-mp-2*s) * std::pow(std::sin(0.5*theta),mp-m+2*s);
        double den = factorial(j+m-s)*factorial(s)*factorial(mp-m+s)*factorial(j-mp-s);

        tmp_sum += num/den;
      }

      return fac1*tmp_sum;
    }

    //correction factor. Needed as we are using the Wigner coefficients based on the normal spherical harmonics
    double correction_factor(const int n, const int m){

        double ret_val = std::sqrt(4.*BEMBEL_PI/(2.*n + 1));

        if(m < 0) ret_val *= std::pow(-1.,m);

        return ret_val;
    }

    Eigen::SparseMatrix<std::complex<double>> make_rotation_matrix(const int num_multipoles,const double theta,const double phi){

      typedef Eigen::Triplet<std::complex<double>> T;
      std::vector<T> tripletList;
      tripletList.reserve(num_multipoles*num_multipoles);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);

      int l_cum = 0;

      for(int l = 0; l <= num_multipoles; ++l){
        for(int k = -l; k <= l; ++k){
          for(int m = -l; m <= l; ++m){

            double d = wigner_d_small(theta,l,m,k);
            std::complex<double> e = std::cos(k*phi) + I*std::sin(k*phi);

            //Needed as we are using the Wigner coefficients based on the normal spherical harmonics
            double f = correction_factor(l,m)/correction_factor(l,k);


            tripletList.push_back(T(l_cum + m + l, l_cum + k + l , f*d*e ));

          }
        }

        l_cum += 2*l + 1;
      }

      int num_coeffs = (num_multipoles+1)*(num_multipoles+1);

      Eigen::SparseMatrix<std::complex<double>> mat(num_coeffs,num_coeffs);

      mat.setFromTriplets(tripletList.begin(), tripletList.end());

      return mat;
    }

    Eigen::VectorXd evaluate_multipole_expansion(const Eigen::MatrixXd &pos, const Eigen::Vector3d &center, const Eigen::VectorXcd &multipoles){

      //number of positions to evaluate
      int num_pos = pos.rows();

      //number of multipoles
      int L = std::sqrt(multipoles.rows())-1;

      //potential container
      Eigen::VectorXd pot(num_pos);

      //number of positions to evaluate
      for(int i = 0 ; i < num_pos ; ++i){

        Eigen::VectorXcd S_lm = Slm_alt(L, pos.row(i), center);

        std::complex<double> tmp = S_lm.transpose() * multipoles;

        pot(i) = tmp.real();

      }

      return pot;

    }

    Eigen::VectorXd evaluate_local_expansion(const Eigen::MatrixXd &pos, const Eigen::Vector3d &center, const Eigen::VectorXcd &locals, bool print = false){

      //number of positions to evaluate
      int num_pos = pos.rows();

      //number of multipoles
      int L = std::sqrt(locals.rows())-1;

      //potential container
      Eigen::VectorXd pot(num_pos);

      //number of positions to evaluate
      for(int i = 0 ; i < num_pos ; ++i){

        Eigen::VectorXcd R_lm;
        if((i == 0) && (print == true)){

          R_lm = Rlm_alt(L, pos.row(i).transpose() - center, print);
        }
        else{
          R_lm = Rlm_alt(L, pos.row(i).transpose() - center);
        }


        std::complex<double> tmp = R_lm.transpose() * locals;

        pot(i) = tmp.real();

      }

      return pot;

    }

    //Spatial derivatives of the irregular solid harmonics S_lm in spherical coordinates
    /*
    Eigen::MatrixXcd S_lm_p_rpt(const int L, const double r, const double theta, const double phi){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);

      //imaginary number as helper
      std::complex<double> I(0.0,1.0);


      Eigen::MatrixXcd sph_harms =  evaluate_harmonic_expansion_fast(L,theta,phi);

      Eigen::MatrixXcd S_rtp(num_coeffs,3);
      S_rtp.setZero();


      //denumerator for S_phi. Again we need to take care of singulatities.
      // From Hospitals rule, we can compute the limit for theta -> 0
      // lim = d_Y_lm_dt / cos(theta) |_theta = 0
      double den_Sp = std::sin(theta);
      int index_Sp = 0;

      if(fabs(theta) < 1e-8){
        den_Sp = 1.;
        index_Sp = 1;

      }
      else if(fabs(theta)  > BEMBEL_PI - 1e-8){
        den_Sp = -1.;
        index_Sp = 1;
      }

      int k = 0;

      for(int l = 0; l <= L; ++l){
        for(int m = -l; m <= l ; ++m){

          S_rtp(k,0) = -(l+1.)*sph_harms(k,0)/std::pow(r,l+2);
          S_rtp(k,1) = sph_harms(k,1)/std::pow(r,l+2);
          S_rtp(k,2) = 1.*m*I*sph_harms(k,index_Sp)/std::pow(r,l+2)/den_Sp;

          ++k;

        }
      }

      return S_rtp;


    }


    //Spatial derivatives of the solid harmonics R_lm
    Eigen::MatrixXcd S_lm_p(const int L, const Eigen::Vector3d &x, const Eigen::Vector3d &y){

      //number of coefficients
      int num_coeffs = (L+1)*(L+1);


      Eigen::VectorXd d = x-y;

      double r = d.norm();
      double theta = std::acos(d(2)/r);
      double phi = std::atan2(d(1),d(0));

      //Derivatives in spherical coordinates
      Eigen::MatrixXcd S_rtp = S_lm_p_rpt(L,r,theta,phi);

      Eigen::MatrixXcd ret_mat(num_coeffs,3);

      ret_mat.col(0) = std::sin(theta)*std::cos(phi)*S_rtp.col(0)
                            + std::cos(theta)*std::cos(phi)*S_rtp.col(1)
                            - std::sin(phi)*S_rtp.col(2);

      ret_mat.col(1) = std::sin(theta)*std::sin(phi)*S_rtp.col(0)
                            + std::cos(theta)*std::sin(phi)*S_rtp.col(1)
                            + std::cos(phi)*S_rtp.col(2);

      ret_mat.col(2) = std::cos(theta)*S_rtp.col(0) - std::sin(theta)*S_rtp.col(1);


      return ret_mat;


    }
    */
    /*
    Eigen::MatrixXcd legendre_recursive_evaluation(const int m, const double x){


      Eigen::MatrixXd ret_mat;
      ret_mat.resize(l+1,2);

      if(fabs(x) < 0){

        for(int i = 0; (i < 2) && (i =< l) ; ++i){
          ret_mat(i,0) = associated_legendre(l , 0, x);

        }



      }


      if(l < 2){


      }
      else{


      }

      for(int m = 0; m <= l; ++l){
-

      }

      return ret_mat;


    }
    */

  }  // namespace Bembel



#endif
