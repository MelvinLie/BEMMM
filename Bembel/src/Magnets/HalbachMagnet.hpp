// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
//
#ifndef BEMBEL_HALBACH_MAGNET_H_
#define BEMBEL_HALBACH_MAGNET_H_

#include <fstream>

#include <Eigen/Dense>



namespace Bembel {

/**
 * \ingroup Laplace
 */
class HalbachMagnet{
  // implementation of the kernel evaluation, which may be based on the
  // information available from the superSpace
 public:
  HalbachMagnet() {

    quad_deg_ = 50;

  }

  void read_geometry_file(const std::string &file_name, const char del = ',') {

    std::ifstream file;
    int numLines = 0;
    std::string current_line;
    std::string current_element;

    std::vector<std::vector<double>> tmp_data;

    file.open(file_name);

    if (!file) {
      std::cerr << "File " << file_name << " doesn't exist!";
      exit(1);
    }

    // first row is irrelevant
    getline(file, current_line);

    while ( std::getline(file, current_line) ){

      std::stringstream current_data(current_line);
      std::vector<double> tmp_vector;


      while(std::getline(current_data,current_element,del)){
        tmp_vector.push_back(atof(current_element.c_str()));

      }
      tmp_data.push_back(tmp_vector);



    }

    file.close();

    num_magnets_ = tmp_data.size();

    p0_.resize(num_magnets_,3);
    p1_.resize(num_magnets_,3);
    p2_.resize(num_magnets_,3);

    eu_.resize(num_magnets_,3);
    ev_.resize(num_magnets_,3);
    ew_.resize(num_magnets_,3);

    h_.resize(num_magnets_);
    a_.resize(num_magnets_);
    b_.resize(num_magnets_);
    l_.resize(num_magnets_);

    H_.resize(num_magnets_,3);

    for (int m = 0; m < num_magnets_; ++m){

        p0_(m,0) = tmp_data[m][0];
        p0_(m,1) = tmp_data[m][1];
        p0_(m,2) = tmp_data[m][2];

        p1_(m,0) = tmp_data[m][3];
        p1_(m,1) = tmp_data[m][4];
        p1_(m,2) = tmp_data[m][5];

        p2_(m,0) = tmp_data[m][6];
        p2_(m,1) = tmp_data[m][7];
        p2_(m,2) = tmp_data[m][8];


        b_(m) = (p1_.row(m) - p0_.row(m)).norm();
        h_(m) = (p2_.row(m) - p0_.row(m)).norm();
        a_(m) = tmp_data[m][9] * b_(m);
        l_(m) = tmp_data[m][10];


        eu_.row(m) = (p1_.row(m) - p0_.row(m))/b_(m);
        ev_.row(m) = (p2_.row(m) - p0_.row(m))/h_(m);

        //cross product
        ew_(m,0) = eu_(m,1)*ev_(m,2) - eu_(m,2)*ev_(m,1);
        ew_(m,1) = eu_(m,2)*ev_(m,0) - eu_(m,0)*ev_(m,2);
        ew_(m,2) = eu_(m,0)*ev_(m,1) - eu_(m,1)*ev_(m,0);

        H_(m,0) = tmp_data[m][11];
        H_(m,1) = tmp_data[m][12];
        H_(m,2) = tmp_data[m][13];



      }
    }

   void print_magnet_info(){

     std::cout << "Halbach Magnet:" << std::endl;
     std::cout << "\tP0_X\tP0_Y\tP0_Z\teu_x\teu_y\teu_z\tev_x\tev_y\tev_z\tew_x\tew_y\tew_z\ta\tb\th\tl \tH_X\tH_Y\tH_Z" << std::endl;
     for (int m = 0; m < num_magnets_; ++m){

        std::cout <<"\t" << p0_(m,0) << "\t";
        std::cout << p0_(m,1) << "\t";
        std::cout << p0_(m,2) << "\t";

        std::cout << eu_(m,0) << "\t";
        std::cout << eu_(m,1) << "\t";
        std::cout << eu_(m,2) << "\t";

        std::cout << ev_(m,0) << "\t";
        std::cout << ev_(m,1) << "\t";
        std::cout << ev_(m,2) << "\t";

        std::cout << ew_(m,0) << "\t";
        std::cout << ew_(m,1) << "\t";
        std::cout << ew_(m,2) << "\t";

        std::cout << a_(m) << "\t";
        std::cout << b_(m) << "\t";
        std::cout << h_(m) << "\t";
        std::cout << l_(m) << "\t";

        std::cout << H_(m,0) << "\t";
        std::cout << H_(m,1) << "\t";
        std::cout << H_(m,2) << std::endl;

       }

   }


   Eigen::MatrixXd assemble_B_matrix(Eigen::MatrixXd r){

     Eigen::MatrixXd B;
     double A;
     Eigen::Vector3d n;
     Eigen::Vector3d tmp_B, tmp_r, pos_0, pos_1;
     Eigen::Vector3d e_j, e_k;
     Eigen::Vector3d face_ctr;
     Eigen::Vector2d xi;
     double w;
     double xp,yp;
     double gamma;
     double d_k,d_j;
     int quad_deg;

     double mu_0 = 4*BEMBEL_PI*1e-7;

     GaussSquare<Constants::maximum_quadrature_degree> GS;

     //std::cout << "Max quad degree = " << Constants::maximum_quadrature_degree << std::endl;


     int num_meas = r.rows();


     B.resize(3*num_meas,3*num_magnets_);
     B.setZero();

     for(int im = 0; im < num_magnets_; ++im){


       //we integrate over the six faces of the REC magnet

       //----------------------------------
       //face 1
       n = -1*ew_.row(im);
       pos_0 = p0_.row(im);

       face_ctr = pos_0.transpose() +  0.5*h_(im)*ev_.row(im);



       for(int ir = 0; ir < num_meas; ++ir){


         tmp_B(0) = 0;
         tmp_B(1) = 0;
         tmp_B(2) = 0;

         quad_deg = get_quad_degree(r.row(ir), face_ctr, 2*b_(im));
         auto Q = GS[quad_deg];

         //quadrature
         for (auto i = 0; i < Q.w_.size(); ++i) {

           //xi = (u , v)
           xi = Q.xi_.col(i);

           //coordinates on trapez
           xp = (2*xi(0)-1)*((a_(im)-b_(im))*xi(1)+b_(im));
           yp = h_(im)*xi(1);

           //surface element
           gamma = 2*h_(im)*((a_(im)-b_(im))*xi(1) + b_(im));

           //integration point
           tmp_r = pos_0.transpose() + xp*eu_.row(im) + yp*ev_.row(im);

           //B-field
           tmp_B += gamma*Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

         }

         //Bx

         B(ir,im)                            += n(0)*tmp_B(0);  //Bx
         B(ir,im+num_magnets_)               += n(1)*tmp_B(0);  //By
         B(ir,im+2*num_magnets_)             += n(2)*tmp_B(0);  //Bz

         //By
         B(ir+num_meas,im)                  += n(0)*tmp_B(1);  //Bx
         B(ir+num_meas,im+num_magnets_)     += n(1)*tmp_B(1);  //By
         B(ir+num_meas,im+2*num_magnets_)   += n(2)*tmp_B(1);  //Bz

         //Bz
         B(ir+2*num_meas,im)                += n(0)*tmp_B(2);  //Bx
         B(ir+2*num_meas,im+num_magnets_)   += n(1)*tmp_B(2);  //By
         B(ir+2*num_meas,im+2*num_magnets_) += n(2)*tmp_B(2);  //Bz


       }

         //----------------------------------
         //face 2

         n = ew_.row(im);
         pos_0 = p0_.row(im) + l_(im)*ew_.row(im);

         face_ctr = pos_0.transpose() +  0.5*h_(im)*ev_.row(im);

         for(int ir = 0; ir < num_meas; ++ir){

           tmp_B(0) = 0;
           tmp_B(1) = 0;
           tmp_B(2) = 0;

           quad_deg = get_quad_degree(r.row(ir), face_ctr, 2*b_(im));
           auto Q = GS[quad_deg];

           //quadrature
           for (auto i = 0; i < Q.w_.size(); ++i) {

             //xi = (u , v)
             xi = Q.xi_.col(i);

             //coordinates on trapez
             xp = (2*xi(0)-1)*((a_(im)-b_(im))*xi(1)+b_(im));
             yp = h_(im)*xi(1);

             //surface element
             gamma = 2*h_(im)*((a_(im)-b_(im))*xi(1) + b_(im));

             //integration point
             tmp_r = pos_0.transpose() + xp*eu_.row(im) + yp*ev_.row(im);

             //B-field
             tmp_B += gamma*Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

           }



           //Bx

           B(ir,im)                            += n(0)*tmp_B(0);  //Bx
           B(ir,im+num_magnets_)               += n(1)*tmp_B(0);  //By
           B(ir,im+2*num_magnets_)             += n(2)*tmp_B(0);  //Bz

           //By
           B(ir+num_meas,im)                  += n(0)*tmp_B(1);  //Bx
           B(ir+num_meas,im+num_magnets_)     += n(1)*tmp_B(1);  //By
           B(ir+num_meas,im+2*num_magnets_)   += n(2)*tmp_B(1);  //Bz

           //Bz
           B(ir+2*num_meas,im)                += n(0)*tmp_B(2);  //Bx
           B(ir+2*num_meas,im+num_magnets_)   += n(1)*tmp_B(2);  //By
           B(ir+2*num_meas,im+2*num_magnets_) += n(2)*tmp_B(2);  //Bz

         }

       //----------------------------------
       //face 3
       pos_0 = p0_.row(im) - b_(im)*eu_.row(im);
       pos_1 = p2_.row(im) - a_(im)*eu_.row(im);

       e_j = (pos_1 - pos_0);
       d_j = e_j.norm();

       A = d_j*l_(im);

       e_j /= d_j;

       face_ctr(0) = pos_0(0) +  0.5*d_j*e_j(0) + 0.5*l_(im) * ew_(im,0);
       face_ctr(1) = pos_0(1) +  0.5*d_j*e_j(1) + 0.5*l_(im) * ew_(im,1);
       face_ctr(2) = pos_0(2) +  0.5*d_j*e_j(2) + 0.5*l_(im) * ew_(im,2);

       //cross product
       n(0) = -e_j(1)*ew_(im,2) + e_j(2)*ew_(im,1);
       n(1) = -e_j(2)*ew_(im,0) + e_j(0)*ew_(im,2);
       n(2) = -e_j(0)*ew_(im,1) + e_j(1)*ew_(im,0);


       for(int ir = 0; ir < num_meas; ++ir){

         tmp_B(0) = 0;
         tmp_B(1) = 0;
         tmp_B(2) = 0;

         quad_deg = get_quad_degree(r.row(ir), face_ctr, std::max(d_j,l_(im)));
         auto Q = GS[quad_deg];

         //quadrature
         for (auto i = 0; i < Q.w_.size(); ++i) {

           //xi = (u , v)
           xi = Q.xi_.col(i);

           //integration point
           tmp_r = pos_0 + l_(im)*xi(0)*ew_.row(im).transpose()
                      + d_j*xi(1)*e_j ;

           //B-field
           tmp_B += Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

         }


         //Bx
         B(ir,im)                            += A*n(0)*tmp_B(0);  //Bx
         B(ir,im+num_magnets_)               += A*n(1)*tmp_B(0);  //By
         B(ir,im+2*num_magnets_)             += A*n(2)*tmp_B(0);  //Bz

         //By
         B(ir+num_meas,im)                  += A*n(0)*tmp_B(1);  //Bx
         B(ir+num_meas,im+num_magnets_)     += A*n(1)*tmp_B(1);  //By
         B(ir+num_meas,im+2*num_magnets_)   += A*n(2)*tmp_B(1);  //Bz

         //Bz
         B(ir+2*num_meas,im)                += A*n(0)*tmp_B(2);  //Bx
         B(ir+2*num_meas,im+num_magnets_)   += A*n(1)*tmp_B(2);  //By
         B(ir+2*num_meas,im+2*num_magnets_) += A*n(2)*tmp_B(2);  //Bz

       }


       //----------------------------------
       //face 4
       pos_0 = p1_.row(im);
       pos_1 = p2_.row(im) + a_(im)*eu_.row(im);

       e_k = (pos_1 - pos_0);
       d_k = e_k.norm();

       A = d_k*l_(im);

       e_k /= d_k;

       face_ctr(0) = pos_0(0) +  0.5*d_k*e_k(0) + 0.5*l_(im) * ew_(im,0);
       face_ctr(1) = pos_0(1) +  0.5*d_k*e_k(1) + 0.5*l_(im) * ew_(im,1);
       face_ctr(2) = pos_0(2) +  0.5*d_k*e_k(2) + 0.5*l_(im) * ew_(im,2);

       //cross product
       n(0) = e_k(1)*ew_(im,2) - e_k(2)*ew_(im,1);
       n(1) = e_k(2)*ew_(im,0) - e_k(0)*ew_(im,2);
       n(2) = e_k(0)*ew_(im,1) - e_k(1)*ew_(im,0);


       for(int ir = 0; ir < num_meas; ++ir){

         tmp_B(0) = 0;
         tmp_B(1) = 0;
         tmp_B(2) = 0;

         quad_deg = get_quad_degree(r.row(ir), face_ctr, std::max(d_k,l_(im)));
         auto Q = GS[quad_deg];

         //quadrature
         for (auto i = 0; i < Q.w_.size(); ++i) {

           //xi = (u , v)
           xi = Q.xi_.col(i);

           //integration point
           tmp_r = pos_0 + d_k*xi(0)*e_k
                    + l_(im)*xi(1)*ew_.row(im).transpose();

           //B-field
           tmp_B += Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

         }

         //Bx

         B(ir,im)                            += A*n(0)*tmp_B(0);  //Bx
         B(ir,im+num_magnets_)               += A*n(1)*tmp_B(0);  //By
         B(ir,im+2*num_magnets_)             += A*n(2)*tmp_B(0);  //Bz

         //By
         B(ir+num_meas,im)                  += A*n(0)*tmp_B(1);  //Bx
         B(ir+num_meas,im+num_magnets_)     += A*n(1)*tmp_B(1);  //By
         B(ir+num_meas,im+2*num_magnets_)   += A*n(2)*tmp_B(1);  //Bz

         //Bz
         B(ir+2*num_meas,im)                += A*n(0)*tmp_B(2);  //Bx
         B(ir+2*num_meas,im+num_magnets_)   += A*n(1)*tmp_B(2);  //By
         B(ir+2*num_meas,im+2*num_magnets_) += A*n(2)*tmp_B(2);  //Bz

       }


       //----------------------------------
       //face 5
       pos_0 = p0_.row(im) - b_(im)*eu_.row(im);



       face_ctr = p0_.row(im) + 0.5*l_(im) * ew_.row(im);

       A = 2*l_(im)*b_(im);


       n = -1*ev_.row(im);


       for(int ir = 0; ir < num_meas; ++ir){

         tmp_B(0) = 0;
         tmp_B(1) = 0;
         tmp_B(2) = 0;

         quad_deg = get_quad_degree(r.row(ir), face_ctr, std::max(2*b_(im),l_(im)));
         auto Q = GS[quad_deg];

         //quadrature
         for (auto i = 0; i < Q.w_.size(); ++i) {

           //xi = (u , v)
           xi = Q.xi_.col(i);

           //integration point
           tmp_r = pos_0 + 2*b_(im)*xi(0)*eu_.row(im).transpose()
                          +  l_(im)*xi(1)*ew_.row(im).transpose();

           //B-field
           tmp_B += Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

         }

         //Bx
         B(ir,im)                            += A*n(0)*tmp_B(0);  //Bx
         B(ir,im+num_magnets_)               += A*n(1)*tmp_B(0);  //By
         B(ir,im+2*num_magnets_)             += A*n(2)*tmp_B(0);  //Bz

         //By
         B(ir+num_meas,im)                  += A*n(0)*tmp_B(1);  //Bx
         B(ir+num_meas,im+num_magnets_)     += A*n(1)*tmp_B(1);  //By
         B(ir+num_meas,im+2*num_magnets_)   += A*n(2)*tmp_B(1);  //Bz

         //Bz
         B(ir+2*num_meas,im)                += A*n(0)*tmp_B(2);  //Bx
         B(ir+2*num_meas,im+num_magnets_)   += A*n(1)*tmp_B(2);  //By
         B(ir+2*num_meas,im+2*num_magnets_) += A*n(2)*tmp_B(2);  //Bz


       }

       //----------------------------------
       //face 6
       pos_0 = p2_.row(im) - a_(im)*eu_.row(im);


       A = 2*l_(im)*a_(im);

       //cross product
       n = ev_.row(im);

       face_ctr = p2_.row(im) + 0.5*l_(im) * ew_.row(im);

       for(int ir = 0; ir < num_meas; ++ir){

         tmp_B(0) = 0;
         tmp_B(1) = 0;
         tmp_B(2) = 0;

         quad_deg = get_quad_degree(r.row(ir), face_ctr, std::max(2*a_(im),l_(im)));
         auto Q = GS[quad_deg];

         //quadrature
         for (auto i = 0; i < Q.w_.size(); ++i) {

           //xi = (u , v)
           xi = Q.xi_.col(i);

           //integration point
           tmp_r = pos_0 + 2*a_(im)*xi(1)*eu_.row(im).transpose()
                        + l_(im)*xi(0)*ew_.row(im).transpose();


           //B-field
           tmp_B += Q.w_(i)*grad_1_over_r(r.row(ir),tmp_r);

         }

         //Bx
         B(ir,im)                            += A*n(0)*tmp_B(0);  //Bx
         B(ir,im+num_magnets_)               += A*n(1)*tmp_B(0);  //By
         B(ir,im+2*num_magnets_)             += A*n(2)*tmp_B(0);  //Bz

         //By
         B(ir+num_meas,im)                  += A*n(0)*tmp_B(1);  //Bx
         B(ir+num_meas,im+num_magnets_)     += A*n(1)*tmp_B(1);  //By
         B(ir+num_meas,im+2*num_magnets_)   += A*n(2)*tmp_B(1);  //Bz

         //Bz
         B(ir+2*num_meas,im)                += A*n(0)*tmp_B(2);  //Bx
         B(ir+2*num_meas,im+num_magnets_)   += A*n(1)*tmp_B(2);  //By
         B(ir+2*num_meas,im+2*num_magnets_) += A*n(2)*tmp_B(2);  //Bz

       }


     }

    return B/4./BEMBEL_PI;

   }

   Eigen::VectorXd flatten_H_matrix(Eigen::MatrixXd H_mat){

     Eigen::VectorXd ret_vec;

     int num_mags = H_mat.rows();

     ret_vec.resize(3*num_mags);

     ret_vec.segment(0,num_mags)          = H_mat.col(0);
     ret_vec.segment(num_mags,num_mags)   = H_mat.col(1);
     ret_vec.segment(2*num_mags,num_mags) = H_mat.col(2);

     return ret_vec;

   }

   Eigen::MatrixXd get_H_matrix(){

     return H_;
   }

private:
  Eigen::MatrixXd p0_,p1_,p2_;
  Eigen::MatrixXd eu_,ev_,ew_;
  Eigen::VectorXd a_,b_,h_,l_;
  Eigen::MatrixXd H_;
  int num_magnets_;
  int quad_deg_;
  int max_quad_degree_ = 12;

  Eigen::Vector3d grad_1_over_r(const Eigen::Vector3d &r,const Eigen::Vector3d &r_p){

    Eigen::Vector3d d = r - r_p;

    return d/std::pow(d.norm(),3);

  }

  int get_quad_degree(Eigen::Vector3d r, Eigen::Vector3d face_ctr, double diam){

    int quad_deg = 10;

    double dist = (r-face_ctr).norm();

    if(dist/diam < 0.75) quad_deg = 20;
    if(dist/diam < 0.5) quad_deg = 30;
    if(dist/diam < 0.25) quad_deg = Constants::maximum_quadrature_degree;


    //return quad_deg;
    return 50;

  }


};

}  // namespace Bembel

#endif
