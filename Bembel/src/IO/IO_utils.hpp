// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universitaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef BEMBEL_IO_UTILS_H_
#define BEMBEL_IO_UTILS_H_

#include <iostream>
#include <fstream>

#include <sys/sysinfo.h>
#include <unistd.h>


namespace Bembel {


  Eigen::MatrixXd read_csv(const std::string &file_name,bool skip_header = false, const char del = ','){

    std::ifstream file;
    int numLines = 0;
    std::string current_line;
    std::string current_element;
    Eigen::MatrixXd ret_val;


    std::vector<std::vector<double>> tmp_data;

    file.open(file_name);

    if (!file) {
      std::cerr << "File " << file_name << " doesn't exist!";
      exit(1);
    }

    if (skip_header)  getline(file, current_line);

    while ( std::getline(file, current_line) ){

      std::stringstream current_data(current_line);
      std::vector<double> tmp_vector;


      while(std::getline(current_data,current_element,del)){
        tmp_vector.push_back(atof(current_element.c_str()));

      }
      tmp_data.push_back(tmp_vector);



    }

    numLines = tmp_data.size();
    int num_cols = tmp_data[0].size();

    ret_val.resize(numLines,num_cols);


    for (int l = 0; l < numLines; ++l){
        for(int j = 0 ; j < num_cols ; ++j){
            ret_val(l,j) = tmp_data[l][j];
        }

      }

      return ret_val;

  }

  std::vector<std::string> read_input_ouput_file(const std::string &file_name){

    std::string current_line;

    std::vector<std::string> ret_list;
    ret_list.resize(2);

    std::ifstream file;

    file.open(file_name);

    if (!file) {
      std::cerr << "File " << file_name << " doesn't exist!";
      exit(1);
    }

    // first row is input directory
    getline(file, current_line);

    ret_list[0] = current_line;

    // second row is output directory
    getline(file, current_line);

    ret_list[1] = current_line;

    file.close();

    return ret_list;
  }

  Eigen::MatrixXd read_evaluation_positions(const std::string &file_name, const char del = ','){

    std::ifstream file;
    int numLines = 0;
    std::string current_line;
    std::string current_element;
    Eigen::MatrixXd ret_val;


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

    numLines = tmp_data.size();

    ret_val.resize(numLines,3);


    for (int l = 0; l < numLines; ++l){

        ret_val(l,0) = tmp_data[l][0];
        ret_val(l,1) = tmp_data[l][1];
        ret_val(l,2) = tmp_data[l][2];

      }

      return ret_val;

  }

  template<class Matrix>
  void write_binary(const char* filename, const Matrix& matrix){
      std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
      typename Matrix::Index rows=matrix.rows(), cols=matrix.cols();
      out.write((char*) (&rows), sizeof(typename Matrix::Index));
      out.write((char*) (&cols), sizeof(typename Matrix::Index));
      out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
      out.close();
  }

  template<class Matrix>
  void read_binary(const char* filename, Matrix& matrix){
      std::ifstream in(filename, std::ios::in | std::ios::binary);
      typename Matrix::Index rows=0, cols=0;
      in.read((char*) (&rows),sizeof(typename Matrix::Index));
      in.read((char*) (&cols),sizeof(typename Matrix::Index));
      matrix.resize(rows, cols);
      in.read( (char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar) );
      in.close();
  }

  void memory_report(){

    unsigned long mem_tot,mem_tot_av,mem_phys,mem_phys_av;
    unsigned long pagesize = (unsigned long) sysconf(_SC_PAGESIZE);

    //unsigned long phys_pages = (unsigned long) sysconf(_SC_PHYS_PAGES);
    //unsigned long av_pages = (unsigned long) sysconf(_SC_AVPHYS_PAGES);
    unsigned long sys_av_pages = (unsigned long) get_phys_pages();
    unsigned long curr_av_pages = (unsigned long) get_avphys_pages();

    //mem_tot = phys_pages*pagesize/1e9;
    //mem_tot_av = av_pages*pagesize/1e9;
    mem_phys = sys_av_pages*pagesize/1e9;
    mem_phys_av = curr_av_pages*pagesize/1e9;

    //std::cout << "Total System RAM \t= " << mem_tot << " GB" << std::endl;
    //std::cout << "Total Available RAM \t= " << mem_tot_av << " GB" << std::endl;
    std::cout << "RAM \t\t= " << mem_phys << " GB" << std::endl;
    std::cout << "Available RAM \t\t= " << mem_phys_av << " GB" << std::endl;


  }

  double get_occupied_ram(){

    double mem_phys,mem_phys_av;
    unsigned long pagesize = (unsigned long) sysconf(_SC_PAGESIZE);

    //unsigned long phys_pages = (unsigned long) sysconf(_SC_PHYS_PAGES);
    //unsigned long av_pages = (unsigned long) sysconf(_SC_AVPHYS_PAGES);
    unsigned long sys_av_pages = (unsigned long) get_phys_pages();
    unsigned long curr_av_pages = (unsigned long) get_avphys_pages();

    //mem_tot = phys_pages*pagesize/1e9;
    //mem_tot_av = av_pages*pagesize/1e9;
    mem_phys = sys_av_pages*pagesize/1e9;
    mem_phys_av = curr_av_pages*pagesize/1e9;

    return mem_phys - mem_phys_av;

  }

}  // namespace Bembel

#endif
