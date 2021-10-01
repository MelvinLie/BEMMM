// This file is part of Bembel, the higher order C++ boundary element library.
// It was written as part of a cooperation of J. Doelz, H. Harbrecht, S. Kurz,
// M. Multerer, S. Schoeps, and F. Wolf at Technische Universtaet Darmstadt,
// Universitaet Basel, and Universita della Svizzera italiana, Lugano. This
// source code is subject to the GNU General Public License version 3 and
// provided WITHOUT ANY WARRANTY, see <http://www.bembel.eu> for further
// information.
#ifndef _BEMBEL_INPUTPROPERTIES_
#define _BEMBEL_INPUTPROPERTIES_


namespace Bembel {
namespace IO {

/**
 * \ingroup IO
 * \brief A simple class
 *
 */

class InputParameter {

 public:

   InputParameter(const std::string parameter_file) {

     std::ifstream file;
     std::string current_line;
     std::string current_element;

     file.open(parameter_file);

     if (!file) {
       std::cerr << "File " << parameter_file << " doesn't exist!";
       exit(1);
     }

     //first two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //on cluster
     getline(file, current_line);
     if (std::atoi( current_line.c_str() ) == 0){
       on_cluster_ = false;
     }
     else{
       on_cluster_ = true;
     }

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //mesh filename
     getline(file, current_line);
     mesh_filename_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //meas filename
     getline(file, current_line);
     meas_filename_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //sensor props filename
     getline(file, current_line);
     sensor_props_filename_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //output directory
     getline(file, current_line);
     output_directory_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //polynomial degrees
     p_.clear();
     getline(file, current_line);
     std::stringstream poly_data(current_line);

     while(std::getline(poly_data,current_element,',')){
       p_.push_back(std::atoi(current_element.c_str()));
     }

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //refinement levels
     h_.clear();
     getline(file, current_line);
     std::stringstream refine_data(current_line);

     while(std::getline(refine_data,current_element,',')){
       h_.push_back(std::atoi(current_element.c_str()));
     }

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //output directory
     getline(file, current_line);
     validation_filename_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //solution directory
     getline(file, current_line);
     solution_filename_ = current_line;

     //next two lines are not informative
     getline(file, current_line);
     getline(file, current_line);

     //input directory
     getline(file, current_line);
     input_directory_ = current_line;


   }

   void print_parameters(){

     std::cout << "*********************************************" << std::endl;
     std::cout << "Input Parameters" << std::endl;
     std::cout << "*********************************************" << std::endl;

     std::cout << "on cluster:" << std::endl;
     if(on_cluster_) std::cout << "\tyes"  << std::endl;
     else std::cout << "\tno"  << std::endl;

     std::cout << "mesh filename:" << std::endl;
     std::cout << "\t" << mesh_filename_ << std::endl;

     std::cout << "measurement filename:" << std::endl;
     std::cout << "\t" << meas_filename_ << std::endl;

     std::cout << "sensor properties filename:" << std::endl;
     std::cout << "\t" << sensor_props_filename_ << std::endl;

     std::cout << "input directory:" << std::endl;
     std::cout << "\t" << input_directory_ << std::endl;

     std::cout << "output directory:" << std::endl;
     std::cout << "\t" << output_directory_ << std::endl;

     std::cout << "validation filename:" << std::endl;
     std::cout << "\t" << validation_filename_ << std::endl;

     std::cout << "solution filename:" << std::endl;
     std::cout << "\t" << solution_filename_ << std::endl;

     std::cout << "polynomial degrees:" << std::endl;
     for(int i = 0 ; i < p_.size() ; ++i) std::cout << "\t" <<  p_[i];
     std::cout << std::endl;

     std::cout << "refinement levels:" << std::endl;
     for(int i = 0 ; i < h_.size() ; ++i) std::cout << "\t" <<  h_[i];
     std::cout << std::endl;
     std::cout << "*********************************************" << std::endl;

   }

   //the parameters of this class are public
   std::string mesh_filename_;
   std::string meas_filename_;
   std::string sensor_props_filename_;
   std::string output_directory_;
   std::string validation_filename_;
   std::string solution_filename_;
   std::string input_directory_;
   std::vector<int> p_,h_;
   bool on_cluster_;


};
}  // namespace Util
}  // namespace Bembel
#endif
