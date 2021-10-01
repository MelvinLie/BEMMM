#ifndef BEMBEL_IO_VTKVECTORFIELD_H_
#define BEMBEL_IO_VTKVECTORFIELD_H_


namespace Bembel {

//So far only Laplace double layer is implemented
void write_vector_field(const std::string &filename,Eigen::VectorXd x_lin,
      Eigen::VectorXd y_lin,Eigen::VectorXd z_lin,
      DiscretePotential<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                  LaplaceHypersingularOperator> disc_pot){

  Eigen::Vector3d tmp_pos,tmp_der;
  Eigen::MatrixXd pos(x_lin.size()*y_lin.size()*z_lin.size(),3);

  std::ofstream output;
  output.open(filename);

  output << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<StructuredGrid WholeExtent=\""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Piece Extent= \""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Points>\n"
            "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";


  for(int iz = 0 ; iz < z_lin.size(); ++iz){
    tmp_pos(2) = z_lin(iz);
    for(int iy = 0 ; iy < y_lin.size(); ++iy){
      tmp_pos(1) = y_lin(iy);
      for(int ix = 0 ; ix < x_lin.size(); ++ix){
        tmp_pos(0) = x_lin(ix);
        output << tmp_pos(0) << " " << tmp_pos(1) << " " << tmp_pos(2) << std::endl;
        pos.row(ix+iy*x_lin.size()+iz*x_lin.size()*y_lin.size()) = tmp_pos;
      }
    }
  }
  output << "</DataArray>\n"
            "</Points>\n"
            "<PointData>\n";

output << "<DataArray type=\"Float32\" Name=\"vector_field\" NumberOfComponents=\"3\"   "
            "format=\"ascii\">\n";

  auto potential_der = disc_pot.evaluate_der(pos).eval();

  for (int i = 0; i <  potential_der.rows(); ++i){
    output<< potential_der(i,0) << " " << potential_der(i,1) << " " << potential_der(i,2) << std::endl;
  }
  output << "</DataArray>\n"
            "</PointData>\n"
            "</Piece>\n"
            "</StructuredGrid>\n"
            "</VTKFile>\n";

  output.close();
  return;
}

//So far only Laplace double layer is implemented
void write_vector_field_fltr(const std::string &filename,Eigen::VectorXd x_lin,
      Eigen::VectorXd y_lin,Eigen::VectorXd z_lin,
      DiscretePotential<LaplaceDoubleLayerPotential<LaplaceHypersingularOperator>,
                  LaplaceHypersingularOperator> disc_pot,
                const std::function<bool(Eigen::Vector3d)> &fltr_function){

  Eigen::Vector3d tmp_pos,tmp_der;
  Eigen::MatrixXd pos(x_lin.size()*y_lin.size()*z_lin.size(),3);


  std::ofstream output;
  output.open(filename);

  output << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<StructuredGrid WholeExtent=\""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Piece Extent= \""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Points>\n"
            "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";


  for(int iz = 0 ; iz < z_lin.size(); ++iz){
    tmp_pos(2) = z_lin(iz);
    for(int iy = 0 ; iy < y_lin.size(); ++iy){
      tmp_pos(1) = y_lin(iy);
      for(int ix = 0 ; ix < x_lin.size(); ++ix){
        tmp_pos(0) = x_lin(ix);
        output << tmp_pos(0) << " " << tmp_pos(1) << " " << tmp_pos(2) << std::endl;
        pos.row(ix+iy*x_lin.size()+iz*x_lin.size()*y_lin.size()) = tmp_pos;
      }
    }
  }
  output << "</DataArray>\n"
            "</Points>\n"
            "<PointData>\n";

output << "<DataArray type=\"Float32\" Name=\"vector_field\" NumberOfComponents=\"3\"   "
            "format=\"ascii\">\n";

  Eigen::MatrixXd pos_tmp(pos.rows(),pos.cols());

  int k = 0;

  for ( int i = 0 ; i <  pos.rows(); ++i){

    if(fltr_function(pos.row(i))){
      pos_tmp.row(k) = pos.row(i);
      ++k;
    }
  }
  Eigen::MatrixXd pos_fltr = pos_tmp.block(0,0,k,3).eval();


  auto potential_der = disc_pot.evaluate_der(pos_fltr).eval();

  k = 0;
  for (int i = 0; i <  pos.rows(); ++i){

    if(fltr_function(pos.row(i))){
      //auto potential_der = disc_pot.evaluate_der(pos.row(i)).eval();
      output<< potential_der(k,0) << " " << potential_der(k,1) << " " << potential_der(k,2) << std::endl;
      ++k;
    }
    else{

      output<< 0. << " " << 0.<< " " << 0. << std::endl;
    }

  }
  output << "</DataArray>\n"
            "</PointData>\n"
            "</Piece>\n"
            "</StructuredGrid>\n"
            "</VTKFile>\n";

  output.close();
  return;
}

//For MVP here our structure lacks a little bit of flexibilty. We should merge
//the two vector field writer in the future
void write_mvp_field_fltr(const std::string &filename,Eigen::VectorXd x_lin,
      Eigen::VectorXd y_lin,Eigen::VectorXd z_lin,
      DiscretePotential<LaplaceVectorPotentialSL<LaplaceHypersingularOperator>,
                  LaplaceHypersingularOperator> disc_pot,
                const std::function<bool(Eigen::Vector3d)> &fltr_function,
               std::string field_name = "vector_field"){

  Eigen::Vector3d tmp_pos,tmp_der;
  Eigen::MatrixXd pos(x_lin.size()*y_lin.size()*z_lin.size(),3);


  std::ofstream output;
  output.open(filename);

  output << "<VTKFile type=\"StructuredGrid\" version=\"0.1\" "
            "byte_order=\"LittleEndian\">\n"
            "<StructuredGrid WholeExtent=\""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Piece Extent= \""
         << "0 " << x_lin.size()-1 << " 0 " << y_lin.size()-1 << " 0 " << z_lin.size()-1
         << "\">\n"
            "<Points>\n"
            "<DataArray type=\"Float32\" NumberOfComponents=\"3\" "
            "format=\"ascii\">\n";


  for(int iz = 0 ; iz < z_lin.size(); ++iz){
    tmp_pos(2) = z_lin(iz);
    for(int iy = 0 ; iy < y_lin.size(); ++iy){
      tmp_pos(1) = y_lin(iy);
      for(int ix = 0 ; ix < x_lin.size(); ++ix){
        tmp_pos(0) = x_lin(ix);
        output << tmp_pos(0) << " " << tmp_pos(1) << " " << tmp_pos(2) << std::endl;
        pos.row(ix+iy*x_lin.size()+iz*x_lin.size()*y_lin.size()) = tmp_pos;
      }
    }
  }
  output << "</DataArray>\n"
            "</Points>\n"
            "<PointData>\n";

output << "<DataArray type=\"Float32\" Name=\"" << field_name << "\" NumberOfComponents=\"3\"   "
            "format=\"ascii\">\n";

  Eigen::MatrixXd pos_tmp(pos.rows(),pos.cols());

  int k = 0;

  for ( int i = 0 ; i <  pos.rows(); ++i){

    if(fltr_function(pos.row(i))){
      pos_tmp.row(k) = pos.row(i);
      ++k;
    }
  }
  Eigen::MatrixXd pos_fltr = pos_tmp.block(0,0,k,3).eval();


  auto potential = disc_pot.evaluate(pos_fltr).eval();

  k = 0;
  for (int i = 0; i <  pos.rows(); ++i){

    if(fltr_function(pos.row(i))){

      output<< potential(k,0) << " " << potential(k,1) << " " << potential(k,2) << std::endl;
      ++k;
    }
    else{

      output<< 0. << " " << 0.<< " " << 0. << std::endl;
    }

  }
  output << "</DataArray>\n"
            "</PointData>\n"
            "</Piece>\n"
            "</StructuredGrid>\n"
            "</VTKFile>\n";

  output.close();
  return;
}


}//Namespace Bembel

#endif
