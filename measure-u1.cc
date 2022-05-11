// measure-u1.cc
/**
 * @file measure-u1.cc
 * @author Carsten Urbach (urbach@hiskp.uni-bonn.de)
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief offline measurements of observables over previously generate gauge configurations
 * @version 0.1
 * @date 2022-05-11
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"omeasurements.hpp"
#include"polyakov_loop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"energy_density.hh"
#include"parse_input_file.hh"
#include"version.hh"
#include"smearape.hh"
#include"output.hh"

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>
#include<algorithm>

namespace po = boost::program_options;

#include <boost/filesystem.hpp>

int main(int ac, char* av[]) {
  std::cout << "## Measuring Tool for U(1) gauge theory" << std::endl;
  std::cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << std::endl;
  std::cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << std::endl;

  namespace gp = global_parameters;
  gp::physics pparams; // physics parameters
  gp::measure_u1 mparams; // measure parameters

  std::string input_file; // yaml input file path
  int err = input_file_parsing::parse_command_line(ac, av, input_file);
  if (err > 0) { return err; }

  namespace in_meas = input_file_parsing::u1::measure;
  err = in_meas::parse_input_file(input_file, pparams, mparams);
  if (err > 0) { return err; }
  
  boost::filesystem::create_directories(boost::filesystem::absolute(mparams.conf_dir));
  boost::filesystem::create_directories(boost::filesystem::absolute(mparams.resdir));


  gaugeconfig<_u1> U(pparams.Lx, pparams.Ly, pparams.Lz, pparams.Lt, pparams.ndims, pparams.beta);

  // get basename for configs
  std::string conf_path_basename = output::get_conf_path_basename(pparams, mparams);
  
  //filename needed for saving results from potential and potentialsmall
  const std::string filename_fine = output::get_filename_fine(pparams, mparams);
  const std::string filename_coarse = output::get_filename_coarse(pparams, mparams);
  const std::string filename_nonplanar = output::get_filename_nonplanar(pparams, mparams);  

  //needed for measuring potential
  std::ofstream resultfile;
  size_t maxsizenonplanar = (pparams.Lx < 4) ? pparams.Lx : 4;
  
  // write explanatory headers into result-files
  if(mparams.potential) {
      //~ open file for saving results
    if(pparams.ndims == 2){
      std::cout << "Currently not working for dim = 2, no measurements for the potential will be made" << std::endl;
      mparams.potential = false;
    }
    
    //~ print heads of columns: W(r, t), W(x, y)
    if(!mparams.append && (pparams.ndims == 3 || pparams.ndims == 4)){
      resultfile.open(filename_fine, std::ios::out);
      resultfile << "##";
      for (size_t t = 1 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
          resultfile << "W(x=" << x << ",t=" << t << ",y=" << 0 << ")  " ;
          }
      }
      resultfile << "counter";
      resultfile << std::endl; 
      resultfile.close();
      
      resultfile.open(filename_coarse, std::ios::out);
      resultfile << "##";
      for (size_t y = 1 ; y <= pparams.Ly*mparams.sizeWloops ; y++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
          resultfile << "W(x=" << x << ",t=" << 0 << ",y=" << y << ")  " ;
          }
      }
      resultfile << "counter";
      resultfile << std::endl; 
      resultfile.close();
        
    }
  }


  if(mparams.potentialsmall) {
      //~ open file for saving results
    
    if(pparams.ndims == 2 || pparams.ndims == 4){
      std::cout << "Currently not working for dim = 2 and dim = 4, no nonplanar measurements will be made" << std::endl;
      mparams.potentialsmall = false;
    }
    
    //~ print heads of columns
    if(!mparams.append && (pparams.ndims == 3)){
      resultfile.open(filename_nonplanar, std::ios::out);
      resultfile << "##";
      for (size_t t = 0 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
        for (size_t x = 0 ; x <= maxsizenonplanar ; x++){
          for (size_t y = 0 ; y <= maxsizenonplanar ; y++){
            resultfile << "W(x=" << x << ",t=" << t << ",y=" << y << ")  " ;
          }
        }
      }
      resultfile << "counter";
      resultfile << std::endl; 
      resultfile.close();        
    }
  }

/** 
 * do the measurements themselves: 
 * load each configuration, check for gauge invariance
 * if selected, measure Wilson-Loop and gradient flow
 * if chosen, do APE-smearing
 * if selected, then measure potential and small potential
 * the "small potential" are the nonplanar loops, but only with small extent in x, y
 * */

  const size_t istart = mparams.icounter==0? mparams.icounter + mparams.nstep : mparams.icounter; 
  for(size_t i = istart; i < mparams.n_meas*mparams.nstep+mparams.icounter; i+=mparams.nstep) {
    std::string path_i = conf_path_basename + "." + std::to_string(i);
    int ierrU =  U.load(path_i);
    if(ierrU == 1){ // cannot load gauge config
      continue;
    }
    
    double plaquette = gauge_energy(U);
    double density = 0., Q=0.;
    energy_density(U, density, Q);
    std::cout << "## Initial Plaquette: " << plaquette/U.getVolume()/double(U.getNc())/6. << std::endl; 
    std::cout << "## Initial Energy density: " << density << std::endl;
    
    random_gauge_trafo(U, mparams.seed);
    plaquette = gauge_energy(U);
    energy_density(U, density, Q);
    std::cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/double(U.getNc())/6. << std::endl; 
    std::cout << "## Energy density: " << density << std::endl;
        
    if(mparams.Wloop) {
      omeasurements::meas_wilson_loop<_u1>(U, i, mparams.conf_dir);
    }

    if(mparams.gradient) {
      omeasurements::meas_gradient_flow<_u1>(U, i, mparams.conf_dir, mparams.tmax);
    }
    
    if(mparams.potential || mparams.potentialsmall){
        //smear lattice
      for (size_t smears = 0 ; smears < mparams.n_apesmear ; smears +=1){
        smearlatticeape(U, mparams.alpha, mparams.smear_spatial_only);
      }
      double loop;
    if(mparams.potential) {
      //~ //calculate wilsonloops for potential
      if(pparams.ndims == 4){
        resultfile.open(filename_fine, std::ios::app);
        for (size_t t = 1 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
            loop  = wilsonloop_non_planar(U, {t, x, 0, 0});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
        
        resultfile.open(filename_coarse, std::ios::app);
        for (size_t y = 1 ; y <= pparams.Ly*mparams.sizeWloops ; y++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
            loop  = wilsonloop_non_planar(U, {0, x, y, 0});
            loop += wilsonloop_non_planar(U, {0, x, 0, y});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/2.0 << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
      }
      if(pparams.ndims == 3){
        resultfile.open(filename_fine, std::ios::app);
        for (size_t t = 1 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
            loop  = wilsonloop_non_planar(U, {t, x, 0});
            //~ loop  += wilsonloop_non_planar(U, {t, 0, x});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
        
        resultfile.open(filename_coarse, std::ios::app);
        for (size_t y = 1 ; y <= pparams.Ly*mparams.sizeWloops ; y++){
          for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
            loop  = wilsonloop_non_planar(U, {0, x, y});
            //~ loop += wilsonloop_non_planar(U, {0, y, x});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
      }
    }
    
    if(mparams.potentialsmall) {
      resultfile.open(filename_nonplanar, std::ios::app);
      for (size_t t = 0 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
        for (size_t x = 0 ; x <= maxsizenonplanar ; x++){
          for (size_t y = 0 ; y <= maxsizenonplanar ; y++){
            loop  = wilsonloop_non_planar(U, {t, x, y});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
      }
      resultfile << i;
      resultfile << std::endl; 
      resultfile.close();  
  }
  }

  }

  return(0);
}
