#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"wilsonloop.hh"
#include"md_update.hh"
#include"monomial.hh"
#include"gradient_flow.hh"
#include"energy_density.hh"
#include"parse_input_file.hh"
#include"version.hh"
#include"smearape.hh"

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>
#include<algorithm>

namespace po = boost::program_options;

#include <boost/filesystem.hpp>

using std::cout;
using std::endl;


int main(int ac, char* av[]) {
  cout << "## Measuring Tool for U(1) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;

  namespace gp = global_parameters;
  gp::general gparams; // general parameters
  gp::measure_u1 mparams; // measure parameters

  std::string input_file; // yaml input file path
  po::options_description desc("Allowed options");
  desc.add_options()
  ("help,h", "produce this help message")
  ("file,f", po::value<std::string>(&input_file)->default_value("NONE"), "yaml input file");

  po::variables_map vm;
  po::store(po::parse_command_line(ac, av, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << desc << "\n";
    return 0;
  }

  namespace in_meas = input_file_parsing::u1::measure;
  int err = in_meas::parse_input_file(input_file, gparams, mparams);
  if (err > 0) {
    return err;
  }
  
  boost::filesystem::create_directories(boost::filesystem::absolute(mparams.confdir));


  gaugeconfig<_u1> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
  
  //needed for measuring potential
  std::ofstream resultfile;
  char filename[200];
  double loop;
  size_t maxsizenonplanar = (gparams.Lx < 4) ? gparams.Lx : 4;
  
  // write explanatory headers into result-files
  if(mparams.potential) {
      //~ open file for saving results
    
    if(gparams.ndims == 2){
      std::cout << "Currently not working for dim = 2, aborting" << std::endl;
      return 0;
    }
    
    //~ print heads of columns: W(r, t), W(x, y)
    if(!mparams.append && (gparams.ndims == 3 || gparams.ndims == 4)){
      if(gparams.ndims == 3){
        sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      if(gparams.ndims == 4){
        sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      resultfile.open(filename, std::ios::out);
      resultfile << "##";
      for (size_t t = 1 ; t <= gparams.Lt*mparams.sizeloops ; t++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
          resultfile << "W(x=" << x << ",t=" << t << ",y=" << 0 << ")  " ;
          }
      }
      resultfile << "counter";
      resultfile << std::endl; 
      resultfile.close();
      if(gparams.ndims == 3){
        sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      if(gparams.ndims == 4){
        sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      resultfile.open(filename, std::ios::out);
      resultfile << "##";
      for (size_t y = 1 ; y <= gparams.Ly*mparams.sizeloops ; y++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
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
    
    if(gparams.ndims == 2 || gparams.ndims == 4){
      std::cout << "Currently not working for dim = 2 and dim = 4, aborting" << std::endl;
      return 0;
    }
    
    //~ print heads of columns
    if(!mparams.append && (gparams.ndims == 3)){
      if(gparams.ndims == 3){
        sprintf(filename, "result2p1d.u1potential.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fnonplanar",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      resultfile.open(filename, std::ios::out);
      resultfile << "##";
      for (size_t t = 0 ; t <= gparams.Lt*mparams.sizeloops ; t++){
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
  for(size_t i = mparams.icounter; i < mparams.nmeas*mparams.nstep+mparams.icounter; i+=mparams.nstep) {
    std::ostringstream os;
    os << mparams.confdir + "/config_u1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt 
      << ".b" << std::fixed << U.getBeta();
    if(gparams.anisotropic){
      os << ".x" << std::fixed << gparams.xi;
    }
    os << "." << i << std::ends;
    int ierrU =  U.load(os.str());
    if(ierrU == 1){ // cannot load gauge config
      continue;
    }
    
    double plaquette = gauge_energy(U);
    double density = 0., Q=0.;
    energy_density(U, density, Q);
    cout << "## Initial Plaquette: " << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 
    cout << "## Initial Energy density: " << density << endl;
    
    random_gauge_trafo(U, mparams.seed);
    plaquette = gauge_energy(U);
    energy_density(U, density, Q);
    cout << "## Plaquette after rnd trafo: " << std::scientific << std::setw(15) << plaquette/U.getVolume()/double(U.getNc())/6. << endl; 
    cout << "## Energy density: " << density << endl;
    
    if(mparams.Wloop) {
      std::ostringstream os;
      os << mparams.confdir + "/wilsonloop.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      os << ".dat" << std::ends;
      compute_all_loops(U, os.str());
    }
    if(mparams.gradient) {
      std::ostringstream os;
      os << mparams.confdir + "/gradient_flow.";
      auto prevw = os.width(6);
      auto prevf = os.fill('0');
      os << i;
      os.width(prevw);
      os.fill(prevf);
      gradient_flow(U, os.str(), mparams.tmax);
    }
    
    if(mparams.potential || mparams.potentialsmall){
        //smear lattice
      for (size_t smears = 0 ; smears < mparams.n_apesmear ; smears +=1){
        smearlatticeape(U, mparams.alpha, mparams.smearspacial);
      }
    if(mparams.potential) {
      //~ //calculate wilsonloops for potential
      if(gparams.ndims == 4){
        sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
        resultfile.open(filename, std::ios::app);
        for (size_t t = 1 ; t <= gparams.Lt*mparams.sizeloops ; t++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
            loop  = wilsonloop_non_planar(U, {t, x, 0, 0});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
        sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
        resultfile.open(filename, std::ios::app);
        for (size_t y = 1 ; y <= gparams.Ly*mparams.sizeloops ; y++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
            loop  = wilsonloop_non_planar(U, {0, x, y, 0});
            loop += wilsonloop_non_planar(U, {0, x, 0, y});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/2.0 << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
      }
      if(gparams.ndims == 3){
        sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
        resultfile.open(filename, std::ios::app);
        for (size_t t = 1 ; t <= gparams.Lt*mparams.sizeloops ; t++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
            loop  = wilsonloop_non_planar(U, {t, x, 0});
            //~ loop  += wilsonloop_non_planar(U, {t, 0, x});
            resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
          }
        }
        resultfile << i;
        resultfile << std::endl; 
        resultfile.close();
        sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
        resultfile.open(filename, std::ios::app);
        for (size_t y = 1 ; y <= gparams.Ly*mparams.sizeloops ; y++){
          for (size_t x = 1 ; x <= gparams.Lx*mparams.sizeloops ; x++){
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
      if(gparams.ndims == 3){
        sprintf(filename, "result2p1d.u1potential.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fnonplanar",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, mparams.n_apesmear, mparams.alpha);
      }
      resultfile.open(filename, std::ios::app);
      for (size_t t = 0 ; t <= gparams.Lt*mparams.sizeloops ; t++){
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
