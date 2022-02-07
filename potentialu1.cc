#include"su2.hh"
#include"u1.hh"
#include"gaugeconfig.hh"
#include"gauge_energy.hh"
#include"random_gauge_trafo.hh"
#include"wilsonloop.hh"
#include"energy_density.hh"
#include"parse_commandline.hh"
#include"version.hh"
#include"smearape.hh"

#include<iostream>
#include<iomanip>
#include<sstream>
#include<vector>
#include<random>
#include<boost/program_options.hpp>

namespace po = boost::program_options;

using std::cout;
using std::endl;


int main(int ac, char* av[]) {
  general_params gparams;
  size_t nstep;
  bool append;
  size_t n_apesmear;
  double alpha;
  bool smearspacial;

  cout << "## Measuring Tool for U(1) gauge theory" << endl;
  cout << "## (C) Carsten Urbach <urbach@hiskp.uni-bonn.de> (2017)" << endl;
  cout << "## GIT branch " << GIT_BRANCH << " on commit " << GIT_COMMIT_HASH << endl << endl;  

  po::options_description desc("Allowed options");
  add_general_options(desc, gparams);
  // add measure specific options
  desc.add_options()
    ("append", po::value<bool>(&append)->default_value(false), "are measurements appended to an existing file, or should it be overwritten?")
    ("nstep", po::value<size_t>(&nstep)->default_value(1), "measure each nstep config")
    ("napesmears", po::value<size_t>(&n_apesmear)->default_value(0), "number of APE smearings done on the lattice before measurement")
    ("apealpha", po::value<double>(&alpha)->default_value(1.0), "parameter alpha for APE smearings")
    ("spacialsmear", po::value<bool>(&smearspacial)->default_value(false), "should smearing be done only for spacial links?")
    ;

  int err = parse_commandline(ac, av, desc, gparams);
  if(err > 0) {
    return err;
  }

  gaugeconfig<_u1> U(gparams.Lx, gparams.Ly, gparams.Lz, gparams.Lt, gparams.ndims, gparams.beta);
  
  //~ open file for saving results
  std::ofstream resultfile;
  char filename[200];
  
  if(gparams.ndims == 2){
    std::cout << "Currently not working for dim = 2, aborting" << std::endl;
    return 0;
  }
  
  //~ print heads of columns: W(r, t), W(x, y)
  if(!append && (gparams.ndims == 3 || gparams.ndims == 4)){
    if(gparams.ndims == 3){
      sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
    }
    if(gparams.ndims == 4){
      sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
    }
    resultfile.open(filename, std::ios::out);
    resultfile << "##";
    for (size_t t = 1 ; t <= gparams.Lt/2 ; t++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
        resultfile << std::setw(5) << "W(x=" << std::setw(2) << x << ", t=" << std::setw(2) << t << ", y=" << std::setw(2) << 0 << ")  " ;
        }
    }
    resultfile << "counter";
    resultfile << std::endl; 
    resultfile.close();
    if(gparams.ndims == 3){
      sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
    }
    if(gparams.ndims == 4){
      sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
    }
    resultfile.open(filename, std::ios::out);
    resultfile << "##";
    for (size_t y = 1 ; y <= gparams.Ly/2 ; y++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
        resultfile << std::setw(5) << "W(x=" << std::setw(2) << x << ", t=" << std::setw(2) << 0 << ", y=" << std::setw(2) << y << ")  " ;
        }
    }
    resultfile << "counter";
    resultfile << std::endl; 
    resultfile.close();
      
  }
  
  double loop;
  for(size_t i = gparams.icounter; i < gparams.N_meas*nstep+gparams.icounter; i += nstep) {
    std::ostringstream os; 
    os << "configu1." << gparams.Lx << "." << gparams.Ly << "." << gparams.Lz << "." << gparams.Lt << ".b" << std::fixed << U.getBeta() << ".x" << gparams.xi << "." << i << std::ends;
    err = U.load(os.str());
    if(err != 0) {
      return err;
    }
    //smear lattice
    for (size_t smears = 0 ; smears < n_apesmear ; smears +=1){
      smearlatticeape(U, alpha, smearspacial);
    }
    //~ //calculate wilsonloops for potential
    if(gparams.ndims == 4){
      sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
      resultfile.open(filename, std::ios::app);
      for (size_t t = 1 ; t <= gparams.Lt/2 ; t++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
          loop  = wilsonloop_non_planar(U, {t, x, 0, 0});
          resultfile << std::setw(14) << std::scientific << loop/U.getVolume() << "  " ;
        }
      }
      resultfile << i;
      resultfile << std::endl; 
      resultfile.close();
      sprintf(filename, "result3p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
      resultfile.open(filename, std::ios::app);
      for (size_t y = 1 ; y <= gparams.Ly/2 ; y++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
          loop  = wilsonloop_non_planar(U, {0, x, y, 0});
          loop += wilsonloop_non_planar(U, {0, x, 0, y});
          //~ loop += wilsonloop_non_planar(U, {0, y, x});
          resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/2.0 << "  " ;
        }
      }
      resultfile << i;
      resultfile << std::endl; 
      resultfile.close();
    }
    if(gparams.ndims == 3){
      sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%ffinedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
      resultfile.open(filename, std::ios::app);
      for (size_t t = 1 ; t <= gparams.Lt/2 ; t++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
          loop  = wilsonloop_non_planar(U, {t, x, 0});
          resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/1.0 << "  " ;
        }
      }
      resultfile << i;
      resultfile << std::endl; 
      resultfile.close();
      sprintf(filename, "result2p1d.u1potential.rotated.Nt%lu.Ns%lu.b%f.xi%f.nape%lu.alpha%fcoarsedistance",gparams.Lt, gparams.Lx,U.getBeta(), gparams.xi, n_apesmear, alpha);
      resultfile.open(filename, std::ios::app);
      for (size_t y = 1 ; y <= gparams.Ly/2 ; y++){
        for (size_t x = 1 ; x <= gparams.Lx/2 ; x++){
          loop  = wilsonloop_non_planar(U, {0, x, y});
          //~ loop += wilsonloop_non_planar(U, {0, y, x});
          resultfile << std::setw(14) << std::scientific << loop/U.getVolume()/1.0 << "  " ;
        }
      }
      resultfile << i;
      resultfile << std::endl; 
      resultfile.close();
    }
  }
  return(0);
}

