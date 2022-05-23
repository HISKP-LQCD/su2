// output.hh
/**
 * @file output.hh
 * @author Simone Romiti (simone.romiti@uni-bonn.de)
 * @brief routines for outputting data
 * @version 0.1
 * @date 2022-05-09
 *
 * @copyright Copyright (c) 2022
 *
 */
#pragma once

#include <iomanip>
#include <sstream>
#include <string>

#include "parameters.hh"

namespace output {
  namespace gp = global_parameters;

  /**
   * @brief Get the beginning of the configuration path's string
   * This function returns the sub-string containing the path to the configuration,
   * without the information about the index of the sweep(Metropolis)/trajectory(HMC),etc.
   * It is assumed that 'sparams' has the attributes:
   * - conf_dir
   * - conf_basename
   * - beta_str_width
   * @tparam S type of structure containing the conf_basename attribute
   * @param pparams physical parameters
   * @param sparams parameters specific to the given
   * @return std::string
   */
  template <class S>
  std::string get_conf_path_basename(const gp::physics &pparams, const S &sparams) {
    // set up name for configs
    const std::string conf_basename = sparams.conf_basename;
    std::stringstream ss;
    ss << sparams.conf_basename << ".";
    ss << pparams.Lx << "." << pparams.Ly << "." << pparams.Lz << "." << pparams.Lt;
    if (pparams.rotating_frame) {
      ss << ".Omega_" << pparams.Omega;
    }
    ss << ".b" << std::fixed << std::setprecision(sparams.beta_str_width) << pparams.beta;
    if (pparams.anisotropic) {
      ss << ".x" << std::fixed << std::setprecision(sparams.beta_str_width) << pparams.xi;
    }

    return sparams.conf_dir+ss.str();
  }
  
  namespace measure {
    
    std::string get_filename_fine(const gp::physics &pparams,
                                  const gp::measure_u1 &mparams) {
      std::ostringstream filename_fine;
    
      filename_fine << mparams.resdir << "/"
                    << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt"
                    << pparams.Lt << ".Ns" << pparams.Lx << ".b" << std::fixed
                    << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi"
                    << std::fixed << std::setprecision(mparams.beta_str_width) << pparams.xi
                    << ".nape" << mparams.n_apesmear << ".alpha" << std::fixed
                    << mparams.alpha << "finedistance" << std::ends;
    
      return filename_fine.str();
    }
    
    std::string get_filename_coarse(const gp::physics &pparams,
                                    const gp::measure_u1 &mparams) {
      std::ostringstream f;
    
      f << mparams.resdir << "/"
        << "result" << pparams.ndims - 1 << "p1d.u1potential.rotated.Nt" << pparams.Lt
        << ".Ns" << pparams.Lx << ".b" << std::fixed
        << std::setprecision(mparams.beta_str_width) << pparams.beta << ".xi" << std::fixed
        << std::setprecision(mparams.beta_str_width) << pparams.xi << ".nape"
        << mparams.n_apesmear << ".alpha" << std::fixed << mparams.alpha << "coarsedistance"
        << std::ends;
    
      return f.str();
    }
    
    std::string get_filename_nonplanar(const gp::physics &pparams,
                                       const gp::measure_u1 &mparams) {
      std::ostringstream f;
    
      f << mparams.resdir << "/"
        << "result" << pparams.ndims - 1 << "p1d.u1potential.Nt" << pparams.Lt << ".Ns"
        << pparams.Lx << ".b" << std::fixed << std::setprecision(mparams.beta_str_width)
        << pparams.beta << ".xi" << std::fixed << std::setprecision(mparams.beta_str_width)
        << pparams.xi << ".nape" << mparams.n_apesmear << ".alpha" << std::fixed
        << mparams.alpha << "nonplanar" << std::ends;
    
      return f.str();
    }
    
    /**
     * @brief writes the headers for the results of the planar Wilson-Loops
     * in the files. coarse is the spacial direction, fine the temporal direction. a_t<=a_s.
     * @param pparams physical parameters
     * @param mparams parameters specific to measurement. May change if a non-valid dimension is given, so cannot be const.
     * @param filenamecoarse, filenamefine: names of the files where results are stored.
     * **/
    void set_header_planar(const gp::physics &pparams,
                           gp::measure_u1 &mparams,
                           const std::string &filename_coarse,
                           const std::string &filename_fine){
        //~ open file for saving results
      std::ofstream resultfile;
      if(pparams.ndims == 2){
        std::cerr << "Currently not working for dim = 2, no measurements for the potential will be made" << std::endl;
        mparams.potential = false;
      }
      
      //~ print heads of columns: W(r, t), W(x, y)
      if(!mparams.append && (pparams.ndims == 3 || pparams.ndims == 4)){
        resultfile.open(filename_fine, std::ios::out);
        resultfile << "## ";
        for (size_t t = 1 ; t <= pparams.Lt*mparams.sizeWloops ; t++){
            for (size_t x = 1 ; x <= pparams.Lx*mparams.sizeWloops ; x++){
            resultfile << "W(x=" << x << ",t=" << t << ",y=" << 0 << ")  " ;
            }
        }
        resultfile << "counter";
        resultfile << std::endl; 
        resultfile.close();
        
        resultfile.open(filename_coarse, std::ios::out);
        resultfile << "## ";
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
    
    /**
     * @brief writes the headers for the results of the non-planar Wilson-Loops in the file.
     * @param pparams physical parameters
     * @param mparams parameters specific to measurement. May change if a non-valid dimension is given, so cannot be const.
     * @param filenamenonplanar: name of the files where results are stored.
     * **/
    void set_header_nonplanar(const gp::physics &pparams,
                           gp::measure_u1 &mparams,
                           const std::string &filename_nonplanar){
        //~ open file for saving results
      std::ofstream resultfile;
      size_t maxsizenonplanar = (pparams.Lx < 4) ? pparams.Lx : 4;
    
      if(pparams.ndims == 2 || pparams.ndims == 4){
        std::cerr << "Currently not working for dim = 2 and dim = 4, no nonplanar measurements will be made" << std::endl;
        mparams.potentialsmall = false;
      }
      
      //~ print heads of columns
      if(!mparams.append && (pparams.ndims == 3)){
        resultfile.open(filename_nonplanar, std::ios::out);
        resultfile << "## ";
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
  } //namespace measure

  namespace hmc {

    std::string get_header(const std::string &sep = " ") {
      std::stringstream ss; // header: column names in the output
      ss << "i" << sep << "getaccept" << sep << "E*A" << sep << "dH" << sep << "rho"
         << sep << "ddH" << sep << "Q"
         << "\n";
      return ss.str();
    }

  } // namespace hmc

} // namespace output
