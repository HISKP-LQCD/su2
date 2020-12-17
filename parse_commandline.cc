#include<iostream>
#include<string>
#include<boost/program_options.hpp>

#include"parse_commandline.hh"
#include"print_program_options.hh"

namespace po = boost::program_options;


void add_general_options(po::options_description &desc, general_params &params) {
  desc.add_options()
    ("help,h", "produce this help message")
    ("spatialsizex,X", po::value<size_t>(&params.Lx), "spatial lattice size x")
    ("spatialsizey,Y", po::value<size_t>(&params.Ly), "spatial lattice size y")
    ("spatialsizez,Z", po::value<size_t>(&params.Lz), "spatial lattice size z")
    ("temporalsize,T", po::value<size_t>(&params.Lt), "temporal lattice size")
    ("ndims", po::value<size_t>(&params.ndims), "number of dimensions")
    ("nsave", po::value<size_t>(&params.N_save)->default_value(1000), "N_save")
    ("nmeas,n", po::value<size_t>(&params.N_meas)->default_value(10), "total number of sweeps")
    ("counter", po::value<size_t>(&params.icounter)->default_value(0), "initial counter for updates")
    ("beta,b", po::value<double>(&params.beta), "beta value")
    ("seed,s", po::value<size_t>(&params.seed)->default_value(13526463), "PRNG seed")
    ("heat", po::value<double>(&params.heat)->default_value(1.), "randomness of the initial config, 1: hot, 0: cold")
    ("restart", "restart from an existing configuration")
    ("configname", po::value< std::string >(&params.configfilename), "configuration filename used in case of restart")
    ;
  return;
}

int parse_commandline(int ac, char * av[], po::options_description &desc, general_params &params) {
  try {
    po::variables_map vm;
    po::store(po::parse_command_line(ac, av, desc), vm);
    po::notify(vm);

    params.restart = false;
    params.acceptreject = true;
    if (vm.count("help")) {
      std::cout << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizex") && !vm.count("help")) {
      std::cerr << "spatial lattice x-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizey") && !vm.count("help")) {
      std::cerr << "spatial lattice y-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("spatialsizez") && !vm.count("help")) {
      std::cerr << "spatial lattice z-size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("temporalsize") && !vm.count("help")) {
      std::cerr << "temporal lattice size must be given!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (!vm.count("beta") && !vm.count("help")) {
      std::cerr << "beta value must be specified!" << std::endl;
      std::cout << std::endl << desc << std::endl;
      return 1;
    }
    if (vm.count("restart")) {
      if(!vm.count("configname")) {
        std::cerr << "for a restart the configuration filename must be specified!" << std::endl;
        std::cout << std::endl << desc << std::endl;
        return 1;
      }
      params.restart = true;
    }
    if (vm.count("no-accept-reject")) {
      params.acceptreject = false;
    }
    if (!vm.count("ndims")) {
      params.ndims = 4;
    }
    PrintVariableMap(vm);
  }
  catch(std::exception& e) {
    std::cerr << "error: " << e.what() << "\n";
    return 1;
  }
  catch(...) {
    std::cerr << "Exception of unknown type!\n";
  }
  return(0);
}
