#include <fstream>
#include <iostream>
#include <string>

double get_available_RAM_Gb() {
  // ACHTUNG: it is assumed you are running on a Linux machine!
  std::ifstream meminfo("/proc/meminfo");

  if (!meminfo.is_open()) {
    std::cerr << "Error opening /proc/meminfo" << std::endl;
    std::abort();
  }

  std::string line;
  double g = 0;

  // first check from the physical RAM
  while (getline(meminfo, line)) {
    if (line.find("MemTotal:") != std::string::npos) {
      // Found the line containing total memory information
      const std::string totalMemoryString = line.substr(10);
      const double k = std::stoi(totalMemoryString); // RAM in kb
      g += k * (1e-6);
      break;
    }
  }

  // now check for the possible SWAP partition
  while (getline(meminfo, line)) {
    if (line.find("SwapTotal:") != std::string::npos) {
      // Found the line containing total memory information
      const std::string swapMemoryString = line.substr(10);
      const double k = std::stoi(swapMemoryString); // RAM in kb
      g += k * (1e-6);
      break;
    }
  }

  meminfo.close();

  return g;
}

void check_RAM_allocation(const size_t &size_in_Gb) {
  if (size_in_Gb > get_available_RAM_Gb()) {
    std::cerr << "## Available RAM " << get_available_RAM_Gb() << " Gb\n";
    std::cerr << "## Requested RAM " << size_in_Gb << " Gb\n";
    std::cerr << "## Error from " << __func__ << ":\n";
    std::cerr << "## Aborting. \n";
    std::abort();
  }
}
