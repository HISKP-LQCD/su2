/**
 * @file error.hpp
 * @brief error messages
 */
#pragma once

#include <iostream>
#include <boost/filesystem.hpp>

inline void fatal_error(const std::string msg, char const *function_name) {
  std::cerr << "# FATAL ERROR.  " << function_name << "(): " << msg << "\nAborting.";
  std::cerr << "# Aborting\n";
  std::abort();
  return;
}

inline void check_file_exists(const std::string &f, const std::string& func_name) {
  namespace fsys = boost::filesystem;
  if (!fsys::exists(fsys::absolute(f))) {
    std::cerr << "Error from " << func_name << "\n";
    std::cerr << f << ": no such file or directory. Aborting.\n";
    std::abort();
  }
}