/**
 * @file error.hpp
 * @brief error messages
 */
#pragma once

#include <iostream>

inline void fatal_error(const std::string msg, char const *function_name) {
  std::cerr << "# FATAL ERROR.  " << function_name << "(): " << msg << "\nAborting.";
  std::cerr << "# Aborting\n";
  std::abort();
  return;
}
