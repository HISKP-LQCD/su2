// parse_input_file.hh
/**
 * @brief parsing input file routines
 * This file contains the following namespaces:
 * - YAML_parsing : parsing parameters from YAML nodes
 * - input_file_parsing : parsing sets of parameters from files
 */

#pragma once

#include <boost/type_index.hpp>
#include <set>

#include "parameters.hh"
#include "yaml.h"

/**
 * @brief namespace of yaml parsing functions
 * This namespace contains functions that read parameters specified in the input file from
 * which YAML nodes are generated
 */
namespace YAML_parsing {

  inline void print_all_nodes(const YAML::Node& node){
    for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
      std::string k = it->first.as<std::string>(); // key
      std::cout << "check " << k << " " << node[k] << "\n";
    }
  }


  /**
   * @brief inspect the input file
   * Reading and initializing input parameters from the YAML file,
   * and checking if the provided ones coincide with the ones initialized by the program.
   * This is done via std::set<std::string> objects, to be checked at the end of the
   * parsing.
   *
   * At the beginning, all keys are stored in a std::set G (given). Nesting is rendered
   * naming the std::strings in it using ":". Each time a read function is called, a
   * parameter string is added to the first set U (used) If in the input file more
   * parameters have been provided than the ones initialized, or some's names have been
   * miswritten, an error is reported and the program aborts.
   */
  class inspect_node {
  private:
    YAML::Node N; // node analyzed
    std::set<std::string> G, U; // set of given (G) and used (U) parameters

    /**
     * @brief fills the std::set (*this).G
     * Finds all nodes and subnodes of (*this).N and stores the keys in (*this).G
     * @param beg beginning of the string key of the node analyzed
     */
    void find_all(const YAML::Node &node, const std::string &beg) {
      for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
        std::string k = it->first.as<std::string>(); // key
        G.insert(beg + k);
        this->find_all(node[k], beg + k + ":");
      }
    }

    /**
     * @brief Get the node str object
     * Example: tree = {"a", "b", "c"} -> "a:b:c"
     * @param tree vector of node names: innermost to outermost
     * @return std::string
     */
    std::string get_node_str(const std::vector<std::string> &tree) const {
      std::string g_str = "";

      const size_t nt = tree.size();
      for (int i = 0; i < nt; ++i) {
        g_str += tree[i] + ":";
      }
      g_str.pop_back(); // remove last character ":"
      return g_str;
    }

  public:

    inspect_node(const YAML::Node& n)
    {
      N = n;
      this->find_all((*this).N, "");
    }

    YAML::Node get_outer_node(const std::vector<std::string> &tree) const {
      YAML::Node node_i = (*this).N;      
      const size_t nt = tree.size();
      for (int i = 0; i < nt; ++i) {
        node_i = node_i[tree[i]];
      }
      return node_i;
    }

    YAML::Node get_root(){ return (*this).N; }

    /**
     * @brief saving value from YAML node
     * Saving in 'x' the value specified in the YAML node under the key string 'name'.
     * If the key doesn't exist, nothing is done
     * @tparam T cast type of the parameter
     * @param x address of the value
     * @param nd YAML node
     * @param name string name of the parameter
     */
    template <class T> void read(T &x, const std::vector<std::string> &tree) {
      const YAML::Node node_i = this->get_outer_node(tree);
      const std::string g_str = this->get_node_str(tree);

      try {
        x = node_i.as<T>();
        G.insert(g_str);
      } catch (...) {
        std::cerr << "Error: check \"" << g_str << "\" in your YAML input file. ";
        std::cerr << boost::typeindex::type_id<T>() << " type was expected. \n";
        abort();
      }
      return;
    }

    // read() and output on std::cout
    template <class T> void read_verb(T &x, const std::vector<std::string> &tree) {
      this->read<T>(x, tree);
      std::cout << "## " << this->get_node_str(tree) << "=" << x << "\n";
      return;
    }

    /*
    Reading optional argument: same as read(), but if the key doesn't exist, nothing is
    done This function should be used with parameters that have a default argument.
    */
    template <class T> void read_optional(T &x, const std::vector<std::string> &tree) {
      if (this->get_outer_node(tree)) {
        this->read<T>(x, tree);
      }
      return;
    }

    // read_optional() and output on std::cout
    template <class T> void read_opt_verb(T &x, const std::vector<std::string> &tree) {
      if (this->get_outer_node(tree)) {
        this->read_verb<T>(x, tree);
      } else {
        std::cout << "## " << this->get_node_str(tree) << "=" << x << " (default)\n";
      }
      return;
    }
  };

} // namespace YAML_parsing

/**
 * @brief parsing from input files
 * This namespace defines functions that, starting from a given input file, initialize the
 * attributes of structures needed for the various programs. The structures are passed as
 * reference to the functions. Their structure is:
 * @param file input file path
 * @param pparams physics parameter structure
 * @param xxx_params xxx program-specific parameter structure (e.g. xxx=hmc_params)
 * @return int error value
 */
namespace input_file_parsing {

  namespace gp = global_parameters;

  /**
   * @brief check geometry parameters
   * Check if the lattice spacetime extensions are sensible and consistent with the number
   * of dimensions (d). It is checked 2<=d. For 2<d<4, spatial dimensions of 'pparams' are
   * flattened according to the following priority order: Lz, Ly, Lx
   * @param file
   * @param pparams
   * @return int
   */
  int validate_geometry(gp::physics &pparams);

  namespace u1 {

    /**
     * @brief parsing geometrical parameters
     * Parsing the spacetime lattice extensions and number of dimensions from the input
     * file. If invalid, the program aborts.
     *
     * The provided node shoud be already the one named "geometry" in the input file
     * @param pparams physics parameter structure
     */
    void parse_geometry(const YAML::Node &nd, gp::physics &pparams);

    namespace hmc {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::hmc_u1 &hmc_params);
    }

    namespace measure {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::measure_u1 &mparams);

    } // namespace measure

    namespace metropolis {
      int parse_input_file(const std::string &file,
                           gp::physics &pparams,
                           gp::metropolis_u1 &mcparams);

    } // namespace metropolis

  } // namespace u1

} // namespace input_file_parsing
