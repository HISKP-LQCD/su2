// parse_input_file.hh
/**
 * @brief parsing input file routines
 * This file contains the following namespaces:
 * - YAML_parsing : parsing parameters from YAML nodes
 * - input_file_parsing : parsing sets of parameters from files
 */

#pragma once

#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/type_index.hpp>
#include <set>
#include <string>

#include "parameters.hh"
#include "yaml.h"

/**
 * @brief namespace of yaml parsing functions
 * This namespace contains functions that read parameters specified in the input file from
 * which YAML nodes are generated
 */
namespace YAML_parsing {

  /**
   * @brief inspect the input file
   * Reading and initializing input parameters from the YAML file,
   * and checking if the provided ones coincide with the ones initialized by the program.
   * This is done via std::set<std::string> objects, to be checked at the end of the
   * parsing.
   *
   * At the beginning, all keys are stored in a std::set G (given). Nesting is rendered
   * naming the std::strings in it using ":". Each time a read function is called, a
   * parameter string is added to the first set U (used). If in the input file more
   * parameters have been provided than the ones initialized, or some's names have been
   * miswritten, an error is reported and the program aborts.
   *
   * The attribute InnerTree is a std::vector<std::string> containing the list of
   * nested nodes identifiers, such that only those submodes will be parsed and the
   * structure of the given (G) and used (U) parameters is preserved.
   * Example: {"a","b","c"} --> only nodes out of "a:b:c" will we parsed by the read()
   * function. This is useful when we one wnats to write functions specific to a given
   * block structure in the YAML input file. In that case, one should set the
   * InnerTree variable before the parsing (and pass the 'inspect_node' object by
   * reference)
   */
  class inspect_node {
  private:
    YAML::Node InnerNode; // node analyzed
    std::set<std::string> G, U; // set of given (G) and used (U) parameters
    std::vector<std::string> InnerTree = {}; // node from which we start reading
    YAML::Node OuterNode; // node analyzed, already prepended according to InnerTree

    /**
     * @brief fills the std::set (*this).G
     * Finds all nodes and subnodes of (*this).InnerNode and stores the keys in (*this).G
     * @param beg beginning of the string key of the node analyzed
     */
    void find_all(const YAML::Node &node, const std::string &beg) {
      for (YAML::const_iterator it = node.begin(); it != node.end(); ++it) {
        std::string k = it->first.as<std::string>(); // key

        if (G.count(beg + k) == 1) {
          std::cerr
            << "Error: The input file contains 2 or more identical nodes with name "
            << beg + k << "\n";
          std::abort();
        }

        G.insert(beg + k);
        // This does the recursion: if there is no further subnode to beg+k, the function
        // does nothing
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
    inspect_node(const YAML::Node &n) {
      (*this).InnerNode = YAML::Clone(n);
      this->find_all((*this).InnerNode, "");
    }

    YAML::Node get_outer_node(const std::vector<std::string> &tree) const {
      YAML::Node node_i = YAML::Clone((*this).InnerNode);
      const size_t nt = tree.size();
      for (int i = 0; i < nt; ++i) {
        node_i = YAML::Clone(node_i[tree[i]]);
      }
      return node_i;
    }

    YAML::Node get_outer_node() const { return this->get_outer_node((*this).InnerTree); }

    std::vector<std::string> get_InnerTree() const { return (*this).InnerTree; }

    YAML::Node get_root() { return YAML::Clone((*this).InnerNode); }

    void set_InnerTree(const std::vector<std::string> &pt) { (*this).InnerTree = pt; }

    /**
     * @brief reads the current state (branch of the YAML node) and goes inside the nested
     * node specified by the relative path pt1
     *
     */
    void dig_deeper(const std::vector<std::string> &pt1) {
      std::vector<std::string> pt2 = this->get_InnerTree();
      pt2.insert(pt2.end(), pt1.begin(), pt1.end());
      this->set_InnerTree(pt2);
    }

    /**
     * @brief reading a nested YAML node
     * @param outer_tree vector of strings determining the path to the input value
     */
    YAML::Node read_node(const std::vector<std::string> &outer_tree) {
      std::vector<std::string> tree = (*this).InnerTree;
      tree.insert(tree.end(), outer_tree.begin(), outer_tree.end());
      const YAML::Node node_i = YAML::Clone(this->get_outer_node(tree));
      const std::string g_str = this->get_node_str(tree);

      if (node_i.IsNull()) {
        std::cerr << "Error: \"" << g_str << "\" not found in your YAML. ";
        std::abort();
      }
      // adding the node string identifiers to the std::set (*this).U
      std::vector<std::string> t1 = {};
      const size_t n = tree.size();
      for (size_t i = 0; i < n; ++i) {
        t1.push_back(tree[i]);
        U.insert(this->get_node_str(t1));
      }

      return node_i;
    }

    /**
     * @brief saving value from YAML node
     * Saving in 'x' the value specified in the YAML node under the key string 'name'.
     * If the key doesn't exist, nothing is done
     * @tparam T cast type of the parameter
     * @param x address of the value
     * @param outer_tree vector of strings determining the path to the input value
     */
    template <class T> void read(T &x, const std::vector<std::string> &outer_tree) {
      YAML::Node node_i = YAML::Clone(this->read_node(outer_tree));

      std::vector<std::string> tree = (*this).InnerTree;
      tree.insert(tree.end(), outer_tree.begin(), outer_tree.end());
      const std::string g_str = this->get_node_str(tree);

      try {
        x = node_i.as<T>();
      } catch (...) {
        std::cerr << "Error: check \"" << g_str << "\" in your YAML input file. ";
        std::cerr << boost::typeindex::type_id<T>() << " type was expected. \n";
        std::abort();
      }
    }

    std::vector<std::string> get_full_tree(const std::vector<std::string> &tree) {
      std::vector<std::string> full_tree = (*this).InnerTree;
      full_tree.insert(full_tree.end(), tree.begin(), tree.end());
      return full_tree;
    }

    // read() and output on std::cout
    template <class T> void read_verb(T &x, const std::vector<std::string> &tree) {
      this->read<T>(x, tree);
      std::vector<std::string> tree2 = this->get_full_tree(tree);
      std::cout << "## " << this->get_node_str(tree2) << "=" << x << "\n";
      return;
    }

    // read() and output on std::cout each component of the list (passed in yaml format)
    template <class T>
    void read_sequence_verb(std::vector<T> &x, const std::vector<std::string> &tree) {
      const YAML::Node node = YAML::Clone(this->read_node(tree));
      std::string s1 = node.as<std::string>(); // = "n1, n2, ..."
      boost::replace_all(s1, " ", ""); // "n1,n2,n3,..." : no spaces
      std::vector<std::string> vs;
      boost::split(vs, s1, boost::is_any_of(",")); // vs = {"n1", "n2", ...}

      const size_t N = vs.size(); // number of elements in the sequence
      x.resize(N); // resizing the container

      std::vector<std::string> tree2 = this->get_full_tree(tree);
      std::cout << "## " << this->get_node_str(tree2) << "={";
      for (size_t i = 0; i < N; i++) {
        if (i > 0) {
          std::cout << ", ";
        }
        x[i] = boost::lexical_cast<T>(vs[i]);
        std::cout << x[i];
      }
      std::cout << "}\n";

      return;
    }

    /*
    Reading optional argument: same as read(), but if the key doesn't exist, nothing is
    done This function should be used with parameters that have a default argument.
    */
    template <class T> void read_opt_verb(T &x, const std::vector<std::string> &tree) {
      std::vector<std::string> tree1 = this->get_full_tree(tree);
      std::vector<std::string> tree2 = tree1;
      tree2.pop_back();
      if (this->get_outer_node(tree2)[tree.back()]) {
        this->read_verb<T>(x, tree);
      } else {
        std::cout << "## " << this->get_node_str(tree1) << "=" << x << " (default)\n";
      }
      return;
    }

    // read_verb() with `optional` as a parameter
    template <class T>
    void read_verb(const bool &optional, T &x, const std::vector<std::string> &tree) {
      if (optional) {
        this->read_opt_verb<T>(x, tree);

      } else {
        this->read_verb<T>(x, tree);
      }
    }

    /**
     * @brief finalize the parsing
     * Check if after all the read functions calls the std::set G and U are the same
     */
    void finalize() {
      if (G != U) {
        size_t ierr = 0;
        for (auto it = G.begin(); it != G.end(); ++it) {
          const size_t c = U.count(*it);
          if (c == 0) {
            std::cerr << "Error: You gave " << *it
                      << " in the input file, but is never used.\n";
            ierr = 1;
          }
        }

        if (ierr == 1) {
          std::abort();
        }
      }
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
   * @brief command line parsing
   * check if "help" has been called or if input file hasn't been specified correctly
   * @param ac argc from main()
   * @param av argv from main()
   * @param input_file input file reference
   */
  void parse_command_line(int ac, char *av[], std::string &input_file);

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
  void validate_beta_str_width(const size_t &n);

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
      void parse_input_file(const YAML::Node &nd,
                            gp::physics &pparams,
                            gp::hmc_u1 &hmc_params);
    }

    namespace measure {
      void parse_input_file(const YAML::Node &nd,
                            gp::physics &pparams,
                            gp::measure_u1 &mparams);

    } // namespace measure

    namespace metropolis {
      void validate_N_hit(const size_t &n);
      void parse_input_file(const YAML::Node &nd,
                            gp::physics &pparams,
                            gp::metropolis_u1 &mcparams);

    } // namespace metropolis

  } // namespace u1

} // namespace input_file_parsing
