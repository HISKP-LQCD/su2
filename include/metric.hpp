/**
 * @file metric.hpp
 * @author Simone Romiti (simone.romiti.1994@gmail.com)
 * @brief classes defining metric types: Minkowski, rotating space, FLRW, etc.
 * @version 0.1
 * @date 2022-05-20
 *
 * @copyright Copyright (c) 2022
 *
 */

#pragma once

namespace metric {

  /**
   * @brief metric_type structure
   * structure allowing to specify the metric type template "at run time"
   * @tparam T
   */
  template <class T> struct metric_type { typedef T type; };

  class flat; // flat metric
  class rotating_z; // uniform rotation around the 'z' axis

} // namespace metric
