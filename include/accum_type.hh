#pragma once

/**
 * @brief generic accumulation type for elements of the group
 *
 * accumulations of group elements (e.g. sums of plaquettes) don't belong to the group
 * anymore. This the "type" attribute of this class defines a templated class representing
 * these quantities.
 *
 * By default, the same class used for the elements of the group is used. It is the case
 * of SU(2). However it becomes more expensive with U(1) and other SU(N) in general.
 *
 * The specialization of this templated struct have to be written
 * explicitly for these other groups.
 * See the implementation of U(1) and SU(3) accum types for an example.
 *
 * @tparam Group
 */
template <class Group> struct accum_type {
  typedef Group type;
};


