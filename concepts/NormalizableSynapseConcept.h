/*************************************************************

*************************************************************/

#ifndef NORMALIZABLESYNAPSECONCEPT_H_
#define NORMALIZABLESYNAPSECONCEPT_H_

#include <concepts>
#include <type_traits>

/*
 *  \class NormalizableSynapseConcept
 *
 *  A normalizable synapse is a synapse that can be normalized.
 *
 *  A model of this concept must meet the requirements for SynapseConcept plus:
 *
 *  The following methods
 *  \li precission_t get_w_max() const
 *  \li precission_t get_weight() const
 *  \li void set_weight(precission_t value)
 */
template <typename T>
concept NormalizableSynapseConcept = requires(T s) {
    typename T::precission_t;
    { s.get_w_max() } -> std::convertible_to<typename T::precission_t>;
    { s.get_weight() } -> std::convertible_to<typename T::precission_t>;
    { s.set_weight(typename T::precission_t{}) } -> std::same_as<void>;
};

#endif /*NORMALIZABLESYNAPSECONCEPT_H_*/
