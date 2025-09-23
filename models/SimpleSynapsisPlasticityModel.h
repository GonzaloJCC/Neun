/*************************************************************


*************************************************************/

#ifndef SIMPLE_SYNAPSIS_PLASTICITY_MODEL_H_
#define SIMPLE_SYNAPSIS_PLASTICITY_MODEL_H_

#ifndef __AVR_ARCH__
#include <type_traits>
#endif  //__AVR_ARCH__
#include<cmath>

/**
 * @brief Implements a simple synapsis with plasticity
 */
template <typename precission = double>
class SimpleSynapsisPlasticityModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

 public:
    enum variable {mslow, n_variables};
    enum parameter {

      v_pre,
      previous_vpre,
      gfast,
      Esyn,
      sfast,
      Vfast,
      ifast,
      i,
      steps_for_degradation, /* Number of steps without synapsis needed for degradation*/
      steps_remaining, /* Current number of steps without synapsis remaining until degradation*/
      n_parameters
    };

    typedef precission precission_t;

 public:
    SimpleSynapsisPlasticityModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {
        incs[mslow] = 0; // this model does not use the slow components
    }
};

#endif  /*SIMPLE_SYNAPSIS_PLASTICITY_MODEL_H_*/
