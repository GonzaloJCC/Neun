/*************************************************************


*************************************************************/

#ifndef LINSKER_SYNAPSE_MODEL_H_
#define LINSKER_SYNAPSE_MODEL_H_

#include "SynapseWeightNormalizer.h"

/**
 * @brief Implements a synapsis based on (Linsker 1986)
 */
template <typename precission = double>
class LinskerSynapseModel {
#ifndef __AVR_ARCH__
  static_assert(std::is_floating_point<precission>::value);
#endif  //__AVR_ARCH__

 public:
    enum variable {
      w,        // Synaptic efficacy
      n_variables
    };

    enum parameter {
      v_pre,    // Presynaptic neuron's voltage (membrane potential)
      v_post,   // Postsynaptic neuron's voltage (membrane potential)
      xo,       // Constant for presynaptic neuron
      yo,       // Constant for postsynaptic neuron
      eta,      // Learning rate
      k1,       // Constant
      w_max,    // Maximum allowed synaptic weight (w_min = -wmax)
      i,        // Synaptic intensity
      n_parameters
    };

    typedef precission precission_t;

 public:
    LinskerSynapseModel() {}

    void eval(const precission* const vars, const precission* const params,
              precission* const incs) const {
        incs[w] = params[eta] * ((params[v_pre] - params[xo]) * (params[v_post] - params[yo]) + params[k1]);
    }
};

#endif /*LINSKER_SYNAPSE_MODEL_H_*/
